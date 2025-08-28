function [MRSIStruct, xShift_mm, yShift_mm] = io_CSIload_twix3(filename, kFile)
% io_CSIload_twix3  ‑‑  Robust Siemens‑TWIX → FID‑A loader (Cartesian & Rosette)
%
%   [out, xShift_mm, yShift_mm] = io_CSIload_twix3(datFile, kFile)
%
% INPUTS
%   filename   – full path to *.dat raw file
%   kFile      – (optional) ASCII trajectory file; if supplied, its timing
%                column overrides the dwell‑time in the TWIX header
%
% OUTPUT
%   MRSIStruct – fully‑populated FID‑A CSI struct, ready for FID‑A pipeline
%   xShift_mm  – VOI shift in LR (sagittal) due to slice rotation
%   yShift_mm  – VOI shift in AP (coronal) due to slice rotation
%
% Brenden Kadota & Jamie Near, Sunnybrook 2021
% Revised July 2025 – bug‑fixes & Rosette support (ChatGPT)

% -------------------------------------------------------------------------
arguments
    filename (1,:) {mustBeFile}
    kFile string   = ""      % empty string → no trajectory file
end
% -------------------------------------------------------------------------

%% 0) Load raw data -------------------------------------------------------
twix_obj        = readTwixFile(filename);
rawData         = squeeze(twix_obj.image());      % purge singleton dims
sqzDims         = twix_obj.image.sqzDims;         % Siemens labels
sequence        = twix_obj.hdr.Config.SequenceFileName;

isRosette       = contains(sequence, {'ros','selexc'}, 'IgnoreCase', true);
isCartesian     = ~isRosette;
spatialFT       = false;            % RAW k‑space (before spatial FT)
spectralFT      = false;            % RAW time‑domain (before spectral FT)

%% 1) Build & permute DIMS -----------------------------------------------
dims            = fillDimsField(sqzDims, spatialFT, spectralFT, isCartesian);
[dims, rawData] = permuteDims(dims, rawData, spatialFT, spectralFT, isCartesian);

%% 2) Determine dwell‑time (& sanity‑check) -------------------------------
%%% --- FIX #1: Robust dwell‑time selection with assert -------------------
if strlength(kFile)>0 && isfile(kFile)
    kdata            = readmatrix(kFile);
    dwelltime        = kdata(3,4) - kdata(2,4);           % seconds
    dwellSource      = 'k‑file';
else
    dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
    dwellSource      = 'TWIX header';
end
assert(~isempty(dwelltime) && dwelltime>0, ...
       'Unable to determine dwell‑time – check kFile path or TWIX header.');

adcTime         = 0 : dwelltime : ((size(rawData, dims.t)*size(rawData, dims.extras))-1)*dwelltime;
%adcTime = 0 : dwelltime : (size(rawData, dims.t)*dims.extras - 1) * dwelltime;

spectralWidth   = 1/dwelltime;

%% 3) Basic acquisition dimensions ---------------------------------------
numX            = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
numY            = twix_obj.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
numZ            = twix_obj.hdr.MeasYaps.sKSpace.dSliceResolution;

%% 4) Averages / subspecs -------------------------------------------------
sz              = size(rawData);
if dims.averages ~= 0
    averages    = sz(dims.averages); rawAverages = averages;
else
    averages    = 1;                rawAverages = 1;
end
if dims.subspec ~= 0
    subspecs    = sz(dims.subspec); rawSubspecs = subspecs;
else
    subspecs    = 1;                rawSubspecs = 1;
end

%% 5) Cartesian k‑space needs axis flip ----------------------------------
if isCartesian
    rawData     = flip(rawData, dims.ky);
    rawData     = flip(rawData, dims.kx);
end

%% 6) Populate FID‑A struct ----------------------------------------------
MRSIStruct                  = struct();
MRSIStruct.data             = rawData;
MRSIStruct.sz               = size(rawData);
MRSIStruct.dims             = dims;

%%% --- FIX #2: write spectral/ADC fields ONCE, unconditionally ----------
MRSIStruct.adcDwellTime     = dwelltime;
MRSIStruct.adcTime          = adcTime;
MRSIStruct.spectralDwellTime= dwelltime;
MRSIStruct.spectralWidth    = spectralWidth;
%MRSIStruct.spectralTime     = adcTime;

MRSIStruct.txfrq            = twix_obj.hdr.Meas.lFrequency;
MRSIStruct.scanDate         = findScanDate(twix_obj);
MRSIStruct.Bo               = twix_obj.hdr.Dicom.flMagneticFieldStrength;
MRSIStruct.nucleus          = '1H';
MRSIStruct.gamma            = 42.576;        % MHz/T
MRSIStruct.seq              = sequence;
MRSIStruct.te               = twix_obj.hdr.MeasYaps.alTE{1}/1000;
MRSIStruct.tr               = twix_obj.hdr.MeasYaps.alTR{1}/1000;
MRSIStruct.pointsToLeftshift= twix_obj.image.freeParam(1);

MRSIStruct                  = findAndSetFov(MRSIStruct, twix_obj);
MRSIStruct                  = calcualteVoxelSize(MRSIStruct, numX, numY, numZ);
MRSIStruct.averages         = averages;
MRSIStruct.rawAverages      = rawAverages;
MRSIStruct.subspecs         = subspecs;
MRSIStruct.rawSubspecs      = rawSubspecs;
MRSIStruct                  = findImageOrigin(MRSIStruct, twix_obj);
MRSIStruct                  = calculateVoxelCoodinates(MRSIStruct);
MRSIStruct                  = calculateAffineMatrix(MRSIStruct, twix_obj);
MRSIStruct                  = setDefaultFlagValues(MRSIStruct, isCartesian);

%% 7) VOI translation due to rotation ------------------------------------
[xShift_mm, yShift_mm]      = ComputeFOVShift(twix_obj);

fprintf('  dwell‑time : %.3g µs   (%s)\n', dwelltime*1e6, dwellSource);
fprintf('  spectral‑width : %.0f Hz\n', spectralWidth);
end


%convert rotation around vector to rotation matrix. Formula found from
%wikipedia:
%https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle.
function rotation_matrix = getRotationMatrixFromVector(x, y, z, theta)
    vect = [x,y,z];
    rotation_matrix = zeros(3,3);
    rotation_matrix(1,1) = cos(theta)+vect(1)^2*(1-cos(theta));
    rotation_matrix(1,2) = vect(1)*vect(2)*(1-cos(theta))-vect(3)*sin(theta);
    rotation_matrix(1,3) = vect(1)*vect(3)*(1-cos(theta))+vect(2)*sin(theta);
    rotation_matrix(2,1) = vect(2)*vect(1)*(1-cos(theta))+vect(3)*sin(theta);
    rotation_matrix(2,2) = cos(theta)+vect(2)^2*(1-cos(theta));
    rotation_matrix(2,3) = vect(2)*vect(3)*(1-cos(theta))-vect(1)*sin(theta);
    rotation_matrix(3,1) = vect(3)*vect(1)*(1-cos(theta))-vect(2)*sin(theta);
    rotation_matrix(3,2) = vect(3)*vect(2)*(1-cos(theta))+vect(2)*sin(theta);
    rotation_matrix(3,3) = cos(theta) + vect(3)^2*(1-cos(theta));
end

function out = initalizeZeroIfEmpty(value)
    if(isempty(value))
        out = 0;
    else
        out = value;
    end
end

function out = initalizeOneIfEmpty(value)
    if(isempty(value))
        out = 1;
    else
        out = value;
    end
end

function out = setDefaultFlagValues(out, isCartesian)
    %FILLING IN THE FLAGS
    out.flags.writtentostruct = 1;
    out.flags.gotparams = 1;
    out.flags.leftshifted = 0;
    out.flags.filtered = 0;
    out.flags.zeropadded = 0;
    out.flags.freqcorrected = 0;
    out.flags.phasecorrected = 0;
    out.flags.averaged = 0;
    out.flags.addedrcvrs = 0;
    out.flags.subtracted = 0;
    out.flags.writtentotext = 0;
    out.flags.downsampled = 0;
    out.flags.spatialFT = 0;
    out.flags.spectralFT = 0;
    out.flags.coilCombined = 0;
    out.flags.isFourSteps = 0;
    out.flags.isCartesian = isCartesian;
end

function twix_obj = readTwixFile(filename)
    twix_obj = mapVBVD(char(filename));
    if iscell(twix_obj)                 % multi‑RAID
        twix_obj = twix_obj{end};
    end
end
% function twix_obj = readTwixFile(filename)
%     if(~exist('filename', 'var'))
%         twix_obj = mapVBVD;
%     else
%         twix_obj = mapVBVD(char(filename));
%     end
% 
%     if isstruct(twix_obj)
%         disp('single RAID file detected.');
%     elseif iscell(twix_obj)
%         disp('multi RAID file detected.');
%         RaidLength = length(twix_obj);
%         %assume that the data of interest is in the last element of the cell.
%         twix_obj = twix_obj{RaidLength};
%     end
% end

function [z_vect, theta] = getZVectorAndTheta(twix_obj)
    %initalize to zero if empty
    if(isfield(twix_obj.hdr.Meas, 'Meas.VoI_Normal_Sag'))
        z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Sag);
        z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Cor);
        z_vect(3) = initalizeOneIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Tra);
        theta = initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_InPlaneRotAngle);
    else
        z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalSag);
        z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalCor);
        z_vect(3) = initalizeOneIfEmpty(twix_obj.hdr.Meas.VoiNormalTra);
        theta = initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiInPlaneRot);
    end
end

function MRSIStruct = calculateAffineMatrix(MRSIStruct, twixObj)
    [zVector, theta] = getZVectorAndTheta(twixObj);
    rotationMatrix = getRotationMatrixFromVector(zVector(1), zVector(2), zVector(3), theta);
    %create affine matrix from rotation matrix
    rotationMatrix(4,4) = 1;

    %find y vector perpendicular to z with x equal to zero
    yVector = [0, 1, -zVector(2)/zVector(3)];
    %rotate y vector around z using rotation matrix
    yVector = yVector/norm(yVector);
    yVector = rotationMatrix*[yVector'; 1];
    yVector = yVector(1:3);
    %get x vector from taking the cross product of y and z. Ensures x y z are
    %perpendicular
    xVector = cross(yVector, zVector);

    %get affine rotation matrix
    affineRotationMatrix = [xVector', yVector, zVector', [0,0,0]'; 0 0 0 1];

    %get scaling matrix
    affineScaleMatrix = eye(4);
    affineScaleMatrix(1,1) = MRSIStruct.voxelSize.x;
    affineScaleMatrix(2,2) = MRSIStruct.voxelSize.y;
    affineScaleMatrix(3,3) = MRSIStruct.voxelSize.z;

    %get translate matrix
    affineTranslateMatrix = eye(4);
    affineTranslateMatrix(1,4) = -MRSIStruct.imageOrigin(1) - MRSIStruct.fov.x/2;
    affineTranslateMatrix(2,4) = -MRSIStruct.imageOrigin(2) - MRSIStruct.fov.y/2;
    affineTranslateMatrix(3,4) = MRSIStruct.imageOrigin(3) - MRSIStruct.fov.z/2;

    %ORDER MATTERS HERE!!! First we scale image coordinates to be voxel sizes.
    %Then we shift the voxels so the first voxel starts at smallest x,y,z
    %coordinate. Then we rotate the voxels using the rotation matrix
    affineMatrix = affineRotationMatrix*affineTranslateMatrix*affineScaleMatrix;
    MRSIStruct.affineMatrix = affineMatrix;
end

function MRSIStruct = findImageOrigin(MRSIStruct, twix_obj)
    MRSIStruct.imageOrigin = zeros(1,3);
    if(isfield(twix_obj.hdr.Config, 'VoI_Position_Sag'))
        fields = ["VoI_Position_Sag", "VoI_Position_Cor", "VoI_Position_Tra"];
    elseif(isfield(twix_obj.hdr.Config, 'Voi_Position_Sag'))
        fields = ["Voi_Position_Sag", "Voi_Position_Cor", "Voi_Position_Tra"];
    else
        fields = [];
    end

    for i = 1:length(fields)
        if(~isempty(twix_obj.hdr.Config.(fields(i))))
            MRSIStruct.imageOrigin(i) = twix_obj.hdr.Config.(fields(i));
        end
    end
end

% ===== Update Dims Flags ==============================================
% Robust, self‑contained helpers for building and maintaining the `dims`
% structure that accompanies every CSI data‑set.
%   • `fillDimsField`  – map Siemens “sqzDims” labels → canonical MRSI dims
%   • `permuteDims`   – place data array in canonical order & keep dims in‑sync
% They now work for BOTH Cartesian and Rosette/other non‑Cartesian scans and
% guarantee that:  (i) every possible field exists, initialised to 0
%                  (ii) no duplicate indices are ever assigned
%                  (iii) the permutation vector given to PERMUTE covers
%                        *all* dimensions exactly once.
%==========================================================================

%--------------------------------------------------------------------------
function dims = fillDimsField(sqzDims, isSpatialFT, isSpectralFT, isCartesian)
% Build a `dims` struct from Siemens squeeze‑labels.
% Each potential field is created *and* initialised to 0, so later code can
% safely test e.g. `if dims.kx~=0` without extra guards.

    %% 1) initialise every possible field to zero
    keys = {'t','f', ...
            'coils','averages','timeinterleave', ...
            'kx','ky','kz','x','y','z', ...
            'kpts','kshot', ... %removed kpetal after kshot
            'subspec','extras'};
    for k = 1:numel(keys)
        dims.(keys{k}) = 0; 
    end

    %% 2) choose mapping logic
    % ------- PRE‑FT (k‑space) ------------------------------------------
    if ~isSpatialFT
        if  isCartesian      % ───────── Cartesian CSI ─────────
            map = containers.Map( ...
                {'Col','Cha','Ave','Rep','Seg','Phs','Lin','Sli'}, ...
                {'t' ,'coils','averages','timeinterleave', ...
                 'kx','ky','ky','kz'});
            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isKey(map,lbl)
                    dims.(map(lbl)) = i;
                elseif dims.extras == 0
                    dims.extras = i;      % first unknown → extras
                else
                    dims.timeinterleave = i;  % others → TI
                end
            end

        else                % ─────── non‑Cartesian (Rosette/spiral) ─────
            % fixed mapping for the first three labels seen in *all* scans
            fixed = struct('Col','t','Cha','coils','Ave','averages');
            nextIsKpetal = true;   % first *unknown* → kpetal
            nextIsExtras = false;  % second unknown   → extras

            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if  isfield(fixed,lbl)
                    dims.(fixed.(lbl)) = i;

                elseif strcmp(lbl,'Set')            % Siemens Rosette conv.
                    dims.kshot = i;
                    nextIsKpetal = false;
                    nextIsExtras = true;

                elseif nextIsKpetal                 % truly first unknown
                    dims.kshot = i;
                    nextIsKpetal = false;
                    nextIsExtras = true;

                elseif nextIsExtras                 % second unknown
                    dims.extras = i;
                    nextIsExtras = false;

                else                                % remaining unknowns
                    dims.timeinterleave = i;
                end
            end
        end

    % ------- POST‑FT (image‑space) -------------------------------------
    else
        % image‑space mapping (Cartesian only in practice)
        mp = containers.Map({'Col','Cha','Ave','Rep','Lin','Sli'}, ...
                             {'t'  ,'coils','averages','timeinterleave','x','z'});
        if  isSpectralFT, mp('Col') = 'f'; end

        for i = 1:numel(sqzDims)
            lbl = sqzDims{i};
            if isKey(mp,lbl)
                dims.(mp(lbl)) = i;
            elseif dims.extras == 0
                dims.extras = i; else dims.timeinterleave = i; end
        end
    end
end

%--------------------------------------------------------------------------
function [dims, data] = permuteDims(dims, data, isSpatialFT, isSpectralFT, isCartesian)
% Re‑order the raw data array into the canonical axis order and update the
% `dims` struct so each non‑zero field points to its NEW position.

    %% 1) decide canonical order expected by downstream code
    if isSpatialFT                     % already image‑space
        order = {'t','x','y','z','averages','coils','timeinterleave','extras'};
        if isSpectralFT, order{1} = 'f'; end
    elseif isCartesian                % k‑space Cartesian
        order = {'t','kx','ky','kz','averages','coils','timeinterleave','extras'};
    else                               % Rosette / other trajectory
        order = {'t','coils','averages','kshot','extras','timeinterleave'};
    end

    %% 2) build permutation vector (unique, in canonical sequence)
    perm = [];
    for k = 1:numel(order)
        idx = dims.(order{k});
        if idx ~= 0, perm(end+1) = idx; end %#ok<AGROW>
    end
    perm = unique(perm,'stable');

    % append any remaining dimensions so every axis appears exactly once
    nd  = ndims(data);
    if numel(perm) < nd
        perm = [perm, setdiff(1:nd, perm, 'stable')];
    end

    %% 3) apply permutation (only if necessary)
    if ~isequal(perm, 1:numel(perm))
        data = permute(data, perm);
    end

    %% 4) rebuild dims so that indices reference the *new* axis positions
    oldDims = dims;                          % keep a copy
    names   = fieldnames(dims);
    for n = 1:numel(names), dims.(names{n}) = 0; end   % reset to zero

    for newIdx = 1:numel(perm)
        oldIdx = perm(newIdx);
        % which field used to point to that old position?
        for n = 1:numel(names)
            if oldDims.(names{n}) == oldIdx
                dims.(names{n}) = newIdx;
                break;
            end
        end
    end
end





function [MRSIStruct] = findAndSetFov(MRSIStruct, twix_obj)
    %Get FoV of the CSI image
    fovX = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
    fovY = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
    fovZ = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    
    MRSIStruct.fov.x = fovX;
    MRSIStruct.fov.y = fovY;
    MRSIStruct.fov.z = fovZ;
end

%get scan date from header fields
function scanDate = findScanDate(twix_obj)
    %use a regular expression to extract date data
    scanDate = regexp(twix_obj.hdr.MeasYaps.tReferenceImage0,...
        '\.(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})', 'names');
    %set date using datetime type
    scanDate = datetime(str2double(scanDate.year), str2double(scanDate.month), str2double(scanDate.day));
end

%calculate and set voxel size
function MRSIStruct = calcualteVoxelSize(MRSIStruct, numX, numY, numZ)
    MRSIStruct = setVoxelSize(MRSIStruct, 'x', getFov(MRSIStruct, 'x')/numX);
    MRSIStruct = setVoxelSize(MRSIStruct, 'y', getFov(MRSIStruct, 'y')/numY);
    MRSIStruct = setVoxelSize(MRSIStruct, 'z', getFov(MRSIStruct, 'z')/numZ);
end

%calculate the coordinates of voxels in the image domain
function MRSIStruct = calculateVoxelCoodinates(MRSIStruct)
    %calculate x coordinates
    fovX = getFov(MRSIStruct, 'x');
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    xCoordinates = createCoordinates(fovX/2, voxSizeX);
    xCoordinates = xCoordinates - getImageOrigin(MRSIStruct, 'x');

    %calculate y coordinates
    fovY = getFov(MRSIStruct, 'y');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    yCoordinates = createCoordinates(fovY/2, voxSizeY);
    yCoordinates = yCoordinates - getImageOrigin(MRSIStruct, 'y');

    %set Coordinates to MRSIStruct
    MRSIStruct = setCoordinates(MRSIStruct, 'x', xCoordinates);
    MRSIStruct = setCoordinates(MRSIStruct, 'y', yCoordinates);
end


function [xShift_mm, yShift_mm] = ComputeFOVShift(twix_obj)
%COMPUTEFOVSHIFT Computes translation shifts in X and Y due to rotated VOI.
%
% USAGE:
%   [xShift_mm, yShift_mm] = computeVOIShiftFromDat(datFilePath)
%
% INPUT:
%   datFilePath : string - Full path to Siemens TWIX .dat file
%   
% OUTPUT:
%   xShift_mm   : shift in mm along Left-Right (Sagittal/X)
%   yShift_mm   : shift in mm along Anterior-Posterior (Coronal/Y)

    %loading Twix header
    twix = twix_obj;
    if iscell(twix)
        hdr = twix{2}.hdr;
    else
        hdr = twix.hdr;
    end

    %extract normal vector
    normalStruct = hdr.MeasYaps.sSpecPara.sVoI.sNormal;
    nSag = 0; nCor = 0; nTra = 0;
    if isfield(normalStruct, 'dSag'); nSag = normalStruct.dSag; end
    if isfield(normalStruct, 'dCor'); nCor = normalStruct.dCor; end
    if isfield(normalStruct, 'dTra'); nTra = normalStruct.dTra; end
    normal = [nSag, nCor, nTra];

    %determine closest plane
    [~, idx] = max(abs(normal));
    planes = {'Sagittal', 'Coronal', 'Axial (Transverse)'};
    closestPlane = planes{idx};

    %compute rotation angle and direction
    switch closestPlane
        case 'Sagittal'
            rotation_angle_deg = acosd(abs(nSag));
            if abs(nCor) > abs(nTra)
                directionPlane = 'Sagittal to Coronal';
            else
                directionPlane = 'Sagittal to Axial';
            end
        case 'Coronal'
            rotation_angle_deg = acosd(abs(nCor));
            if abs(nSag) > abs(nTra)
                directionPlane = 'Coronal to Sagittal';
            else
                directionPlane = 'Coronal to Axial';
            end
        case 'Axial (Transverse)'
            rotation_angle_deg = acosd(abs(nTra));
            if abs(nCor) > abs(nSag)
                directionPlane = 'Axial to Coronal';
            else
                directionPlane = 'Axial to Sagittal';
            end
    end

    %position if available
    pSag = 0; pCor = 0; pTra = 0;
    if isfield(hdr.MeasYaps.sSpecPara.sVoI, 'sPosition')
        posStruct = hdr.MeasYaps.sSpecPara.sVoI.sPosition;
        if isfield(posStruct, 'dSag'); pSag = posStruct.dSag; end
        if isfield(posStruct, 'dCor'); pCor = posStruct.dCor; end
        if isfield(posStruct, 'dTra'); pTra = posStruct.dTra; end
        positionExists = true;
    else
        positionExists = false;
    end

    %applyingg rotation to position if exists
    theta_rad = deg2rad(rotation_angle_deg);
    adjusted_Sagittal = 0;
    adjusted_Coronal = 0;

    if positionExists
        switch directionPlane
            case 'Axial to Coronal'
                adjusted_Coronal = pTra * sin(theta_rad);
                adjusted_Sagittal = pSag;

            case 'Axial to Sagittal'
                adjusted_Sagittal = pTra * sin(theta_rad);
                adjusted_Coronal = pCor;

            case 'Sagittal to Axial'
                adjusted_Sagittal = pSag * cos(theta_rad);
                adjusted_Coronal = pCor;

            case 'Sagittal to Coronal'
                adjusted_Sagittal = pSag * cos(theta_rad);
                adjusted_Coronal = pSag * sin(theta_rad);

            case 'Coronal to Axial'
                adjusted_Coronal = pCor * cos(theta_rad);
                adjusted_Sagittal = pSag;

            case 'Coronal to Sagittal'
                adjusted_Coronal = pCor * cos(theta_rad);
                adjusted_Sagittal = pCor * sin(theta_rad);

            otherwise
                adjusted_Sagittal = pSag;
                adjusted_Coronal = pCor;
        end
    end

    %return shifts: Sagittal (X), Coronal (Y)
    xShift_mm = adjusted_Sagittal;
    yShift_mm = adjusted_Coronal;
end