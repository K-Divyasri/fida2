function [MRSIStruct, phaseMap, freqMap] = op_CSIB0Correction(MRSIStruct, phaseMap, freqMap, plottingArguments)
    arguments
        MRSIStruct (1, 1) struct
        phaseMap double = []
        freqMap double = []
        plottingArguments.isPlot = false;
    end

    if ~isfield(MRSIStruct.flags, 'addedrcvrs') || MRSIStruct.flags.addedrcvrs == 0
        error("Error: Coils must be combined first (addedrcvrs = 1).");
    end
    if ~MRSIStruct.flags.spatialFT
        error("Please Fourier transform along the spatial dimension.");
    end
    if ~MRSIStruct.flags.spectralFT
        error("Please Fourier transform along the spectral dimension.");
    end

    % Find voxel with highest intensity
    [x, y] = findMaxCoord(MRSIStruct);
    referenceMRS = op_CSItoMRS(MRSIStruct, x, y);
    referenceMRS = setFlags(referenceMRS, 'averaged', true);

    % Use provided maps or compute them
    if ~isempty(phaseMap) && ~isempty(freqMap)
        MRSIStruct = applyFreqMap(MRSIStruct, freqMap);
        MRSIStruct = applyPhaseMap(MRSIStruct, phaseMap);
    else
        sz = MRSIStruct.sz;
        phaseMap = zeros(sz(3), sz(2));  % [y, x]
        freqMap  = zeros(sz(3), sz(2));
        data = MRSIStruct.data;

        for y = 1:sz(3)
            for x = 1:sz(2)
                alignMRS = op_CSItoMRS(MRSIStruct, x, y);
                alignMRS = setFlags(alignMRS, 'averaged', true);
                [alignedMRS, phase, freq] = op_alignScans(alignMRS, referenceMRS, referenceMRS.t(end));

                data(:, y, x) = alignedMRS.specs;
                phaseMap(y, x) = phase;
                freqMap(y, x) = freq;
            end
        end

        MRSIStruct.data = data;
        if plottingArguments.isPlot
            plotFreqAndPhaseMap(phaseMap, MRSIStruct, freqMap);
        end
        MRSIStruct = setFlags(MRSIStruct, 'phasecorrected', 1);
        MRSIStruct = setFlags(MRSIStruct, 'freqcorrected', 1);
    end
end

function [x, y] = findMaxCoord(in)
    [~, i] = max(in.data(:));
    [~, y, x] = ind2sub(size(in.data), i);
end

function MRSIStruct = applyFreqMap(MRSIStruct, freqMap)
    for y = 1:size(freqMap, 1)
        for x = 1:size(freqMap, 2)
            alignMRS = op_CSItoMRS(MRSIStruct, x, y);
            alignMRS = op_freqshift(alignMRS, freqMap(y, x));
            MRSIStruct.data(:, y, x) = alignMRS.specs;
        end
    end
end

function MRSIStruct = applyPhaseMap(MRSIStruct, phaseMap)
    for y = 1:size(phaseMap, 1)
        for x = 1:size(phaseMap, 2)
            alignMRS = op_CSItoMRS(MRSIStruct, x, y);
            alignMRS = op_addphase(alignMRS, phaseMap(y, x));
            MRSIStruct.data(:, y, x) = alignMRS.specs;
        end
    end
end

function plotFreqAndPhaseMap(phaseMap, MRSIStruct, freqMap)
    figure;
    subplot(2, 1, 1);
    imagesc(MRSIStruct.coordinates.x, MRSIStruct.coordinates.y, phaseMap');
    title('Phase Map'); xlabel('x'); ylabel('y'); colorbar;

    subplot(2, 1, 2);
    imagesc(MRSIStruct.coordinates.x, MRSIStruct.coordinates.y, freqMap');
    title('Frequency Map'); xlabel('x'); ylabel('y'); colorbar;
end

% function [MRSIStruct, phaseMap, freqMap] = op_CSIB0Correction(MRSIStruct, varargin)
% % op_CSIB0Correction: B0 field inhomogeneity correction for CSI data
% %
% % USAGE:
% % [correctedStruct, phaseMap, freqMap] = op_CSIB0Correction(MRSIStruct)
% % [correctedStruct] = op_CSIB0Correction(MRSIStruct, phaseMap, freqMap)
% %
% % INPUTS:
% % MRSIStruct    = CSI data structure
% % phaseMap      = (optional) pre-computed phase map
% % freqMap       = (optional) pre-computed frequency map
% %
% % OUTPUTS:
% % MRSIStruct    = B0-corrected CSI data structure
% % phaseMap      = phase correction map
% % freqMap       = frequency correction map
% 
% % Check input arguments
% if nargin < 1
%     error('op_CSIB0Correction requires at least one input argument');
% end
% 
% % Initialize output variables
% phaseMap = [];
% freqMap = [];
% 
% % Check if pre-computed maps are provided
% if nargin >= 3
%     phaseMap = varargin{1};
%     freqMap = varargin{2};
%     applyCorrection = true;
%     computeMaps = false;
% elseif nargin == 2
%     error('Both phaseMap and freqMap must be provided together');
% else
%     applyCorrection = false;
%     computeMaps = true;
% end
% 
% % Validate input structure
% if ~isstruct(MRSIStruct) || ~isfield(MRSIStruct, 'sz') || ~isfield(MRSIStruct, 'dims')
%     error('Invalid MRSIStruct input. Must contain sz and dims fields.');
% end
% 
% % Check for data field (can be either 'data' or 'fids')
% if isfield(MRSIStruct, 'data')
%     dataField = 'data';
% elseif isfield(MRSIStruct, 'fids')
%     dataField = 'fids';
% else
%     error('MRSIStruct must contain either data or fids field.');
% end
% 
% % Get data dimensions
% dataSize = MRSIStruct.sz;
% dims = MRSIStruct.dims;
% 
% % Validate dimensions
% if length(dataSize) < 3
%     error('Data must have at least 3 dimensions (spectral, x, y)');
% end
% 
% % Find dimension indices
% spectralDim = find([dims.f, dims.t] > 0, 1);
% if isempty(spectralDim)
%     spectralDim = 1; % Default to first dimension
% end
% 
% xDim = dims.x;
% yDim = dims.y;
% 
% if xDim == 0 || yDim == 0
%     error('Spatial dimensions (x, y) must be defined and non-zero');
% end
% 
% % Get data dimensions
% nSpectral = dataSize(spectralDim);
% nX = dataSize(xDim);
% nY = dataSize(yDim);
% 
% if nX == 0 || nY == 0
%     error('Spatial dimensions cannot be zero');
% end
% 
% % Reshape data for processing
% try
%     % Get the data array
%     dataArray = MRSIStruct.(dataField);
% 
%     % Ensure data is properly sized
%     if numel(dataArray) ~= prod(dataSize)
%         warning('Data size does not match expected dimensions. Attempting to reshape...');
%         dataArray = reshape(dataArray, dataSize);
%     end
% 
%     % Reshape to [spectral, spatial] format
%     if spectralDim == 1
%         fidsMatrix = reshape(dataArray, nSpectral, []);
%     else
%         % Move spectral dimension to first position
%         permOrder = 1:length(dataSize);
%         permOrder(1) = spectralDim;
%         permOrder(spectralDim) = 1;
%         fidsPermuted = permute(dataArray, permOrder);
%         fidsMatrix = reshape(fidsPermuted, nSpectral, []);
%     end
% 
% catch ME
%     error('Failed to reshape data: %s', ME.message);
% end
% 
% % Initialize correction maps if computing them
% if computeMaps
%     phaseMap = zeros(nX, nY);
%     freqMap = zeros(nX, nY);
% 
%     % Find reference spectrum (center of k-space or highest SNR)
%     centerX = round(nX/2);
%     centerY = round(nY/2);
%     refIdx = (centerY-1)*nX + centerX;
% 
%     if refIdx > size(fidsMatrix, 2)
%         refIdx = 1; % Fallback to first voxel
%     end
% 
%     refSpectrum = fidsMatrix(:, refIdx);
% 
%     % Set up time vector for fitting
%     if isfield(MRSIStruct, 'adcDwellTime') && MRSIStruct.adcDwellTime > 0
%         dt = MRSIStruct.adcDwellTime;
%     elseif isfield(MRSIStruct, 'spectralDwellTime') && MRSIStruct.spectralDwellTime > 0
%         dt = MRSIStruct.spectralDwellTime;
%     elseif isfield(MRSIStruct, 'dwelltime') && MRSIStruct.dwelltime > 0
%         dt = MRSIStruct.dwelltime;
%     elseif isfield(MRSIStruct, 'spectralwidth') && MRSIStruct.spectralwidth > 0
%         dt = 1 / MRSIStruct.spectralwidth;
%     else
%         dt = 1e-4; % Default 0.1 ms
%         warning('Using default dwell time of 0.1 ms');
%     end
% 
%     t = (0:nSpectral-1)' * dt;
% 
%     % Process each voxel
%     fprintf('Computing B0 correction maps...\n');
% 
%     for voxelIdx = 1:size(fidsMatrix, 2)
%         if mod(voxelIdx, 100) == 0
%             fprintf('Processing voxel %d/%d\n', voxelIdx, size(fidsMatrix, 2));
%         end
% 
%         targetSpectrum = fidsMatrix(:, voxelIdx);
% 
%         % Skip if spectrum is essentially zero
%         if max(abs(targetSpectrum)) < 1e-10
%             continue;
%         end
% 
%         try
%             [phase, freq] = align_LM_robust(targetSpectrum, refSpectrum, t);
% 
%             % Convert linear index to subscripts
%             [yIdx, xIdx] = ind2sub([nY, nX], voxelIdx);
%             phaseMap(yIdx, xIdx) = phase;
%             freqMap(yIdx, xIdx) = freq;
% 
%         catch ME
%             % Skip problematic voxels
%             continue;
%         end
%     end
% 
%     fprintf('B0 correction maps computed.\n');
% end
% 
% % Apply corrections
% fprintf('Applying B0 corrections...\n');
% 
% % Set up time vector for corrections (needed for both compute and apply modes)
% if isfield(MRSIStruct, 'adcDwellTime') && MRSIStruct.adcDwellTime > 0
%     dt = MRSIStruct.adcDwellTime;
% elseif isfield(MRSIStruct, 'spectralDwellTime') && MRSIStruct.spectralDwellTime > 0
%     dt = MRSIStruct.spectralDwellTime;
% elseif isfield(MRSIStruct, 'dwelltime') && MRSIStruct.dwelltime > 0
%     dt = MRSIStruct.dwelltime;
% elseif isfield(MRSIStruct, 'spectralwidth') && MRSIStruct.spectralwidth > 0
%     dt = 1 / MRSIStruct.spectralwidth;
% else
%     dt = 1e-4; % Default 0.1 ms
%     warning('Using default dwell time of 0.1 ms for corrections');
% end
% 
% t = (0:nSpectral-1)' * dt;
% 
% correctedFids = fidsMatrix;
% for voxelIdx = 1:size(fidsMatrix, 2)
%     % Convert linear index to subscripts
%     [yIdx, xIdx] = ind2sub([nY, nX], voxelIdx);
% 
%     if yIdx <= size(phaseMap, 1) && xIdx <= size(phaseMap, 2)
%         phase = phaseMap(yIdx, xIdx);
%         freq = freqMap(yIdx, xIdx);
% 
%         % Apply correction
%         correctionFactor = exp(1j * (phase + 2*pi*freq*t));
%         correctedFids(:, voxelIdx) = fidsMatrix(:, voxelIdx) .* correctionFactor;
%     end
% end
% 
% % Reshape back to original format
% if spectralDim == 1
%     MRSIStruct.(dataField) = reshape(correctedFids, dataSize);
% else
%     % Reshape and permute back
%     correctedReshaped = reshape(correctedFids, [nSpectral, dataSize(2:end)]);
%     permOrder = 1:length(dataSize);
%     permOrder(1) = spectralDim;
%     permOrder(spectralDim) = 1;
%     MRSIStruct.(dataField) = ipermute(correctedReshaped, permOrder);
% end
% 
% % Update flags
% if isfield(MRSIStruct, 'flags')
%     MRSIStruct.flags.freqcorrected = 1;
%     MRSIStruct.flags.phasecorrected = 1;
% else
%     MRSIStruct.flags = struct('freqcorrected', 1, 'phasecorrected', 1);
% end
% 
% fprintf('B0 correction completed.\n');
% 
% end
% 
% %% Helper function for robust alignment
% function [phase, freq] = align_LM_robust(target, reference, t)
%     % Robust Levenberg-Marquardt alignment with better error handling
% 
%     % Ensure both spectra are column vectors
%     target = target(:);
%     reference = reference(:);
%     t = t(:);
% 
%     % Validate inputs
%     if length(target) ~= length(reference) || length(target) ~= length(t)
%         error('Target, reference, and time vectors must have the same length');
%     end
% 
%     % Skip if either spectrum is essentially zero
%     if max(abs(target)) < 1e-10 || max(abs(reference)) < 1e-10
%         phase = 0;
%         freq = 0;
%         return;
%     end
% 
%     % Use real parts for fitting to avoid complex arithmetic issues
%     yData = double(real(target));
% 
%     % Initial parameter guess
%     beta0 = [0, 0]; % [phase, frequency]
% 
%     % Parameter bounds
%     lb = [-pi, -500]; % Lower bounds
%     ub = [pi, 500];   % Upper bounds
% 
%     % Optimization options
%     opts = optimoptions('lsqnonlin', ...
%         'Display', 'off', ...
%         'MaxIterations', 100, ...
%         'FunctionTolerance', 1e-8, ...
%         'OptimalityTolerance', 1e-8);
% 
%     % Model function
%     model = @(b) double(real(reference .* exp(1j*(b(1) + 2*pi*b(2)*t))));
% 
%     % Objective function
%     objective = @(b) yData - model(b);
% 
%     try
%         % Perform optimization
%         beta = lsqnonlin(objective, beta0, lb, ub, opts);
%         phase = beta(1);
%         freq = beta(2);
% 
%     catch ME
%         % Fallback to simple cross-correlation if optimization fails
%         warning('Optimization failed, using fallback method: %s', ME.message);
% 
%         % Simple phase-only correction using cross-correlation
%         crossCorr = sum(target .* conj(reference));
%         phase = angle(crossCorr);
%         freq = 0; % No frequency correction in fallback
%     end
% end