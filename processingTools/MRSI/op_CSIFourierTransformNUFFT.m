%function [mrsi, mrsi_w,mrsi_sr,mrsi_sr_w,f,p,f_w,p_w,ccav,ccav_w,ccav_sr,ccav_sr_w,timeCombined,timeCombined_w,sr,sr_w] = rosettePipeline(fileName, fileNameWaterSuppressed, kFile)
%     arguments 
%         fileName (1, :) char {mustBeFile}
%         fileNameWaterSuppressed (1, :) char {mustBeFile}
%         kFile (1, :) char {mustBeFile}
%     end
    %In case k-space file needs to be generated:  
    %obj=RosetteMoreParameters();
    %obj.
    
    fileName='meas_MID01415_FID95676_lb_RosetteSpinEcho_ek_v1a.dat';
    fileName_w='meas_MID01416_FID95677_lb_RosetteSpinEcho_ek_v1a_w.dat';
    kFile='RosetteTraj_40x40_sw1587p3.txt';
    kFile_path = 'C:\Users\divya\Downloads\contents\DataForDivya\invivo_10Sep2024_40x40_rose\invivo_10Sep2024_40x40_rose\RosetteTraj_40x40_sw1587p3.txt';
    %Load data
    rose = io_CSIload_twix1(fileName,kFile_path);
    rose_w = io_CSIload_twix1(fileName_w,kFile_path);
    % rose = roio_CSIload_twix(fileName);
    % rose_w = io_CSIload_twix(fileName_w);

    %Do the complex conj right off the bat to improve the rest of the
    %pipeline:
    rose_cc=rose;
    rose_cc_w=rose_w;
    
    % 
    % %Input correct ADC dwell time (**This needs to be automated):
    % rose_cc.adcDwellTime = 5e-6; 
    % rose_cc_w.adcDwellTime = 5e-6;

    %Combine the ADC blocks
    timeCombined = op_CSICombineTime(rose_cc, 'extras');
    timeCombined_w = op_CSICombineTime(rose_cc_w, 'extras');
   

    % %Re-calculate the adcTime vector, using the correct dwell time (**This
    % %needs to be automated). 
    % timeCombined.adcTime = timeCombined.adcDwellTime*[0:(timeCombined.sz(1)-1)];
    % timeCombined_w.adcTime = timeCombined_w.adcDwellTime*[0:(timeCombined_w.sz(1)-1)];
   
%     %Don't do spectral registration
%     [out_sr,f,p] = op_CSISpecReg(timeCombined,'shortTerm',126);
%     [out_sr_w,f_w,p_w] = op_CSISpecReg(timeCombined_w,'shortTerm',126);

    %Do density compensation
    dComp = op_CSIPSFCorrection_jn(timeCombined,kFile);
    dComp_w = op_CSIPSFCorrection_jn(timeCombined_w,kFile);
    data = dComp.data;
    %% 
    data_resh = permute(data,[2,3,4,1]);
    data_resh=reshape(data_resh,[16,4,63,126,576]);
    data_resh = permute(data_resh,[5,1,2,3,4]);
    data_resh = reshape(data_resh,[576,16,4,7938]);
    dComp.data = data_resh;
    dComp.sz = size(data_resh);
    data = dComp_w.data;
    %% Debugging: Verify reshaping
    figure; imagesc(abs(squeeze(dComp.data(1,:,1,:)))); colorbar;
    title('Debug: Reshaped Data (coil 1, avg 1)');
    data_resh = permute(data,[2,3,1]);
    data_resh=reshape(data_resh,[16,63,126,576]);
    data_resh = permute(data_resh,[4,1,2,3]);
    data_resh = reshape(data_resh,[576,16,7938]);
    dComp_w.data = data_resh;
    dComp_w.sz = size(data_resh);
    dComp = op_CSIFourierTransformNUFFT1(dComp, 'RosetteTraj_40x40_sw1587p3.txt',kFile_path);
    dComp_w = op_CSIFourierTransformNUFFT1(dComp_w, 'RosetteTraj_40x40_sw1587p3.txt',kFile_path);
    % dComp = op_CSIFourierTransform4_nufft(dComp, 'RosetteTraj_40x40_sw1587p3.txt');
    % dComp.dims.ky=0;
    % dComp.dims.x=4;
    % dComp.dims.y=5;
    % dComp_w = op_CSIFourierTransform4_nufft(dComp_w, 'RosetteTraj_40x40_sw1587p3.txt');
    % dComp_w.data = squeeze(dComp_w.data); % Removes singleton dimensions
    % 
    % % Update the size field in the structure
    % dComp_w.sz = size(dComp_w.data);
    % dComp_w.dims.ky=0;
    % dComp_w.dims.x=3;
    % dComp_w.dims.y=4;
    %dComp = op_CSIFourierTransform_nufft(dComp, 'RosetteTraj_40x40_sw1587p3.txt');
    %Do the k->spatial fourier transform
    % ftSpatial = op_CSIFourierTrans(dComp, kFile, 'spatial', true, 'spectral', false);
    % ftSpatial_w = op_CSIFourierTrans(dComp_w, kFile, 'spatial', true, 'spectral', false);

    % ftSpatial = op_CSIFourierTransform1(dComp, kFile, 'spatial', true, 'spectral', false);
    % ftSpatial_w = op_CSIFourierTransform1(dComp_w, kFile, 'spatial', true, 'spectral', false);
    ftSpatial = dComp;
    ftSpatial_w = dComp_w;
    %Combine the RF channels
    [coilCombined_w, phase, weights] = op_CSICombineCoils(ftSpatial_w);
    [coilCombined] = op_CSICombineCoils(ftSpatial, 1, phase, weights);
  
    %Combine the averages
    ccav = op_CSIAverage(coilCombined);
    ccav_w = op_CSIAverage(coilCombined_w);
    figure;
    imagesc(abs(squeeze(ccav_w.data(1,:,:)))); colormap('gray');
    title('Final Combined Brain Image (Time 1)');
    axis equal; axis tight;
    
    %Do the time->spectral fourier transform
    % ftSpec = op_CSIFourierTrans(ccav);
    % ftSpec_w = op_CSIFourierTrans(ccav_w);
    ftSpec = op_CSIFourierTrans(ccav);
    ftSpec_w = op_CSIFourierTrans(ccav_w);
    elapsedTime = toc;
    fprintf('Time taken to run the rosette pipeline: %.2f seconds\n', elapsedTime);
    save('ftSpatialData.mat', 'ftSpatial', 'ftSpatial_w', '-v7.3');
    %load('ftSpatialData.mat', 'ftSpatial', 'ftSpatial_w');

    % Print the contents to the MATLAB command window
    %disp('Contents of ftSpatial:');
    %disp(ftSpatial);

    %disp('Contents of ftSpatial_w:');
    %disp(ftSpatial_w);
    %Remove the residual lipids from the water suppressed data:
    ftSpec_rmlip = op_CSIssp(ftSpec,1.0,1.88);

    %Remove  residual water from the water suppressed data:
    ftSpec_rmw = op_CSIRemoveLipids(ftSpec_rmlip,lipidPPMRange=[4.4 5.0],linewidthRange=[0.5 20]);

    %Do a B0 correction:
    [ftSpec_B0corr_w,phaseMap,freqMap] = op_CSIB0Correction(ftSpec_w);
    [ftSpec_B0corr] = op_CSIB0Correction(ftSpec_rmw,phaseMap,freqMap);
    save('ftSpatial_B0corr.mat', 'ftSpec_B0corr', 'ftSpec_B0corr_w', '-v7.3');

    %Finally, do spatial smoothing
    %% 
    [ftSpec_smooth] = op_CSIApodize(ftSpec_B0corr,'functionType','gaussian','fullWidthHalfMax',30);
    [ftSpec_smooth_w] = op_CSIApodize(ftSpec_B0corr_w,'functionType','gaussian','fullWidthHalfMax',30);
    save('ftSpatial_Smooth.mat', 'ftSpec_smooth', 'ftSpec_smooth_w', '-v7.3');
    op_CSIPlot(ftSpec_smooth);



function MRSIStruct = op_CSIFourierTransformNUFFT1(MRSIStruct, kFile, irtPath)
    arguments
        MRSIStruct (1,1) struct
        kFile (1,:) char
        irtPath (1,:) char = ""
    end

    disp('=== op_CSIFourierTransformNUFFT1_vec_fixed START ===');

    % Add path to IRT if needed
    if ~isempty(irtPath)
        addpath(genpath(irtPath));
    end

    % Read k-space trajectory correctly
    [kTable, kCoords] = readKFile1(kFile);
    figure; plot(kCoords(:,1), kCoords(:,2), '.-'); axis equal;
    title('Debug: Rosette Trajectory'); xlabel('kx'); ylabel('ky');
    Nk = size(kCoords, 1);

    % Verify Nk matches data
    dataIn = MRSIStruct.data;
    szIn = size(dataIn);
    assert(szIn(end) == Nk, ...
        'Mismatch in k-points: data (%d) vs trajectory (%d)', szIn(end), Nk);

    % Set NUFFT parameters
    Nx = 40; Ny = 40;
    Nd = [Nx, Ny];
    Jd = [6, 6];     % typical choice
    Kd = [2*Nx, 2*Ny];
    n_shift = Nd/2;

    % Initialize NUFFT structure
    st = nufft_init(kCoords, Nd, Jd, Kd, n_shift, 'kaiser');

    % Vectorize the data dimensions
    sz = size(MRSIStruct.data);
    Nt = sz(1);
    Ncoils = sz(2);
    Navg = 1;
    if numel(sz) == 4
        Navg = sz(3);
        Nk = sz(4);
        dataIn = reshape(MRSIStruct.data, [Nt*Ncoils*Navg, Nk]).';
    elseif numel(sz) == 3
        Nk = sz(3);
        dataIn = reshape(MRSIStruct.data, [Nt*Ncoils, Nk]);
    else
        error('Data dimensions not supported.');
    end

    % NUFFT Adjoint (vectorized, efficient)
    imgData = nufft_adj(dataIn.', st); % [Nk,L] → [Nx,Ny,L]

    % reshape imgData back to [Nx,Ny,Nt,Ncoils,Navg]
    imgData = reshape(imgData, Nx, Ny, Nt, Ncoils, Navg);
    imgData = permute(imgData, [3,4,5,1,2]); % [Nt,Ncoils,Navg,Nx,Ny]

    if Navg == 1
        imgData = squeeze(imgData); % [Nt,Ncoils,Nx,Ny]
        MRSIStruct.data = imgData;
        MRSIStruct.sz = [Nt, Ncoils, Nx, Ny];
        MRSIStruct.dims = struct('t',1,'coils',2,'x',3,'y',4,'ky',0);
    else
        MRSIStruct.data = imgData;
        MRSIStruct.sz = [Nt, Ncoils, Navg, Nx, Ny];
        MRSIStruct.dims = struct('t',1,'coils',2,'averages',3,'x',4,'y',5,'ky',0);
    end
    MRSIStruct.flags.spatialFT = true;
    [kTable, kArray]  = readKFile(kFile);
    kPtsPerCycle      = getKPtsPerCycle(kTable);
    NPtemporal        = getTemporalPts(kTable, MRSIStruct);
    disp(NPtemporal);
    disp(kPtsPerCycle);
    disp(kTable);
    MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal);
    disp('=== NUFFT Finished ===');
end

% Helper function to correctly read your k-space trajectory file
function [kTable, kArray] = readKFile1(fname)
    opts = detectImportOptions(fname);
    opts.DataLines = [2, Inf]; % Skip header
    kTable = readtable(fname, opts);
    kArray = [kTable{:,2}, kTable{:,3}];
    kMax = max(abs(kArray(:)));
    kArray = (kArray / kMax) * pi;
end

function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
    spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, kPtsPerCycle);
    spectralWidth     = 1 / spectralDwellTime;
    spectralTime      = calculateSpectralTime(spectralDwellTime, NPtemporal);

    MRSIStruct = setSpectralWidth(MRSIStruct, spectralWidth);    
    MRSIStruct = setSpectralDwellTime(MRSIStruct, spectralDwellTime);
    MRSIStruct = setSpectralTime(MRSIStruct, spectralTime);
end


function dt = calculateSpectralDwellTime(MRSIStruct, spatialPoints)
    adcDt = getAdcDwellTime(MRSIStruct);
    dt    = spatialPoints * adcDt;
end

function tAxis = calculateSpectralTime(spectralDwellTime, spatialPoints)
    tAxis = 0:spectralDwellTime:spectralDwellTime*575;
end
