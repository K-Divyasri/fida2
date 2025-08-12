%% === INPUTS ===
fileName = 'C:\Users\divya\Downloads\Lexar\DataForDivya\invivo_10Sep2024_40x40_rose\invivo_10Sep2024_40x40_rose\meas_MID01415_FID95676_lb_RosetteSpinEcho_ek_v1a.dat';
fileName_w = 'C:\Users\divya\Downloads\Lexar\DataForDivya\invivo_10Sep2024_40x40_rose\invivo_10Sep2024_40x40_rose\meas_MID01416_FID95677_lb_RosetteSpinEcho_ek_v1a_w.dat';
kFile = 'C:\Users\divya\Downloads\Lexar\DataForDivya\invivo_10Sep2024_40x40_rose\invivo_10Sep2024_40x40_rose\RosetteTraj_40x40_sw1587p3.txt';

fprintf('=== TWO FILE NON-CARTESIAN MODE ===\n');

%% === LOAD TWIX FILES ===
[timeCombined_sh, timeCombined_sh_w] = load_twix2(fileName, fileName_w,kFile);

nComp = 8;  % choose, e.g., 8 virtual channels out of your 32
[timeCombined,  W] = op_CSICoilCompression_SVD(timeCombined_sh,  nComp);
[timeCombined_w, ~] = op_CSICoilCompression_SVD(timeCombined_sh_w, nComp, false);
timeCombined_w.coilCompression.W = W;    % overwrite voxel weights
 % overwrite voxel weights

%% === RESHAPE ===
timeCombined_rs = reshape_twix_data(timeCombined_sh,kFile);
timeCombined_rs_w = reshape_twix_data(timeCombined_sh_w,kFile);

%% === K-SPACE COtRRECTION AND NUFFT ===
fprintf('Applying k-space correction and NUFFT...\n');
dComp = op_CSIPSFCorrection_jn(timeCombined_rs, kFile);
dComp_w = op_CSIPSFCorrection_jn(timeCombined_rs_w, kFile);


ftSpatial = op_NUFFTSpatial(dComp, kFile);
ftSpatial_w = op_NUFFTSpatial(dComp_w, kFile);

%% === COIL COMBINATION ===
fprintf('Combining coils...\n');
[coilCombined_w, phase, weights] = op_CSICombineCoils1(ftSpatial_w);
coilCombined = op_CSICombineCoils1(ftSpatial, 1, phase, weights);

%% === COMBINE AVERAGES ===
fprintf('Combining averages...\n');
ccav = op_CSIAverage(coilCombined);
ccav_w = op_CSIAverage(coilCombined_w);

%% === SPECTRAL FT ===
fprintf('Performing spectral Fourier transform...\n');
ftSpec = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);

%% === LIPID REMOVAL ===
fprintf('Removing lipids...\n');
ftSpec_rmlip = op_CSIssp(ftSpec, 1.0, 1.88);
ftSpec_rmw = op_CSIRemoveLipids(ftSpec_rmlip, lipidPPMRange=[4.4 5.0], linewidthRange=[0.5 20]);


%Do a B0 correction:
[ftSpec_B0corr_w,phaseMap,freqMap] = op_CSIB0Correction(ftSpec_w);
[ftSpec_B0corr] = op_CSIB0Correction(ftSpec_rmw,phaseMap,freqMap);
   

%% === APODIZATION ===
fprintf('Applying Gaussian apodization...\n');
ftSpec_smooth = op_CSIApodize(ftSpec_B0corr, 'functionType', 'gaussian', 'fullWidthHalfMax', 25);
ftSpec_smooth_w = op_CSIApodize(ftSpec_B0corr_w, 'functionType', 'gaussian', 'fullWidthHalfMax', 25);

%% === PLOT ===
op_CSIPlot(ftSpec_smooth);

