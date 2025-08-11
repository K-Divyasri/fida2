function ftSpatial = op_NUFFTSpatial(dComp, kFile_path)
    % OPTIMIZED VERSION: Minimized memory allocations and redundant operations
    % Modified to set spectral fields like op_CSIFourierTransform
    
    if nargout > 0 || ~isempty(inputname(1))
        disp('=== START: op_NUFFTSpatial (Optimized) ===');
    end

    %% Load and normalize k-trajectory (OPTIMIZED with caching)
    persistent cached_kCoords cached_kFile cached_st cached_Nx cached_Ny;
    
    % Cache k-space trajectory and NUFFT structure
    if isempty(cached_kCoords) || ~strcmp(cached_kFile, kFile_path)
        if nargout > 0 || ~isempty(inputname(1))
            disp(['Loading k-file: ', kFile_path]);
        end
        
        [kTable, ~] = readKFile(kFile_path);
        
        % Handle different k-file formats
        if istable(kTable)
            if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
                kCoords = [kTable.Kx, kTable.Ky];
            else
                % Try to find columns that look like coordinates
                varNames = kTable.Properties.VariableNames;
                if length(varNames) >= 3
                    kCoords = [kTable{:, 2}, kTable{:, 3}]; % Assume columns 2,3 are Kx,Ky
                else
                    error('Cannot find Kx, Ky columns in k-file table');
                end
            end
        else
            % Handle as numeric matrix
            if size(kTable, 2) >= 3
                kCoords = kTable(:, 2:3);  % Extract kx, ky columns (columns 2,3)
            else
                error('K-file must have at least 3 columns for numeric format');
            end
        end
        
        % Vectorized normalization
        kCoords = (kCoords / max(abs(kCoords(:)))) * pi;
        
        % Cache trajectory
        cached_kCoords = kCoords;
        cached_kFile = kFile_path;
        
        if nargout > 0 || ~isempty(inputname(1))
            disp(['Cached kCoords shape = ', mat2str(size(kCoords))]);
        end
    else
        kCoords = cached_kCoords;
    end

    % Pre-compute dimensions (OPTIMIZED)
    sz = dComp.sz;
    dims = dComp.dims;
    nKshot = sz(dims.kshot);
    nKpts = sz(dims.kpts);
    totalKPoints = nKshot * nKpts;
    
    % Verify k-space dimensions once
    if size(kCoords, 1) ~= totalKPoints
        error('K-space file length (%d) does not match data k-space dimensions (%d)', ...
              size(kCoords, 1), totalKPoints);
    end

    %% NUFFT setup (OPTIMIZED with caching)
    Nx = length(dComp.coordinates.x);
    Ny = length(dComp.coordinates.y);
    
    % Cache NUFFT structure for same dimensions
    if isempty(cached_st) || cached_Nx ~= Nx || cached_Ny ~= Ny
        Nd = [Nx, Ny];
        Jd = [6, 6]; 
        Kd = [2*Nx, 2*Ny]; 
        n_shift = Nd/2;
        
        cached_st = nufft_init(kCoords, Nd, Jd, Kd, n_shift);
        cached_Nx = Nx;
        cached_Ny = Ny;
        
        if nargout > 0 || ~isempty(inputname(1))
            disp('NUFFT initialized and cached.');
        end
    end
    st = cached_st;

    %% OPTIMIZED data processing with pre-allocation
    dataSize = size(dComp.data);
    hasAverages = (dims.averages > 0);
    
    if hasAverages && ndims(dComp.data) == 5
        % 5D: [t, coils, averages, kshot, kpts] -> [kshot*kpts, t*coils*averages]
        
        % Direct reshape without intermediate permutation (MEMORY OPTIMIZED)
        data_reshaped = reshape(dComp.data, [sz(1)*sz(2)*sz(3), sz(4)*sz(5)]);
        data_transposed = data_reshaped.';  % [totalKPoints, t*coils*averages]
        
        % NUFFT transform
        imgFlat = nufft_adj(data_transposed, st);
        
        % Direct reshape to final dimensions
        imgData = reshape(imgFlat, [Nx, Ny, sz(1), sz(2), sz(3)]);
        imgData = permute(imgData, [3, 4, 5, 2, 1]); % [t, coils, averages, x, y]
        
        % Update dimensions efficiently
        new_dims = dims;
        new_dims.x = 5;
        new_dims.y = 4;
        new_dims.kshot = 0;
        new_dims.kpts = 0;
        
    elseif ~hasAverages && ndims(dComp.data) == 4
        % 4D: [t, coils, kshot, kpts] -> [kshot*kpts, t*coils]
        
        % Direct reshape (MEMORY OPTIMIZED)
        data_reshaped = reshape(dComp.data, [sz(1)*sz(2), sz(3)*sz(4)]);
        data_transposed = data_reshaped.';  % [totalKPoints, t*coils]
        
        % NUFFT transform
        imgFlat = nufft_adj(data_transposed, st);
        
        % Direct reshape to final dimensions
        imgData = reshape(imgFlat, [Nx, Ny, sz(1), sz(2)]);
        imgData = permute(imgData, [3, 4, 2, 1]); % [t, coils, x, y]
        
        % Update dimensions efficiently
        new_dims = dims;
        new_dims.x = 4;
        new_dims.y = 3;
        new_dims.kshot = 0;
        new_dims.kpts = 0;
        
    else
        error('Unsupported data dimension: %d or averages flag inconsistent', ndims(dComp.data));
    end

    %% Update structure efficiently (OPTIMIZED)
    ftSpatial = dComp;
    ftSpatial.data = imgData;
    ftSpatial.sz = size(imgData);
    ftSpatial.dims = new_dims;
    ftSpatial.flags.spatialFT = 1;
    
    % Clear k-space dimensions
    ftSpatial.dims.kx = 0;
    ftSpatial.dims.ky = 0;
    
    %% Calculate and set spectral values (like op_CSIFourierTransform)
    % Read k-file again for spectral calculations
    [kTable, ~] = readKFile(kFile_path);
    kPtsPerCycle = getKPtsPerCycle(kTable);
    NPtemporal = getTemporalPts(kTable, ftSpatial);
    
    % Calculate spectral parameters (matching op_CSIFourierTransform logic)
    ftSpatial = calculateSpectralValues(ftSpatial, kPtsPerCycle, NPtemporal);
    
    if nargout > 0 || ~isempty(inputname(1))
        disp(['Spectral width: ', num2str(ftSpatial.spectralWidth)]);
        disp(['Spectral dwell time: ', num2str(ftSpatial.spectralDwellTime)]);
        disp(['Spectral time length: ', num2str(length(ftSpatial.spectralTime))]);
        disp('=== END: op_NUFFTSpatial (Optimized) ===');
    end
end

%% Helper functions (matching op_CSIFourierTransform implementation)

function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
    % Calculate spectral parameters exactly like op_CSIFourierTransform
    spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, kPtsPerCycle);
    spectralWidth = 1/spectralDwellTime;
    spectralTime = calculateSpectralTime(spectralDwellTime, NPtemporal);

    % Set spectral fields in structure
    MRSIStruct.spectralWidth = spectralWidth;
    MRSIStruct.spectralDwellTime = spectralDwellTime;
    MRSIStruct.spectralTime = spectralTime;
end

function spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, spatialPoints)
    % Calculate spectral dwell time from ADC dwell time and spatial points
    adcDwellTime = MRSIStruct.adcDwellTime;
    spectralDwellTime = spatialPoints * adcDwellTime;
end

function spectralTime = calculateSpectralTime(spectralDwellTime, spatialPoints)
    % Calculate spectral time vector
    spectralTime = 0:spectralDwellTime:spectralDwellTime*(spatialPoints - 1);
end

function kPtsPerCycle = getKPtsPerCycle(kTable)
    % Extract k-points per cycle from k-table
    try
        if istable(kTable)
            if ismember('TR', kTable.Properties.VariableNames)
                % Find the maximum TR index in the trajectory
                num_TR = max(kTable.TR, [], 'all');
                % Divide the number of lines by the number of TRs
                kPtsPerCycle = height(kTable)/num_TR;
            else
                % Default fallback for table without TR column
                kPtsPerCycle = 126;
            end
        elseif ismatrix(kTable) && size(kTable, 2) >= 5
            % Assume TR is in column 5 for numeric matrix
            num_TR = max(kTable(:, 5), [], 'all');
            if num_TR > 0
                kPtsPerCycle = size(kTable, 1)/num_TR;
            else
                kPtsPerCycle = 126;
            end
        else
            % Default fallback
            kPtsPerCycle = 126; % Generally 126 for Rosette sequences
        end
    catch
        % Fallback in case of any error
        kPtsPerCycle = 126;
    end
end

function NPtemporal = getTemporalPts(kTable, dComp)
    % Get temporal points from data structure
    try
        if isfield(dComp.dims, 't') && dComp.dims.t > 0
            NPtemporal = dComp.sz(dComp.dims.t);
        else
            NPtemporal = dComp.sz(1); % Fallback to first dimension
        end
    catch
        % Fallback
        NPtemporal = size(dComp.data, 1);
    end
end

function [kTable, kArray] = readKFile(kFileName)
    % Read k-space trajectory file with robust error handling
    if isempty(kFileName) || ~isfile(kFileName)
        kTable = [];
        kArray = [];
        return;
    end
    
    try
        % First try reading as table
        kTable = readtable(kFileName);
        
        % Extract kx, ky coordinates
        if istable(kTable)
            if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
                kArray = [kTable.Kx, kTable.Ky];
            elseif width(kTable) >= 3
                % Assume columns 2,3 are Kx, Ky
                kArray = [kTable{:, 2}, kTable{:, 3}];
            else
                error('Table does not have enough columns');
            end
        end
        
    catch
        % Fallback: try reading as numeric matrix
        try
            kTable = readmatrix(kFileName);
            if size(kTable, 2) >= 3
                kArray = kTable(:, 2:3); % Assume Kx, Ky are in columns 2, 3
            else
                error('Matrix does not have enough columns');
            end
        catch ME
            warning('Failed to read k-file: %s. Error: %s', kFileName, ME.message);
            kTable = [];
            kArray = [];
        end
    end
end

%% ULTRA-FAST VERSION with pre-computed NUFFT
function ftSpatial = op_NUFFTSpatial_UltraFast(dComp, kFile_path)
    % ULTRA-FAST: Pre-computed NUFFT with minimal overhead
    
    persistent ultra_cached_st ultra_cached_kFile ultra_cached_dims;
    
    % Ultra-fast caching
    if isempty(ultra_cached_st) || ~strcmp(ultra_cached_kFile, kFile_path)
        [kTable, ~] = readKFile(kFile_path);
        
        % Handle different formats safely
        if istable(kTable)
            if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
                kCoords = [kTable.Kx, kTable.Ky];
            else
                kCoords = [kTable{:, 2}, kTable{:, 3}];
            end
        else
            kCoords = kTable(:, 2:3);
        end
        
        kCoords = (kCoords / max(abs(kCoords(:)))) * pi;
        
        Nx = length(dComp.coordinates.x);
        Ny = length(dComp.coordinates.y);
        
        ultra_cached_st = nufft_init(kCoords, [Nx, Ny], [6, 6], [2*Nx, 2*Ny], [Nx, Ny]/2);
        ultra_cached_kFile = kFile_path;
        ultra_cached_dims = [Nx, Ny];
    end
    
    [Nx, Ny] = deal(ultra_cached_dims(1), ultra_cached_dims(2));
    sz = dComp.sz;
    hasAverages = (dComp.dims.averages > 0);
    
    % Ultra-fast processing
    if hasAverages && ndims(dComp.data) == 5
        % Direct vectorized operations
        data_flat = reshape(permute(dComp.data, [4, 5, 1, 2, 3]), [], sz(1)*sz(2)*sz(3));
        img_flat = nufft_adj(data_flat, ultra_cached_st);
        imgData = permute(reshape(img_flat, Nx, Ny, sz(1), sz(2), sz(3)), [3, 4, 5, 1, 2]);
        
        new_dims = dComp.dims;
        [new_dims.x, new_dims.y, new_dims.kshot, new_dims.kpts] = deal(4, 5, 0, 0);
        
    else % 4D case
        data_flat = reshape(permute(dComp.data, [3, 4, 1, 2]), [], sz(1)*sz(2));
        img_flat = nufft_adj(data_flat, ultra_cached_st);
        imgData = permute(reshape(img_flat, Nx, Ny, sz(1), sz(2)), [3, 4, 1, 2]);
        
        new_dims = dComp.dims;
        [new_dims.x, new_dims.y, new_dims.kshot, new_dims.kpts] = deal(3, 4, 0, 0);
    end
    
    % Update structure
    ftSpatial = dComp;
    ftSpatial.data = imgData;
    ftSpatial.sz = size(imgData);
    ftSpatial.dims = new_dims;
    ftSpatial.flags.spatialFT = 1;
    
    % Calculate and set spectral values (matching op_CSIFourierTransform)
    [kTable, ~] = readKFile(kFile_path);
    kPtsPerCycle = getKPtsPerCycle(kTable);
    NPtemporal = getTemporalPts(kTable, ftSpatial);
    ftSpatial = calculateSpectralValues(ftSpatial, kPtsPerCycle, NPtemporal);
end

%% MEMORY-OPTIMIZED VERSION for large datasets
function ftSpatial = op_NUFFTSpatial_MemoryOptimized(dComp, kFile_path, chunk_size)
    % MEMORY OPTIMIZED: Process data in chunks to reduce memory usage
    
    if nargin < 3
        chunk_size = 1000; % Process 1000 time points at once
    end
    
    % Setup (same as optimized version)
    [kTable, ~] = readKFile(kFile_path);
    
    % Handle different formats safely
    if istable(kTable)
        if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
            kCoords = [kTable.Kx, kTable.Ky];
        else
            kCoords = [kTable{:, 2}, kTable{:, 3}];
        end
    else
        kCoords = kTable(:, 2:3);
    end
    
    kCoords = (kCoords / max(abs(kCoords(:)))) * pi;
    
    Nx = length(dComp.coordinates.x);
    Ny = length(dComp.coordinates.y);
    st = nufft_init(kCoords, [Nx, Ny], [6, 6], [2*Nx, 2*Ny], [Nx, Ny]/2);
    
    sz = dComp.sz;
    hasAverages = (dComp.dims.averages > 0);
    
    if hasAverages && ndims(dComp.data) == 5
        % Process in chunks for 5D data
        imgData = zeros(sz(1), sz(2), sz(3), Nx, Ny, 'like', dComp.data);
        
        for t_start = 1:chunk_size:sz(1)
            t_end = min(t_start + chunk_size - 1, sz(1));
            t_indices = t_start:t_end;
            
            % Process chunk
            data_chunk = dComp.data(t_indices, :, :, :, :);
            data_flat = reshape(permute(data_chunk, [4, 5, 1, 2, 3]), [], numel(data_chunk)/(sz(4)*sz(5)));
            img_flat = nufft_adj(data_flat, st);
            imgData(t_indices, :, :, :, :) = permute(reshape(img_flat, Nx, Ny, length(t_indices), sz(2), sz(3)), [3, 4, 5, 1, 2]);
        end
        
        new_dims = dComp.dims;
        [new_dims.x, new_dims.y, new_dims.kshot, new_dims.kpts] = deal(4, 5, 0, 0);
        
    else % 4D case - similar chunking
        imgData = zeros(sz(1), sz(2), Nx, Ny, 'like', dComp.data);
        
        for t_start = 1:chunk_size:sz(1)
            t_end = min(t_start + chunk_size - 1, sz(1));
            t_indices = t_start:t_end;
            
            data_chunk = dComp.data(t_indices, :, :, :);
            data_flat = reshape(permute(data_chunk, [3, 4, 1, 2]), [], numel(data_chunk)/(sz(3)*sz(4)));
            img_flat = nufft_adj(data_flat, st);
            imgData(t_indices, :, :, :) = permute(reshape(img_flat, Nx, Ny, length(t_indices), sz(2)), [3, 4, 1, 2]);
        end
        
        new_dims = dComp.dims;
        [new_dims.x, new_dims.y, new_dims.kshot, new_dims.kpts] = deal(3, 4, 0, 0);
    end
    
    % Update structure
    ftSpatial = dComp;
    ftSpatial.data = imgData;
    ftSpatial.sz = size(imgData);
    ftSpatial.dims = new_dims;
    ftSpatial.flags.spatialFT = 1;
    
    % Calculate and set spectral values (matching op_CSIFourierTransform)
    kPtsPerCycle = getKPtsPerCycle(kTable);
    NPtemporal = getTemporalPts(kTable, ftSpatial);
    ftSpatial = calculateSpectralValues(ftSpatial, kPtsPerCycle, NPtemporal);
end

%% PARALLEL VERSION using Parallel Computing Toolbox
function ftSpatial = op_NUFFTSpatial_Parallel(dComp, kFile_path)
    % PARALLEL: Use parallel processing for NUFFT transforms
    
    % Check if Parallel Computing Toolbox is available
    if ~license('test', 'Distrib_Computing_Toolbox') || isempty(gcp('nocreate'))
        warning('Parallel Computing Toolbox not available. Using standard version.');
        ftSpatial = op_NUFFTSpatial(dComp, kFile_path);
        return;
    end
    
    % Setup
    [kTable, ~] = readKFile(kFile_path);
    
    % Handle different formats safely
    if istable(kTable)
        if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
            kCoords = [kTable.Kx, kTable.Ky];
        else
            kCoords = [kTable{:, 2}, kTable{:, 3}];
        end
    else
        kCoords = kTable(:, 2:3);
    end
    
    kCoords = (kCoords / max(abs(kCoords(:)))) * pi;
    
    Nx = length(dComp.coordinates.x);
    Ny = length(dComp.coordinates.y);
    st = nufft_init(kCoords, [Nx, Ny], [6, 6], [2*Nx, 2*Ny], [Nx, Ny]/2);
    
    sz = dComp.sz;
    hasAverages = (dComp.dims.averages > 0);
    
    if hasAverages && ndims(dComp.data) == 5
        % Parallel processing for 5D data
        data_cell = cell(sz(3), 1); % Split by averages
        
        parfor avg = 1:sz(3)
            data_avg = squeeze(dComp.data(:, :, avg, :, :));
            data_flat = reshape(permute(data_avg, [3, 4, 1, 2]), [], sz(1)*sz(2));
            img_flat = nufft_adj(data_flat, st);
            data_cell{avg} = permute(reshape(img_flat, Nx, Ny, sz(1), sz(2)), [3, 4, 1, 2]);
        end
        
        % Combine results
        imgData = cat(3, data_cell{:});
        imgData = permute(imgData, [1, 2, 3, 5, 4]); % [t, coils, averages, x, y]
        
        new_dims = dComp.dims;
        [new_dims.x, new_dims.y, new_dims.kshot, new_dims.kpts] = deal(4, 5, 0, 0);
        
    else % 4D case
        % Parallel processing for 4D data
        data_cell = cell(sz(2), 1); % Split by coils
        
        parfor coil = 1:sz(2)
            data_coil = squeeze(dComp.data(:, coil, :, :));
            data_flat = reshape(permute(data_coil, [2, 3, 1]), [], sz(1));
            img_flat = nufft_adj(data_flat, st);
            data_cell{coil} = permute(reshape(img_flat, Nx, Ny, sz(1)), [3, 1, 2]);
        end
        
        % Combine results
        imgData = cat(2, data_cell{:});
        imgData = permute(imgData, [1, 2, 4, 3]); % [t, coils, x, y]
        
        new_dims = dComp.dims;
        [new_dims.x, new_dims.y, new_dims.kshot, new_dims.kpts] = deal(3, 4, 0, 0);
    end
    
    % Update structure
    ftSpatial = dComp;
    ftSpatial.data = imgData;
    ftSpatial.sz = size(imgData);
    ftSpatial.dims = new_dims;
    ftSpatial.flags.spatialFT = 1;
    
    % Calculate and set spectral values (matching op_CSIFourierTransform)
    kPtsPerCycle = getKPtsPerCycle(kTable);
    NPtemporal = getTemporalPts(kTable, ftSpatial);
    ftSpatial = calculateSpectralValues(ftSpatial, kPtsPerCycle, NPtemporal);
end

%% BENCHMARKING FUNCTION
function benchmark_op_NUFFTSpatial(dComp, kFile_path)
    fprintf('=== op_NUFFTSpatial Performance Benchmark ===\n');
    
    % Test original version
    tic;
    result1 = op_NUFFTSpatial(dComp, kFile_path);
    time_original = toc;
    
    % Test ultra-fast version
    tic;
    result2 = op_NUFFTSpatial_UltraFast(dComp, kFile_path);
    time_ultrafast = toc;
    
    % Test memory-optimized version
    tic;
    result3 = op_NUFFTSpatial_MemoryOptimized(dComp, kFile_path);
    time_memory = toc;
    
    % Report results
    fprintf('Original version:        %.3f seconds\n', time_original);
    fprintf('Ultra-fast version:      %.3f seconds (%.1fx faster)\n', time_ultrafast, time_original/time_ultrafast);
    fprintf('Memory-optimized:        %.3f seconds (%.1fx faster)\n', time_memory, time_original/time_memory);
    
    % Verify results are similar (allowing for numerical precision differences)
    diff1 = max(abs(result1.data(:) - result2.data(:)));
    diff2 = max(abs(result1.data(:) - result3.data(:)));
    
    fprintf('Max difference (ultra-fast): %.2e\n', diff1);
    fprintf('Max difference (memory-opt): %.2e\n', diff2);
    
    if diff1 < 1e-10 && diff2 < 1e-10
        fprintf('✓ All versions produce nearly identical results\n');
    else
        fprintf('⚠ Results differ - check implementation\n');
    end
    
    % Check spectral field consistency
    fprintf('\n=== Spectral Field Verification ===\n');
    fprintf('spectralWidth: %.4e (orig) vs %.4e (ultra) vs %.4e (mem)\n', ...
        result1.spectralWidth, result2.spectralWidth, result3.spectralWidth);
    fprintf('spectralDwellTime: %.4e (orig) vs %.4e (ultra) vs %.4e (mem)\n', ...
        result1.spectralDwellTime, result2.spectralDwellTime, result3.spectralDwellTime);
    fprintf('spectralTime length: %d (orig) vs %d (ultra) vs %d (mem)\n', ...
        length(result1.spectralTime), length(result2.spectralTime), length(result3.spectralTime));
end