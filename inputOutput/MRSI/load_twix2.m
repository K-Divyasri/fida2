
function [timeCombined_sh, timeCombined_sh_w] = load_twix2(fileName, fileName_w, kFile)
    % Load both files completely in one read each
    fprintf('Loading primary file: %s\n', fileName);
    [rose, x, y] = io_CSIload_twix3(fileName,kFile);
    
    fprintf('Loading water-suppressed file: %s\n', fileName_w);
    [rose_w, ~, ~] =io_CSIload_twix3(fileName_w,kFile); % x,y should be same for both files
    
    % Combine time dimensions
    fprintf('Combining time dimensions...\n');
    timeCombined   = op_CSICombineTime1(rose, 'extras');
    timeCombined_w = op_CSICombineTime1(rose_w, 'extras');
    
    % Set adcTime
    %timeCombined.adcTime   = timeCombined.adcDwellTime * (0:(timeCombined.sz(1)-1));
    %timeCombined_w.adcTime = timeCombined_w.adcDwellTime * (0:(timeCombined_w.sz(1)-1));
    
    % Shift
    fprintf('Applying spatial shifts...\n');
    if isempty(kFile)
        % Cartesian shift using default or dummy inputs
        fprintf('Using Cartesian mode (no k-space file)\n');
        timeCombined_sh   = op_CSIShift(timeCombined, x, y, "");
        timeCombined_sh_w = op_CSIShift(timeCombined_w, x, y, "");
    else
        % Non-Cartesian shift using k-space trajectory file
        fprintf('Using non-Cartesian mode with k-file: %s\n', kFile);
        timeCombined_sh   = op_CSIShift(timeCombined, x, y, kFile);
        timeCombined_sh_w = op_CSIShift(timeCombined_w, x, y, kFile);
    end
    
    fprintf('File loading and initial processing completed.\n');
end


% function [mrsistruct, x, y] = load_twix_complete(filename,kFile)
%     % Single complete read of the twix file
%     % This replaces the separate header and data loading steps
% 
%     fprintf('Reading twix file: %s\n', filename);
%     [mrsistruct, x, y] = io_CSIload_twix3(filename,kFile);
% 
%     % Set the ADC dwell time (this was previously done in load_twix_header)
% 
%     fprintf('Completed reading file. Data size: %s\n', mat2str(mrsistruct.sz));
% end