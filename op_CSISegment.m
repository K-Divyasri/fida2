function MRSIStruct = op_CSISegment(MRSIStruct, blob, stepsize)
% op_CSISegment - Segments lipid and brain regions from CSI magnitude image.
%
% Inputs:
%   MRSIStruct : CSI data structure (e.g., ccav)
%   blob       : Max size of small components to remove (default = 50)
%   stepsize   : Threshold increment step for adaptive thresholding (default = 0.01)
%
% Output:
%   MRSIStruct.mask : Structure containing 'lipmasks', 'brainmasks', and thresholds

arguments
    MRSIStruct (1,1) struct
    blob (1,1) double = 50
    stepsize (1,1) double = 0.01
end

%% Step 1: Extract the 1st timepoint (magnitude image)
ccav = abs(squeeze(MRSIStruct.data(1,:,:)));

% Apply a Gaussian filter for smoothing
I = imgaussfilt(ccav);

% Initialize mask storage
counter = 0;

%% Step 2: Loop over thresholds
for sens = 0:stepsize:1
    % Adaptive threshold
    thresh = adaptthresh(I, sens);
    BW = imbinarize(I, thresh);

    % Morphological operations
    SE = strel('disk', 1);
    BW = imdilate(BW, SE);
    BW = imerode(BW, SE);

    % Remove small blobs
    clean_mask = bwareaopen(BW, blob);

    % Analyze connected components
    CC = bwconncomp(clean_mask);
    stats = regionprops(CC, 'Centroid','ConvexArea','Solidity', ...
        'FilledImage','FilledArea');

    % Reject masks that don't meet structural criteria
    if length(stats) > 1 || isempty(stats)
        % disp(['Rejected threshold ', num2str(sens), ': multiple or no regions']);
        continue;
    end

    s = stats(1);
    if s.Solidity > 0.5
        % disp(['Rejected threshold ', num2str(sens), ': solidity too high']);
        continue;
    elseif size(s.FilledImage,1) < size(s.FilledImage,2)
        % disp(['Rejected threshold ', num2str(sens), ': not ring-like']);
        continue;
    elseif s.FilledArea < s.ConvexArea/2
        % disp(['Rejected threshold ', num2str(sens), ': not complete ring']);
        continue;
    end

    % Passed all conditions → save masks
    counter = counter + 1;
    lipmasks(:,:,counter)   = clean_mask;
    brainmasks(:,:,counter) = imfill(clean_mask, 'holes');
    threshes(counter,1)     = sens;
end

%% Step 3: Handle case of no valid masks
if counter == 0
    warning('No valid masks found. Returning empty mask structure.');
    MRSIStruct.mask = struct( ...
        'lipmasks', [], ...
        'brainmasks', [], ...
        'thresh', [], ...
        'numthresh', 0 ...
    );
    return;
end

%% Step 4: Visualize results
figure;
subplot(2,2,1), imshow(ccav, []), title('Original Image');
subplot(2,2,2), imshow(lipmasks(:,:,1), []), title('Clean Segmented Image (lipid)');
subplot(2,2,3), imshow(brainmasks(:,:,1)-lipmasks(:,:,1), []), title('Brain Mask (outer minus lipid)');
subplot(2,2,4), imshow(ccav - lipmasks(:,:,1), []), title('Subtracted Image');

%% Step 5: Package results
mask.lipmasks    = lipmasks;
mask.brainmasks  = logical(brainmasks(:,:,1) - lipmasks(:,:,1));
mask.thresh      = threshes;
mask.numthresh   = counter;

MRSIStruct.mask = mask;
end

% function MRSIStruct = op_CSISegment(MRSIStruct, blob, stepsize)
% % op_CSISegment.m
% % Kaito Hara-Lee, SunnyBrook Hospital 2024.
% %
% % Description: Creates a 2D mask of the skull using adaptive thresholding.
% %              All viable masks generated and their respective thresholds
% %              are saved to a field called "Masks".
% %
% % USAGE:
% % out = op_CSISegment(in, blob, stepsize)
% %
% % input:  in        = MRSI Structure (ccav)
% %         blob      = Max pixel group size to remove (default = 50)
% %         stepsize  = Increment step for threshold scan (default = 0.01)
% %
% % output: out       = Structure with added .mask field
% 
% arguments
%     MRSIStruct (1,1) struct
%     blob (1,1) double = 50
%     stepsize (1,1) double = 0.01
% end
% 
% %% === Step 1: Extract central slice ===
% ccav = abs(squeeze(MRSIStruct.data(1,:,:)));
% I = imgaussfilt(ccav);  % smooth image
% 
% counter = 0;
% 
% %% === Step 2: Thresholding loop ===
% for sens = 0:stepsize:1
%     thresh = adaptthresh(I, sens);
%     BW = imbinarize(I, thresh);
% 
%     SE = strel('disk', 1);    % smooth binary edges
%     BW = imdilate(BW, SE);
%     BW = imerode(BW, SE);
% 
%     clean_mask = bwareaopen(BW, blob);  % remove small blobs
% 
%     CC = bwconncomp(clean_mask);
%     stats = regionprops(CC, 'Centroid', 'ConvexArea', ...
%         'Solidity', 'FilledImage', 'FilledArea');
% 
%     % Debug why rejections happen (optional)
%     if length(stats) > 1 || isempty(stats)
%         continue;
%     elseif stats.Solidity > 0.5
%         continue;
%     elseif size(stats.FilledImage,1) < size(stats.FilledImage,2)
%         continue;
%     elseif stats.FilledArea < stats.ConvexArea/2
%         continue;
%     end
% 
%     % Passed all checks
%     counter = counter + 1;
%     lipmasks(:, :, counter) = clean_mask;
%     brainmasks(:, :, counter) = imfill(clean_mask, 'holes');
%     threshes(counter, :) = sens;
% end
% 
% %% === Step 3: Finalize mask and display ===
% 
% if counter == 0
%     warning('No valid masks found. Returning empty mask structure.');
%     MRSIStruct.mask = struct('lipmasks', [], ...
%                              'brainmasks', [], ...
%                              'thresh', [], ...
%                              'numthresh', 0);
%     return;
% end
% 
% % Visualize results
% figure;
% subplot(2,2,1), imshow(ccav, []), title('Original Image');
% subplot(2,2,2), imshow(lipmasks(:,:,1), []), title('Clean Segmented Image (lipid)');
% subplot(2,2,3), imshow(brainmasks(:,:,1) - lipmasks(:,:,1), []), title('Clean Segmented Image (brain)');
% subplot(2,2,4), imshow(ccav - lipmasks(:,:,1), []), title('Subtracted Image');
% 
% % Output mask structure
% mask.lipmasks = lipmasks;
% mask.brainmasks = logical(brainmasks(:,:,1) - lipmasks(:,:,1));
% mask.thresh = threshes;
% mask.numthresh = counter;
% 
% MRSIStruct.mask = mask;
% 
% end
