function MRS = op_CSItoMRS(MRSIStruct, x, y)
    % Get data and dimension indices
    data = MRSIStruct.data;
    dim_f = MRSIStruct.dims.f;
    dim_x = MRSIStruct.dims.x;
    dim_y = MRSIStruct.dims.y;

    % Build indexing cell array
    idx = repmat({':'}, 1, numel(MRSIStruct.sz));
    idx{dim_x} = x;
    idx{dim_y} = y;

    % Extract spectrum at (x,y)
    specs = squeeze(data(idx{:}));

    % Create output MRS structure
    MRS = struct();
    MRS.specs = specs;
    MRS.fids = fft(ifftshift(specs, dim_f), [], dim_f);
    MRS.sz = size(MRS.specs);
    MRS.ppm = MRSIStruct.ppm;
    MRS.t = MRSIStruct.spectralTime;
    MRS.dwelltime = MRSIStruct.spectralDwellTime;
    MRS.spectralwidth = MRSIStruct.spectralWidth;
    MRS.txfrq = MRSIStruct.txfrq;
    MRS.Bo = MRSIStruct.Bo;
    MRS.te = MRSIStruct.te;
    MRS.tr = MRSIStruct.tr;

    % === Add dims field ===
    MRS.dims = struct();
    MRS.dims.t = 1;  % always true after squeeze

    % === Copy flags including isFourSteps ===
    MRS.flags = struct();
    flagFields = fieldnames(MRSIStruct.flags);
    for i = 1:length(flagFields)
        MRS.flags.(flagFields{i}) = MRSIStruct.flags.(flagFields{i});
    end
end

% function MRS = op_CSItoMRS(MRSIStruct, x, y)
%     % Get data and dimension indices
%     data = MRSIStruct.data;
%     dim_f = MRSIStruct.dims.f;
%     dim_x = MRSIStruct.dims.x;
%     dim_y = MRSIStruct.dims.y;
% 
%     % Build indexing cell array
%     idx = repmat({':'}, 1, numel(MRSIStruct.sz));
%     idx{dim_x} = x;
%     idx{dim_y} = y;
% 
%     % Extract spectrum at (x,y)
%     specs = squeeze(data(idx{:}));
% 
%     % Create output MRS structure
%     MRS = struct();
%     MRS.specs = specs;
%     MRS.fids = fft(ifftshift(specs, dim_f), [], dim_f);
%     MRS.sz = size(MRS.specs);
%     MRS.ppm = MRSIStruct.ppm;
%     MRS.t = MRSIStruct.spectralTime;
%     MRS.dwelltime = MRSIStruct.spectralDwellTime;
%     MRS.spectralwidth = MRSIStruct.spectralWidth;
%     MRS.txfrq = MRSIStruct.txfrq;
%     MRS.Bo = MRSIStruct.Bo;
%     MRS.te = MRSIStruct.te;
%     MRS.tr = MRSIStruct.tr;
%     MRS.flags = struct();
%     MRS = setFlags(MRS, 'averaged', false);
%     MRS = setFlags(MRS, 'addedrcvrs', true); % Assume already combined
%     MRS = setFlags(MRS, 'spectralft', true);
% end
% 
% % % CSItoTwix.m
% % % converts a voxel from CSI to an MRS structure for FID-A. The CSI voxel
% % % selected is based on MRSIndex variable.
% % % 
% % % USAGE:
% % % out = op_CSItoMRS(in, xCoordinate, yCoordinate)
% % %
% % % INPUT: 
% % % in = FID-A CSI object
% % % xCoordinate = X coordinate wanted to be converted
% % % yCoordinate = Y coorindate wanted to be converted
% % %
% % % OUTPUT:
% % % out = single voxel FID-A object
% % 
% % function MRS = op_CSItoMRS(MRSIStruct, xCoordinate, yCoordinate, index)
% %     arguments
% %         MRSIStruct (1, 1) struct
% %         xCoordinate (1, 1) double
% %         yCoordinate (1, 1) double
% %         index.averageIndex (1, 1) double = 0
% %         index.coilIndex (1, 1) double = 0
% %         index.subSpecIndex (1, 1) double = 0
% %         index.extraIndex (1, 1) double = 0
% %         index.linearIndex (1, 1) double = 0
% %     end
% %     checkArguments(MRSIStruct, index);
% %     MRS = MRSIStruct;
% % 
% %     MRSIndex = buildIndex(MRSIStruct, xCoordinate, yCoordinate, index);
% % 
% %     data = getData(MRSIStruct);
% %     timeDimension = getDimension(MRSIStruct, 't');
% %     if(getFlags(MRSIStruct, 'spectralft') == true)
% %         MRS.specs = data(MRSIndex{:});
% %         %MRS.fids = ifft(fftshift(MRS.specs, timeDimension), [],  timeDimension);
% %         MRS.fids = fft(fftshift(MRS.specs, timeDimension), [],  timeDimension);
% %     else
% %         MRS.fids = data(MRSIndex{:});
% %         MRS.specs = fftshift(ifft(MRS.fids, [], timeDimension), timeDimension);
% %     end
% %     MRS = rmfield(MRS, 'data');
% %     MRS.specs = squeeze(MRS.specs);
% % 
% %     %adjust structure to match FID-A object
% %     MRS.sz = size(MRS.specs);
% %     MRS.t = MRSIStruct.spectralTime;
% %     MRS.dwelltime = MRSIStruct.spectralDwellTime;
% %     MRS.spectralwidth = MRSIStruct.spectralWidth;
% %     if(isfield(MRSIStruct, 'ppm'))
% %         MRS.ppm = MRSIStruct.ppm;
% %     else
% %         tSize = getSizeFromDimensions(MRSIStruct, {'t'});
% %         step = MRSIStruct.spectralWidth/tSize;
% %         freqBounds = MRSIStruct.spectralWidth/2 - step/2;
% %         frequency = -freqBounds:step:freqBounds;
% % 
% %         MRS.ppm = frequency/(MRSIStruct.Bo*MRSIStruct.gamma);
% %     end
% % 
% % 
% %     MRS = removeDimension(MRS, 'x');
% %     MRS = removeDimension(MRS, 'y');
% %     %remove dimensions if they are being indexed
% %     if(index.averageIndex >= 1); MRS = removeDimension(MRS, 'averages'); end
% %     if(index.coilIndex >= 1); MRS = removeDimension(MRS, 'coils'); end
% %     if(index.subSpecIndex >= 1); MRS = removeDimension(MRS, 'subspec'); end
% %     if(index.extraIndex >= 1); MRS = removeDimension(MRS, 'extras'); end
% % end
% % 
% % function checkArguments(in, index)
% %     if(getFlags(in, 'spatialft') == false)
% %         error('please fourier transform along the spatial and spectral dimension before using this function');
% %     end
% %     if(index.linearIndex > 0 && (index.averageIndex || index.coilIndex || index.subSpecIndex || index.extraIndex))
% %     end
% % 
% % end
% % 
% % function MRSIndex = buildIndex(MRSIStruct, xCoordinate, yCoordinate, index)
% %     %initalize index variable
% %     MRSIndex = cell(length(MRSIStruct.sz), 1);
% %     %all index dimensions are default to colon. (ie. all values)
% %     for i = 1:length(MRSIndex)
% %         MRSIndex{i} = ':';
% %     end
% % 
% %     %set x and y indexes
% %     MRSIndex{MRSIStruct.dims.x} = xCoordinate;
% %     MRSIndex{MRSIStruct.dims.y} = yCoordinate;
% % 
% %     averageDimension = getDimension(MRSIStruct, 'averages');
% %     %if averages exist set average index
% %     if(averageDimension ~= 0 && index.averageIndex > 0)
% %         MRSIndex{averageDimension} = index.averageIndex;
% %     end
% % 
% %     coilDimension = getDimension(MRSIStruct, 'coils');
% %     %if averages exist set average index
% %     if(coilDimension ~= 0 && index.coilIndex > 0)
% %         MRSIndex{coilDimension} = index.coilIndex;
% %     end
% % 
% %     subSpecDimension = getDimension(MRSIStruct, 'subspecs');
% %     %if averages exist set average index
% %     if(subSpecDimension ~= 0 && index.subSpecIndex > 0)
% %         MRSIndex{subSpecDimension} = index.subSpecIndex;
% %     end
% % 
% %     extraDimension = getDimension(MRSIStruct, 'extras');
% %     %if averages exist set average index
% %     if(extraDimension ~= 0 && index.extraIndex > 0)
% %         MRSIndex{extraDimension} = index.extraIndex;
% %     end
% % end
