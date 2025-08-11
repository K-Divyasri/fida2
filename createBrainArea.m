function [brain_mask, lipid_mask] = createBrainArea(MRSIStructSpectral, origin_x, origin_y, varargin)
% CREATEBRAINAREA  Build logical brain & lipid masks from a CSI data set.
%
%   [brain_mask, lipid_mask] = createBrainArea(specStruct, x0, y0)
%   [...] = createBrainArea(specStruct, x0, y0, 'Name',Value, ...)
%
% INPUTS
%   specStruct   : FID‑A CSI structure *already* Fourier‑transformed in space & freq.
%   x0 , y0      : Floating‑point coordinates (in voxel units) of the brain centre.
%
% NAME–VALUE OPTIONS
%   'Plot'       : true  → show QC figure (default = false)
%   'SpatialStruct': additional CSI struct in the spatial domain (for nicer QC plots)
%   'ZThresh'    : Z‑score threshold for lipid detection (default = 3)
%
% OUTPUTS
%   brain_mask   : logical matrix (Ny×Nx)  — voxels inside the convex hull.
%   lipid_mask   : logical matrix (Ny×Nx)  — voxels along the outer lipid ring.
%
% NOTES / IMPROVEMENTS OVER ORIGINAL VERSION
%   • Accepts Name–Value syntax instead of an options struct.
%   • Works whether or not an extras/averages dimension exists.
%   • Uses built‑in *bwconvhull* to obtain the convex hull (faster & safer).
%   • Ensures lipid ring is exactly a 1‑pixel perimeter via *bwperim*.
%   • All intermediate indices are bounds‑checked.
% -------------------------------------------------------------------------

%% ------------------ 0. Parse inputs -----------------------------------------
    p = inputParser;
    p.addRequired ( 'MRSIStructSpectral', @(s)isstruct(s)&&isfield(s,'data') );
    p.addRequired ( 'origin_x', @isnumeric );
    p.addRequired ( 'origin_y', @isnumeric );
    p.addParameter( 'Plot',          false, @islogical );
    p.addParameter( 'SpatialStruct', struct(), @(s)isstruct(s) );
    p.addParameter( 'ZThresh',       3,     @(x)isnumeric(x)&&x>0 );
    p.parse(MRSIStructSpectral, origin_x, origin_y, varargin{:});
    args = p.Results;

%% ------------------ 1. 2‑D lipid integral map -------------------------------
    lipid3D = op_CSIintegrate(args.MRSIStructSpectral, 0.9, 1.8, 'mag');  % y×x×extras
    lipid2D = squeeze(lipid3D(:,:,1));                                    % force 2‑D
    [Ny, Nx] = size(lipid2D);

%% ------------------ 2. Set up ray geometry ----------------------------------
    maxR     = hypot(Nx/2, Ny/2);               % radius to an image corner
    nSlices  = ceil( 2*pi / atan(1/(max(Nx,Ny)/2)) );
    angles   = linspace(0, 2*pi, nSlices+1); angles(end) = [];

    rayMask  = false(Ny, Nx);

%% ------------------ 3. Ray‑wise lipid voxel detection ------------------------
    for ang = angles
        % tip point in 2‑D image coordinates
        x_tip =  origin_x + maxR*cos(-ang);
        y_tip =  origin_y + maxR*sin(-ang);

        [cx, cy, c] = improfile(abs(lipid2D), [origin_x, x_tip], [origin_y, y_tip]);
        if isempty(c), continue, end
        valid = isfinite(c(:));  c=c(valid);  cx=cx(valid);  cy=cy(valid);
        if isempty(c), continue, end

        zc = zscore(double(c));
        candidate = zc > args.ZThresh;
        if ~any(candidate)          % fallback: pick global max
            [~, idx] = max(c);
            candidate = false(size(c)); candidate(idx)=true;
        end
        % write into mask (bounds‑checked)
        for k = find(candidate).'
            xi = min(max(round(cx(k)),1), Nx);
            yi = min(max(round(cy(k)),1), Ny);
            rayMask(yi,xi) = true;
        end
    end

%% ------------------ 4. Build lipid & brain masks ----------------------------
    if ~any(rayMask(:))
        brain_mask = false(Ny,Nx);
        lipid_mask = false(Ny,Nx);
        return
    end

    lipid_mask = bwperim(bwconvhull(rayMask));   % 1‑pixel lipid perimeter
    brain_mask = imfill(~lipid_mask,'holes');    % everything inside hull
    brain_mask(lipid_mask) = false;             % ensure disjoint sets

%% ------------------ 5. Optional QC plot -------------------------------------
    if args.Plot
        figure('Color','k'); hold on
        if ~isempty(fieldnames(args.SpatialStruct))
            imagesc(abs(squeeze(args.SpatialStruct.data(1,:,:))))
        else
            imagesc(abs(lipid2D))
        end
        colormap(gray); axis image ij off
        scatter(find(brain_mask), brain_mask(brain_mask), 15, 'g', 'filled');
        scatter(find(lipid_mask), lipid_mask(lipid_mask), 15, 'r', 'filled');
        title('createBrainArea – QC','Color','w')
        legend({'Brain','Lipid'},'TextColor','w');
    end
end
