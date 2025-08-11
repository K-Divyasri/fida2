function MRSIStruct = op_CSIPSFCorrection_jn(MRSIStruct, k_space_file, NameValueArgs)
% Density compensation for Rosette CSI with three methods:
%   'nearest'     : k-NN area estimate (fast)
%   'voronoi'     : Voronoi cell areas with model target (Uniform/Gaussian/FlatEdge)
%   'pipe_menon'  : Fixed-point NUFFT (Fessler) as in Pipe & Menon, Mark Chiew
%
% Stores weights in MRSIStruct.densityComp, and applies them to data.

arguments
    MRSIStruct (1,1) struct
    k_space_file (1,:) char {mustBeFile}
    NameValueArgs.method (1,:) char {mustBeMember(NameValueArgs.method,{'nearest','voronoi','pipe_menon'})} = 'nearest'
    NameValueArgs.modelType (1,:) char {mustBeMember(NameValueArgs.modelType,{'Uniform','Gaussian','FlatEdge'})} = 'Uniform'
    NameValueArgs.sigma (1,1) double = sqrt(-(0.1^2)/(2*log(0.01)))
    NameValueArgs.steep (1,1) double = -100
    NameValueArgs.numNeighbors (1,1) double = 10
    NameValueArgs.numIterations (1,1) double = 25
    % plots (kept from your version)
    NameValueArgs.isPlotRosette (1,1) logical = false
    NameValueArgs.isPlotVoronoi (1,1) logical = false
    NameValueArgs.isPlotDensity (1,1) logical = false
    NameValueArgs.isPlotDensityPSF (1,1) logical = false
    NameValueArgs.isPlotDesiredDensity (1,1) logical = false
    NameValueArgs.isPlotDesiredPSF (1,1) logical = false
    NameValueArgs.isPlotWeights (1,1) logical = false
    NameValueArgs.isPlotWeightsPSF (1,1) logical = false
end

% --- guards ---
if isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'spatialft') && MRSIStruct.flags.spatialft
    error('Density Compensation must be performed prior to spatial FT.');
end
if ~(isfield(MRSIStruct,'seq') && contains(MRSIStruct.seq,["rst","ros"],'IgnoreCase',true))
    error('Apply only to Rosette sequences.');
end

% --- k-space and dims ---
[kx, ky] = readK(k_space_file);
Nt     = MRSIStruct.sz(MRSIStruct.dims.t);
Ncoils = MRSIStruct.sz(MRSIStruct.dims.coils);
Nkpts  = MRSIStruct.sz(MRSIStruct.dims.kpts);
Nkshot = MRSIStruct.sz(MRSIStruct.dims.kshot);
hasAverages = isfield(MRSIStruct.dims,'averages') && MRSIStruct.dims.averages>0;
if hasAverages, Nave = MRSIStruct.sz(MRSIStruct.dims.averages); else, Nave=1; end
if numel(kx) ~= Nkpts*Nkshot
    error('K-file length (%d) != Nkpts*Nkshot (%d*%d).', numel(kx), Nkpts, Nkshot);
end
if NameValueArgs.isPlotRosette, figure; plot(kx,ky,'.-'); axis equal; title('Rosette trajectory'); end

% --- compute weights per method ---
switch lower(NameValueArgs.method)
    case 'nearest'
        [rho, kRep] = density_nearest(kx, ky, NameValueArgs.numNeighbors);
        rho = rho .* kRep;
        rho = rho / max(min(rho(rho>0)), eps);
        w = 1 ./ max(rho, eps);

    case 'voronoi'
        [rho, w] = density_voronoi(kx, ky, MRSIStruct, ...
            NameValueArgs.modelType, NameValueArgs.sigma, NameValueArgs.steep);

    case 'pipe_menon'
        [w, rho] = density_pipe_menon(kx, ky, MRSIStruct, NameValueArgs.numIterations);
end

% optional plots
if NameValueArgs.isPlotDensity
    figure; scatter3(kx,ky,rho,12,rho,'filled'); title('Sampling density'); xlabel Kx; ylabel Ky; zlabel \rho;
end
if NameValueArgs.isPlotWeights
    figure; scatter3(kx,ky,w,12,w,'filled'); title('DCF weights'); xlabel Kx; ylabel Ky; zlabel w;
end

% --- apply weights to data ---
W2 = reshape(w,[Nkpts Nkshot]);                 % [kpts x kshot]
if hasAverages && ndims(MRSIStruct.data)==5
    fullW = ones(Nt,Ncoils,Nave,Nkpts,Nkshot,'like',MRSIStruct.data);
    for t=1:Nt, for c=1:Ncoils, for a=1:Nave, fullW(t,c,a,:,:) = W2; end, end, end
elseif ~hasAverages && ndims(MRSIStruct.data)==4
    fullW = ones(Nt,Ncoils,Nkpts,Nkshot,'like',MRSIStruct.data);
    for t=1:Nt, for c=1:Ncoils, fullW(t,c,:,:) = W2; end, end
else
    error('Unexpected data dims: %s', mat2str(size(MRSIStruct.data)));
end
MRSIStruct.data = MRSIStruct.data .* fullW;

% --- stash for PSF ---
MRSIStruct.densityComp = struct( ...
    'method', lower(NameValueArgs.method), ...
    'numNeighbors', NameValueArgs.numNeighbors, ...
    'numIterations', NameValueArgs.numIterations, ...
    'modelType', NameValueArgs.modelType, ...
    'sigma', NameValueArgs.sigma, ...
    'steep', NameValueArgs.steep, ...
    'rho', rho(:), 'weights', w(:), 'kx', kx(:), 'ky', ky(:));

end

% ================= helpers =================
function [kx,ky] = readK(f)
    T = readtable(f);
    if any(strcmpi(T.Properties.VariableNames,'Kx')), kx=T.Kx; ky=T.Ky;
    elseif any(strcmpi(T.Properties.VariableNames,'kx')), kx=T.kx; ky=T.ky;
    else, error('Kx/Ky columns not found in %s', f);
    end
    kx = double(kx(:)); ky = double(ky(:));
end

function [rho, kRep] = density_nearest(kx, ky, K)
    N = numel(kx); rho = zeros(N,1); kRep=rho;
    rmax = max(hypot(kx,ky));
    for i=1:N
        d = hypot(kx-kx(i), ky-ky(i));
        kRep(i) = sum(d==0);
        sd = sort(d);
        kdist = sd(min(K+kRep(i), numel(sd)));
        db = rmax - hypot(kx(i), ky(i));
        cf = 1 + max(0, 1 - db/max(kdist,eps));
        rho(i) = cf * K / (pi*max(kdist,eps)^2);
    end
end

function [rho, w] = density_voronoi(kx, ky, S, modelType, sigma, steep)
    % Voronoi cell areas -> rho = 1/area. Then map to desired model (Uniform/Gaussian/FlatEdge)
    [xc, yc, z0, ~] = cleanK(kx, ky);   % keep one origin sample
    [v,c] = voronoin([xc yc]);
    FOV = getFOV(S); Nx = getNx(S);
    Kmax = 1 / (2*(FOV/Nx));
    areas = zeros(size(c,1),1);
    for i=1:size(c,1)
        idx = c{i};
        if any(idx==1) || numel(idx)<3, areas(i)=NaN; continue; end
        vx = v(idx,1); vy = v(idx,2);
        r = hypot(vx,vy); o = r > Kmax;     % clip to circle
        vx(o) = Kmax*vx(o)./r(o); vy(o)=Kmax*vy(o)./r(o);
        areas(i) = polyarea(vx,vy);
    end
    % fill NaNs fairly
    areas = fillmissing(areas,'constant',median(areas(~isnan(areas))));
    rho_c = 1./max(areas,eps);

    % desired rho model on cleaned set
    R = hypot(xc,yc);
    switch lower(modelType)
        case 'uniform'
            rho_des = double(R < Kmax);
        case 'gaussian'
            rho_des = max(rho_c)*exp(-(xc.^2/(2*sigma^2) + yc.^2/(2*sigma^2)));
        case 'flatedge'
            rho_des = -1./(1+exp(steep*(R - 0.64*Kmax))) + 1;
    end
    w_c = rho_des ./ max(rho_c,eps);   % weights on cleaned points

    % map back to full ordering (replicate origin count)
    [w, rho] = mapBack(kx,ky,xc,yc,w_c,rho_c,z0);
    rho = rho / max(min(rho(rho>0)), eps);
    w   = 1 ./ max(rho,eps);           % ensure consistent definition
end

function [w, rho] = density_pipe_menon(kx, ky, S, iters)
    % Mark Chiew / Fessler fixed-point scheme
    assert(exist('nufft_init','file')==2,'Fessler IRT not found.');
    % normalise to [-pi,pi]
    kmax = max(hypot(kx,ky));
    om = [ (kx./max(kmax,eps))*pi, (ky./max(kmax,eps))*pi ];
    Nx = getNx(S); Nd=[Nx Nx]; Jd=[6 6]; Kd=[2*Nx 2*Nx]; nshift=Nd/2;
    st = nufft_init(om, Nd, Jd, Kd, nshift);    % table-less ok

    w = ones(size(kx));
    for ii=1:iters
        % w = w ./ real( A * (A^H * w) )
        ktmp = st.p' * w;                % A^H w  -> image
        ktmp = st.p  * ktmp;             % A (A^H w) -> back to k
        den  = real(ktmp);  den(abs(den)<1e-12) = 1e-12;
        w    = w ./ den;
    end
    % normalise like Chiew’s code
    w = w * st.sn(ceil(end/2),ceil(end/2))^(-2)/prod(st.Kd);
    w = sqrt(abs(w));
    rho = 1./max(w,eps);                 % pseudo-density for reporting
    rho = rho / max(min(rho(rho>0)), eps);
end

function [xc,yc,zerosCount,keepOrigin] = cleanK(kx,ky)
    idx0 = (kx==0 & ky==0);
    zerosCount = sum(idx0);
    keep = true(size(kx));
    % keep exactly one origin
    if zerosCount>0
        first = find(idx0,1,'first');
        keep(idx0) = false; keep(first)=true;
        keepOrigin = first;
    else
        keepOrigin = [];
    end
    xc = kx(keep); yc = ky(keep);
end

function [w, rho] = mapBack(kx,ky,xc,yc,wc,rc,zeroIdx)
    N = numel(kx); w = zeros(N,1); rho = zeros(N,1);
    for i=1:N
        d = hypot(xc-kx(i), yc-ky(i));
        [m, j] = min(d);
        if m<1e-12, w(i)=wc(j); rho(i)=rc(j); else, w(i)=wc(j); rho(i)=rc(j); end
    end
    % replicate origin weight if there were repeats
    % (already handled because all identical points map to same j)
end

function FOV = getFOV(S)
    if isfield(S,'fov') && isfield(S.fov,'x'), FOV = S.fov.x; else, FOV = 240; end
end
function Nx = getNx(S)
    if isfield(S,'voxelSize') && isfield(S.voxelSize,'x') && S.voxelSize.x>0
        Nx = round(getFOV(S) / S.voxelSize.x);
    else
        Nx = 64;
    end
end
