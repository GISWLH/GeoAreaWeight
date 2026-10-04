function out = gaw_area_mean(x, lat, lon, varargin)
%GAW_AREA_MEAN  Area-weighted spatial mean; other dimensions (e.g. time) are kept.
%   m = gaw_area_mean(x, lat, lon)                 x is lon x lat (x time ...) as from ncread
%   m = gaw_area_mean(x, lat, lon, 'Method','cos', 'Mask',M, 'LatDim',2, 'LonDim',1, ...)
%   Name-value options:
%     'Method'    'band' (default) | 'cos' | 'ellipsoid' | 'none'
%     'Weights'   nlat x nlon cell-area matrix (e.g. areacella); overrides Method
%     'Mask'      nlat x nlon logical (true = include): regional / land-only mean
%     'LatBounds','LonBounds'   (n+1) edges or n x 2
%     'LatDim' (default 2), 'LonDim' (default 1)
%     'SkipNaN'   true (default): weights renormalised over finite cells
    o = struct('Method', 'band', 'Weights', [], 'Mask', [], 'LatBounds', [], 'LonBounds', [], ...
               'LatDim', 2, 'LonDim', 1, 'SkipNaN', true);
    for k = 1:2:numel(varargin)
        if ~isfield(o, varargin{k}), error('gaw:option', 'unknown option %s', varargin{k}); end
        o.(varargin{k}) = varargin{k+1};
    end
    if isempty(o.Weights)
        W = gaw_weights_2d(lat, lon, o.LatBounds, o.LonBounds, o.Method, false);
    else
        W = double(o.Weights);
    end
    out = gaw_core(x, W, o.LatDim, o.LonDim, o.Mask, o.SkipNaN);
end
