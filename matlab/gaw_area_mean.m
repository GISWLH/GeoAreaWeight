function out = gaw_area_mean(x, lat, lon, varargin)
%GAW_AREA_MEAN  Area-weighted spatial mean; other dimensions (e.g. time) are kept.
%   m = gaw_area_mean(x, lat, lon)        x is lon x lat (x time ...) as returned by ncread
%   m = gaw_area_mean(x, lat, lon, 'Method', 'cos', 'Mask', M, 'LatDim', 2, 'LonDim', 1, ...)
%   m = gaw_area_mean(x, [], [], 'Weights', areacella)
%
%   Name-value options (case-insensitive):
%     'Method'     'band' (default, exact sphere) | 'cos' | 'ellipsoid' (exact WGS84) |
%                  'none' (plain unweighted mean, for comparison)
%     'Weights'    nlat x nlon cell areas (e.g. CMIP areacella); overrides Method. NaN -> 0.
%     'Mask'       nlat x nlon logical (true = include); 0 / NaN exclude
%     'LatBounds', 'LonBounds'   (n+1) edges or n x 2 bounds; inferred if omitted
%     'LatDim' (default 2), 'LonDim' (default 1)
%     'SkipNaN'    true (default): skip NaN and renormalise the weights over valid
%                  cells; false: any NaN inside the (masked) domain gives NaN
%   2-D fields (Weights, Mask) are nlat x nlon; an nlon x nlat field is
%   transposed automatically when nlat ~= nlon.
%
%   Example
%     lat = (-89.5:89.5)'; lon = 0.5:359.5;
%     x = repmat(abs(lat'), numel(lon), 1);          % lon x lat
%     gaw_area_mean(x, lat, lon)                     % 32.71 (naive mean(x(:)) = 45)
%
%   See also GAW_REGION_MEAN, GAW_AREA_INTEGRAL, GAW_AREA_FRACTION, GAW_WEIGHTS_2D.
    o = gaw_parse_options(struct('Method', 'band', 'Weights', [], 'Mask', [], 'LatBounds', [], ...
                                 'LonBounds', [], 'LatDim', 2, 'LonDim', 1, 'SkipNaN', true), ...
                          varargin, 'gaw_area_mean');
    d = size(x); d(end+1:max(o.LatDim, o.LonDim)) = 1;
    if isempty(o.Weights)
        if isempty(lat) || isempty(lon), error('gaw:args', 'lat and lon (or ''Weights'') are required'); end
        W = gaw_weights_2d(lat, lon, o.LatBounds, o.LonBounds, o.Method, false);
    else
        W = gaw_grid_field(o.Weights, d(o.LatDim), d(o.LonDim), 'Weights', 'weights');
    end
    out = gaw_core(x, W, o.LatDim, o.LonDim, o.Mask, o.SkipNaN, 'mean');
end
