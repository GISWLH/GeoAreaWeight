function out = gaw_area_integral(x, lat, lon, varargin)
%GAW_AREA_INTEGRAL  Area integral sum(x .* cellArea) over the grid (per-area flux -> total).
%   t = gaw_area_integral(x, lat, lon)             x is lon x lat (x time ...); result in x-units * m2
%   t = gaw_area_integral(x, lat, lon, 'Units', 'km2', 'Mask', M, ...)
%
%   Name-value options (case-insensitive):
%     'Ellipsoid' false (default) | true: WGS84 cell areas
%     'Radius'    sphere radius [m], default 6371008.8
%     'Units'     'm2' (default) | 'km2'
%     'Weights'   nlat x nlon cell areas to use instead (any unit)
%     'Mask', 'LatBounds', 'LonBounds', 'LatDim' (2), 'LonDim' (1)
%     'SkipNaN'   true (default): NaN cells contribute 0; NaN only if no valid cell
%
%   See also GAW_CELL_AREA, GAW_AREA_MEAN.
    o = gaw_parse_options(struct('Ellipsoid', false, 'Radius', [], 'Units', 'm2', 'Weights', [], ...
                                 'Mask', [], 'LatBounds', [], 'LonBounds', [], 'LatDim', 2, ...
                                 'LonDim', 1, 'SkipNaN', true), varargin, 'gaw_area_integral');
    d = size(x); d(end+1:max(o.LatDim, o.LonDim)) = 1;
    if isempty(o.Weights)
        if isempty(lat) || isempty(lon), error('gaw:args', 'lat and lon (or ''Weights'') are required'); end
        A = gaw_cell_area(lat, lon, o.LatBounds, o.LonBounds, o.Ellipsoid, o.Radius, o.Units);
    else
        A = gaw_grid_field(o.Weights, d(o.LatDim), d(o.LonDim), 'Weights', 'weights');
    end
    out = gaw_core(x, A, o.LatDim, o.LonDim, o.Mask, o.SkipNaN, 'sum');
end
