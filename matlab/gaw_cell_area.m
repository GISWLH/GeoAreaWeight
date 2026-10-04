function A = gaw_cell_area(lat, lon, latBounds, lonBounds, ellipsoid, radius, units)
%GAW_CELL_AREA  Area of lat-lon cells [m2 or km2]; (nlat x nlon) if lon given,
%   otherwise a column vector per radian of longitude.
%   ellipsoid=true -> WGS84 (authalic radius); radius default 6371008.8 m. units 'm2'|'km2'.
    if nargin < 2, lon = []; end
    if nargin < 3, latBounds = []; end
    if nargin < 4, lonBounds = []; end
    if nargin < 5 || isempty(ellipsoid), ellipsoid = false; end
    if nargin < 6 || isempty(radius), radius = 6371008.8; end
    if nargin < 7 || isempty(units), units = 'm2'; end
    if ellipsoid, R = gaw_authalic_radius(); else, R = radius; end
    f = gaw_band_weights(lat, latBounds, ellipsoid, false);
    if isempty(lon)
        A = R^2 * f;
    else
        lon = double(lon(:)); n = numel(lon);
        if isempty(lonBounds)
            b = gaw_infer_bounds(lon); ed = [b(1:n), b(2:n+1)];
        else
            ed = gaw_edges(lonBounds, n, 'lonBounds');
        end
        dlon = abs(ed(:, 2) - ed(:, 1)) * pi / 180;
        A = R^2 * (f * dlon');
    end
    if strcmp(units, 'km2'), A = A / 1e6; elseif ~strcmp(units, 'm2'), error('gaw:units', 'units must be m2 or km2'); end
end
