function A = gaw_cell_area(lat, lon, latBounds, lonBounds, ellipsoid, radius, units)
%GAW_CELL_AREA  Area of lat-lon cells [m2 or km2], exact for graticule cells.
%   A = gaw_cell_area(lat, lon)                      nlat x nlon, m2, sphere R = 6371008.8 m
%   A = gaw_cell_area(lat, lon, latBounds, lonBounds, ellipsoid, radius, units)
%   A = gaw_cell_area(lat)                           column vector: area per radian of longitude
%   ellipsoid=true -> WGS84 (authalic radius 6371007.181 m); units 'm2' (default) | 'km2'.
%   Longitudes may cross the 0/360 or +-180 meridian.
%
%   See also GAW_BAND_WEIGHTS, GAW_AREA_INTEGRAL.
    if nargin < 2, lon = []; end
    if nargin < 3, latBounds = []; end
    if nargin < 4, lonBounds = []; end
    if nargin < 5 || isempty(ellipsoid), ellipsoid = false; end
    if nargin < 6 || isempty(radius), radius = 6371008.8; end
    if nargin < 7 || isempty(units), units = 'm2'; end
    if ~any(strcmp(units, {'m2', 'km2'})), error('gaw:units', 'units must be ''m2'' or ''km2'''); end
    if ellipsoid, R = gaw_authalic_radius(); else, R = radius; end
    f = gaw_band_weights(lat, latBounds, ellipsoid, false);
    if isempty(lon)
        A = R^2 * f;
    else
        dlon = gaw_lon_widths(lon, lonBounds, false) * pi / 180;
        A = R^2 * (f * dlon');
    end
    if strcmp(units, 'km2'), A = A / 1e6; end
end
