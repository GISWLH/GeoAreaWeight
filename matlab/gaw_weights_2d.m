function W = gaw_weights_2d(lat, lon, latBounds, lonBounds, method, normalize)
%GAW_WEIGHTS_2D  2-D weight field (nlat x nlon).
%   W = gaw_weights_2d(lat, lon)                                   exact band weights, sum(W(:)) == 1
%   W = gaw_weights_2d(lat, lon, latBounds, lonBounds, method, normalize)
%   A 2-D (curvilinear) lat gives cos(lat) weights; 'band'/'ellipsoid' then warn.
%   Pass real cell areas (e.g. areacella) to GAW_AREA_MEAN via 'Weights' instead.
%
%   See also GAW_AREA_WEIGHTS, GAW_CELL_AREA, GAW_AREA_MEAN.
    if nargin < 2, lon = []; end
    if nargin < 3, latBounds = []; end
    if nargin < 4, lonBounds = []; end
    if nargin < 5 || isempty(method), method = 'band'; end
    if nargin < 6 || isempty(normalize), normalize = true; end
    m = gaw_canonical_method(method);
    if ~isvector(lat)
        if strcmp(m, 'none')
            W = ones(size(lat));
        else
            if ~strcmp(m, 'cos')
                warning('gaw:curvilinear', ['method ''%s'' needs 1-D latitudes; using cos(lat) for the ' ...
                        '2-D grid. Pass cell areas via ''Weights'' for exact results.'], method);
            end
            W = gaw_cos_lat_weights(lat);
        end
    else
        if isempty(lon), error('gaw:lon', 'lon is required to build a 2-D weight field'); end
        w1 = gaw_area_weights(lat, latBounds, m, false);
        if strcmp(m, 'none')
            W = ones(numel(lat), numel(lon));
        else
            W = w1(:) * gaw_lon_widths(lon, lonBounds, true)';
        end
    end
    if normalize, W = W / sum(W(:)); end
end
