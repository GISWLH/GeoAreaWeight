function w = gaw_area_weights(lat, latBounds, method, normalize)
%GAW_AREA_WEIGHTS  1-D latitude weights.
%   method: 'cos' | 'band' (default, exact spherical band) | 'ellipsoid' (WGS84) | 'none'
    if nargin < 2, latBounds = []; end
    if nargin < 3 || isempty(method), method = 'band'; end
    if nargin < 4 || isempty(normalize), normalize = true; end
    m = lower(method);
    switch m
        case 'cos'
            w = gaw_cos_lat_weights(lat(:));
        case {'band', 'sphere', 'sphere_band'}
            w = gaw_band_weights(lat, latBounds, false, normalize); return
        case {'ellipsoid', 'wgs84'}
            w = gaw_band_weights(lat, latBounds, true, normalize); return
        case {'none', 'equal', 'unweighted'}
            w = ones(numel(lat), 1);
        otherwise
            error('gaw:method', 'unknown method: %s', method);
    end
    if normalize, w = w / sum(w); end
end
