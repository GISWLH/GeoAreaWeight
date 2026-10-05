function w = gaw_area_weights(lat, latBounds, method, normalize)
%GAW_AREA_WEIGHTS  1-D latitude weights (column vector, one per latitude).
%   w = gaw_area_weights(lat)                               exact band weights, sum(w) == 1
%   w = gaw_area_weights(lat, latBounds, method, normalize)
%   method: 'band' (default, exact sphere) | 'cos' | 'ellipsoid' (WGS84) | 'none'
%
%   See also GAW_BAND_WEIGHTS, GAW_COS_LAT_WEIGHTS, GAW_WEIGHTS_2D.
    if nargin < 2, latBounds = []; end
    if nargin < 3 || isempty(method), method = 'band'; end
    if nargin < 4 || isempty(normalize), normalize = true; end
    switch gaw_canonical_method(method)
        case 'band'
            w = gaw_band_weights(lat, latBounds, false, normalize); return
        case 'ellipsoid'
            w = gaw_band_weights(lat, latBounds, true, normalize); return
        case 'cos'
            w = gaw_cos_lat_weights(lat(:));
        case 'none'
            gaw_check_lat(lat);
            w = ones(numel(lat), 1);
    end
    if normalize, w = w / sum(w); end
end
