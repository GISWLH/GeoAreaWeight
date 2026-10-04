function w = gaw_band_weights(lat, latBounds, ellipsoid, normalize)
%GAW_BAND_WEIGHTS  Exact area weight of each latitude band |sin(phi_n) - sin(phi_s)|.
%   w = gaw_band_weights(lat)                         bounds inferred from centres (clipped +-90)
%   w = gaw_band_weights(lat, latBounds, ellipsoid, normalize)
%   latBounds: [] | (n+1) edges | (n x 2) [lower upper];  ellipsoid: WGS84 authalic (default false);
%   normalize: sum(w)==1 (default true).  Returns a column vector.
    if nargin < 2, latBounds = []; end
    if nargin < 3 || isempty(ellipsoid), ellipsoid = false; end
    if nargin < 4 || isempty(normalize), normalize = true; end
    lat = double(lat(:)); n = numel(lat);
    if isempty(latBounds)
        b = gaw_infer_bounds(lat, -90, 90); ed = [b(1:n), b(2:n+1)];
    else
        ed = gaw_edges(latBounds, n, 'latBounds');
    end
    if ellipsoid
        s = gaw_authalic_sin_lat(ed);
    else
        s = sin(ed * pi / 180);
    end
    w = abs(s(:, 2) - s(:, 1));
    if normalize, w = w / sum(w); end
end
