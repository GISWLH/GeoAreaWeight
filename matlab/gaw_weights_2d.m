function W = gaw_weights_2d(lat, lon, latBounds, lonBounds, method, normalize)
%GAW_WEIGHTS_2D  (nlat x nlon) weight field.  2-D lat (curvilinear): cos(lat) only.
    if nargin < 2, lon = []; end
    if nargin < 3, latBounds = []; end
    if nargin < 4, lonBounds = []; end
    if nargin < 5 || isempty(method), method = 'band'; end
    if nargin < 6 || isempty(normalize), normalize = true; end
    m = lower(method); isnone = any(strcmp(m, {'none', 'equal', 'unweighted'}));
    if ~isvector(lat)
        if isnone, W = ones(size(lat)); else, W = max(cos(double(lat) * pi / 180), 0); end
    else
        if isempty(lon), error('gaw:lon', 'lon is required to build a 2-D weight field'); end
        if isnone
            W = ones(numel(lat), numel(lon));
        else
            w1 = gaw_area_weights(lat, latBounds, method, false);
            n = numel(lon);
            if isempty(lonBounds)
                b = gaw_infer_bounds(lon); ed = [b(1:n), b(2:n+1)];
            else
                ed = gaw_edges(lonBounds, n, 'lonBounds');
            end
            W = w1(:) * abs(ed(:, 2) - ed(:, 1))';
        end
    end
    if normalize, W = W / sum(W(:)); end
end
