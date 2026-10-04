function out = gaw_area_integral(x, lat, lon, varargin)
%GAW_AREA_INTEGRAL  sum(x .* cellArea) over the grid (flux per area -> total). NaN skipped.
%   options: 'Ellipsoid' (false), 'Mask', 'LatBounds', 'LonBounds', 'LatDim' (2), 'LonDim' (1), 'Units' ('m2')
    o = struct('Ellipsoid', false, 'Mask', [], 'LatBounds', [], 'LonBounds', [], 'LatDim', 2, 'LonDim', 1, 'Units', 'm2');
    for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
    A = gaw_cell_area(lat, lon, o.LatBounds, o.LonBounds, o.Ellipsoid, [], o.Units);
    if ~isempty(o.Mask), A = A .* double(o.Mask ~= 0); end
    d = size(x); nd = max(numel(d), max(o.LatDim, o.LonDim)); d(end+1:nd) = 1;
    rest = setdiff(1:nd, [o.LatDim o.LonDim]);
    xa = permute(double(x), [o.LatDim o.LonDim rest]);
    nrest = prod(d(rest)); if isempty(rest), nrest = 1; end
    xa = reshape(xa, d(o.LatDim) * d(o.LonDim), nrest);
    xa(~isfinite(xa)) = 0;
    res = sum(xa .* repmat(A(:), 1, nrest), 1);
    if isempty(rest), out = res; elseif numel(rest) == 1, out = res(:); else, out = reshape(res, d(rest)); end
end
