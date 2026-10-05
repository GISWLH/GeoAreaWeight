function b = gaw_infer_bounds(centers, lo, hi)
%GAW_INFER_BOUNDS  Cell edges (n+1 x 1) from n strictly monotonic cell centres.
%   b = gaw_infer_bounds(centers)            mid-points; outer edges extrapolated by half a step
%   b = gaw_infer_bounds(centers, lo, hi)    additionally clipped to [lo, hi]
    c = double(centers(:)); n = numel(c);
    if n < 2
        error('gaw:bounds', 'need >= 2 grid points to infer cell bounds; pass explicit bounds');
    end
    if any(~isfinite(c)), error('gaw:bounds', 'centers contain NaN or infinite values'); end
    d = diff(c);
    if ~(all(d > 0) || all(d < 0))
        error('gaw:bounds', ['centers must be strictly increasing or strictly decreasing to ' ...
              'infer cell bounds; sort the grid or pass explicit bounds']);
    end
    b = zeros(n + 1, 1);
    b(2:n) = 0.5 * (c(1:end-1) + c(2:end));
    b(1)   = c(1) - 0.5 * (c(2) - c(1));
    b(n+1) = c(n) + 0.5 * (c(n) - c(n-1));
    if nargin >= 2 && ~isempty(lo), b = max(b, lo); end
    if nargin >= 3 && ~isempty(hi), b = min(b, hi); end
end
