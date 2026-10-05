function ed = gaw_edges(bounds, n, name)
%GAW_EDGES  (n+1) edge vector or (n x 2) bounds -> n x 2 matrix. (internal)
    if isvector(bounds) && numel(bounds) == n + 1
        b = double(bounds(:)); ed = [b(1:n), b(2:n+1)];
    elseif ismatrix(bounds) && isequal(size(bounds), [n 2])
        ed = double(bounds);
    else
        error('gaw:bounds', '%s must have n+1 elements or be n x 2 (n = %d)', name, n);
    end
    if any(~isfinite(ed(:))), error('gaw:bounds', '%s contains NaN or infinite values', name); end
end
