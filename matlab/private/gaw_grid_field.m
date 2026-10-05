function f = gaw_grid_field(f, nlat, nlon, name, kind)
%GAW_GRID_FIELD  Validate a mask / weight field against the nlat x nlon grid. (internal)
%   An nlon x nlat field is transposed when nlat ~= nlon. Masks: false/0/NaN exclude.
%   Weights: NaN -> 0, negative -> error.
    if ~isequal(size(f), [nlat nlon])
        if nlat ~= nlon && isequal(size(f), [nlon nlat])
            f = f.';
        elseif isscalar(f)
            f = repmat(f, nlat, nlon);
        else
            error('gaw:size', '%s is %s but the grid is nlat x nlon = %d x %d', ...
                  name, mat2str(size(f)), nlat, nlon);
        end
    end
    if strcmp(kind, 'mask')
        if islogical(f), return; end
        f = double(f);
        f = isfinite(f) & f ~= 0;
    else
        f = double(f);
        if any(f(isfinite(f)) < 0), error('gaw:weights', '%s must be non-negative', name); end
        f(~isfinite(f)) = 0;
    end
end
