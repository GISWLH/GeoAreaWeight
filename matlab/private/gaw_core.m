function out = gaw_core(x, W, latDim, lonDim, mask, skipnan, how)
%GAW_CORE  Weighted mean or sum over (latDim, lonDim); other dimensions are kept. (internal)
    if nargin < 7, how = 'mean'; end
    if latDim == lonDim, error('gaw:dims', 'LatDim and LonDim must differ'); end
    d = size(x); nd = max(numel(d), max(latDim, lonDim));
    d(end+1:nd) = 1;
    nlat = d(latDim); nlon = d(lonDim);
    if ~isequal(size(W), [nlat nlon])
        error('gaw:size', ['lat/lon give a %s grid but the data has %d x %d along ' ...
              '(LatDim=%d, LonDim=%d); check LatDim/LonDim'], mat2str(size(W)), nlat, nlon, latDim, lonDim);
    end
    rest = setdiff(1:nd, [latDim lonDim]);
    nrest = prod(d(rest));
    xa = reshape(permute(double(x), [latDim lonDim rest]), nlat * nlon, nrest);
    w = W(:);
    if ~isempty(mask), w = w .* double(reshape(gaw_grid_field(mask, nlat, nlon, 'Mask', 'mask'), [], 1)); end
    finite = isfinite(xa);
    used = w > 0;
    xa(~finite) = 0;
    num = w' * xa;
    nvalid = double(used') * finite;
    if strcmp(how, 'mean')
        if skipnan, den = w' * double(finite); else, den = repmat(sum(w), 1, nrest); end
        res = num ./ den; res(den <= 0) = NaN;
    else
        res = num;
    end
    res(nvalid == 0) = NaN;
    if ~skipnan, res(nvalid < sum(used)) = NaN; end
    if isempty(rest)
        out = res;
    elseif numel(rest) == 1
        out = res(:);
    else
        out = reshape(res, d(rest));
    end
end
