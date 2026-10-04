function out = gaw_core(x, W, latDim, lonDim, mask, skipnan)
%GAW_CORE  weighted mean over (latDim, lonDim) with renormalised weights. (internal)
    d = size(x); nd = max(numel(d), max(latDim, lonDim));
    d(end+1:nd) = 1;
    rest = setdiff(1:nd, [latDim lonDim]);
    xa = permute(double(x), [latDim lonDim rest]);
    nspace = d(latDim) * d(lonDim);
    nrest = prod(d(rest)); if isempty(rest), nrest = 1; end
    xa = reshape(xa, nspace, nrest);
    if ~isequal(size(W), [d(latDim) d(lonDim)]), error('gaw:size', 'weights do not match data grid'); end
    w = W(:);
    if ~isempty(mask), w = w .* double(mask(:) ~= 0); end
    Wm = repmat(w, 1, nrest);
    if skipnan, valid = isfinite(xa); else, valid = true(nspace, nrest); end
    Wv = Wm .* valid;
    xa(~valid) = 0;
    num = sum(xa .* Wv, 1); den = sum(Wv, 1);
    res = num ./ den; res(den <= 0) = NaN;
    if ~skipnan
        bad = sum(~isfinite(reshape(permute(double(x), [latDim lonDim rest]), nspace, nrest)) & Wm > 0, 1) > 0;
        res(bad) = NaN;
    end
    if isempty(rest)
        out = res;
    elseif numel(rest) == 1
        out = res(:);
    else
        out = reshape(res, d(rest));
    end
end
