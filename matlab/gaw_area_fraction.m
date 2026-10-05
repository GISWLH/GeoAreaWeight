function out = gaw_area_fraction(cond, lat, lon, varargin)
%GAW_AREA_FRACTION  Fraction (0-1) of the valid area where cond is true.
%   f = gaw_area_fraction(spei <= -1, lat, lon, 'Valid', isfinite(spei))
%   cond: logical, or double 0/1 with NaN for invalid cells.
%   Note that (spei <= -1) is false, not NaN, where spei is NaN: pass 'Valid'
%   so that missing cells are excluded from the denominator.
%   'Valid' is broadcast along trailing dimensions (e.g. a lon x lat field for
%   a lon x lat x time cond). Other options are passed to GAW_AREA_MEAN.
%
%   See also GAW_AREA_MEAN.
    c = double(cond);
    k = find(strcmpi('Valid', varargin(1:2:end)), 1);
    if ~isempty(k)
        v = logical(varargin{2*k});
        varargin(2*k-1:2*k) = [];
        v = repmat(v, size(c) ./ [size(v) ones(1, ndims(c) - ndims(v))]);
        c(~v) = NaN;
    end
    out = gaw_area_mean(c, lat, lon, varargin{:});
end
