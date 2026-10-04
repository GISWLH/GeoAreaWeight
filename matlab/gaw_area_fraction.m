function out = gaw_area_fraction(cond, lat, lon, varargin)
%GAW_AREA_FRACTION  Area fraction (0-1) of cells where cond is true, among valid cells.
%   Use NaN in a double cond to exclude cells, e.g. gaw_area_fraction(double(spei < -1) + 0*spei, lat, lon).
    out = gaw_area_mean(double(cond), lat, lon, varargin{:});
end
