function out = gaw_region_mean(x, lat, lon, mask, varargin)
%GAW_REGION_MEAN  Area-weighted mean inside mask (nlat x nlon, true = in region).
%   m = gaw_region_mean(x, lat, lon, mask, ...)  is gaw_area_mean(x, lat, lon, 'Mask', mask, ...)
%
%   See also GAW_AREA_MEAN.
    if isempty(mask), error('gaw:args', 'mask is required'); end
    out = gaw_area_mean(x, lat, lon, 'Mask', mask, varargin{:});
end
