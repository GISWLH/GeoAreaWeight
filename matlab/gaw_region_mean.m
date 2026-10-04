function out = gaw_region_mean(x, lat, lon, mask, varargin)
%GAW_REGION_MEAN  Area-weighted mean inside mask (nlat x nlon, true = in region).
    out = gaw_area_mean(x, lat, lon, 'Mask', mask, varargin{:});
end
