function lat = gaw_check_lat(lat, name)
%GAW_CHECK_LAT  Check that latitudes are within [-90, 90] and clip round-off. (internal)
    if nargin < 2, name = 'lat'; end
    lat = double(lat);
    f = lat(isfinite(lat));
    if ~isempty(f) && max(abs(f)) > 90 + 1e-6
        error('gaw:lat', '%s must be within [-90, 90] degrees (max |value| = %g); did you swap lat and lon?', ...
              name, max(abs(f)));
    end
    lat = min(max(lat, -90), 90);
end
