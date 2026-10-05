function w = gaw_cos_lat_weights(lat)
%GAW_COS_LAT_WEIGHTS  cos(latitude) weights (same shape as lat; 0 at the poles).
%   Approximate (O(dlat^2)); prefer GAW_BAND_WEIGHTS for exact cell areas.
    w = max(cos(gaw_check_lat(lat) * pi / 180), 0);
end
