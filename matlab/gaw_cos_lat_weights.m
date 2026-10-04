function w = gaw_cos_lat_weights(lat)
%GAW_COS_LAT_WEIGHTS  cos(latitude) weights, clipped at 0.
    w = max(cos(double(lat) * pi / 180), 0);
end
