function r = gaw_authalic_sin_lat(latDeg)
%GAW_AUTHALIC_SIN_LAT  sin(authalic latitude) for WGS84 geodetic latitude(s) in degrees.
%   Cell area on the ellipsoid = Rq^2 * dlon * (r(north) - r(south)), Rq = GAW_AUTHALIC_RADIUS.
    f = 1 / 298.257223563; e2 = f * (2 - f); e = sqrt(e2);
    s = sin(double(latDeg) * pi / 180);
    q = @(sn) (1 - e2) .* (sn ./ (1 - e2 .* sn.^2) - log((1 - e .* sn) ./ (1 + e .* sn)) ./ (2 * e));
    r = q(s) ./ q(1);
end
