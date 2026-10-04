function R = gaw_authalic_radius()
%GAW_AUTHALIC_RADIUS  Radius [m] of the sphere with the same area as the WGS84 ellipsoid.
    a = 6378137.0; f = 1 / 298.257223563; e2 = f * (2 - f); e = sqrt(e2);
    qp = 1 + (1 - e2) / e * atanh(e);
    R = a * sqrt(qp / 2);
end
