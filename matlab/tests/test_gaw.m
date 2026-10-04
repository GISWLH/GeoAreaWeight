% Basic self-test for GeoAreaWeight (MATLAB / Octave).  Run:  cd matlab; addpath('.'); test_gaw
R = 6371008.8;
fprintf('Running GeoAreaWeight MATLAB/Octave tests ...\n');

% 1. areas sum to sphere / ellipsoid
for cfg = [1 0; 2.5 1; 5 0]'
    res = cfg(1);
    if cfg(2) == 1, lat = (-90:res:90)'; lon = 0:res:360-res; else, lat = (-90+res/2:res:90)'; lon = (res/2:res:360); end
    A = gaw_cell_area(lat, lon);
    assert(abs(sum(A(:)) / (4*pi*R^2) - 1) < 1e-12, 'sphere area');
end
lat = (-89.5:1:89.5)'; lon = (0.5:1:359.5);
Ae = gaw_cell_area(lat, lon, [], [], true);
assert(abs(sum(Ae(:)) / 1e6 / 510065621.7 - 1) < 1e-7, 'ellipsoid area');
assert(abs(gaw_authalic_radius() - 6371007.1809) < 1e-3, 'authalic radius');

% 2. descending / irregular latitude
lat2 = [80 60 30 10 -5 -40 -75]';
w = gaw_band_weights(lat2);
assert(abs(sum(w) - 1) < 1e-14);
assert(max(abs(w - flipud(gaw_band_weights(flipud(lat2))))) < 1e-14);
assert(abs(2*pi*sum(gaw_cell_area(lat2)) / (4*pi*R^2) - 1) < 1e-12);

% 3. constant field, analytic |lat| mean, unweighted mean
nl = numel(lon); nt = numel(lat);
x = 3.7 * ones(nl, nt);
ms = {'cos', 'band', 'ellipsoid', 'none'};
for k = 1:4, assert(abs(gaw_area_mean(x, lat, lon, 'Method', ms{k}) - 3.7) < 1e-10); end
y = repmat(abs(lat'), nl, 1);
assert(abs(gaw_area_mean(y, lat, lon) - (pi/2 - 1)*180/pi) < 3e-3);
assert(abs(gaw_area_mean(y, lat, lon, 'Method', 'none') - 45) < 1e-9);

% 4. NaN renormalisation, time dimension, mask
lat5 = (-87.5:5:87.5)'; lon5 = (2.5:5:357.5); n1 = numel(lon5); n2 = numel(lat5);
rng_x = reshape(sin((1:n1*n2*6) * 0.37), n1, n2, 6);
rng_x(:, 1:5, :) = NaN; rng_x(10:30, 10:20, 4) = NaN;
m = gaw_area_mean(rng_x, lat5, lon5);
assert(isequal(size(m), [6 1]));
W = gaw_weights_2d(lat5, lon5, [], [], 'band', false);   % lat x lon
Wt = W';
for t = 1:6
    xt = rng_x(:, :, t); v = isfinite(xt);
    ref = sum(xt(v) .* Wt(v)) / sum(Wt(v));
    assert(abs(m(t) - ref) < 1e-12, 'NaN renorm');
end
assert(all(isnan(gaw_area_mean(rng_x, lat5, lon5, 'SkipNaN', false))));
mk = false(n2, n1); mk(20:30, 5:40) = true;
mm = gaw_area_mean(rng_x, lat5, lon5, 'Mask', mk);
xt = rng_x(:, :, 1); Wm = Wt .* mk'; v = isfinite(xt);
assert(abs(mm(1) - sum(xt(v) .* Wm(v)) / sum(Wm(v))) < 1e-12, 'mask');

% 5. area fraction & integral
cond = repmat((lat' >= 30), nl, 1);
assert(abs(gaw_area_fraction(cond, lat, lon) - 0.25) < 1e-12);
assert(abs(gaw_area_integral(ones(nl, nt), lat, lon) / (4*pi*R^2) - 1) < 1e-12);

% 6. user-supplied weights equal built-in
A = gaw_cell_area(lat5, lon5);
assert(abs(gaw_area_mean(rng_x(:,:,1), [], [], 'Weights', A) - gaw_area_mean(rng_x(:,:,1), lat5, lon5)) < 1e-12);
fprintf('All GeoAreaWeight MATLAB/Octave tests passed.\n');
