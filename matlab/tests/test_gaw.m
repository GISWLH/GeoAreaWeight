function test_gaw()
%TEST_GAW  Self-test for GeoAreaWeight (MATLAB / Octave).
% Run from the repository root:  addpath('matlab'); addpath('matlab/tests'); test_gaw
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

    % 7. options are case-insensitive; unknown options, methods and bad input raise
    assert(abs(gaw_area_mean(x, lat, lon, 'method', 'COS') - 3.7) < 1e-10);
    expect_error(@() gaw_area_mean(x, lat, lon, 'Methd', 'cos'), 'gaw:option');
    expect_error(@() gaw_area_weights(lat, [], 'nonsense'), 'gaw:method');
    expect_error(@() gaw_area_weights((0:10:350)'), 'gaw:lat');
    expect_error(@() gaw_band_weights([0 10 5]), 'gaw:bounds');
    expect_error(@() gaw_cell_area([0 10], 5), 'gaw:lon');

    % 8. area integral: NaN skipped, NaN when no valid cell
    xi = ones(n1, n2, 2); xi(:, :, 2) = NaN; xi(1, 1, 1) = NaN;
    A5 = gaw_cell_area(lat5, lon5);
    ri = gaw_area_integral(xi, lat5, lon5);
    assert(abs(ri(1) / (sum(A5(:)) - A5(1, 1)) - 1) < 1e-12 && isnan(ri(2)));
    ri = gaw_area_integral(xi, lat5, lon5, 'SkipNaN', false);
    assert(isnan(ri(1)));

    % 9. mask given as lon x lat (transposed automatically); NaN mask entries exclude
    mm2 = gaw_area_mean(rng_x, lat5, lon5, 'Mask', mk');
    assert(abs(mm2(1) - mm(1)) < 1e-12);
    mf = double(mk); mf(~mk) = NaN;
    mm3 = gaw_area_mean(rng_x, lat5, lon5, 'Mask', mf);
    assert(abs(mm3(1) - mm(1)) < 1e-12);

    % 10. longitudes crossing the 0/360 meridian
    lw = [182.5:5:357.5, 2.5:5:177.5];
    Aw = gaw_cell_area(lat5, lw); Ar = gaw_cell_area(lat5, 2.5:5:357.5);
    assert(max(abs(Aw(:) - Ar(:)) ./ Ar(:)) < 1e-12);

    % 11. area fraction with 'Valid'
    sp = sin((1:n1*n2) * 0.71); sp = reshape(sp, n1, n2); sp(:, 1:4) = NaN;
    fr = gaw_area_fraction(sp <= -0.5, lat5, lon5, 'Valid', isfinite(sp));
    c = double(sp <= -0.5); c(~isfinite(sp)) = NaN;
    assert(abs(fr - gaw_area_mean(c, lat5, lon5)) < 1e-12);

    % 12. 2-D latitude falls back to cos(lat) with a warning
    [LA, ~] = ndgrid(lat5, lon5);
    warning('off', 'gaw:curvilinear');
    W2 = gaw_weights_2d(LA, [], [], [], 'band', false);
    warning('on', 'gaw:curvilinear');
    assert(max(abs(W2(:) - cosd(LA(:)))) < 1e-14);

    % 13. SkipNaN = false on complete data equals the default
    xc = reshape(sin((1:n1*n2*3) * 0.13), n1, n2, 3);
    assert(max(abs(gaw_area_mean(xc, lat5, lon5, 'SkipNaN', false) - gaw_area_mean(xc, lat5, lon5))) < 1e-14);
    xc(1, 1, 2) = NaN;
    r = gaw_area_mean(xc, lat5, lon5, 'SkipNaN', false);
    assert(isnan(r(2)) && all(isfinite(r([1 3]))));

    fprintf('All GeoAreaWeight MATLAB/Octave tests passed.\n');
end

function expect_error(f, id)
    try
        f();
    catch err
        assert(strcmp(err.identifier, id), 'expected error %s, got %s', id, err.identifier);
        return
    end
    error('test:noerror', 'expected error %s was not raised', id);
end
