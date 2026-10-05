function ref_octave(fname)
%REF_OCTAVE  Reference values from the MATLAB/Octave implementation (see compare.py).
%   octave-cli --eval "addpath('../matlab'); ref_octave('ref_octave.csv')"
lat = (-90:2.5:90)'; lon = 0:2.5:357.5; nt = 5; nl = numel(lon); nla = numel(lat);
x = zeros(nl, nla, nt);
for i = 1:nl, for j = 1:nla, for t = 1:nt
    if mod(i*7 + j*3 + t, 11) == 0, x(i,j,t) = NaN; else, x(i,j,t) = sin(0.3*i + 0.17*j + 1.1*t) + 0.01*j; end
end, end, end
fid = fopen(fname, 'w');
put = @(k, a) arrayfun(@(i) fprintf(fid, '%s[%d],%.15e\n', k, i-1, a(i)), 1:numel(a));
put('w_cos', gaw_area_weights(lat, [], 'cos')); put('w_band', gaw_area_weights(lat)); put('w_ell', gaw_area_weights(lat, [], 'ellipsoid'));
A = gaw_cell_area(lat, lon); put('area_sum_sphere', sum(A(:)));
A = gaw_cell_area(lat, lon, [], [], true); put('area_sum_ell', sum(A(:)));
ms = {'cos','band','ellipsoid','none'};
for k = 1:4, put(['mean_' ms{k}], gaw_area_mean(x, lat, lon, 'Method', ms{k})); end
mk = false(nla, nl); mk(21:50, 11:90) = true;
put('mean_mask', gaw_area_mean(x, lat, lon, 'Mask', mk)); put('mean_skipna_false', gaw_area_mean(x, lat, lon, 'SkipNaN', false));
x0 = x(:,:,1); x0(isnan(x0)) = 0;
put('integral', gaw_area_integral(x0, lat, lon) / 1e12);
lat2 = [85 60 33 10 -12 -41 -70 -88]'; lon2 = [0 7 30 31 100 200 300];
x2 = zeros(7, 8); for i = 1:7, for j = 1:8, x2(i,j) = cos(0.5*i + 0.9*j) + 0.1*j; end, end
W = gaw_weights_2d(lat2, lon2); put('w2d_irreg', W');   % python row-major (lat, lon)  -> print W transposed col-major
put('mean_irreg_band', gaw_area_mean(x2, lat2, lon2)); put('mean_irreg_cos', gaw_area_mean(x2, lat2, lon2, 'Method','cos'));
put('mean_irreg_ell', gaw_area_mean(x2, lat2, lon2, 'Method','ellipsoid'));
c = double(x(:,:,1) > 0.5); c(isnan(x(:,:,1))) = NaN; put('frac', gaw_area_fraction(c, lat, lon));
fclose(fid);
end
