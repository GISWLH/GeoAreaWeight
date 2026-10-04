% needs: addpath('../matlab')
f = 'air.mon.mean.nc';
x = ncread(f, 'air'); lat = ncread(f, 'lat'); lon = ncread(f, 'lon');   % lon x lat x time
gm    = gaw_area_mean(x, lat, lon);                                      % area-weighted
naive = squeeze(mean(mean(x, 1), 2));                                    % the mistake
fprintf('mean: unweighted %.2f  weighted %.2f\n', mean(naive), mean(gm));
