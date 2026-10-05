% Unweighted vs area-weighted global mean of NCEP/NCAR R1 near-surface air temperature.
% Data: air.mon.mean.nc from https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.surface.html
% Run from the examples folder after addpath('../matlab').
f = 'air.mon.mean.nc';
x = ncread(f, 'air'); lat = ncread(f, 'lat'); lon = ncread(f, 'lon');   % lon x lat x time
gm    = gaw_area_mean(x, lat, lon);                                      % area-weighted
naive = squeeze(mean(mean(x, 1), 2));                                    % plain mean: the mistake
fprintf('mean: unweighted %.2f  weighted %.2f degC\n', mean(naive), mean(gm));
