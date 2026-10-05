% GeoAreaWeight - area weights and area-weighted statistics for lat-lon grids
% Version 0.1.0
%
% Statistics (x is lon x lat x ... as returned by ncread; see 'LatDim'/'LonDim')
%   gaw_area_mean       - Area-weighted mean (intensive quantities), other dims kept
%   gaw_region_mean     - Area-weighted mean inside a mask
%   gaw_area_integral   - sum(x .* cellArea): per-area flux -> total
%   gaw_area_fraction   - Fraction of the valid area where a condition holds
%
% Weights and areas (2-D fields are nlat x nlon)
%   gaw_cell_area       - Exact cell areas [m2 | km2], sphere or WGS84
%   gaw_weights_2d      - 2-D weight field
%   gaw_area_weights    - 1-D latitude weights ('band' | 'cos' | 'ellipsoid' | 'none')
%   gaw_band_weights    - Exact band weights |sin(phi_n) - sin(phi_s)|
%   gaw_cos_lat_weights - cos(latitude) weights
%   gaw_infer_bounds    - Cell edges from cell centres
%   gaw_authalic_sin_lat - sin(authalic latitude), WGS84
%   gaw_authalic_radius - Authalic radius of WGS84 [m]
%
% Tested with GNU Octave; run tests/test_gaw.m.
