\name{geoareaweight}
\alias{geoareaweight}
\alias{area_mean}
\alias{region_mean}
\alias{area_integral}
\alias{area_fraction}
\alias{area_weights}
\alias{cos_lat_weights}
\alias{band_weights}
\alias{cell_area}
\alias{weights_2d}
\alias{infer_bounds}
\alias{authalic_sin_lat}
\alias{authalic_radius}
\title{Area weights and area-weighted means for latitude-longitude grids}
\description{
Cosine-latitude, exact spherical band and WGS84 authalic cell-area weights, and
NaN-aware area-weighted means, regional means, integrals and area fractions.
Arrays follow the ncdf4 convention (lon x lat x time) by default (\code{lon_dim = 1}, \code{lat_dim = 2}).
}
\usage{
area_mean(x, lat = NULL, lon = NULL, method = "band", weights = NULL, mask = NULL,
          lat_bounds = NULL, lon_bounds = NULL, lat_dim = 2, lon_dim = 1, skipna = TRUE)
region_mean(x, lat = NULL, lon = NULL, mask, ...)
area_integral(x, lat, lon, ellipsoid = FALSE, mask = NULL, lat_bounds = NULL,
              lon_bounds = NULL, lat_dim = 2, lon_dim = 1, units = "m2")
area_fraction(cond, lat = NULL, lon = NULL, valid = NULL, ...)
area_weights(lat, lat_bounds = NULL, method = "band", normalize = TRUE)
cos_lat_weights(lat)
band_weights(lat, lat_bounds = NULL, ellipsoid = FALSE, normalize = TRUE)
cell_area(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL, ellipsoid = FALSE,
          radius = NULL, units = "m2")
weights_2d(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL, method = "band", normalize = TRUE)
infer_bounds(centers, lo = -Inf, hi = Inf)
authalic_sin_lat(lat_deg, e2 = WGS84_E2)
authalic_radius(a = WGS84_A, e2 = WGS84_E2)
}
\arguments{
\item{x}{numeric array with a latitude and a longitude dimension.}
\item{lat, lon}{latitude / longitude in degrees (1-D; 2-D lat falls back to cos(lat)).}
\item{method}{"band" (exact spherical band area, default), "cos", "ellipsoid" (WGS84), "none".}
\item{weights}{optional nlat x nlon cell-area matrix (e.g. CMIP areacella); overrides method.}
\item{mask}{optional nlat x nlon logical matrix (TRUE = include).}
\item{skipna}{weights are renormalised over finite cells.}
\item{lat_bounds, lon_bounds}{cell edges (n+1) or n x 2 matrix; inferred from centres if NULL.}
\item{lat_dim, lon_dim}{integer positions of the latitude / longitude dimensions in x.}
\item{ellipsoid}{use WGS84 authalic areas instead of a sphere.}
\item{cond}{logical array; \code{valid} optional logical array of valid cells.}
\item{valid}{optional logical array: cells to include (others set to NA).}
\item{...}{further arguments passed to \code{area_mean}.}
\item{units}{"m2" or "km2".}
\item{normalize}{scale weights to sum to 1.}
\item{radius}{sphere radius in metres (default 6371008.8).}
\item{centers, lo, hi}{cell centres and optional clipping limits for \code{infer_bounds}.}
\item{lat_deg, e2}{geodetic latitude (degrees) and eccentricity squared.}
\item{a}{semi-major axis in metres.}
}
\value{An array with the lat/lon dimensions removed (area_mean), or weights / areas.}
\examples{
lat <- seq(-89.5, 89.5, 1); lon <- seq(0.5, 359.5, 1)
x <- matrix(rep(abs(lat), each = length(lon)), length(lon))
area_mean(x, lat, lon)                    # ~32.7  (analytic: (pi/2 - 1) rad = 32.70 deg)
area_mean(x, lat, lon, method = "none")   # 45, the naive mean
}
