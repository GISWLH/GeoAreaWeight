\name{geoareaweight}
\alias{geoareaweight}
\alias{geoareaweight-package}
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
\title{Area Weights and Area-Weighted Statistics for Latitude-Longitude Grids}
\description{
On a regular latitude-longitude grid every row has the same number of cells,
but a cell's area shrinks with the cosine of latitude. A plain \code{mean()}
therefore over-weights high latitudes. These functions compute correct weights
and area-weighted statistics.

\itemize{
  \item \code{area_mean}: area-weighted mean (intensive quantities such as
    temperature or precipitation rate); other dimensions (time, ...) are kept.
  \item \code{region_mean}: \code{area_mean} inside a mask.
  \item \code{area_integral}: \code{sum(x * cell_area)}, e.g. a flux per unit
    area to a total.
  \item \code{area_fraction}: fraction of the valid area where a condition
    holds (e.g. drought area).
  \item \code{cell_area}, \code{weights_2d}, \code{area_weights},
    \code{band_weights}, \code{cos_lat_weights}: the weights themselves.
}

Array convention: data arrays follow the \pkg{ncdf4} order
\code{lon x lat x time} by default (\code{lon_dim = 1}, \code{lat_dim = 2}).
2-D weight, area and mask fields are \code{nlat x nlon}; an \code{nlon x nlat}
field is transposed automatically when \code{nlat != nlon}.
}
\usage{
area_mean(x, lat = NULL, lon = NULL, method = "band", weights = NULL, mask = NULL,
          lat_bounds = NULL, lon_bounds = NULL, lat_dim = 2, lon_dim = 1, skipna = TRUE)
region_mean(x, lat = NULL, lon = NULL, mask, ...)
area_integral(x, lat = NULL, lon = NULL, ellipsoid = FALSE, mask = NULL, lat_bounds = NULL,
              lon_bounds = NULL, lat_dim = 2, lon_dim = 1, units = "m2",
              radius = NULL, weights = NULL, skipna = TRUE)
area_fraction(cond, lat = NULL, lon = NULL, valid = NULL, ...)
area_weights(lat, lat_bounds = NULL, method = "band", normalize = TRUE)
cos_lat_weights(lat)
band_weights(lat, lat_bounds = NULL, ellipsoid = FALSE, normalize = TRUE)
cell_area(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL, ellipsoid = FALSE,
          radius = NULL, units = "m2")
weights_2d(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL, method = "band",
           normalize = TRUE)
infer_bounds(centers, lo = -Inf, hi = Inf)
authalic_sin_lat(lat_deg, e2 = WGS84_E2)
authalic_radius(a = WGS84_A, e2 = WGS84_E2)
}
\arguments{
\item{x}{numeric array or matrix with a latitude and a longitude dimension.}
\item{lat, lon}{cell-centre latitude / longitude in degrees. 1-D vectors
  (ascending, descending or irregular; longitudes may cross 0/360 or +-180).
  A matrix \code{lat} (curvilinear grid) falls back to cos(lat) with a warning.}
\item{method}{\code{"band"} (exact spherical cell area, default), \code{"cos"}
  (cos of the centre latitude), \code{"ellipsoid"} (exact WGS84 cell area) or
  \code{"none"} (plain unweighted mean, for comparison). Aliases:
  \code{"sphere"}, \code{"wgs84"}, \code{"unweighted"}.}
\item{weights}{optional \code{nlat x nlon} matrix of your own cell areas (e.g.
  CMIP \code{areacella}); overrides \code{method} (and the computed areas in
  \code{area_integral}). NA weights count as 0.}
\item{mask}{optional \code{nlat x nlon} logical or numeric matrix;
  \code{TRUE} / non-zero = include, \code{FALSE} / 0 / NA = exclude.}
\item{skipna}{if \code{TRUE}, NA cells are skipped and the weights are
  renormalised over the valid cells of every time step. If \code{FALSE}, any
  NA inside the (masked) domain gives NA.}
\item{lat_bounds, lon_bounds}{cell edges (length n+1) or an n x 2 matrix of
  bounds; inferred from the centres (mid-points, clipped to +-90) if NULL.}
\item{lat_dim, lon_dim}{positions of the latitude / longitude dimensions in \code{x}.}
\item{ellipsoid}{use WGS84 cell areas instead of a sphere.}
\item{radius}{sphere radius in metres (default 6371008.8, the IUGG mean radius).}
\item{units}{area unit of computed cell areas: \code{"m2"} or \code{"km2"}.}
\item{cond}{logical array (or numeric 0/1 with NA for invalid cells).}
\item{valid}{optional logical array of cells that count in the denominator;
  recycled along trailing dimensions. Note that \code{spei <= -1} is NA (not
  FALSE) where \code{spei} is NA, so in R missing cells are excluded anyway.}
\item{...}{further arguments passed to \code{area_mean}.}
\item{normalize}{scale the weights to sum to 1.}
\item{centers, lo, hi}{cell centres and optional clipping limits for \code{infer_bounds}.}
\item{lat_deg, e2}{geodetic latitude (degrees) and eccentricity squared.}
\item{a}{semi-major axis in metres.}
}
\value{
\code{area_mean}, \code{region_mean}, \code{area_integral} and
\code{area_fraction} return \code{x} with the latitude and longitude
dimensions removed (a number, a vector, or an array); NA where no valid cell
remains. \code{cell_area} and \code{weights_2d} return an \code{nlat x nlon}
matrix (\code{cell_area} without \code{lon}: area per radian of longitude);
\code{area_weights}, \code{band_weights}, \code{cos_lat_weights} return one
weight per latitude; \code{infer_bounds} returns n+1 edges.
}
\details{
\code{"band"} weights are \code{|sin(lat_north) - sin(lat_south)|}, the exact
area of a spherical latitude band; \code{"ellipsoid"} replaces the sine by the
sine of the authalic latitude, which gives exact WGS84 areas with the authalic
radius 6371007.181 m. On fine regular grids \code{"cos"} and \code{"band"}
agree closely; \code{"band"} is also correct for coarse, irregular and pole-point
grids (a 90N row becomes a half-height polar cap).
}
\examples{
lat <- seq(-89.5, 89.5, 1); lon <- seq(0.5, 359.5, 1)
x <- matrix(rep(abs(lat), each = length(lon)), length(lon))   # lon x lat
area_mean(x, lat, lon)                    # 32.71 (sphere, analytic 90 - 180/pi = 32.70)
area_mean(x, lat, lon, method = "none")   # 45, the naive mean
sum(cell_area(lat, lon, units = "km2"))   # 510065881 km2 (sphere)
area_fraction(x >= 30, lat, lon)          # 0.5 = area poleward of 30 degrees
}
\keyword{spatial}
