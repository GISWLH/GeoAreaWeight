# GeoAreaWeight (R) -- mirrors the Python API. Angles in degrees.
# Array convention: ncdf4 style, defaults lon_dim = 1, lat_dim = 2.

WGS84_A  <- 6378137.0
WGS84_F  <- 1 / 298.257223563
WGS84_E2 <- WGS84_F * (2 - WGS84_F)
EARTH_RADIUS <- 6371008.8   # IUGG mean radius R1 [m]

.deg2rad <- function(x) x * pi / 180

#' Cell edges (n+1) from cell centres; mid-points, outer edges extrapolated by half a step.
infer_bounds <- function(centers, lo = -Inf, hi = Inf) {
  c <- as.numeric(centers); n <- length(c)
  if (n < 2) stop("need >= 2 grid points to infer cell bounds; pass explicit bounds")
  b <- numeric(n + 1)
  b[2:n] <- 0.5 * (c[-n] + c[-1])
  b[1] <- c[1] - 0.5 * (c[2] - c[1])
  b[n + 1] <- c[n] + 0.5 * (c[n] - c[n - 1])
  pmin(pmax(b, lo), hi)
}

.edges <- function(bounds, n, name) {
  b <- bounds
  if (is.matrix(b) && all(dim(b) == c(n, 2))) return(b)
  if (is.null(dim(b)) && length(b) == n + 1) return(cbind(b[-(n + 1)], b[-1]))
  stop(name, " must have length n+1 or be an n x 2 matrix")
}

#' cos(lat) weights (clipped at 0)
cos_lat_weights <- function(lat) pmax(cos(.deg2rad(as.numeric(lat))), 0)

#' sin(authalic latitude) of geodetic latitude(s)
authalic_sin_lat <- function(lat_deg, e2 = WGS84_E2) {
  s <- sin(.deg2rad(lat_deg)); e <- sqrt(e2)
  q <- function(sn) (1 - e2) * (sn / (1 - e2 * sn^2) - log((1 - e * sn) / (1 + e * sn)) / (2 * e))
  r <- q(s) / q(1)
  if (!is.null(dim(lat_deg))) dim(r) <- dim(lat_deg)
  r
}

#' Authalic radius of the WGS84 ellipsoid [m]
authalic_radius <- function(a = WGS84_A, e2 = WGS84_E2) {
  e <- sqrt(e2)
  qp <- 1 + (1 - e2) / e * atanh(e)
  a * sqrt(qp / 2)
}

#' Exact area weight of each latitude band |sin(phi_n) - sin(phi_s)|
band_weights <- function(lat, lat_bounds = NULL, ellipsoid = FALSE, normalize = TRUE) {
  lat <- as.numeric(lat); n <- length(lat)
  ed <- if (is.null(lat_bounds)) {
    b <- infer_bounds(lat, -90, 90); cbind(b[-(n + 1)], b[-1])
  } else .edges(lat_bounds, n, "lat_bounds")
  s <- if (ellipsoid) authalic_sin_lat(ed) else sin(.deg2rad(ed))
  w <- abs(s[, 2] - s[, 1])
  if (normalize) w <- w / sum(w)
  w
}

#' 1-D latitude weights: method = "cos", "band" (default), "ellipsoid", "none"
area_weights <- function(lat, lat_bounds = NULL, method = "band", normalize = TRUE) {
  m <- tolower(method)
  if (m == "cos") w <- cos_lat_weights(lat)
  else if (m %in% c("band", "sphere", "sphere_band")) return(band_weights(lat, lat_bounds, FALSE, normalize))
  else if (m %in% c("ellipsoid", "wgs84")) return(band_weights(lat, lat_bounds, TRUE, normalize))
  else if (m %in% c("none", "equal", "unweighted")) w <- rep(1, length(lat))
  else stop("unknown method: ", method)
  if (normalize) w / sum(w) else w
}

#' Cell area [m2 or km2]; matrix (nlat x nlon) if lon given, else vector per radian of longitude
cell_area <- function(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL,
                      ellipsoid = FALSE, radius = NULL, units = "m2") {
  R <- if (ellipsoid) authalic_radius() else if (is.null(radius)) EARTH_RADIUS else radius
  f <- band_weights(lat, lat_bounds, ellipsoid, normalize = FALSE)
  if (is.null(lon)) {
    area <- R^2 * f
  } else {
    n <- length(lon)
    ed <- if (is.null(lon_bounds)) { b <- infer_bounds(lon); cbind(b[-(n + 1)], b[-1]) }
          else .edges(lon_bounds, n, "lon_bounds")
    dlon <- .deg2rad(abs(ed[, 2] - ed[, 1]))
    area <- R^2 * outer(f, dlon)
  }
  if (units == "km2") area <- area / 1e6 else if (units != "m2") stop("units must be 'm2' or 'km2'")
  area
}

#' 2-D weight field (nlat x nlon). 2-D `lat` (curvilinear): cos(lat) only.
weights_2d <- function(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL,
                       method = "band", normalize = TRUE) {
  m <- tolower(method)
  if (is.matrix(lat)) {
    w <- if (m %in% c("none", "equal", "unweighted")) matrix(1, nrow(lat), ncol(lat)) else pmax(cos(.deg2rad(lat)), 0)
  } else {
    if (is.null(lon)) stop("lon is required to build a 2-D weight field")
    if (m %in% c("none", "equal", "unweighted")) {
      w <- matrix(1, length(lat), length(lon))
    } else {
      w1 <- area_weights(lat, lat_bounds, method, normalize = FALSE)
      n <- length(lon)
      ed <- if (is.null(lon_bounds)) { b <- infer_bounds(lon); cbind(b[-(n + 1)], b[-1]) }
            else .edges(lon_bounds, n, "lon_bounds")
      w <- outer(w1, abs(ed[, 2] - ed[, 1]))
    }
  }
  if (normalize) w <- w / sum(w)
  w
}

# core: x has dims; returns weighted mean over (lat_dim, lon_dim), weights renormalised over valid cells
.core <- function(x, w2d, lat_dim, lon_dim, mask, skipna) {
  d <- dim(x); if (is.null(d)) stop("x must be an array/matrix")
  nd <- length(d)
  rest <- setdiff(seq_len(nd), c(lat_dim, lon_dim))
  xa <- aperm(x, c(lat_dim, lon_dim, rest))
  nspace <- d[lat_dim] * d[lon_dim]
  nrest <- if (length(rest)) prod(d[rest]) else 1
  dim(xa) <- c(nspace, nrest)
  if (!all(dim(w2d) == c(d[lat_dim], d[lon_dim]))) stop("weights do not match the data grid")
  w <- as.vector(w2d)
  if (!is.null(mask)) w <- w * as.vector(mask != 0)
  W <- matrix(w, nspace, nrest)
  valid <- if (skipna) is.finite(xa) else matrix(TRUE, nspace, nrest)
  Wv <- W * valid
  xa0 <- xa; xa0[!valid] <- 0
  num <- colSums(xa0 * Wv); den <- colSums(Wv)
  out <- ifelse(den > 0, num / den, NA_real_)
  if (!skipna) out[colSums(!is.finite(xa) & W > 0) > 0] <- NA_real_
  if (length(rest)) array(out, dim = d[rest]) else as.numeric(out)
}

.make_w <- function(lat, lon, method, weights, lat_bounds, lon_bounds) {
  if (!is.null(weights)) return(as.matrix(weights))
  weights_2d(lat, lon, lat_bounds, lon_bounds, method, normalize = FALSE)
}

#' Area-weighted spatial mean of an array; time/other dims are kept.
#' @param x array (default layout lon x lat x time, as read by ncdf4)
#' @param weights optional nlat x nlon cell-area matrix (e.g. areacella)
#' @param mask optional nlat x nlon logical/0-1 matrix (TRUE = include)
area_mean <- function(x, lat = NULL, lon = NULL, method = "band", weights = NULL, mask = NULL,
                      lat_bounds = NULL, lon_bounds = NULL, lat_dim = 2, lon_dim = 1, skipna = TRUE) {
  if (is.null(lat) && is.null(weights)) stop("lat (or weights) is required")
  w <- .make_w(lat, lon, method, weights, lat_bounds, lon_bounds)
  .core(x, w, lat_dim, lon_dim, mask, skipna)
}

region_mean <- function(x, lat = NULL, lon = NULL, mask, ...) area_mean(x, lat, lon, mask = mask, ...)

#' Sum of x * cell_area (e.g. flux per area -> total)
area_integral <- function(x, lat, lon, ellipsoid = FALSE, mask = NULL, lat_bounds = NULL,
                          lon_bounds = NULL, lat_dim = 2, lon_dim = 1, units = "m2") {
  a <- cell_area(lat, lon, lat_bounds, lon_bounds, ellipsoid = ellipsoid, units = units)
  if (!is.null(mask)) a <- a * (mask != 0)
  d <- dim(x); nd <- length(d); rest <- setdiff(seq_len(nd), c(lat_dim, lon_dim))
  xa <- aperm(x, c(lat_dim, lon_dim, rest)); dim(xa) <- c(d[lat_dim] * d[lon_dim], if (length(rest)) prod(d[rest]) else 1)
  xa[!is.finite(xa)] <- 0
  out <- colSums(xa * as.vector(a))
  if (length(rest)) array(out, dim = d[rest]) else as.numeric(out)
}

#' Area fraction (0-1) of cells where `cond` is TRUE among valid cells
area_fraction <- function(cond, lat = NULL, lon = NULL, valid = NULL, ...) {
  c0 <- array(as.numeric(cond), dim = dim(cond))
  if (!is.null(valid)) c0[!valid] <- NA_real_
  area_mean(c0, lat, lon, ...)
}
