# GeoAreaWeight (R): area weights and area-weighted statistics for lat-lon grids.
# Mirrors the Python API. Angles are in degrees.
# Array convention: ncdf4 order, i.e. defaults lon_dim = 1, lat_dim = 2.
# 2-D weight / area / mask fields are nlat x nlon (an nlon x nlat field is
# transposed automatically when nlat != nlon).

WGS84_A  <- 6378137.0
WGS84_F  <- 1 / 298.257223563
WGS84_E2 <- WGS84_F * (2 - WGS84_F)
EARTH_RADIUS <- 6371008.8   # IUGG mean radius R1 [m]

.deg2rad <- function(x) x * pi / 180

.canonical_method <- function(method) {
  m <- tolower(as.character(method))
  alias <- c(band = "band", sphere = "band", sphere_band = "band", cos = "cos",
             ellipsoid = "ellipsoid", wgs84 = "ellipsoid",
             none = "none", equal = "none", unweighted = "none")
  if (length(m) != 1 || !(m %in% names(alias)))
    stop("unknown method '", method, "'; expected one of 'band', 'cos', 'ellipsoid', 'none'")
  unname(alias[m])
}

.check_lat <- function(lat, name = "lat") {
  lat <- as.numeric(lat)
  f <- lat[is.finite(lat)]
  if (length(f) && max(abs(f)) > 90 + 1e-6)
    stop(name, " must be within [-90, 90] degrees (max |value| = ", max(abs(f)),
         "); did you swap lat and lon?")
  pmin(pmax(lat, -90), 90)
}

.check_finite <- function(x, name) {
  if (any(!is.finite(x))) stop(name, " contains NA, NaN or infinite values")
  x
}

.unwrap_lon <- function(lon) {
  if (length(lon) < 2) return(lon)
  d <- ((diff(lon) + 180) %% 360) - 180
  c(lon[1], lon[1] + cumsum(d))
}

infer_bounds <- function(centers, lo = -Inf, hi = Inf) {
  c <- .check_finite(as.numeric(centers), "centers"); n <- length(c)
  if (n < 2) stop("need >= 2 grid points to infer cell bounds; pass explicit bounds")
  d <- diff(c)
  if (!(all(d > 0) || all(d < 0)))
    stop("centers must be strictly increasing or strictly decreasing to infer cell bounds; ",
         "sort the grid or pass explicit bounds")
  b <- numeric(n + 1)
  b[2:n] <- 0.5 * (c[-n] + c[-1])
  b[1] <- c[1] - 0.5 * (c[2] - c[1])
  b[n + 1] <- c[n] + 0.5 * (c[n] - c[n - 1])
  pmin(pmax(b, lo), hi)
}

.edges <- function(bounds, n, name) {
  b <- bounds
  out <- if (is.matrix(b) && all(dim(b) == c(n, 2))) {
    b
  } else if (is.null(dim(b)) && length(b) == n + 1) {
    cbind(b[-(n + 1)], b[-1])
  } else {
    stop(name, " must have length n+1 or be an n x 2 matrix (n = ", n, ")")
  }
  storage.mode(out) <- "double"
  .check_finite(out, name)
}

.lon_widths <- function(lon, lon_bounds, allow_single) {
  lon <- .check_finite(as.numeric(lon), "lon"); n <- length(lon)
  if (!is.null(lon_bounds)) {
    ed <- .edges(lon_bounds, n, "lon_bounds")
    w <- abs(ed[, 2] - ed[, 1])
    return(ifelse(w > 180, 360 - w, w))   # a pair such as [359, 1] crosses the wrap point
  }
  if (n == 1) {
    if (allow_single) return(1)
    stop("cannot infer the width of a single longitude; pass lon_bounds")
  }
  abs(diff(infer_bounds(.unwrap_lon(lon))))
}

cos_lat_weights <- function(lat) {
  r <- pmax(cos(.deg2rad(.check_lat(lat))), 0)
  if (!is.null(dim(lat))) dim(r) <- dim(lat)
  r
}

authalic_sin_lat <- function(lat_deg, e2 = WGS84_E2) {
  s <- sin(.deg2rad(lat_deg)); e <- sqrt(e2)
  q <- function(sn) (1 - e2) * (sn / (1 - e2 * sn^2) - log((1 - e * sn) / (1 + e * sn)) / (2 * e))
  r <- q(s) / q(1)
  if (!is.null(dim(lat_deg))) dim(r) <- dim(lat_deg)
  r
}

authalic_radius <- function(a = WGS84_A, e2 = WGS84_E2) {
  e <- sqrt(e2)
  qp <- 1 + (1 - e2) / e * atanh(e)
  a * sqrt(qp / 2)
}

band_weights <- function(lat, lat_bounds = NULL, ellipsoid = FALSE, normalize = TRUE) {
  lat <- .check_lat(.check_finite(as.numeric(lat), "lat")); n <- length(lat)
  ed <- if (is.null(lat_bounds)) {
    b <- infer_bounds(lat, -90, 90); cbind(b[-(n + 1)], b[-1])
  } else {
    matrix(.check_lat(.edges(lat_bounds, n, "lat_bounds"), "lat_bounds"), n, 2)
  }
  s <- if (ellipsoid) authalic_sin_lat(ed) else sin(.deg2rad(ed))
  w <- abs(s[, 2] - s[, 1])
  if (normalize) w <- w / sum(w)
  w
}

area_weights <- function(lat, lat_bounds = NULL, method = "band", normalize = TRUE) {
  m <- .canonical_method(method)
  if (m == "band") return(band_weights(lat, lat_bounds, FALSE, normalize))
  if (m == "ellipsoid") return(band_weights(lat, lat_bounds, TRUE, normalize))
  lat <- .check_lat(.check_finite(as.numeric(lat), "lat"))
  w <- if (m == "cos") cos_lat_weights(lat) else rep(1, length(lat))
  if (normalize) w / sum(w) else w
}

cell_area <- function(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL,
                      ellipsoid = FALSE, radius = NULL, units = "m2") {
  if (!(units %in% c("m2", "km2"))) stop("units must be 'm2' or 'km2'")
  R <- if (ellipsoid) authalic_radius() else if (is.null(radius)) EARTH_RADIUS else radius
  f <- band_weights(lat, lat_bounds, ellipsoid, normalize = FALSE)
  area <- if (is.null(lon)) R^2 * f else R^2 * outer(f, .deg2rad(.lon_widths(lon, lon_bounds, FALSE)))
  if (units == "km2") area / 1e6 else area
}

weights_2d <- function(lat, lon = NULL, lat_bounds = NULL, lon_bounds = NULL,
                       method = "band", normalize = TRUE) {
  m <- .canonical_method(method)
  if (is.matrix(lat)) {
    if (m == "none") {
      w <- matrix(1, nrow(lat), ncol(lat))
    } else {
      if (m != "cos")
        warning("method '", method, "' needs 1-D latitudes; using cos(lat) for the 2-D ",
                "(curvilinear) grid. Pass cell areas via weights= for exact results.")
      w <- cos_lat_weights(lat)
    }
  } else {
    if (is.null(lon)) stop("lon is required to build a 2-D weight field")
    lw <- area_weights(lat, lat_bounds, m, normalize = FALSE)
    w <- if (m == "none") matrix(1, length(lw), length(lon))
         else outer(lw, .lon_widths(lon, lon_bounds, TRUE))
  }
  if (normalize) w <- w / sum(w)
  w
}

# validate a mask / weight field against the nlat x nlon grid
.grid_field <- function(f, nlat, nlon, name, kind) {
  f <- as.matrix(f)
  if (!all(dim(f) == c(nlat, nlon))) {
    if (nlat != nlon && all(dim(f) == c(nlon, nlat))) f <- t(f)
    else if (length(f) == 1) f <- matrix(f, nlat, nlon)
    else stop(name, " has dimensions ", paste(dim(f), collapse = " x "),
              " but the grid is nlat x nlon = ", nlat, " x ", nlon)
  }
  if (kind == "mask") {
    if (is.logical(f)) return(!is.na(f) & f)
    f <- suppressWarnings(as.numeric(f))
    return(matrix(is.finite(f) & f != 0, nlat, nlon))
  }
  f <- matrix(as.numeric(f), nlat, nlon)
  if (any(f[is.finite(f)] < 0)) stop(name, " must be non-negative")
  f[!is.finite(f)] <- 0
  f
}

# weighted mean or sum over (lat_dim, lon_dim); other dimensions are kept
.reduce <- function(x, w2d, lat_dim, lon_dim, mask, skipna, how) {
  d <- dim(x)
  nd <- length(d)
  if (lat_dim == lon_dim || max(lat_dim, lon_dim) > nd) stop("invalid lat_dim / lon_dim")
  nlat <- d[lat_dim]; nlon <- d[lon_dim]
  if (!all(dim(w2d) == c(nlat, nlon)))
    stop("lat/lon give a ", paste(dim(w2d), collapse = " x "), " grid but the data has ",
         nlat, " x ", nlon, " along (lat_dim = ", lat_dim, ", lon_dim = ", lon_dim,
         "); check the dimension arguments")
  rest <- setdiff(seq_len(nd), c(lat_dim, lon_dim))
  xa <- aperm(x, c(lat_dim, lon_dim, rest))
  nspace <- nlat * nlon
  nrest <- if (length(rest)) prod(d[rest]) else 1
  storage.mode(xa) <- "double"
  dim(xa) <- c(nspace, nrest)
  w <- as.vector(w2d)
  if (!is.null(mask)) w <- w * as.vector(.grid_field(mask, nlat, nlon, "mask", "mask"))
  finite <- is.finite(xa)
  used <- w > 0
  x0 <- xa; x0[!finite] <- 0
  num <- colSums(x0 * w)
  n_valid <- colSums(finite & used)
  out <- if (how == "mean") {
    den <- if (skipna) colSums(finite * w) else rep(sum(w), nrest)
    ifelse(den > 0, num / den, NA_real_)
  } else num
  out[n_valid == 0] <- NA_real_
  if (!skipna) out[n_valid < sum(used)] <- NA_real_
  if (length(rest) > 1) array(out, dim = d[rest]) else as.numeric(out)
}

.data_grid <- function(x, lat_dim, lon_dim) {
  d <- dim(x)
  if (is.null(d)) stop("x must be an array or matrix with a lat and a lon dimension")
  if (max(lat_dim, lon_dim) > length(d)) stop("invalid lat_dim / lon_dim")
  c(d[lat_dim], d[lon_dim])
}

area_mean <- function(x, lat = NULL, lon = NULL, method = "band", weights = NULL, mask = NULL,
                      lat_bounds = NULL, lon_bounds = NULL, lat_dim = 2, lon_dim = 1, skipna = TRUE) {
  g <- .data_grid(x, lat_dim, lon_dim)
  w <- if (!is.null(weights)) {
    .grid_field(weights, g[1], g[2], "weights", "weights")
  } else {
    if (is.null(lat) || is.null(lon)) stop("lat and lon (or weights) are required")
    weights_2d(lat, lon, lat_bounds, lon_bounds, method, normalize = FALSE)
  }
  .reduce(x, w, lat_dim, lon_dim, mask, skipna, "mean")
}

region_mean <- function(x, lat = NULL, lon = NULL, mask, ...) {
  if (missing(mask) || is.null(mask)) stop("mask is required")
  area_mean(x, lat, lon, mask = mask, ...)
}

area_integral <- function(x, lat = NULL, lon = NULL, ellipsoid = FALSE, mask = NULL, lat_bounds = NULL,
                          lon_bounds = NULL, lat_dim = 2, lon_dim = 1, units = "m2",
                          radius = NULL, weights = NULL, skipna = TRUE) {
  g <- .data_grid(x, lat_dim, lon_dim)
  a <- if (!is.null(weights)) {
    .grid_field(weights, g[1], g[2], "weights", "weights")
  } else {
    if (is.null(lat) || is.null(lon)) stop("lat and lon (or weights) are required")
    cell_area(lat, lon, lat_bounds, lon_bounds, ellipsoid = ellipsoid, radius = radius, units = units)
  }
  .reduce(x, a, lat_dim, lon_dim, mask, skipna, "sum")
}

area_fraction <- function(cond, lat = NULL, lon = NULL, valid = NULL, ...) {
  d <- dim(cond)
  c0 <- as.numeric(cond)
  if (!is.null(valid)) c0[!(as.logical(valid) %in% TRUE)] <- NA_real_
  dim(c0) <- d
  area_mean(c0, lat, lon, ...)
}
