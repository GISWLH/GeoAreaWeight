grid <- function(res, pole = FALSE) {
  lat <- if (pole) seq(-90, 90, by = res) else seq(-90 + res / 2, 90 - res / 2, by = res)
  lon <- seq(0, 360 - res, by = res) + if (pole) 0 else res / 2
  list(lat = lat, lon = lon)
}
R <- 6371008.8
test_that("areas sum to the sphere / WGS84 ellipsoid", {
  for (cfg in list(c(1, 0), c(2.5, 1), c(5, 0))) {
    g <- grid(cfg[1], cfg[2] == 1)
    expect_equal(sum(cell_area(g$lat, g$lon)), 4 * pi * R^2, tolerance = 1e-12)
  }
  g <- grid(1)
  expect_equal(sum(cell_area(g$lat, g$lon, ellipsoid = TRUE)) / 1e6, 510065621.7, tolerance = 1e-7)
  expect_equal(authalic_radius(), 6371007.1809, tolerance = 1e-9)
})
test_that("descending/irregular latitude", {
  lat <- c(80, 60, 30, 10, -5, -40, -75)
  expect_equal(sum(band_weights(lat)), 1)
  expect_equal(band_weights(lat), rev(band_weights(rev(lat))))
  expect_equal(2 * pi * sum(cell_area(lat)), 4 * pi * R^2, tolerance = 1e-12)
})
test_that("constant field and analytic |lat| mean", {
  g <- grid(1); nl <- length(g$lon); nt <- length(g$lat)
  x <- matrix(3.7, nl, nt)
  for (m in c("cos", "band", "ellipsoid", "none")) expect_equal(area_mean(x, g$lat, g$lon, method = m), 3.7)
  y <- matrix(rep(abs(g$lat), each = nl), nl, nt)
  expect_equal(area_mean(y, g$lat, g$lon), (pi / 2 - 1) * 180 / pi, tolerance = 1e-4)
  expect_equal(area_mean(y, g$lat, g$lon, method = "none"), 45)
})
test_that("NaN renormalisation, mask, time dim, fraction", {
  g <- grid(5); nl <- length(g$lon); nt <- length(g$lat)
  set.seed(1); x <- array(rnorm(nl * nt * 6), c(nl, nt, 6))
  x[, 1:5, ] <- NA; x[10:30, 10:20, 4] <- NA
  m <- area_mean(x, g$lat, g$lon)
  expect_length(m, 6)
  w <- t(weights_2d(g$lat, g$lon, normalize = FALSE))    # lon x lat
  for (t in 1:6) { v <- is.finite(x[, , t]); expect_equal(m[t], sum(x[, , t][v] * w[v]) / sum(w[v])) }
  expect_true(all(is.na(area_mean(x, g$lat, g$lon, skipna = FALSE))))
  mk <- matrix(FALSE, nt, nl); mk[20:30, 5:40] <- TRUE
  mm <- area_mean(x, g$lat, g$lon, mask = mk)
  wm <- w * t(mk); v <- is.finite(x[, , 1])
  expect_equal(mm[1], sum(x[, , 1][v] * wm[v]) / sum(wm[v]))
  g1 <- grid(1); cond <- matrix(rep(g1$lat >= 30, each = length(g1$lon)), length(g1$lon))
  expect_equal(area_fraction(cond, g1$lat, g1$lon), 0.25, tolerance = 1e-12)
  expect_equal(area_integral(matrix(1, length(g1$lon), length(g1$lat)), g1$lat, g1$lon), 4 * pi * R^2, tolerance = 1e-12)
})
test_that("integral: NaN skipped, NA when no valid cell, skipna = FALSE propagates", {
  g <- grid(10); nl <- length(g$lon); nt <- length(g$lat)
  x <- array(1, c(nl, nt, 2)); x[, , 2] <- NA; x[1, 1, 1] <- NA
  a <- cell_area(g$lat, g$lon)                 # nlat x nlon
  r <- area_integral(x, g$lat, g$lon)
  expect_equal(r[1], sum(a) - a[1, 1])
  expect_true(is.na(r[2]))
  expect_true(is.na(area_integral(x, g$lat, g$lon, skipna = FALSE)[1]))
  expect_equal(area_integral(x[, , 1], weights = a), r[1])
})
test_that("method validation, 2-D lat fallback warns", {
  g <- grid(10)
  expect_error(area_weights(g$lat, method = "nonsense"), "unknown method")
  lat2 <- outer(g$lat, rep(1, length(g$lon)))
  expect_error(weights_2d(lat2, method = "nonsense"), "unknown method")
  expect_warning(w <- weights_2d(lat2, normalize = FALSE), "cos")
  expect_equal(w, cos(lat2 * pi / 180))
  expect_equal(area_weights(g$lat, method = "WGS84"), area_weights(g$lat, method = "ellipsoid"))
})
test_that("input validation", {
  expect_error(area_weights(seq(0, 350, 10)), "swap lat and lon")
  expect_error(band_weights(c(0, 10, 5)), "strictly")
  expect_error(cell_area(c(0, 10), 5), "single longitude")
})
test_that("masks: NA/0 exclude, lon x lat orientation accepted, weights validated", {
  g <- grid(10); nl <- length(g$lon); nt <- length(g$lat)
  set.seed(2); x <- matrix(rnorm(nl * nt), nl, nt)
  sel <- outer(g$lat > 30, rep(TRUE, nl))      # nlat x nlon
  ref <- area_mean(x, g$lat, g$lon, mask = sel)
  expect_equal(area_mean(x, g$lat, g$lon, mask = ifelse(sel, 1, NA)), ref)
  expect_equal(area_mean(x, g$lat, g$lon, mask = t(sel)), ref)
  a <- cell_area(g$lat, g$lon)
  expect_equal(area_mean(x, weights = t(a)), area_mean(x, g$lat, g$lon))
  expect_error(area_mean(x, weights = -a), "non-negative")
})
test_that("longitudes crossing the wrap point", {
  lat <- c(-45, 45)
  lw <- c(seq(182.5, 357.5, 5), seq(2.5, 177.5, 5))
  a <- cell_area(lat, lw)
  expect_equal(a, cell_area(lat, seq(2.5, 357.5, 5)))
})
test_that("skipna = FALSE on complete data equals the default", {
  g <- grid(10); nl <- length(g$lon); nt <- length(g$lat)
  set.seed(3); x <- array(rnorm(nl * nt * 3), c(nl, nt, 3))
  expect_equal(area_mean(x, g$lat, g$lon, skipna = FALSE), area_mean(x, g$lat, g$lon))
  x[1, 1, 2] <- NA
  r <- area_mean(x, g$lat, g$lon, skipna = FALSE)
  expect_true(is.na(r[2]) && all(is.finite(r[c(1, 3)])))
})
