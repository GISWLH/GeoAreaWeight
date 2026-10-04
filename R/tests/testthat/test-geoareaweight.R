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
