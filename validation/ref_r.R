suppressMessages(library(geoareaweight)); args <- commandArgs(TRUE)
out <- character(0)
put <- function(k, a) { a <- as.numeric(a); out <<- c(out, sprintf("%s[%d],%.15e", k, seq_along(a) - 1, a)) }
lat <- seq(-90, 90, by = 2.5); lon <- seq(0, 357.5, by = 2.5); nt <- 5
nl <- length(lon); nla <- length(lat)
x <- array(NA_real_, c(nl, nla, nt))
for (i in 1:nl) for (j in 1:nla) for (t in 1:nt) x[i, j, t] <- if ((i * 7 + j * 3 + t) %% 11 == 0) NA else sin(0.3 * i + 0.17 * j + 1.1 * t) + 0.01 * j
put("w_cos", area_weights(lat, method = "cos")); put("w_band", area_weights(lat)); put("w_ell", area_weights(lat, method = "ellipsoid"))
put("area_sum_sphere", sum(cell_area(lat, lon))); put("area_sum_ell", sum(cell_area(lat, lon, ellipsoid = TRUE)))
for (m in c("cos", "band", "ellipsoid", "none")) put(paste0("mean_", m), area_mean(x, lat, lon, method = m))
mk <- matrix(FALSE, nla, nl); mk[21:50, 11:90] <- TRUE
put("mean_mask", area_mean(x, lat, lon, mask = mk)); put("mean_skipna_false", area_mean(x, lat, lon, skipna = FALSE))
x0 <- x[, , 1]; x0[is.na(x0)] <- 0
put("integral", area_integral(x0, lat, lon) / 1e12)
lat2 <- c(85, 60, 33, 10, -12, -41, -70, -88); lon2 <- c(0, 7, 30, 31, 100, 200, 300)
x2 <- outer(seq_along(lon2), seq_along(lat2), function(i, j) cos(0.5 * i + 0.9 * j) + 0.1 * j)   # lon x lat
put("w2d_irreg", t(weights_2d(lat2, lon2)))   # python order is (lat, lon) row-major
put("mean_irreg_band", area_mean(x2, lat2, lon2)); put("mean_irreg_cos", area_mean(x2, lat2, lon2, method = "cos"))
put("mean_irreg_ell", area_mean(x2, lat2, lon2, method = "ellipsoid"))
cond <- ifelse(is.na(x[, , 1]), NA_real_, as.numeric(x[, , 1] > 0.5))
put("frac", area_fraction(cond, lat, lon))
writeLines(out, args[1])
