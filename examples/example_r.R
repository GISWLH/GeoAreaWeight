library(geoareaweight); library(ncdf4)
nc  <- nc_open("air.mon.mean.nc")
lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
x   <- ncvar_get(nc, "air")                    # lon x lat x time
gm  <- area_mean(x, lat, lon)                  # area-weighted global mean (time series)
naive <- apply(x, 3, mean)                     # the mistake
cat(sprintf("mean: unweighted %.2f  weighted %.2f\n", mean(naive), mean(gm)))
