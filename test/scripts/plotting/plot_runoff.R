library(ncdf4)

quartz(width = 10, height = 6)

# ----- read in mizuRoute simulations -----

mizu_path <- 'test/CAN_05BB001/distributed/work/'
mizu_name <- 'CAN_05BB001_lumped_route.h.1989-01-01-43200.nc'

mizu_file <- file.path(mizu_path, mizu_name)

nc <- nc_open(mizu_file)

  segid_mizu    <- ncvar_get(nc, "reachID")
  bounds_mizu   <- ncvar_get(nc, "time_bounds")
  qrunoff_mizu  <- ncvar_get(nc, "basRunoff")
  qrouted_mizu  <- ncvar_get(nc, "KWroutedRunoff")

nc_close(nc)

# ----- read in mizuRoute simulations -----

fuse_path <- 'test/CAN_05BB001/lumped/work/'
fuse_name <- 'CAN_05BB001_2__runs_opt.nc'

fuse_file <- file.path(fuse_path, fuse_name)

nc <- nc_open(fuse_file)

  segid_fuse    <- ncvar_get(nc, "seg")
  bounds_fuse   <- ncvar_get(nc, "time_bnds")
  qrunoff_fuse  <- ncvar_get(nc, "q_routed")
  qrouted_fuse  <- ncvar_get(nc, "q_reach")
  qobs_fuse     <- ncvar_get(nc, "q_obs")

nc_close(nc)

# ----- read in observations -----

obs_path <- 'test/CAN_05BB001/lumped/input/'
obs_name <- 'CAN_05BB001_daymet_qobs_merged.nc'

obs_file <- file.path(obs_path, obs_name)

nc <- nc_open(obs_file)

  time_obs   <- ncvar_get(nc, "time")
  runoff_obs <- ncvar_get(nc, "runoff_obs")
  basin_area <- ncvar_get(nc, "basin_area")

  fill <- ncatt_get(nc, "runoff_obs", "_FillValue")$value

nc_close(nc)

runoff_obs[runoff_obs == fill] <- NA

# ----- extract information -----

time_mizu <- colMeans(bounds_mizu)
time_fuse <- colMeans(bounds_fuse)

date_mizu <- as.POSIXct("1950-01-01", tz = "UTC") + time_mizu * 86400
date_fuse <- as.POSIXct("1950-01-01", tz = "UTC") + time_fuse * 86400
date_obs  <- as.POSIXct("1950-01-01", tz = "UTC") + time_obs * 86400

# choose reach
iRch_mizu <- which(segid_mizu == 710285850)
iRch_fuse <- which(segid_fuse == 710285850)

# convert routed streamflow from m3/s to mm/day
routed_mizu <- qrouted_mizu[iRch_mizu, ] * 86400 * 1000 / (basin_area * 1e6)
routed_fuse <- qrouted_fuse[iRch_fuse, ] * 86400 * 1000 / (basin_area * 1e6)

# extract basin runoff (spatially constant)
runoff_mizu <- qrunoff_mizu[iRch_mizu, ]
runoff_fuse <- qrunoff_fuse # lumped -- singleton dimensions removed

# ----- plot results -----

# define xlim and ylim
xlim <- as.POSIXct(c("1990-05-01", "1990-07-01"), tz = "UTC")
ylim <- c(0,10)

plot(date_obs, runoff_obs, type = "n",
     xlab = "Date", ylab = "Runoff (mm/day)",
     xlim = xlim,
     ylim = ylim,
     main = paste("Reach", segid_mizu[iRch_mizu]))

lines(date_obs,  runoff_obs,  col = "black")
lines(date_mizu, runoff_mizu, col = "lightblue")
lines(date_mizu, routed_mizu, col = "darkblue")

lines(date_fuse, runoff_fuse, col = "magenta", lty=2)
lines(date_fuse, routed_fuse, col = "magenta")


points(date_mizu, runoff_mizu, col = "lightblue", pch = 16, cex = 0.5)
points(date_mizu, routed_mizu, col = "darkblue",  pch = 16, cex = 0.5)

legend("topleft",
       legend = c("Observed", "Runoff", "Routed"),
       col = c("black", "lightblue", "darkblue"),
       lty = 1)
