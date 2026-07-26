#!/usr/bin/env Rscript

# ============================================================
# compute_oudin_pet.R
#
# Compute daily Oudin potential evapotranspiration from daily
# mean air temperature using airGR and store the results in the
# merged FUSE NetCDF file.
#
# Usage:
#
#   Rscript compute_oudin_pet.R merged_file.nc
#
# Required packages:
#
#   airGR
#   ncdf4
# ============================================================

library(airGR)
library(ncdf4)

# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------

debug <- FALSE

if (debug) {

    nc_file <- "../CAN_05BB001/lumped/input/CAN_05BB001_daymet_qobs_merged.nc"

} else {

    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) != 1L) {
        stop(
            "Usage: Rscript compare_oudin_pet.R <merged_netcdf_file>",
            call. = FALSE
        )
    }

    nc_file <- args[1]

}

if (!file.exists(nc_file)) {
    stop("NetCDF file not found: ", nc_file, call. = FALSE)
}

# ------------------------------------------------------------
# Read NetCDF file
# ------------------------------------------------------------

# Open NetCDF file
nc <- ncdf4::nc_open(nc_file)

# Read time information
time <-  ncdf4::ncvar_get(nc, "time")
units <- ncdf4::ncatt_get(nc, "time", "units")$value

# Extract the reference date from:
# "days since 1950-01-01 00:00:00"
origin <- as.Date(sub(".*since\\s+([0-9]{4}-[0-9]{2}-[0-9]{2}).*", "\\1", units))

# compute the julian day from the time coordinate
date <- origin + floor(time)
jday <- as.integer(format(date, "%j"))

# Read temperature
temperature <- ncdf4::ncvar_get(nc, "temp", collapse_degen = TRUE)
temp_fill   <- ncdf4::ncatt_get(nc, "temp", "_FillValue")$value

# Read latitude
latitude <- ncdf4::ncvar_get(nc, "latitude", collapse_degen = TRUE)

# Close NetCDF file
ncdf4::nc_close(nc)

# ------------------------------------------------------------
# Compute Oudin PET
# ------------------------------------------------------------

pet_airgr <- rep(NA_real_, length(temperature))

temperature[temperature == temp_fill] <- NA
valid_temperature <- !is.na(temperature)

pet_airgr[valid_temperature] <- airGR::PE_Oudin(
    JD = jday[valid_temperature],
    Temp = as.numeric(temperature[valid_temperature]),
    Lat = as.numeric(latitude),
    LatUnit = "deg",
    TimeStepIn = "daily",
    TimeStepOut = "daily",
    RunFortran = FALSE
)

# ------------------------------------------------------------
# Add the Oudin PET variable to the NetCDF file
# ------------------------------------------------------------

nc <- ncdf4::nc_open(nc_file, write = TRUE)

if (!"pet_oudin" %in% names(nc$var)) {

    pet_oudin_def <- ncdf4::ncvar_def(
        name = "pet_oudin",
        units = "mm/day",
        dim = list(
            nc$dim$hru,
            nc$dim$time
        ),
        missval = -9999.0,
        longname = paste(
            "potential evapotranspiration estimated",
            "with the Oudin method"
        ),
        prec = "float"
    )

    nc <- ncdf4::ncvar_add(
        nc,
        pet_oudin_def
    )
}

ncdf4::ncvar_put(
    nc,
    varid = "pet_oudin",
    vals = matrix(pet_airgr, nrow = 1)
)

ncdf4::nc_close(nc)
