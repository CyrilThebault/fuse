###############################################################################
# Create a hydrologic fabric
#
# This script creates a hydrologic fabric NetCDF file from MERIT-Basins river
# and catchment shapefiles. Variable definitions, metadata, and NetCDF names
# are obtained from an external metadata table, allowing the script to remain
# independent of the underlying shapefile attribute names.
#
# Processing steps:
#
#   1. Read river and catchment shapefiles.
#   2. Read metadata describing variables and global attributes.
#   3. Define NetCDF dimensions and variables.
#   4. Create the NetCDF file.
#   5. Write variable data and global metadata.
#   6. Perform basic consistency checks.
#
# Usage:
#
#   Rscript create_hydrofabric.R \
#       catchments.shp \
#       rivers.shp \
#       metadata.csv \
#       hydrofabric.nc
#
# Arguments:
#
#   catchments.shp
#       Catchment polygon shapefile.
#
#   rivers.shp
#       River network shapefile.
#
#   metadata.csv
#       CSV file describing NetCDF variables and global attributes.
#
#   hydrofabric.nc
#       Output hydrologic fabric.
#
# Example:
#
#   catchment_shp=test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_basin.shp
#   river_shp=test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_river.shp
#   metadata=test/metadata/merit_basins_shapefile_metadata.csv
#   hydrofabric_merit=test/CAN_05BB001/distributed/input/hydrofabric_merit.nc
#
#   Rscript test/scripts/distributed/create_hydrofabric.R \
#       "${catchment_shp}" \
#       "${river_shp}" \
#       "${metadata}" \
#       "${hydrofabric_merit}"
#
# Dependencies:
#
#   sf
#       Reading shapefiles.
#
#   ncdf4
#       Creating and writing NetCDF files.
#
###############################################################################

library(sf)
library(ncdf4)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage: Rscript create_hydrofabric.R ",
    "<catchment_shapefile> <river_shapefile> ",
    "<metadata_csv> <output_netcdf>"
  )
}

# Example (Bow River distributed test case):
#
# cd path/to/fuse
#
# Rscript test/scripts/distributed/create_hydrofabric.R \
#   test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_basin.shp \
#   test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_river.shp \
#   test/metadata/merit_basins_shapefile_metadata.csv \
#   test/CAN_05BB001/distributed/input/hydrofabric.nc

catchment_shapefile <- args[1]
river_shapefile     <- args[2]
metadata_file       <- args[3]
output_file         <- args[4]


###############################################################################
# Read shapefiles
###############################################################################

# read shapefiles
catchment <- st_read(catchment_shapefile, quiet = TRUE)
river     <- st_read(river_shapefile, quiet = TRUE)

# ensure both shapefiles contain the same number of features
stopifnot(nrow(river) == nrow(catchment))
n_features <- nrow(river)

###############################################################################
# Read metadata
###############################################################################

# Read metadata
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# split metadata
global_metadata   <- subset(metadata, scope == "global")
variable_metadata <- subset(metadata, scope == "variable")

river_metadata <- subset(
  variable_metadata,
  shapefile %in% c("riv", "all")
)

catchment_metadata <- subset(
  variable_metadata,
  shapefile %in% c("cat", "all")
)

cat("Number of metadata entries :", nrow(metadata), "\n")
cat("River variables            :", nrow(river_metadata), "\n")
cat("Catchment variables        :", nrow(catchment_metadata), "\n")

###############################################################################
# Define NetCDF dimensions
###############################################################################

hru_dim <- ncdim_def(
  name = "hru",
  vals = seq_len(nrow(catchment)),
  units = "", create_dimvar = FALSE
)

seg_dim <- ncdim_def(
  name = "seg",
  vals = seq_len(nrow(river)),
  units = "", create_dimvar = FALSE
)

###############################################################################
# Define NetCDF variables
###############################################################################

var_defs <- list()

for (i in seq_len(nrow(variable_metadata))) {

  # select the NetCDF dimension for this variable
  if (variable_metadata$shapefile[i] == "cat") {
    dim_use <- list(hru_dim)
  } else {
    dim_use <- list(seg_dim)
  }

  # define the variable
  var_defs[[length(var_defs) + 1]] <- ncvar_def(
    name     = variable_metadata$nc_name[i],
    units    = variable_metadata$units[i],
    dim      = dim_use,
    missval  = switch(variable_metadata$type[i], integer = -9999L, double  = -9999.0),
    longname = variable_metadata$long_name[i],
    prec     = variable_metadata$type[i]
  )

}

###############################################################################
# Create and write NetCDF file
###############################################################################

# Set outlet(s): downstream ID not present in the network
river$NextDownID[!(river$NextDownID %in% river$COMID)] <- -999

nc <- nc_create(
  filename = output_file,
  vars = var_defs
)

for (i in seq_len(nrow(variable_metadata))) {

  # select the source shapefile
  if (variable_metadata$shapefile[i] == "cat") {
    values <- catchment[[variable_metadata$name[i]]]
  } else {
    values <- river[[variable_metadata$name[i]]]
  }

  # Wouter's preprocessing for CAMELS-SPAT distinguishes split reaches by 
  # appending decimal suffixes to COMIDs (e.g., 71029071.1, 71029071.2).
  # Convert these to unique integer identifiers for mizuRoute.
  if (variable_metadata$name[i] %in% c("COMID", "hruSegId", "NextDownID")) {
    values <- as.integer(round(values * 10))
  }

  # write the variable
  ncvar_put(
    nc,
    variable_metadata$nc_name[i],
    values
  )

}

# write global metadata
for (i in seq_len(nrow(global_metadata))) {
    ncatt_put(
        nc,
        varid   = 0,
        attname = global_metadata$name[i],
        attval  = global_metadata$value[i]
    )
}

nc_close(nc)

###############################################################################
# Checks
###############################################################################

area_diff <- sum(catchment$unitarea) - max(river$uparea)

if (abs(area_diff) > 1e-6) {
    stop(sprintf(
        "Catchment areas inconsistent with upstream area (difference = %.6f km²)",
        area_diff
    ))
}

cat("\n")
cat("Hydrologic fabric creation completed successfully.\n")
cat("Created:", output_file, "\n")
