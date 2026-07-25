# FUSE preprocessing scripts

This directory contains the scripts used to prepare meteorological forcing,
streamflow observations, and supporting geospatial data for FUSE from the
public CAMELS-SPAT dataset.

Most users do **not** need to run these scripts. The FUSE repository already
includes the processed input files required to run the example application.
These scripts are provided to ensure complete reproducibility and to
facilitate preparation of FUSE applications for additional CAMELS-SPAT
catchments.

## Workflows

The preprocessing scripts support two independent workflows:

- **Hydrometeorological preprocessing** prepares meteorological forcing and
  streamflow observations for use with FUSE.

- **Geospatial preprocessing** prepares supporting geospatial datasets used
  by distributed FUSE model configurations.

## Hydrometeorological preprocessing

Meteorological forcing and streamflow observations are prepared using the
master preprocessing script

```bash
./prepare_fuse_input_data.sh \
    $forcing_file \
    $streamflow_file \
    $output_file
```

where

| Argument | Description |
|----------|-------------|
| `$forcing_file`    | NetCDF file containing the meteorological forcing variables. |
| `$streamflow_file` | NetCDF file containing the observed streamflow time series for the catchment. |
| `$output_file`     | Name of the output NetCDF file that will contain the processed FUSE input data. |

This workflow

- standardizes time-coordinate conventions;
- simplifies forcing coordinate variables;
- computes daily mean air temperature;
- merges meteorological forcing with streamflow observations;
- estimates potential evapotranspiration using the Oudin method;
- adds basin area and runoff depth;
- simplifies metadata; and
- creates a legacy-format FUSE input file.

For the Bow River test case, these shell variables can be defined as
```bash
forcing_file=test/CAN_05BB001/input/forcing/CAN_05BB001_daymet_lumped.nc
streamflow_file=test/CAN_05BB001/input/q_obs/CAN_05BB001_daily_flow_observations.nc
output_file=test/CAN_05BB001/input/CAN_05BB001_daymet_qobs_merged.nc
```

Intermediate files are written to a `work/` directory beneath the output
directory and are retained to facilitate inspection and debugging.

### `prepare_fuse_input_data.sh`

Master preprocessing script that coordinates the complete workflow. It calls
the individual preprocessing scripts in sequence and produces the final FUSE
input file together with a legacy-format version compatible with older FUSE
applications.

### `standardize_time_coordinates.sh`

Standardizes the time coordinate conventions used in CAMELS-SPAT forcing,
streamflow observation, and legacy FUSE input files. The script converts all
datasets to a common representation using CF-compliant time coordinates with
units of `"days since 1950-01-01"` and reconstructs the time coordinate from
the interval bounds. Common processing routines are implemented in
`common_time_coordinates.sh`. 

### `fix_coords.sh`

Simplifies CAMELS-SPAT forcing files by removing the unnecessary time
dimension from the latitude, longitude, and HRU identifier variables while
preserving their values.

### `compute_mean_temperature.sh`

Computes daily mean air temperature from the daily minimum and maximum air
temperature fields and adds the resulting variable to the forcing dataset.

### `merge_netcdf_files.sh`

Creates a single FUSE input file by merging meteorological forcing with
streamflow observations. The forcing record defines the master time axis,
while streamflow observations are inserted over the period where the two
datasets overlap. Missing values are assigned outside the observation
period.

### `compute_oudin_pet.R`

Computes daily potential evapotranspiration using the Oudin method
(implemented in the `airGR` package) from daily mean air temperature and
catchment latitude, and appends the resulting variable to the merged NetCDF
file.

### `common_time_coordinates.sh`

Provides shared routines used to standardize time coordinates across forcing,
streamflow, and legacy FUSE datasets. This script is called internally by
`standardize_time_coordinates.sh` and is not intended to be run directly.

## Geospatial preprocessing

Distributed FUSE model configurations require an elevation-band
description derived from a digital elevation model (DEM) and a catchment
boundary.


### `make_elev_bands.R`

Creates an elevation-band description for distributed FUSE applications from
a digital elevation model (DEM) and a catchment boundary. The R script
partitions the catchment into fixed-width elevation bands and computes the
area fraction and mean elevation of each band. The resulting NetCDF file is
used by FUSE to distribute meteorological forcing and model states across
elevation bands.

This preprocessing step is performed using

```bash
Rscript make_elev_bands.R \
    $dem_file \
    $catchment_file \
    $band_width_m \
    $output_file
```

where

| Argument | Description |
|----------|-------------|
| `$dem_file` | Digital elevation model (DEM) covering the catchment. |
| `$catchment_file` | Catchment boundary polygon (ESRI Shapefile). |
| `$band_width_m` | Width of each elevation band, in metres. |
| `$output_file` | Output NetCDF file containing the elevation-band description for FUSE. |

This workflow

- clips the DEM to the catchment boundary;
- partitions the catchment into elevation bands;
- computes the area and mean elevation of each band; and
- writes the elevation-band description to a NetCDF file compatible with FUSE.

For the Bow River test case, these shell variables can be defined as

```bash
dem_file=test/CAN_05BB001/input/geospatial/tif/CAN_05BB001_merit_hydro_elv.tif
catchment_file=test/CAN_05BB001/input/geospatial/shp/CAN_05BB001_lumped.shp
band_width_m=100
output_file=test/CAN_05BB001/input/CAN_05BB001_elev_bands.nc
```

## Software requirements

The shell scripts require the NetCDF Operators (NCO) which can be installed
using homebrew

```bash
brew install nco
```

The R scripts require the following packages:

- `airGR`
- `ncdf4`
- `terra`

