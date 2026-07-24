# FUSE preprocessing scripts

This directory contains the scripts used to prepare meteorological forcing,
streamflow observations, and supporting geospatial data for FUSE from the
public CAMELS-SPAT dataset.

Most users do **not** need to run these scripts. The FUSE repository already
includes the processed input files required to run the example application.
These scripts are provided to ensure complete reproducibility and to
facilitate preparation of FUSE applications for additional CAMELS-SPAT
catchments. :contentReference[oaicite:0]{index=0}

## Workflows

The preprocessing scripts support two independent workflows:

- **Hydrometeorological preprocessing** prepares meteorological forcing and
  streamflow observations for use with FUSE.

- **Geospatial preprocessing** prepares supporting geospatial datasets used
  by distributed FUSE model configurations.

## Hydrometeorological preprocessing

The complete hydrometeorological preprocessing workflow is executed using

```bash
./prepare_fuse_input_data.sh \
    forcing.nc \
    streamflow.nc \
    output.nc
```

This driver script performs the following steps:

1. Standardize the time coordinates in the forcing and streamflow files.
2. Remove unnecessary dimensions from forcing coordinate variables.
3. Compute daily mean air temperature.
4. Merge meteorological forcing and streamflow observations.
5. Compute Oudin potential evapotranspiration.
6. Add basin area and observed runoff depth.
7. Simplify global metadata.
8. Create a legacy-format FUSE input file. :contentReference[oaicite:1]{index=1}

Intermediate files are written to a `work/` directory beneath the output
directory and are retained to facilitate inspection and debugging. :contentReference[oaicite:2]{index=2}

### `prepare_fuse_input_data.sh`

Master preprocessing script that coordinates the complete workflow. It calls
the individual preprocessing scripts in sequence and produces the final FUSE
input file together with a legacy-format version compatible with older FUSE
applications. :contentReference[oaicite:3]{index=3}

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
preserving their values. :contentReference[oaicite:5]{index=5}

### `compute_mean_temperature.sh`

Computes daily mean air temperature from the daily minimum and maximum air
temperature fields and adds the resulting variable to the forcing dataset. :contentReference[oaicite:6]{index=6}

### `merge_netcdf_files.sh`

Creates a single FUSE input file by merging meteorological forcing with
streamflow observations. The forcing record defines the master time axis,
while streamflow observations are inserted over the period where the two
datasets overlap. Missing values are assigned outside the observation
period. :contentReference[oaicite:7]{index=7}

### `compute_oudin_pet.R`

Computes daily potential evapotranspiration using the Oudin method
(implemented in the `airGR` package) from daily mean air temperature and
catchment latitude, and appends the resulting variable to the merged NetCDF
file. :contentReference[oaicite:8]{index=8}

### `common_time_coordinates.sh`

Provides shared routines used to standardize time coordinates across forcing,
streamflow, and legacy FUSE datasets. This script is called internally by
`standardize_time_coordinates.sh` and is not intended to be run directly. :contentReference[oaicite:9]{index=9}

## Geospatial preprocessing

Distributed FUSE model configurations require a description of the
catchment elevation bands. This workflow derives the required information
from a digital elevation model (DEM) and a catchment boundary.

The workflow is executed using

```bash
Rscript make_elev_bands.R \
    dem.tif \
    catchment.shp \
    band_width_m \
    elevation_bands.nc
```

### `make_elev_bands.R`

Creates an elevation-band description for distributed FUSE applications from
a digital elevation model (DEM) and a catchment boundary. The script
partitions the catchment into fixed-width elevation bands and computes the
area fraction and mean elevation of each band. The resulting NetCDF file is
used by FUSE to distribute meteorological forcing and model states across
elevation bands. :contentReference[oaicite:10]{index=10}

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

