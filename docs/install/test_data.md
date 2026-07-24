# Preparing input data for FUSE

FUSE requires two types of input data:

- hydrometeorological forcing and streamflow observations; and
- geospatial data describing elevation bands in each catchment.

The workflows described on this page use the CAMELS-SPAT dataset as an
illustrative example of how FUSE input files can be prepared. The preprocessing
scripts included in the repository are provided primarily to document the
processing steps used to generate the example datasets distributed with FUSE.
They are intended as example implementations rather than general-purpose
preprocessing tools.

Users preparing FUSE applications with other datasets will typically need
to develop preprocessing workflows appropriate for their own data sources.

Further details for the example preprocessing scripts are documented
in `scripts/README.md`.

## Example dataset: CAMELS-SPAT

CAMELS-SPAT provides meteorological forcing, streamflow observations,
catchment attributes, and geospatial data for 1,426 catchments across
Canada and the United States. The dataset and its accompanying
documentation are described in

> Knoben, W. J. M., et al. (2025). *CAMELS-SPAT: a dataset of distributed
> hydro-meteorological variables and geospatial attributes for large-sample
> hydrology in Canada and the United States*. Hydrology and Earth System
> Sciences, 29, 5791–5824.
> https://doi.org/10.5194/hess-29-5791-2025

The complete dataset is publicly available from the Federated Research
Data Repository (FRDR):

https://doi.org/10.20383/103.01306

## Hydrometeorological preprocessing

Meteorological forcing and streamflow observations are prepared using the
master preprocessing script

```bash
./prepare_fuse_input_data.sh \
    forcing.nc \
    streamflow.nc \
    output.nc
```

where

| Argument | Description |
|----------|-------------|
| `forcing.nc` | NetCDF file containing the meteorological forcing variables. |
| `streamflow.nc` | NetCDF file containing the observed streamflow time series for the catchment. |
| `output.nc` | Name of the output NetCDF file that will contain the processed FUSE input data. |


This workflow

- standardizes time-coordinate conventions;
- simplifies forcing coordinate variables;
- computes daily mean air temperature;
- merges meteorological forcing with streamflow observations;
- estimates potential evapotranspiration using the Oudin method;
- adds basin area and runoff depth;
- simplifies metadata; and
- creates a legacy-format FUSE input file.

Intermediate files are written to a `work/` directory beneath the output
directory to facilitate inspection and debugging.

Detailed descriptions of the individual preprocessing scripts are provided
in `scripts/README.md`.

## Geospatial preprocessing

Distributed FUSE model configurations require an elevation-band
description derived from a digital elevation model (DEM) and a catchment
boundary. This preprocessing step is performed using

```bash
Rscript make_elev_bands.R \
    dem.tif \
    catchment.shp \
    band_width_m \
    elevation_bands.nc
```

where

| Argument | Description |
|----------|-------------|
| `dem.tif` | Digital elevation model (DEM) covering the catchment. |
| `catchment.shp` | Catchment boundary polygon (ESRI Shapefile). |
| `band_width_m` | Width of each elevation band, in metres. |
| `elevation_bands.nc` | Output NetCDF file containing the elevation-band description for FUSE. |

Additional information on this workflow is provided in
`scripts/README.md`.
