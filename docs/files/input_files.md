# Input data

The input data provide the hydrometeorological time series and ancillary
datasets required to run a FUSE simulation. These files are stored in
the input directory specified by the `input_dir` entry in the control
file.

For the Bow River test case, the control file specifies

```toml
[filepaths]
input_dir = "test/CAN_05BB001/input/"
```

The input directory contains the following files.

## Hydrometeorological time series

The hydrometeorological input file provides the time-varying data
required to drive and evaluate the hydrologic model. The file is
identified by appending the `forcing_suffix` specified in the control
file to the basin identifier.

For the Bow River test case, the control file specifies

```toml
[input]
forcing_suffix = "_daymet_qobs_merged_legacy.nc"
```

which identifies the input file

```text
test/CAN_05BB001/input/CAN_05BB001_daymet_qobs_merged_legacy.nc
```

The input file is stored in NetCDF format and contains the
hydrometeorological time series required by FUSE, including
meteorological forcing variables, streamflow observations, and derived
variables such as potential evapotranspiration. These datasets are
prepared from the CAMELS-SPAT dataset using the preprocessing scripts
distributed with the FUSE repository.

The hydrometeorological variables use the same NetCDF structure for
both catchment and gridded simulations. FUSE determines the spatial
configuration from the dimensions of the latitude and longitude
coordinates. A lumped catchment simulation is performed when both
dimensions have length one. When either spatial dimension contains
multiple values, FUSE applies the model over the corresponding grid.

The names of the coordinate and data variables are specified in the
`forcing_coords` and `forcing_vars` sections of the control file. For
the Bow River test case, these entries are

```toml
[forcing_coords]
time      = "time"
latitude  = "latitude"
longitude = "longitude"

[forcing_vars]
precip = "prcp"
temp   = "temp"
pet    = "pet_oudin"
qobs   = "runoff_obs"
```

## Elevation bands

The elevation-band file defines the elevation bands used by the snow
module. The file is identified by appending the `elevbands_suffix`
specified in the control file to the basin identifier.

For the Bow River test case, the control file specifies

```toml
[input]
elevbands_suffix = "_elev_bands.nc"
```

which identifies the elevation-band file

```text
test/CAN_05BB001/input/CAN_05BB001_elev_bands.nc
```

The elevation-band file is stored in NetCDF format and is required when
the snow module is active. It defines the elevation-band
characteristics used to distribute meteorological forcing and model
states within each spatial unit. The elevation bands are derived from a
digital elevation model using the preprocessing scripts distributed
with the FUSE repository.

The spatial dimensions of the elevation-band file must be consistent
with those of the hydrometeorological input file.
