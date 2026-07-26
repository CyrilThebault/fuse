# Output data

FUSE writes model results to the output directory specified by the
`output_dir` entry in the control file.

For the Bow River test case, the control file specifies

```toml
[filepaths]
output_dir = "test/CAN_05BB001/lumped/output/"
```

The output directory contains the following files.

- **Model output file** (`*_runs.nc`) — simulated hydrologic time series
  together with coordinate variables and metadata.

- **Parameter file** (`*_para.nc`) — parameter values used during the
  simulation or generated during calibration.

Both the model output and parameter files include software provenance
metadata, including the FUSE version, build time, Git branch, and Git
commit hash, to support reproducibility.

## Model output time series

The model output file contains the simulated time series produced during
a FUSE run. The filename depends on the simulation mode (evaluation,
calibration, sensitivity analysis, etc.), with the final component of
the filename identifying the FUSE mode used to generate the output.

The output file is stored in NetCDF format and contains the simulated
model variables together with the coordinate variables, model metadata,
and software version information.

The output file always contains

- simulation time;
- latitude and longitude coordinates;
- the simulated streamflow variables (`q_instnt` and
  `q_routed`); and

When elevation bands are enabled, variables associated with snow
processes include an additional elevation-band dimension.

## Parameter values

The parameter file contains the parameter values used during the model
run. Depending on the simulation mode, this may represent a single
parameter set or multiple parameter sets generated during calibration.

The parameter file is stored in NetCDF format. Each parameter is written
as a separate variable together with its descriptive metadata,
including the parameter name, description, units, and missing-value
definition.

Parameters associated with elevation bands are stored with an
additional elevation-band dimension, while scalar parameters are stored
as one-dimensional arrays indexed by parameter set.
