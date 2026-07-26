# Simulation setup

The simulation setup defines how the hydrologic model is applied and
evaluated. This includes the simulation period, the evaluation period,
the objective function used for model assessment, calibration options,
the spatial configuration, and the locations of the input, output, and
settings files.

These settings are specified in a TOML control file, which serves as the
primary entry point for a FUSE simulation. An example control file is
provided with the Bow River test case:

```text
test/CAN_05BB001/lumped/settings/fuse_v2/fuse_control_CAN_05BB001.toml
```

The control file is organized into several sections, each defining a
different aspect of the simulation setup.

| Section | Purpose |
|---------|---------|
| `filepaths` | Defines the locations of the input, output, and settings directories. |
| `input` | Specifies the forcing and elevation-band input files. |
| `model` | Identifies the settings files that define the model structure, parameter definitions, and numerical solution methods. |
| `spatial` | Defines the spatial configuration, including catchment mode or grid mode and the use of elevation bands. |
| `forcing_coords` | Specifies the names of the coordinate variables in the forcing dataset. |
| `forcing_vars` | Specifies the names of the meteorological forcing and observed streamflow variables in the forcing dataset. |
| `output` | Defines output options, including the model identifier and variables to write. |
| `run_periods` | Specifies the simulation and evaluation periods. |
| `calibration` | Defines the objective function and any transformations used during parameter estimation. |
| `sce` | Specifies the control parameters for the Shuffled Complex Evolution (SCE) optimization algorithm. |
