# Model parameters

The model parameters define the values of the coefficients used in the
hydrologic model equations. For each parameter, FUSE specifies a default
value together with the allowable range of parameter values and whether
the parameter may be adjusted during calibration.

The model parameters are specified through the parameter constraints
file, which is identified by the `constraints_file` entry in the control
file.

The parameter definitions correspond to those described by
[Clark et al. (2008)](https://doi.org/10.1029/2007WR006735) and
[Henn et al. (2015)](https://doi.org/10.1002/2014WR016736).

The Bow River test case uses the parameter constraints file

```text
test/CAN_05BB001/settings/fuse_v2/fuse_zConstraints_snow.txt
```

Each entry in the parameter constraints file specifies

- the parameter name;
- the default parameter value;
- the lower and upper bounds; and
- whether the parameter is included in model calibration.

The default parameter values are used when running FUSE in its default
mode. During calibration, parameter values are constrained to remain
within the specified lower and upper bounds.

The calibration procedure itself is configured separately through the
control file, which specifies the calibration and evaluation periods,
the objective function, and the optimization settings. See
**[Simulation setup](simulation_setup.md)** for details.
