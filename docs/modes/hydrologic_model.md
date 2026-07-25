# Hydrologic model definition

The hydrologic model used in a FUSE simulation is defined through three
complementary components:

- the **model structure**, which specifies the process
  representations included in the model;
- the **model parameters**, which specify the parameter values,
  allowable ranges, and calibration flags; and
- the **numerical solution methods**, which specify how the governing
  equations are integrated.

These components are specified through three settings files.

The settings files reside in the directory specified by the
`settings_dir` entry in the control file. Individual settings files are
identified through entries in the `model` section of the control file.

Example settings files are provided with the Bow River test case included
in the FUSE repository. Specifically, the settings directory

```text
test/CAN_05BB001/settings/fuse_v2/
```

contains the following files:

- **[Model structure](structure.md)** (`fuse_zDecisions_2.txt`,
  referenced by `decisions_file`) defines the model structure by
  selecting one process representation for each modelling decision. The
  available modelling decisions are described by
  [Clark et al. (2008)](https://doi.org/10.1029/2007WR006735), with the
  snow model options (Decision 9) described by
  [Henn et al. (2015)](https://doi.org/10.1002/2014WR016736).

- **[Model parameters](parameters.md)** (`fuse_zConstraints_snow.txt`,
  referenced by `constraints_file`) defines the model parameters,
  including their default values, lower and upper bounds, and
  calibration flags. The parameter definitions correspond to those
  described by
  [Clark et al. (2008)](https://doi.org/10.1029/2007WR006735) and
  [Henn et al. (2015)](https://doi.org/10.1002/2014WR016736).

- **[Numerical solution methods](numerical_methods.md)**
  (`fuse_zNumerix.txt`, referenced by `numerics_file`) specifies the
  numerical methods and solver options used to integrate the model
  equations. The influence of these numerical methods on model
  behaviour is discussed by
  [Clark and Kavetski (2010)](https://doi.org/10.1029/2009WR008894) and
  [Kavetski and Clark (2010)](https://doi.org/10.1029/2009WR008896).

The following pages describe each component of the hydrologic model in
more detail.
