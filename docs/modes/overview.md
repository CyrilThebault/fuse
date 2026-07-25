# Overview

FUSE (Framework for Understanding Structural Errors) is a flexible
hydrologic modelling framework for constructing, calibrating, and
evaluating conceptual rainfall–runoff models. Rather than providing a
single fixed model, FUSE allows users to construct alternative model
instantiations within a common software framework.

A central concept in FUSE is that constructing a hydrologic model
requires a series of subjective modelling decisions. These decisions
include, for example, how to represent soil moisture storage, runoff
generation, evaporation, routing, and snow processes. Rather than
embedding these choices within a fixed model structure, FUSE represents
each decision explicitly, allowing alternative model structures to be
assembled by selecting among alternative process representations.

The model structure is therefore defined independently of the model
parameters and numerical solution methods. This separation makes it
possible to investigate the influence of structural, parametric, and
numerical choices on model behaviour within a consistent modelling
framework.

A FUSE simulation is configured by defining five complementary aspects
of the simulation:

**[Simulation setup](simulation_setup.md).** Define the simulation
period, evaluation period, calibration options, objective function, the
spatial configuration, and the locations of the input, output, and
settings files. These settings are specified in the control file. The
Bow River test case uses

```text
test/CAN_05BB001/settings/fuse_v2/fuse_control_CAN_05BB001.toml
```

**[Spatial configuration](spatial_configuration.md).** Define whether
the hydrologic model is applied to a single catchment or a gridded
domain, together with the use of elevation bands to improve the
representation of snow accumulation and melt. These settings are
specified through entries in the toml control file.

**[Hydrologic model definition](hydrologic_model.md).** Define the
hydrologic model by specifying the model structure, parameter
definitions, and numerical solution methods. These settings are
specified through three settings files identified by entries in the
control file.

**[Input data files](../files/input_files.md).** Provide the meteorological forcing, catchment
characteristics, observations, and other information required for a
simulation. The Bow River test case stores these files in

```text
test/CAN_05BB001/input/
```

**[Output files](../files/output_files.md).** FUSE simulations of streamflow, model
state variables, fluxes, and diagnostic information are written during
model execution. The output directory is specified by the `output_dir`
entry in the control file.

The remainder of this User Guide describes each of these aspects of the
simulation in detail.

## Running FUSE

FUSE is controlled through a command-line interface that supports
multiple run modes for model evaluation and parameter calibration.

For details of the command-line interface, required arguments, and
example simulations, see
**[Running simulations](../install/test_cases.md)**.
