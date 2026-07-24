# Running simulations

The FUSE repository includes a complete example application for the Bow
River basin (station CAN_05BB001), including all processed input files
required to run simulations immediately. The following sections describe
the FUSE command-line interface and illustrate the most common execution
modes.

Users interested in preparing new applications from the public
CAMELS-SPAT dataset should see **Preparing input data from CAMELS-SPAT**.

## Command-line interface

FUSE is controlled from the command line. To view the available command-line
arguments and execution modes, run

```bash
bin/fuse.exe --help
```

(or equivalently

```bash
bin/fuse.exe -h
```

).

The command-line help provides a concise summary of the required arguments,
available execution modes, and optional runtime settings. Rather than
duplicating the complete help message here, the sections below describe the
most commonly used execution modes and illustrate their use.

## Required arguments

All FUSE simulations require three command-line arguments:

- `--domid` specifies the domain (or catchment) identifier.
- `--control` specifies the TOML control file.
- `--runmode` specifies how FUSE should execute.

Since all of the examples below use the same control file, it is
convenient to define a shell variable:

```bash
control_file=test/CAN_05BB001/settings/fuse_control_CAN_05BB001.toml
```

For the Bow River example included with FUSE, the required arguments are

```bash
-d CAN_05BB001 \
-c $control_file 
-m <run mode>
```

In the examples below, the domain identifier and control file remain the same; only the run mode and optional arguments change.

The following sections describe the different run modes.

## Run with the default parameter values (`def`)

This mode runs FUSE using the default parameter values specified in the
parameter definition file.

```bash
bin/fuse.exe \
    -d CAN_05BB001 \
    -c $control_file 
    -m def
```

## Run a specified parameter set (`idx`)

This mode runs FUSE using a specified parameter set stored in a NetCDF
parameter file. The `--index` argument identifies the parameter set to use.

```bash
bin/fuse.exe \
    -d CAN_05BB001 \
    -c $control_file 
    -m idx \
    -s calibration.nc \
    -i 17
```

## Run the best parameter set (`opt`)

This mode runs FUSE using the optimal parameter set contained in a NetCDF
parameter file.

```bash
bin/fuse.exe \
    -d CAN_05BB001 \
    -c $control_file 
    -m opt \
    -s calibration.nc
```

## Calibrate model parameters (`sce`)

This mode estimates model parameters using the Shuffled Complex Evolution
(SCE) optimization algorithm.

```bash
bin/fuse.exe \
    -d CAN_05BB001 \
    -c $control_file 
    -m sce
```

## Optional runtime arguments

Several optional command-line arguments modify model execution without changing
the control file.

- `--param` overrides one or more parameter values.
- `--restart` specifies how frequently restart files are written.
- `--tag` appends a string to the output filenames.

For example, the following command overrides two parameter values while using
the default execution mode.

```bash
bin/fuse.exe \
    -d CAN_05BB001 \
    -c $control_file 
    -m def \
    -p MAXWATR_1=1000 \
    -p PERCRTE=0.25
```

For a complete description of all command-line arguments and options, consult
the built-in help:

```bash
bin/fuse.exe --help
```
