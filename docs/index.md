## FUSE in a nutshell

The Framework for Understanding Structural Errors or [FUSE](https://github.com/CH-Earth/fuse) is modular modelling framework which enables generating a myriad of conceptual hydrological models by recombining elements from commonly-used models. This modular architecture makes it possible to (i) systematically understand the effects of model structural choices, (ii) generate large ensembles of alternative model structures, and (iii) represent spatial variability by combining process representations appropriate for different landscapes.

## Description and credits

The original implementation of the Framework for Understanding Structural Errors (FUSE1) is described in
[Clark et al. (2008)](https://doi.org/10.1029/2007WR006735). FUSE1 introduced a modular framework for constructing conceptual hydrologic models by recombining alternative representations of individual hydrologic processes. Subsequent developments extended the framework to support multiple numerical solution methods
([Clark and Kavetski, 2010](https://doi.org/10.1029/2009WR008894);
[Kavetski and Clark, 2010](https://doi.org/10.1029/2009WR008896))
and an optional temperature-index snow model
([Henn et al., 2015](https://doi.org/10.1002/2014WR016736)).

FUSE2 builds on these scientific foundations with a redesigned software architecture that supports reproducible workflows, scalable model execution, and differentiable hydrologic modeling.

Major additions in FUSE2 include

- an improved command-line interface supporting runtime configuration, parameter overrides, and multiple execution modes
- structured configuration files based on TOML that simplify experiment specification and future software extensions
- self-describing NetCDF-based input, output, and parameter files that facilitate visualization, analysis, and integration with external software
- extension from lumped to spatially distributed model configurations
- MPI-based parallel execution across large model domains
- a differentiable physics architecture for gradient-based calibration and machine learning;
- a modular software architecture that improves maintainability, extensibility, and reproducibility.

## License

FUSE is distributed under the GNU Public License Version 3. For details see the file `LICENSE` in the FUSE root directory or visit the [online version](https://www.gnu.org/licenses/gpl-3.0.html).
