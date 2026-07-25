# Build dependencies

The required compiler, libraries, and build tools must be installed in order
for FUSE to compile successfully, 

## GNU Fortran

FUSE has been successfully compiled using

- GNU Fortran (`gfortran`, version 6 or later)
- Intel Fortran (`ifort`)

Other modern Fortran compilers should also work.

On macOS, install the GNU Compiler Collection (GCC), which includes
`gfortran`:

```bash
brew install gcc
```

## Required libraries and build tools

FUSE requires

- CMake (used to build the bundled `toml-f` library);
- the NetCDF-C and NetCDF-Fortran libraries;
- `pkg-config` (used to locate installed libraries during compilation).

On macOS, these can be installed with Homebrew:

```bash
brew install cmake
brew install netcdf
brew install netcdf-fortran
brew install pkg-config
```

Homebrew automatically installs HDF5 as a dependency of NetCDF.

The `toml-f` library is included with the FUSE source code and is built
automatically during compilation.
