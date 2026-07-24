# Installing FUSE

FUSE has been successfully compiled and tested on Linux and macOS using a range
of modern Fortran compilers.

## Download FUSE

The source code is available from

https://github.com/CH-EARTH/fuse

There are several ways to obtain FUSE:

- download the latest stable release from the **Releases** page;
- clone the `develop` branch to obtain the latest development version; or
- fork the repository if you plan to contribute to FUSE development.

## Compile FUSE

From the top-level FUSE directory, run

```bash
cd build
make
```

If compilation completes successfully, the executable will be written to

```text
bin/fuse.exe
```

If the build fails because required libraries or tools cannot be found, see
**Build dependencies** below.

## Verify the installation

Run

```bash
bin/fuse.exe
```

When executed without command-line arguments, FUSE prints a short usage message
listing the required command-line arguments. This confirms that the executable
was built successfully.

To display the complete command-line interface, run

```bash
bin/fuse.exe --help
```

which lists all execution modes, required arguments, and optional settings.

## Run the example application

The `test/` directory contains a single example application that demonstrates a
lumped FUSE simulation for the Bow River above Banff, Alberta, Canada. This example
provides a simple way to verify that FUSE has been installed correctly and that
the complete modeling workflow functions as expected.

The next section of this documentation, **[Quick Start](quick-start/)**,
provides a step-by-step guide for running the Bow River example and
interpreting the model output.

Additional example applications, including distributed and gridded model
configurations, are described in later sections of this documentation.

Once the Bow River example runs successfully, your FUSE installation is ready
for use.

---

# Build dependencies

If FUSE does not compile successfully, ensure the required compiler, libraries,
and build tools are installed.

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
