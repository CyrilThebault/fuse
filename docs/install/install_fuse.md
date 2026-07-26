# Installing FUSE

FUSE has been successfully compiled and tested on Linux and macOS using a range
of modern Fortran compilers.

## Download FUSE

Clone the repository and initialize the required Git submodules:

```bash
git clone https://github.com/CH-Earth/fuse.git
cd fuse
git submodule update --init --recursive
```

The default branch is suitable for most users. If you wish to work with a
different development branch, switch branches after cloning. For example,

```bash
git checkout develop
```

or

```bash
git checkout staging
```

where

- **develop** contains the latest development version and is the primary
  branch for ongoing development; and
- **staging** contains code that is being prepared for the next release.

Alternatively, you may download the latest stable release from the
**Releases** page if you do not require the development version.

## Build dependencies

Before compiling FUSE, ensure the required compiler, libraries, and
build tools are installed. The required software is described in
**[Build dependencies](build_dependencies.md)**.

## Compile FUSE

From the top-level FUSE directory, run

```bash
cd build
make all
```

If compilation completes successfully, the executable will be written to

```text
bin/fuse.exe
```

If the build fails because required libraries or tools cannot be found, see
**[Build dependencies](build_dependencies.md)**.

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

The `test/CAN_05BB001/` directory contains examples that demonstrate
FUSE simulations for the Bow River above Banff, Alberta, Canada. These examples
provides a simple way to verify that FUSE has been installed correctly and that
the complete modeling workflow functions as expected.

The next section of this documentation, **[Running simulations](test_cases.md/)**,
provides a step-by-step guide for running the Bow River example and
interpreting the model output.

Additional example applications, including distributed and gridded model
configurations, are described in later sections of this documentation.

Once the Bow River example runs successfully, your FUSE installation is ready
for use.
