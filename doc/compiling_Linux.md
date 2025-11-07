# Compiling Delft3D on Linux

## Prerequisites

## Build steps
- build.sh
  Execute "./build.sh --help" to show the usage
  Currently used as default build process: "./build.sh fm-suite --compiler intel21"
  This will execute "src/setenv.sh" on Deltares systems. On other systems, the environment must be prepared upfront.
  For instructions, see [Setup your own Linux environment](Linux_setup.md).

## Alternative: without build-script
Refer to the [README](src\cmake\README) file in the CMake folder.
WARNING: When building without our build-script, the collection of the resulting binaries will need attention
