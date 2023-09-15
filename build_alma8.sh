#!/bin/bash
#. /usr/share/Modules/init/bash

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:'./src/third_party_open/netcdf/netCDF 4.6.1/lib/pkgconfig/'
module load intel/2023.2.1
module load intelmpi/2021.10.0
module load petsc/3.19.0_gcc12.2.0


cmake ./src/cmake -DCMAKE_Fortran_COMPILER=mpiifort -D CONFIGURATION_TYPE:STRING=DflowFM -D CMAKE_BUILD_TYPE=Debug -B build_all

cd build_all
cmake --build . -j --target install --config Debug