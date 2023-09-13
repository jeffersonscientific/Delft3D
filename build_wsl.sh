#!/bin/bash
#. /usr/share/Modules/init/bash

source /opt/intel/oneapi/setvars.sh
export PKG_CONFIG_PATH=/mnt/c/Development/petsc-3.19.4/arch-linux-c-opt/lib/pkgconfig/
cmake ./src/cmake -DCMAKE_Fortran_COMPILER=mpiifort -D CONFIGURATION_TYPE:STRING=DflowFM -D CMAKE_BUILD_TYPE=Debug -B build_all

cd build_all
cmake --build . -j --target install --config Debug