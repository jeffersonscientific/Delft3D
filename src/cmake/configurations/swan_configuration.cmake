# Specify the modules to be included
add_subdirectory(${CMAKE_SOURCE_DIR}/src/${swan_mpi_lib_module} swan_mpi_lib)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/${swan_mpi_module} swan_mpi)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/${swan_omp_module} swan_omp)

# netcdf
if(NOT TARGET netcdff)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${netcdf_module} netcdff)
endif()

# Project name must be at the end of the configuration: it might get a name when including other configurations and needs to overwrite that
project(swan)
