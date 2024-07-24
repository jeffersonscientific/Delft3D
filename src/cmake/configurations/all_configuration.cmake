# Specify the modules to be included

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/include/dflowfm_configuration_basic.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dwaq_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dwaves_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dimr_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/include/tools_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/include/windows_postbuild_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/drr_configuration.cmake)

# Project name must be at the end of the configuration: it might get a name when including other configurations and needs to overwrite that
project(all)
