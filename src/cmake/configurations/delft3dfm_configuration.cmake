# Specify the modules to be included

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dflowfm_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dwaq_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dwaves_configuration.cmake)

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/dimr_configuration.cmake)

# Additional includes

include(${CMAKE_SOURCE_DIR}/src/cmake/configurations/include/tools_delft3dfm_configuration.cmake)

# Project name must be at the end of the configuration: it might get a name when including other configurations and needs to overwrite that
project(delft3dfm)
