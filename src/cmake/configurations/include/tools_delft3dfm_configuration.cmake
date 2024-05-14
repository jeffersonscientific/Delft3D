# Tools for Delft3D-FM

# Tools_gpl

# Mormerge
if(NOT TARGET mormerge)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${mormerge_module} mormerge)
endif()

# dfmoutput
if(NOT TARGET dfmoutput)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${dfmoutput_module} dfmoutput)
endif()

# dfm_volume_tool
if(NOT TARGET dfm_volume_tool)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${dfm_volume_tool_module} dfm_volume_tool)
endif()

# dfm_api_access
if(NOT TARGET dfm_api_access)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${dfm_api_access_module} dfm_api_access)
endif()

# cosumo_bmi
if(NOT TARGET cosumo_bmi)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${cosumo_bmi_module} cosumo_bmi)
endif()

# D-Waq tools
include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/include/dwaq/dwaq_tools.cmake)