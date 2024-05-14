# RTC Real Time Control
# =====================
add_subdirectory(${CMAKE_SOURCE_DIR}/src/${rtc_module} rtc)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/${rtc_plugin_c_module} plugin_rtc_c)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/${rtc_kernel_module} rtc_kernel)


# Utils
# =====
if(NOT TARGET rr_rtc_tools)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${rr_rtc_tools_module} rr_rtc_tools)
endif()

if(NOT TARGET wl_openmi_support)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${wl_openmi_support_module} wl_openmi_support)
endif()


# Utils LGPL
# =====

if(NOT TARGET delftio_shm)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${delftio_shm_module} delftio_shm)
endif()

if(NOT TARGET delftio)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${delftio_module} delftio)
endif()

if(NOT TARGET deltares_common)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${deltares_common_module} deltares_common)
endif()

if(NOT TARGET deltares_common_c)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${deltares_common_c_module} deltares_common_c)
endif()

if(NOT TARGET io_netcdf)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${io_netcdf_module} io_netcdf)
endif()

if(NOT TARGET io_netcdf_data)
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/${io_netcdf_data_module} io_netcdf_data)
endif()

# Third party
# ===========

if(WIN32)
    if(NOT TARGET netcdff)
        add_subdirectory(${CMAKE_SOURCE_DIR}/src/${netcdf_module} netcdff)
    endif()
endif()

# Project name must be at the end of the configuration: it might get a name when including other configurations and needs to overwrite that
project(rtc)
