set(executable_name test_io_netcdf)

file(GLOB source src/*.f90)

add_executable(${executable_name} ${source})

set(oss_dependencies    fortrangis
                        ftnunit
                        io_netcdf
                        netcdf-fortran::netcdf-fortran)

if (UNIX)
    target_link_libraries(${executable_name} ${oss_dependencies})
    set_property(TARGET ${executable_name} PROPERTY LINKER_LANGUAGE Fortran)
endif(UNIX)

if (WIN32)
    oss_include_libraries(${executable_name} oss_dependencies)
    target_link_libraries(${executable_name} ${oss_dependencies})
endif(WIN32)

# Set additional compilation properties, specific for Debug mode
target_compile_options(${executable_name} PRIVATE "$<$<CONFIG:Debug>:${check_pointers_flag}>")
target_compile_options(${executable_name} PRIVATE "$<$<CONFIG:Debug>:${floating_point_exception_flag}>")

# Create the folder structure in vfproj
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${source})
set_target_properties (${executable_name} PROPERTIES FOLDER tests/utils_lgpl/io_netcdf)

install(TARGETS ${executable_name} RUNTIME DESTINATION bin)