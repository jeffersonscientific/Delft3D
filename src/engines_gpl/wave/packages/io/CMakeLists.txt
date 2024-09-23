# Set directory which contains all source files
set(src_path src)

# Gather source files
get_fortran_source_files(${src_path} source_files)

# Define library
set(library_name wave_io)
add_library(${library_name} ${source_files})

# Set dependencies
if (WIN32)
    set(oss_dependencies wave_data
                         deltares_common)
    oss_include_libraries(${library_name} oss_dependencies)
endif(WIN32)

if (UNIX)
    set(oss_dependencies wave_data
                         deltares_common)
    oss_include_libraries(${library_name} oss_dependencies)
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
endif(UNIX)

target_link_libraries(${library_name} PRIVATE netcdf-fortran::netcdf-fortran)

# Define how the files should be structured within Visual Studio
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${source_files})
set_target_properties (${library_name} PROPERTIES FOLDER engines_gpl/wave)
