# Note that the CMAKE_CURRENT_SOURCE_DIR has a path to THIS CMakeLists location
set(wave_root_path ${CMAKE_CURRENT_SOURCE_DIR}/../../) # Trace back to the D-Waves root

# Set directory which contains all source files
set(src_path src)

# Gather source files
get_fortran_source_files(${src_path} source_files)

# Define library
set(library_name wave_manager)
add_library(${library_name} ${source_files})

# Set dependencies
set(oss_dependencies    wave_data
                        delftio
                        deltares_common
                        ec_module
                        wave_io
                        wave_kernel) 

target_link_libraries(${library_name} PUBLIC ${oss_dependencies})

# Add dependencies
if (UNIX)
    # the `pkg_check_modules` function is created with this call
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(precice REQUIRED IMPORTED_TARGET libprecice)

    target_link_libraries(${library_name} PRIVATE PkgConfig::precice)
endif(UNIX)

include_directories(${mpi_module_path})

# Define how the files should be structured within Visual Studio
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${source_files})
set_target_properties (${library_name} PROPERTIES FOLDER engines_gpl/wave)
