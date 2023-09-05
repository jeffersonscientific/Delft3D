set(library_name blas)

# Define include directories
# get_module_include_path(${checkout_src_root}/${deltares_common_module} 
#                         ${library_name}
#                         deltares_common_include_path)

# message(STATUS "Configuring '${library_name}' with the following deltares_common include path: '${deltares_common_include_path}'")


# Note that the CMAKE_CURRENT_SOURCE_DIR has a path to THIS CMakeLists location
# set(include_path ${CMAKE_CURRENT_SOURCE_DIR}/include)
# include_directories(${deltares_common_include_path})

# Set directory which contains all source files
set(src_path "${CMAKE_CURRENT_SOURCE_DIR}/../src")

# Gather source files
get_fortran_source_files(${src_path}/ source)

#file(GLOB preprocessor_files ${src_path}/general/pldep.F90)

# Define library
add_library(${library_name} ${source})

# Set additional compilation properties
#set_source_files_properties(${preprocessor_files}
#                            PROPERTIES COMPILE_OPTIONS "${file_preprocessor_flag}")

# Set dependencies
#set(oss_dependencies    deltares_common
#                        deltares_common_mpi
#                        delftio
#                        ec_module
#                        flow2d3d_data
#                        flow2d3d_io_dol_f
#                        flow2d3d_kernel_dd_f
#                        morphology_data
#                        morphology_kernel) 
#oss_include_libraries(${library_name} oss_dependencies)

#if(UNIX)
#    target_link_libraries(${library_name} ${oss_dependencies})
#endif(UNIX)

# Define how the files should be structured within Visual Studio
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/.. FILES ${source})
set_target_properties (${library_name} PROPERTIES FOLDER third_party_open/blas)

# Change project flags at the end of this file
# target_compile_options(${library_name} PRIVATE "${qauto_threaded_flags}")
# target_compile_options(${library_name} PRIVATE "$<$<CONFIG:Debug>:${check_nobounds_flag}>")
