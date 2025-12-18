# Shared definitions
set(wave_root_path ${CMAKE_CURRENT_SOURCE_DIR}/../../) # Trace back to the D-Waves root
set(wave_version_path ${wave_root_path}/version)
set(rc_version_file ${wave_version_path}/wave_version.rc)
set(version_include_dir ${CMAKE_SOURCE_DIR}/../version_includes)

# Set directory which contains all source files
set(src_path src)

# Wave Library definition
include(wave_dll.cmake)

# Wave Executable definition
include(wave_exe.cmake)

# yoder:
# OpenMP should be all in the include() bits, above
#project(dwaves LANGUAGES C CXX Fortran)
#find_package(OpenMP REQUIRED)
#find_package(OpenMP REQUIRED COMPONENTS Fortran)
#target_link_libraries(${library_name} PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran OpenMP::OPenMP_C)
#target_link_libraries(${executable_name} PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran OpenMP::OPenMP_C)


#
#if (OPENMP_FOUND)
#    #
#    # older syntax:
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

