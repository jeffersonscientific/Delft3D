# Set Intel compiler specific flags:
enable_language (Fortran)
set(src_root_dir ${CMAKE_SOURCE_DIR}/..)

if (WIN32)
    message(FATAL_ERROR "GNU compilers are not supported on Windows. CMake will exit.")
endif(WIN32)


if (UNIX)
    # Set optional flags:
    message(STATUS "Setting Fortran compiler flags in Unix")

    # yoder: remove -fopenmp from all the things.
    #   As it is, it looks like swan  we can eitehr compile with MPI or OpenMP... or at lest we cannot compile with OpenMP (have not actually tried building without MPI)
    #   Having both triggers, "Error: invalid branch to/from OpenMP structured block"
    #   Unfortunately, there appear to be other packages that explicitly require OpenMP, so let's see if we can take out -pthread here, and add it into the CmakeLists.txt
    #   for those modules.
    #
    #set(CMAKE_CXX_FLAGS_RELEASE      "-O2 -fPIC -pthread")
    # -lm really does not belong here, but there are about 1000 CMakeLists.txt files, and so far most of them need target_link_libraries(${executable_name} PRIVATE m) (or similar...)
    #   so let's just try it and see how we go.
    set(CMAKE_CXX_FLAGS_RELEASE      "-O2 -fPIC -lm ")
    set(CMAKE_C_FLAGS_RELEASE        "-O2 -fPIC -lm ")
    set(CMAKE_Fortran_FLAGS          "-O2 -fPIC -lm -ffixed-line-length-132 -ffree-line-length-512 -fallow-argument-mismatch -cpp")
    set(CMAKE_CXX_FLAGS_DEBUG        "-g -O0 -fPIC ")
    set(CMAKE_C_FLAGS_DEBUG          "-g -O0 -fPIC ")
    
    set(cpp_compiler_flags           "-std=c++17")
    set(dialect                      "-std=f2008")
    set(bounds                       "-fbounds-check")

    set(file_preprocessor_flag       "-cpp")
    set(traceback_flag               "-fbacktrace")

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif(UNIX)

set(waq_default_flags ${file_preprocessor_flag} ${traceback_flag})
