# Set Intel compiler specific flags:
enable_language (Fortran)
set(src_root_dir ${CMAKE_SOURCE_DIR}/..)

if (WIN32)
    message(FATAL_ERROR "GNU compilers are not supported on Windows. CMake will exit.")
endif(WIN32)

add_library(all_compiler_warnings INTERFACE)
set(all_warning_flags "-Wall" "-pedantic")
target_compile_options(all_compiler_warnings INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${all_warning_flags}>")

add_library(compiler_warnings_as_errors INTERFACE)
set(warning_error_flag "-Werror")
target_compile_options(compiler_warnings_as_errors INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${warning_error_flag}>")

add_library(limit_compiler_warnings INTERFACE)
set(disabled_warning_flags "")
target_compile_options(limit_compiler_warnings INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${disabled_warning_flags}>")

add_library(no_compiler_warnings INTERFACE)
set(no_warning_flags "-w")
target_compile_options(no_compiler_warnings INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${no_warning_flags}>")

if (UNIX)
    # Set optional flags:
    message(STATUS "Setting Fortran compiler flags in Unix")

    set(CMAKE_CXX_FLAGS_RELEASE      "-O2 -fPIC -fopenmp")
    set(CMAKE_C_FLAGS_RELEASE        "-O2 -fPIC -fopenmp")
    set(CMAKE_Fortran_FLAGS          "-O2 -fPIC -fopenmp -ffixed-line-length-132 -ffree-line-length-512 -fallow-argument-mismatch")
    set(CMAKE_CXX_FLAGS_DEBUG        "-g -O0 -fPIC -fopenmp")
    set(CMAKE_C_FLAGS_DEBUG          "-g -O0 -fPIC -fopenmp")
    
    set(cpp_compiler_flags           "-std=c++17")
    set(dialect                      "-std=f2008")
    set(bounds                       "-fbounds-check")

    set(file_preprocessor_flag       "-cpp")
    set(traceback_flag               "-fbacktrace")

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif(UNIX)

set(waq_default_flags ${file_preprocessor_flag} ${traceback_flag})
