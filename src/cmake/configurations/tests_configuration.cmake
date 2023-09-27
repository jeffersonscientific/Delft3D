# Specify the modules to be included

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/dflowfm_configuration.cmake)

# Test binaries
add_subdirectory(${checkout_src_root}/${test_deltares_common_module} test_deltares_common)
add_subdirectory(${checkout_src_root}/${test_ec_module}              test_ec_module)
add_subdirectory(${checkout_src_root}/${test_waq_utils_f}            test_waq_utils_f)
add_subdirectory(${checkout_src_root}/${test_dflowfm_kernel}         test_dflowfm_kernel)

if(UNIX)
    # install
    add_subdirectory(${checkout_src_root}/${install_tests_module} install_tests)
endif(UNIX)

# Project name must be at the end of the configuration: it might get a name when including other configurations and needs to overwrite that
project(tests)
