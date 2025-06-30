if(WIN32)
    # On Windows, build Google Test from source
    set(BUILD_GMOCK OFF)
    set(INSTALL_GTEST OFF)
    set(gtest_force_shared_crt ON)

    add_subdirectory(googletest-main)

    set_target_properties(gtest PROPERTIES FOLDER third_party_open/googletest)
    set_target_properties(gtest_main PROPERTIES FOLDER third_party_open/googletest)
else()
    # On Linux, use pkg-config to find sytem (docker) installed Google Test
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GTEST REQUIRED gtest)
    find_library(GTEST_STATIC_LIB NAMES libgtest.a gtest REQUIRED)

    if(NOT TARGET gtest)
        add_library(gtest STATIC IMPORTED)
        set_target_properties(gtest PROPERTIES
            IMPORTED_LOCATION ${GTEST_STATIC_LIB}
            INTERFACE_INCLUDE_DIRECTORIES "${GTEST_INCLUDE_DIRS}"
            INTERFACE_COMPILE_OPTIONS "${GTEST_CFLAGS_OTHER}"
        )
    endif()
endif()
