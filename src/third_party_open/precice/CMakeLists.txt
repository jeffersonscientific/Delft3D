add_library(precice_fortran precice.F90)

if(UNIX)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(precice REQUIRED IMPORTED_TARGET libprecice)
    target_link_libraries(precice_fortran PRIVATE PkgConfig::precice)
endif(UNIX)
