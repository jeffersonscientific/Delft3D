install(SCRIPT "${CMAKE_SOURCE_DIR}/src/cmake/install_tests/install_and_bundle.cmake")
if(NOT DEFINED COMMIT_VERSION)
  set(COMMIT_VERSION "0")
endif()
configure_package_installer("Tests" "${CMAKE_SOURCE_DIR}/src/cmake/install_tests/ReadMe.txt" "1" "0" "${COMMIT_VERSION}" "STGZ;TGZ")

