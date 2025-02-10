message(STATUS "Installing Intel Redistributable Libraries")

if(WIN32)
  include(intel_redist_windows.cmake)
endif(WIN32)

if(UNIX)
  include(intel_redist_linux.cmake)
endif(UNIX)
