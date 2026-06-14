# This CMake snippet is adapted from the AMReX Guided-Tutorial.
# Original source: https://github.com/AMReX-Codes/amrex-tutorials
# Copyright (c) 2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory.
# Licensed under the BSD-3-Clause License.
# See LICENSE-AMReX.txt for full license details.

# To use a pre-installed AMReX build, run:
#    cmake -DAMReX_ROOT=/path/to/installdir
# Otherwise cmake will download AMReX from GitHub
if(NOT DEFINED AMReX_ROOT)
  message("-- Download and configure AMReX from GitHub")

  # Download AMReX from GitHub
  include(FetchContent)
  # set(FETCHCONTENT_QUIET OFF) # for more verbosity

  FetchContent_Declare(
    amrex_code
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
    GIT_TAG        origin/development
    )
    set(AMReX_FORTRAN            OFF)
    set(AMReX_FORTRAN_INTERFACES OFF)

  FetchContent_MakeAvailable(amrex_code)

else()

  # Add AMReX
  message("-- Searching for AMReX install directory at ${AMReX_ROOT}")
  message("${AMReX_ROOT}/lib/cmake/AMReX/AMReXConfig.cmake")
  find_package(AMReX
               PATHS ${AMReX_ROOT}/lib/cmake/AMReX/)

endif()
