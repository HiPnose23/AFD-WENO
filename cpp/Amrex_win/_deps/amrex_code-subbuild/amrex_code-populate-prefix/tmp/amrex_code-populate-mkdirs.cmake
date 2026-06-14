# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-src")
  file(MAKE_DIRECTORY "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-src")
endif()
file(MAKE_DIRECTORY
  "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-build"
  "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix"
  "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix/tmp"
  "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp"
  "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src"
  "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/max/Desktop/Thesis/code/cpp/Amrex/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
