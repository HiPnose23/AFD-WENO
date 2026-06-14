# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitclone-lastrun.txt" AND EXISTS "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitinfo.txt" AND
  "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitclone-lastrun.txt" IS_NEWER_THAN "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/cluster/software/stacks/2024-06/sfos/git"
            clone --no-checkout --config "advice.detachedHead=false" "https://github.com/AMReX-Codes/amrex.git/" "amrex_code-src"
    WORKING_DIRECTORY "/cluster/home/mkryukov/Thesis/Amrex/build/_deps"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/AMReX-Codes/amrex.git/'")
endif()

execute_process(
  COMMAND "/cluster/software/stacks/2024-06/sfos/git"
          checkout "origin/development" --
  WORKING_DIRECTORY "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'origin/development'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/cluster/software/stacks/2024-06/sfos/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-src"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitinfo.txt" "/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/cluster/home/mkryukov/Thesis/Amrex/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/amrex_code-populate-gitclone-lastrun.txt'")
endif()
