cmake_minimum_required(VERSION 3.15)

# include AMReX setup
include(cmake/SetupAMReX.cmake)

# project name and source file languages
project(weno_amrex LANGUAGES C CXX)

# add files
add_executable(weno_3d_amrex weno_3d_amrex.cpp)
add_executable(radial weno_3d_radial.cpp)
add_executable(test_radial test_radial.cpp)
add_executable(euler euler_v2.cpp)

# require C++20 for std::format
target_compile_features(weno_3d_amrex PUBLIC cxx_std_20)
set_target_properties(weno_3d_amrex PROPERTIES CXX_EXTENSIONS OFF)

# link to AMReX
target_link_libraries(weno_3d_amrex PRIVATE AMReX::amrex)

# require C++20 for std::format
target_compile_features(radial PUBLIC cxx_std_20)
set_target_properties(radial PROPERTIES CXX_EXTENSIONS OFF)

target_compile_features(test_radial PUBLIC cxx_std_20)
set_target_properties(test_radial PROPERTIES CXX_EXTENSIONS OFF)

# link to AMReX
target_link_libraries(test_radial PRIVATE AMReX::amrex)
target_link_libraries(radial PRIVATE AMReX::amrex)
target_compile_features(euler PUBLIC cxx_std_20)
set_target_properties(euler PROPERTIES CXX_EXTENSIONS OFF)

# link to AMReX
target_link_libraries(euler PRIVATE AMReX::amrex)

