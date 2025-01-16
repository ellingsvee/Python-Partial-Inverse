cmake_minimum_required(VERSION 3.4...3.18)
project(partial_inverse)

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

add_subdirectory(pybind11)
pybind11_add_module(partial_inverse partial_inverse.cpp)