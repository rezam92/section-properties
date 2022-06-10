cmake_minimum_required(VERSION 3.22)
project(section_properties_cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(section_properties_cpp main.cpp pre/geometry.cpp pre/geometry.h pre/pre.cpp pre/pre.h)
