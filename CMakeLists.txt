cmake_minimum_required(VERSION 3.22)
project(section_properties_cpp)

#set(QHULL_DIR /usr/local/Cellar/qhull/2020.2) # set the lib path
#include_directories(${QHULL_DIR}/include/libqhull) # include the lib
#include_directories(${QHULL_DIR}/include/libqhull_r) # include the lib
#include_directories(${QHULL_DIR}/include/libqhullcpp) # include the lib
#link_directories(${QHULL_DIR}/lib/) # link the lib

include_directories(/usr/local/include/libqhull) # include the lib
include_directories(/usr/local/include/libqhull_r) # include the lib
include_directories(/usr/local/include/libqhullcpp) # include the lib
link_directories(/usr/local/lib) # link the lib

set(CMAKE_CXX_STANDARD 17)

add_library(section_properties_cpp main.cpp pre/geometry.cpp pre/geometry.h pre/pre.cpp pre/pre.h
        libdistmesh/include/distmesh/constants.h libdistmesh/include/distmesh/distance_function.h
        libdistmesh/include/distmesh/distmesh.h libdistmesh/include/distmesh/functional.h
        libdistmesh/include/distmesh/triangulation.h libdistmesh/include/distmesh/utils.h
        libdistmesh/src/distance_function.cpp libdistmesh/src/distmesh.cpp libdistmesh/src/functional.cpp
        libdistmesh/src/triangulation.cpp libdistmesh/src/utils.cpp
    )

add_executable(section_properties_run main.cpp)

target_link_libraries(section_properties_cpp libqhull_r)
