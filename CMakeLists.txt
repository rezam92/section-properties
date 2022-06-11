cmake_minimum_required(VERSION 3.22)
project(section_properties_cpp)

set(CMAKE_CXX_STANDARD 17)

include_directories("qhull")

#add_executable(qhull qhull/src/libqhull/libqhull.c qhull/src/libqhull/libqhull.h)

add_executable(section_properties_cpp main.cpp pre/geometry.cpp pre/geometry.h pre/pre.cpp pre/pre.h
        libdistmesh/include/distmesh/constants.h libdistmesh/include/distmesh/distance_function.h
        libdistmesh/include/distmesh/distmesh.h libdistmesh/include/distmesh/functional.h
        libdistmesh/include/distmesh/triangulation.h libdistmesh/include/distmesh/utils.h
        libdistmesh/src/distance_function.cpp libdistmesh/src/distmesh.cpp libdistmesh/src/functional.cpp
        libdistmesh/src/triangulation.cpp libdistmesh/src/utils.cpp

    )

target_link_libraries(section_properties_cpp qhull/libqhullstatic.a)
target_link_libraries(section_properties_cpp qhull/libqhull_r.so)
target_link_libraries(section_properties_cpp qhull/libqhullcpp.a)