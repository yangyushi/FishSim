cmake_minimum_required(VERSION 3.10)

# set the project name
project(
    cfish_sim VERSION 0.0
    DESCRIPTION "simulate the collective behaviour of fish as active matter"
    )
find_package(pybind11 REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

# setup usual compile flag
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# build the cpp targets
add_library(neighbour_list STATIC src/neighbour_list.cpp src/neighbour_list.h)
target_link_libraries(neighbour_list Eigen3::Eigen)

add_library(network STATIC src/network.cpp src/network.h)
target_link_libraries(network Eigen3::Eigen)

add_library(vicsek STATIC src/vicsek.cpp src/vicsek.h)
target_link_libraries(vicsek Eigen3::Eigen)

# build the python module
add_library(cfish_sim MODULE src/csimulate.cpp)

target_link_libraries(cfish_sim PRIVATE pybind11::module pybind11::lto Eigen3::Eigen neighbour_list network vicsek)
pybind11_extension(cfish_sim)
set_target_properties(cfish_sim PROPERTIES CXX_VISIBILITY_PRESET "hidden")

# move the libraries to lib folder
SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/lib/fish_sim)
install(TARGETS cfish_sim DESTINATION ${CMAKE_INSTALL_PREFIX})
