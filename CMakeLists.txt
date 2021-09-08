cmake_minimum_required(VERSION 3.10)

# set the project name
project(
    cfish_sim VERSION 0.0
    DESCRIPTION "simulate the collective behaviour of fish as active matter"
)

find_package(OpenMP)
add_subdirectory(extern/pybind11)
add_subdirectory(extern/eigen3 EXCLUDE_FROM_ALL)

# setup usual compile flag
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -O3 -fPIC)

set(CMAKE_CXX_STANDARD 14 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# build the cpp targets
add_library(core STATIC src/core.cpp src/core.hpp)
target_link_libraries(core eigen)

add_library(neighbour_list STATIC src/neighbour_list.cpp src/neighbour_list.hpp)
target_link_libraries(neighbour_list core eigen)

add_library(network STATIC src/network.cpp src/network.hpp)
target_link_libraries(network core eigen)

add_library(vicsek STATIC src/vicsek.cpp src/vicsek.hpp)
target_link_libraries(vicsek core neighbour_list eigen)

if(OpenMP_CXX_FOUND)
    target_link_libraries(neighbour_list OpenMP::OpenMP_CXX)
    target_link_libraries(core OpenMP::OpenMP_CXX)
endif()

# build the python module
add_library(cmodel MODULE src/cmodel.cpp)
target_link_libraries(
    cmodel PRIVATE pybind11::module pybind11::lto
    core neighbour_list network vicsek eigen
    )
pybind11_extension(cmodel)
set_target_properties(cmodel PROPERTIES CXX_VISIBILITY_PRESET "hidden")

# move the libraries to lib folder
SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/lib/fish_sim)
install(TARGETS cmodel DESTINATION ${CMAKE_INSTALL_PREFIX})
