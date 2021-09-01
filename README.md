# FishSim

This is a collection of simulation code I wrote for my PhD project in the university of Bristol.

## Dependency

You need to have a C++ compiler and CMake to build the C++ source code.

The Python module requires the following packages: `numpy`, `matplotlib`, `scipy`, and `numba`.

## Install

the following command will build the package.

```sh
git submodule update --init  # download the Eigen & pybind11

mkdir build
cd build
cmake ..
make
make install
```

## Use the module

1. Copy `lib/fish_sim` into your project folder
2. Add the path to `lib` to $PYTHONPATH
