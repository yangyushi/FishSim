# FishSim

This is a collection of simulation code I wrote for my PhD project in the university of Bristol.


## Use the module

1. Copy `lib/fish_sim` into your project folder.
2. Add the path to `lib` to `$PYTHONPATH`. (For instance, you can edit the `~/.bashrc` file.)

## Dependency

The Python module requires the following packages: `numpy`, `matplotlib`, `scipy`, and `numba`.

They can be installed via command `python -m pip install numpy matplotlib scipy numba`

## \[Optional\] Compile the C++ Code

The following command will build a python package, binded by `pybind11`. The product will be a extra
Python model, `cmodel`. The simulation of models in `cmodel` is much faster, because they are written
in C++.

(You will need a C++ compiler and `CMake` to build the `cmodel` module.)

```sh
git submodule update --init  # download the Eigen & pybind11

mkdir build
cd build
cmake ..
make
make install
```

