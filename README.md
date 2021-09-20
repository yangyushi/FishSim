# FishSim

This is a collection of simulation code I wrote for my PhD project in the university of Bristol.

![](couzin.gif)



I implemented the following models for the collective behaviour.

- The voter model (JAP, MA, CH, HL, *PRE*, **2008**)
- The Vectorial Network Model (MA, CH, *J. Stat. Phys.*, **2003**)
- The Vicsek Model (TV, AC, EBJ, IC, OS, *PRL*, **1995**)
- The Couzin Model (IDC, JK, GDR, NRF, *J. theor. Biol.*, **2002**)


## Use the module

1. Copy `lib/fish_sim` into your project folder.
2. Add the path to `lib` to `$PYTHONPATH`. (For instance, you can edit the `~/.bashrc` file.)

(Some examples are in the `tests` folder.)

## Dependency

The Python module requires the following packages: `numpy`, `matplotlib`, `scipy`, and `numba`.

They can be installed via command `python -m pip install numpy matplotlib scipy numba`

## \[Optional\] Fast Models

The following command will build an extra python package `cmodel`. It is written in C++, and binded by `pybind11`. The simulation of models in `cmodel` is faster, and the way to use these models is consistent with the Python models.

(You will need a C++ compiler and `CMake` to build the `cmodel` module.)

```sh
git submodule update --init  # download the Eigen & pybind11

mkdir build
cd build
cmake ..
make
make install
```

