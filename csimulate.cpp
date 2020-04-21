#include "vicsek.h"
#include "neighbour_list.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <string>

namespace py = pybind11;

const char* docstring_vicsek_3d_pbc =\
    "3D Vicsek simulation with scarlar (intrinsic) noise"\
    " and periodic boundary condition\n";

py::array_t<double> vicsek_3d_pbc(
        int n, double box, double r, double eta, double v0,
        int pre_steps, int run_steps, int jump
        ){
    /*
     * Run a vicsek simulation and return the positions & velocities as a (T, n, 6) array
     *
     * The system is inside a cubic periodic boundary with following parameters
     *     - n: particle number
     *     - box: size of box
     *     - r: interaction radius
     *     - eta: noise level
     *     - vo: speed
     *
     * The simulation setup is specified by following arguments
     *     - pre_steps: number of simulation steps to get into steady state
     *     - run_steps: the number of simulation steps performed after pre_steps
     *     - jump: every jump steps were recorded in the output numpy array
     */
    int update_step = floor(r / v0); // the frequency to update the cell list
    const int offset_frame = 6 * n;  // size of each frame
    int offset = 0;

    Vicsek3DPBC system{n, r, eta, box, v0};
    system.dump("test.xyz");

    py::array_t<double> result = py::array_t<double>(run_steps * n * 6);
    auto buffer_result = result.request();  // result buffer
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
        system.dump("test.xyz");
    }

    system.dump("test.xyz");

    for (int step = 0; step < run_steps * jump; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
        if (step % jump == 0){
            for (int i = 0; i < n; i++){
                for (int d = 0; d < 3; d++){
                    offset = (offset_frame * step) + (i * 6) + d;
                    ptr_result[offset] = system.positions(d, i);
                }
                for (int d = 0; d < 3; d++){
                    offset = (offset_frame * step) + (i * 6) + d + 3;
                    ptr_result[offset] = system.velocities(d, i);
                }
            }
        }
    }
    result.resize({run_steps, n, 6});
    return result;
}


PYBIND11_MODULE(csimulate, m){
    m.doc() = "common simulations for collective behaviours";
    m.def("vicsek_3d_pbc", &vicsek_3d_pbc, docstring_vicsek_3d_pbc,
           py::return_value_policy::copy);
}
