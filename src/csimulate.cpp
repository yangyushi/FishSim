#include "vicsek.h"
#include "neighbour_list.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <string>

namespace py = pybind11;


const char* docstring_vicsek_3d_pbc =\
    "3D Vicsek simulation with scarlar (intrinsic) noise"\
    " and periodic boundary condition\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
    "   eta (:obj:`float`): the magnitude of noise, from 0 to 1\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 6)\n"\
    ;


const char* docstring_attanasi2014pcb =\
    "3D Vicsek simulation in paper 10.1371/journal.pcbi.1003697\n"\
    "(scalar noise + harmonic motion in 3D)\n"\
    " \n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   eta (:obj:`float`): the magnitude of noise, from 0 to 1\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   b (:obj:`float`): the beta parameter controlling the"\
    "   strength of the resotring force of the harmonic oscillator\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n"\
    "   output_file (:obj:`str`) = "": output running frames into an xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 6)\n"\
    ;


const char* docstring_vicsek_3d_pbc_inertia =\
    "3D Vicsek simulation in PBC with a inertia term\n"\
    " \n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
    "   eta (:obj:`float`): the magnitude of noise, from 0 to 1\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   alpha (:obj:`float`): larger alpha corresponds to larger mass, ranging from 0 to 1\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n"\
    "   output_file (:obj:`str`) = "": output running frames into an xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 6)\n"\
    ;


const char* docstring_vicsek_3d_pbc_inertia_af =\
    "3D Vicsek simulation in PBC with a inertia term with J_ij = -1 (J_ii = 1)\n"\
    " \n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
    "   eta (:obj:`float`): the magnitude of noise, from 0 to 1\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   alpha (:obj:`float`): larger alpha corresponds to larger mass, ranging from 0 to 1\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n"\
    "   output_file (:obj:`str`) = "": output running frames into an xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 6)\n"\
    ;



const char* docstring_vicsek_2d_pbc =\
    "2D Vicsek simulation with scarlar (intrinsic) noise"\
    " and periodic boundary condition"\
    "\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
    "   eta (:obj:`float`): the magnitude of noise, from 0 to 1\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 4)\n"\
    ;


const char* docstring_vicsek_2d_pbc_vn =\
    "2D Vicsek simulation with vectorial (extrinsic) noise"\
    " and periodic boundary condition"\
    "\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
    "   eta (:obj:`float`): the magnitude of noise, from 0 to 1\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 4)\n"\
    ;


py::array_t<double> vicsek_3d_pbc(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1, string load_file=""
        ){
    int update_step = floor(r / v0); // the frequency to update the cell list
    if (update_step < 1){ update_step = 1; }
    const int offset_frame = 6 * n;  // size of each frame
    int offset = 0;

    Vicsek3DPBC system{n, r, eta, box, v0};

    if (load_file.size() > 0){
        system.load(load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * n * 6);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
    }

    int cursor = 0;
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
                    offset = (offset_frame * cursor) + (i * 6) + d;
                    ptr_result[offset] = system.positions_(d, i);
                }
                for (int d = 0; d < 3; d++){
                    offset = (offset_frame * cursor) + (i * 6) + d + 3;
                    ptr_result[offset] = system.velocities_(d, i);
                }
            }
            cursor++;
        }
    }
    result.resize({run_steps, n, 6});
    return result;
}


py::array_t<double> vicsek_3d_pbc_inertia(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1,
        string load_file="", string output_file=""
        ){
    int update_step = floor(r / v0); // the frequency to update the cell list
    if (update_step < 1){ update_step = 1; }
    const int offset_frame = 6 * n;  // size of each frame
    int offset = 0;

    bool output = output_file.size() > 0;

    Vicsek3DPBCInertia system{n, r, eta, box, v0, alpha};

    if (load_file.size() > 0){
        system.load(load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * n * 6);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
    }

    int cursor = 0;
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
                    offset = (offset_frame * cursor) + (i * 6) + d;
                    ptr_result[offset] = system.positions_(d, i);
                }
                for (int d = 0; d < 3; d++){
                    offset = (offset_frame * cursor) + (i * 6) + d + 3;
                    ptr_result[offset] = system.velocities_(d, i);
                }
            }
            cursor++;
            if (output){
                system.dump(output_file);
            }
        }
    }
    result.resize({run_steps, n, 6});
    return result;
}


py::array_t<double> vicsek_3d_pbc_inertia_af(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1,
        string load_file="", string output_file=""
        ){
    int update_step = floor(r / v0); // the frequency to update the cell list
    if (update_step < 1){ update_step = 1; }
    const int offset_frame = 6 * n;  // size of each frame
    int offset = 0;

    bool output = output_file.size() > 0;

    Vicsek3DPBCInertiaAF system{n, r, eta, box, v0, alpha};

    if (load_file.size() > 0){
        system.load(load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * n * 6);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
    }

    int cursor = 0;
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
                    offset = (offset_frame * cursor) + (i * 6) + d;
                    ptr_result[offset] = system.positions_(d, i);
                }
                for (int d = 0; d < 3; d++){
                    offset = (offset_frame * cursor) + (i * 6) + d + 3;
                    ptr_result[offset] = system.velocities_(d, i);
                }
            }
            cursor++;
            if (output){
                system.dump(output_file);
            }
        }
    }
    result.resize({run_steps, n, 6});
    return result;
}


py::array_t<double> attanasi2014pcb(
        int n, double eta, double v0, double b, double r,
        int pre_steps, int run_steps, int jump=1,
        string load_file="", string output_file=""
        ){
    int update_step = floor(r / v0); // the frequency to update the cell list
    if (update_step < 1){ update_step = 1; }
    const int offset_frame = 6 * n;  // size of each frame
    int offset = 0;
    bool output = output_file.size() > 0;

    Attanasi2014PCB system{n, r, eta, v0, b};

    if (load_file.size() > 0){
        system.load(load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * n * 6);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
    }
    int cursor = 0;
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
                    offset = (offset_frame * cursor) + (i * 6) + d;
                    ptr_result[offset] = system.positions_(d, i);
                }
                for (int d = 0; d < 3; d++){
                    offset = (offset_frame * cursor) + (i * 6) + d + 3;
                    ptr_result[offset] = system.velocities_(d, i);
                }
            }
            cursor++;
            if (output){
                system.dump(output_file);
            }
        }
    }
    result.resize({run_steps, n, 6});
    return result;
}


py::array_t<double> vicsek_2d_pbc(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1, string load_file=""
        ){
    int update_step = floor(r / v0); // the frequency to update the cell list
    if (update_step < 1){ update_step = 1; }
    const int offset_frame = 4 * n;  // size of each frame
    int offset = 0;

    Vicsek2DPBC system{n, r, eta, box, v0};

    if (load_file.size() > 0){
        system.load(load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * n * 4);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
    }

    int cursor = 0;
    for (int step = 0; step < run_steps * jump; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
        if (step % jump == 0){
            for (int i = 0; i < n; i++){
                for (int d = 0; d < 2; d++){
                    offset = (offset_frame * cursor) + (i * 4) + d;
                    ptr_result[offset] = system.positions(d, i);
                }
                for (int d = 0; d < 2; d++){
                    offset = (offset_frame * cursor) + (i * 4) + d + 2;
                    ptr_result[offset] = system.velocities(d, i);
                }
            }
            cursor++;
        }
    }
    result.resize({run_steps, n, 4});
    return result;
}


py::array_t<double> vicsek_2d_pbc_vn(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1, string load_file=""
        ){
    int update_step = floor(r / v0); // the frequency to update the cell list
    if (update_step < 1){ update_step = 1; }
    const int offset_frame = 4 * n;  // size of each frame
    int offset = 0;

    Vicsek2DPBCVN system{n, r, eta, box, v0};

    if (load_file.size() > 0){
        system.load(load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * n * 4);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
    }

    int cursor = 0;
    for (int step = 0; step < run_steps * jump; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else{
            system.move(false);
        }
        if (step % jump == 0){
            for (int i = 0; i < n; i++){
                for (int d = 0; d < 2; d++){
                    offset = (offset_frame * cursor) + (i * 4) + d;
                    ptr_result[offset] = system.positions(d, i);
                }
                for (int d = 0; d < 2; d++){
                    offset = (offset_frame * cursor) + (i * 4) + d + 2;
                    ptr_result[offset] = system.velocities(d, i);
                }
            }
            cursor++;
        }
    }
    result.resize({run_steps, n, 4});
    return result;
}


PYBIND11_MODULE(csimulate, m){
    m.doc() = "common simulations for collective behaviours";

    m.def(
            "vicsek_3d_pbc", &vicsek_3d_pbc, docstring_vicsek_3d_pbc,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"), py::arg("jump")=1,
            py::arg("load_file")=""
           );

    m.def(
            "vicsek_3d_pbc_inertia", &vicsek_3d_pbc_inertia, docstring_vicsek_3d_pbc_inertia,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"), py::arg("alpha"),
            py::arg("pre_steps"), py::arg("run_steps"), py::arg("jump")=1,
            py::arg("load_file")="", py::arg("output_file")
           );

    m.def(
            "vicsek_3d_pbc_inertia_af", &vicsek_3d_pbc_inertia_af, docstring_vicsek_3d_pbc_inertia_af,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"), py::arg("alpha"),
            py::arg("pre_steps"), py::arg("run_steps"), py::arg("jump")=1,
            py::arg("load_file")="", py::arg("output_file")
           );

    m.def(
            "attanasi2014pcb", &attanasi2014pcb, docstring_attanasi2014pcb,
            py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("b"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"), py::arg("jump")=1,
            py::arg("load_file")="", py::arg("output_file")=""
           );

    m.def(
            "vicsek_2d_pbc", &vicsek_2d_pbc, docstring_vicsek_2d_pbc,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")=""
           );

    m.def(
            "vicsek_2d_pbc_vn", &vicsek_2d_pbc_vn, docstring_vicsek_2d_pbc_vn,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")=""
           );
}
