#include "vicsek.h"
#include "spin.h"
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


const char* docstring_ism_3d =\
    "3D Inertial spin model simulation"\
    " with open space periodic boundary condition\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   T (:obj:`float`): the magnitude of noise, noted as the temperature\n"\
    "   j (:obj:`float`): the coupling constant\n"\
    "   m (:obj:`float`): the strength of the inertial effect like the mass\n"\
    "   f (:obj:`float`): the friction coefficent\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 6)\n"\
    ;


const char* docstring_ism_3d_pbc =\
    "3D Inertial spin model simulation"\
    " with periodic boundary condition\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
    "   r (:obj:`float`): the interaction range\n"\
    "   v0 (:obj:`float`): the speed in the simulation\n"\
    "   T (:obj:`float`): the magnitude of noise, noted as the temperature\n"\
    "   j (:obj:`float`): the coupling constant\n"\
    "   m (:obj:`float`): the strength of the inertial effect like the mass\n"\
    "   f (:obj:`float`): the friction coefficent\n"\
    "   pre_steps (:obj:`int`): the simulation steps running without writing to the result\n"\
    "                    typically this is the simulation steps to reach a steady state\n"\
    "   run_steps (:obj:`int`): the simulation steps that write to the result\n"\
    "   jump (:obj:`int`) = 1: every ``jump`` step is writting into the result\n"\
    "   load_file (:obj:`str`) = "": load initial configuration from a xyz file\n\n"\
    "Return:\n"\
    "   :obj:`numnp.ndarray`: the positions and velocities of particles"\
    " in each frame, shape (run_steps, n, 6)\n"\
    ;



/*
 * get the frequency to update the neighbour list
 */
int get_update_step(double r, double v0){
    int update_step = floor(r / v0);
    if (update_step < 1){
        update_step = 1;
    }
    return update_step;
}

template<class T>
py::array_t<double> simulate(
        T system, int update_step, int pre_steps, int run_steps, int jump,
        string load_file, string dump_file
        ){
    int dim = system.positions_.rows();
    const int offset_frame = 2 * dim * system.n_;  // size of each frame
    int offset = 0;
    bool should_dump = dump_file.size() > 0;

    if (load_file.size() > 0){
        load(system, load_file);
    }

    py::array_t<double> result = py::array_t<double>(run_steps * system.n_ * dim * 2);
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
            for (int i = 0; i < system.n_; i++){
                for (int d = 0; d < dim; d++){
                    offset = (offset_frame * cursor) + (i * 2 * dim) + d;
                    ptr_result[offset] = system.positions_(d, i);
                }
                for (int d = 0; d < dim; d++){
                    offset = (offset_frame * cursor) + (i * 2 * dim) + d + dim;
                    ptr_result[offset] = system.velocities_(d, i);
                }
            }
            cursor++;
            if (should_dump){
                dump(system, dump_file);
            }
        }
    }
    result.resize({run_steps, system.n_, 2 * dim});

    return result;
}

py::array_t<double> vicsek_3d_pbc(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    Vicsek3DPBC system{n, r, eta, box, v0};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> vicsek_3d_pbc_inertia(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    Vicsek3DPBCInertia system{n, r, eta, box, v0, alpha};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> vicsek_3d_pbc_inertia_af(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    Vicsek3DPBCInertiaAF system{n, r, eta, box, v0, alpha};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> attanasi2014pcb(
        int n, double eta, double v0, double b, double r,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    Attanasi2014PCB system{n, r, eta, v0, b};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> vicsek_2d_pbc(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    Vicsek2DPBC system{n, r, eta, box, v0};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> vicsek_2d_pbc_vn(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    Vicsek2DPBCVN system{n, r, eta, box, v0};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> ism_3d(
        int n, double r, double v0,
        double T, double j, double m, double f,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    InertialSpin3D system{n, r, v0, T, j, m, f};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


py::array_t<double> ism_3d_pbc(
        int n, double box, double r, double v0,
        double T, double j, double m, double f,
        int pre_steps, int run_steps, int jump=1, string load_file="", string dump_file=""
        ){
    InertialSpin3DPBC system{n, box, r, v0, T, j, m, f};
    int update_step = get_update_step(r, v0);
    return simulate(system, update_step, pre_steps, run_steps, jump, load_file, dump_file);
}


PYBIND11_MODULE(csimulate, m){
    m.doc() = "common simulations for collective behaviours";

    m.def(
            "vicsek_3d_pbc", &vicsek_3d_pbc, docstring_vicsek_3d_pbc,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")=""
           );

    m.def(
            "vicsek_3d_pbc_inertia", &vicsek_3d_pbc_inertia, docstring_vicsek_3d_pbc_inertia,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"), py::arg("alpha"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")=""
           );

    m.def(
            "vicsek_3d_pbc_inertia_af", &vicsek_3d_pbc_inertia_af, docstring_vicsek_3d_pbc_inertia_af,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"), py::arg("alpha"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")=""
           );

    m.def(
            "attanasi2014pcb", &attanasi2014pcb, docstring_attanasi2014pcb,
            py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("b"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")=""
           );

    m.def(
            "vicsek_2d_pbc", &vicsek_2d_pbc, docstring_vicsek_2d_pbc,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")=""
           );

    m.def(
            "vicsek_2d_pbc_vn", &vicsek_2d_pbc_vn, docstring_vicsek_2d_pbc_vn,
            py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file=")=""
           );
    m.def(
            "ism_3d", &ism_3d, docstring_ism_3d,
            py::arg("n"), py::arg("r"), py::arg("v0"),
            py::arg("T"), py::arg("j"), py::arg("m"), py::arg("f"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file=")=""
           );
    m.def(
            "ism_3d_pbc", &ism_3d_pbc, docstring_ism_3d_pbc,
            py::arg("n"), py::arg("box"), py::arg("r"), py::arg("v0"),
            py::arg("T"), py::arg("j"), py::arg("m"), py::arg("f"),
            py::arg("pre_steps"), py::arg("run_steps"),
            py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file=")=""
           );
}
