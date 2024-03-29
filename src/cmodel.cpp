#include "core.hpp"
#include "vicsek.hpp"
#include "couzin.hpp"
#include "network.hpp"
#include "boundary.hpp"
#include "montecarlo.hpp"
#include "neighbour_list.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <string>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


namespace py = pybind11;


const char* docstring_vicsek_2d =\
    "2D Vicsek simulation with scarlar (intrinsic) noise"\
    " and periodic boundary condition"\
    "\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
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


const char* docstring_vicsek_3d =\
    "3D Vicsek simulation with scarlar (intrinsic) noise"\
    "\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
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


const char* docstring_vicsek_2d_pbc_inertia =\
    "2D Vicsek simulation in PBC with a inertia term\n"\
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


const char* docstring_vicsek_2d_pbc_vn =\
    "2D Vicsek simulation with vectorial (extrinsic) noise"\
    " and periodic boundary condition"\
    "\n\n\n"\
    "Args:\n\n"\
    "   n (:obj:`int`): number of particles in the system\n"\
    "   box (:obj:`float`): the size of the cubic box\n"\
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
        std::string load_file, std::string dump_file
        ){
    int dim = system.positions_.rows();
    const int offset_frame = 2 * dim * system.n_;  // size of each frame
    int offset = 0;
    bool should_dump = dump_file.size() > 0;

    if (load_file.size() > 0){
        load(system, load_file);
    }

    py::array_t<double> result = py::array_t<double>(
            run_steps * system.n_ * dim * 2
    );
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


/*
 * simulate vicsek-like systems without using the neighbour list
 */
template<class T>
py::array_t<double> simulate_no_nl(
        T system, int pre_steps, int run_steps, int jump,
        std::string load_file, std::string dump_file
        ){
    int dim = system.positions_.rows();
    const int offset_frame = 2 * dim * system.n_;  // size of each frame
    int offset = 0;
    bool should_dump = dump_file.size() > 0;

    if (load_file.size() > 0){
        load(system, load_file);
    }

    py::array_t<double> result = py::array_t<double>(
            run_steps * system.n_ * dim * 2
    );
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;

    for (int step = 0; step < pre_steps; step++){
        system.move_no_nl();
    }

    int cursor = 0;
    for (int step = 0; step < run_steps * jump; step++){
        system.move_no_nl();
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

/*
 * Perform a simulation and return a numpy.ndarray
 */
template<class T>
py::array_t<double> run(
        T system, int update_step, int run_steps, int jump
        ){
    int dim = system.positions_.rows();
    const int offset_frame = 2 * dim * system.n_;  // size of each frame
    int offset = 0;

    py::array_t<double> result = py::array_t<double>(run_steps * system.n_ * dim * 2);
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;
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
        }
    }
    result.resize({run_steps, system.n_, 2 * dim});

    return result;
}


/*
 * Perform a simulation and return a numpy.ndarray
 * without using neighbour list
 */
template<class T>
py::array_t<double> run_no_nl(T system, int run_steps, int jump){
    int dim = system.positions_.rows();
    const int offset_frame = 2 * dim * system.n_;  // size of each frame
    int offset = 0;
    py::array_t<double> result = py::array_t<double>(
            run_steps * system.n_ * dim * 2
            );
    auto buffer_result = result.request();
    auto *ptr_result = (double *) buffer_result.ptr;
    int cursor = 0;
    for (int step = 0; step < run_steps * jump; step++){
        system.move_no_nl();
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
        }
    }
    result.resize({run_steps, system.n_, 2 * dim});
    return result;
}


py::array_t<double> vicsek_3d_pbc(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek3DPBC system{n, r, eta, box, v0};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump, load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump, load_file, dump_file
        );
    }
}


py::array_t<double> continue_vicsek_3d_pbc(
        int n, double box, double eta, double v0, double r,
        int run_steps, Coord3D positions, Coord3D velocities,
        int jump=1, bool use_nl=false
        ){
    Vicsek3DPBC system{n, r, eta, box, v0};
    system.positions_ = positions;
    system.velocities_ = velocities;
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return run(system, update_step, run_steps, jump);
    } else {
        return run_no_nl(system, run_steps, jump);
    }
}


py::array_t<double> vicsek_3d(
        int n, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek3D system{n, r, eta, v0};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump, load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump, load_file, dump_file
        );
    }
}


py::array_t<double> continue_vicsek_3d(
        int n, double eta, double v0, double r,
        int run_steps, Coord3D positions, Coord3D velocities,
        int jump=1, bool use_nl=false
        ){
    Vicsek3D system{n, r, eta, v0};
    system.positions_ = positions;
    system.velocities_ = velocities;
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return run(system, update_step, run_steps, jump);
    } else {
        return run_no_nl(system, run_steps, jump);
    }
}


py::array_t<double> vicsek_2d(
        int n, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek2D system{n, r, eta, v0};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump, load_file, dump_file
        );
    }
}


py::array_t<double> continue_vicsek_2d(
        int n, double eta, double v0, double r,
        int run_steps, Coord2D positions, Coord2D velocities,
        int jump=1, bool use_nl=false
        ){
    Vicsek2D system{n, r, eta, v0};
    system.positions_ = positions;
    system.velocities_ = velocities;
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return run(system, update_step, run_steps, jump);
    } else {
        return run_no_nl(system, run_steps, jump);
    }
}


py::array_t<double> vicsek_3d_pbc_inertia(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek3DPBCInertia system{n, r, eta, box, v0, alpha};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}


py::array_t<double> vicsek_2d_pbc_inertia(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek2DPBCInertia system{n, r, eta, box, v0, alpha};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}



py::array_t<double> vicsek_3d_pbc_inertia_af(
        int n, double box, double eta, double v0, double r, double alpha,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek3DPBCInertiaAF system{n, r, eta, box, v0, alpha};
    if (use_nl) {
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}


py::array_t<double> attanasi2014pcb(
        int n, double eta, double v0, double b, double r,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Attanasi2014PCB system{n, r, eta, v0, b};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}


py::array_t<double> vicsek_2d_pbc(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek2DPBC system{n, r, eta, box, v0};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}


py::array_t<double> continue_vicsek_2d_pbc(
        int n, double box, double eta, double v0, double r,
        int run_steps, Coord2D positions, Coord2D velocities,
        int jump=1, bool use_nl=false
        ){
    Vicsek2DPBC system{n, r, eta, box, v0};
    system.positions_ = positions;
    system.velocities_ = velocities;
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return run(system, update_step, run_steps, jump);
    } else {
        return run_no_nl(system, run_steps, jump);
    }
}


py::array_t<double> vicsek_2d_pbc_vn(
        int n, double box, double eta, double v0, double r,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    Vicsek2DPBCVN system{n, r, eta, box, v0};
    int update_step = get_update_step(r, v0);
    if (use_nl){
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}


py::array_t<double> continue_vicsek_2d_pbc_vn(
        int n, double box, double eta, double v0, double r,
        int run_steps, Coord2D positions, Coord2D velocities,
        int jump=1, bool use_nl=false
        ){
    Vicsek2DPBCVN system{n, r, eta, box, v0};
    system.positions_ = positions;
    system.velocities_ = velocities;
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return run(system, update_step, run_steps, jump);
    } else {
        return run_no_nl(system, run_steps, jump);
    }
}


py::array_t<double> ism_3d(
        int n, double r, double v0,
        double T, double j, double m, double f,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    InertialSpin3D system{n, r, v0, T, j, m, f};
    if (use_nl){
    int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
        );
    }
}


py::array_t<double> ism_3d_pbc(
        int n, double box, double r, double v0,
        double T, double j, double m, double f,
        int pre_steps, int run_steps, int jump=1,
        std::string load_file="", std::string dump_file="",
        bool use_nl=false
        ){
    InertialSpin3DPBC system{n, box, r, v0, T, j, m, f};
    if (use_nl){
        int update_step = get_update_step(r, v0);
        return simulate(
            system, update_step, pre_steps, run_steps, jump,
            load_file, dump_file
            );
    } else {
        return simulate_no_nl(
            system, pre_steps, run_steps, jump,
            load_file, dump_file
            );
    }
}


ConnMat get_random_vnm_adjecancy_matrix(int d, int n, bool include_self){
    Graph graph;
    if (include_self){
        graph = random_vnm_graph_force_self(d, n);
    } else {
        graph = random_vnm_graph(d, n);
    }
    return graph.as_matrix();
}


ConnMat get_random_regular_adjecancy_matrix(int d, int n, bool include_self){
    Graph graph;
    if (include_self){
        graph = random_regular_graph(d - 1, n);
        for (auto node : graph.nodes_){
            graph.edges_.emplace(Index2D{node, node});
        }
    } else {
        graph = random_regular_graph(d, n);
    }
    return graph.as_matrix();
}


Index2D _unravel_index_2d(int index, Index2D shape){
    return unravel_index(index, shape);
}

Index3D _unravel_index_3d(int index, Index3D shape){
    return unravel_index(index, shape);
}


RotMat _get_rotation_matrix_from_001(double x, double y, double z){
    return get_rotation_matrix(x, y, z);
}

RotMat _get_rotation_matrix(Vec3D v1, Vec3D v2){
    return get_rotation_matrix(v1, v2);
}

RotMat _get_rotation_matrix_limited(Vec3D v1, Vec3D v2, double theta){
    return get_rotation_matrix(v1, v2, theta);
}

std::vector<int> _argsort(const std::vector<double> &values){
    return argsort(values);
}


Indices3D _product_3d(Indices a1, Indices& a2, Indices& a3){
    return product_3d(a1, a2, a3);
}

Indices2D _product_2d(Indices a1, Indices& a2){
    return product_2d(a1, a2);
}


PYBIND11_MODULE(cmodel, m){
    m.doc() = "Simulating different models to understand the collective behaviours";

    m.def(
        "_argsort", &_argsort, "c++ version of np.argsort", py::arg("array")
    );

    m.def(
        "_product_2d", &_product_2d, "2D Cartesian product of two arrays"
    );

    m.def(
        "_product_3d", &_product_3d, "3D Cartesian product of two arrays"
    );

    m.def(
        "_unravel_index_2d", &_unravel_index_2d,
        "c++ version of np.unravel_index for 2D array",
        py::arg("index"), py::arg("shape")
    );

    m.def(
        "_unravel_index_3d", &_unravel_index_3d,
        "c++ version of np.unravel_index for 3D array",
        py::arg("index"), py::arg("shape")
    );

    m.def(
        "_get_rotation_matrix_from_001", &_get_rotation_matrix_from_001,
        "get the rotation matrix from 001 to x, y, z",
        py::arg("x"), py::arg("y"), py::arg("z")
    );

    m.def(
        "_get_rotation_matrix", &_get_rotation_matrix,
        "get the rotation matrix from v1 to v2",
        py::arg("v1"), py::arg("v2")
    );

    m.def(
        "_get_rotation_matrix_limited", &_get_rotation_matrix_limited,
        "get the rotation matrix from v1 to v2, with the maximum rotation value of theta",
        py::arg("v1"), py::arg("v2"), py::arg("theta")
    );

    m.def(
        "vicsek_3d_pbc", &vicsek_3d_pbc, docstring_vicsek_3d_pbc,
        py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "continue_vicsek_3d_pbc",
        &continue_vicsek_3d_pbc, docstring_vicsek_3d_pbc,
        py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("run_steps"), py::arg("positions"), py::arg("velocities"),
        py::arg("jump")=1, py::arg("use_nl")=false
    );

    m.def(
        "vicsek_3d", &vicsek_3d, docstring_vicsek_3d,
        py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "continue_vicsek_3d",
        &continue_vicsek_3d, docstring_vicsek_3d,
        py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("run_steps"), py::arg("positions"), py::arg("velocities"),
        py::arg("jump")=1, py::arg("use_nl")=false
    );

    m.def(
        "vicsek_2d", &vicsek_2d, docstring_vicsek_2d,
        py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "continue_vicsek_2d",
        &continue_vicsek_2d, docstring_vicsek_2d,
        py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("run_steps"), py::arg("positions"), py::arg("velocities"),
        py::arg("jump")=1, py::arg("use_nl")=false
    );

    m.def(
        "vicsek_3d_pbc_inertia",
        &vicsek_3d_pbc_inertia, docstring_vicsek_3d_pbc_inertia,
        py::arg("n"), py::arg("box"), py::arg("eta"),
        py::arg("v0"), py::arg("r"), py::arg("alpha"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "vicsek_2d_pbc_inertia",
        &vicsek_2d_pbc_inertia, docstring_vicsek_2d_pbc_inertia,
        py::arg("n"), py::arg("box"), py::arg("eta"),
        py::arg("v0"), py::arg("r"), py::arg("alpha"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "vicsek_3d_pbc_inertia_af",
        &vicsek_3d_pbc_inertia_af, docstring_vicsek_3d_pbc_inertia_af,
        py::arg("n"), py::arg("box"), py::arg("eta"),
        py::arg("v0"), py::arg("r"), py::arg("alpha"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "attanasi2014pcb", &attanasi2014pcb, docstring_attanasi2014pcb,
        py::arg("n"), py::arg("eta"), py::arg("v0"), py::arg("b"), py::arg("r"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "vicsek_2d_pbc", &vicsek_2d_pbc, docstring_vicsek_2d_pbc,
        py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file")="",
        py::arg("use_nl")=false
    );

    m.def(
        "continue_vicsek_2d_pbc", &continue_vicsek_2d_pbc, docstring_vicsek_2d_pbc,
        py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("run_steps"), py::arg("positions"), py::arg("velocities"),
        py::arg("jump")=1, py::arg("use_nl")=false
    );

    m.def(
        "vicsek_2d_pbc_vn", &vicsek_2d_pbc_vn, docstring_vicsek_2d_pbc_vn,
        py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file=")="",
        py::arg("use_nl")=false
    );

    m.def(
        "continue_vicsek_2d_pbc_vn",
        &continue_vicsek_2d_pbc_vn, docstring_vicsek_2d_pbc_vn,
        py::arg("n"), py::arg("box"), py::arg("eta"), py::arg("v0"), py::arg("r"),
        py::arg("run_steps"), py::arg("positions"), py::arg("velocities"),
        py::arg("jump")=1, py::arg("use_nl")=false
    );

    m.def(
        "ism_3d", &ism_3d, docstring_ism_3d,
        py::arg("n"), py::arg("r"), py::arg("v0"),
        py::arg("T"), py::arg("j"), py::arg("m"), py::arg("f"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file=")="",
        py::arg("use_nl")=false
    );

    m.def(
        "ism_3d_pbc", &ism_3d_pbc, docstring_ism_3d_pbc,
        py::arg("n"), py::arg("box"), py::arg("r"), py::arg("v0"),
        py::arg("T"), py::arg("j"), py::arg("m"), py::arg("f"),
        py::arg("pre_steps"), py::arg("run_steps"),
        py::arg("jump")=1, py::arg("load_file")="", py::arg("dump_file=")="",
        py::arg("use_nl")=false
    );

    m.def(
        "get_random_vnm_adjecancy_matrix", &get_random_vnm_adjecancy_matrix,
        R"---(
        Generate a random adjecancy matrix for the vertorial network

        Args:
            k (int): the number of neighbours of each node.
            n (int): the total number of nodes.
            include_self (bool): if true, all nodes are forced to
                have links to themselves.
        Return:
            np.ndarray: the adjecancy matrix of the graph.
        )---",
        py::arg("k"), py::arg("n"), py::arg("include_self")=false
    );

    m.def(
        "get_random_regular_adjecancy_matrix", &get_random_regular_adjecancy_matrix,
        R"---(
        Generate a random adjecancy matrix for a random graph

        Args:
            k (int): the number of neighbours of each node.
            n (int): the total number of nodes.
            include_self (bool): if true, all nodes are forced to
                have links to themselves. If false, no self-loop
                is allowed. (The graph is a simple graph.)
        Return:
            np.ndarray: the adjecancy matrix of the graph.
        )---",
        py::arg("k"), py::arg("n"), py::arg("include_self")=false
    );

    py::class_<Vicsek3D>(m, "Vicsek3D")
        .def(
            py::init<
                int, double, double, double
            >(),
            pybind11::arg("n"),
            pybind11::arg("r"),
            pybind11::arg("eta"),
            pybind11::arg("v0")
        )
        .def("move", &Vicsek3D::move, py::arg("rebuild")=true)
        .def("evolve", &Vicsek3D::evolve, py::arg("steps"), py::arg("rebuild")=true)
        .def("get_polarisation", &Vicsek3D::get_polarisation)
        .def("get_velocities", &Vicsek3D::get_velocities)
        .def("get_positions", &Vicsek3D::get_positions)
        .def("load_velocities", &Vicsek3D::load_velocities)
        .def_readonly("dim", &Vicsek3D::dim_)
        .def_property("positions",
                &Vicsek3D::get_positions, &Vicsek3D::load_positions,
                py::return_value_policy::copy
        )
        .def_property("velocities",
            &Vicsek3D::get_velocities, &Vicsek3D::load_velocities,
            py::return_value_policy::copy
        );

    py::class_<AVS3D>(m, "AVS3D")
        .def(
            py::init< int, double >(),
            pybind11::arg("n"),
            pybind11::arg("eta")
        )
        .def("move", &AVS3D::add_noise)
        .def("get_polarisation", &AVS3D::get_polarisation)
        .def_readonly("dim", &AVS3D::dim_)
        .def_property("orientations",
                &AVS3D::get_orientations, &AVS3D::load_orientations,
                py::return_value_policy::copy
        );

    py::class_<AVS2D>(m, "AVS2D")
        .def(
            py::init< int, double >(),
            pybind11::arg("n"),
            pybind11::arg("eta")
        )
        .def("move", &AVS2D::add_noise)
        .def("get_polarisation", &AVS2D::get_polarisation)
        .def_readonly("dim", &AVS2D::dim_)
        .def_property("orientations",
                &AVS2D::get_orientations, &AVS2D::load_orientations,
                py::return_value_policy::copy
        );

    py::class_<Vicsek3DPBC>(m, "Vicsek3DPBC")
        .def(
            py::init<
                int, double, double, double, double
            >(),
            pybind11::arg("n"),
            pybind11::arg("r"),
            pybind11::arg("eta"),
            pybind11::arg("box"),
            pybind11::arg("v0")
        )
        .def("move", &Vicsek3DPBC::move, py::arg("rebuild")=true)
        .def("evolve", &Vicsek3DPBC::evolve, py::arg("steps"), py::arg("rebuild")=true)
        .def("get_polarisation", &Vicsek3DPBC::get_polarisation)
        .def("get_velocities", &Vicsek3DPBC::get_velocities)
        .def("get_positions", &Vicsek3DPBC::get_positions)
        .def("load_velocities", &Vicsek3DPBC::load_velocities)
        .def_property("positions",
                &Vicsek3DPBC::get_positions, &Vicsek3DPBC::load_positions,
                py::return_value_policy::copy
        )
        .def_property("velocities",
                &Vicsek3DPBC::get_velocities, &Vicsek3DPBC::load_velocities,
                py::return_value_policy::copy
        )
        .def_readonly("dim", &Vicsek3DPBC::dim_);

    py::class_<Vicsek3DPBCInertia>(m, "Vicsek3DPBCInertia")
        .def(
            py::init<
                int, double, double, double, double, double
            >(),
            pybind11::arg("n"),
            pybind11::arg("r"),
            pybind11::arg("eta"),
            pybind11::arg("box"),
            pybind11::arg("v0"),
            pybind11::arg("alpha")
        )
        .def(
            "move", &Vicsek3DPBCInertia::move,
            R"---(
            Advance the dynamic, move one Monte-Carlo step

            Args:
                rebuild (bool): if true, the neighobur will be\
                    used to accelerate the calculation.
            )---",
            py::arg("rebuild")=true
        )
        .def(
            "move_no_nl", &Vicsek3DPBCInertia::move_no_nl,
            R"---(
            Advance the dynamic, move one Monte-Carlo step
            )---"
        )
       .def(
            "evolve", &Vicsek3DPBCInertia::evolve,
            R"---(
            Advance the dynamic, move multiple Monte-Carlo steps

            Args:
                steps (int): the number of Monte-Carlo steps to move.
                rebuild (bool): if true, the neighobur will be\
                    used to accelerate the calculation.

            )---",
            py::arg("steps"), py::arg("rebuild")=true
        )
        .def(
            "get_polarisation", &Vicsek3DPBCInertia::get_polarisation,
            "Calculate the polarisation of current state of the system",
            py::return_value_policy::copy
        )
        .def(
            "get_velocities", &Vicsek3DPBCInertia::get_velocities,
            "Retrieve the current velocities as numpy array of the system"
        )
        .def(
            "get_positions", &Vicsek3DPBCInertia::get_positions,
            "Retrieve the current positions as numpy array of the system"
        )
        .def(
            "load_velocities", &Vicsek3DPBCInertia::load_velocities,
            "Set the current velocities of the system",
            py::arg("velocities")
        )
        .def_property("positions",
            &Vicsek3DPBCInertia::get_positions, &Vicsek3DPBCInertia::load_positions,
            py::return_value_policy::copy
        )
        .def_property("velocities",
            &Vicsek3DPBCInertia::get_velocities, &Vicsek3DPBCInertia::load_velocities,
            py::return_value_policy::copy
        )
        .def_readonly("dim", &Vicsek3DPBCInertia::dim_);

    py::class_<Vicsek2DPBCInertia>(m, "Vicsek2DPBCInertia")
        .def(
            py::init<
                int, double, double, double, double, double
            >(),
            pybind11::arg("n"),
            pybind11::arg("r"),
            pybind11::arg("eta"),
            pybind11::arg("box"),
            pybind11::arg("v0"),
            pybind11::arg("alpha")
        )
        .def(
            "move", &Vicsek2DPBCInertia::move,
            R"---(
            Advance the dynamic, move one Monte-Carlo step

            Args:
                rebuild (bool): if true, the neighobur will be\
                    used to accelerate the calculation.
            )---",
            py::arg("rebuild")=true
        )
        .def(
            "move_no_nl", &Vicsek2DPBCInertia::move_no_nl,
            R"---(
            Advance the dynamic, move one Monte-Carlo step
            )---"
        )
        .def(
            "get_polarisation", &Vicsek2DPBCInertia::get_polarisation,
            "Calculate the polarisation of current state of the system",
            py::return_value_policy::copy
        )
        .def(
            "get_velocities", &Vicsek2DPBCInertia::get_velocities,
            "Retrieve the current velocities as numpy array of the system"
        )
        .def(
            "get_positions", &Vicsek2DPBCInertia::get_positions,
            "Retrieve the current positions as numpy array of the system"
        )
        .def(
            "load_velocities", &Vicsek2DPBCInertia::load_velocities,
            "Set the current velocities of the system",
            py::arg("velocities")
        )
        .def_property("positions",
            &Vicsek2DPBCInertia::get_positions, &Vicsek2DPBCInertia::load_positions,
            py::return_value_policy::copy
        )
        .def_property("velocities",
            &Vicsek2DPBCInertia::get_velocities, &Vicsek2DPBCInertia::load_velocities,
            py::return_value_policy::copy
        )
        .def_readonly("dim", &Vicsek2DPBCInertia::dim_);

    py::class_<Network3D>(m, "Network3D")
        .def(
            py::init<int, int, double>(),
            pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("eta")
        )
        .def(
            "move", &Network3D::move,
            R"---(
            Advance the dynamic, move one Monte-Carlo step

            Args:
                new_graph (bool): if true, a new graph will be
                    generated.
            )---",
            py::arg("new_graph")=true
        )
        .def(
            "evolve", &Network3D::evolve,
            R"---(
            Advance the dynamic, move multiple Monte-Carlo steps

            Args:
                steps (int): the number of Monte-Carlo steps to move.\
                dynamic (int): 0 = quenched dynamic, no graph update.\
                    1 = anneleda dynamic, graph update every step.

            )---",
            py::arg("steps"), py::arg("rebuild")=true
        )
        .def(
            "get_polarisation", &Network3D::get_polarisation,
            "Calculate the polarisation of current state of the system",
            py::return_value_policy::copy
        )
        .def(
            "get_orientations", &Network3D::get_orientations,
            "Retrieve the current orientations as numpy array of the system"
        )
        .def_property("orientations",
            &Network3D::get_orientations, &Network3D::load_orientations,
            py::return_value_policy::copy
        )
        .def_property("adj_mat",
            &Network3D::get_adj_mat, &Network3D::set_adj_mat,
            py::return_value_policy::copy
        )
        .def_readonly("dim", &Network3D::dim_);

    py::class_<Network3DRG>(m, "Network3DRG")
        .def(
            py::init<int, int, double>(),
            pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("eta")
        )
        .def(
            "move", &Network3DRG::move,
            R"---(
            Advance the dynamic, move one Monte-Carlo step

            Args:
                new_graph (bool): if true, a new regular graph will be
                    generated.
            )---",
            py::arg("new_graph")=true
        )
        .def(
            "evolve", &Network3DRG::evolve,
            R"---(
            Advance the dynamic, move multiple Monte-Carlo steps

            Args:
                steps (int): the number of Monte-Carlo steps to move.\
                dynamic (int): 0 = quenched dynamic, no graph update.\
                    1 = anneleda dynamic, graph update every step.

            )---",
            py::arg("steps"), py::arg("rebuild")=true
        )
        .def(
            "get_polarisation", &Network3DRG::get_polarisation,
            "Calculate the polarisation of current state of the system",
            py::return_value_policy::copy
        )
        .def(
            "get_orientations", &Network3DRG::get_orientations,
            "Retrieve the current velocities as numpy array of the system"
        )
        .def_property("velocities",
            &Network3DRG::get_orientations, &Network3DRG::load_orientations,
            py::return_value_policy::copy
        )
        .def_property("adj_mat",
            &Network3DRG::get_adj_mat, &Network3DRG::set_adj_mat,
            py::return_value_policy::copy
        )
        .def_readonly("dim", &Network3DRG::dim_);

    py::class_<Voter>(m, "Voter")
        .def(
            py::init<int, int, double>(),
            pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("eta")
        )
        .def(
            "move", &Voter::move,
            R"---(
            Advance the dynamic, move one Monte-Carlo step

            Args:
                new_graph (bool): if true, a new regular graph will be
                    generated.
            )---",
            py::arg("new_graph")=true
        )
        .def(
            "evolve", &Voter::evolve,
            R"---(
            Advance the dynamic, move multiple Monte-Carlo steps

            Args:
                steps (int): the number of Monte-Carlo steps to move.\
                dynamic (int): 0 = quenched dynamic, no graph update.\
                    1 = anneleda dynamic, graph update every step.

            )---",
            py::arg("steps"), py::arg("rebuild")=true
        )
        .def(
            "get_polarisation", &Voter::get_polarisation,
            "Calculate the get_polarisation per spin of current state of the system",
            py::return_value_policy::copy
        )
        .def_property("spins",
            &Voter::get_spins, &Voter::load_spins,
            py::return_value_policy::copy
        )
        .def_property("adj_mat",
            &Voter::get_adj_mat, &Voter::set_adj_mat,
            py::return_value_policy::copy
        )
        .def_readonly("dim", &Voter::dim_);

    py::class_<Couzin3D>(m, "Couzin3D")
        .def(
            py::init<
                int, double, double, double,
                double, double, double, double, double
            >(),
            pybind11::arg("n"),
            pybind11::arg("rr"),
            pybind11::arg("ro"),
            pybind11::arg("ra"),
            pybind11::arg("perception"),
            pybind11::arg("noise"),
            pybind11::arg("speed"),
            pybind11::arg("turn_rate"),
            pybind11::arg("dt")
        )
        .def("move", &Couzin3D::move, py::arg("rebuild")=true)
        .def("evolve", &Couzin3D::evolve, py::arg("steps"), py::arg("rebuild")=true)
        .def("get_polarisation", &Couzin3D::get_polarisation)
        .def("get_mill", &Couzin3D::get_mill)
        .def("get_velocities", &Couzin3D::get_velocities)
        .def("get_positions", &Couzin3D::get_positions)
        .def("load_velocities", &Couzin3D::load_velocities)
        .def_readonly("dim", &Couzin3D::dim_)
        .def_property("positions",
                &Couzin3D::get_positions, &Couzin3D::load_positions,
                py::return_value_policy::copy
        )
        .def_property("velocities",
            &Couzin3D::get_velocities, &Couzin3D::load_velocities,
            py::return_value_policy::copy
        );

    py::class_<CouzinTank3D>(m, "CouzinTank3D")
        .def(
            py::init<
                int, double, double, double,
                double, double, double, double, double,
                double, double, double, bool, bool,
                double
            >(),
            pybind11::arg("n"),
            pybind11::arg("rr"),
            pybind11::arg("ro"),
            pybind11::arg("ra"),
            pybind11::arg("perception"),
            pybind11::arg("noise"),
            pybind11::arg("speed"),
            pybind11::arg("turn_rate"),
            pybind11::arg("dt"),
            pybind11::arg("c"),
            pybind11::arg("h"),
            pybind11::arg("kw")=0,
            pybind11::arg("align_cap")=false,
            pybind11::arg("align_base")=false,
            pybind11::arg("g")=0
        )
        .def("move", &CouzinTank3D::move_in_tank, py::arg("rebuild")=true)
        .def("evolve", &CouzinTank3D::evolve_in_tank, py::arg("steps"), py::arg("rebuild")=true)
        .def("get_polarisation", &CouzinTank3D::get_polarisation)
        .def("get_velocities", &CouzinTank3D::get_velocities)
        .def("get_positions", &CouzinTank3D::get_positions)
        .def("load_velocities", &CouzinTank3D::load_velocities)
        .def_readonly("dim", &CouzinTank3D::dim_)
        .def_property("positions",
            &CouzinTank3D::get_positions, &CouzinTank3D::load_positions,
            py::return_value_policy::copy
        )
        .def_property("velocities",
            &CouzinTank3D::get_velocities, &CouzinTank3D::load_velocities,
            py::return_value_policy::copy
        );

    py::class_<Tank3D>(m, "Tank3D")
        .def(
            py::init< double, double >(),
            pybind11::arg("c"),
            pybind11::arg("z_max")
        )
        .def(
            "get_random_positions", &Tank3D::get_random_positions,
            R"---(
            Generate random points inside the boundary

            Args:
                n (int): the number of particles inside the fish bowl
            )---",
            py::arg("n")
        )
        .def(
            "get_random_positions_reject", &Tank3D::get_random_positions_reject,
            R"---(
            Generate random points inside the boundary using rejection method

            Args:
                n (int): the number of particles inside the fish bowl
            )---",
            py::arg("n")
        )
        .def(
            "project", &Tank3D::project,
            R"---(
            Generate random points inside the boundary.

            Args:
                xyz (numpy.ndarray): the coordinates, shape (3, n).

            Return:
               numpy.ndarray: the points projected on the bottom of the bowl.
            )---",
            py::arg("xyz")
        )
        .def(
            "is_inside", &Tank3D::is_inside_single,
            R"---(
            Check if a point is inside the tank.

            Args:
                xyz (numpy.ndarray): the coordinates, shape (3,).

            Return:
               bool: if True the particle is in the tank
            )---",
            py::arg("xyz")
        )
        .def(
            "get_dist_to_tank_single", &Tank3D::get_dist_to_tank_single,
            R"---(
            Check if a point is inside the tank.

            Args:
                xyz (numpy.ndarray): the coordinates, shape (3,).

            Return:
               float: the distance to the tank
            )---",
            py::arg("xyz")
        )
        .def(
            "project_single", &Tank3D::project_single,
            R"---(
            Generate random points inside the boundary.

            Args:
                x (float): the x coordinate
                y (float): the y coordinate
                z (float): the z coordinate

            Return:
               numpy.ndarray: the points projected on the bottom of the bowl.
            )---",
            py::arg("x"), py::arg("y"), py::arg("z")
        );

    py::class_<FishMCMC>(m, "FishMCMC")
        .def(
            py::init<
                size_t, double, double, double, double
            >(),
            pybind11::arg("n"),
            pybind11::arg("beta"),
            pybind11::arg("si"),
            pybind11::arg("sh"),
            pybind11::arg("dx")
        )
        .def("sweep", &FishMCMC::sweep)
        .def("evolve", &FishMCMC::evolve, py::arg("steps"))
        .def("get_energy", &FishMCMC::get_total_energy)
        .def("get_positions", &FishMCMC::get_positions)
        .def("load_positions", &FishMCMC::load_positions)
        .def_property("positions",
                &FishMCMC::get_positions, &FishMCMC::load_positions,
                py::return_value_policy::copy
        );
}
