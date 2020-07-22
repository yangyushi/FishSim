#ifndef SPIN
#define SPIN

#include <vector>
#include <array>
#include <map>
#include <Eigen/Dense>
#include "vicsek.h"
#include "neighbour_list.h"

/*
 *  * please refer to 10.1007/s10955-014-1119-3 for the origin of the model
 *   */
class InertialSpin3D{
    private:
        double speed_;
        double mass_;
        double imass_;
        double friction_;
        double dt_;
        double T_; // temperature
        double J_; // coupling constant 
        double v0_sq_inv_;
        Conn connections_;
        VerletList3D verlet_list_;
        void add_alignment();
        void add_noise();
        void update_velocity();
        void update_spin();

    public:
        int n_;
        Coord3D spins_;
        Coord3D positions_;
        Coord3D velocities_;
        InertialSpin3D(
            int n, double r, double eta, double v0, double j, double m, double f
        );
        void move(bool rebuid);
};

#endif
