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
        double T_; // temperature
        double J_; // coupling constant 
        double v0_sq_inv_;
        VerletList3D verlet_list_;

    public:
        int n_;
        Coord3D spins_;
        Coord3D positions_;
        Coord3D velocities_;
        Conn connections_;
        double dt_;
        InertialSpin3D(
            int n, double r, double eta, double v0, double j, double m, double f
        );
        void move(bool rebuid);
        void add_alignment();
        void add_noise();
        void update_velocity();
        void update_spin();
};


class InertialSpin3DPBC : public InertialSpin3D {
    protected:
        double box_;
        CellList3D cell_list_;
        inline void fix_positions(){
            positions_.array() -= box_ * (positions_.array() / box_).floor();
        }

    public:
        InertialSpin3DPBC(
            int n, double box, double r, double eta, double v0,
            double j, double m, double f
        );
        void move(bool rebuid);
};

#endif
