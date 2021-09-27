#ifndef COUZIN
#define COUZIN

#include "core.hpp"
#include "boundary.hpp"
#include "neighbour_list.hpp"


/*
 * Implementing the "couzin model"
 * see IDC JK RJ GDR NRF, J. theor. Biol., 2002 (10.1006/yjtbi.3065)
 */
class Couzin3D : public AVS3D {
    protected:
        VerletList<Coord3D> verlet_list_;
        double speed_,
               r_repel_,
               r_align_,
               r_attract_,
               a_percept_,
               v_turn_,
               dt_;
        Conn conn_repel_;
        Conn conn_align_;
        Conn conn_attract_;
        bool is_visible(int i, int j); // agent i can see agent j
        inline void update_velocity(){
            velocities_ = orientations_ * speed_;
        };

    public:
        Coord3D positions_;
        Coord3D velocities_;  // velocities_ represent the target moving direction
        Couzin3D(
            int n,
            double rr, double ro, double ra,  // repel, align, attract ranges
            double ap,  // angle of perception
            double eta, // noise
            double v0,  // speed
            double vt,  // turning rate
            double dt   // delta time
        );
        void move(bool rebuild);
        void evolve(int steps, bool rebuild);
        inline Coord3D& get_positions() { return positions_; };
        void load_positions(Coord3D);
        inline Coord3D& get_velocities() { return velocities_; };
        inline void load_velocities(Coord3D v) { velocities_ = v ;};
};


class CouzinTank3D : public Couzin3D {
    Tank3D tank_;

    public:
        CouzinTank3D(
            int n,
            double rr, double ro, double ra,  // repel, align, attract ranges
            double ap,  // angle of perception
            double eta, // noise
            double v0,  // speed
            double vt,  // turning rate
            double dt,  // delta time
            double c,   // shape parameter for tank
            double h,   // height of tank
            double kw,   // strength of wall interaction
            bool align   // true: align with wall; false: reflect from wall
        );
        void move_in_tank(bool rebuild);
        inline void evolve_in_tank(int steps, bool rebuild){
            for (int i = 0; i < steps; i++){
                this->move_in_tank(rebuild);
            }
        };
};

#endif
