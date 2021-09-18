#ifndef COUZIN
#define COUZIN

#include "core.hpp"
#include "neighbour_list.hpp"


/*
 * Implementing the "couzin model"
 * see IDC JK RJ GDR NRF, J. theor. Biol., 2002 (10.1006/yjtbi.3065)
 */
class Couzin3D : public AVS3D {
    protected:
        VerletList<Coord3D> verlet_list_;
        double r_repel_,
               r_align_,
               r_attract_,
               a_percept_,
               v_turn_,
               dt_;
        Conn conn_repel_;
        Conn conn_align_;
        Conn conn_attract_;
        bool is_visible(int i, int j); // agent i can see agent j

    public:
        Coord3D positions_;
        Coord3D velocities_real_;  // velocities_ represent the target moving direction
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
};


class Couzin3DPBC : public Couzin3D {
    protected:
        double box_;
        CellList3D cell_list;
};


#endif
