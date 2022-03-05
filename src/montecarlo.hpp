#ifndef MONTECARLO
#define MONTECARLO

#include "core.hpp"
#include "boundary.hpp"


class FishMCMC{
    private:
        Tank3D tank_{0.74, 0.35};
        Coord3D positions_;
        std::array<double, 3> r_holes {0.096, 0.213, 0.327};
        //r1_ = 0.096, r2_=0.213, r3_=0.327;
        // lognormal fit of gr
        double h_=4.636158, c_=0.042691, w_=0.113731, a_=1.302245;
        double get_dE(Vec3D p_new, Vec3D p_old, size_t i);
        double get_energy(Vec3D pos, size_t i);

    public:
        /*
         *  sh -- strength of holes
         *  si -- strength of pairwise interaction
         */
        FishMCMC(size_t n, double beta, double si, double sh, double dx);
        size_t n_;
        double beta_, si_, sh_, dx_;
        void sweep();
        void evolve(size_t steps);
        inline Coord3D get_positions() { return positions_; };
        inline void load_positions(Coord3D pos) {
            this->positions_ = pos;
        };
        double get_total_energy();
};

#endif
