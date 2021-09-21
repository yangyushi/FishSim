#ifndef VICSEK
#define VICSEK

#include "core.hpp"
#include "neighbour_list.hpp"


class Vicsek3D : public AVS3D{
    protected:
        double speed_;
        double rc_;
        inline void update_velocity(){ velocities_ = orientations_ * speed_; }
        Conn connections_;
        VerletList<Coord3D> verlet_list_;

    public:
        Coord3D positions_;
        Coord3D velocities_;
        Vicsek3D(int n, double r, double eta, double v0);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
        inline void evolve(int steps, bool rebuild){
            for (int s=0; s < steps; s++){
                this->move(rebuild);
            }
        };
        // for interact with python
        inline Coord3D& get_positions() { return positions_; };
        void load_positions(Coord3D);
        inline Coord3D& get_velocities() { return velocities_; }
        void load_velocities(Coord3D v);
};


class Vicsek3DPBC : public Vicsek3D{
    protected:
        double box_;
        CellList3D cell_list_;
        inline void fix_positions(){
            positions_.array() -= box_ * (positions_.array() / box_).floor();
        }

    public:
        Vicsek3DPBC(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


class Attanasi2014PCB : public Vicsek3D{
    protected:
        void harmonic_align();
        double beta_;
    public:
        Attanasi2014PCB(int n, double r, double eta, double v0, double beta);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


class Vicsek3DPBCVN : public Vicsek3DPBC{
    public:
        Vicsek3DPBCVN(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


class Vicsek3DPBCInertia : public Vicsek3DPBC{
    public:
        Coord3D old_orientations_;
        double alpha_;
        Vicsek3DPBCInertia(int n, double r, double eta, double box, double v0, double alpha);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


/*
 * Anti-ferromagnetic version of Vicsek model, J = -1 not 1
 */
class Vicsek3DPBCInertiaAF : public Vicsek3DPBCInertia{
    private:
        void align_af();
    public:
        Vicsek3DPBCInertiaAF(int n, double r, double eta, double box, double v0, double alpha);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


class Vicsek2D : public AVS2D{
    protected:
        double rc_;
        double speed_;
        Conn connections_;
        VerletList<Coord2D> verlet_list_;
        void update_velocity();

    public:
        Coord2D positions_;
        Coord2D velocities_;
        Vicsek2D(int n, double r, double eta, double v0);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list

        // for interact with python
        inline Coord2D& get_positions() { return positions_; };
        inline void load_positions(Coord2D p) {positions_ = p;} ;
        inline Coord2D& get_velocities() { return velocities_; }
        void load_velocities(Coord2D v);
};


class Vicsek2DPBC : public Vicsek2D{
    protected:
        double box_;
        Conn connections_;
        inline void fix_positions() {
            positions_.array() -= box_ * (positions_.array() / box_).floor();
        }

    public:
        CellList2D cell_list_;
        Vicsek2DPBC(int n, double r, double eta, double box, double v0);
        Vec2D get_shift(Vec2D p1, Vec2D p2);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


class Vicsek2DPBCVN : public Vicsek2DPBC{
    public:
        Vicsek2DPBCVN(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


class Vicsek2DPBCVNCO : public Vicsek2DPBC{
    private:
        void apply_interaction();
        double alpha_;
        double beta_;
        double ra_;
        double re_;
        double rc_;
        double c_;
    public:
        Vicsek2DPBCVNCO(
                int n, double r, double eta, double box, double v0,  // vicsek model
                double a, double b, double ra, double re, double rc  // for cohesion
                );
        // default cohesion parameter
        Vicsek2DPBCVNCO(int n, double r, double eta, double box, double v0);  
        void move(bool rebuild);
        void move_no_nl();  // without neighbour list
};


/*
 *  * please refer to 10.1007/s10955-014-1119-3 for the origin of the model
 */
class InertialSpin3D{
    protected:
        double mass_;
        double imass_;
        double friction_;
        double T_; // temperature
        double J_; // coupling constant 
        double rc_;
        double v0_sq_inv_;
        VerletList<Coord3D> verlet_list_;

    public:
        int n_;
        Coord3D spins_;
        Coord3D positions_;
        Coord3D velocities_;
        Conn connections_;
        double dt_;
        double speed_;
        InertialSpin3D(
            int n, double r, double v0, double T, double j, double m, double f
        );
        void move(bool rebuid);
        void move_no_nl();  // without neighbour list
        void add_alignment();
        void add_noise();
        void update_velocity_half();
        void update_velocity_full();
        void update_spin();
};


/*
 * use topological distance to find neighbours instead of metric
 */
class InertialSpin3DTP : public InertialSpin3D {
    private:
        int nc_;  // nc_ nearest neighbours were considered
    public:
        InertialSpin3DTP(
            int n, int nc, double v0,
            double T, double j, double m, double f
        );
        void move();  // without neighbour list
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
            int n, double box, double r, double v0,
            double T, double j, double m, double f
        );
        void move(bool rebuid);
        void move_no_nl();  // without neighbour list
};

#endif
