#ifndef VICSEK
#define VICSEK

#include <string>
#include <iostream>
#include <random>
#include <cmath>
#include <regex>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include "neighbour_list.h"

using namespace std;

using Property = Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)
using PropertyInt = Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)

using Coord3D = Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor>;  // (3, n)
using Vec3D = Eigen::Array3d;  // (3, 1)
using RotMat = Eigen::Matrix3d;

using Coord2D = Eigen::Array<double, 2, Eigen::Dynamic, Eigen::RowMajor>;  // (2, n)
using Vec2D = Eigen::Array2d;  // (2, )


class Vicsek3D{
    protected:
        double noise_;
        double speed_;
        int n_;
        Conn connections_;
        void align();
        void add_noise();
        void rotate_noise(Coord3D& noise_xyz);
        void rotate_noise_fast(Coord3D& noise_xyz);
        void rotate_noise_xyz(Coord3D& noise_xyz);
        VerletList3D verlet_list_;

    public:
        Coord3D positions_;
        Coord3D velocities_;

        Vicsek3D(int n, double r, double eta, double v0);
        void move(bool rebuild);
        void dump(string filename);
        void load(string filename);
};


class Vicsek3DPBC : public Vicsek3D{
    protected:
        double box_;
        CellList3D cell_list_;
        inline void fix_positions(){
            positions_ -= box_ * (positions_ / box_).floor();
        }

    public:
        Vicsek3DPBC(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
};


class Attanasi2014PCB : public Vicsek3D{
    protected:
        void harmonic_align();
        double beta_;
    public:
        Attanasi2014PCB(int n, double r, double eta, double v0, double beta);
        void move(bool rebuild);
};


class Vicsek3DPBCVN : public Vicsek3DPBC{
    private:
        void noisy_align();
    public:
        Vicsek3DPBCVN(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
};


class Vicsek3DPBCInertia : public Vicsek3DPBC{
    public:
        double alpha_;
        Coord3D old_velocities_;
        Vicsek3DPBCInertia(int n, double r, double eta, double box, double v0, double alpha);
        void move(bool rebuild);
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
};


class Vicsek2DPBC{
    protected:
        double noise;
        double speed;
        double box;
        int n;
        ConnMat conn_mat;

        void align();
        void add_noise();
        inline void fix_positions(){positions -= box * (positions / box).floor();}

    public:
        Coord2D positions;
        Coord2D velocities;
        CellList2D cell_list;

        Vicsek2DPBC(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void dump(string filename);
        void load(string filename);
        Vec2D get_shift(Vec2D p1, Vec2D p2);
};


class Vicsek2DPBCVN : public Vicsek2DPBC{
    private:
        void noisy_align();
    public:
        Vicsek2DPBCVN(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
};


class Vicsek2DPBCVNCO : public Vicsek2DPBC{
    private:
        void update_velocity();
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
};



Coord3D xyz_to_sphere(Coord3D& xyz);
Coord3D sphere_to_xyz(Coord3D& sphere);

void normalise(Coord3D& xyz);
void normalise(Coord3D& xyz, double spd);

void normalise(Coord2D& xy);
void normalise(Coord2D& xy, double spd);

#endif
