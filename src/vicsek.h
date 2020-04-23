#ifndef VICSEK
#define VICSEK

#include <string>
#include <iostream>
#include <cmath>
#include <regex>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include "neighbour_list.h"

using namespace std;

using Property = Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)

using Coord3D = Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor>;  // (3, n)
using Vec3D = Eigen::Array3d;  // (3, )
using RotMat = Eigen::Matrix3d;

using Coord2D = Eigen::Array<double, 2, Eigen::Dynamic, Eigen::RowMajor>;  // (2, n)
using Vec2D = Eigen::Array2d;  // (2, )


class Vicsek3DPBC{
    private:
        double noise;
        double speed;
        double box;
        int n;
        DistMat dist_mat;

        void align();
        void add_noise();
        void rotate_noise(Coord3D& noise_xyz);
        inline void fix_positions(){positions -= box * (positions / box).floor();}

    public:
        Coord3D positions;
        Coord3D velocities;
        CellList3D cell_list;

        Vicsek3DPBC(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void dump(string filename);
        void load(string filename);
};


class Vicsek2DPBC{
    private:
        double noise;
        double speed;
        double box;
        int n;
        DistMat dist_mat;

        void align();
        void add_noise();
        void rotate_noise(Coord2D& noise_xyz);
        inline void fix_positions(){positions -= box * (positions / box).floor();}

    public:
        Coord2D positions;
        Coord2D velocities;
        CellList2D cell_list;

        Vicsek2DPBC(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void dump(string filename);
        void load(string filename);
};


Coord3D xyz_to_sphere(Coord3D& xyz);
Coord3D sphere_to_xyz(Coord3D& sphere);

void normalise(Coord3D& xyz);
void normalise(Coord3D& xyz, double spd);

void normalise(Coord2D& xy);
void normalise(Coord2D& xy, double spd);

#endif
