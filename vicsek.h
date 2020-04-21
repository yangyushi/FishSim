#ifndef VICSEK
#define VICSEK

#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <cmath>
#include <fstream>
#include "neighbour_list.h"

using namespace std;

using Coord3D = Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor>;  // (3, n)
using Property = Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)
using Vec3D = Eigen::Array3d;  // (3, )
using RotMat = Eigen::Matrix3d;


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
};


Coord3D xyz_to_sphere(Coord3D& xyz);

Coord3D sphere_to_xyz(Coord3D& sphere);

void normalise(Coord3D& xyz);

void normalise(Coord3D& xyz, double spd);

#endif
