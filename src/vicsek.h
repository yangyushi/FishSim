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

class Vicsek3D{
    protected:
        double noise_;
        double speed_;
        Conn connections_;
        void align();
        void add_noise();
        void rotate_noise(Coord3D& noise_xyz);
        void rotate_noise_fast(Coord3D& noise_xyz);
        void rotate_noise_xyz(Coord3D& noise_xyz);
        VerletList3D verlet_list_;

    public:
        int n_;
        Coord3D positions_;
        Coord3D velocities_;

        Vicsek3D(int n, double r, double eta, double v0);
        void move(bool rebuild);
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
        double noise_;
        double speed_;
        double box_;
        int n_;
        ConnMat conn_mat_;

        void align();
        void add_noise();
        inline void fix_positions(){positions_ -= box_ * (positions_ / box_).floor();}

    public:
        Coord2D positions_;
        Coord2D velocities_;
        CellList2D cell_list_;

        Vicsek2DPBC(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
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

/*
 * dump the current phase point to an xyz file
 * it works with both 2D and 3D system
 */
template<class T>
void dump(T system, string filename){
    ofstream f;
    f.open(filename, ios::out | ios::app);
    f << system.n_ << endl;
    if (system.positions_.rows() == 3){
        f << "id, x, y, z, vx, vy, vz" << endl;
        for (int i = 0; i < system.n_; i++ ) {
            f << i << " "
                << system.positions_(0, i)  << " "
                << system.positions_(1, i)  << " "
                << system.positions_(2, i)  << " "
                << system.velocities_(0, i) << " "
                << system.velocities_(1, i) << " "
                << system.velocities_(2, i) << endl;
        }
    } else if (system.positions_.rows() == 2){
        f << "id, x, y, vx, vy" << endl;
        for (int i = 0; i < system.n_; i++ ) {
            f << i << " "
              << system.positions_(0, i)  << " "
              << system.positions_(1, i)  << " "
              << system.velocities_(0, i) << " "
              << system.velocities_(1, i) << " " << endl;
        }
    } else {
        throw("invalid dimension");
    }
    f.close();
}


/*
 * load the phase point from the *last* frame of an xyz file
 * it works with both 2D and 3D system
 * The xyz file should be generated with the `dump` function
 */
template<class T>
void load(T system, string filename){
    ifstream f;
    string line;
    regex head_pattern{"\\d+"};
    smatch matched;
    int head_lines = 2;
    string num;
    int N = 0;
    int total_frame = 0;
    
    // find total number of frames
    f.open(filename, ios::in);
    while (f) {
        getline(f, line);
        if (regex_match(line, matched, head_pattern)){
            N = stoi(line);
            total_frame += 1;
            for (int i=0; i<N; i++) getline(f, line);
        }
    }
    f.close();
    
    // jump to the last frame 
    f.open(filename, ios::in);
    for (int i = 0; i < total_frame - 1; i++){
        for (int j = 0; j < N + head_lines; j++){
        getline(f, line);
        }
    }
    
    // load the data
    if (system.positions_.rows() == 3){
        for (int i = 0; i < N + head_lines; i++){
            getline(f, line);
            if (i > 1) {
                istringstream ss(line);
                ss >> num;
                for (int j = 0; j < 3; j++){
                    ss >> system.positions_(j, i - head_lines);
                }
                for (int j = 0; j < 3; j++){
                    ss >> system.velocities_(j, i - head_lines);
                }
            }
        }
    } else if (system.positions_.rows() == 2) {
        for (int i = 0; i < N + head_lines; i++){
            getline(f, line);
            
            if (i > 1) {
                istringstream ss(line);
                ss >> num;
                for (int j = 0; j < 2; j++){
                    ss >> system.positions_(j, i - head_lines);
                }
                for (int j = 0; j < 2; j++){
                    ss >> system.velocities_(j, i - head_lines);
                }
            }
        }
    } else {
        throw("invalid dimension");
    }
    f.close();
}

#endif
