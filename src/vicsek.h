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

const double PI = 3.141592653589793238463;

using namespace std;

using Property = Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)
using PropertyInt = Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)


/*
 * calculate the vicsek alignment of velocities
 */
template<class T>
void vicsek_align(T& velocities, const Conn& connections){
    int dim = velocities.rows();
    int n = velocities.cols();
    T new_velocities{dim, n};
    new_velocities.setZero();
    for (int i = 0; i < n; i++){
        for (auto j : connections[i]){
            new_velocities.col(i) += velocities.col(j);
        }
    }
    velocities = new_velocities;
}


class Vicsek3D{
    protected:
        double noise_;
        double speed_;
        double rc_;
        Conn connections_;
        void align(){
            vicsek_align(velocities_, connections_);
        }
        void add_noise();
        void rotate_noise(Coord3D& noise_xyz);
        void rotate_noise_xyz(Coord3D& noise_xyz);
        void rotate_noise_fast(Coord3D& noise_xyz);
        VerletList<Coord3D> verlet_list_;

    public:
        int n_;
        Coord3D positions_;
        Coord3D velocities_;

        Vicsek3D(int n, double r, double eta, double v0);
        void move(bool rebuild);
        void move();  // without neighbour list
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
        void move();  // without neighbour list
};


class Attanasi2014PCB : public Vicsek3D{
    protected:
        void harmonic_align();
        double beta_;
    public:
        Attanasi2014PCB(int n, double r, double eta, double v0, double beta);
        void move(bool rebuild);
        void move();  // without neighbour list
};


class Vicsek3DPBCVN : public Vicsek3DPBC{
    private:
        void noisy_align();
    public:
        Vicsek3DPBCVN(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void move();  // without neighbour list
};


class Vicsek3DPBCInertia : public Vicsek3DPBC{
    public:
        double alpha_;
        Coord3D old_velocities_;
        Vicsek3DPBCInertia(int n, double r, double eta, double box, double v0, double alpha);
        void move(bool rebuild);
        void move();  // without neighbour list
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
        void move();  // without neighbour list
};


class Vicsek2D{
    protected:
        double noise_;
        double speed_;
        double rc_;
        Conn connections_;
        void align(){
            vicsek_align(velocities_, connections_);
        }
        void add_noise();
        VerletList<Coord2D> verlet_list_;

    public:
        int n_;
        Coord2D positions_;
        Coord2D velocities_;
        Vicsek2D(int n, double r, double eta, double v0);
        void move(bool rebuild);
        void move();  // without neighbour list
};



class Vicsek2DPBC{
    protected:
        double noise_;
        double speed_;
        double box_;
        ConnMat conn_mat_;

        void align();
        void add_noise();
        inline void fix_positions() {
            positions_.array() -= box_ * (positions_.array() / box_).floor();
            }

    public:
        int n_;
        Coord2D positions_;
        Coord2D velocities_;
        CellList2D cell_list_;

        Vicsek2DPBC(int n, double r, double eta, double box, double v0);
        Vec2D get_shift(Vec2D p1, Vec2D p2);
        void move(bool rebuild);
        void move();  // without neighbour list
};


class Vicsek2DPBCVN : public Vicsek2DPBC{
    private:
        void noisy_align();
    public:
        Vicsek2DPBCVN(int n, double r, double eta, double box, double v0);
        void move(bool rebuild);
        void move();  // without neighbour list
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
        void move();  // without neighbour list
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


/*
 *  * please refer to 10.1007/s10955-014-1119-3 for the origin of the model
 */
class InertialSpin3D{
    private:
        double mass_;
        double imass_;
        double friction_;
        double T_; // temperature
        double J_; // coupling constant 
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
        void move();
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
            int n, double box, double r, double v0, double T, double j, double m, double f
        );
        void move(bool rebuid);
};


template <typename T>
vector<int> argsort(const vector<T> &v) {
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    stable_sort(
            idx.begin(), idx.end(),
            [&v](int i1, int i2) {return v[i1] < v[i2];}
            );
    return idx;
}


/*
 * get the n nearest neighbours, should work with any dimension
 */
template<class T>
Conn get_topology_connections(T positions, int nc){
    Conn connections;
    vector<int> sorted_indices;
    int n_total = positions.cols();
    for (int i = 0; i < n_total; i++){
        vector<double> squared_distances {};
        for (int j = 0; j < n_total; j++){
            squared_distances.push_back(
                   (positions.col(i) - positions.col(j)).squaredNorm()
                   );
        }
        sorted_indices = argsort(squared_distances);
        connections.push_back(vector<int> {
                sorted_indices.begin(), sorted_indices.begin() + nc
                });
    }
    return connections;
}


/*
 * get the neighbours within a metric range
 */
template<class T>
Conn get_connections(T& positions, double rc){
    Conn connections;
    double rc2 = rc * rc;
    for (int i=0; i<positions.cols(); i++){
        connections.push_back(vector<int>{});
        for (int j=0; j<positions.cols(); j++){
            if ((positions.col(i) - positions.col(j)).squaredNorm() <= rc2){
                connections[i].push_back(j);
            };
        }
    }
    return connections;
}


/*
 * get the neighbours within a metric range in a cubic PBC box
 */
template<class T>
Conn get_connections_pbc(T& positions, double rc, double box){
    Conn connections;
    double rc2 = rc * rc;
    double half_box_2 = box * box / 4;
    for (int i=0; i<positions.cols(); i++){
        connections.push_back(vector<int>{});
        for (int j=0; j<positions.cols(); j++){
            double dist_nd_2 = 0;
            for (auto shift_1d : (positions.col(i) - positions.col(j))){
                double dist_1d = abs(shift_1d);
                if (dist_1d > half_box_2){
                    dist_nd_2 += pow(box - dist_1d, 2);
                } else {
                    dist_nd_2 += pow(dist_1d, 2);
                }
            }
            if (dist_nd_2 <= rc2){
                connections[i].push_back(j);
            };
        }
    }
    return connections;
}


#endif