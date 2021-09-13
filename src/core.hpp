#ifndef CORE
#define CORE

#include <map>
#include <cmath>
#include <array>
#include <regex>
#include <random>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>

using Property = Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)
using PropertyInt = Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)

using Conn = std::vector< std::vector <int> >;
using DistMat = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; // (n, n)
using ConnMat = Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; // (n, n)
using Indices = std::vector<int>;  // size is n
using Neighbours = std::vector< Indices >;

using Coord3D = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;  // (3, n)
using CellIndices3D = Eigen::Array<int, 3, Eigen::Dynamic, Eigen::RowMajor>; // (3, n)
using Index3D = std::array<int, 3>;  // size is 3
using Indices3D = std::vector< Index3D >;
using Head3D = std::map<Index3D, int>;

using Coord2D = Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>;  // (2, n)
using CellIndices2D = Eigen::Array<int, 2, Eigen::Dynamic, Eigen::RowMajor>; // (2, n)
using Index2D = std::array<int, 2>;  // size is 3
using Indices2D = std::vector< Index2D >;
using Head2D = std::map<Index2D, int>;

using RotMat = Eigen::Matrix3d;
using Vec2D = Eigen::Matrix<double, 2, 1>;  // (2, 1)
using Vec3D = Eigen::Matrix<double, 3, 1>;  // (3, 1)


// Same behaviour as ``numpy.unravel_index``
Index3D unravel_index(int index, Index3D shape);
Index2D unravel_index(int index, Index2D shape);

// Calculate Cartesian product of arr repeated 3 times
Indices3D product_3d(Indices& arr);
Indices3D product_3d(Indices& arr_1, Indices& arr_2, Indices& arr_3);

// Calculate Cartesian product of arr repeated 2 times
Indices2D product_2d(Indices& arr);
Indices2D product_2d(Indices& arr_1, Indices& arr_2);

// coordinate transform
Coord3D xyz_to_sphere(Coord3D& xyz);
Coord3D sphere_to_xyz(Coord3D& sphere);
Property xy_to_phi(Coord2D& xy);
Coord2D phi_to_xy(Property& phi);
Coord2D phi_to_xy(Property& phi, double spd);


/*
 * get the rotation matrix to rotate vector (0, 0, 1) to (x, y, z)
 */
RotMat get_rotation_matrix(double x, double y, double z);


template <typename T>
std::vector<int> argsort(const std::vector<T> &v) {
    std::vector<int> idx(v.size());
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
    std::vector<int> sorted_indices;
    int n_total = positions.cols();
    for (int i = 0; i < n_total; i++){
        std::vector<double> squared_distances {};
        for (int j = 0; j < n_total; j++){
            squared_distances.push_back(
                   (positions.col(i) - positions.col(j)).squaredNorm()
                   );
        }
        sorted_indices = argsort(squared_distances);
        connections.push_back(std::vector<int> {
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
        connections.push_back(std::vector<int>{});
    }
    #pragma omp parallel for
    for (int i=0; i<positions.cols(); i++){
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
    double half_box = box / 2;
    for (int i=0; i<positions.cols(); i++){
        connections.push_back(std::vector<int>{});
    }
    #pragma omp parallel for
    for (int i=0; i<positions.cols(); i++){
        for (int j=0; j<positions.cols(); j++){
            double dist_nd_2 = 0;
            for (auto shift_1d : (positions.col(i) - positions.col(j))){
                double dist_1d = abs(shift_1d);
                if (dist_1d > half_box){
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


template<class T>
void normalise(T& xyz){
    int i = 0;
    for (auto L : xyz.colwise().norm()){
        xyz.col(i) /= L;
        i++;
    }
}


template<class T>
void normalise(T& xyz, double norm_value){
    int i = 0;
    for (auto L : xyz.colwise().norm()){
        xyz.col(i) = xyz.col(i) / L * norm_value;
        i++;
    }
}

/*
 * calculate the vicsek alignment of velocities with [v]ectorial [n]oise
 *     it should work with any dimension
 * See ginellEPJ2016 for equations
 */
template<class T>
void vicsek_align_vn(T& velocities, const Conn& connections, double noise){
    int dim = velocities.rows();
    int n = velocities.cols();
    T new_velocities{dim, n};
    T noise_nd{dim, n};
    new_velocities.setZero();
    noise_nd.setRandom();
    normalise(noise_nd, noise);
    for (int i = 0; i < n; i++){
        for (auto j : connections[i]){
            new_velocities.col(i) += velocities.col(j);
        }
        new_velocities.col(i) += noise_nd.col(i) * connections[i].size();
    }
    velocities = new_velocities;
}


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


/*
 * dump the current phase point to an xyz file
 * it works with both 2D and 3D system
 */
template<class T>
void dump(T system, std::string filename){
    std::ofstream f;
    f.open(filename, std::ios::out | std::ios::app);
    f << system.n_ << std::endl;
    if (system.positions_.rows() == 3){
        f << "id, x, y, z, vx, vy, vz" << std::endl;
        for (int i = 0; i < system.n_; i++ ) {
            f << i << " "
                << system.positions_(0, i)  << " "
                << system.positions_(1, i)  << " "
                << system.positions_(2, i)  << " "
                << system.velocities_(0, i) << " "
                << system.velocities_(1, i) << " "
                << system.velocities_(2, i) << std::endl;
        }
    } else if (system.positions_.rows() == 2){
        f << "id, x, y, vx, vy" << std::endl;
        for (int i = 0; i < system.n_; i++ ) {
            f << i << " "
              << system.positions_(0, i)  << " "
              << system.positions_(1, i)  << " "
              << system.velocities_(0, i) << " "
              << system.velocities_(1, i) << " " << std::endl;
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
void load(T system, std::string filename){
    std::ifstream f;
    std::string line;
    std::regex head_pattern{"\\d+"};
    std::smatch matched;
    int head_lines = 2;
    std::string num;
    int N = 0;
    int total_frame = 0;
    
    // find total number of frames
    f.open(filename, std::ios::in);
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
    f.open(filename, std::ios::in);
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
                std::istringstream ss(line);
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
                std::istringstream ss(line);
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


class AVS2D{  // [A]ligning [V]ectors with [S]carlar noise in [2D]
    public:
        int n_;
        double speed_ = 1;
        const int dim_ = 2;
        Coord2D velocities_; 

        AVS2D(int n, double noise, double speed);
        void add_noise();

        Coord2D& get_velocities() { return velocities_; };
        void load_velocities(Coord2D);
        inline double get_polarisation() {
            return velocities_.rowwise().sum().norm() / n_ / speed_;
        };

    protected:
        double noise_;
};


class AVS3D{  // [A]ligning [V]ectors with [S]carlar noise in [3D]
    public:
        int n_;
        double speed_ = 1;
        const int dim_ = 3;
        Coord3D velocities_; 

        AVS3D(int n, double noise, double speed);
        Coord3D& get_velocities() { return velocities_; };
        void load_velocities(Coord3D);
        inline double get_polarisation() {
            return velocities_.rowwise().sum().norm() / n_ / speed_;
        };
        void add_noise();

    protected:
        double noise_;
};


#endif
