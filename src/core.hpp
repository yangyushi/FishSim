#ifndef CORE
#define CORE

#include <map>
#include <cmath>
#include <array>
#include <regex>
#include <random>
#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>

using PropertyComplex = Eigen::Array<
    std::complex<double>, 1, Eigen::Dynamic, Eigen::RowMajor
    >;  // (1, n)
using Property = Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)
using Spins = Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>;  // (1, n)

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

const double _EPS = 100 * std::numeric_limits<double>::epsilon();


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


// get the rotation matrix to rotate vector (0, 0, 1) to (x, y, z)
RotMat get_rotation_matrix(double x, double y, double z);

// get the rotation matrix to rotate vector v1 to v2
RotMat get_rotation_matrix(Vec3D v1, Vec3D v2);

/* 
 * get the rotation matrix to rotate vector v1 to v2, with the maximum rotation
 * angle of theta.
 * If the angle between v1 and v2 is greater than theta, then only rotate by
 * the angle of theta, in the correct direction
*/
RotMat get_rotation_matrix(Vec3D v1, Vec3D v2, double theta);


template<class T>
void check_nan(T numbers, std::string info){
    for (auto num : numbers){
        if (std::isnan(num)){
            std::cout << info << std::endl;
            throw std::runtime_error("null value in ");
        }
    }
}


template <typename T>
double get_angle(T v1, T v2){
    double tmp = v1.transpose() * v2;
    tmp = tmp / v1.norm();
    tmp = tmp / v2.norm();
    return acos(tmp);
}

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
    //#pragma omp parallel for
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
 * calculate the vicsek alignment of orientations with [v]ectorial [n]oise
 *     it should work with any dimension
 * See ginellEPJ2016 for equations
 */
template<class T>
void vicsek_align_vn(T& orientations, const Conn& connections, double noise){
    int dim = orientations.rows();
    int n = orientations.cols();
    T new_orientations{dim, n};
    T noise_nd{dim, n};
    new_orientations.setZero();
    noise_nd.setRandom();
    normalise(noise_nd, noise);
    for (int i = 0; i < n; i++){
        for (auto j : connections[i]){
            new_orientations.col(i) += orientations.col(j);
        }
        new_orientations.col(i) += noise_nd.col(i) * connections[i].size();
    }
    orientations = new_orientations;
    normalise(orientations);
}


/*
 * calculate the vicsek alignment of orientations
 */
template<class T>
void vicsek_align(T& orientations, const Conn& connections){
    size_t dim = orientations.rows();
    size_t n = orientations.cols();
    T new_orientations{dim, n};
    new_orientations.setZero();
    for (size_t i = 0; i < n; i++){
        for (auto j : connections[i]){
            new_orientations.col(i) += orientations.col(j);
        }
    }
    orientations = new_orientations;
    normalise(orientations);
}


template<class T>
double get_mill_order(const T& orientations, const T& positions){
    Vec3D rot {0, 0, 0};
    Vec3D r, v, m;
    Vec3D group_centre = positions.rowwise().mean();
    size_t n = orientations.cols();
    for (size_t i = 0; i < n; i++){
        r << positions.col(i).array() - group_centre.array();
        v << orientations.col(i).array();
        r = r / r.norm();
        m << v.cross(r);
        rot += m;
    }
    return rot.norm() / (float) n;
}

void voter_align(Spins& spins, Property& states, const Conn& connections);


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
        const int dim_ = 2;
        Coord2D orientations_; 

        AVS2D(int n, double noise);
        void add_noise();

        Coord2D& get_orientations() { return orientations_; };
        void load_orientations(Coord2D);
        inline double get_polarisation() {
            return orientations_.rowwise().sum().norm() / n_;
        };

    protected:
        double noise_;
};


class AVS3D{  // [A]ligning [V]ectors with [S]carlar noise in [3D]
    public:
        int n_;
        const int dim_ = 3;
        Coord3D orientations_; 

        AVS3D(int n, double noise);
        virtual ~AVS3D(){};
        Coord3D& get_orientations() { return orientations_; };
        void load_orientations(Coord3D);
        inline double get_polarisation() {
            return orientations_.rowwise().sum().norm() / n_;
        };
        void add_noise();

    protected:
        double noise_;
};


class ABS{  // [A]ligning [B]inary [S]pins
    public:
        int n_;
        int dim_ = 1;
        Spins spins_;  // takse value of -1 or +1
        ABS(int n, double noise);
        void add_noise();
        inline Spins get_spins(){ return spins_; };
        inline void load_spins(Spins to_load) { spins_ << to_load; };
        inline double get_polarisation() {
            return abs(static_cast<float>(spins_.sum()) / n_);
        };

    protected:
        double noise_;
        Property states_;  // helper array
};

#endif
