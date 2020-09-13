#ifndef NETWORK
#define NETWORK

#include <vector>
#include <array>
#include <map>
#include <Eigen/Dense>
#include "vicsek.h"
#include "neighbour_list.h"

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
        VerletList3D verlet_list_;

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


#endif
