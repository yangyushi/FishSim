#ifndef NEIGHBOUR_LIST
#define NEIGHBOUR_LIST

#include "core.hpp"


template<class T>
class VerletList{
    /*
     * Using the Verlet list to accelerate distance calculation with a
     * cutoff for 3D simulation
     * This is suitable for accelerating the simulation without a box
     */
    private:
        double rc_;
        double rc2_;
        double rl_;
        double rl2_;
        std::vector<int> nlist_;
        std::vector<int> point_;
        int point_size_ = 0;  // particle number + 1, point_.size()
        int size_ = 0;  // particle number
        Conn get_blank_connections(){
            Conn connections{};
            for (int i = 0; i < size_; i++){
                connections.push_back(std::vector<int>{});
            }
            return connections;
        }

    public:
        VerletList(double r_cut, double r_skin);
        void build(const T& positoins);
        void get_dmat(const T& positoins, DistMat& dist_mat);
        void get_cmat(const T& positoins, ConnMat& conn_mat);
        void get_cmat_slow(const T& positoins, ConnMat& conn_mat);
        Conn get_conn(const T& positions);
        /*
         * find many different connections:
         *     {dij < r1}, {r1 <= dij < r2} , ... , {rc_ <= dij}
         */
        std::vector<Conn> get_conn(const T& positions, std::vector<double> r_vals);

        Conn get_neighbour_conn(const T& positions);
        std::vector<Conn> get_neighbour_conn(const T& positions, std::vector<double> r_vals);
        std::vector<Conn> get_neighbour_conn_slow(const T& positions, std::vector<double> r_vals);
        Conn get_conn_slow(const T& positions);
};


template<class T>
VerletList<T>::VerletList(double r_cut, double r_skin)
    : rc_{r_cut}, rl_{r_skin} {
        rl2_ = rl_ * rl_;
        rc2_ = rc_ * rc_;
}


template<class T>
void VerletList<T>::build(const T& positions){
    double d2{0};
    size_ = positions.cols();
    point_.clear();
    nlist_.clear();

    point_.push_back(0);
    for (int i = 0; i < size_; i++){
        for (int j = 0; j < size_; j++){
            d2 = (positions.col(i) - positions.col(j)).array().pow(2).sum();
            if (d2 < rl2_){
                nlist_.push_back(j);
            }
        }
        point_.push_back(nlist_.size());
    }
    point_size_ = point_.size(); // particle number + 1
}


template<class T>
void VerletList<T>::get_dmat(const T& positions, DistMat& dist_mat){
    dist_mat.setConstant(-1);
    #pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i+1]; j++){
            idx_j = nlist_[j];
            dist2 = (positions.col(i) - positions.col(idx_j)).array().pow(2).sum();
            if (dist2 < rc2_){
                dist_mat(i, idx_j) = sqrt(dist2);
            }
        }
    }
}


template<class T>
void VerletList<T>::get_cmat(const T& positions, ConnMat& conn_mat){
    conn_mat.setZero();
    #pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i + 1]; j++){
            idx_j = nlist_[j];
            dist2 = (positions.col(i) - positions.col(idx_j)).array().pow(2).sum();
            if (dist2 < rc2_){
                conn_mat(i, idx_j) = 1;
            }
        }
    }
}


template<class T>
void VerletList<T>::get_cmat_slow(const T& positions, ConnMat& conn_mat){
    conn_mat.setZero();
    #pragma omp parallel for
    for (int i = 0; i < size_; i++){
        for (int j = 0; j < size_; j++){
            double dist2 = (positions.col(i) - positions.col(j)).squaredNorm();
            if (dist2 < rc2_){
                conn_mat(i, j) = 1;
            }
        }
    }
}


template<class T>
Conn VerletList<T>::get_conn(const T& positions){
    Conn connections = get_blank_connections();
    #pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i + 1]; j++){
            idx_j = nlist_[j];
            dist2 = (positions.col(i) - positions.col(idx_j)).array().pow(2).sum();
            if (dist2 < rc2_){
                connections[i].push_back(idx_j);
            }
        }
    }
    return connections;
}

template<class T>
Conn VerletList<T>::get_neighbour_conn(const T& positions){
    Conn connections = get_blank_connections();
    #pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i + 1]; j++){
            idx_j = nlist_[j];
            if (idx_j == i){
                continue;
            }
            dist2 = (positions.col(i) - positions.col(idx_j)).array().pow(2).sum();
            if (dist2 < rc2_){
                connections[i].push_back(idx_j);
            }
        }
    }
    return connections;
}



/*
 * get the indices of neighbours in different regions 
 *
 *  Args:
 *      r_bins: [r1, r2, ..., rk], where r1 < r2 < ... rk < rc_
 *
 *  Return:
 *     [c1, c2, ..., cn]  where ci is the connection in region i.
 *     and ci counts the region where r(i) < dij <= r(i+1)
 */
template<class T>
std::vector<Conn> VerletList<T>::get_conn(
        const T& positions, std::vector<double> r_bins
    ){
    size_t conn_num = r_bins.size() - 1;
    std::vector<Conn> conn_list;
    std::vector<double> r2_bins;
    for (int i = 0; i < conn_num; i++){
        conn_list.push_back(get_blank_connections());
        r2_bins.push_back(r_bins[i] * r_bins[i]);
    }
    r2_bins.push_back(r_bins[conn_num] * r_bins[conn_num]);
    #pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i + 1]; j++){
            idx_j = nlist_[j];
            dist2 = (positions.col(i) - positions.col(idx_j)).array().pow(2).sum();
            for (int k=0; k < conn_num; k++){
                if ((dist2 > r2_bins[k]) and (dist2 <= r2_bins[k+1])){
                    conn_list[k][i].push_back(idx_j);
                    break;
                }
            }
        }
    }
    return conn_list;
}


/*
 * not including self-loop
 */
template<class T>
std::vector<Conn> VerletList<T>::get_neighbour_conn(
        const T& positions, std::vector<double> r_bins
    ){
    size_t conn_num = r_bins.size() - 1;
    std::vector<Conn> conn_list;
    std::vector<double> r2_bins;
    for (int i = 0; i < conn_num; i++){
        conn_list.push_back(get_blank_connections());
        r2_bins.push_back(r_bins[i] * r_bins[i]);
    }
    r2_bins.push_back(r_bins[conn_num] * r_bins[conn_num]);

    //#pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i + 1]; j++){
            idx_j = nlist_[j];
            if (idx_j == i){
                continue;
            }
            dist2 = (positions.col(i) - positions.col(idx_j)).array().pow(2).sum();
            std::cout << dist2 << " -- " << r2_bins[1] << std::endl;
            for (int k=0; k < conn_num; k++){ 
                if ((dist2 > r2_bins[k]) and (dist2 <= r2_bins[k+1])){
                    conn_list[k][i].push_back(idx_j);
                    break;
                }
            }
        }
    }
    return conn_list;
}


template<class T>
std::vector<Conn> VerletList<T>::get_neighbour_conn_slow(
        const T& positions, std::vector<double> r_bins
    ){
    size_t conn_num = r_bins.size() - 1;
    std::vector<Conn> conn_list;
    std::vector<double> r2_bins;
    for (int i = 0; i < conn_num; i++){
        conn_list.push_back(get_blank_connections());
        r2_bins.push_back(r_bins[i] * r_bins[i]);
    }
    r2_bins.push_back(r_bins[conn_num] * r_bins[conn_num]);

    //#pragma omp parallel for
    for (int i = 0; i < size_; i++){
        double dist2 = 0;
        for (int j = 0; j < size_; j++){
            if (j == i){
                continue;
            }
            dist2 = (positions.col(i) - positions.col(j)).array().pow(2).sum();
            for (int k=0; k < conn_num; k++){ 
                if ((dist2 > r2_bins[k]) and (dist2 <= r2_bins[k+1])){
                    conn_list[k][i].push_back(j);
                    break;
                }
            }
        }
    }
    return conn_list;
}



template<class T>
Conn VerletList<T>::get_conn_slow(const T& positions){
    Conn connections = get_blank_connections();
    #pragma omp parallel for
    for (int i = 0; i < size_; i++){
        for (int j = 0; j < size_; j++){
            double dist2 = (positions.col(i) - positions.col(j)).array().pow(2).sum();
            if (dist2 < rc2_){
                connections[i].push_back(j);
            }
        }
    }
    return connections;
}


class CellList3D{
    /*
     * Using cell list to accelerate distance calculation with a cutoff for 3D simulation
     */
    private:
        double rc_;
        double box_;
        bool pbc_;
        int size_; // number of particles
        const int ndim_ = 3;
        int sc_;  // box is divided into sc * sc * sc cells
        Indices clist_;  // cell list
        Head3D chead_;  // cell head
        Index3D head_shape_;
        void refill();
        void fill_neighbour_indices_1d(Neighbours& neighbours, int cell_idx_1d);
        Indices3D get_neighbours_indices(const Index3D& cell_idx);
        CellIndices3D get_cell_indices(const Coord3D& positions);

    public:
        CellList3D(double r_cut, double box, bool pbc);
        void update_sc(int new_sc);
        void build(Coord3D& positions);
        void get_dmat(Coord3D& positions, DistMat& dist_mat);
        void get_cmat(Coord3D& positions, ConnMat& conn_mat);
        Conn get_conn(Coord3D& positions);
};


class CellList2D{
    /*
     * Using cell list to accelerate distance calculation with a cutoff for 2D simulation
     */
    private:
        double rc;
        double box;
        int sc;
        int size; // number of particles
        int ndim = 2;
        bool pbc;
        Indices clist;  // cell list
        Head2D chead;  // cell head
        Index2D head_shape_;

        void refill();
        Indices2D get_neighbours_indices(Index2D cell_idx);

    public:
        CellList2D(double r_cut, double box, bool pbc);
        void update_sc(int new_sc);
        void build(Coord2D& positions);
        void get_dmat(Coord2D& positions, DistMat& dist_mat);
        void get_cmat(Coord2D& positions, ConnMat& conn_mat);
        Conn get_conn(Coord2D& positions);
};


#endif
