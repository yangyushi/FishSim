#ifndef NEIGHBOUR_LIST
#define NEIGHBOUR_LIST
#include <vector>
#include <array>
#include <map>
#include <Eigen/Dense>

using namespace std;

using Coord3D = Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor>;  // (3, n)
using CellIndices = Eigen::Array<int, 3, Eigen::Dynamic, Eigen::RowMajor>; // (3, n)
using DistMat = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; // (n, n)

using Index3D = array<int, 3>;  // size is 3
using Index2D = array<int, 2>;  // size is 3

using Indices3D = vector< Index3D >;
using Indices2D = vector< Index2D >;

using Indices = vector<int>;  // size is n
using Neighbours = vector<Indices>;
using Head = map<Index3D, int>;

class CellList3D{
    /*
     * Using cell list to accelerate distance matrix calculation with a cutoff
     * It works with N-dimensional data in a *(hyper) cube*
     */
    private:
        double rc;
        double box;
        int sc;
        int size; // number of particles
        int ndim;
        bool pbc;
        Indices clist;  // cell list
        Head chead;  // cell head

        inline vector<int> chead_shape() {
            vector<int> cs;
            for (int i = 0; i < ndim; i++){
                cs.push_back(sc);
            }
            return cs;
        }
        void refill();
        Indices3D get_neighbours_indices(Index3D cell_idx);

    public:
        CellList3D(double r_cut, double box, bool pbc);
        void build(Coord3D& positions);
        DistMat get(Coord3D positions);
};

Index3D unravel_index(int index, vector<int> shape);

Indices3D product_3d(Indices& arr);
Indices3D product_3d(Indices& arr_1, Indices& arr_2, Indices& arr_3);

#endif
