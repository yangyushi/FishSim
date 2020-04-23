#ifndef NEIGHBOUR_LIST
#define NEIGHBOUR_LIST
#include <vector>
#include <array>
#include <map>
#include <Eigen/Dense>

using namespace std;

using DistMat = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; // (n, n)
using Indices = vector<int>;  // size is n
using Neighbours = vector<Indices>;

using Coord3D = Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor>;  // (3, n)
using CellIndices3D = Eigen::Array<int, 3, Eigen::Dynamic, Eigen::RowMajor>; // (3, n)
using Index3D = array<int, 3>;  // size is 3
using Indices3D = vector< Index3D >;
using Head3D = map<Index3D, int>;

using Coord2D = Eigen::Array<double, 2, Eigen::Dynamic, Eigen::RowMajor>;  // (2, n)
using CellIndices2D = Eigen::Array<int, 2, Eigen::Dynamic, Eigen::RowMajor>; // (2, n)
using Index2D = array<int, 2>;  // size is 3
using Indices2D = vector< Index2D >;
using Head2D = map<Index2D, int>;



Index3D unravel_index(int index, Index3D shape);
Index2D unravel_index(int index, Index2D shape);

Indices3D product_3d(Indices& arr);
Indices3D product_3d(Indices& arr_1, Indices& arr_2, Indices& arr_3);

Indices2D product_2d(Indices& arr);
Indices2D product_2d(Indices& arr_1, Indices& arr_2);


class CellList3D{
    /*
     * Using cell list to accelerate distance calculation with a cutoff for 3D simulation
     */
    private:
        double rc;
        double box;
        int sc;
        int size; // number of particles
        int ndim = 3;
        bool pbc;
        Indices clist;  // cell list
        Head3D chead;  // cell head

        inline Index3D chead_shape() {
            Index3D cs;
            for (int d = 0; d < ndim; d++){
                cs[d] = sc;
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

        inline Index2D chead_shape() {
            Index2D cs;
            for (int d = 0; d < ndim; d++){
                cs[d] = sc;
            }
            return cs;
        }
        void refill();
        Indices2D get_neighbours_indices(Index2D cell_idx);

    public:
        CellList2D(double r_cut, double box, bool pbc);
        void build(Coord2D& positions);
        DistMat get(Coord2D positions);
};


#endif
