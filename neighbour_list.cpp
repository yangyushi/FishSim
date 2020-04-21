#include "neighbour_list.h"


CellList3D::CellList3D(double r_cut, double box, bool pbc) :
    rc(r_cut), box(box), pbc(pbc) {
    ndim = 0;
    size = 0;
    sc = floor(box / rc / 2);
}


void CellList3D::refill(){
    /*
    * Filling chead and clist with zeros
    */
    chead.clear();
    vector<int> sc_range;
    for (int i=0; i < sc; i++){
        sc_range.push_back(i);
    }

    for (auto item : product_3d(sc_range)){
        chead[idx] = 0;
    }

    clist.clear();
    for (int i = 0; i < size + 1; i++){
        clist.push_back(0);
    }
}


Indices3D CellList3D::get_neighbours_indices(vector<int> cell_idx){
    /*
     * For a cell indexed as (a, b, c), find all its neighbours
     * The cell itself is included as a "neighbour"
     */
    Indices3D neighbours; // (dimension, neighbour_num)

    for (int d = 0; d < ndim; d++){
        neighbours.push_back(Index3D{cell_idx[d]});
        if (cell_idx[d] == 0){
            if (pbc) {
                neighbours.back().push_back(this->sc - 1);
                neighbours.back().push_back(1);
            } else {
                neighbours.back().push_back(1);
            }
        } else if (cell_idx[d] == this->sc - 1){
            if (pbc) {
                neighbours.back().push_back(cell_idx[d] - 1);
                neighbours.back().push_back(0);
            } else {
                neighbours.back().push_back(cell_idx[d] - 1);
            }
        } else {
            neighbours.back().push_back(cell_idx[d] + 1);
            neighbours.back().push_back(cell_idx[d] - 1);
        }
    }

    return product_3d(neighbours[0], neighbours[1], neighbours[2]);
}


void CellList3D::build(Coord3D& positions){
    ndim = positions.rows();
    size = positions.cols();
    CellIndices ci(ndim, size);  // discretise positions into cell indices
    ci << floor(positions / box * sc).cast<int>();

    refill();

    Index3D cell_idx;
    for (int i=1; i < size + 1; i++){
        for (int j = 0; j < ndim; j++){
            cell_idx.push_back(ci(j, i-1));
        }
        clist[i] = chead[cell_idx];
        chead[cell_idx] = i;
        cell_idx.clear();
    }
}


DistMat CellList3D::get(Coord3D positions){
    DistMat dist_mat(size, size);
    dist_mat.setConstant(-1);

    double rc2 = rc * rc;
    double dist2 = 0;
    double dist_1d = 0;

    CellIndices ci(ndim, size);
    ci << floor(positions / box * sc).cast<int>();

    Indices3D neighbours;
    Indices in_neighbour;
    Indices in_cell;
    Index3D cell_idx;
    int cursor;
    // itering over all cells
    for (int x=0; x < pow(sc, ndim); x++){ 
        in_cell.clear();
        in_neighbour.clear();
        cell_idx = unravel_index(x, chead_shape());

        // ignore empty cells
        if (chead[cell_idx] == 0) continue;

        // Collecting particle indices in current cell
        in_cell.push_back(chead[cell_idx]);
        cursor = chead[cell_idx];
        while (clist[cursor] > 0) {
            cursor = clist[cursor];
            in_cell.push_back(cursor);
        }

        // Collecting particle indices in neighbour cells
        neighbours = get_neighbours_indices(cell_idx);
        for (auto nc : neighbours){ // nc -> neighbour_cell
            if (chead[nc] == 0) continue;
            in_neighbour.push_back(chead[nc]);
            cursor = chead[nc];
            while (clist[cursor] > 0){
                cursor = clist[cursor];
                in_neighbour.push_back(cursor);
            }
        }

        for (auto i : in_cell){
        for (auto j : in_neighbour){
            dist2 = 0;
            for (int d = 0; d < ndim; d++){
                dist_1d = abs(positions(d, i-1) - positions(d, j-1));
                if (pbc and dist_1d > box / 2){
                    dist_1d = box - dist_1d;
                }
                dist2 += dist_1d * dist_1d;
            }
            if (dist2 < rc2){
                dist_mat(i-1, j-1) = dist2;
            }
        }}
    }

    dist_mat = (dist_mat > 0).select(dist_mat.sqrt(), dist_mat);

    return dist_mat;
}


Index3D unravel_index(int index, vector<int> shape){
    /*
     * Same behaviour as ``numpy.unravel_index``
     */
    Index3D result;
    int dim = shape.size();
    int size;
    int tmp;
    for (int d1 = 0; d1 < dim; d1++){
        size = 1;
        for (int d2 = d1 + 1; d2 < dim; d2++){
            size *= shape[d2];
        }
        tmp = floor(index / size);
        result.push_back(tmp);
        index -= tmp * size;
    }
    return result;
}


Indices3D product_3d(Indices& arr){
    /*
     * Calculate Cartesian product of indices in different dimensions
     */
    Indices3D result;
    Index3D temp{0, 0, 0};
    int size = arr.size();
    for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
    for (int k = 0; k < size; k++){
        temp[0] = arr[i];
        temp[1] = arr[j];
        temp[2] = arr[k];
        result.push_back(temp);
    }}}
    return result;
}


Indices3D product_3d(Indices& arr_1, Indices& arr_2, Indices& arr_3){
    /*
     * Calculate Cartesian product of arr repeated 3 times
     */
    Indices3D result;
    Index3D temp{0, 0, 0};
    int size = arr_1.size();
    for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
    for (int k = 0; k < size; k++){
        temp[0] = arr_1[i];
        temp[1] = arr_2[j];
        temp[2] = arr_3[k];
        result.push_back(temp);
    }}}
    return result;
}
