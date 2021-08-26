#include "neighbour_list.h"


Index3D unravel_index(int index, Index3D shape){
    /*
     * Same behaviour as ``numpy.unravel_index``
     */
    int dim = shape.size();
    Index3D result;
    int size;
    int tmp;
    for (int d1 = 0; d1 < dim; d1++){
        size = 1;
        for (int d2 = d1 + 1; d2 < dim; d2++){
            size *= shape[d2];
        }
        tmp = floor(index / size);
        result[d1] = tmp;
        index -= tmp * size;
    }
    return result;
}


Index2D unravel_index(int index, Index2D shape){
    /*
     * Same behaviour as ``numpy.unravel_index``
     */
    int dim = shape.size();
    Index2D result;
    int size;
    int tmp;
    for (int d1 = 0; d1 < dim; d1++){
        size = 1;
        for (int d2 = d1 + 1; d2 < dim; d2++){
            size *= shape[d2];
        }
        tmp = floor(index / size);
        result[d1] = tmp;
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


Indices2D product_2d(Indices& arr){
    /*
     * Calculate Cartesian product of indices in different dimensions
     */

    Indices2D result;
    Index2D temp{0, 0};

    int size = arr.size();
    for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
        temp[0] = arr[i];
        temp[1] = arr[j];
        result.push_back(temp);
    }}
    return result;
}


Indices2D product_2d(Indices& arr_1, Indices& arr_2){
    /*
     * Calculate Cartesian product of arr repeated 3 times
     */
    Indices2D result;
    Index2D temp{0, 0};

    int size = arr_1.size();

    for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
        temp[0] = arr_1[i];
        temp[1] = arr_2[j];
        result.push_back(temp);
    }}
    return result;
}


CellList3D::CellList3D(double r_cut, double box, bool pbc)
    : rc_(r_cut), box_(box), pbc_(pbc) {
    ndim_ = 3;
    size_ = 0;
    sc_ = floor(box_ / rc_ / 2);
    if (sc_ < 1) {
        sc_ = 1;
    }
    for (int d = 0; d < ndim_; d++){
        head_shape_[d] = sc_;
    }
}


void CellList3D::update_sc(int new_sc){
    if (new_sc < 1){
        new_sc = 1;  // at least split into 1 x 1 x 1
    }
    sc_ = new_sc;
    for (int d = 0; d < ndim_; d++){
        head_shape_[d] = sc_;
    }
}


void CellList3D::refill(){
    /*
    * Filling chead and clist with zeros
    */

    Indices sc_range;
    for (int i=0; i < sc_; i++){
        sc_range.push_back(i);
    }

    chead_.clear();
    for (auto item : product_3d(sc_range)){
        chead_[item] = 0;
    }

    clist_.clear();
    for (int i = 0; i < size_ + 1; i++){
        clist_.push_back(0);
    }
}


CellIndices3D CellList3D::get_cell_indices(const Coord3D& positions){
    /*
     * Discretise positions into cell indices
     */
    CellIndices3D ci(ndim_, size_);
    ci << floor(positions.array() / box_ * sc_).cast<int>();
    return ci;
}


void CellList3D::build(Coord3D& positions){
    size_ = positions.cols();
    CellIndices3D ci = get_cell_indices(positions);
    refill();
    Index3D cell_idx;
    for (int i = 1; i < size_ + 1; i++){
        for (int d = 0; d < ndim_; d++){
            cell_idx[d] = ci(d, i-1);
        }
        clist_[i] = chead_[cell_idx];
        chead_[cell_idx] = i;
    }
}


void CellList3D::get_dmat(Coord3D& positions, DistMat& dist_mat){
    dist_mat.setConstant(-1);
    CellIndices3D ci = get_cell_indices(positions);
    double rc2 = rc_ * rc_;
    int max_cell_idx = pow(sc_, ndim_);

    #pragma omp parallel for
    for (int x = 0; x < max_cell_idx; x++){ 

        int cursor = 0;
        double dist_1d = 0;
        Indices3D neighbours;
        Indices in_neighbour;
        Indices in_cell;
        Index3D cell_idx = unravel_index(x, head_shape_);

        // ignore empty cells
        if (chead_[cell_idx] == 0) continue;

        // Collecting particle indices in current cell
        in_cell.push_back(chead_[cell_idx]);
        cursor = chead_[cell_idx];
        while (clist_[cursor] > 0) {
            cursor = clist_[cursor];
            in_cell.push_back(cursor);
        }

        // Collecting particle indices in neighbour cells
        neighbours = get_neighbours_indices(cell_idx);

        for (auto nc : neighbours){ // nc -> neighbour_cell
            if (chead_[nc] == 0) continue;
            in_neighbour.push_back(chead_[nc]);
            cursor = chead_[nc];
            while (clist_[cursor] > 0){
                cursor = clist_[cursor];
                in_neighbour.push_back(cursor);
            }
        }

        for (int i : in_cell){
        for (int j : in_neighbour){
            double dist2 = 0;
            for (int d = 0; d < ndim_; d++){
                dist_1d = abs(positions(d, i-1) - positions(d, j-1));
                if (pbc_ and dist_1d > box_ / 2){
                    dist_1d = box_ - dist_1d;
                }
                dist2 += dist_1d * dist_1d;
            }
            if (dist2 < rc2){
                dist_mat(i-1, j-1) = dist2;
            }
        }}
    }

    dist_mat = (dist_mat > 0).select(dist_mat.sqrt(), dist_mat);
}



void CellList3D::get_cmat(Coord3D& positions, ConnMat& conn_mat){
    conn_mat.setZero();

    CellIndices3D ci = get_cell_indices(positions);
    double rc2 = rc_ * rc_;
    int max_cell_idx = pow(sc_, ndim_);

    #pragma omp parallel for
    for (int x = 0; x < max_cell_idx; x++){ 

        int cursor = 0;
        double dist_1d = 0;
        Indices3D neighbours;
        Indices in_neighbour;
        Indices in_cell;
        Index3D cell_idx = unravel_index(x, head_shape_);

        // ignore empty cells
        if (chead_[cell_idx] == 0) continue;

        // Collecting particle indices in current cell
        in_cell.push_back(chead_[cell_idx]);
        cursor = chead_[cell_idx];
        while (clist_[cursor] > 0) {
            cursor = clist_[cursor];
            in_cell.push_back(cursor);
        }

        // Collecting particle indices in neighbour cells
        neighbours = get_neighbours_indices(cell_idx);

        for (auto nc : neighbours){ // nc -> neighbour_cell
            if (chead_[nc] == 0) continue;
            in_neighbour.push_back(chead_[nc]);
            cursor = chead_[nc];
            while (clist_[cursor] > 0){
                cursor = clist_[cursor];
                in_neighbour.push_back(cursor);
            }
        }

        for (int i : in_cell){
        for (int j : in_neighbour){
            double dist2 = 0;
            for (int d = 0; d < ndim_; d++){
                dist_1d = abs(positions(d, i-1) - positions(d, j-1));
                if (pbc_ and dist_1d > box_ / 2){
                    dist_1d = box_ - dist_1d;
                }
                dist2 += dist_1d * dist_1d;
            }
            if (dist2 < rc2){
                conn_mat(i-1, j-1) = 1;
            }
        }}
    }
}


void CellList3D::fill_neighbour_indices_1d(Neighbours& neighbours, int cell_idx_1d){
    neighbours.push_back( Indices{ cell_idx_1d } );
    if (sc_ == 1){
        return;
    } else if (sc_ == 2){
        neighbours.back().push_back(1 - cell_idx_1d);
    } else {
        if (cell_idx_1d == 0){
            if (pbc_) {
                neighbours.back().push_back(sc_ - 1);
                neighbours.back().push_back(1);
            } else {
                neighbours.back().push_back(1);
            }
        } else if (cell_idx_1d == sc_ - 1){
            if (pbc_) {
                neighbours.back().push_back(cell_idx_1d - 1);
                neighbours.back().push_back(0);
            } else {
                neighbours.back().push_back(cell_idx_1d - 1);
            }
        } else {
            neighbours.back().push_back(cell_idx_1d + 1);
            neighbours.back().push_back(cell_idx_1d - 1);
        }
    }
}


Indices3D CellList3D::get_neighbours_indices(const Index3D& cell_idx){
    /*
     * For a cell indexed as (a, b, c), find all its neighbours
     * The cell itself is included as a "neighbour"
     */
    Neighbours neighbours; // (dimension, neighbour_num)
    for (int d = 0; d < ndim_; d++){
        fill_neighbour_indices_1d(neighbours, cell_idx[d]);
    }
    return product_3d(neighbours[0], neighbours[1], neighbours[2]);
}


Conn CellList3D::get_conn(Coord3D& positions){
    Conn connections{};
    for (int i = 0; i < positions.cols(); i++){
        connections.push_back(vector<int>{});
    }

    CellIndices3D ci = get_cell_indices(positions);
    double rc2 = rc_ * rc_;
    int max_cell_idx = pow(sc_, ndim_);

    #pragma omp parallel for
    for (int x = 0; x < max_cell_idx; x++){ 
        int cursor = 0;
        double dist_1d = 0;
        Indices3D neighbours;
        Indices in_neighbour;
        Indices in_cell;
        Index3D cell_idx = unravel_index(x, head_shape_);

        // ignore empty cells
        if (chead_[cell_idx] == 0) continue;

        // Collecting particle indices in current cell
        in_cell.push_back(chead_[cell_idx]);
        cursor = chead_[cell_idx];
        while (clist_[cursor] > 0) {
            cursor = clist_[cursor];
            in_cell.push_back(cursor);
        }

        // Collecting particle indices in neighbour cells
        neighbours = get_neighbours_indices(cell_idx);

        for (auto& nc : neighbours){ // nc -> neighbour_cell
            if (chead_[nc] == 0) continue;
            in_neighbour.push_back(chead_[nc]);
            cursor = chead_[nc];
            while (clist_[cursor] > 0){
                cursor = clist_[cursor];
                in_neighbour.push_back(cursor);
            }
        }

        for (int i : in_cell){
        for (int j : in_neighbour){
            double dist2 = 0;
            for (int d = 0; d < ndim_; d++){
                dist_1d = abs(positions(d, i - 1) - positions(d, j - 1));
                if (pbc_ and (dist_1d > (box_ / 2))){
                    dist_1d = box_ - dist_1d;
                }
                dist2 += dist_1d * dist_1d;
            }
            if (dist2 < rc2){
                connections[i - 1].push_back(j - 1);
            }
        }}
    }
    return connections;
}


CellList2D::CellList2D(double r_cut, double box, bool pbc):
    rc(r_cut), box(box), pbc(pbc) {
    ndim = 2;
    size = 0;
    sc = floor(box / rc / 2);
    for (int d = 0; d < ndim; d++){
        head_shape_[d] = sc;
    }
}


void CellList2D::refill(){
    /*
    * Filling chead and clist with zeros
    */
    chead.clear();
    Indices sc_range;

    for (int i=0; i < sc; i++){
        sc_range.push_back(i);
    }

    for (auto item : product_2d(sc_range)){
        chead[item] = 0;
    }

    clist.clear();
    for (int i = 0; i < size + 1; i++){
        clist.push_back(0);
    }
}


Indices2D CellList2D::get_neighbours_indices(Index2D cell_idx){
    /*
     * For a cell indexed as (a, b, c), find all its neighbours
     * The cell itself is included as a "neighbour"
     * TODO: FIX The Cases where Thera are only 1 and 2 cells
     */
    Neighbours neighbours; // (dimension, neighbour_num)

    for (int d = 0; d < ndim; d++){

        neighbours.push_back( Indices{cell_idx[d]} );

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

    return product_2d(neighbours[0], neighbours[1]);
}


void CellList2D::build(Coord2D& positions){
    size = positions.cols();
    CellIndices2D ci(ndim, size);  // discretise positions into cell indices
    ci << floor(positions.array() / box * sc).cast<int>();

    refill();

    Index2D cell_idx;

    for (int i=1; i < size + 1; i++){
        for (int d = 0; d < ndim; d++){
            cell_idx[d] = ci(d, i-1);
        }
        clist[i] = chead[cell_idx];
        chead[cell_idx] = i;
    }
}


void CellList2D::update_sc(int new_sc){
    this->sc = new_sc;
    for (int d = 0; d < ndim; d++){
        head_shape_[d] = sc;
    }
}


void CellList2D::get_dmat(Coord2D& positions, DistMat& dist_mat){
    dist_mat.setConstant(-1);

    double rc2 = rc * rc;

    CellIndices2D ci(ndim, size);
    ci << floor(positions.array() / box * sc).cast<int>();


    int max_cell_idx = pow(sc, ndim);
    // itering over all cells

    #pragma omp parallel for
    for (int x=0; x < max_cell_idx; x++){ 

        int cursor = 0;

        Indices in_neighbour;
        Indices in_cell;
        Index2D cell_idx;
        Indices2D neighbours;

        in_cell.clear();
        in_neighbour.clear();
        cell_idx = unravel_index(x, head_shape_);

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

        for (int i : in_cell){
        for (int j : in_neighbour){
            double dist2 = 0;
            double dist_1d = 0;
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
}


void CellList2D::get_cmat(Coord2D& positions, ConnMat& conn_mat){
    conn_mat.setZero();

    double rc2 = rc * rc;

    CellIndices2D ci(ndim, size);
    ci << floor(positions.array() / box * sc).cast<int>();


    int max_cell_idx = pow(sc, ndim);
    // itering over all cells

    #pragma omp parallel for
    for (int x=0; x < max_cell_idx; x++){ 

        int cursor = 0;

        Indices in_neighbour;
        Indices in_cell;
        Index2D cell_idx;
        Indices2D neighbours;

        in_cell.clear();
        in_neighbour.clear();
        cell_idx = unravel_index(x, head_shape_);

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

        for (int i : in_cell){
        for (int j : in_neighbour){
            double dist2 = 0;
            double dist_1d = 0;
            for (int d = 0; d < ndim; d++){
                dist_1d = abs(positions(d, i-1) - positions(d, j-1));
                if (pbc and dist_1d > box / 2){
                    dist_1d = box - dist_1d;
                }
                dist2 += dist_1d * dist_1d;
            }
            if (dist2 < rc2){
                conn_mat(i-1, j-1) = 1;
            }
        }}
    }
}

Conn CellList2D::get_conn(Coord2D& positions){
    Conn connections;
    for (int i = 0; i < positions.cols(); i++){
        connections.push_back(vector<int>{});
    }

    double rc2 = rc * rc;

    CellIndices2D ci(ndim, size);
    ci << floor(positions.array() / box * sc).cast<int>();


    int max_cell_idx = pow(sc, ndim);
    // itering over all cells

    #pragma omp parallel for
    for (int x=0; x < max_cell_idx; x++){ 

        int cursor = 0;

        Indices in_neighbour;
        Indices in_cell;
        Index2D cell_idx;
        Indices2D neighbours;

        in_cell.clear();
        in_neighbour.clear();
        cell_idx = unravel_index(x, head_shape_);

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

        for (int i : in_cell){
        for (int j : in_neighbour){
            double dist2 = 0;
            double dist_1d = 0;
            for (int d = 0; d < ndim; d++){
                dist_1d = abs(positions(d, i-1) - positions(d, j-1));
                if (pbc and dist_1d > box / 2){
                    dist_1d = box - dist_1d;
                }
                dist2 += dist_1d * dist_1d;
            }
            if (dist2 < rc2){
                connections[i - 1].push_back(j - 1);
            }
        }}
    }
    return connections;
}
