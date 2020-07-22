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


VerletList3D::VerletList3D(double r_cut, double r_skin)
    : rc_{r_cut}, rl_{r_skin} {
        rl2_ = rl_ * rl_;
        rc2_ = rc_ * rc_;
    }


void VerletList3D::build(Coord3D& positions){
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
    point_size_ = point_.size();
}


void VerletList3D::get_dmat(Coord3D& positions, DistMat& dist_mat){
    dist_mat.setConstant(-1);
    #pragma omp parallel for
    for (int i = 0; i < point_size_ - 1; i++){
        int p0 = point_[i];
        int idx_j = 0;
        double dist2 = 0;
        for (int j = p0; j < point_[i+1]; j++){
            idx_j = nlist_[j];
            dist2 = (positions.col(i) - positions.col(j)).array().pow(2).sum();
            if (dist2 < rc2_){
                dist_mat(i, idx_j) = sqrt(dist2);
            }
        }
    }
}


void VerletList3D::get_cmat(Coord3D& positions, ConnMat& conn_mat){
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


Conn VerletList3D::get_conn(Coord3D& positions){
    Conn connections;
    for (int i = 0; i < point_size_ - 1; i++){
        connections.push_back(vector<int>{});
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


void VerletList3D::get_cmat_slow(Coord3D& positions, ConnMat& conn_mat){
    conn_mat.setZero();
    int n = positions.cols();
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double dist2 = (positions.col(i) - positions.col(j)).array().pow(2).sum();
            if (dist2 < rc2_){
                conn_mat(i, j) = 1;
            }
        }
    }
}


CellList3D::CellList3D(double r_cut, double box, bool pbc)
    : rc(r_cut), box(box), pbc(pbc) {
    ndim = 3;
    size = 0;
    sc = floor(box / rc / 2);
    for (int d = 0; d < ndim; d++){
        head_shape[d] = sc;
    }
}


void CellList3D::update_sc(int new_sc){
    this->sc = new_sc;
    for (int d = 0; d < ndim; d++){
        head_shape[d] = sc;
    }
}


void CellList3D::refill(){
    /*
    * Filling chead and clist with zeros
    */
    chead.clear();
    Indices sc_range;

    for (int i=0; i < sc; i++){
        sc_range.push_back(i);
    }

    for (auto item : product_3d(sc_range)){
        chead[item] = 0;
    }

    clist.clear();
    for (int i = 0; i < size + 1; i++){
        clist.push_back(0);
    }
}


Indices3D CellList3D::get_neighbours_indices(Index3D cell_idx){
    /*
     * For a cell indexed as (a, b, c), find all its neighbours
     * The cell itself is included as a "neighbour"
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

    return product_3d(neighbours[0], neighbours[1], neighbours[2]);
}


void CellList3D::build(Coord3D& positions){
    size = positions.cols();
    CellIndices3D ci(ndim, size);  // discretise positions into cell indices
    ci << floor(positions.array() / box * sc).cast<int>();

    refill();

    Index3D cell_idx;

    for (int i=1; i < size + 1; i++){
        for (int d = 0; d < ndim; d++){
            cell_idx[d] = ci(d, i-1);
        }
        clist[i] = chead[cell_idx];
        chead[cell_idx] = i;
    }
}


void CellList3D::get_dmat(Coord3D& positions, DistMat& dist_mat){
    dist_mat.setConstant(-1);

    CellIndices3D ci(ndim, size);
    ci << floor(positions.array() / box * sc).cast<int>();

    double rc2 = rc * rc;
    int max_cell_idx = pow(sc, ndim);

    #pragma omp parallel for
    for (int x = 0; x < max_cell_idx; x++){ 

        int cursor = 0;
        double dist_1d = 0;
        Indices3D neighbours;
        Indices in_neighbour;
        Indices in_cell;
        Index3D cell_idx = unravel_index(x, head_shape);

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


void CellList3D::get_cmat(Coord3D& positions, ConnMat& conn_mat){
    conn_mat.setZero();

    CellIndices3D ci(ndim, size);
    ci << floor(positions.array() / box * sc).cast<int>();

    double rc2 = rc * rc;
    int max_cell_idx = pow(sc, ndim);

    #pragma omp parallel for
    for (int x = 0; x < max_cell_idx; x++){ 

        int cursor = 0;
        double dist_1d = 0;
        Indices3D neighbours;
        Indices in_neighbour;
        Indices in_cell;
        Index3D cell_idx = unravel_index(x, head_shape);

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


Conn CellList3D::get_conn(Coord3D& positions){
    Conn connections;
    for (int i = 0; i < positions.cols(); i++){
        connections.push_back(vector<int>{});
    }

    CellIndices3D ci(ndim, size);
    ci << floor(positions.array() / box * sc).cast<int>();

    double rc2 = rc * rc;
    int max_cell_idx = pow(sc, ndim);

    #pragma omp parallel for
    for (int x = 0; x < max_cell_idx; x++){ 
        int cursor = 0;
        double dist_1d = 0;
        Indices3D neighbours;
        Indices in_neighbour;
        Indices in_cell;
        Index3D cell_idx = unravel_index(x, head_shape);

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

        for (auto& nc : neighbours){ // nc -> neighbour_cell
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
            for (int d = 0; d < ndim; d++){
                dist_1d = abs(positions(d, i - 1) - positions(d, j - 1));
                if (pbc and (dist_1d > (box / 2))){
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


CellList2D::CellList2D(double r_cut, double box, bool pbc) :
    rc(r_cut), box(box), pbc(pbc) {
    ndim = 2;
    size = 0;
    sc = floor(box / rc / 2);
    for (int d = 0; d < ndim; d++){
        head_shape[d] = sc;
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
    ci << floor(positions.array() / box * sc);

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
        head_shape[d] = sc;
    }
}


void CellList2D::get_dmat(Coord2D& positions, DistMat& dist_mat){
    dist_mat.setConstant(-1);

    double rc2 = rc * rc;

    CellIndices2D ci(ndim, size);
    ci << floor(positions.array() / box * sc);


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
        cell_idx = unravel_index(x, head_shape);

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
    ci << floor(positions.array() / box * sc);


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
        cell_idx = unravel_index(x, head_shape);

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
    ci << floor(positions.array() / box * sc);


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
        cell_idx = unravel_index(x, head_shape);

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
