#include "couzin.hpp"


Coord3D random_pos_no_overlap(double d, int n){
    std::vector<Vec3D> positions;
    Coord3D coordinates_no_overlap {3, n};
    double D = d * pow((double) 10.0 / 1.0 * n, 1.0 / 3.0);
    double pz, pxy, theta, r;
    Vec3D p;
    for (int i = 0; i < n; i++){
        bool overlap = true;
        while (overlap){
            overlap = false;
            pz = (double) 2 * rand() / (RAND_MAX) - 1;  // U(-1, 1)
            theta = (double) rand() / (RAND_MAX) * M_PI * 2;  // U(0, 2*pi)
            r = (double) rand() / (RAND_MAX) * D / 2;  // U(0, D)
            pxy = sqrt(1 - pow(pz, 2));
            p << pxy * cos(theta) * r, pxy * sin(theta) * r , pz * r;
            for (auto p2 : positions){
                double dist = (p2 - p).norm();
                if (dist < d){
                    overlap = true;
                }
            }
        }
        positions.push_back(p);
    }

    for (int i = 0; i < n; i++){
        coordinates_no_overlap.col(i) = positions[i];
    }

    return coordinates_no_overlap;
}


Couzin3D::Couzin3D(
    int n, double rr, double ro, double ra, double ap,
    double eta, double v0, double vt, double dt
) :
    AVS3D{n, eta}, verlet_list_{ra, ra * 3}, speed_(v0),
    r_repel_(rr), r_align_(ro), r_attract_(ra),
    a_percept_(ap), v_turn_(vt), dt_(dt),
    positions_{3, n}, velocities_{3, n}
{
    velocities_ = orientations_;
    positions_ = random_pos_no_overlap(rr, n);
}


void Couzin3D::load_positions(Coord3D positions){
    this->positions_ = positions;
    verlet_list_.build(this->positions_);
}


bool Couzin3D::is_visible(int i, int j){
    Vec3D v1, v2;
    v1 = velocities_.col(i);
    v2 = positions_.col(j);
    double aij = abs(get_angle(v1, v2));
    if (aij < a_percept_){
        return true;
    }
    return false;
}


void Couzin3D::move(bool rebuild){
    if (rebuild) {
        verlet_list_.build(positions_);
        std::vector<Conn> conn_list = verlet_list_.get_neighbour_conn_slow(
            positions_, std::vector<double>{0, r_repel_, r_align_, r_attract_}
        );
        conn_repel_ = conn_list[0];
        conn_align_ = conn_list[1];
        conn_attract_ = conn_list[2];
    }

    Vec3D shift_ij, oi, vj, oi_attr_, oi_align_;
    double dij = 0;

    for (int i = 0; i < n_; i++){  // determine the target orientations
        if (conn_repel_[i].size() > 0){  // repel dominate
            oi << 0, 0, 0;
            for (auto j : conn_repel_[i]){
                shift_ij = positions_.col(j) - positions_.col(i);
                dij = shift_ij.norm();
                oi -= shift_ij / dij;
            }
            orientations_.col(i) = oi / oi.norm();
            continue;
        }

        oi_align_ << velocities_.col(i);
        for (auto j : conn_align_[i]){
            if (is_visible(i, j)) {
                vj = velocities_.col(j);
                oi_align_ += vj;
            }
        }
        oi_align_ = oi_align_ / oi_align_.norm();

        bool see_neighbour = false;
        oi_attr_ << 0, 0, 0;
        for (auto j : conn_attract_[i]){
            if (is_visible(i, j)) {
                shift_ij = positions_.col(j) - positions_.col(i);
                dij = shift_ij.norm();
                oi_attr_ += shift_ij / dij;
                see_neighbour = true;
            }
        }
        if (see_neighbour){
            oi_attr_ = oi_attr_ / oi_attr_.norm();
        }

        oi = oi_align_ + oi_attr_;
        orientations_.col(i) = oi / oi.norm();
    }

    this->add_noise();  // rotate orientations_ randomly

    // try to align velocities to gargeted directions
    for (int i = 0; i < n_; i++){
        Vec3D vi_old = velocities_.col(i);
        Vec3D vi_target = orientations_.col(i);
        RotMat R = get_rotation_matrix(vi_old, vi_target, v_turn_ * dt_);
        velocities_.col(i) << R * vi_old;
    }
    positions_ += velocities_ * speed_ * dt_;
}


void Couzin3D::evolve(int steps, bool rebuild){
    for (int i = 0; i < steps; i++){
        this->move(rebuild);
    }
}


CouzinTank3D::CouzinTank3D(
            int n,
            double rr, double ro, double ra,  // repel, align, attract ranges
            double ap,  // angle of perception
            double eta, // noise
            double v0,  // speed
            double vt,  // turning rate
            double dt,  // delta time
            double c,   // shape parameter for tank
            double h,   // height of tank
            double kw,   // strength of wall interaction
            bool align_cap,   // true: align with wall; false: reflect from wall
            bool align_base,   // true: align with wall; false: reflect from wall
            double g   // gravity strength
        ) :
    Couzin3D(n, rr, ro, ra, ap, eta, v0, vt, dt),
    tank_{c, h, kw, align_cap, align_base}, gravity_{g}
{
    positions_ = tank_.get_random_positions(n);
    tank_.fix_orientations(positions_, orientations_, dt_);
}


void CouzinTank3D::move_in_tank(bool rebuild){
    if (rebuild) {
        verlet_list_.build(positions_);
        std::vector<Conn> conn_list = verlet_list_.get_neighbour_conn_slow(
            positions_, std::vector<double>{0, r_repel_, r_align_, r_attract_}
        );
        conn_repel_ = conn_list[0];
        conn_align_ = conn_list[1];
        conn_attract_ = conn_list[2];
    }

    Vec3D shift_ij, oi, vj, oi_attr_, oi_align_;
    double dij = 0;

    for (int i = 0; i < n_; i++){  // determine the target orientations
        if (conn_repel_[i].size() > 0){  // repel dominate
            oi << 0, 0, 0;
            for (auto j : conn_repel_[i]){
                shift_ij = positions_.col(j) - positions_.col(i);
                dij = shift_ij.norm();
                oi -= shift_ij / dij;
            }
            orientations_.col(i) = oi / oi.norm();
            continue;
        }

        oi_align_ << velocities_.col(i);
        for (auto j : conn_align_[i]){
            if (is_visible(i, j)) {
                vj = velocities_.col(j);
                oi_align_ += vj;
            }
        }
        oi_align_ = oi_align_ / oi_align_.norm();

        bool see_neighbour = false;
        oi_attr_ << 0, 0, 0;
        for (auto j : conn_attract_[i]){
            if (is_visible(i, j)) {
                shift_ij = positions_.col(j) - positions_.col(i);
                dij = shift_ij.norm();
                oi_attr_ += shift_ij / dij;
                see_neighbour = true;
            }
        }
        if (see_neighbour){
            oi_attr_ = oi_attr_ / oi_attr_.norm();
        }

        oi = oi_align_ + oi_attr_;
        orientations_.col(i) = oi / oi.norm();
    }

    this->add_noise();  // rotate orientations_ randomly


    // try to align velocities to gargeted directions
    for (int i = 0; i < n_; i++){
        Vec3D vi_old = velocities_.col(i);
        Vec3D vi_target = orientations_.col(i);
        RotMat R = get_rotation_matrix(vi_old, vi_target, v_turn_ * dt_);
        velocities_.col(i) << R * vi_old;
    }

    gravity_.fix_orientations(positions_, velocities_, dt_);
    tank_.fix_orientations(positions_, velocities_, dt_);

    positions_ += velocities_ * speed_ * dt_;
    tank_.fix_positions(positions_, dt_);
}

