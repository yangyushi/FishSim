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
    AVS3D{n, eta, v0}, verlet_list_{ra, ra * 3},
    r_repel_(rr), r_align_(ro), r_attract_(ra),
    a_percept_(ap), v_turn_(vt), dt_(dt),
    positions_{3, n},
    velocities_real_{3, n}
{
    Property vz{n}, vphi{n}, vrxy{n};
    vz.setRandom(); // vz ~ U(-1, 1)
    vrxy = sqrt(1 - vz.pow(2));
    vphi.setRandom();
    vphi *= M_PI;  // vphi ~ U(-pi, pi)
    velocities_real_.row(0) << vrxy * cos(vphi);
    velocities_real_.row(1) << vrxy * sin(vphi);
    velocities_real_.row(2) << vz;
    normalise(velocities_real_);
    normalise(velocities_);

    positions_ = random_pos_no_overlap(rr, n);
}


void Couzin3D::load_positions(Coord3D positions){
    this->positions_ = positions;
    verlet_list_.build(this->positions_);
}


bool Couzin3D::is_visible(int i, int j){
    Vec3D v1, v2;
    v1 = velocities_real_.col(i);
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
            velocities_.col(i) = oi / oi.norm();
            continue;
        }

        oi_align_ << velocities_real_.col(i);
        for (auto j : conn_align_[i]){
            if (is_visible(i, j)) {
                vj = velocities_real_.col(j);
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
        velocities_.col(i) = oi / oi.norm();
    }

    this->add_noise();  // rotate velocities_ randomly according to the noise
    normalise(velocities_);

    // try to align velocities to gargeted directions
    for (int i = 0; i < n_; i++){
        Vec3D vi_old = velocities_real_.col(i);
        Vec3D vi_target = velocities_.col(i);
        RotMat R = get_rotation_matrix(vi_old, vi_target, v_turn_ * dt_);
        velocities_real_.col(i) << R * vi_old;
    }
    positions_ += velocities_real_ * speed_ * dt_;
}


void Couzin3D::evolve(int steps, bool rebuild){
    for (int i = 0; i < steps; i++){
        this->move(rebuild);
    }
}



