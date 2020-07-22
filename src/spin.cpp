#include "spin.h"


InertialSpin3D::InertialSpin3D(int n, double r, double v0,
        double T, double j, double m, double f)
    :speed_{v0}, mass_{m}, friction_{f}, T_{T}, J_{j}, verlet_list_{r, r + v0},
     n_{n}, spins_{3, n}, positions_{3, n}, velocities_{3, n} {
        Property vz{n}, vphi{n}, vrxy{n};
        RotMat R;
        positions_.setRandom(3, n);
        positions_ /= v0;
        velocities_.setZero();
        velocities_.row(0).array() += 1;
        normalise(velocities_, speed_);
        spins_.setZero(); // particles do not intend to change direction
        imass_ = 1 / mass_;
        v0_sq_inv_ = 1 / v0 / v0;
        dt_ = 0.1;
}

void InertialSpin3D::move(bool rebuild){
    if (rebuild) {
        verlet_list_.build(positions_);
    }
    connections_ = verlet_list_.get_conn(positions_);
    add_alignment();
    update_velocity();
    add_noise();
    add_alignment();
    update_velocity();
    positions_ += velocities_ * dt_;
}

void InertialSpin3D::update_velocity(){
    for (int i = 0; i < n_; i++){
        velocities_.col(i) += 0.5 * dt_ * imass_ * spins_.col(i).cross(velocities_.col(i));
    }
}

void InertialSpin3D::add_alignment(){
    Vec3D velocity_sum;
    for (int i = 0; i < n_; i++){
        velocity_sum.setZero();
        for (int j : connections_[i]){
            velocity_sum += velocities_.col(j);
        }
        spins_.col(i) += 0.5 * dt_ * velocities_.col(i).cross(velocity_sum * J_ * v0_sq_inv_);
    }
}

void InertialSpin3D::add_noise(){
    Vec3D noise;
    random_device rd{};
    mt19937 gen{rd()};
    normal_distribution<float> d{0, 1};
    for (int i = 0; i < n_; i++){
        noise << d(gen), d(gen), d(gen);
        spins_.col(i) = exp(-friction_ / mass_) * spins_.col(i) +\
                        sqrt((1 - exp(-2 * friction_ / mass_)) * mass_ * T_) *\
                        velocities_.col(i).cross(noise);
    }
}

void InertialSpin3D::update_spin(){
    Vec3D noise;
    Vec3D velocity_sum;
    double sigma = 2 * 3 * friction_ * T_;
    random_device rd{};
    mt19937 gen{rd()};
    normal_distribution<double> d{0, sigma};
    for (int i = 0; i < n_; i++){
        velocity_sum.setZero();
        noise << d(gen), d(gen), d(gen);
        for (int j : connections_[i]){
            velocity_sum += velocities_.col(j);
        }
        spins_.col(i) = spins_.col(i) \
                      + velocities_.col(i).cross(
                              J_ * v0_sq_inv_ * velocity_sum
                        ) \
                      - friction_ * imass_ * spins_.col(i) \
                      + velocities_.col(i).cross(noise);
    }
}

InertialSpin3DPBC::InertialSpin3DPBC(
        int n, double box, double r, double eta, double v0,
        double j, double m, double f
        )
    : InertialSpin3D{n, r, eta, v0, j, m, f}, cell_list_{r, box, true}{
        double r_box = pow(box * box * box / n, 1.0/3.0);
        if (r_box > r){
            int sc = floor(box / r_box / 2);
            cell_list_.update_sc(sc);
        }
        positions_.setRandom(3, n);
        positions_ = (positions_.array() + 1) / 2 * box_;
    }

void InertialSpin3DPBC::move(bool rebuild){
    if (rebuild) {
        cell_list_.build(positions_);
    }
    connections_ = cell_list_.get_conn(positions_);
    add_alignment();
    update_velocity();
    add_noise();
    add_alignment();
    update_velocity();
    positions_ += velocities_ * dt_;
}
