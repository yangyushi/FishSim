#include "vicsek.hpp"


Vicsek3D::Vicsek3D(int n, double r, double eta, double v0)
    : AVS3D{n, eta}, speed_(v0), rc_(r),
    verlet_list_{r, r + v0}, positions_(3, n), velocities_(3, n) {
        positions_.setRandom(3, n);
        this->update_velocity();
    }


void Vicsek3D::load_positions(Coord3D positions) {
    this->positions_ = positions;
    verlet_list_.build(this->positions_);
};

void Vicsek3D::load_velocities(Coord3D velocities) {
    this->velocities_ = velocities;
    this->orientations_ = velocities_;
    normalise(this->orientations_);
}

void Vicsek3D::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);
    connections_ = verlet_list_.get_conn(positions_);
    vicsek_align(orientations_, connections_);
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
}

void Vicsek3D::move_no_nl(){
    connections_ = get_connections(positions_, rc_);
    vicsek_align(velocities_, connections_);
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
}


Vicsek3DPBC::Vicsek3DPBC(int n, double r, double eta, double box, double v0)
    : Vicsek3D{n, r, eta, v0}, box_{box}, cell_list_{r, box, true}
   {
        double L_per_particle = pow(box * box * box / n, 1.0/3.0);
        if (L_per_particle > r){
            int sc = floor(box / L_per_particle / 2);
            cell_list_.update_sc(sc);
        }
        positions_.setRandom(3, n);
        positions_ = (positions_.array() + 1) / 2 * box_;
    }


void Vicsek3DPBC::move(bool rebuild){
    if (rebuild) {
        cell_list_.build(positions_);
    }
    connections_.clear();
    connections_ = cell_list_.get_conn(positions_);
    vicsek_align(velocities_, connections_);
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBC::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    vicsek_align(velocities_, connections_);
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


Vicsek3DPBCInertia::Vicsek3DPBCInertia(
    int n, double r, double eta, double box, double v0, double alpha
    ) : Vicsek3DPBC{n, r, eta, box, v0}, old_orientations_{3, n}, alpha_{alpha}
    {}


void Vicsek3DPBCInertia::move(bool rebuild){
    if (rebuild){
        cell_list_.build(positions_);
    }
    connections_ = cell_list_.get_conn(positions_);
    old_orientations_ << orientations_;
    vicsek_align(orientations_, connections_);
    this->add_noise();
    orientations_ = (old_orientations_ * alpha_ + orientations_ * (1 - alpha_));
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBCInertia::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    old_orientations_ << orientations_;
    vicsek_align(orientations_, connections_);
    add_noise();
    orientations_ = (old_orientations_ * alpha_ + orientations_ * (1 - alpha_));
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


Vicsek3DPBCInertiaAF::Vicsek3DPBCInertiaAF(
    int n, double r, double eta, double box, double v0, double alpha
    ) : Vicsek3DPBCInertia{n, r, eta, box, v0, alpha} {}


void Vicsek3DPBCInertiaAF::align_af(){
    Coord3D new_orientations{3, n_};
    new_orientations.setZero();
    for (int i = 0; i < n_; i++){
        for (auto j : connections_[i]){
            if (i == j){
                new_orientations.col(i) += orientations_.col(j);
            } else {
                new_orientations.col(i) -= orientations_.col(j);
            }
        }
    }
    orientations_ = new_orientations;
}


void Vicsek3DPBCInertiaAF::move(bool rebuild){
    if (rebuild){
        cell_list_.build(positions_);
    }
    connections_ = cell_list_.get_conn(positions_);
    old_orientations_ << orientations_;
    this->align_af();
    this->add_noise();
    orientations_ = (old_orientations_ * alpha_ + orientations_ * (1 - alpha_));
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBCInertiaAF::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    old_orientations_ << orientations_;
    this->align_af();
    this->add_noise();
    orientations_ = (old_orientations_ * alpha_ + orientations_ * (1 - alpha_));
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


Attanasi2014PCB::Attanasi2014PCB(
    int n, double r, double eta, double v0, double beta
    ) : Vicsek3D{n, r, eta, v0}, beta_{beta} {}


void Attanasi2014PCB::harmonic_align(){
    Coord3D new_orientations{3, n_};
    new_orientations.setZero();
    for (int i = 0; i < n_; i++){
        for (auto j : connections_[i]){
            new_orientations.col(i) += orientations_.col(j);
        }
        new_orientations.col(i) -= beta_ * positions_.col(i);
    }
    orientations_ = new_orientations;
}


void Attanasi2014PCB::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);
    connections_ = verlet_list_.get_conn(positions_);
    this->harmonic_align();
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
}


void Attanasi2014PCB::move_no_nl(){
    connections_ = get_connections(positions_, rc_);
    this->harmonic_align();
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
}


Vicsek3DPBCVN::Vicsek3DPBCVN(
        int n, double r, double eta, double box, double v0
    ) : Vicsek3DPBC(n, r, eta, box, v0) {}


void Vicsek3DPBCVN::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);
    vicsek_align_vn(orientations_, connections_, noise_);
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBCVN::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    vicsek_align_vn(orientations_, connections_, noise_);
    this->update_velocity();
    positions_ += velocities_;
    fix_positions();
}


Vicsek2D::Vicsek2D(int n, double r, double eta, double v0)
    : AVS2D{n, eta}, rc_(r), speed_(v0), verlet_list_(r, r + v0),
      positions_(2, n), velocities_(2, n){
          positions_.setRandom();
          this->update_velocity();
    }


void Vicsek2D::update_velocity(){
    velocities_ = orientations_ * speed_;
}


void Vicsek2D::load_velocities(Coord2D v){
    velocities_ = v;
    orientations_ = v;
    normalise(orientations_);
}


void Vicsek2D::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);
    connections_ = verlet_list_.get_conn(positions_);
    vicsek_align(orientations_, connections_);
    this->add_noise();
    this->update_velocity();
    positions_ += velocities_;
}
 

void Vicsek2D::move_no_nl(){
    connections_ = get_connections(positions_, rc_);

    vicsek_align(orientations_, connections_);
    this->add_noise();
    this->update_velocity();

    positions_ += velocities_;
}


Vicsek2DPBC::Vicsek2DPBC(int n, double r, double eta, double box, double v0)
    : Vicsek2D{n, r, eta, v0}, box_(box), cell_list_(r, box, true){
        double L_per_particle = pow(box * box/ n, 0.5);
        if (L_per_particle > r){
            int sc = floor(box / L_per_particle / 2);
            cell_list_.update_sc(sc);
        }
        positions_.setRandom(2, n_);
        positions_ = (positions_.array() + 1) / 2 * box;
        this->update_velocity();
    }


void Vicsek2DPBC::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);

    vicsek_align(orientations_, connections_);
    this->add_noise();
    this->update_velocity();

    positions_ += velocities_;
    fix_positions();
}


void Vicsek2DPBC::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);

    vicsek_align(orientations_, connections_);
    this->add_noise();
    this->update_velocity();

    positions_ += velocities_;
    fix_positions();
}


Vec2D Vicsek2DPBC::get_shift(Vec2D p1, Vec2D p2){
    double shift_1d{0};
    Vec2D shift;
    for (int d = 0; d < 2; d++){
        shift_1d = p1(d) - p2(d);
        if (abs(shift_1d) > box_ / 2){
            if (shift_1d > 0){
                shift_1d = box_ - shift_1d;
            }
            else {
                shift_1d += box_;
            }
        }
        shift(d) = shift_1d;
    }
    return shift;
}


Vicsek2DPBCVN::Vicsek2DPBCVN(int n, double r, double eta, double box, double v0)
    : Vicsek2DPBC(n, r, eta, box, v0) {}


void Vicsek2DPBCVN::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);

    vicsek_align_vn(orientations_, connections_, noise_);
    this->update_velocity();

    positions_ += velocities_;
    fix_positions();
}


void Vicsek2DPBCVN::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);

    vicsek_align_vn(orientations_, connections_, noise_);
    this->update_velocity();

    positions_ += velocities_;
    fix_positions();
}


Vicsek2DPBCVNCO::Vicsek2DPBCVNCO(
        int n, double r, double eta, double box, double v0,
        double a, double b, double ra, double re, double rc)
    : Vicsek2DPBC(n, r, eta, box, v0), alpha_{a}, beta_{b}, ra_{ra}, re_{re}, rc_{rc} {
        c_ = 1 / (ra - re) / 4;  // fij = c * (rij - re)
    }


Vicsek2DPBCVNCO::Vicsek2DPBCVNCO(int n, double r, double eta, double box, double v0)
    : Vicsek2DPBC(n, r, eta, box, v0){
        alpha_ = 1.5;
        beta_ = 1;
        ra_ = 0.8;
        re_ = 0.5; 
        rc_ = 0.2;
        c_ = 1 / (ra_ - re_) / 4;  // fij = c * (rij - re)
    }


void Vicsek2DPBCVNCO::apply_interaction(){
    /*
     * See chatelEPJ2008 for equations
     */
    Property noise_phi{1, n_};
    Coord2D noise_xy{2, n_};
    Coord2D new_orientations{2, n_};
    new_orientations.setZero();
    Vec2D shift;
    Vec2D v_force;
    Vec2D v_avoid;
    double dij{0};
    double n_avoid{0};
    double a = alpha_ / speed_;

    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= M_PI; // phi ~ U(-PI, PI)

    noise_xy.row(0) = cos(noise_phi) * noise_; // -> eta * xi_i^t
    noise_xy.row(1) = sin(noise_phi) * noise_;

    for (int i = 0; i < n_; i++){
        v_force << 0, 0;
        v_avoid << 0, 0;
        n_avoid = 0;
        for (auto j : connections_[i]){
            if (i == j){
                v_force += a * orientations_.col(j);
            } else {
                shift = get_shift(positions_.col(j), positions_.col(i));
                dij = shift.matrix().norm();
                shift /= dij;
                if (dij < rc_) {
                    n_avoid += 1;  // fij = infty
                    v_avoid += shift * -1;
                } else if (dij < ra_) {
                    v_force += a * orientations_.col(j) + beta_ * shift * c_ * (dij - re_);
                } else {
                    v_force += a * orientations_.col(j) + beta_ * shift;  // fij = 1
                }
            }
        }
        if (n_avoid > 0){
            new_orientations.col(i) = v_avoid;
        } else {
            new_orientations.col(i) = v_force + noise_xy.col(i) * connections_[i].size();
        }
    }
    orientations_ = new_orientations;
    this->update_velocity();
}


void Vicsek2DPBCVNCO::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);
    this->apply_interaction();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek2DPBCVNCO::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    this->apply_interaction();
    positions_ += velocities_;
    fix_positions();
}


InertialSpin3D::InertialSpin3D(
        int n, double r, double v0, double T, double j, double m, double f
        )
    :mass_{m}, friction_{f}, T_{T}, J_{j}, rc_{r}, verlet_list_{r, r + v0},
     n_{n}, spins_{3, n}, positions_{3, n}, velocities_{3, n}, speed_{v0} {
        positions_.setRandom(3, n);
        velocities_.setZero();
        velocities_.row(0).array() += speed_;
        spins_.setZero(); // particles do not intend to change direction
        imass_ = 1 / mass_;
        v0_sq_inv_ = 1 / v0 / v0;
        dt_ = 0.1 * sqrt(j / m);
}


void InertialSpin3D::move(bool rebuild){
    if (rebuild) {
        verlet_list_.build(positions_);
    }
    connections_ = verlet_list_.get_conn(positions_);
    update_velocity_full();
    update_spin();
    positions_ += velocities_ * dt_;
    normalise(velocities_, speed_);
}


void InertialSpin3D::move_no_nl(){
    connections_ = get_connections(positions_, rc_);
    update_velocity_full();
    update_spin();
    positions_ += velocities_ * dt_;
    normalise(velocities_, speed_);
}


void InertialSpin3D::update_velocity_half(){
    for (int i = 0; i < n_; i++){
        velocities_.col(i) += 0.5 * dt_ * imass_ * spins_.col(i).cross(velocities_.col(i));
    }
}


void InertialSpin3D::update_velocity_full(){
    for (int i = 0; i < n_; i++){
        velocities_.col(i) += dt_ * imass_ * spins_.col(i).cross(velocities_.col(i));
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
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<float> d{0, 1};
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
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d{0, sigma};
    for (int i = 0; i < n_; i++){
        velocity_sum.setZero();
        noise << d(gen), d(gen), d(gen);
        for (int j : connections_[i]){
            velocity_sum += velocities_.col(j);
        }
        spins_.col(i) += velocities_.col(i).cross(
                              J_ * v0_sq_inv_ * velocity_sum
                        ) * dt_\
                      - friction_ * imass_ * spins_.col(i) \
                      + velocities_.col(i).cross(noise) * sqrt(dt_) / speed_;
    }
}


InertialSpin3DTP::InertialSpin3DTP(
        int n, int nc, double v0,
        double T, double j, double m, double f
        )
    : InertialSpin3D{n, 1, v0, T, j, m, f}, nc_{nc} { }


void InertialSpin3DTP::move(){
    connections_ = get_topology_connections(positions_, nc_);
    update_velocity_full();
    update_spin();
    positions_ += velocities_ * dt_;
    normalise(velocities_, speed_);
}


InertialSpin3DPBC::InertialSpin3DPBC(
        int n, double box, double r, double v0,
        double T, double j, double m, double f
        )
    : InertialSpin3D{n, r, v0, T, j, m, f}, box_{box}, cell_list_{r, box, true}{
        double r_box = pow(box * box * box / n, 1.0/3.0);
        if (r_box > r){
            int sc = floor(box / r_box / 2);
            cell_list_.update_sc(sc);
        }
        positions_.setRandom();
        positions_ = (positions_.array() + 1) / 2 * box_;
        velocities_.setRandom();
        normalise(velocities_, speed_);
    }


void InertialSpin3DPBC::move(bool rebuild){
    if (rebuild) {
        cell_list_.build(positions_);
    }
    connections_ = cell_list_.get_conn(positions_);
    update_velocity_full();
    update_spin();
    positions_ += velocities_ * dt_;
    fix_positions();
    normalise(velocities_, speed_);
}


void InertialSpin3DPBC::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    update_velocity_full();
    update_spin();
    positions_ += velocities_ * dt_;
    fix_positions();
    normalise(velocities_, speed_);
}
