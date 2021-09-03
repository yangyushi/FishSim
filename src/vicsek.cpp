#include "vicsek.hpp"

const double PI = 3.141592653589793238463;
std::default_random_engine generator;


Coord3D xyz_to_sphere(Coord3D& xyz){
    int n = xyz.cols();

    Coord3D sphere(3, xyz.cols());
    Property r(n);

    r = xyz.colwise().norm();

    sphere.row(0) << r;
    sphere.row(1) << xyz.row(1).binaryExpr(
            xyz.row(0), [] (double a, double b) { return atan2(a,b);}
            ); // azimuth, phi

    sphere.row(2) << asin(xyz.row(2).array() / r);

    return sphere;
}


Coord3D sphere_to_xyz(Coord3D& sphere){
    int n = sphere.cols();
    Coord3D xyz(3, n);
    Property r_xy(n);

    r_xy << sphere.array().row(0) * cos(sphere.row(2).array());

    xyz.row(0) << r_xy * cos(sphere.row(1).array());
    xyz.row(1) << r_xy * sin(sphere.row(1).array());
    xyz.row(2) << sphere.array().row(0) * sin(sphere.array().row(2));

    return xyz;
}


Property xy_to_phi(Coord2D& xy){
    Property phi{1, xy.cols()};

    for (int i = 0; i < phi.cols(); i++){
        phi(0, i) = atan2(xy(1, i), xy(0, i));
    }

    return phi;
}


Coord2D phi_to_xy(Property& phi){
    Coord2D xy{2, phi.cols()};
    xy.row(0) = cos(phi);
    xy.row(1) = sin(phi);
    return xy;
}


Coord2D phi_to_xy(Property& phi, double spd){
    Coord2D xy{2, phi.cols()};
    xy.row(0) = cos(phi) * spd;
    xy.row(1) = sin(phi) * spd;
    return xy;
}


Vicsek3D::Vicsek3D(int n, double r, double eta, double v0)
    : noise_(eta), speed_(v0), rc_(r), verlet_list_(r, r + v0), n_(n),
    positions_(3, n), velocities_(3, n) {
        positions_.setRandom(3, n);
        Property vz{n}, vphi{n}, vrxy{n};
        vz.setRandom(); // vz ~ U(-1, 1)
        vrxy = sqrt(1 - vz.pow(2));
        vphi.setRandom();
        vphi *= PI;  // vphi ~ U(-pi, pi)
        velocities_.row(0) << vrxy * cos(vphi);
        velocities_.row(1) << vrxy * sin(vphi);
        velocities_.row(2) << vz;
        normalise(velocities_, speed_);
    }


void Vicsek3D::load_positions(Coord3D positions) {
    this->positions_ = positions;
    verlet_list_.build(this->positions_);
};


void Vicsek3D::load_velocities(Coord3D velocities) {
    this->velocities_ = velocities;
    normalise(this->velocities_, speed_);
};


void Vicsek3D::rotate_noise(Coord3D& noise_xyz){
    /*
    * Rotate nosie pointing at (0, 0, 1) to direction xyz 
    * and then add noise_xyz to velocities
    * here velocities have unit norm
    */
    RotMat F, G, R;
    Vec3D A, B, u, v, w;
    double x{0}, y{0}, z{0}, rxy{0};
    for (int i = 0; i < n_; i++){
        x = velocities_(0, i);
        y = velocities_(1, i);
        z = velocities_(2, i);
        rxy = sqrt(x * x + y * y);
        A << 0, 0, 1;
        B << x, y, z;
        if ((A - B).norm() < 1e-10){
            velocities_.col(i) = noise_xyz.col(i);
        } else if ((A + B).norm() < 1e-10){  // B close to (0, 0, -1)
            velocities_.col(i) << noise_xyz(0, i), noise_xyz(1, i), -noise_xyz(2, i);
        } else {
            u = A;
            v << x / rxy, y / rxy, 0;
            w = B.cross(A);  // w = B x A
            G <<   z, -rxy, 0,
                 rxy,    z, 0,
                   0,    0, 1;
            F << u.transpose(), v.transpose(), w.transpose();
            R = F.inverse() * (G * F);
            velocities_.col(i) = (R * noise_xyz.col(i)) * speed_;
        }
    }
}


void Vicsek3D::rotate_noise_xyz(Coord3D& noise_xyz){
    /*
    * Rotate nosie pointing at (0, 0, 1) to direction xyz 
    * and then add noise_xyz to velocities
    * here velocities have unit norm
    *   R1 -- rotate (0, 0, 1) to (0, 0, 1)
    *   R2 -- rotate (0, 0, 1) to (0, 0, 1) in uvw coordinate
    *   R3 -- elevate (0, 0, 1) in uvw coordinate
    *   uvw -- basis for uvw coordiante
    */
    RotMat R, R1, R2, R3, uvw;
    Vec3D v;
    R1 <<  cos(PI/2),   0, sin(PI/2),
                   0,   1,         0,
          -sin(PI/2),   0, cos(PI/2);
    double phi, theta;
    for (int i = 0; i < n_; i++){
        v = velocities_.col(i);
        phi = atan2(v[1], v[0]);
        theta = asin(v[2] / v.norm());
        uvw << cos(phi), -sin(phi), 0,
               sin(phi),  cos(phi), 0,
                      0,         0, 1;
        R2 << cos(phi), -sin(phi), 0,
              sin(phi), cos(phi), 0,
                      0,        0, 1;
        R3 << cos(-theta), 0, sin(-theta),
                        0, 1,           0,
             -sin(-theta), 0, cos(-theta);
        R = uvw * R3 * R2 * uvw.inverse() * R1;
        velocities_.col(i) = R * noise_xyz.col(i) * speed_;
    }
}


void Vicsek3D::rotate_noise_fast(Coord3D& noise_xyz){
    /*
    * Rotate nosie pointing at (0, 0, 1) to direction xyz 
    * and then add noise_xyz to velocities
    * here velocities have unit norm
    */

    RotMat R1, R2, R;
    Vec3D A, B, v;
    double c;

    for (int i = 0; i < n_; i++){
        A << 0, 0, 1;
        B << velocities_.col(i);
        v = A.cross(B);  // w = A cross B
        c = A.dot(B);  // A dot B
        R1 << 0, -v[2], v[1], // the skew-symmetric matrix
              v[2], 0, -v[0],
             -v[1], v[0], 0;
        R << 1, 0, 0,
             0, 1, 0,
             0, 0, 1; 
        R = R + R1 + (R1 * R1) / (1 + c);
        velocities_.col(i) = (R * noise_xyz.col(i)) * speed_;
    }
}


void Vicsek3D::add_noise(){
    Property noise_rxy{1, n_}, noise_phi{1, n_};
    Coord3D noise_xyz{3, n_};

    std::uniform_real_distribution<> dist_phi(-PI, PI);
    std::uniform_real_distribution<> dist_z(1 - 2 * noise_, 1);

    auto rand_phi = [&] (double) {return dist_phi(generator);};
    auto rand_z = [&] (double) {return dist_z(generator);};

    noise_phi.row(0) = Eigen::ArrayXd::NullaryExpr(n_, rand_phi); // phi ~ U(-PI, PI)
    noise_xyz.row(2) = Eigen::ArrayXd::NullaryExpr(n_, rand_z); //  z ~ U(1 - 2 * noise, 1)

    noise_rxy = sqrt(1 - noise_xyz.array().row(2).pow(2));
    noise_xyz.row(0) = noise_rxy * cos(noise_phi);
    noise_xyz.row(1) = noise_rxy * sin(noise_phi);

    rotate_noise(noise_xyz);
}


void Vicsek3D::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);
    connections_ = verlet_list_.get_conn(positions_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
}


void Vicsek3D::move_no_nl(){
    connections_ = get_connections(positions_, rc_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
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
        velocities_.setRandom(3, n);
        normalise(velocities_, speed_);
    }


void Vicsek3DPBC::move(bool rebuild){
    if (rebuild) {
        cell_list_.build(positions_);
    }
    connections_.clear();
    connections_ = cell_list_.get_conn(positions_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBC::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
    fix_positions();
}


Vicsek3DPBCInertia::Vicsek3DPBCInertia(
    int n, double r, double eta, double box, double v0, double alpha
    ) : Vicsek3DPBC{n, r, eta, box, v0}, alpha_{alpha}, old_velocities_{3, n}
    {}


void Vicsek3DPBCInertia::move(bool rebuild){
    if (rebuild){
        cell_list_.build(positions_);
    }
    connections_ = cell_list_.get_conn(positions_);
    old_velocities_ << velocities_;
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
    velocities_ = (old_velocities_ * alpha_ + velocities_ * (1 - alpha_));
    normalise(velocities_, speed_);
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBCInertia::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    old_velocities_ << velocities_;
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
    velocities_ = (old_velocities_ * alpha_ + velocities_ * (1 - alpha_));
    normalise(velocities_, speed_);
    positions_ += velocities_;
    fix_positions();
}


Vicsek3DPBCInertiaAF::Vicsek3DPBCInertiaAF(
    int n, double r, double eta, double box, double v0, double alpha
    ) : Vicsek3DPBCInertia{n, r, eta, box, v0, alpha} {}


void Vicsek3DPBCInertiaAF::align_af(){
    Coord3D new_velocities{3, n_};
    new_velocities.setZero();
    for (int i = 0; i < n_; i++){
        for (auto j : connections_[i]){
            if (i == j){
                new_velocities.col(i) += velocities_.col(j);
            } else {
                new_velocities.col(i) -= velocities_.col(j);
            }
        }
    }
    velocities_ = new_velocities;
}


void Vicsek3DPBCInertiaAF::move(bool rebuild){
    if (rebuild){
        cell_list_.build(positions_);
    }
    connections_ = cell_list_.get_conn(positions_);
    old_velocities_ << velocities_;
    align_af();
    normalise(velocities_);
    add_noise();
    velocities_ = (old_velocities_ * alpha_ + velocities_ * (1 - alpha_));
    normalise(velocities_, speed_);
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBCInertiaAF::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    old_velocities_ << velocities_;
    align_af();
    normalise(velocities_);
    add_noise();
    velocities_ = (old_velocities_ * alpha_ + velocities_ * (1 - alpha_));
    normalise(velocities_, speed_);
    positions_ += velocities_;
    fix_positions();
}


Attanasi2014PCB::Attanasi2014PCB(
    int n, double r, double eta, double v0, double beta
    ) : Vicsek3D{n, r, eta, v0}, beta_{beta} {}


void Attanasi2014PCB::harmonic_align(){
    Coord3D new_velocities{3, n_};
    new_velocities.setZero();
    for (int i = 0; i < n_; i++){
        for (auto j : connections_[i]){
            new_velocities.col(i) += velocities_.col(j);
        }
        new_velocities.col(i) -= beta_ * positions_.col(i);
    }
    velocities_ = new_velocities;
}


void Attanasi2014PCB::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);
    connections_ = verlet_list_.get_conn(positions_);
    harmonic_align();
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
}


void Attanasi2014PCB::move_no_nl(){
    connections_ = get_connections(positions_, rc_);
    harmonic_align();
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
}


Vicsek3DPBCVN::Vicsek3DPBCVN(
        int n, double r, double eta, double box, double v0
    ) : Vicsek3DPBC(n, r, eta, box, v0) {}


void Vicsek3DPBCVN::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);
    vicsek_align_vn(velocities_, connections_, noise_);
    positions_ += velocities_;
    fix_positions();
}


void Vicsek3DPBCVN::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    vicsek_align_vn(velocities_, connections_, noise_);
    positions_ += velocities_;
    fix_positions();
}


Vicsek2D::Vicsek2D(int n, double r, double eta, double v0)
    : noise_(eta), speed_(v0), rc_(r), verlet_list_(r, r + v0), n_(n),
      positions_(2, n), velocities_(2, n) {
          positions_.setRandom();
          Property vr{n}, vphi{n};
          vr.setConstant(speed_);
          vphi.setRandom();
          vphi *= PI;  // vphi ~ U(-pi, pi)
          velocities_.row(0) << vr * cos(vphi);
          velocities_.row(1) << vr * sin(vphi);
    }


void Vicsek2D::add_noise(){
    Property noise_phi{1, n_};
    Property phi = xy_to_phi(velocities_);
    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= PI * noise_; // phi ~ U(-PI * eta, PI * eta)
    phi += noise_phi;
    velocities_.row(0) = speed_ * cos(phi);
    velocities_.row(1) = speed_ * sin(phi);
}


void Vicsek2D::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);
    connections_ = verlet_list_.get_conn(positions_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);  // speed = 1
    add_noise();  // speed = speed
    positions_ += velocities_;
}
 

void Vicsek2D::move_no_nl(){
    connections_ = get_connections(positions_, rc_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
}


Vicsek2DPBC::Vicsek2DPBC(int n, double r, double eta, double box, double v0)
    : Vicsek2D{n, r, eta, v0}, box_(box), cell_list_(r, box, true){
        double L_per_particle = pow(box * box/ n, 0.5);
        if (L_per_particle > r){
            int sc = floor(box / L_per_particle / 2);
            cell_list_.update_sc(sc);
        }
        positions_.setRandom(2, n);
        positions_ = (positions_.array() + 1) / 2 * box;
        velocities_.setRandom(2, n);
        normalise(velocities_, v0);
    }


void Vicsek2DPBC::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);  // speed = 1
    add_noise();  // speed = speed
    positions_ += velocities_;
    fix_positions();
}


void Vicsek2DPBC::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
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
    vicsek_align_vn(velocities_, connections_, noise_);
    normalise(velocities_, speed_);
    positions_ += velocities_;
    fix_positions();
}


void Vicsek2DPBCVN::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    vicsek_align_vn(velocities_, connections_, noise_);
    normalise(velocities_, speed_);
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


void Vicsek2DPBCVNCO::update_velocity(){
    /*
     * See chatelEPJ2008 for equations
     */
    Property noise_phi{1, n_};
    Coord2D noise_xy{2, n_};
    Coord2D new_velocities{2, n_};
    new_velocities.setZero();
    Vec2D shift;
    Vec2D v_force;
    Vec2D v_avoid;
    double dij{0};
    double n_avoid{0};
    double a = alpha_ / speed_;

    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= PI; // phi ~ U(-PI, PI)

    noise_xy.row(0) = cos(noise_phi) * noise_; // -> eta * xi_i^t
    noise_xy.row(1) = sin(noise_phi) * noise_;

    for (int i = 0; i < n_; i++){
        v_force << 0, 0;
        v_avoid << 0, 0;
        n_avoid = 0;
        for (auto j : connections_[i]){
            if (i == j){
                v_force += a * velocities_.col(j);
            } else {
                shift = get_shift(positions_.col(j), positions_.col(i));
                dij = shift.matrix().norm();
                shift /= dij;
                if (dij < rc_) {
                    n_avoid += 1;  // fij = infty
                    v_avoid += shift * -1;
                } else if (dij < ra_) {
                    v_force += a * velocities_.col(j) + beta_ * shift * c_ * (dij - re_);
                } else {
                    v_force += a * velocities_.col(j) + beta_ * shift;  // fij = 1
                }
            }
        }
        if (n_avoid > 0){
            new_velocities.col(i) = v_avoid;
        } else {
            new_velocities.col(i) = v_force + noise_xy.col(i) * connections_[i].size();
        }
    }
    normalise(new_velocities, speed_);
    velocities_ = new_velocities;
}


void Vicsek2DPBCVNCO::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);
    connections_ = cell_list_.get_conn(positions_);
    update_velocity();
    positions_ += velocities_;
    fix_positions();
}


void Vicsek2DPBCVNCO::move_no_nl(){
    connections_ = get_connections_pbc(positions_, rc_, box_);
    update_velocity();
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

