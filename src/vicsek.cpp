#include "vicsek.h"

default_random_engine generator;

const double PI = 3.141592653589793238463;


Coord3D xyz_to_sphere(Coord3D& xyz){
    int n = xyz.cols();

    Coord3D sphere(3, xyz.cols());
    Property r(n);

    r << sqrt( xyz.row(0).pow(2) +
               xyz.row(1).pow(2) +
               xyz.row(2).pow(2) );

    sphere.row(0) << r;
    sphere.row(1) << xyz.row(1).binaryExpr(
            xyz.row(0), [] (double a, double b) { return atan2(a,b);}
            ); // azimuth, phi

    sphere.row(2) << asin(xyz.row(2) / r);

    return sphere;
}


Coord3D sphere_to_xyz(Coord3D& sphere){
    int n = sphere.cols();
    Coord3D xyz(3, n);
    Property r_xy(n);

    r_xy << sphere.row(0) * cos(sphere.row(2));

    xyz.row(0) << r_xy * cos(sphere.row(1));
    xyz.row(1) << r_xy * sin(sphere.row(1));
    xyz.row(2) << sphere.row(0) * sin(sphere.row(2));

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


void normalise(Coord3D& xyz){
    int i = 0;
    for (auto L : xyz.colwise().norm()){
        xyz.col(i) /= L;
        i++;
    }
}


void normalise(Coord3D& xyz, double spd){
    int i = 0;
    for (auto L : xyz.colwise().norm()){
        xyz.col(i) = xyz.col(i) / L * spd;
        i++;
    }
}


void normalise(Coord2D& xy){
    int i = 0;
    for (auto L : xy.colwise().norm()){
        xy.col(i) /= L;
        i++;
    }
}


void normalise(Coord2D& xy, double spd){
    int i = 0;
    for (auto L : xy.colwise().norm()){
        xy.col(i) = xy.col(i) / L * spd;
        i++;
    }
}


Vicsek3D::Vicsek3D(int n, double r, double eta, double v0)
    : noise_(eta), speed_(v0), n_(n), conn_mat_(n, n),
    verlet_list_(r, r*3), positions_(3, n), velocities_(3, n) {
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


void Vicsek3D::align(){
    Coord3D new_velocities{3, n_};
    new_velocities.setZero();
    for (int i = 0; i < n_; i++){
        for (int j = 0; j < n_; j++){
            if (conn_mat_(i, j) > 0){
                new_velocities.col(i) += velocities_.col(j);
            }
        }
    }
    velocities_ = new_velocities;
}


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
        if ((A - B).matrix().norm() < 1e-10){
            velocities_.col(i) = noise_xyz.col(i);
        } else if ((A + B).matrix().norm() < 1e-10){  // B close to (0, 0, -1)
            velocities_.col(i) << noise_xyz(0, i), noise_xyz(1, i), -noise_xyz(2, i);
        } else {
            u = A;
            v << x / rxy, y / rxy, 0;
            w = B.matrix().cross(A.matrix());  // w = B x A
            G <<   z, -rxy, 0,
                 rxy,    z, 0,
                   0,    0, 1;
            F << u.transpose(), v.transpose(), w.transpose();
            R = F.inverse() * (G * F);
            velocities_.col(i) = (R * noise_xyz.col(i).matrix()) * speed_;
        }
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
        if ((A - B).matrix().norm() < 1e-10){ // B close to (0, 0, 1)
            velocities_.col(i) = noise_xyz.col(i);  // R -> I
        } else if ((A + B).matrix().norm() < 1e-10){  // B close to (0, 0, -1)
            velocities_.col(i) << noise_xyz(0, i),  // R -> rotate to (0, 0, -1)
                                  noise_xyz(1, i),
                                - noise_xyz(2, i);
        } else {
            v = A.matrix().cross(B.matrix());  // w = A cross B
            c = (A * B).sum();  // A dot B
            R1 << 0, -v[2], v[1], // the skew-symmetric matrix
                  v[2], 0, -v[0],
                 -v[1], v[0], 0;
            R << 1, 0, 0,
                 0, 1, 0,
                 0, 0, 1; 
            R = R + R1 + (R1 * R1) / (1 + c);
            velocities_.col(i) = (R * noise_xyz.col(i).matrix()) * speed_;
        }
    }
}


void Vicsek3D::add_noise(){
    uniform_real_distribution<> dist_phi(-PI, PI);
    uniform_real_distribution<> dist_z(1 - 2 * noise_, 1);
    auto rand_phi = [&] (double) {return dist_phi(generator);};
    auto rand_z = [&] (double) {return dist_z(generator);};
    Property noise_rxy{1, n_}, noise_phi{1, n_};
    Coord3D noise_xyz{3, n_};

    //noise_phi.setRandom(); // ~ U(-1, 1)
    //noise_phi *= PI; // phi ~ U(-PI, PI)

    noise_phi.row(0) = Eigen::ArrayXd::NullaryExpr(n_, rand_phi); // phi ~ U(-PI, PI)
    noise_xyz.row(2) = Eigen::ArrayXd::NullaryExpr(n_, rand_z); //  z ~ U(1 - 2 * noise, 1)

    //noise_xyz.row(2).setRandom(); // ~ U(-1, 1)
    //noise_xyz.row(2) = noise_xyz.row(2) * noise_ - noise_ + 1; // z ~ U(1 - 2 * noise, 1)

    noise_rxy = sqrt(1 - noise_xyz.row(2).pow(2));
    noise_xyz.row(0) = noise_rxy * cos(noise_phi);
    noise_xyz.row(1) = noise_rxy * sin(noise_phi);

    rotate_noise_fast(noise_xyz);
}


void Vicsek3D::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);

    verlet_list_.get_cmat(positions_, conn_mat_);

    align();

    normalise(velocities_);

    add_noise();

    positions_ += velocities_;
}


void Vicsek3D::dump(string filename){
    ofstream f;
    f.open(filename, ios::out | ios::app);
    f << n_ << endl;
    f << "id, x, y, z, vx, vy, vz" << endl;
    for (int i = 0; i < n_; i++ ) {
        f << i << " "
            << positions_(0, i)  << " "
            << positions_(1, i)  << " "
            << positions_(2, i)  << " "
            << velocities_(0, i) << " "
            << velocities_(1, i) << " "
            << velocities_(2, i) << endl;
    }
    f.close();
}


void Vicsek3D::load(string filename){
    /*
     * Load the configuration from the last frame in a xyz file
     * The xyz file should be produced by this->dump
     */
    ifstream f;
    string line;
    regex head_pattern{"\\d+"};
    smatch matched;
    int head_lines = 2;
    string num;
    int N = 0;
    int total_frame = 0;
    
    f.open(filename, ios::in);
    while (f) {
        getline(f, line);
        if (regex_match(line, matched, head_pattern)){
            N = stoi(line);
            total_frame += 1;
            for (int i=0; i<N; i++) getline(f, line);
        }
    }
    f.close();
    
    f.open(filename, ios::in);
    for (int i = 0; i < total_frame - 1; i++){
        for (int j = 0; j < N + head_lines; j++){
        getline(f, line);
        }
    }
    
    for (int i = 0; i < N + head_lines; i++){
        getline(f, line);
        
        if (i > 1) {
            istringstream ss(line);
            ss >> num;
            for (int j = 0; j < 3; j++){
                ss >> positions_(j, i - head_lines);
            }
            for (int j = 0; j < 3; j++){
                ss >> velocities_(j, i - head_lines);
            }
        }
    }
    f.close();
}


Vicsek3DPBC::Vicsek3DPBC(int n, double r, double eta, double box, double v0)
    : Vicsek3D{n, r, eta, v0}, box_{box}, cell_list_{r, box, true} {
        positions_.setRandom(3, n);
        positions_ = (positions_ + 1) / 2 * box_;
        velocities_.setRandom(3, n);
        normalise(velocities_, speed_);
    }


void Vicsek3DPBC::move(bool rebuild){
    if (rebuild) {
        cell_list_.build(positions_);
    }
    cell_list_.get_cmat(positions_, conn_mat_);
    align();
    normalise(velocities_);
    add_noise();
    positions_ += velocities_;
    fix_positions();
}


Vicsek3DPBCInertia::Vicsek3DPBCInertia(int n, double r, double eta, double box, double v0, double alpha)
    : Vicsek3DPBC{n, r, eta, box, v0}, alpha_{alpha}, old_velocities_{3, n} {}


void Vicsek3DPBCInertia::move(bool rebuild){
    if (rebuild){
        cell_list_.build(positions_);
    }
    cell_list_.get_cmat(positions_, conn_mat_);
    old_velocities_ << velocities_;
    align();
    normalise(velocities_);
    add_noise();
    velocities_ = (old_velocities_ * alpha_ + velocities_ * (1 - alpha_));
    normalise(velocities_, speed_);
    positions_ += velocities_;
    fix_positions();
}


Attanasi2014PCB::Attanasi2014PCB(int n, double r, double eta, double v0, double beta)
    : Vicsek3D{n, r, eta, v0}, beta_{beta} {}


void Attanasi2014PCB::apply_harmonic_force(){
    for (int i = 0; i < n_; i++){
        velocities_.col(i) -= beta_ * positions_.col(i);
    }
}


void Attanasi2014PCB::move(bool rebuild){
    if (rebuild) verlet_list_.build(positions_);

    verlet_list_.get_cmat(positions_, conn_mat_);

    align();

    apply_harmonic_force();

    normalise(velocities_);

    add_noise();

    positions_ += velocities_;
}


Vicsek3DPBCVN::Vicsek3DPBCVN(int n, double r, double eta, double box, double v0)
    : Vicsek3DPBC(n, r, eta, box, v0) {}


void Vicsek3DPBCVN::noisy_align(){
    /*
     * See ginellEPJ2016 for equations
     */
    Property noise_phi{1, n_};
    Property noise_theta{1, n_};
    PropertyInt neighbour_nums{1, n_};

    Coord3D noise_xyz{3, n_};
    Coord3D new_velocities{3, n_};
    Vec3D v_align;

    neighbour_nums = conn_mat_.colwise().sum(); // -> m_i

    noise_phi.setRandom();   noise_phi   = (noise_phi + 1) * PI; // phi ~ U(0, 2PI)
    noise_theta.setRandom(); noise_theta = asin(noise_theta);    // theta ~ U(0, PI)

    noise_xyz.row(0) = cos(noise_theta) * cos(noise_phi) * speed_ * noise_; // -> eta * xi_i^t
    noise_xyz.row(1) = cos(noise_theta) * sin(noise_phi) * speed_ * noise_;
    noise_xyz.row(2) = sin(noise_theta) * speed_ * noise_;

    for (int i = 0; i < n_; i++){
        v_align << 0, 0, 0;

        for (int j = 0; j < conn_mat_.cols(); j++){  // \sum n_ij^t s_j^t
            if (conn_mat_(i, j) > 0){
                v_align += velocities_.col(j);
            }
        }
        new_velocities.col(i) = v_align + noise_xyz.col(i) * neighbour_nums(0, i);
    }

    normalise(new_velocities, speed_);
    velocities_ = new_velocities;
}


void Vicsek3DPBCVN::move(bool rebuild){
    if (rebuild) cell_list_.build(positions_);

    cell_list_.get_cmat(positions_, conn_mat_);

    noisy_align();

    positions_ += velocities_;

    fix_positions();
}


Vicsek2DPBC::Vicsek2DPBC(int n, double r, double eta, double box, double v0)
    : noise(eta), speed(v0), box(box), n(n), conn_mat(n, n),
      positions(2, n), velocities(2, n), cell_list(r, box, true){
        positions.setRandom(2, n);
        positions = (positions + 1) / 2 * box;
        velocities.setRandom(2, n);
        normalise(velocities, v0);
    }


void Vicsek2DPBC::align(){
    Coord2D new_velocities{2, n};
    Vec2D v_align;

    for (int i = 0; i < this->n; i++){
        v_align << 0, 0;

        for (int j = 0; j < conn_mat.rows(); j++){
            if (conn_mat(i, j) > 0){
                v_align += velocities.col(j);
            }
        }
        new_velocities.col(i) << v_align;
    }

    velocities = new_velocities;
}


void Vicsek2DPBC::add_noise(){
    Property noise_phi{1, this->n};
    Property phi = xy_to_phi(velocities);

    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= PI * this->noise; // phi ~ U(-PI * eta, PI * eta)
    phi += noise_phi;

    velocities.row(0) = speed * cos(phi);
    velocities.row(1) = speed * sin(phi);
}


void Vicsek2DPBC::move(bool rebuild){
    if (rebuild) cell_list.build(positions);
    cell_list.get_cmat(positions, conn_mat);

    align();
    normalise(velocities);  // speed = 1
    add_noise();  // speed = speed

    for (int d = 0; d < 2; d++){
        positions.row(d) += velocities.row(d);
    }
    fix_positions();
}


void Vicsek2DPBC::dump(string filename){
    ofstream f;
    f.open(filename, ios::out | ios::app);
    f << this->n << endl;
    f << "id, x, y, vx, vy" << endl;
    for (int i = 0; i < this->n; i++ ) {
        f << i << " "
          << positions(0, i)  << " " <<  positions(1, i)  << " "
          << velocities(0, i) << " " <<  velocities(1, i) << " " << endl;
    }
    f.close();
}


void Vicsek2DPBC::load(string filename){
    /*
     * Load the configuration from the last frame in a xyz file
     * The xyz file should be produced by this->dump
     */
    ifstream f;
    string line;
    regex head_pattern{"\\d+"};
    smatch matched;
    int head_lines = 2;
    string num;
    int N = 0;
    int total_frame = 0;
    
    f.open(filename, ios::in);
    while (f) {
        getline(f, line);
        if (regex_match(line, matched, head_pattern)){
            N = stoi(line);
            total_frame += 1;
            for (int i=0; i<N; i++) getline(f, line);
        }
    }
    f.close();
    
    f.open(filename, ios::in);
    for (int i = 0; i < total_frame - 1; i++){
        for (int j = 0; j < N + head_lines; j++){
            getline(f, line);
        }
    }
    
    for (int i = 0; i < N + head_lines; i++){
        getline(f, line);
        
        if (i > 1) {
            istringstream ss(line);
            ss >> num;
            for (int j = 0; j < 2; j++){
                ss >> positions(j, i - head_lines);
            }
            for (int j = 0; j < 2; j++){
                ss >> velocities(j, i - head_lines);
            }
        }
    }
    f.close();
}


Vec2D Vicsek2DPBC::get_shift(Vec2D p1, Vec2D p2){
    double shift_1d{0};
    Vec2D shift;

    for (int d = 0; d < 2; d++){
        shift_1d = p1(d) - p2(d);
        if (abs(shift_1d) > this->box / 2){
            if (shift_1d > 0){
                shift_1d = this->box - shift_1d;
            }
            else {
                shift_1d += this->box;
            }
        }
        shift(d) = shift_1d;
    }
    return shift;
}


Vicsek2DPBCVN::Vicsek2DPBCVN(int n, double r, double eta, double box, double v0)
    : Vicsek2DPBC(n, r, eta, box, v0) {}


void Vicsek2DPBCVN::noisy_align(){
    /*
     * See ginellEPJ2016 for equations
     */
    Property noise_phi{1, this->n};
    Coord2D noise_xy{2, this->n};
    PropertyInt neighbour_nums{1, this->n};
    Coord2D new_velocities{2, this->n};
    Vec2D v_align;

    neighbour_nums = conn_mat.colwise().sum(); // -> m_i

    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= PI; // phi ~ U(-PI, PI)

    noise_xy.row(0) = cos(noise_phi) * speed * noise; // -> eta * xi_i^t
    noise_xy.row(1) = sin(noise_phi) * speed * noise;

    for (int i = 0; i < this->n; i++){
        v_align << 0, 0;

        for (int j = 0; j < conn_mat.cols(); j++){  // \sum n_ij^t s_j^t
            if (conn_mat(i, j) > 0){
                v_align += velocities.col(j);
            }
        }
        new_velocities.col(i) = v_align + noise_xy.col(i) * neighbour_nums(0, i);
    }

    normalise(new_velocities, speed);
    velocities = new_velocities;
}


void Vicsek2DPBCVN::move(bool rebuild){
    if (rebuild) cell_list.build(positions);

    cell_list.get_cmat(positions, conn_mat);
    noisy_align();

    for (int d = 0; d < 2; d++){ positions.row(d) += velocities.row(d); }
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
    Property noise_phi{1, this->n};
    Coord2D noise_xy{2, this->n};
    PropertyInt neighbour_nums{1, this->n};
    Coord2D new_velocities{2, this->n};
    new_velocities.setZero();
    Vec2D shift;
    Vec2D v_force;
    Vec2D v_avoid;
    double dij{0};
    double n_avoid{0};
    double a = alpha_ / speed;

    neighbour_nums = conn_mat.colwise().sum(); // -> m_i

    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= PI; // phi ~ U(-PI, PI)

    noise_xy.row(0) = cos(noise_phi) * noise; // -> eta * xi_i^t
    noise_xy.row(1) = sin(noise_phi) * noise;

    for (int i = 0; i < this->n; i++){
        v_force << 0, 0;
        v_avoid << 0, 0;
        n_avoid = 0;
        for (int j = 0; j < this->n; j++){  // \sum n_ij^t s_j^t
            if (conn_mat(i, j) > 0){
                if (i == j){
                    v_force += a * velocities.col(j);
                } else {
                    shift = get_shift(positions.col(j), positions.col(i));
                    dij = shift.matrix().norm();
                    shift /= dij;
                    if (dij < rc_) {
                        n_avoid += 1;  // fij = infty
                        v_avoid += shift * -1;
                    } else if (dij < ra_) {
                        v_force += a * velocities.col(j) + beta_ * shift * c_ * (dij - re_);
                    } else {
                        v_force += a * velocities.col(j) + beta_ * shift;  // fij = 1
                    }
                }
            }
        }
        if (n_avoid > 0){
            new_velocities.col(i) = v_avoid;
        } else {
            new_velocities.col(i) = v_force + noise_xy.col(i) * neighbour_nums(0, i);
        }
    }

    normalise(new_velocities, speed);
    velocities = new_velocities;
}


void Vicsek2DPBCVNCO::move(bool rebuild){
    if (rebuild) cell_list.build(positions);

    cell_list.get_cmat(positions, conn_mat);
    update_velocity();

    for (int d = 0; d < 2; d++){ positions.row(d) += velocities.row(d); }
    fix_positions();
}

