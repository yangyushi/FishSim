#include "vicsek.h"

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


Vicsek3DPBC::Vicsek3DPBC(int n, double r, double eta, double box, double v0) :
        noise(eta), speed(v0), box(box), n(n), conn_mat(n, n),
        positions(3, n), velocities(3, n), cell_list(r, box, true){
        positions.setRandom(3, n);
        positions = (positions + 1) / 2 * box;
        velocities.setRandom(3, n);
        normalise(velocities, v0);
    }


void Vicsek3DPBC::align(){
    Coord3D new_velocities{3, this->n};
    new_velocities.setZero();

    for (int i = 0; i < this->n; i++){

        for (int j = 0; j < conn_mat.rows(); j++){
            if (conn_mat(i, j) > 0){
                new_velocities.col(i) += velocities.col(j);
            }
        }

    }
    velocities << new_velocities;
}


void Vicsek3DPBC::rotate_noise(Coord3D& noise_xyz){
    /*
    * Rotate nosie pointing at (0, 0, 1) to direction xyz 
    * and then add noise_xyz to velocities
    * here velocities have unit norm
    */
    RotMat F, G, R;
    Vec3D A, B, u, v, w;
    double x{0}, y{0}, z{0}, rxy{0};
    for (int i = 0; i < this->n; i++){
        x = velocities(0, i);
        y = velocities(1, i);
        z = velocities(2, i);
        rxy = sqrt(x * x + y * y);
        A << 0, 0, 1;
        B << x, y, z;
        if ((A - B).matrix().norm() < 1e-10){
            velocities.col(i) = noise_xyz.col(i);
        } else {
            u = A;
            v << x / rxy, y / rxy, 0;
            w = B.matrix().cross(A.matrix());  // w = B x A
            G <<   z, -rxy, 0,
                 rxy,    z, 0,
                   0,    0, 1;
            F << u.transpose(), v.transpose(), w.transpose();
            R = F.inverse() * (G * F);
            velocities.col(i) = (R * noise_xyz.col(i).matrix()) * speed;
        }
    }
}


void Vicsek3DPBC::add_noise(){
    Property noise_rxy{1, this->n}, noise_phi{1, this->n};
    Coord3D noise_xyz{3, this->n};

    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= PI; // phi ~ U(-PI, PI)

    noise_xyz.row(2).setRandom(); // ~ U(-1, 1)
    noise_xyz.row(2) = noise_xyz.row(2) * noise - noise + 1; // z ~ U(1 - 2 * noise, 1)
    noise_rxy = sqrt(1 - noise_xyz.row(2).pow(2));
    noise_xyz.row(0) = noise_rxy * cos(noise_phi);
    noise_xyz.row(1) = noise_rxy * sin(noise_phi);
    rotate_noise(noise_xyz);
}


void Vicsek3DPBC::move(bool rebuild){
    if (rebuild) cell_list.build(positions);

    cell_list.get_cmat(positions, conn_mat);
    align();
    normalise(velocities);
    add_noise();

    for (int d = 0; d < 3; d++){
        positions.row(d) += velocities.row(d);
    }

    fix_positions();
}


void Vicsek3DPBC::dump(string filename){
    ofstream f;
    f.open(filename, ios::out | ios::app);
    f << this->n << endl;
    f << "id, x, y, z, vx, vy, vz" << endl;
    for (int i = 0; i < this->n; i++ ) {
        f << i << " "
          << positions(0, i)  << " " <<  positions(1, i)  << " " << positions(2, i)  << " "
          << velocities(0, i) << " " <<  velocities(1, i) << " " << velocities(2, i) << endl;
    }
    f.close();
}


void Vicsek3DPBC::load(string filename){
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
                ss >> positions(j, i - head_lines);
            }
            for (int j = 0; j < 3; j++){
                ss >> velocities(j, i - head_lines);
            }
        }
    }
    f.close();
}


Vicsek3DPBCVN::Vicsek3DPBCVN(int n, double r, double eta, double box, double v0)
    : Vicsek3DPBC(n, r, eta, box, v0) {}


void Vicsek3DPBCVN::noisy_align(){
    /*
     * See ginellEPJ2016 for equations
     */
    Property noise_phi{1, this->n};
    Property noise_theta{1, this->n};
    PropertyInt neighbour_nums{1, this->n};

    Coord3D noise_xyz{3, this->n};
    Coord3D new_velocities{3, this->n};
    Vec3D v_align;

    neighbour_nums = conn_mat.colwise().sum(); // -> m_i

    noise_phi.setRandom();   noise_phi   = (noise_phi + 1) * PI; // phi ~ U(0, 2PI)
    noise_theta.setRandom(); noise_theta = asin(noise_theta);    // theta ~ U(0, PI)

    noise_xyz.row(0) = cos(noise_theta) * cos(noise_phi) * speed * noise; // -> eta * xi_i^t
    noise_xyz.row(1) = cos(noise_theta) * sin(noise_phi) * speed * noise;
    noise_xyz.row(2) = sin(noise_theta) * speed * noise;

    for (int i = 0; i < this->n; i++){
        v_align << 0, 0, 0;

        for (int j = 0; j < conn_mat.cols(); j++){  // \sum n_ij^t s_j^t
            if (conn_mat(i, j) > 0){
                v_align += velocities.col(j);
            }
        }
        new_velocities.col(i) = v_align + noise_xyz.col(i) * neighbour_nums(0, i);
    }

    normalise(new_velocities, speed);
    velocities = new_velocities;
}


void Vicsek3DPBCVN::move(bool rebuild){
    if (rebuild) cell_list.build(positions);

    cell_list.get_cmat(positions, conn_mat);
    noisy_align();

    positions += velocities;
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

