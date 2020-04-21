#include "vicsek.h"

const double PI  =3.141592653589793238463;


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
        xyz.col(i) /= L;
        i++;
    }
    xyz *= spd;
}


Vicsek3DPBC::Vicsek3DPBC(int n, double r, double eta, double box, double v0) :
        noise(eta), speed(v0), box(box), n(n),
        positions(3, n), velocities(3, n), cell_list(r, box, true) {
        positions.setRandom(3, n);
        positions = (positions + 1) / 2 * box;
        velocities.setRandom(3, n);
        normalise(velocities, v0);
    }


void Vicsek3DPBC::align(){
    int neighbour_num = 0;
    Coord3D new_velocities{3, n};
    Vec3D v_align;

    for (int i = 0; i < this->n; i++){
        v_align << 0, 0, 0;
        neighbour_num = 0;

        for (int j = 0; j < dist_mat.rows(); j++){
            if (dist_mat(i, j) >= 0){
                v_align += velocities.col(j);
                neighbour_num += 1;
            }
        }
        new_velocities.col(i) << v_align / neighbour_num;
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
        rxy = sqrt(1 - z * z);
        A << 0, 0, 1;
        B << x, y, z;
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
    this->dist_mat = cell_list.get(positions);

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
