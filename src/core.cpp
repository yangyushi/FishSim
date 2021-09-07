#include "core.hpp"

std::default_random_engine generator;


Index3D unravel_index(int index, Index3D shape){
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


RotMat get_rotation_matrix(double x, double y, double z){
    double rxy = sqrt(x * x + y * y);
    RotMat F, G, R;
    Vec3D A, B, u, v, w;
    A << 0, 0, 1;
    B << x, y, z;
    u = A;
    v << x / rxy, y / rxy, 0;
    w = B.cross(A);  // w = B x A
    G <<   z, -rxy, 0,
         rxy,    z, 0,
           0,    0, 1;
    F << u.transpose(), v.transpose(), w.transpose();
    R = F.inverse() * (G * F);
    return R;
}



AVS2D::AVS2D(int n, double noise, double speed)
    : n_(n), speed_(speed), velocities_(dim_, n), noise_(noise) {
    // generate random vector on a unit circle
    Property vr{n_}, vphi{n_};
    vr.setConstant(speed_);
    vphi.setRandom();
    vphi *= M_PI;  // vphi ~ U(-pi, pi)
    velocities_.row(0) << vr * cos(vphi);
    velocities_.row(1) << vr * sin(vphi);
}

void AVS2D::load_velocities(Coord2D velocities) {
    this->velocities_ = velocities;
    normalise(this->velocities_, speed_);
};

void AVS2D::add_noise(){
    Property noise_phi{1, n_};
    Property phi = xy_to_phi(velocities_);
    noise_phi.setRandom(); // ~ U(-1, 1)
    noise_phi *= M_PI * noise_; // phi ~ U(-PI * eta, PI * eta)
    phi += noise_phi;
    velocities_.row(0) = speed_ * cos(phi);
    velocities_.row(1) = speed_ * sin(phi);
}


AVS3D::AVS3D(int n, double noise, double speed)
    : n_(n), speed_(speed), velocities_(dim_, n), noise_(noise){
    Property vz{n}, vphi{n}, vrxy{n};
    vz.setRandom(); // vz ~ U(-1, 1)
    vrxy = sqrt(1 - vz.pow(2));
    vphi.setRandom();
    vphi *= M_PI;  // vphi ~ U(-pi, pi)
    velocities_.row(0) << vrxy * cos(vphi);
    velocities_.row(1) << vrxy * sin(vphi);
    velocities_.row(2) << vz;
    normalise(velocities_, speed_);
}

void AVS3D::load_velocities(Coord3D velocities) {
    this->velocities_ = velocities;
    normalise(this->velocities_, speed_);
};

void AVS3D::add_noise(){
    Property noise_rxy{1, n_}, noise_phi{1, n_};
    Coord3D noise_xyz{3, n_};

    std::uniform_real_distribution<> dist_phi(-M_PI, M_PI);
    std::uniform_real_distribution<> dist_z(1 - 2 * noise_, 1);

    auto rand_phi = [&] (double) {return dist_phi(generator);};
    auto rand_z = [&] (double) {return dist_z(generator);};

    noise_phi.row(0) = Eigen::ArrayXd::NullaryExpr(n_, rand_phi); // phi ~ U(-PI, PI)
    noise_xyz.row(2) = Eigen::ArrayXd::NullaryExpr(n_, rand_z); //  z ~ U(1 - 2 * noise, 1)

    noise_rxy = sqrt(1 - noise_xyz.array().row(2).pow(2));
    noise_xyz.row(0) = noise_rxy * cos(noise_phi);
    noise_xyz.row(1) = noise_rxy * sin(noise_phi);

    for (int i = 0; i < n_; i++){
        RotMat R  = get_rotation_matrix(
            velocities_(0, i), velocities_(1, i), velocities_(2, i)
        );
        velocities_.col(i) = (R * noise_xyz.col(i)) * speed_;
    }
}
