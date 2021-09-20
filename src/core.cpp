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

RotMat get_rotation_matrix(Vec3D v1, Vec3D v2){
    RotMat F, G, R;
    R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if ((v1 - v2).norm() <= 100 * std::numeric_limits<double>::epsilon()){
        return R;
    }
    double n1 = (double) v1.norm(),
           n2 = (double) v2.norm();
    double tmp = v1.transpose() * v2;
    double t12 = acos(tmp / (n1 * n2));
    Vec3D u, v, w;
    if (cos(t12) > 0){
        u << v1 / n1 * n2 * cos(t12);
        v << v2 - u;
        u = u / u.norm();
        v = v / v.norm();
    } else {
        u << v1 / n1 * n2 * -cos(t12);
        v << v2 + u;
        u = u / u.norm();
        v = v / v.norm();
    }
    w << u.cross(v);
    G << cos(t12), -sin(t12), 0,
         sin(t12), cos(t12), 0,
         0, 0, 1;
    F << u.transpose(), v.transpose(), w.transpose();
    R = F.inverse() * G * F;

    return R;
}

RotMat get_rotation_matrix(Vec3D v1, Vec3D v2, double theta){
    RotMat F, G, R;
    R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if (abs((v1 - v2).norm()) <= 100 * std::numeric_limits<double>::epsilon()){
        return R;
    }
    double n1 = v1.norm(),
           n2 = v2.norm();
    double tmp = v1.transpose() * v2;
    double t12 = acos(tmp / (n1 * n2));
    Vec3D u, v, w;
    if (cos(t12) > 0){
        u << v1 / n1 * n2 * cos(t12);
        v << v2 - u;
        u = u / u.norm();
        v = v / v.norm();
    } else {
        u << v1 / n1 * n2 * -cos(t12);
        v << v2 + u;
        u = u / u.norm();
        v = v / v.norm();
    }
    w = u.cross(v);
    if (abs(t12) > theta){
        if (t12 < 0){
            t12 = -theta;
        } else {
            t12 = theta;
        }
    }
    G << cos(t12), -sin(t12), 0,
         sin(t12), cos(t12), 0,
         0, 0, 1;
    F << u.transpose(), v.transpose(), w.transpose();
    R = F.inverse() * (G * F);
    return R;
}


void voter_align(Spins& spins, Property& states, const Conn& connections){
    size_t n = spins.cols();
    states.setZero();
    for (size_t i = 0; i < n; i++){
        for (auto j : connections[i]){
            states(0, i) += spins(0, j);
        }
    }
    for (size_t i = 0; i < n; i++){
        if (states(0, i) <= 0){
            spins(0, i) = -1;
        } else {
            spins(0, i) = 1;
        }
    }
}


AVS2D::AVS2D(int n, double noise, double speed)
    : n_(n), speed_(speed), velocities_(dim_, n), noise_(noise) {
    // generate random vector on a unit circle

    if (noise > 1){
        throw std::invalid_argument("noise should in range (0, 1)");
    } else if (noise < 0){
        throw std::invalid_argument("noise should in range (0, 1)");
    }
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


AVS3D::AVS3D(int n, double noise)
    : n_(n), orientations_(dim_, n), noise_(noise){

    if (noise > 1){
        throw std::invalid_argument("noise should in range (0, 1)");
    } else if (noise < 0){
        throw std::invalid_argument("noise should in range (0, 1)");
    }
    Property vz{n}, vphi{n}, vrxy{n};
    vz.setRandom(); // vz ~ U(-1, 1)
    vrxy = sqrt(1 - vz.pow(2));
    vphi.setRandom();
    vphi *= M_PI;  // vphi ~ U(-pi, pi)
    orientations_.row(0) << vrxy * cos(vphi);
    orientations_.row(1) << vrxy * sin(vphi);
    orientations_.row(2) << vz;
    normalise(orientations_);
}

void AVS3D::load_orientations(Coord3D velocities) {
    orientations_ = velocities;
    normalise(orientations_);
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
            orientations_(0, i), orientations_(1, i), orientations_(2, i)
        );
        orientations_.col(i) = R * noise_xyz.col(i);
    }
}


ABS::ABS(int n, double noise)
    : n_(n), spins_{1, n}, noise_(noise), states_{1, n}
{
    if (noise > 1){
        throw std::invalid_argument("noise should in range (0, 1)");
    } else if (noise < 0){
        throw std::invalid_argument("noise should in range (0, 1)");
    }

    states_.setRandom();  // ~U(-1, 1)
    size_t i = 0;
    for (const auto& num : states_) {
        if (num > 0) {
            spins_(0, i) = 1;
        } else {
            spins_(0, i) = -1;
        }
        i++;
    }
}


void ABS::add_noise(){
    states_.setRandom();  // ~U(0, 1)
    size_t i = 0;
    double threshold = noise_ - 1;
    for (const auto& num : states_) {
        if (num < threshold) {
            spins_(0, i) = -spins_(0, i);
        }
        i++;
    }
}
