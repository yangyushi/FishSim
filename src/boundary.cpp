#include "boundary.hpp"

Tank3D::Tank3D(double c, double z_max, double kw) :
    c_(c), z_max_(z_max), kw_(kw),
    r_max_(sqrt(z_max / c)), volume_(M_PI/2/c * z_max*z_max) {
}


Coord3D Tank3D::get_random_positions(size_t n){
    Property rand_z2{n}, rand_r2{n}, rand_r{n}, rand_phi{n}, rand_z{n};
    rand_z2.setRandom();
    rand_phi.setRandom();

    rand_z2 = (rand_z2 + 1) / 2;  // ~ U(0, 1)
    rand_z = sqrt(rand_z2) * z_max_;

    rand_r2.setRandom();  // ~U(-1, 1)
    rand_r2 = (rand_r2 + 1) / 2;  // ~U(0, 1)

    rand_r = sqrt(rand_r2 * rand_z / c_);
    rand_phi *= M_PI; // ~U(-PI, PI)

    Coord3D rand_xyz{3, n};

    rand_xyz.row(0) << rand_r * cos(rand_phi);
    rand_xyz.row(1) << rand_r * sin(rand_phi);
    rand_xyz.row(2) << rand_z;

    return rand_xyz;
}

bool Tank3D::is_inside(double x, double y, double z){
    if (z > z_max_) {
        return false;
    }
    double r2 = x * x + y * y;
    if (c_ * r2 < z) {
        return true;
    }
    return false;
}

bool Tank3D::is_inside(const Vec3D position){
    double x = position(0, 0);
    double y = position(1, 0);
    double z = position(2, 0);
    double r2 = x * x + y * y;

    if (z > z_max_) {
        return false;
    }

    if (c_ * r2 < z) {
        return true;
    }
    return false;
}


Coord3D Tank3D::get_random_positions_reject(size_t n){
    Coord3D result{3, n};
    std::uniform_real_distribution<> dist_r(-r_max_, r_max_);
    std::uniform_real_distribution<> dist_z(0, z_max_);
    std::default_random_engine generator;
    size_t count = 0;
    double rand_x, rand_y, rand_z;
    while (count < n){
        rand_x = dist_r(generator);
        rand_y = dist_r(generator);
        rand_z = dist_z(generator);
        if (this->is_inside(rand_x, rand_y, rand_z)){
            result.col(count) << rand_x, rand_y, rand_z;
            count += 1;
        }
    }
    return result;
}


Coord3D Tank3D::project(const Coord3D& xyz){
    using namespace std::complex_literals;

    Property rsq = pow(xyz.row(0).array(), 2) +
                  pow(xyz.row(1).array(), 2);

    Property phi = xyz.row(1).binaryExpr(
        xyz.row(0), [](double y, double x) {return atan2(y, x);}
    );

    // different constants for the roots
    double k0 = pow(2, 1.0 / 3.0);
    double k3 = 12 * k0 * c_ * c_;
    double r, r0, r1, r2, z, z0, z1, z2, dsq0, dsq1, dsq2, r_proj, z_proj;
    std::complex<double> k1 = 1.0 + sqrt(3) * 1i;
    std::complex<double> k2 = 1.0 - sqrt(3) * 1i;

    PropertyComplex _t1 = 2 * c_ * xyz.row(2).array() - 1; // 2cz - 1
    PropertyComplex _t2 = 864 * pow(c_, 6) * pow(_t1, 3);  // 864 c^6 (2cz - 1)^3
    PropertyComplex _t3 = 11664 * pow(c_, 8) * rsq;         // 11664 c^8 r^2
    PropertyComplex _t4 = 108 * pow(c_, 4) * sqrt(rsq);     // 108 c^4 r^2
    PropertyComplex _t5 = pow(-_t4 + sqrt(_t3 - _t2), 1.0/3.0);  // (-108 xxx)^(1/3)
    PropertyComplex _t5_inv = 1.0 / _t5;

    Property r_proj_0 = (- k0 * _t1 * _t5_inv - _t5 / (6 * k0 * c_ * c_)).real();
    Property r_proj_1 = ((k1 / k0 / k0) * _t1 * _t5_inv +  k2 / k3 * _t5).real();
    Property r_proj_2 = ((k2 / k0 / k0) * _t1 * _t5_inv +  k1 / k3 * _t5).real();

    Coord3D xyz_proj {3, xyz.cols()};
    for (size_t idx = 0; idx < xyz.cols(); idx++){
        r = sqrt(rsq(idx)); z = xyz(2, idx);
        r0 = r_proj_0(idx); z0 = c_ * r0 * r0;
        r1 = r_proj_1(idx); z1 = c_ * r1 * r1;
        r2 = r_proj_2(idx); z2 = c_ * r2 * r2;
        dsq0 = pow(r - r0, 2) + pow(z - z0, 2);
        dsq1 = pow(r - r1, 2) + pow(z - z1, 2);
        dsq2 = pow(r - r2, 2) + pow(z - z2, 2);
        std::vector<double> dists_sq {dsq0, dsq1, dsq2};
        std::vector<double> r_vals {r0, r1, r2};
        size_t idx_min = std::min_element(dists_sq.begin(), dists_sq.end()) - dists_sq.begin();
        r_proj = r_vals[idx_min];
        z_proj = c_ * r_proj * r_proj;
        xyz_proj.col(idx) << r_proj * cos(phi(idx)), r_proj * sin(phi(idx)), z_proj;
    }

    return xyz_proj;
}


Coord3D Tank3D::project_single(const Vec3D& xyz){
    using namespace std::complex_literals;

    double x, y, z;
    x = xyz(0, 0); y = xyz(1, 0); z = xyz(2, 0);
    double rsq = pow(x, 2) + pow(y, 2);
    double phi = atan2(y, x);

    // different constants for the roots
    double k0 = pow(2, 1.0 / 3.0);
    double k3 = 12 * k0 * c_ * c_;
    double r, r0, r1, r2, z0, z1, z2, dsq0, dsq1, dsq2, r_proj, z_proj;
    std::complex<double> k1 = 1.0 + sqrt(3) * 1i;
    std::complex<double> k2 = 1.0 - sqrt(3) * 1i;

    double _t1 = 2 * c_ * z - 1; // 2cz - 1
    double _t2 = 864 * pow(c_, 6) * pow(_t1, 3);  // 864 c^6 (2cz - 1)^3
    std::complex<double> _t3 = 11664 * pow(c_, 8) * rsq;         // 11664 c^8 r^2
    std::complex<double> _t4 = 108 * pow(c_, 4) * sqrt(rsq);     // 108 c^4 r^2
    std::complex<double> _t5 = pow(-_t4 + sqrt(_t3 - _t2), 1.0/3.0);  // (-108 xxx)^(1/3)
    std::complex<double> _t5_inv = 1.0 / _t5;

    r0 = (- k0 * _t1 * _t5_inv - _t5 / (6 * k0 * c_ * c_)).real();
    r1 = ((k1 / k0 / k0) * _t1 * _t5_inv +  k2 / k3 * _t5).real();
    r2 = ((k2 / k0 / k0) * _t1 * _t5_inv +  k1 / k3 * _t5).real();

    z0 = c_ * r0 * r0;
    z1 = c_ * r1 * r1;
    z2 = c_ * r2 * r2;

    Vec3D xyz_proj {3, xyz.cols()};
    r = sqrt(rsq);

    dsq0 = pow(r - r0, 2) + pow(z - z0, 2);
    dsq1 = pow(r - r1, 2) + pow(z - z1, 2);
    dsq2 = pow(r - r2, 2) + pow(z - z2, 2);
    std::vector<double> dists_sq {dsq0, dsq1, dsq2};
    std::vector<double> r_vals {r0, r1, r2};
    size_t idx_min = std::min_element(dists_sq.begin(), dists_sq.end()) - dists_sq.begin();
    r_proj = r_vals[idx_min];
    z_proj = c_ * r_proj * r_proj;

    xyz_proj << r_proj * cos(phi), r_proj * sin(phi), z_proj;

    return xyz_proj;

}

Vec3D Tank3D::to_rzt(const Vec3D& xyz){
    Vec3D rzt {3, 1};
    rzt.row(0) << sqrt(xyz(0, 0) * xyz(0, 0) + xyz(1, 0) * xyz(1, 0));
    rzt.row(1) << xyz(2, 0);
    rzt.row(2) << atan2(xyz(1, 0), xyz(0, 0));
    return rzt;
}


void Tank3D::fix_orientations(const Coord3D& positions, Coord3D& orientations, double dt){
    size_t n = positions.cols();
    Vec3D pos, orient, o_proj, o_target_cap, o_target_base, o_target;
    double dist_to_base, dist_to_cap;

    for (size_t i = 0; i < n; i++){
        pos << positions.col(i);
        orient << orientations.col(i);
        o_proj << pos - this->project_single(pos);
        dist_to_base = o_proj.norm();
        dist_to_cap = abs(z_max_ - pos(2, 0));

        o_target_base = get_rotation_matrix(
            orient, o_proj, kw_ * pow(dist_to_base, -2) * dt
        ) * orient;

        o_target_cap = get_rotation_matrix(
            orient, Vec3D{0, 0, -1}, kw_ * pow(dist_to_cap, -2) * dt
        ) * orient;

        o_target = o_target_cap + o_target_base;
        normalise(o_target);
        orientations.col(i) = o_target;
    }
}


void Tank3D::fix_positions(Coord3D& positions){
    size_t n = positions.cols();
    Vec3D pos, pos_proj;
    double r2, theta;
    double r2_max = r_max_ * r_max_;

    for (size_t i = 0; i < n; i++){
        pos << positions.col(i);
        r2 = pos(0, 0) * pos(0, 0) + pos(1, 0) * pos(1, 0);
        if (not this->is_inside(pos)){
            if (pos(2, 0) > z_max_){
                pos(2, 0) = z_max_;
                if (r2 > r2_max){
                    theta = atan2(pos(1, 0), pos(0, 0));
                    pos(0, 0) = r_max_ * cos(theta);
                    pos(1, 0) = r_max_ * sin(theta);
                }
            } else {
                pos << this->project_single(pos);
            }
            positions.col(i) = pos;
        }
    }
}
