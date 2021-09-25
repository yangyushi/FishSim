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
    /*
     * it is not numerically stable
     * good region: c ~ 0.7; z_max ~ 0.5
     * TODO: implement a numerical root finding projection
     */
    double p = pow(2, 1.0 / 3.0);

    Property r2 = pow(xyz.row(0).array(), 2) +
                  pow(xyz.row(1).array(), 2);

    Property phi = xyz.row(1).binaryExpr(
        xyz.row(0), [](double y, double x) {return atan2(y, x);}
    );

    Property _t1 = 2 * c_ * xyz.row(2).array() - 1;
    Property _t2 = 55296 * pow(c_, 6) * pow(_t1, 3);
    Property _t3 = 746496 * pow(c_, 8) * r2;
    Property _t4 = 864 * pow(c_, 4) * sqrt(r2);
    Property _t5 = pow(_t4 + sqrt(_t3 - _t2), 1.0/3.0);

    Property _term_1 = (2 * p * _t1) / _t5;
    Property _term_2 = _t5 / (12 * p * c_*c_);

    Property r_proj = _term_1 + _term_2;
    Property z_proj = c_ * pow(r_proj, 2);

    Coord3D xyz_proj {3, xyz.cols()};
    xyz_proj.row(0) << r_proj * cos(phi);
    xyz_proj.row(1) << r_proj * sin(phi);
    xyz_proj.row(2) << z_proj;

    return xyz_proj;
}


Coord3D Tank3D::project_single(const Vec3D& xyz){
    /*
     * it is not numerically stable
     * good region: c ~ 0.7; z_max ~ 0.5
     * TODO: implement a numerical root finding projection
     */
    double p = pow(2, 1.0 / 3.0);

    Property r2 = pow(xyz.row(0).array(), 2) +
                  pow(xyz.row(1).array(), 2);

    Property phi = xyz.row(1).binaryExpr(
        xyz.row(0), [](double y, double x) {return atan2(y, x);}
    );

    Property _t1 = 2 * c_ * xyz.row(2).array() - 1;
    Property _t2 = 55296 * pow(c_, 6) * pow(_t1, 3);
    Property _t3 = 746496 * pow(c_, 8) * r2;
    Property _t4 = 864 * pow(c_, 4) * sqrt(r2);
    Property _t5 = pow(_t4 + sqrt(_t3 - _t2), 1.0/3.0);

    Property _term_1 = (2 * p * _t1) / _t5;
    Property _term_2 = _t5 / (12 * p * c_*c_);

    Property r_proj = _term_1 + _term_2;
    Property z_proj = c_ * pow(r_proj, 2);

    Vec3D xyz_proj {3, 1};
    xyz_proj.row(0) << r_proj * cos(phi);
    xyz_proj.row(1) << r_proj * sin(phi);
    xyz_proj.row(2) << z_proj;

    return xyz_proj;
}

Vec3D Tank3D::to_rzt(const Vec3D& xyz){
    Vec3D rzt {3, 1};
    rzt.row(0) << sqrt(xyz(0, 0) * xyz(0, 0) + xyz(1, 0) * xyz(1, 0));
    rzt.row(1) << xyz(2, 0);
    rzt.row(2) << atan2(xyz(1, 0), xyz(0, 0));
    return rzt;
}


void Tank3D::fix_orientations(const Coord3D& positions, Coord3D& orientations){
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
            orient, o_proj, kw_ * pow(dist_to_base, -2)
        ) * orient;

        o_target_cap = get_rotation_matrix(
                orient, Vec3D{0, 0, -1}, kw_ * pow(dist_to_cap, -2)
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
