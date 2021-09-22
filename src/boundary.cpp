#include "boundary.hpp"

FishBowl3D::FishBowl3D(double c, double z_max) :
    c_(c), z_max_(z_max),
    r_max_(sqrt(z_max / c)), volume_(M_PI/2/c * z_max*z_max) {
}


Coord3D FishBowl3D::get_random_positions(size_t n){
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

bool FishBowl3D::is_inside(double x, double y, double z){
    if (z > z_max_) {
        return false;
    }
    double r2 = x * x + y * y;
    if (c_ * r2 < z) {
        return true;
    }
    return false;
}

bool FishBowl3D::is_inside(const Vec3D position){
    double x = position(0, 0);
    double y = position(1, 0);
    double z = position(2, 0);

    if (z > z_max_) {
        return false;
    }
    double r2 = x * x + y * y;
    if (c_ * r2 < z) {
        return true;
    }
    return false;
}


Coord3D FishBowl3D::get_random_positions_reject(size_t n){
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


Coord3D FishBowl3D::project(const Coord3D& xyz){
    double p = pow(2, 1.0 / 3.0);
    double sx = 1 / xyz.row(0).maxCoeff();
    double sz = 1 / xyz.row(2).maxCoeff();
    double c_rescale = c_ / sz * sx * sx;
    sx = 1;
    sz = 1;

    std::cout << "scale: "  << sx << ", " << sz
        << "; c_rescale: " << c_rescale << std::endl;

    Coord3D xyz_rescale {3, xyz.cols()};
    xyz_rescale << xyz;
    xyz_rescale.row(0) *= sx;
    xyz_rescale.row(1) *= sx;
    xyz_rescale.row(2) *= sz;

    Property r2 = pow(xyz_rescale.row(0).array(), 2) +
                  pow(xyz_rescale.row(1).array(), 2);

    Property phi = xyz.row(1).binaryExpr(
        xyz.row(0), [](double y, double x) {return atan2(y, x);}
    );

    Property term = 108 * pow(c_rescale, 4) * sqrt(r2) + sqrt(
        11664 * pow(c_rescale, 8) * r2 +
        864 * pow(c_rescale, 6) * pow(
                (1 - 2 * c_rescale * xyz_rescale.row(2).array()), 3
            )
    );
    term = pow(term, 1.0 / 3.0);
    size_t count = 0;
    for (auto num : term){
        std::cout << count << ": " << num
            << "("
            << xyz_rescale(0, count) << ", "
            << xyz_rescale(1, count) << ", "
            << xyz_rescale(2, count) << ", "
            << "), r2 = "
            << r2(0, count) << ", "
            << "; phi = "
            << phi(0, count) << ", "
            << std::endl;
        count ++;
    }
    Property r_proj = -(p * (
                1 - 2 * c_rescale * xyz_rescale.row(2).array()
                )) / term +
        term / (6 * p * c_rescale * c_rescale);
    Property z_proj = c_rescale * pow(r_proj, 2);

    Coord3D xyz_proj {3, xyz.cols()};
    xyz_proj.row(0) << r_proj * cos(phi) / sx;
    xyz_proj.row(1) << r_proj * sin(phi) / sx;
    xyz_proj.row(2) << z_proj / sz;

    return xyz_proj;
}

void FishBowl3D::fix_orientations(const Coord3D& positions, Coord3D& orientations){
    size_t n = positions.cols();
    Vec3D pos;
    for (size_t i = 0; i < n; i++){
        pos << positions.col(i);
        if (not this->is_inside(pos)){
        }
    }
}
