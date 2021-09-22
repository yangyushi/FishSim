#ifndef BOUNDARY
#define BOUNDARY
#include "core.hpp"


struct Boundary{
    virtual void fix_positions(Coord3D& positions){};
    virtual void fix_orientations(
        const Coord3D& positions, Coord3D& orientations
    ){};
    Boundary(){};
    virtual ~Boundary(){};
};


/*
 * z = c * r^2
 * r = sqrt(z / c)
 */
struct FishBowl3D : public Boundary {
    // generate uniform random points inside the boundary 
    Coord3D get_random_positions(size_t n);
    Coord3D get_random_positions_reject(size_t n);
    // for all 3D points, calculate the cloest point on the bottom of the bowl
    Coord3D project(const Coord3D& xyz);
    const double c_, z_max_, r_max_, volume_;

    FishBowl3D(double c, double z_max);

    bool is_inside(double x, double y, double z) ;
    bool is_inside(Vec3D position) ;
    inline void fix_positions(Coord3D& positions) {} ; // do nothing
    void fix_orientations(
        const Coord3D& positions, Coord3D& orientations
    );
};

#endif
