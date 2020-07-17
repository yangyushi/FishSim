#ifndef SPIN
#define SPIN

#include <vector>
#include <array>
#include <map>
#include <Eigen/Dense>
#include "neighbour_list.h"

/*
 *  * please refer to 10.1007/s10955-014-1119-3 for the origin of the model
 *   */
class InertialSpin {
    private:
        double friction_;
        double dt_;
};

#endif
