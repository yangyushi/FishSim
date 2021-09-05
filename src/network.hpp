#ifndef NETWORK
#define NETWORK

#include <vector>
#include <Eigen/Dense>
#include "core.hpp"

/*
 * The vectorial network model introduced by Aldana and Huepe in 2003
 */
class Network3D : public AVS3D{
    protected:
        int k_; // connection number 
        std::vector<int> indices_;
        int dynamic_ = 0;  // 0 - quenched dynamic; 1 - annealed dynamic

    public:
        Network3D(int n, int k, double eta);
        void random_align();
        void move();
        void evolve(int steps);

};

class TCNetwork3D : public Network3D{
};

#endif
