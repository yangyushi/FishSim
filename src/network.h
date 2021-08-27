#ifndef NETWORK
#define NETWORK

#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <Eigen/Dense>
#include "vicsek.h"
#include "neighbour_list.h"

class Network3D : public Vicsek3D{
    protected:
        int k_; // connection number 
        vector<int> indices_;

    public:
        Network3D(int n, int k, double eta);
        void random_align();
        void move();
        void evolve(int steps);
};


#endif
