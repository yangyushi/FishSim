#ifndef NETWORK
#define NETWORK

#include <algorithm>
#include "core.hpp"

using Edges = Indices2D;  // vector < array< int, 2 > >
using Nodes = std::vector< int >;


struct Graph {
    Edges edges_;
    Nodes nodes_;
    Graph(int n);
    Conn as_connections();
    ConnMat as_matrix();
    inline int size() {return *std::max_element(nodes_.begin(), nodes_.end()) + 1;}
};

Graph random_regular_graph(int d, int n);


/*
 * The vectorial network model introduced by Aldana and Huepe in 2003
 */
class Network3D : public AVS3D{
    protected:
        int k_; // connection number 
        std::vector<int> indices_;
        int dynamic_ = 0;  // 0 - quenched dynamic; 1 - annealed dynamic
        Graph graph_;

    public:
        Network3D(int n, int k, double eta);
        void random_align();
        void move();
        void evolve(int steps);
};

#endif
