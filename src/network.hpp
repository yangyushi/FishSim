#ifndef NETWORK
#define NETWORK

#include <algorithm>
#include <numeric>
#include "core.hpp"

using Edges = Indices2D;  // vector < array< int, 2 > >
using Nodes = std::vector< int >;


struct Graph {
    Nodes nodes_;
    Edges edges_;
    Graph();
    Graph(int n);
    Graph(ConnMat adj_mat);
    Graph(Nodes nodes, Edges edges);
    Conn as_connections();
    ConnMat as_matrix();
    inline int size() { return nodes_.size(); }
};

Graph random_regular_graph(int d, int n);

// random graph where all nodes have d random neighbours
Graph random_vnm_graph(int d, int n);

// random graph where all nodes have d random neighbours, must
// include itself
Graph random_vnm_graph_force_self(int d, int n);

/*
 * The vectorial network model introduced by Aldana and Huepe in 2003
 */
class Network3D : public AVS3D{
    protected:
        int k_; // connection number 
        int dynamic_ = 0;  // 0 - quenched dynamic; 1 - annealed dynamic
        Graph graph_;
        Conn connections_;

    public:
        Network3D(int n, int k, double eta);
        void random_align();
        void move(bool new_graph);
        void evolve(int steps, int dynamic);
};

#endif
