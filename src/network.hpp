#ifndef NETWORK
#define NETWORK

#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include "core.hpp"


struct PairHasher{
    std::size_t operator()(const Index2D&) const;
};

using Nodes = std::vector< int >;
using Edges = std::unordered_set<Index2D, PairHasher>;


struct Graph {
    Nodes nodes_;
    Edges edges_;
    Graph();
    Graph(int n);
    Graph(ConnMat adj_mat);
    Graph(Nodes nodes, Edges edges, bool directional);
    Conn as_connections();
    ConnMat as_matrix();
    inline int size() { return nodes_.size(); }
};


// random simple regular graph
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
    
    private:
        virtual void update_graph();

    public:
        Network3D(int n, int k, double eta);
        virtual ~Network3D(){};
        void move(bool new_graph);
        void evolve(int steps, int dynamic);
        ConnMat get_adj_mat();
        void set_adj_mat(ConnMat adj_mat);
};


class Network3DRG : public Network3D {
    public:
        Network3DRG(int n, int k, double eta);

    private:
        void update_graph();
};


class Voter : public ABS {
    protected:
        int k_;
        Graph graph_;
        Conn connections_;

    private:
        void update_graph();

    public:
        Voter(int n, int k, double eta);
        void move(bool new_graph);
        void evolve(int steps, int dynamic);
        inline ConnMat get_adj_mat(){ return graph_.as_matrix(); };
        void set_adj_mat(ConnMat);
};

#endif
