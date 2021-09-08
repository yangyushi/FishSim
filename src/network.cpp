#include "network.hpp"

std::random_device rd;
std::mt19937 g(rd());

Graph::Graph() : nodes_{}, edges_{} {;}
Graph::Graph(Nodes n, Edges e) : nodes_{n}, edges_{e} {;}
Graph::Graph(int n) : nodes_{n}{
    std::iota(nodes_.begin(), nodes_.end(), 0);
}

Conn Graph::as_connections(){
    Conn connections {};
    for (int i = 0; i < this->size(); i++){
        connections.push_back(Nodes{});
    }
    for (auto edge : edges_){
        connections[edge[0]].push_back(edge[1]);
    }
    return connections;
}

ConnMat Graph::as_matrix(){
    int n = this->size();
    ConnMat mat{n, n};
    mat.setZero();
    for (auto edge : edges_){
        mat(edge[0], edge[1]) += 1;
    }
    return mat;
}


Graph random_regular_graph(int d, int n){
    if (n * d % 2 == 0){
        throw("In valid argument, n * d shuold be even");
    }
    if ((d < 0) or (d >= n)) {
        throw("the 0 <= d < n inequality must be satisfied");
    }
    if (d == 0) {return Graph(n);}
    return Graph(n);
}


Graph random_vnc_graph(int d, int n){
    Edges edges{};
    Nodes nodes(n);
    std::iota(nodes.begin(), nodes.end(), 0);

    for (int i = 0; i < n; i++){
        shuffle(nodes.begin(), nodes.end(), g);
        for (int j = 0; j < d; j++){
            edges.push_back(std::array<int, 2> {i, nodes[j]});
        }
    }
    return Graph{nodes, edges};
}

Graph random_vnc_graph_force_self_loop(int d, int n){
    Edges edges ;
    Nodes nodes(n);
    std::iota(nodes.begin(), nodes.end(), 0);

    for (int i1 = 0; i1 < n; i1++){
        bool self_included = false;
        shuffle(nodes.begin(), nodes.end(), g);
        for (int j = 0; j < d - 1; j++){
            int i2 = nodes[j];
            if (i2 == i1){ self_included = true; }
            edges.push_back(std::array<int, 2> {i1, i2});
        }
        if (self_included){
            edges.push_back(std::array<int, 2> {i1, nodes[d - 1]});
        } else {
            edges.push_back(std::array<int, 2> {i1, i1});
        }
    }
    return Graph{nodes, edges};

}


Network3D::Network3D(int n, int k, double eta)
    : AVS3D{n, eta, 1.0}, k_(k) {
        graph_ = random_vnc_graph(k, n);
}


void Network3D::move(bool new_graph){
    if (new_graph){
        graph_ = random_vnc_graph(k_, n_);
    }
    vicsek_align(velocities_, graph_.as_connections());
    normalise(velocities_);
    add_noise();
}


void Network3D::evolve(int steps, int dynamic){
    if (dynamic == 0){
        for (int i = 0; i < steps; i++){
            this->move(false);
        }
    } else if (dynamic == 1){
        for (int i = 0; i < steps; i++){
            this->move(true);
        }
    } else {
        throw(
            "invalid dynamic type, choose between [0]: quenched or [1]: anneled"
        );
    }
}

