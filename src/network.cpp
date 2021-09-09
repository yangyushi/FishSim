#include "network.hpp"

std::random_device rd;
std::mt19937 g(rd());

Graph::Graph() : nodes_{}, edges_{} {;}
Graph::Graph(Nodes n, Edges e) : nodes_{n}, edges_{e} {;}
Graph::Graph(int n) : nodes_(n){
    std::iota(nodes_.begin(), nodes_.end(), 0);
}
Graph::Graph(ConnMat adj_mat) : nodes_{}, edges_{}{
    int n = adj_mat.rows();
    for (int i = 0; i < n; i++){
        nodes_.push_back(i);
    }
    for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
        if (adj_mat(i, j) > 0){
            edges_.push_back(std::array<int, 2> {i, j});
        }
    }
    }
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


bool check_existence(int i, Nodes &collection){
    for (auto element : collection){
        if (i == element){
            return true;
        }
    }
    return false;
}

void collect_vnm_edges(int i, int d, int n, Edges& edges){
    Nodes neighbours {};
    int i2;
    for (int j = 0; j < d; j++){
        i2 = std::rand() % n;
        bool is_duplicate = check_existence(i2, neighbours);
        if (is_duplicate){
            j--;
        } else {
            edges.push_back(std::array<int, 2> {i, i2});
            neighbours.push_back(i2);
        }
    }
}


void collect_vnm_edges_force_diag(int i1, int d, int n, Edges& edges){
    Nodes neighbours {};
    bool self_included = false;
    for (int j = 0; j < d - 1; j++){
        int i2 = std::rand() % n;
        bool is_duplicate = check_existence(i2, neighbours);
        if (is_duplicate){
            j--;
        } else {
            edges.push_back(std::array<int, 2> {i1, i2});
            neighbours.push_back(i2);
            if (i2 == i1){ self_included = true; }
        }
    }
    if (self_included){
        bool should_continue = true;
        int i2;
        while (should_continue){
            should_continue = false;
            i2 = std::rand() % n;
            for (auto neighbour : neighbours) {
                if (i2 == neighbour) {should_continue = true;}
            }
        }
        edges.push_back(std::array<int, 2> {i1, i2});
    } else {
        edges.push_back(std::array<int, 2> {i1, i1});
    }
}


void collect_vnm_edges_no_diag(int i1, int d, int n, Edges& edges){
    Nodes neighbours {};
    for (int j = 0; j < d; j++){
        int i2 = std::rand() % n;
        bool is_duplicate = check_existence(i2, neighbours);
        is_duplicate = is_duplicate or (i1 == i2);
        if (is_duplicate){
            j--;
        } else {
            edges.push_back(std::array<int, 2> {i1, i2});
            neighbours.push_back(i2);
        }
    }
}


Graph random_vnm_graph(int d, int n){
    if (d == 0){
        return Graph(n);
    } else if (d == n) {
        ConnMat adj_mat {n, n};
        adj_mat.setConstant(1);
        return Graph(adj_mat);
    } else {
        Nodes nodes {};
        for (int i = 0; i < n; i++){
            nodes.push_back(i);
        }
        Edges edges{};
        if (d < n / 2){
            for (int i = 0; i < n; i++){
                collect_vnm_edges(i, d, n, edges);
            }
            return Graph(nodes, edges);
        } else {
            for (int i = 0; i < n; i++){
                collect_vnm_edges(i, n - d, n, edges);
            }
            Graph g_inv {nodes, edges};
            ConnMat adj_mat_inv = g_inv.as_matrix();
            return Graph (1 - adj_mat_inv); 
        }
    }
}


Graph random_vnm_graph_force_self(int d, int n){
    if (d == 0) {
        return Graph(n);
    } else if (d == n) {
        ConnMat adj_mat {n, n};
        adj_mat.setConstant(1);
        return Graph(adj_mat);
    } else {
        Nodes nodes {};
        for (int i = 0; i < n; i++){
            nodes.push_back(i);
        }

        Edges edges {};
        if (d < n / 2){
            for (int i = 0; i < n; i++){
                collect_vnm_edges_force_diag(i, d, n, edges);
            }
            return Graph{nodes, edges};
        } else {
            for (int i = 0; i < n; i++){
                collect_vnm_edges_no_diag(i, n - d, n, edges);
            }
            Graph g_inv {nodes, edges};
            ConnMat adj_mat_inv = g_inv.as_matrix();
            return Graph (1 - adj_mat_inv); 
        }
    }
}


Network3D::Network3D(int n, int k, double eta)
    : AVS3D{n, eta, 1.0}, k_(k) {
        graph_ = random_vnm_graph(k, n);
        connections_ = graph_.as_connections();
}


void Network3D::move(bool new_graph){
    if (new_graph){
        graph_ = random_vnm_graph(k_, n_);
        connections_ = graph_.as_connections();
    }
    vicsek_align(velocities_, connections_);
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

