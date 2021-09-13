#include "network.hpp"

std::random_device rd;
std::mt19937 g(rd());

Graph::Graph() : nodes_{}, edges_{} {;}

Graph::Graph(Nodes n, Edges e, bool directional) : nodes_{n}, edges_{e} {
    if (not directional){
        for (auto& edge : e){
            edges_.emplace(Index2D {edge[1], edge[0]});
        }
    }
}

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
            edges_.emplace(Index2D{i, j});
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


std::size_t PairHasher::operator()(const Index2D& a) const {
    std::size_t h = 0;
    for (auto e : a) {
        h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
}


void _sort_pair(int& i1, int& i2){
    if (i1 > i2) {
        int tmp = i1;
        i1 = i2;
        i2 = tmp;
    }
}

/* 
 * pairs in edges are sorted, also i1 shuold < i2
 *   the unordered_set is MUCH faster than std::vector
*/
bool _not_in_edges(int i1, int i2, const Edges& edges){
    auto query = edges.find(Index2D{i1, i2});
    return query == edges.end();
}


bool _is_suitable_for_graph(const Edges& edges, const std::unordered_map<int, int>& potential_edges) {
    if (potential_edges.size() == 0){
        return true;
    }
    int s1, s2;
    for (const auto& [k1, val_1] : potential_edges){
    for (const auto& [k2, val_2] : potential_edges){
        s1 = k1; s2 = k2;
        if (s1 == s2) break;
        _sort_pair(s1, s2);
        if (_not_in_edges(s1, s2, edges)) {
            return true;
        }
    }
    }
    return false;
}


Edges _try_create_edges(int d, int n){
    Edges edges;
    std::vector<int> stubs;
    std::unordered_map<int, int> potential_edges;
    int s1, s2, smax;

    for (int repeat = 0; repeat < d; repeat++){
        for (int i = 0; i < n; i++){
            stubs.push_back(i);
        }
    }

    int repeat_num = 0; 
    while (stubs.size() > 0){
        repeat_num++;
        potential_edges.clear();
        shuffle(stubs.begin(), stubs.end(), g);

        smax = stubs.size();
        for (int count = 0; count < smax; count++){
            s1 = stubs[count++];
            s2 = stubs[count];

            _sort_pair(s1, s2);

            if ((s1 != s2) and _not_in_edges(s1, s2, edges)){
                edges.emplace(Index2D{s1, s2});
            } else {
                potential_edges[s1]++;  //std::map will initialise value of 0
                potential_edges[s2]++;
            }
        }

        if (not _is_suitable_for_graph(edges, potential_edges)){
            return Edges {};
        }

        stubs.clear();
        for (const auto& [node, node_count] : potential_edges){
            for (int j = 0; j < node_count; j++){
                stubs.push_back(node);  // push back twice to to create pairs
            }
        }
    }
    return edges;
}


bool _check_existence(int i, Nodes &collection){
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
        bool is_duplicate = _check_existence(i2, neighbours);
        if (is_duplicate){
            j--;
        } else {
            edges.emplace(Index2D{i, i2});
            neighbours.push_back(i2);
        }
    }
}


void collect_vnm_edges_force_diag(int i1, int d, int n, Edges& edges){
    Nodes neighbours {};
    bool self_included = false;
    for (int j = 0; j < d - 1; j++){
        int i2 = std::rand() % n;
        bool is_duplicate = _check_existence(i2, neighbours);
        if (is_duplicate){
            j--;
        } else {
            edges.emplace(Index2D {i1, i2});
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
        edges.emplace(Index2D {i1, i2});
    } else {
        edges.emplace(Index2D {i1, i1});
    }
}


void collect_vnm_edges_no_diag(int i1, int d, int n, Edges& edges){
    Nodes neighbours {};
    for (int j = 0; j < d; j++){
        int i2 = std::rand() % n;
        bool is_duplicate = _check_existence(i2, neighbours);
        is_duplicate = is_duplicate or (i1 == i2);
        if (is_duplicate){
            j--;
        } else {
            edges.emplace(Index2D {i1, i2});
            neighbours.push_back(i2);
        }
    }
}


Edges get_regular_edges(int d, int n){
    Edges edges = _try_create_edges(d, n);
    while (edges.size() == 0){
        edges = _try_create_edges(d, n);
    }
    return edges;
}

Graph random_regular_graph(int d, int n){
    if (d == 0){
        return Graph{n};
    }
    if (n * d % 2 != 0){
        throw std::invalid_argument("n * d shuold be even");
    }
    if ((d < 0) or (d >= n)) {
        throw std::invalid_argument("the 0 <= d < n inequality must be satisfied");
    }
    Nodes nodes;
    for (int i = 0; i < n; i++){
        nodes.push_back(i);
    }
    if (d < n / 2){
        Edges edges = get_regular_edges(d, n);
        return Graph{nodes, edges, false};
    } else {
        Edges edges_inv = get_regular_edges(n-d, n);
        Graph g_inv{nodes, edges_inv, true};
        ConnMat adj_mat_inv = g_inv.as_matrix();
        return Graph{1 - adj_mat_inv};
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
            return Graph(nodes, edges, true);
        } else {
            for (int i = 0; i < n; i++){
                collect_vnm_edges(i, n - d, n, edges);
            }
            Graph g_inv {nodes, edges, true};
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
            return Graph{nodes, edges, true};
        } else {
            for (int i = 0; i < n; i++){
                collect_vnm_edges_no_diag(i, n - d, n, edges);
            }
            Graph g_inv {nodes, edges, true};
            ConnMat adj_mat_inv = g_inv.as_matrix();
            return Graph (1 - adj_mat_inv); 
        }
    }
}


Network3D::Network3D(int n, int k, double eta)
    : AVS3D{n, eta, 1.0}, k_(k) {
        this->update_graph();
}


void Network3D::move(bool new_graph){
    if (new_graph){
        this->update_graph();
    }
    vicsek_align(velocities_, connections_);
    normalise(velocities_);
    add_noise();
}


void Network3D::update_graph(){
    graph_ = random_vnm_graph(k_, n_);
    connections_ = graph_.as_connections();
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
        throw std::invalid_argument(
            "invalid dynamic type, choose between [0]: quenched or [1]: anneled"
        );
    }
}


ConnMat Network3D::get_adj_mat(){
    return graph_.as_matrix();
}


void Network3D::set_adj_mat(ConnMat adj_mat){
    Graph new_graph{adj_mat};
    graph_ = new_graph;
}


Network3DRG::Network3DRG(int n, int k, double eta) : Network3D(n, k, eta){}

void Network3DRG::update_graph(){
    graph_ = random_regular_graph(k_ - 1, n_);
    for (auto node : graph_.nodes_){
        graph_.edges_.emplace(Index2D{node, node});
    }
    connections_ = graph_.as_connections();
}
