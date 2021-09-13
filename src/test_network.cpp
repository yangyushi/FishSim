#include "network.hpp"


int main(){
    int n = 50;
    int k = 2;
    int total_steps = 100;

    Graph g; 
    g = random_regular_graph(k, n);
    g = random_vnm_graph(k, n);
    g = random_vnm_graph_force_self(k, n);
    std::cout << "graph functions okay" << std::endl;

    Network3DRG system{n, k, 0.3};
    std::cout << "system creation okay" << std::endl;

    for (int step=0; step < total_steps; step++){
        system.move(false);
    }

    std::cout << "quenched movement okay" << std::endl;

    for (int step=0; step < total_steps; step++){
        system.move(true);
    }

    std::cout << "anneled movement okay" << std::endl;

    return 0;
}
