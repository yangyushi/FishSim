#include "network.hpp"


int main(){
    int n = 50;
    int k = 2;
    int steps = 10;

    Graph g; 
    g = random_regular_graph(k, n);
    g = random_vnm_graph(k, n);
    g = random_vnm_graph_force_self(k, n);
    std::cout << "graph functions okay" << std::endl;


    /*
     * Test the vectorial network model
     */
    Network3D s1{n, k, 0.3};
    s1.evolve(steps, 0);

    for (int step=0; step < steps; step++){
        s1.move(true);
    }

    /*
     * Test the vectorial network model on a regular graph
     */
    Network3DRG s2{n, k, 0.3};
    s2.evolve(steps, 0);
    for (int step=0; step < steps; step++){
        s2.move(true);
    }

    /*
     * Test the voter model
     */
    Voter s3{n, k, 0.2};
    s3.evolve(steps, 0);
    for (int step=0; step < steps; step++){
        s3.move(true);
    }
    return 0;
}
