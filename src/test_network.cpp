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
    std::cout << "Network3D system creation okay" << std::endl;

    s1.evolve(steps, 0);
    std::cout << "Network3D quenched movement okay" << std::endl;

    for (int step=0; step < steps; step++){
        s1.move(true);
        std::cout << "Polarisation: " << s1.get_polarisation() << std::endl;
    }
    std::cout << "Network3D anneled movement okay" << std::endl;

    /*
     * Test the vectorial network model on a regular graph
     */
    Network3DRG s2{n, k, 0.3};
    std::cout << "Network3DRG system creation okay" << std::endl;

    s2.evolve(steps, 0);
    std::cout << "Network3DRG quenched movement okay" << std::endl;

    for (int step=0; step < steps; step++){
        s2.move(true);
        std::cout << "Polarisation: " << s2.get_polarisation() << std::endl;
    }
    std::cout << "Network3DRG anneled movement okay" << std::endl;


    /*
     * Test the voter model
     */
    Voter s3{n, k, 0.2};
    std::cout << "Voter system creation okay" << std::endl;

    s3.evolve(steps, 0);
    std::cout << "Voter quenched movement okay" << std::endl;

    for (int step=0; step < steps; step++){
        s3.move(true);
        std::cout << "Polarisation: " << s3.get_polarisation() << std::endl;
    }
    std::cout << "Voter anneled movement okay" << std::endl;


    return 0;
}
