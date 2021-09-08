#include "network.hpp"


int main(){
    int n = 50;
    int k = 10;
    int total_steps = 100;

    Network3D system{n, k, 0.3};

    std::cout << "system created" << std::endl;

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
