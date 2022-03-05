#include "montecarlo.hpp"

int main(){
    FishMCMC system{50, 10, 1e-5, 1e-10, 0.05};
    system.evolve(100);
    for (size_t i = 0; i < 100; i++){
        system.sweep();
        std::cout << "E = " << system.get_total_energy() << std::endl;
    }
}
