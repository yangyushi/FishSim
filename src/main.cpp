#include "network.hpp"
#include "vicsek.hpp"
#include <chrono>

int main(){
    int n = 500;
    int total_steps = 1000;
    int jump = 1;
    int update_step = 100;

    Vicsek2DPBC system{
        n, 1.0, 0.4, 10, 0.1
    };
    dump(system, "test.xyz");

    std::cout << "system created" << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move();
        if (step % jump == 0) dump(system, "test.xyz");
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time spent without neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;

    std::cout << "movement works" << std::endl;

    load(system, "test.xyz");

    begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) dump(system, "test.xyz");
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Time spent without neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;
    return 0;
}
