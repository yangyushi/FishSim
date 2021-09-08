#include "network.hpp"
#include "vicsek.hpp"
#include <chrono>

int main(){
    int n = 1000;
    int total_steps = 100;
    int jump = 1;

    Vicsek2DPBC system{
        n, 1.0, 0.3, 20, 0.05
    };

    std::cout << "system created" << std::endl;

    dump(system, "test.xyz");

    std::cout << "dump okay" << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move_no_nl();
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
        system.move(true);
        if (step % jump == 0) dump(system, "test.xyz");
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Time spent with neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;
    return 0;
}
