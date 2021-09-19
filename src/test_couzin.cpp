#include "couzin.hpp"
#include <chrono>

int main(){
    int n = 50;
    int total_steps = 1000;
    int jump = 1;
    float theta = 3.141592653 / 180 * 40;
    float perception = 270.0 / 360.0 * 3.1415;

    Couzin3D system{
        n,
        1.0, 2.0, 20.0, // repel, align, attract
        perception,  // perception
        0.05,  // eta
        3.0, // v0
        theta, // turning rate
        0.1 // dt
    };

    system.positions_ *= 1;

    std::cout << "system created" << std::endl;

    dump(system, "test_couzin.xyz");

    std::cout << "dump okay" << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move(true);
        if (step % jump == 0) dump(system, "test_couzin.xyz");
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time spent without neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;

    std::cout << "movement works" << std::endl;

    load(system, "test_couzin.xyz");

    begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move(true);
        if (step % jump == 0) dump(system, "test_couzin.xyz");
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Time spent with neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;
    return 0;
}
