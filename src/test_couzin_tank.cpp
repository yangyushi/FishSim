#include "couzin.hpp"
#include <chrono>

int main(){
    int n = 20;
    int total_steps = 10000;
    int jump = 100;
    float theta = 3.141592653 / 180 * 10;
    float perception = 270.0 / 360.0 * 3.1415;
    double c = 0.734, h = 0.35;
    double kw = 0.001;
    bool align_cap = false;
    bool align_base = false;
    double g = 0.01;

    CouzinTank3D system{
        n,
        0.03, 0.05, 0.2, // repel, align, attract
        perception,  // perception
        0.05,  // eta
        0.01, // v0
        theta, // turning rate
        0.05, // dt
        c, h, kw, align_cap, align_base, g
    };

    system.positions_ *= 0.9;

    std::cout << "system created" << std::endl;

    dump(system, "test_couzin_tank.xyz");

    std::cout << "dump okay" << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move_in_tank(true);
        if (step % jump == 0) dump(system, "test_couzin_tank.xyz");
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time spent without neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;

    std::cout << "movement works" << std::endl;

    load(system, "test_couzin_tank.xyz");

    begin = std::chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move_in_tank(true);
        if (step % jump == 0) dump(system, "test_couzin_tank.xyz");
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Time spent with neighbour list: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
         << "[ms]" << std::endl;
    for (int step=0; step < total_steps; step++){
        system.move_in_tank(true);
        if (step % jump == 0) {
            std::cout << "mill: " << system.get_mill()
                      << "; polarisation: " << system.get_polarisation()
                      << std::endl;
        }
    }

    return 0;
}
