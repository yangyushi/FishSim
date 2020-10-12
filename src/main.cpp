#include "network.h"
#include "vicsek.h"
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

    cout << "system created" << endl;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        system.move();
        if (step % jump == 0) dump(system, "test.xyz");
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time spent without neighbour list: "
         << chrono::duration_cast<chrono::milliseconds>(end - begin).count()
         << "[ms]" << endl;

    cout << "movement works" << endl;

    load(system, "test.xyz");

    begin = chrono::steady_clock::now();
    for (int step=0; step < total_steps; step++){
        if (step % update_step == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) dump(system, "test.xyz");
    }
    end = chrono::steady_clock::now();
    cout << "Time spent without neighbour list: "
         << chrono::duration_cast<chrono::milliseconds>(end - begin).count()
         << "[ms]" << endl;
    return 0;
}
