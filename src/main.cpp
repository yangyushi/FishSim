#include "spin.h"
#include "vicsek.h"

int main(){
    int jump = 1;

    InertialSpin3DPBC system{
        50, 5, 1.0, 0.1,
        1e-4, 0.8, 1.25, 1
    };
    dump(system, "test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 500; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) dump(system, "test.xyz");
    }

    cout << "movement works" << endl;

    load(system, "test.xyz");
    for (int step=0; step < 500; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) dump(system, "test.xyz");
    }
    return 0;
}
