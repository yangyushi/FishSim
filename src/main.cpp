#include "spin.h"
#include "vicsek.h"

int main(){
    int jump = 100;

    InertialSpin3D system{
        50, 5, 0.1, 8e-2, 0.8, 1.25, 0.3
    };
    dump(system, "test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 50000; step++){
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
    for (int step=0; step < 50000; step++){
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
