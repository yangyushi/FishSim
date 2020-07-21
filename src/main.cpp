#include "vicsek.h"

int main(){
    int jump = 10;

    Vicsek3DPBCInertia system{
        50, 1.0, 0.4, 3.68, 0.01, 0.74
    };
    dump(system, "test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 1000; step++){
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
    for (int step=0; step < 1000; step++){
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
