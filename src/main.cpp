#include "vicsek.h"

int main(){
    int jump = 10;

    Vicsek3DPBCInertia system{
        50, 1.0, 0.8, 3.68, 5, 0.74
    };
    system.dump("test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 10000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
    }

    cout << "movement works" << endl;

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) system.dump("test.xyz");
    }
    return 0;
}
