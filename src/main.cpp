#include "vicsek.h"

int main(){
    int jump = 1;

    Vicsek3DPBCInertia system{
        250, 1.0, 0.4, 5, 0.05, 0.8,
    };

    cout << "system created" << endl;

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) system.dump("test.xyz");
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
