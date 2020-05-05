#include "vicsek.h"

int main(){
    Attanasi2014PCB system{
        200, 1.0, 0.6, 0.02, 0.001,
    };

    cout << "system created" << endl;

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 20 == 0) system.dump("test.xyz");
    }

    cout << "movement works" << endl;

    for (int step=0; step < 10000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 20 == 0) system.dump("test.xyz");
    }
    return 0;
}
