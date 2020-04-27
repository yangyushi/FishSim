#include "vicsek.h"

int main(){
    Vicsek3DPBC system{1000, 1.0, 0.36, 8.0, 0.05};

    for (int step=0; step < 5000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 5 == 0) system.dump("test.xyz");
    }

    system.load("test.xyz");

    for (int step=0; step < 5000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 5 == 0) system.dump("test.xyz");
    }
    return 0;
}
