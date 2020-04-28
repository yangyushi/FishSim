#include "vicsek.h"

int main(){
    Vicsek2DPBC system{200, 1.0, 0.40, 8.0, 0.1};

    for (int step=0; step < 10000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 5 == 0) system.dump("test.xyz");
    }

    system.load("test.xyz");

    for (int step=0; step < 10000; step++){
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
