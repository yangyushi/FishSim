#include "vicsek.h"

int main(){
    Vicsek2DPBC system{2000, 1.0, 0.05, 30.0, 0.05};

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
    }

    system.load("test.xyz");

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        system.dump("test.xyz");
    }
    return 0;
}
