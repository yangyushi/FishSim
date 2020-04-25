#include "vicsek.h"

int main(){
    Vicsek2DPBCVN system{200, 1.0, 0.62, 10.0, 0.05};

    for (int step=0; step < 1000; step++){
        if (step % 100 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        system.dump("test.xyz");
    }

    system.load("test.xyz");

    for (int step=0; step < 1000; step++){
        if (step % 100 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        system.dump("test.xyz");
    }

    return 0;
}
