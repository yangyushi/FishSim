#include "vicsek.h"

int main(){
    Vicsek2DPBCVNCO system{
        200, 1.0, 1.0, 5.0, 0.05,
        5.0, 3.0, 0.8, 0.5, 0.2
    };

    for (int step=0; step < 100000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
    }

    for (int step=0; step < 100000; step++){
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
