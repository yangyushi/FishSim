#include "vicsek.h"

int main(){
    Vicsek2DPBC system{500, 1.0, 0.1, 10.0, 0.05};

    for (int step=0; step < 2000; step++){
        if (step % 100 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 10 == 0){
            system.dump("test.xyz");
        }
    }

    system.load("test.xyz");

    for (int step=0; step < 2000; step++){
        if (step % 100 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 10 == 0){
            system.dump("test.xyz");
        }
    }

    return 0;
}
