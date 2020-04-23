#include "vicsek.h"

int main(){
    Vicsek3DPBC system{500, 1.0, 0.01, 40.0, 0.01};

    for (int step=0; step < 10000; step++){
        if (step % 100 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 100 == 0){
            system.dump("test.xyz");
        }
    }

    system.load("test.xyz");

    for (int step=0; step < 10000; step++){
        if (step % 100 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % 100 == 0){
            system.dump("test.xyz");
        }
    }

    return 0;
}
