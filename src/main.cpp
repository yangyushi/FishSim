#include "vicsek.h"

int main(){
    Vicsek3DPBC system{1000, 1.0, 0.4, 10.0, 0.1};

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
    }

    //system.load("test.xyz");

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        //system.dump("test.xyz");
    }
    return 0;
}
