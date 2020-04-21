#include "vicsek.h"

int main(){
    Vicsek3DPBC system{50, 1.0, 0.05, 7.0, 0.1};
    system.move(true);
    for (int step=0; step < 20000; step++){
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
