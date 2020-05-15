#include "vicsek.h"

int main(){
    int jump = 10;

    Vicsek3DPBC system{
        200, 1.0, 0.01, 10, 0.05,
    };
    system.dump("test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) system.dump("test.xyz");
    }

    cout << "movement works" << endl;

    for (int step=0; step < 1000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) system.dump("test.xyz");
    }
    return 0;
}
