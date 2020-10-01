#include "network.h"
#include "vicsek.h"

int main(){
    int jump = 1;

    Vicsek3DPBC system{
        3, 1.0, 0.2, 5, 0.5
    };
    dump(system, "test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 2; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) dump(system, "test.xyz");
    }

    cout << "movement works" << endl;

    load(system, "test.xyz");
    for (int step=0; step < 2; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) dump(system, "test.xyz");
    }
    return 0;
}
