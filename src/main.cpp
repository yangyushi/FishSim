#include "vicsek.h"

int main(){
    int jump = 1;

    Attanasi2014PCB system{
        100, 1.0, 0.05, 0.02, 0.03,
    };

    cout << "system created" << endl;

    for (int step=0; step < 10000; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
    }

    cout << "movement works" << endl;

    for (int step=0; step < 100; step++){
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
