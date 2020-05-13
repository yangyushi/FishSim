#include "vicsek.h"

int main(){
    int jump = 1;

    Attanasi2014PCB system{
        100, 1.0, 0.05, 0.02, 0.03,
    };
    system.dump("test.xyz");

    cout << "system created" << endl;

    for (int step=0; step < 500; step++){
        if (step % 20 == 0){
            system.move(true);
        }
        else {
            system.move(false);
        }
        if (step % jump == 0) system.dump("test.xyz");
    }

    cout << "movement works" << endl;

    for (int step=0; step < 500; step++){
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
