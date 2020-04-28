#include "neighbour_list.h"
#include <iostream>


int main(){
    Indices arr{1, 2, 3};
    for (auto item : product_3d(arr)) {
        for (auto num : item){
            cout << num << ", ";
        }
        cout << endl;
    };
}
