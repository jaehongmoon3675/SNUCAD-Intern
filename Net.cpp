#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"

void printNetInfo(Net* NET_array, int N){
    printf("Net Info\n");
    for(int i = 1; i <= N; i++){
        NET_array[i].print_Net();
    }
}