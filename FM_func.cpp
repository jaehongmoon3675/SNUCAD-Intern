#include <iostream>
#include "Cell.h"
#include "Net.h"

void getMax(const int C, const Cell* CELL_array, int &pmax, int &smax){
    pmax = -1;
    smax = -1;

    pmax = CELL_array[1].get_pin();
    smax = CELL_array[1].get_size();

    for(int i = 2; i <= C; i++){    
        if(pmax < CELL_array[i].get_pin())
            pmax = CELL_array[i].get_pin();

        if(smax < CELL_array[i].get_size())
            smax = CELL_array[i].get_size();
    }
}