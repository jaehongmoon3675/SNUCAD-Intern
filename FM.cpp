#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"

int main(){
    int N, C; //the num of net and cell, respectively
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    read_hgr_area(C, CELL_array);

    for(int i = 1; i <= C; i++){
        printf("Cell num: %d, ", i);
        CELL_array[i].print_Cell();
    }

    return 0;
}