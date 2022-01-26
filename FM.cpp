#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    int r = 0.5; //balance factor
    int pass = 5; //how many pass we go through
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array);

    /*
    printf("Cell Info\n");
    for(int i = 1; i <= C; i++){
        printf("Cell num: %d, ", i);
        CELL_array[i].print_Cell();
    }
    */
    
    int pmax = 0, smax = 0;
    int balance_low_bound = 0, balance_up_bound = 0;

    get_max(C, CELL_array, pmax, smax);
    balance_low_bound = r*W + smax;
    balance_up_bound = r*W + smax;

    return 0;
}