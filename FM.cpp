#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    int pass = 5; //how many pass we go through
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;

    get_max(C, CELL_array, pmax, smax);
    balance_low_bound = r*W - smax;
    balance_up_bound = r*W + smax;

    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d\n", pmax, smax);

    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "A");
    Block B(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "B");

    BlockInitialization(A, B, CELL_array, C);
    BlockReinitialization(A, B, CELL_array);

    //printCellInfo(CELL_array, C);

    printf("Block A\n");
    A.print_Block();

    printf("\nBlock B\n");
    B.print_Block();


    delete[] NET_array;
    delete[] CELL_array;

    return 0;
}