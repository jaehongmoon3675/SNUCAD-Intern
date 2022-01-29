#include <iostream>
#include <fstream>
#include <list>
#include <stack>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "WriteFile.h"

extern std::stack<Cell*> FreeCellList;

void read_output_part(Block &A, Block &B, const int C, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("output.part");

    //printf("read_hgr_map\n");

    int block;
    std::string temp_cell_name;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_name >> block;
            
            if(block == 1)
                CELL_array[i].set_current_block(&A);
            else
                CELL_array[i].set_current_block(&B);

            FreeCellList.push(CELL_array + i);
        }

        readFile.close();
    }
    else
        printf("No map file\n");
}

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    int pass = 20; //how many pass we go through
    int min_cutset_num;
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);
    balance_low_bound = r*W - smax * 10;
    balance_up_bound = r*W + smax * 10;

    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d\n", pmax, smax);

    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "A");
    Block B(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "B");

    read_output_part(A, B, C, CELL_array);
    BlockReinitialization(A, B, CELL_array);


    printf("cutnet_num: %d\n", CountCutNet(A, NET_array, N));


    delete[] NET_array;
    delete[] CELL_array;

    return 0;
}