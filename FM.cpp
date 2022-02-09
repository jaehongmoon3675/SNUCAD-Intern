//#define NDEBUG

#include <iostream>
#include <cassert>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "WriteFile.h"

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    int pass = 10; //how many pass we go through
    int k = 10;
    int min_cutnet;
    int cutnet;
    bool balance_option = false; //true면 smax 기반, false면 비율 기반
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);

    if(balance_option){
        printf("balance ratio: %f, num of pass: %d, k: %d\n", r, pass, k);
        balance_low_bound = r*W - smax * k;
        balance_up_bound = r*W + smax * k;
    }
    else{
        printf("balance ratio +- 0.1: %f, num of pass: %d\n", r, pass);
        balance_low_bound = (r - 0.05)*W;
        balance_up_bound = (r + 0.05)*W;
    }

    /*
    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d\n", pmax, smax);
    printf("balance low bound: %f, balance up bound: %f\n", balance_low_bound, balance_up_bound);
    printf("balance ratio: %f, num of pass: %d, k: %d\n", r, pass, k);
    */
    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "A");
    Block B(pmax, W - balance_up_bound, W - balance_low_bound, C, N, W, r, "B");

    //ver1
    //BlockInitialization(A, B, CELL_array, C);
    //ver2
    BlockInitialization(A, B, CELL_array, NET_array, C, N);
    BlockReinitialization(A, B, CELL_array);

    min_cutnet = CountCutNet(A, NET_array, N);
    cutnet = min_cutnet;

    CellDist MinDist(C, W*r, cutnet, A.get_size(), B.get_size(), A.get_cell_count(), B.get_cell_count(), &A, &B);

    printf("initial cutnet num: %d\n", min_cutnet);
    
    //printCellInfo(CELL_array, C);
    //printNetInfo(NET_array, N);

    /*
    printf("start!\n");
    printf("Block A\n");
    A.print_Block(CELL_array);

    printf("\nBlock B\n");
    B.print_Block(CELL_array);
    printf("\n\n");
    */

    Cell* BaseCell = nullptr;
    int temp;
    int move_count = 1;

    for(int i = 0; i < pass; i++){
        BaseCell = ChooseBaseCell(A, B, r);
        //BaseCell->print_Cell();

        while(BaseCell != nullptr){
            if(BaseCell->get_current_block() == &A){
                MoveCell(A, B, BaseCell);
                cutnet -= A.gain[BaseCell->get_cell_num()];
            }
            else{
                MoveCell(B, A, BaseCell);
                cutnet -= A.gain[BaseCell->get_cell_num()];
            }

            /*
            printf("Move %d\n", move_count++);
            printf("\nBlock A\n");
            A.print_Block_short(CELL_array);

            printf("\nBlock B\n");
            B.print_Block_short(CELL_array);
            printf("\n\n\n");
            */

            assert(CountCutNet(A, NET_array, N) == cutnet);
            //printf("num of cutnet: %d\n", temp);


            if(MinDist.update(CELL_array, C, A.get_size(), B.get_size(), A.get_cell_count(), B.get_cell_count(), cutnet))
                min_cutnet = cutnet;

            BaseCell = ChooseBaseCell(A, B, r);
            /*
            if(BaseCell != nullptr)
                BaseCell->print_Cell();
            */
            //printf("go on...!\n");
        }

        
        //printf("Pass %d finished\n", i);
        LoadDistribution(MinDist, CELL_array, C);
        BlockReinitialization(A, B, CELL_array);
        /*
        printf("After Pass %d...", i + 1);
        printf("Block A\n");
        A.print_Block(CELL_array);

        printf("\nBlock B\n");
        B.print_Block(CELL_array);
        printf("\n\n\n\n");
        */
    }
    /*
    printf("\n\n\n\nAfter Pass...\n");
    printf("Block A\n");
    A.print_Block(CELL_array);

    printf("\nBlock B\n");
    B.print_Block(CELL_array);
    */

    MinDist.writeCellDist(CELL_array, C);

    printf("Final min num of cutnet: %d\n", min_cutnet);

    Check(A, B, NET_array, N);


    delete[] NET_array;
    delete[] CELL_array;

    return 0;
}