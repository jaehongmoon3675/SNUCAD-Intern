#include <iostream>
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
    int k = 1;
    int min_cutset_num;
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);
    balance_low_bound = r*W - smax * k;
    balance_up_bound = r*W + smax * k;

    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d\n", pmax, smax);

    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "A");
    Block B(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "B");

    BlockInitialization(A, B, CELL_array, C);
    BlockReinitialization(A, B, CELL_array);

    min_cutset_num = CountCutNet(A, NET_array, N);
    printf("initial cutset num: %d\n", min_cutset_num);
    write_output(A, CELL_array, C);
    //printCellInfo(CELL_array, C);

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

    for(int i = 0; i < pass; i++){
        BaseCell = ChooseBaseCell(A, B, r);
        //BaseCell->print_Cell();

        while(BaseCell != nullptr){
            /*
            printf("before move...\n");
            printf("\nBlock A\n");
            A.print_Block(CELL_array);

            printf("\nBlock B\n");
            B.print_Block(CELL_array);
            printf("\n\n\n\n");
            */
            if(BaseCell->get_current_block() == &A)
                MoveCell(A, B, BaseCell);
            else
                MoveCell(B, A, BaseCell);
            /*
            printf("after move...\n");
            printf("\nBlock A\n");
            A.print_Block(CELL_array);

            printf("\nBlock B\n");
            B.print_Block(CELL_array);
            printf("\n\n\n\n");
            */
            
            temp = CountCutNet(A, NET_array, N);
            //printf("num of cutset: %d\n", temp);

            if(min_cutset_num > temp){
                min_cutset_num = temp;
                //CountCutNetAgain(A, NET_array, N);

                printf("min chagne: %d\n", temp);
                write_output(A, CELL_array, C);
            }
            BaseCell = ChooseBaseCell(A, B, r);
            /*
            if(BaseCell != nullptr)
                BaseCell->print_Cell();
            */
            //printf("go on...!\n");
        }

        
        //printf("Pass %d finished\n", i);
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
    printf("Final min num of cutset: %d\n", min_cutset_num);

    Check(A, B, NET_array, N);


    delete[] NET_array;
    delete[] CELL_array;

    return 0;
}