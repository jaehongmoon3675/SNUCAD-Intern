#define NDEBUG

#include <iostream>
#include <cassert>
#include <list>
#include <ctime>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "CellDist.h"
#include "WriteFile.h"

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    int pass = 5; //how many pass we go through
    int k = 1;
    int InitVer = 2;
    int min_cutnet;
    int cutnet;
    bool balance_option = false; //true면 smax 기반, false면 비율 기반
    bool time_check = false;
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    clock_t start, end, move_time, reinit_time;
    double total_time;

    start = clock();

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);

    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d\n", pmax, smax);
    //printf("balance ratio: %f, num of pass: %d, k: %d\n", r, pass, k);

    if(balance_option){
        printf("balance ratio: %f, num of pass: %d, k: %d\n", r, pass, k);
        balance_low_bound = r*W - smax * k;
        balance_up_bound = r*W + smax * k;
    }
    else{
        printf("balance ratio +- 0.05, num of pass: %d\n", pass);
        balance_low_bound = (r - 0.05)*W;
        balance_up_bound = (r + 0.05)*W;
    }

    printf("balance low bound: %f, balance up bound: %f\n", balance_low_bound, balance_up_bound);
    
    
    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "A");
    Block B(pmax, W - balance_up_bound, W - balance_low_bound, C, N, W, r, "B");


    switch(InitVer){
    case 1:
        BlockInitialization(A, B, CELL_array, C);
        break;
    case 2:
        BlockInitialization(A, B, CELL_array, NET_array, C, N);
        break;
    case 3:
        BlockInitialization(A, B, CELL_array, NET_array, C, N, 0);
        break;
    default:
        BlockInitialization(A, B, CELL_array, NET_array, C, N, 0);
    }

    BlockReinitialization(A, B, CELL_array);

    min_cutnet = CountCutNet(A, NET_array, N);
    cutnet = min_cutnet;

    CellDist GlobalMinDist(C, W*r, cutnet, A.get_size(), B.get_size(), &A, &B);
    CellDist LocalMinDist(C, W*r, cutnet, A.get_size(), B.get_size(), &A, &B);

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
    int move_count = 0;

    bool pass_start = true;
    bool stuck = false;
    printf("FM start\n");

    if(time_check){
        end = clock();
        total_time = (double)end - start;
        printf("time: %fs\n\n", total_time / CLOCKS_PER_SEC);

        move_time = clock();
        reinit_time = clock();
    }

    for(int i = 0; i < pass; i++){
        BaseCell = ChooseBaseCell_gain(A, B, r);
        //printf("Pass %d starts\n", i);
        //BaseCell->print_Cell();

        pass_start = true;
        stuck = false;
        move_count = 0;

        while(BaseCell != nullptr){
            if(BaseCell->get_current_block() == &A){
                cutnet -= A.gain[BaseCell->get_cell_num()];
                MoveCell(A, B, BaseCell);
            }
            else{
                cutnet -= B.gain[BaseCell->get_cell_num()];
                MoveCell(B, A, BaseCell);
            }
            
            if(time_check){
                if(move_count % 10000 == 5000){
                    total_time = clock() - (double)move_time;
                    move_time = clock();
                    printf("move time: %fs\n\n", total_time / CLOCKS_PER_SEC);

                    clock_t move_time = clock();
                }
                move_count++;
            }

            /*
            printf("Move %d\n", move_count++);
            printf("\nBlock A\n");
            A.print_Block_short(CELL_array);

            printf("\nBlock B\n");
            B.print_Block_short(CELL_array);
            printf("\n\n\n");
            */
            //printf("move!\n");
            
            if(GlobalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet)){
                //printf("cutnet update: %d\n", min_cutnet);
                min_cutnet = cutnet;
                //printf("pass: %d\n", i);
                //stuck = false;
            }

            
            
            if(pass_start){
                LocalMinDist.overWrite(CELL_array, C, A.get_size(), B.get_size(), cutnet);    
                pass_start = false;
            }
            else{
                /*
                if((i/5) % 4 != 0 && cutnet >= min_cutnet){
                //if((i/8) % 10 != 0 && cutnet >= min_cutnet){

                }
                else if(LocalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet))
                    stuck = false;
                */

                if(LocalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet))
                    stuck = false;
                
            }
            
            //printf("expected: %d, real: %d\n", cutnet, CountCutNet(A, NET_array, N));
            assert(CountCutNet(A, NET_array, N) == cutnet);
            //printf("num of cutnet: %d\n", temp);

            BaseCell = ChooseBaseCell_gain(A, B, r);
            /*
            if(BaseCell != nullptr)
                BaseCell->print_Cell();
            */
            //printf("go on...!\n");
        }

        if(stuck){
            printf("stuck!\n");
            break;
        }

        if(time_check){  
            printf("Pass %d finished, Total move: %d\n", i, move_count);
            printf("min cutnet: %d\n", min_cutnet);
            total_time = clock() - (double)end;
            end = clock();
            printf("time: %fs\n\n", total_time / CLOCKS_PER_SEC);

            reinit_time = clock();
            LoadDistribution(LocalMinDist, CELL_array, C, cutnet);
            BlockReinitialization(A, B, CELL_array);
            total_time = clock() - (double)reinit_time;
            move_time = clock();
            printf("reinit time: %fs\n", total_time / CLOCKS_PER_SEC);
        }
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


    printf("\n---------- After FM Algorithm----------\n");
    

    GlobalMinDist.writeCellDist(CELL_array, C);
    GlobalMinDist.printCellDist();

    //printf("Final min num of cutnet: %d\n", min_cutnet);

    Check(A, B, NET_array, N);


    delete[] NET_array;
    delete[] CELL_array;

    end = clock();
    total_time = (double)end - start;

    printf("\ntime: %fs\n", total_time / CLOCKS_PER_SEC);


    return 0;
}