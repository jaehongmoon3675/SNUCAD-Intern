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

int main(int argc, char ** argv){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    const int pass = 1; //how many pass we go through
    const int k = 1;
    const int InitVer = 2;
    int global_min_cutnet, local_min_cutnet;
    int cutnet;
    const int stuck_criteria = 4;
    const bool balance_option = false; //true면 smax 기반, false면 비율 기반
    const bool time_option = true;
    const bool stuck_option = true;
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

    global_min_cutnet = CountCutNet(A, NET_array, N);
    local_min_cutnet = global_min_cutnet;
    cutnet = global_min_cutnet;

    CellDist GlobalMinDist(C, W*r, cutnet, A.get_size(), B.get_size(), &A, &B);
    CellDist LocalMinDist(C, W*r, cutnet, A.get_size(), B.get_size(), &A, &B);

    printf("initial cutnet num: %d\n", global_min_cutnet);
    
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
    int stuck_count = 0;
    bool stuck_check = true;
    bool stuck = false;

    bool pass_start = true;
    printf("FM start\n");

    if(time_option){
        end = clock();
        total_time = (double)end - start;
        printf("time: %fs\n\n", total_time / CLOCKS_PER_SEC);

        move_time = clock();
        reinit_time = clock();
    }

    double move_min = 5, move_max = 0;
    int move_max_loc = 0;

    for(int i = 0; i < pass; i++){
        if(!stuck)
            BaseCell = ChooseBaseCell_gain(A, B, r);
        else
            BaseCell = ChooseBaseCell_balance(A, B, r);
        //printf("Pass %d starts\n", i);
        //BaseCell->print_Cell();

        pass_start = true;
        stuck_check = true;
        move_count = 0;

        

        while(BaseCell != nullptr){
            clock_t Move_time = clock();
            if(BaseCell->get_current_block() == &A){
                cutnet -= A.gain[BaseCell->get_cell_num()];
                MoveCell(A, B, BaseCell);
            }
            else{
                cutnet -= B.gain[BaseCell->get_cell_num()];
                MoveCell(B, A, BaseCell);
            }

            double temp = clock() - (double)Move_time;

            if(move_max < temp){
                move_max = temp;
                move_max_loc = move_count;
            }

            if(move_count == 1960)
                printf("hello\n");

            
            if(temp < move_min)
                move_min = temp;

            if(time_option){
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
                //printf("cutnet update: %d\n", global_min_cutnet);
                global_min_cutnet = cutnet;
                stuck_check = false;
                //printf("pass: %d\n", i);
                //stuck = false;
            }

            
            
            if(pass_start){
                LocalMinDist.overWrite(CELL_array, C, A.get_size(), B.get_size(), cutnet);    
                pass_start = false;
            }
            else{
                /*
                if((i/5) % 4 != 0 && cutnet >= global_min_cutnet){
                //if((i/8) % 10 != 0 && cutnet >= global_min_cutnet){

                }
                else if(LocalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet))
                    stuck = false;
                */

                //if(LocalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet))
                //    stuck_check = false;

                LocalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet);
            }
            
            //printf("expected: %d, real: %d\n", cutnet, CountCutNet(A, NET_array, N));
            assert(CountCutNet(A, NET_array, N) == cutnet);
            //printf("num of cutnet: %d\n", temp);

            if(!stuck)
                BaseCell = ChooseBaseCell_gain(A, B, r);
            else
                BaseCell = ChooseBaseCell_balance(A, B, r);
            /*
            if(BaseCell != nullptr)
                BaseCell->print_Cell();
            */
            //printf("go on...!\n");
        }

        if(time_option){  
            printf("Pass %d finished, Total move: %d\n", i, move_count);
            printf("min cutnet: %d\n", global_min_cutnet);
            total_time = clock() - (double)end;
            end = clock();
            printf("time: %fs\n\n", total_time / CLOCKS_PER_SEC);

            reinit_time = clock();
        }

        LoadDistribution(LocalMinDist, CELL_array, C, cutnet);
        BlockReinitialization(A, B, CELL_array);

        if(time_option){
            total_time = clock() - (double)reinit_time;
            move_time = clock();
            printf("reinit time: %fs\n", total_time / CLOCKS_PER_SEC);
        }

        if(!stuck && stuck_check)
            stuck_count++;
        
        if(stuck)
            stuck_count--;

        if(stuck_option){
            if(stuck_count == stuck_criteria){
                stuck = true;
            }
        }

        if(stuck_count == 0)
            stuck = false;

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

    printf("Move time min: %f, Move time max: %f\n", move_min, move_max);
    printf("Move max loc: %d\n", move_max_loc);


    printf("\n---------- After FM Algorithm----------\n");
    

    GlobalMinDist.writeCellDist(CELL_array, C);
    GlobalMinDist.printCellDist();

    //printf("Final min num of cutnet: %d\n", global_min_cutnet);

    Check(A, B, NET_array, N);


    delete[] NET_array;
    delete[] CELL_array;

    end = clock();
    total_time = (double)end - start;

    printf("\ntime: %fs\n", total_time / CLOCKS_PER_SEC);


    return 0;
}