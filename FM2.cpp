#define NDEBUG
#define MAX_MOVE_TIME 10713

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

int gain_update_count = 0;

int main(int argc, char ** argv){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    const int pass = 10; //how many pass we go through
    const int k = 1;
    const int InitVer = 2;
    int global_min_cutnet, local_min_cutnet;
    int cutnet;
    const int stuck_criteria = 5;
    const bool stuck_out_itr = 1;
    const bool balance_option = false; //true면 smax 기반, false면 비율 기반
    const bool time_option = true;
    const bool stuck_option = true;
    const bool no_large_net = false;
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    clock_t start, end, move_time, reinit_time, Move_time;
    double total_time;

    start = clock();

    P = read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    W = read_hgr_area(C, CELL_array); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);

    int max_cell_count = get_max_cell_count(NET_array, N);

    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d, max_cell_count: %d\n", pmax, smax, max_cell_count);
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

    BlockReinitialization(A, B, CELL_array, NET_array, false);

    global_min_cutnet = CountCutNet(A, NET_array, N);
    local_min_cutnet = global_min_cutnet;
    cutnet = global_min_cutnet;

    if(no_large_net){
        A.deactivate_large_net(NET_array);
        B.deactivate_large_net(NET_array);
    }

    CellDist GlobalMinDist(C, N, W*r, cutnet, A.get_size(), B.get_size(), &A, &B);
    CellDist LocalMinDist(C, N, W*r, cutnet, A.get_size(), B.get_size(), &A, &B);

    printf("initial cutnet num: %d\n", global_min_cutnet);

    Cell* BaseCell = nullptr;
    int temp;
    int move_count = 0;
    int stuck_count = 0;
    int gain_update_count_temp = 0;
    bool stuck_check = true; //얼마나 오래동안 stuck 되었는지를 체크하기 위함
    bool stuck = false; //실제 stuck 되어있는지 여부
    bool destroy_balance = false;

    bool pass_start = true;

    if(time_option){
        printf("FM start\n");
        end = clock();
        total_time = (double)end - start;
        printf("time: %fs\n\n", total_time / CLOCKS_PER_SEC);

        move_time = clock();
        reinit_time = clock();
    }

    double move_min = 5, move_max = 0;
    int move_max_loc = 0;

    double record_time = 0;
    clock_t record_time_temp;

    for(int i = 0; i < pass; i++){
        record_time = 0;

        if(!stuck)
            BaseCell = ChooseBaseCell_gain(A, B, r);
        else{
            BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
        }

        pass_start = true;
        stuck_check = true;        

        while(BaseCell != nullptr){
            if(time_option)
                Move_time = clock();            

            if(BaseCell->get_current_block() == &A){
                cutnet -= A.gain[BaseCell->get_cell_num()];
                MoveCell(A, B, BaseCell);
            }
            else{
                cutnet -= B.gain[BaseCell->get_cell_num()];
                MoveCell(B, A, BaseCell);
            }       
            
            if(pass_start){
                LocalMinDist.overWrite(CELL_array, C, A.get_size(), B.get_size(), cutnet);    
                pass_start = false;
            }
            else{
                LocalMinDist.update(CELL_array, C, A.get_size(), B.get_size(), cutnet);
            }
            
            assert(CountCutNet(A, NET_array, N) == cutnet);


            if(!stuck)
                BaseCell = ChooseBaseCell_gain(A, B, r);
            else{
                //destroy_balance = !A.get_balance();
                //destroy_balance = true;
                BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
            }
        }


        LoadDistribution(LocalMinDist, CELL_array, C, cutnet);
        BlockReinitialization(A, B, CELL_array, NET_array, no_large_net);


        if(!stuck && stuck_check)
            stuck_count++;
        
        if(stuck)
            stuck_count--;

        if(stuck_option){
            if(stuck_count == stuck_criteria){
                stuck = true;
                stuck_count *= stuck_out_itr;
                destroy_balance = !A.get_balance();
                destroy_balance = true;
            }
        }

        if(stuck_count == 0)
            stuck = false;

    }


    printf("\n---------- After FM Algorithm----------\n");
    

    GlobalMinDist.writeCellDist(CELL_array, C);
    GlobalMinDist.printCellDist();

    //printf("Final min num of cutnet: %d\n", global_min_cutnet);

    //Check(A, B, NET_array, N);


    delete[] NET_array;
    delete[] CELL_array;

    end = clock();
    total_time = (double)end - start;

    printf("\ntime: %fs\n", total_time / CLOCKS_PER_SEC);


    return 0;
}