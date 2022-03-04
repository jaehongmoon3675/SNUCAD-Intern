#define NDEBUG

#include <iostream>
#include <cassert>
#include <list>
#include <algorithm>
#include <ctime>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "CellDist.h"
#include "WriteFile.h"

//stuck check에 유의

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight

    const int block_num = 3;
    const int InitVer = 1;
    const int pass = 100;
    const double skew = 0.1;
  
    const int file_num = 3;
    const std::string file_name_arr[4] = {"aes_128", "ldpc", "jpeg", "initPlace"};
    std::string file_name = file_name_arr[file_num];
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    clock_t start, end;
    double total_time;

    start = clock();

    P = read_hgr(N, C, NET_array, CELL_array, file_name);
    read_hgr_map(C, CELL_array, file_name);
    W = read_hgr_area(C, CELL_array, file_name); 
    
    int cutnet = FM(InitVer, pass, CELL_array, NET_array, C, N, P, W, block_num, nullptr, skew, 0, true);
    int degree = calculate_degree(CELL_array, C, NET_array, N, block_num, cutnet);

    write_output(CELL_array, C, file_name, block_num, InitVer, pass);

    printf("cutnet: %d, degree: %d\n", cutnet, degree);

    delete[] NET_array;
    delete[] CELL_array;

    end = clock();
    total_time = (double)end - start;

    printf("\ntime: %fs\n", total_time / CLOCKS_PER_SEC);

    return 0;
}

/*
int main(int argc, char ** argv){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    const double r = 0.5; //balance factor
    const double pm_r = 0.05;
    const int pass = 200; //how many pass we go through
    const int k = 5;
    const int InitVer = 2;
    const int file_num = 1;
    int global_min_cutnet, local_min_cutnet;
    int cutnet;
    //const int stuck_criteria = 8;
    //const double stuck_out_itr = 0.5;
    const int stuck_criteria = 10;
    const int stuck_out_itr = 5;
    const bool balance_option = false; //true면 smax 기반, false면 비율 기반
    const bool time_option = false;
    const bool stuck_option = true;
    const std::string file_name_arr[4] = {"aes_128", "ldpc", "jpeg", "initPlace"};
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    int min_pass = 0;

    std::string file_name = file_name_arr[file_num];

    clock_t start, end, move_time, reinit_time, Move_time;
    double total_time;

    start = clock();

    P = read_hgr(N, C, NET_array, CELL_array, file_name);
    read_hgr_map(C, CELL_array, file_name);
    W = read_hgr_area(C, CELL_array, file_name); 
    
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
        printf("balance ratio +- %f, num of pass: %d\n", pm_r, pass);
        balance_low_bound = (r - pm_r)*W;
        balance_up_bound = (r + pm_r)*W;
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
    default:
        read_output_part(A, B, C, CELL_array);
    }

    BlockReinitialization(C, A, B, CELL_array, NET_array, 0);

    global_min_cutnet = CountCutNet(A, NET_array, N);
    local_min_cutnet = global_min_cutnet;
    cutnet = global_min_cutnet;

    CellDist GlobalMinDist(C, N, W*r, cutnet, A.get_size(), B.get_size(), &A, &B, CELL_array);
    CellDist LocalMinDist(C, N, W*r, cutnet, A.get_size(), B.get_size(), &A, &B, CELL_array);

    printf("initial cutnet num: %d\n", global_min_cutnet);

    Cell* BaseCell = nullptr;
    int temp;
    int move_count = 0;
    int stuck_count = 0;
    int stuck_out_count = 0;
    bool big_wave = false;
    const int big_wave_criteria = 8;
    int gain_update_count_temp = 0;
    bool stuck_check = true; //얼마나 오래동안 stuck 되었는지를 체크하기 위함
    bool stuck = false; //실제 stuck 되어있는지 여부
    bool destroy_balance = false;
    int min_cutnet = GlobalMinDist.get_cutnet();

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
    clock_t record_time_temp, pass_time_temp;

    int pass_num = 0;

    for(int i = 0; i < pass; i++){
        //printf("\npass %d\n", i);
        pass_time_temp = clock();
        //LocalMinDist.printCellDist();

        cutnet = LocalMinDist.get_cutnet();

        assert(A.get_size() == LocalMinDist.get_A_size());
        assert(B.get_size() == LocalMinDist.get_B_size());

        FM_pass(C, N, r, i, CELL_array, NET_array, A, B, stuck, LocalMinDist, big_wave);

        if(stuck_option){
            if(!stuck && (min_cutnet <= LocalMinDist.get_cutnet()))
                stuck_count++;

            if(!stuck && min_cutnet > LocalMinDist.get_cutnet()){
                min_cutnet = LocalMinDist.get_cutnet();
                stuck_count = 0;
            }

            if(stuck)
                stuck_count--;

            if(!stuck && stuck_count == stuck_criteria){
                stuck = true;
                stuck_count = stuck_out_itr;
                destroy_balance = !A.get_balance();
                destroy_balance = true;
            }

            if(stuck_count == 0){
                min_cutnet = LocalMinDist.get_cutnet();
                stuck = false;
                stuck_out_count++;
            }

            if(big_wave)
                stuck_out_count++;

            if((stuck_out_count + 1) % big_wave_criteria == 0)
                big_wave = true;
            else
                big_wave = false;

            //printf("\nstuck_count: %d\n", stuck_count);
        }

        
        //if(stuck && 0 <= (i + 1) % 128 && (i + 1) % 128 < stuck_criteria){
        //    printf("out from stuck!\n");
        //    StuckOut(A, B, CELL_array, NET_array, C, N, !A.bigger());
        //    BlockReinitialization(C, A, B, CELL_array, NET_array, i);
        //    LocalMinDist.update(CELL_array, c, A.get_size(), B.get_size(), CountCutNet(A, NET_array, N));
        //    stuck_count = 0;
        //    stuck = false;
        //}
        
                
        if(GlobalMinDist.get_cutnet() > LocalMinDist.get_cutnet()){
            GlobalMinDist = LocalMinDist;
            min_cutnet = GlobalMinDist.get_cutnet();
            stuck_count = 0;
            min_pass = i;
            stuck = false;
            stuck_out_count++;
        }

        if(time_option){
            printf("pass %d end, cutnet: %d, running time: %fs\n", i, LocalMinDist.get_cutnet(), (clock() - (double)pass_time_temp) / CLOCKS_PER_SEC);
            printf("A size: %d, B size: %d\n", LocalMinDist.get_A_size(), LocalMinDist.get_B_size());
        }
        printf("pass %d end, cutnet: %d\n", i, LocalMinDist.get_cutnet());
        //printf("A size: %d, B size: %d\n", LocalMinDist.get_A_size(), LocalMinDist.get_B_size());
    }


    printf("\n---------- After FM Algorithm----------\n");
    
    GlobalMinDist.writeCellDist(CELL_array, C, file_name, InitVer, pass);
    printf("Move: %d ,", min_pass);
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
*/