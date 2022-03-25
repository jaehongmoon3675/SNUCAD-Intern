#define NDEBUG
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <list>
#include <algorithm>
#include <ctime>
#include <vector>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "CellDist.h"
#include "WriteFile.h"

int ALPHA = 100;

//stuck check에 유의

int main(){
    int _N, N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    bool bin_area_print = true;

    int ll_x, ll_y, ur_x, ur_y;
    std::srand(std::time(nullptr));

    //세로가 n, 가로는 m
    const int map_n = 2;
    const int map_m = 2;
    printf("map_n: %d, map_m: %d\n", map_n, map_m);
    int overlap_x, overlap_y;

    std::vector<int> *BIN_array = new std::vector<int>[map_n * map_m];

    //alpha_tune이 커질수록 overlap늘고 cutnet 준다
    int alpha_tune = 100;
    const int block_num = 6;
    const int InitVer = 1;
    const int pass = 10;
    const double skew = 0.05;
    const bool rand_mode = true;
    int tight = 1; //tight이 커질수록 overlap에 대한 정밀도가 낮아짐...

    //printf("map_n: %d, map_m: %d, block_num: %d, pass: %d, skew: %.2f\n", map_n, map_m, block_num, pass, skew);
  
    const int file_num = 3;
    const std::string file_name_arr[8] = {"aes_128", "ldpc", "jpeg", "wb_dma", "ecg", "ac97", "nova", "tate_pairing"};
    std::string file_name = file_name_arr[file_num];
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    clock_t start, end;
    double total_time;

    start = clock();
    
    read_place(C, CELL_array, file_name, map_n, map_m, BIN_array, ll_x, ll_y, ur_x, ur_y);
    overlap_x = (ur_x - ll_x) / (2 * tight);
    overlap_y = ur_y - ll_y + 2; //overlap_y의 scale은 항상 1이 되도록 하자.
    P = read_hgr(_N, C, NET_array, CELL_array, file_name, overlap_x, overlap_y);
    read_hgr_map(C, CELL_array, file_name);
    W = read_hgr_area(C, CELL_array, file_name);
    //read_place(C, CELL_array, file_name, map_n, map_m, BIN_array);
    N = _N + overlap_x * overlap_y;
    read_place(C, CELL_array, _N, N, NET_array, file_name, map_n, map_m, BIN_array, ll_x, ll_y, ur_x, ur_y, overlap_x, overlap_y);
    read_partial_part(C, CELL_array, file_name);
    
    int overlap_count = 0;
    int overlap_max = 0;
    for(int i = _N + 1; i <= N; i++){
        overlap_count += NET_array[i].get_cell_count();
        if(overlap_max <= NET_array[i].get_cell_count())
            overlap_max = NET_array[i].get_cell_count();
    }

    //printf("overlap count: %d\n", overlap_count);

    for(int i = 1; i <= N; i++)
        NET_array[i].adjust_weight(overlap_max, alpha_tune, rand_mode);

    /*
    bool **map = new bool*[ur_x - ll_x + 1];
    
    for(int i = 0; i <= ur_x -ll_x; i++){
        map[i] = new bool[ur_y - ll_y + 1];

        for(int j = 0; j <= ur_y - ll_y + 1; j++)
            map[i][j] = false;
    }

    int overlap_max = 0;
    int overlap_temp = 0;

    for(int i = 1; i <= N; i++){
        overlap_temp = NET_array[i].adjust_overlap(map, ll_x, ll_y);

        if(overlap_temp > overlap_max)
            overlap_max = overlap_temp;
    }

    printf("max overlap: %d\n", overlap_max);

    CountOverlap(map, N, NET_array, ll_x, ll_y, ur_x, ur_y, block_num);

    for(int i = 1; i <= N; i++)
        NET_array[i].adjust_weight(overlap_max, 50);
    */

    bin_based_FM(InitVer, pass, CELL_array, NET_array, C, N, P, W, block_num, skew, map_n, map_m, BIN_array);
    //printf("bin_based finish\n");

    int cutnet, degree;
    //int degree = calculate_degree(CELL_array, C, NET_array, N, block_num, cutnet);
    //printf("before... \ncutnet: %d, degree: %d\n", cutnet, degree);

    //write_output(CELL_array, C, file_name, block_num, InitVer, pass);

    int CB = block_num * map_n * map_m;

    Cell* CELL_BIN_array = new Cell[CB + 1];
    Net* NET_BIN_array = new Net[N + 1];
    int* bin_area = new int[map_m * map_n];

    for(int i = 0; i < map_m * map_n; i++)
        bin_area[i] = 0;

    //printf("cell2 bin: %d\n", CELL_array[2].get_bin());

    for(int i = 1; i <= CB; i++){
        CELL_BIN_array[i].set_cell_num(i);
        CELL_BIN_array[i].set_current_block_num((i - 1) % block_num);
        CELL_BIN_array[i].set_bin((i - 1) / block_num);
    }

    int current_CB;

    for(int i = 1; i <= C; i++){
        current_CB = block_num * CELL_array[i].get_bin() + CELL_array[i].get_current_block_num() + 1;
        CELL_BIN_array[current_CB].Cell_set.push_back(i);

        if(!CELL_BIN_array[current_CB].get_fixed() && CELL_array[i].get_fixed())
            CELL_BIN_array[current_CB].set_fixed(true);

        bin_area[CELL_array[i].get_bin()] += CELL_array[i].get_size();

        //printf("%d, %d, %d\n", i, CELL_BIN_array[current_CB].get_current_block_num(), CELL_array[i].get_current_block_num());
        assert(CELL_BIN_array[current_CB].get_current_block_num() == CELL_array[i].get_current_block_num());
        assert(current_CB <= CB);
    }

    bool *NET_BIN_check = new bool[1 + CB];
    int net_bin_count = 0;
    int NB = 0;
    int exists_count = 0;

    for(int i = 1; i <= N; i++){
        net_bin_count = 0;
        
        
        for(int j = 1; j <= CB; j++)
            NET_BIN_check[j] = true;

        for(auto itr = NET_array[i].cell_list.begin(); itr != NET_array[i].cell_list.end(); itr++){
            current_CB = block_num * (*itr)->get_bin() + (*itr)->get_current_block_num() + 1;

            if(NET_BIN_check[current_CB]){
                net_bin_count++;
                NET_BIN_check[current_CB] = false;
            }
        }

        if(net_bin_count <= 1){
            exists_count++;
            continue;
        }
        
        NB++;
        NET_BIN_array[NB].set_net_num(NB);
        NET_BIN_array[NB].set_weight(NET_array[i].get_weight());

        for(int j = 1; j <= CB; j++)
            NET_BIN_check[j] = true;

        for(auto itr = NET_array[i].cell_list.begin(); itr != NET_array[i].cell_list.end(); itr++){
            current_CB = block_num * (*itr)->get_bin() + (*itr)->get_current_block_num() + 1;

            if(NET_BIN_check[current_CB]){
                NET_BIN_array[NB].push_cell(CELL_BIN_array + current_CB);
                CELL_BIN_array[current_CB].push_net(NET_BIN_array + NB);
                NET_BIN_check[current_CB] = false;
            }
        }
    }
    //printf("exist!: %d\n", exists_count);
    //printf("%d %d\n", N, NB);

    /*
    for(int i = 1; i <= NB; i++){
        if(NET_BIN_array[i].get_cell_count() > CB)
            printf("too many cell!\n");
    }
    */

    /*
    for(int i = 1; i <= CB; i++)
        assert(CELL_BIN_array[i].get_size() == 1);
    */
    //W = CB, P는 유지
    //printf("----------------------------\n");
    //printf("bin FM starts!\n");

    int *bin_block_area = new int[CB + 1];
    int *block_area = new int[block_num];

    for(int i = 1; i <= CB; i++)
        bin_block_area[i] = 0;

    for(int i = 0; i < block_num; i++)
        block_area[i] = 0;
    
    /*
    for(int i = 0; i < map_m * map_n; i++)
        printf("%d ", bin_area[i]);
    */

    for(int i = 1; i <= C; i++){
        current_CB = block_num * CELL_array[i].get_bin() + CELL_array[i].get_current_block_num() + 1;
        bin_block_area[current_CB] += CELL_array[i].get_size();
    }

    if(bin_area_print){
        printf("\t");
        for(int i = 0; i < block_num; i++)
            printf("Block %d\t\t", i);
        printf("skew\n");
    }

    int current_bin = 0;
    int min = W, max = 0;
    for(int i = 1; i <= CB; i++){
        if(i % block_num == 1){
            if(bin_area_print)
                printf("bin %d\t", i / block_num);
            current_bin = i / block_num;
            min = bin_area[current_bin]; max = 0;
        }

        if(min > bin_block_area[i])
            min = bin_block_area[i];
        
        if(max < bin_block_area[i])
            max = bin_block_area[i];

        if(bin_area_print)    
            printf("%d / %.2f\t", bin_block_area[i], (double)bin_block_area[i] / bin_area[current_bin]);

        block_area[(i - 1) % block_num] += bin_block_area[i];

        
        if(bin_area_print && i % block_num == 0){
            printf("%.2f", (double)(max - min) / bin_area[current_bin]);
            printf("\n");
        }
    }

    FM(InitVer, pass, CELL_BIN_array, NET_BIN_array, CB, NB, P, CB, block_num, nullptr, 0, 0, true, 0);

    int current_block_num;
    for(int i = 1; i <= CB; i++){
        current_block_num = CELL_BIN_array[i].get_current_block_num();

        for(int j = 0; j < CELL_BIN_array[i].Cell_set.size(); j++){
            CELL_array[CELL_BIN_array[i].Cell_set[j]].set_current_block_num(current_block_num);
            //printf("%d %d %d %d\n", CELL_array[j].get_bin(), (i - 1) / block_num, i, j);
            assert(CELL_array[CELL_BIN_array[i].Cell_set[j]].get_bin() == (i - 1) / block_num);
        }
    }

    for(int i = 1; i <= CB; i++)
        bin_block_area[i] = 0;

    for(int i = 0; i < block_num; i++)
        block_area[i] = 0;
    
    /*
    for(int i = 0; i < map_m * map_n; i++)
        printf("%d ", bin_area[i]);
    */

    for(int i = 1; i <= C; i++){
        current_CB = block_num * CELL_array[i].get_bin() + CELL_array[i].get_current_block_num() + 1;
        bin_block_area[current_CB] += CELL_array[i].get_size();
    }

    if(bin_area_print){
        printf("\t");
        for(int i = 0; i < block_num; i++)
            printf("Block %d\t\t", i);
        printf("skew\n");
    }

    current_bin = 0;
    min = W, max = 0;
    for(int i = 1; i <= CB; i++){
        if(i % block_num == 1){
            if(bin_area_print)
                printf("bin %d\t", i / block_num);
            current_bin = i / block_num;
            min = bin_area[current_bin]; max = 0;
        }

        if(min > bin_block_area[i])
            min = bin_block_area[i];
        
        if(max < bin_block_area[i])
            max = bin_block_area[i];

        if(bin_area_print)    
            printf("%d / %.2f\t", bin_block_area[i], (double)bin_block_area[i] / bin_area[current_bin]);

        block_area[(i - 1) % block_num] += bin_block_area[i];

        
        if(bin_area_print && i % block_num == 0){
            printf("%.2f", (double)(max - min) / bin_area[current_bin]);
            printf("\n");
        }
    }
    
    if(bin_area_print){
        min = W; max = 0;
        printf("total\t");
        for(int i = 0; i < block_num; i++){
            if(min > block_area[i])
                min = block_area[i];
            
            if(max < block_area[i])
                max = block_area[i];
            
            printf("%d / %.2f\t", block_area[i], (double)block_area[i] / W);
        }
        printf("%.2f\n", (double)(max - min) / W);
    }


    //int cutnet = FM(InitVer, pass, CELL_array, NET_array, C, N, P, W, block_num, nullptr, skew, 0, true);
    cutnet;
    degree = calculate_degree(CELL_array, C, NET_array, _N, block_num, cutnet);
    //printf("\ncutnet: %d, degree: %d\n", cutnet, degree);

    write_output(CELL_array, C, file_name, map_n, map_m, block_num, InitVer, pass);
    check_partial_part(C, CELL_array, file_name);
    overlap_count = CalculateTotalOverlap(NET_array + _N, overlap_x * overlap_y, block_num);
    //check_place(C, CELL_array, file_name, map_n,  map_m);

    printf("cutnet: %d, degree: %d, overlap: %d\n", cutnet, degree, overlap_count);
    //CountOverlap(map, N, NET_array, ll_x, ll_y, ur_x, ur_y, block_num);

    delete[] NET_array;
    delete[] CELL_array;
    delete[] BIN_array;
    delete[] CELL_BIN_array;
    delete[] NET_BIN_array;
    delete[] NET_BIN_check;
    delete[] bin_area;
    delete[] bin_block_area;
    delete[] block_area;

    /*
    for(int i = 0; i <= ur_x -ll_x; i++){
        delete map[i];
    }

    delete[] map;
    */
    end = clock();
    total_time = (double)end - start;


    printf("runtime: %.3fs\n", total_time / CLOCKS_PER_SEC);
    printf("\n");

    system("pause");
    
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