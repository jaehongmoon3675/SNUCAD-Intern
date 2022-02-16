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
extern double get_max_gain_cell_time;

//stuck check에 유의
void FM_pass(int C, int N, double r, Cell* CELL_array, Net* NET_array, Block &A, Block &B, const bool stuck, CellDist& LocalMinDist){
    Cell* BaseCell = nullptr;
    bool review = false;
    int init_pass_start = C / 100;
    int move_count = 0, move_count_at_min = 1;
    int A_size_at_min = A.get_size();
    int min_cutnet, cutnet;

    min_cutnet = LocalMinDist.get_cutnet();

    //Block TempA = A;
    //Block TempB = B;

    if(stuck)
        printf("stuck!\n");
    else
        printf("not stuck!\n");

    const bool destroy_balance = A.get_balance();
    //bool destroy_balance = true;
    
    if(destroy_balance)
        printf("Destroy...!\n");
    else
        printf("balance\n");      
    
    clock_t choosemove_time_temp, reinit_time_temp, choose_time_temp, move_time_temp;
    double choose_time = 0, move_time = 0;
    double get_max_gain_cell_time_temp;

    do{ 
        int pass_start = init_pass_start; //계속해서 초기화 시켜여줘야함에 유의!
        cutnet = LocalMinDist.get_cutnet();
        //printf("----------------------------");

        //assert(TempA == A);
        //assert(TempB == B);
        choose_time = 0, move_time = 0;

        choosemove_time_temp = clock();

        get_max_gain_cell_time_temp = get_max_gain_cell_time;
        
        choose_time_temp = clock();
        if(!stuck)
            BaseCell = ChooseBaseCell_gain(A, B, r);
            //BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
        else{
            BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
        }
        choose_time += (clock() - (double)choose_time_temp);

        while(BaseCell != nullptr){        
            //A.print_Block(CELL_array);
            //B.print_Block(CELL_array);
            move_time_temp = clock();
            if(BaseCell->get_current_block() == &A){
                cutnet -= A.gain[BaseCell->get_cell_num()];
                MoveCell(A, B, BaseCell);
            }
            else{
                cutnet -= B.gain[BaseCell->get_cell_num()];
                MoveCell(B, A, BaseCell);
            }
            move_time += ((clock() - (double)move_time_temp));

            if(review)
                move_count_at_min--;
            else
                move_count++;

            if(move_count_at_min == 0)
                break;       
            
            if(cutnet < 0){
                printf("negative cutnet error! move_count: %d\n", move_count);
                assert(cutnet >= 0);
            }

            if(pass_start > 0){    
                if(!review){              
                    min_cutnet = cutnet;
                    move_count_at_min = move_count;
                    //printf("%d ", move_count_at_min);
                    A_size_at_min = A.get_size();
                }
                pass_start--;
            }
            else if(!review){
                if(min_cutnet > cutnet){
                    min_cutnet = cutnet;
                    move_count_at_min = move_count;
                    A_size_at_min = A.get_size();
                }
                else if(min_cutnet == cutnet){
                    if(std::abs(A.get_size() - A.get_W()) < std::abs(A_size_at_min - A.get_W())){
                        min_cutnet = cutnet;
                        move_count_at_min = move_count;
                        A_size_at_min = A.get_size();
                    }
                }
            }

            choose_time_temp = clock();
            if(!stuck)
                BaseCell = ChooseBaseCell_gain(A, B, r);
                //BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
            else{
                BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
            }
            choose_time += (clock() - (double)choose_time_temp);
            /*
            if(!stuck)
                BaseCell = ChooseBaseCell_gain(A, B, r);
            else{
                //destroy_balance = !A.get_balance();
                //destroy_balance = true;
                BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
            }
            */

           //if(stuck && move_count >= init_pass_start / 2 && A.get_strict_balance() >= 0 && A.get_strict_balance() != destroy_balance){
            if(stuck && move_count >= init_pass_start / 5 && A.get_balance() != destroy_balance){
                min_cutnet = cutnet;
                move_count_at_min = move_count;
                break;
           }

            
           if(move_count > C - C / 4){
               move_count = 0;
               break;
           }
           
        }

        printf("choose + move time: %fs, move_count_at_min: %d\n", (clock() - (double)choosemove_time_temp) / CLOCKS_PER_SEC, move_count_at_min);
        //printf("choose time: %fs, move time: %fs\n", choose_time/ CLOCKS_PER_SEC, move_time / CLOCKS_PER_SEC);
        //printf("get_max_gain_cell_time: %fs\n", (get_max_gain_cell_time - get_max_gain_cell_time_temp)/ CLOCKS_PER_SEC);
        //printf("review finish\n");

        reinit_time_temp = clock();
        if(!review)
            LoadDistribution(LocalMinDist, A, B, CELL_array, C);
        
        BlockReinitialization(C, A, B, CELL_array, NET_array);
        //printf("reinit time: %fs\n\n", (clock() - (double)reinit_time_temp) / CLOCKS_PER_SEC);


        review = !review;
    }while(review);

    //printf("min_cutnet: %d, cutnet: %d\n", min_cutnet, cutnet);
    assert(min_cutnet == cutnet);

    //printf("\n final \n");
    //A.print_Block(CELL_array);
    //B.print_Block(CELL_array);
    

    LocalMinDist.overWrite(CELL_array, C, A.get_size(), B.get_size(), min_cutnet);
}

int main(int argc, char ** argv){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    const double r = 0.5; //balance factor
    const double pm_r = 0.05;
    const int pass = 1000; //how many pass we go through
    const int k = 5;
    const int InitVer = 2;
    int global_min_cutnet, local_min_cutnet;
    int cutnet;
    const int stuck_criteria = 5;
    const double stuck_out_itr = 1;
    const bool balance_option = false; //true면 smax 기반, false면 비율 기반
    const bool time_option = true;
    const bool stuck_option = true;
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
    case 3:
        BlockInitialization(A, B, CELL_array, NET_array, C, N, 0);
        break;
    default:
        BlockInitialization(A, B, CELL_array, NET_array, C, N, 0);
    }

    BlockReinitialization(C, A, B, CELL_array, NET_array);

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
    clock_t record_time_temp, pass_time_temp;

    for(int i = 0; i < pass; i++){
        //printf("\npass %d\n", i);
        pass_time_temp = clock();
        //LocalMinDist.printCellDist();

        cutnet = LocalMinDist.get_cutnet();

        assert(A.get_size() == LocalMinDist.get_A_size());
        assert(B.get_size() == LocalMinDist.get_B_size());

        FM_pass(C, N, r, CELL_array, NET_array, A, B, stuck, LocalMinDist);

        if(stuck_option){
            if(!stuck && (GlobalMinDist.get_cutnet() <= LocalMinDist.get_cutnet()))
                stuck_count++;
        
            if(stuck)
                stuck_count--;

            if(stuck_count == stuck_criteria){
                stuck = true;
                stuck_count *= stuck_out_itr;
                destroy_balance = !A.get_balance();
                destroy_balance = true;
            }

            if(stuck_count == 0)
                stuck = false;
            
            printf("\nstuck_count: %d\n", stuck_count);
        }

        
        if(GlobalMinDist.get_cutnet() > LocalMinDist.get_cutnet()){
            GlobalMinDist = LocalMinDist;
            stuck_count = 0;
            stuck = false;
        }
        //assert(CountCutNet(A, NET_array, N) == cutnet);

        printf("pass %d end, cutnet: %d, running time: %fs\n", i, LocalMinDist.get_cutnet(), (clock() - (double)pass_time_temp) / CLOCKS_PER_SEC);
        printf("A size: %d, B size: %d\n", LocalMinDist.get_A_size(), LocalMinDist.get_B_size());
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