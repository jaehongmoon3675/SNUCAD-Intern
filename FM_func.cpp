#define NDEBUG
#include <iostream>
#include <cassert>
#include <list>
#include <algorithm>
#include <ctime>
#include <string>
#include <vector>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "Block.h"
#include "CellDist.h"
#include "WriteFile.h"
#include "FM_func.h"

extern int ALPHA;

void get_max(const int C, const Cell* CELL_array, int &pmax, int &smax){
    pmax = -1;
    smax = -1;

    pmax = CELL_array[1].get_pin();
    smax = CELL_array[1].get_size();

    for(int i = 2; i <= C; i++){    
        if(pmax < CELL_array[i].get_pin())
            pmax = CELL_array[i].get_pin();

        if(smax < CELL_array[i].get_size())
            smax = CELL_array[i].get_size();
    }
}

int CountCutNet(Block &block, Net *NET_array, int N){
    int cutnet = 0;
    
    for(int i = 1; i <= N; i++){
        if((block.ith_net_distribution(i) < NET_array[i].cell_list.size()) && (block.ith_net_distribution(i) != 0))
            cutnet += NET_array[i].get_weight();
        else if(block.ith_net_distribution(i) > NET_array[i].cell_list.size())
            printf("CountcutNet error\n");
    }

    return cutnet;
}

void CountCutNetAgain(Block &block, Net *NET_array, int N){
    
    for(int i = 1; i <= N; i++){
        if(block.ith_net_distribution(i) < NET_array[i].cell_list.size())
            printf("Net %d\n", i);
        else if(block.ith_net_distribution(i) > NET_array[i].cell_list.size())
            printf("CountCutNetAgain error\n");
    }
}

void Check(Block &A, Block &B, Net* NET_array, int N){
    for(int i = 1; i <= N; i++)
        if(A.ith_net_distribution(i) + B.ith_net_distribution(i) != NET_array[i].cell_list.size())
            printf("%d th net error\n", i);
}

int get_max_cell_count(Net* NET_array, int N){
    int max_cell_count = 0;
    for(int i = 1; i <= N; i++)
        if(max_cell_count < NET_array[i].get_cell_count())
            max_cell_count = NET_array[i].get_cell_count();
    
    return max_cell_count;
}

void FM_pass(int C, int N, double r, int pass_num, Cell* CELL_array, Net* NET_array, Block &A, Block &B, const bool stuck, CellDist& LocalMinDist, bool big_wave, bool alternate){
    Cell* BaseCell = nullptr;
    bool review = false;
    bool how_to_start = true;
    //init_pass_start의 값이 크면 결과는 빨리 나오나, 운적인 요소가 크게 작용
    int init_pass_start = (100 < (C / 50))? 100 : (C / 50) + 2;
    //printf("init_pass: %d\n", init_pass_start);
    int move_count = 0, move_count_at_min = 2;
    int A_size_at_min = A.get_size();
    const int out_from_stuck_citeria = 1;
    int out_from_stuck = out_from_stuck_citeria;
    int min_cutnet, cutnet;

    BlockReinitialization(C, A, B, CELL_array, NET_array, pass_num);

    min_cutnet = LocalMinDist.get_cutnet();
    //Block TempA = A;
    //Block TempB = B;

    /*
    if(stuck)
        printf("stuck!\n");
    else
        printf("not stuck!\n");
    */

    const bool destroy_balance = A.get_balance();
    const bool bigger = A.bigger();
    bool stuck_temp = stuck;
    Block* alternate_block;
    //bool destroy_balance = true;
    
    /*
    if(destroy_balance)
        printf("Destroy...!\n");
    else
        printf("balance\n");      
    */
    double choose_time = 0, move_time = 0;
    double get_max_gain_cell_time_temp;
    int pass_start;

    do{ 
        //int pass_start = init_pass_start; //계속해서 초기화 시켜여줘야함에 유의!
        pass_start = (how_to_start)? init_pass_start : 2;
        cutnet = LocalMinDist.get_cutnet();
        alternate_block = bigger? &A : &B;
        //printf("----------------------------");

        //assert(TempA == A);
        //assert(TempB == B);
        choose_time = 0, move_time = 0;
        if(alternate){
            BaseCell = alternate_block->get_max_gain_cell();

            if(BaseCell != nullptr){
                alternate_block->remove_from_BUCKET(BaseCell);

            }
            /*
            else{
                printf("Something is wrong...\n");
                printf("A: %d, %f\n", A.get_size(), A.get_ubound());
                assert(BaseCell != nullptr);
            }
            */
            alternate_block = (alternate_block == &A)? &B : &A;
        }
        else if(!stuck_temp)
            BaseCell = ChooseBaseCell_gain(A, B, r);
            //BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
        else{
            BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance, bigger);
        }

        while(BaseCell != nullptr){     
            //A.print_Block(CELL_array);
            //B.print_Block(CELL_array);
            if(BaseCell->get_current_block() == &A){
                /*
                if(!A.check_block() || !B.check_block()){
                    printf("hello 3\n");
                    printf("hdlhldddddddh\n");
                }
                */
               //printf("A gain %d\n", A.gain[BaseCell->get_cell_num()]);
                cutnet -= A.gain[BaseCell->get_cell_num()];
                if(!how_to_start && A.gain[BaseCell->get_cell_num()] > 0)
                    pass_start = false;
                MoveCell(A, B, BaseCell);

                /*
                if(!A.check_block() || !B.check_block()){
                    printf("hello 3\n");
                    printf("hdlhlddddh\n");
                }
                */
            }
            else{
                //printf("B gain %d\n", B.gain[BaseCell->get_cell_num()]);
                cutnet -= B.gain[BaseCell->get_cell_num()];
                if(!how_to_start && B.gain[BaseCell->get_cell_num()] > 0)
                    pass_start = false;
                MoveCell(B, A, BaseCell);
            }

            if(review)
                move_count_at_min--;
            else
                move_count++;

            if(move_count_at_min == 0)
                break;       
            
            if(cutnet < 0){
                //printf("negative cutnet error! move_count: %d\n", move_count);
                //assert(cutnet >= 0);
            }

            if(pass_start){    
                if(!review){              
                    min_cutnet = cutnet;
                    move_count_at_min = move_count;
                    //printf("%d ", move_count_at_min);
                    A_size_at_min = A.get_size();

                    if(how_to_start)
                        pass_start--;
                }
            }
            else if(!review){
                if(!alternate){
                    if(min_cutnet > cutnet){
                        min_cutnet = cutnet;
                        move_count_at_min = move_count;
                        A_size_at_min = A.get_size();
                    }
                }
                else{
                    if(min_cutnet > cutnet && move_count %2 == 0){
                        min_cutnet = cutnet;
                        move_count_at_min = move_count;
                        A_size_at_min = A.get_size();
                    }
                }
            }


            if(alternate){
                BaseCell = alternate_block->get_max_gain_cell();

                if(BaseCell != nullptr){
                    alternate_block->remove_from_BUCKET(BaseCell);
                }
                
                alternate_block = (alternate_block == &A)? &B : &A;
            }
            else if(!stuck_temp)
                BaseCell = ChooseBaseCell_gain(A, B, r);
                //BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
            else{
                BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance, bigger);
            }
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
            if(stuck_temp && A.bigger() != bigger && A.get_balance() != destroy_balance && (!alternate || move_count %2 == 0)){
                if(!review){
                    min_cutnet = cutnet;
                    move_count_at_min = move_count;
                }
                stuck_temp = false;
            }

            if(big_wave && stuck && !stuck_temp && (A.bigger() == bigger) && (!alternate || move_count %2 == 0)){
                stuck_temp = true;

                if(move_count > C / 4){
                    move_count = 0;
                    break;
                }
            }
            
        
            if(move_count > C - C / 8 && move_count % 2 == 0){
               move_count = 0;
               break;
            }
           
        }

        //printf("move_count_at_min: %d\n", move_count_at_min);
        assert(!alternate || move_count_at_min % 2 == 0);

        //printf("choose + move time: %fs, move_count_at_min: %d\n", (clock() - (double)choosemove_time_temp) / CLOCKS_PER_SEC, move_count_at_min);
        //printf("choose time: %fs, move time: %fs\n", choose_time/ CLOCKS_PER_SEC, move_time / CLOCKS_PER_SEC);
        //printf("get_max_gain_cell_time: %fs\n", (get_max_gain_cell_time - get_max_gain_cell_time_temp)/ CLOCKS_PER_SEC);
        //printf("review finish\n");

        if(!review){
            LoadDistribution(LocalMinDist, A, B, CELL_array, C);
            if(stuck != stuck_temp)
                stuck_temp = stuck;
            
            BlockReinitialization(C, A, B, CELL_array, NET_array, pass_num);
        }
        //printf("reinit time: %fs\n\n", (clock() - (double)reinit_time_temp) / CLOCKS_PER_SEC);

        out_from_stuck = out_from_stuck_citeria;
        stuck_temp = stuck;
        review = !review;
    }while(review);

    //printf("min_cutnet: %d, cutnet: %d\n", min_cutnet, cutnet);
    assert(min_cutnet == cutnet);

    //printf("\n final \n");
    //A.print_Block(CELL_array);
    //B.print_Block(CELL_array);
    

    LocalMinDist.overWrite(CELL_array, C, NET_array, N, min_cutnet);
}


void read_past_block_record(const int _N, const int _C, const Net* _NET_array, const Cell* _CELL_array, int N, int C, Net* &NET_array, Cell* &CELL_array, const Block* current_block, int &pin_num, int &total_weight, int * current_to_past){
    pin_num = 0;
    total_weight = 0;

    char c;
    int temp_cell_num;

    NET_array = new Net[N + 1];
    CELL_array = new Cell[C + 1];

    int* past_to_current = new int[_C + 1];

    int j = 1;
    for(int i = 1; i <= _C; i++){
        if(_CELL_array[i].get_current_block() != current_block)
            continue;
        
        CELL_array[j].set_cell_num(j);
        CELL_array[j].set_size(_CELL_array[i].get_size());
        total_weight += CELL_array[j].get_size();
        CELL_array[j].set_name(_CELL_array[i].get_cell_name());
        CELL_array[j].set_bin(_CELL_array[i].get_bin());
        CELL_array[j].set_fixed(_CELL_array[i].get_fixed());

        past_to_current[i] = j;
        current_to_past[j] = i;

        j++;
    }

    j = 1;

    for(int i = 1; i <= _N; i++){
        if(_NET_array[i].get_current_block() != current_block)
            continue;

        assert(j <= N);

        for(auto itr = _NET_array[i].cell_list.begin(); itr != _NET_array[i].cell_list.end(); itr++){
            pin_num++;

            assert((*itr)->get_current_block() == current_block);

            temp_cell_num = past_to_current[(*itr)->get_cell_num()];

            assert(temp_cell_num <= C);

            NET_array[j].set_net_num(j);
            NET_array[j].set_weight(_NET_array[i].get_weight());
            NET_array[j].push_cell(CELL_array + temp_cell_num);
            CELL_array[temp_cell_num].push_net(NET_array + j);
        }

        j++;
    }

    delete[] past_to_current;
}

//최종적으로 globalmin을 load하고 reinit 후, 반드시 Block의 set_current_block_of_net을 호출시켜주어야
//alternate시 초기화가 이미 되어있어야.
int FM(const int InitVer, const int pass, Cell* _CELL_array, Net* _NET_array, const int _C, const int _N, const int _P, const int _W, const int block_num, const Block* current_block, double skew, int bias, bool alternate, int bin_num){
    int C, N;
    int P, W; //P: total pin num, W: total weight

    const double r = (double)(block_num / 2) / block_num; //balance factor
    const double pm_r = skew / 2;
    //const int pass = 30; //how many pass we go through
    //const int InitVer = 1;
    int global_min_cutnet, local_min_cutnet;
    int min_cutnet, cutnet, total_cutnet;
    const int stuck_criteria = 10;
    const int stuck_out_itr = 5;
    const bool stuck_option = true;
    
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;
    int* current_to_past;

    if(current_block == nullptr){
        CELL_array = _CELL_array;
        NET_array = _NET_array;
        C = _C;
        N = _N;
        P = _P;
        W = _W;
        current_to_past = new int[C + 1];

        for(int i = 1; i <= C; i++)
            current_to_past[i] = i;
    }
    else{
        C = current_block->get_cell_num(_CELL_array, _C);
        N = current_block->get_uncut_count(_NET_array, _N);
        current_to_past = new int[C + 1];
        read_past_block_record(_N, _C, _NET_array, _CELL_array, N, C, NET_array, CELL_array, current_block, P, W, current_to_past);
    }

    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);

    //printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    //printf("pmax: %d, smax: %d, max_cell_count: %d\n", pmax, smax, max_cell_count);
    //printf("balance ratio: %f, num of pass: %d, k: %d\n", r, pass, k);


    //printf("balance ratio +- %f, num of pass: %d\n", pm_r, pass);
    if(alternate){
        balance_low_bound = (C / block_num) * (block_num / 2) - 1 - 0.5;
        balance_up_bound = (C / block_num) * (block_num / 2) + 1 + 0.5;
    }
    else{
        if(pm_r * W >= smax){
            balance_low_bound = (r - pm_r)*W;
            balance_up_bound = (r + pm_r)*W;
        }
        else{
            balance_low_bound = r*W - smax;
            balance_up_bound = r*W + smax;
        }
    }
    //printf("balance low bound: %f, balance up bound: %f\n", balance_low_bound, balance_up_bound);
    
    
    //Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, bias, bias + block_num/2 - 1, "A");
    //Block B(pmax, W - balance_up_bound, W - balance_low_bound, C, N, W, 1 - r, bias + block_num/2, bias + block_num, "B");
    
    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, bias + (block_num + 1) / 2, bias + block_num, "A");
    Block B(pmax, W - balance_up_bound, W - balance_low_bound, C, N, W, 1 - r, bias, bias + (block_num + 1) / 2 - 1, "B");


    //alternate인 경우 초기화 조건이 주어져야 한다.
    if(alternate){
        BlockInitialization_cell_bin(A, B, CELL_array, C, block_num);
    }
    else{
    switch(InitVer){
        case 1:
            //BlockInitialization(A, B, CELL_array, C);
            BlockInitialization(B, A, CELL_array, C, bin_num);
            break;
        case 2:
            BlockInitialization(B, A, CELL_array, NET_array, C, N);
            break;
        default:
            read_output_part(A, B, C, CELL_array);
        }   
    }

    BlockReinitialization(C, A, B, CELL_array, NET_array, 0);

    /*
    printf("cell 8 info\n");
    std::cout << CELL_array[8].get_current_block()->get_block_name() <<std::endl;
    printf("gain: %d\n", CELL_array[8].get_current_block()->gain[8]);
    */
    min_cutnet = CountCutNet(A, NET_array, N);
    cutnet = min_cutnet;

    //printf("before fm cutnet: %d\n", cutnet);

    CellDist GlobalMinDist(C, N, W*r, cutnet, &A, &B, CELL_array, NET_array);
    CellDist LocalMinDist(C, N, W*r, cutnet, &A, &B, CELL_array, NET_array);

    //printf("initial cutnet num: %d\n", global_min_cutnet);

    Cell* BaseCell = nullptr;
    int move_count = 0;
    int stuck_count = 0;
    int stuck_out_count = 0;
    bool big_wave = false;
    const int big_wave_criteria = 8;
    bool stuck_check = true; //얼마나 오래동안 stuck 되었는지를 체크하기 위함
    bool stuck = false; //실제 stuck 되어있는지 여부
    bool destroy_balance = false;

    bool pass_start = true;

    int pass_num = 0;

    for(int i = 0; i < pass; i++){
        //printf("pass %d\n", i);

        /*
        if(i == 0 || i == pass - 1)
            printf("A size: %d, B size: %d, A ubound: %f, B ubound: %f, A cell count: %d, B cell count: %d\n", A.get_size(), B.get_size(), A.get_ubound(), B.get_ubound(), A.get_cell_num(CELL_array, C), B.get_cell_num(CELL_array, C));
        */
        //LocalMinDist.printCellDist();

        cutnet = LocalMinDist.get_cutnet();

        assert(A.get_size() == LocalMinDist.get_A_size());
        assert(B.get_size() == LocalMinDist.get_B_size());


        if(alternate)
            FM_pass(C, N, r, i, CELL_array, NET_array, A, B, stuck, LocalMinDist, big_wave, alternate, block_num);
        else
            FM_pass(C, N, r, i, CELL_array, NET_array, A, B, stuck, LocalMinDist, big_wave, alternate);

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

        /*
        if(stuck && 0 <= (i + 1) % 128 && (i + 1) % 128 < stuck_criteria){
            printf("out from stuck!\n");
            StuckOut(A, B, CELL_array, NET_array, C, N, !A.bigger());
            BlockReinitialization(C, A, B, CELL_array, NET_array, i);
            LocalMinDist.update(CELL_array, c, A.get_size(), B.get_size(), CountCutNet(A, NET_array, N));
            stuck_count = 0;
            stuck = false;
        }
        */
                
        if(GlobalMinDist.get_cutnet() > LocalMinDist.get_cutnet()){
            GlobalMinDist = LocalMinDist;
            min_cutnet = GlobalMinDist.get_cutnet();
            stuck_count = 0;
            stuck = false;
            stuck_out_count++;
        }
        //printf("A size: %d, B size: %d\n", LocalMinDist.get_A_size(), LocalMinDist.get_B_size());
    }


    LoadDistribution(GlobalMinDist, A, B, CELL_array, C);
    BlockReinitialization(C, A, B, CELL_array, NET_array, -1);
    A.set_current_block_of_net(NET_array, N);
    B.set_current_block_of_net(NET_array, N);

    //printf("after fm cutnet: %d\n", CountCutNet(A, NET_array, N));

    total_cutnet = GlobalMinDist.get_cutnet();
    int next_step_1 = block_num / 2;
    int next_step_2 = block_num - next_step_1;

    if(next_step_1 > 1){
        total_cutnet += FM(InitVer, pass, CELL_array, NET_array, C, N, P, W, next_step_1, &A, skew, bias + next_step_2, alternate);
        for(int i = 1; i <= C; i++){
            if(GlobalMinDist[i] == 1){
                if(!_CELL_array[current_to_past[i]].get_fixed())
                    _CELL_array[current_to_past[i]].set_current_block_num(CELL_array[i].get_current_block_num());
            }
                //_CELL_array[current_to_past[i]].set_current_block_num(bias + CELL_array[i].get_current_block_num());
        }
    }
    else{
        for(int i = 1; i <= C; i++){
            if(GlobalMinDist[i] == 1)
                if(!_CELL_array[current_to_past[i]].get_fixed())
                    _CELL_array[current_to_past[i]].set_current_block_num(bias + next_step_2);
        }
        //printf("Block %d area: %d, cell_num: %d\n", bias + next_step_2, A.get_size(), A.get_cell_num(CELL_array, C));      
    }

    if(next_step_2 > 1){
        total_cutnet += FM(InitVer, pass, CELL_array, NET_array, C, N, P, W, next_step_2, &B, skew, bias, alternate);
        for(int i = 1; i <= C; i++){
            if(GlobalMinDist[i] == 0)
                if(!_CELL_array[current_to_past[i]].get_fixed())
                    _CELL_array[current_to_past[i]].set_current_block_num(CELL_array[i].get_current_block_num());
                //_CELL_array[current_to_past[i]].set_current_block_num(bias + CELL_array[i].get_current_block_num());
        }
    }
    else{
        for(int i = 1; i <= C; i++){
            if(GlobalMinDist[i] == 0)
                if(!_CELL_array[current_to_past[i]].get_fixed())
                    _CELL_array[current_to_past[i]].set_current_block_num(bias);
        }
        //printf("Block %d area: %d, cell_num: %d\n", bias, B.get_size(), B.get_cell_num(CELL_array, C));      
    }


    //printf("Final min num of cutnet: %d\n", global_min_cutnet);

    //Check(A, B, NET_array, N);

    delete[] current_to_past;
    if(current_block != nullptr){
        delete[] NET_array;
        delete[] CELL_array;
    }

    return total_cutnet;
}

int calculate_degree(Cell* &CELL_array, int C, Net* NET_array, int N, int block_num, int &cutnet){

    bool* block_check = new bool[block_num];

    int cutnet_check = 0;
    int degree = 0;
    int degree_temp = 0;

    for(int i = 1; i <= N; i++){
        for(int j = 0; j < block_num; j++)
            block_check[j] = false;

        auto itr = NET_array[i].cell_list.begin();
        block_check[(*itr)->get_current_block_num()] = true;
        degree_temp = 1;
        itr++;


        for(; itr != NET_array[i].cell_list.end(); itr++){
            if(!block_check[(*itr)->get_current_block_num()]){
                if(degree_temp == 1)
                    cutnet_check++;
                
                block_check[(*itr)->get_current_block_num()] = true;

                degree_temp++;
            }

            if(degree_temp == block_num)
                break;
        }

        if(degree_temp != 1)
            degree += degree_temp;
    }


    cutnet = cutnet_check;
    delete[] block_check;

    return degree;
}

void read_bin_record(const int _N, const int _C, const Net* _NET_array, const Cell* _CELL_array, int N, int C, Net* &NET_array, Cell* &CELL_array, const int current_bin, int &pin_num, int &total_weight, int *current_to_past){
    pin_num = 0;
    total_weight = 0;

    char c;
    int temp_cell_num;

    NET_array = new Net[N + 1];
    CELL_array = new Cell[C + 1];

    int* past_to_current = new int[_C + 1];

    /*
    printf("test1\n");
    int c_count = 0;
    for(int i = 1; i <= _C; i++){
        if(_CELL_array[i].get_bin() == current_bin)
            c_count++;
    }
    if(c_count != C)
        printf("fail %d %d\n", c_count, C);
    */
    int j = 1;
    for(int i = 1; i <= _C; i++){
        if(_CELL_array[i].get_bin() != current_bin)
            continue;
        
        CELL_array[j].set_cell_num(j);
        CELL_array[j].set_size(_CELL_array[i].get_size());
        total_weight += CELL_array[j].get_size();
        CELL_array[j].set_name(_CELL_array[i].get_cell_name());
        CELL_array[j].set_fixed(_CELL_array[i].get_fixed());
        CELL_array[j].set_current_block_num(_CELL_array[i].get_current_block_num());
        CELL_array[j].set_bin(current_bin);

        past_to_current[i] = j;
        current_to_past[j] = i;

        j++;
    }

    j = 0;

    for(int i = 1; i <= _N; i++){
        assert(j <= N);

        bool check = true;

        for(auto itr = _NET_array[i].cell_list.begin(); itr != _NET_array[i].cell_list.end(); itr++){
            if((*itr)->get_bin() == current_bin){
                if(check){
                    j++;
                    check = false;
                    NET_array[j].set_net_num(j);
                    NET_array[j].set_weight(_NET_array[i].get_weight());
                }

                temp_cell_num = past_to_current[(*itr)->get_cell_num()];
                NET_array[j].push_cell(CELL_array + temp_cell_num);
                CELL_array[temp_cell_num].push_net(NET_array + j);
                pin_num++;
            }
        }
    }

    /*
    printf("bin_size_test\n");
    int test_size = 0;
    for(int i = 1; i <= C; i++){
        test_size += CELL_array[i].get_size();

        if(CELL_array[i].get_bin() == -1)
            printf("wrong1\n");
    }

    if(test_size != total_weight)
        printf("wrong2!\n");
    */
    delete[] past_to_current;
}

void bin_based_FM(const int InitVer, const int pass, Cell* _CELL_array, Net* _NET_array, const int _C, const int _N, const int _P, const int _W, const int block_num, double skew, int map_n, int map_m, std::vector<int> *BIN_array){
    int C, N;
    int P, W; //P: total pin num, W: total weight
    
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;
    int* current_to_past = nullptr;

    int* net_in_bin = new int[map_n * map_m];

    int cutnet;

    for(int i = 0; i < map_n * map_m; i++)
        net_in_bin[i] = 0;
    
    for(int i = 1; i <= _N; i++){
        for(auto itr = _NET_array[i].cell_list.begin(); itr != _NET_array[i].cell_list.end(); itr++){
            net_in_bin[(*itr)->get_bin()]++;
        }
    }

    for(int i = 0; i < map_n * map_m; i++){
        C = BIN_array[i].size();
        N = net_in_bin[i];
        current_to_past = new int[C + 1];
        read_bin_record(_N, _C, _NET_array, _CELL_array, N, C, NET_array, CELL_array, i, P, W, current_to_past);

        cutnet = FM(InitVer, pass, CELL_array, NET_array, C, N, P * block_num, W, block_num, nullptr, skew, 0, false, i);

        printf("bin %d cutnet: %d\n", i, cutnet);

        for(int i = 1; i <= C; i++){
            _CELL_array[current_to_past[i]].set_current_block_num(CELL_array[i].get_current_block_num());
        }

        //check_place(C, CELL_array, "initPlace", map_n, map_m);

        delete[] NET_array;
        delete[] CELL_array;
        delete[] current_to_past;
    }

    /*
    printf("bin test\n");
    for(int i = 1; i <= C; i++){
        if(_CELL_array[i].get_bin() == -1)
            printf("error");
    }
    */
    delete[] net_in_bin;
}

void FM_pass(int C, int N, double r, int pass_num, Cell* CELL_array, Net* NET_array, Block &A, Block &B, const bool stuck, CellDist& LocalMinDist, bool big_wave, bool alternate, int block_num){
    Cell* BaseCell = nullptr;
    bool review = false;
    bool how_to_start = true;
    //init_pass_start의 값이 크면 결과는 빨리 나오나, 운적인 요소가 크게 작용
    int init_pass_start = (100 < (C / 50))? 100 : C / 50;
    //int init_pass_start = C / 50;
    int move_count = 0, move_count_at_min = 2;
    int A_size_at_min = A.get_size();
    const int out_from_stuck_citeria = 1;
    int out_from_stuck = out_from_stuck_citeria;
    int min_cutnet, cutnet;


    BlockReinitialization(C, A, B, CELL_array, NET_array, pass_num);

    min_cutnet = LocalMinDist.get_cutnet();
    //Block TempA = A;
    //Block TempB = B;

    /*
    if(stuck)
        printf("stuck!\n");
    else
        printf("not stuck!\n");
    */

    const bool destroy_balance = A.get_balance();
    const bool bigger = A.bigger();
    bool stuck_temp = stuck;
    Block* alternate_block;
    //bool destroy_balance = true;
    
    /*
    if(destroy_balance)
        printf("Destroy...!\n");
    else
        printf("balance\n");      
    */
    /*
    printf("cell 8 info\n");
    std::cout << CELL_array[8].get_current_block()->get_block_name() <<std::endl;
    printf("gain: %d\n", CELL_array[8].get_current_block()->gain[8]);
    */

    double choose_time = 0, move_time = 0;
    double get_max_gain_cell_time_temp;
    int pass_start;

    bool pair = false;
    int pair_bin = -1;
    int moved_cell = -1;
    do{ 
        //int pass_start = init_pass_start; //계속해서 초기화 시켜여줘야함에 유의!
        pass_start = (how_to_start && !alternate)? init_pass_start : 2;
        cutnet = LocalMinDist.get_cutnet();
        alternate_block = bigger? &A : &B;
        //printf("----------------------------");

        //assert(TempA == A);
        //assert(TempB == B);
        choose_time = 0, move_time = 0;

        pair = false;
        pair_bin = -1;
        moved_cell = -1;

        if(alternate && !pair){
            BaseCell = alternate_block->get_max_gain_cell();

            if(BaseCell != nullptr){
                alternate_block->remove_from_BUCKET(BaseCell);

                pair_bin = BaseCell->get_bin();
                moved_cell = BaseCell->get_cell_num();
                pair = true;

                alternate_block = (alternate_block == &A)? &B : &A;
            }
            /*
            else{
                printf("Something is wrong...\n");
                printf("A: %d, %f\n", A.get_size(), A.get_ubound());
                assert(BaseCell != nullptr);
            }
            */
        }
        else if(alternate && pair){
            int gain_max = -N * ALPHA;
            for(int i = pair_bin * block_num + 1; i <= (pair_bin + 1) * block_num; i++){
                if(CELL_array[i].get_current_block() == alternate_block){
                    if(alternate_block->gain[i] > gain_max && i != moved_cell && !CELL_array[i].locked){
                        gain_max = alternate_block->gain[i];
                        BaseCell = CELL_array + i;
                    }
                }

                assert(i - moved_cell < block_num);
            }

            pair = false;
            pair_bin = -1;
            moved_cell = -1;

            if(gain_max == -N * ALPHA)
                BaseCell = nullptr;

            //printf("BaseCell gain: %d\n", BaseCell->get_current_block()->gain[BaseCell->get_cell_num()]);

            if(BaseCell != nullptr)
                alternate_block->remove_from_BUCKET(BaseCell);

            alternate_block = (alternate_block == &A)? &B : &A;
        }
        else if(!stuck_temp)
            BaseCell = ChooseBaseCell_gain(A, B, r);
            //BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
        else{
            BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance, bigger);
        }


        while(BaseCell != nullptr){        
            //A.print_Block(CELL_array);
            //B.print_Block(CELL_array);

            /*
            printf("cell 11 gain: %d\n", CELL_array[11].get_current_block()->gain[11]);
            
            std::cout << BaseCell->get_current_block()->get_block_name() << " ";
            printf("BaseCell %d gain: %d\n", BaseCell->get_cell_num(), BaseCell->get_current_block()->gain[BaseCell->get_cell_num()]);
            */
            //printf("before move!\n");
            if(BaseCell->get_current_block() == &A){
                //printf("A gain %d\n", A.gain[BaseCell->get_cell_num()]);
                cutnet -= A.gain[BaseCell->get_cell_num()];
                if(!how_to_start && A.gain[BaseCell->get_cell_num()] > 0)
                    pass_start = false;
                MoveCell(A, B, BaseCell);
            }
            else{
                //printf("B gain %d\n", B.gain[BaseCell->get_cell_num()]);
                cutnet -= B.gain[BaseCell->get_cell_num()];
                if(!how_to_start && B.gain[BaseCell->get_cell_num()] > 0)
                    pass_start = false;
                MoveCell(B, A, BaseCell);
            }
            //printf("move!\n");
            if(review)
                move_count_at_min--;
            else
                move_count++;

            if(move_count_at_min == 0)
                break;       
            
            /*
            printf("%d %d\n", cutnet, CountCutNet(A, NET_array, N));
            printf("move_count: %d\n", move_count);
            */


            //printf("cc: %d %d\n", cutnet, CountCutNet(A, NET_array, N));
            assert(cutnet == CountCutNet(A, NET_array, N));
            if(cutnet < 0){
                //printf("negative cutnet error! move_count: %d\n", move_count);
                assert(cutnet >= 0);
            }

            if(pass_start){    
                if(!review){              
                    min_cutnet = cutnet;
                    move_count_at_min = move_count;
                    //printf("%d ", move_count_at_min);
                    A_size_at_min = A.get_size();

                    if(how_to_start)
                        pass_start--;
                }
            }
            else if(!review){
                if(!alternate){
                    if(min_cutnet > cutnet){
                        min_cutnet = cutnet;
                        move_count_at_min = move_count;
                        A_size_at_min = A.get_size();
                    }
                }
                else{
                    if(min_cutnet > cutnet && move_count %2 == 0){
                        min_cutnet = cutnet;
                        move_count_at_min = move_count;
                        A_size_at_min = A.get_size();
                    }
                }
            }


            if(alternate && !pair){
                BaseCell = alternate_block->get_max_gain_cell();

                if(BaseCell != nullptr){
                    alternate_block->remove_from_BUCKET(BaseCell);

                    pair_bin = BaseCell->get_bin();
                    moved_cell = BaseCell->get_cell_num();
                    pair = true;
                }
                /*
                else{
                    printf("Something is wrong...\n");
                    printf("A: %d, %f\n", A.get_size(), A.get_ubound());
                    assert(BaseCell != nullptr);
                }
                */
                alternate_block = (alternate_block == &A)? &B : &A;
            }
            else if(alternate && pair){
                int gain_max = -N  * ALPHA;
                for(int i = pair_bin * block_num + 1; i <= (pair_bin + 1) * block_num; i++){
                    if(CELL_array[i].get_current_block() == alternate_block){
                        if(alternate_block->gain[i] > gain_max && i != moved_cell && !CELL_array[i].locked){
                            gain_max = alternate_block->gain[i];
                            BaseCell = CELL_array + i;
                        }
                    }

                    //printf("%d, %d, %d, %d\n", pair_bin, i, moved_cell, block_num);
                    assert(i - moved_cell < block_num);
                }

                if(gain_max == -N * ALPHA)
                    BaseCell = nullptr;

                //printf("basecell gain: %d\n", BaseCell->get_current_block()->gain[BaseCell->get_cell_num()]);

                if(BaseCell != nullptr){
                    alternate_block->remove_from_BUCKET(BaseCell);
                    pair = false;
                    pair_bin = -1;
                    moved_cell = -1;
                    alternate_block = (alternate_block == &A)? &B : &A;
                }
            }
            else if(!stuck_temp)
                BaseCell = ChooseBaseCell_gain(A, B, r);
                //BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance);
            else{
                BaseCell = ChooseBaseCell_balance(A, B, r, destroy_balance, bigger);
            }

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
            if(stuck_temp && A.bigger() != bigger && A.get_balance() != destroy_balance && (!alternate || move_count %2 == 0)){
                if(!review){
                    min_cutnet = cutnet;
                    move_count_at_min = move_count;
                }
                stuck_temp = false;
            }

            

            if(big_wave && stuck && !stuck_temp && (A.bigger() == bigger) && (!alternate || move_count %2 == 0)){
                stuck_temp = true;

                if(move_count > C / 4){
                    move_count = 0;
                    break;
                }
            }
            
        
            if(move_count > C - C / 8 && move_count % 2 == 0){
               move_count = 0;
               break;
            }
           
        }

        //printf("move_count_at_min: %d\n", move_count_at_min);
        assert(!alternate || move_count_at_min % 2 == 0);

        //printf("choose + move time: %fs, move_count_at_min: %d\n", (clock() - (double)choosemove_time_temp) / CLOCKS_PER_SEC, move_count_at_min);
        //printf("choose time: %fs, move time: %fs\n", choose_time/ CLOCKS_PER_SEC, move_time / CLOCKS_PER_SEC);
        //printf("get_max_gain_cell_time: %fs\n", (get_max_gain_cell_time - get_max_gain_cell_time_temp)/ CLOCKS_PER_SEC);
        //printf("review\n");

        if(!review){
            LoadDistribution(LocalMinDist, A, B, CELL_array, C);
            if(stuck != stuck_temp)
                stuck_temp = stuck;
            
            BlockReinitialization(C, A, B, CELL_array, NET_array, pass_num);
        }
        //printf("reinit time: %fs\n\n", (clock() - (double)reinit_time_temp) / CLOCKS_PER_SEC);

        out_from_stuck = out_from_stuck_citeria;
        stuck_temp = stuck;
        review = !review;
    }while(review);

    //printf("min_cutnet: %d, cutnet: %d\n", min_cutnet, cutnet);
    assert(min_cutnet == cutnet);

    //printf("\n final \n");
    //A.print_Block(CELL_array);
    //B.print_Block(CELL_array);
    

    LocalMinDist.overWrite(CELL_array, C, NET_array, N, min_cutnet);
}

int CountOverlap(bool **map, int N, Net *NET_array, int ll_x, int ll_y, int ur_x, int ur_y, int block_num){
    int overlap_count = 0;
    
    for(int i = 1; i <= N; i++){
        for(int j = 0; j < block_num; j++)
            overlap_count += NET_array[i].count_overlap(map, ll_x, ll_y, j);
    }

    printf("Total overlap: %d\n", overlap_count);

    return overlap_count;
}

int CalculateTotalOverlap(Net* NET_array, int N, int block_num){
    bool* block_check = new bool[block_num];
    int overlap_count = 0;

    for(int i = 1; i <= N; i++){
        for(int i = 0; i < block_num; i++)
            block_check[i] = false;
        
        for(auto itr = NET_array[i].cell_list.begin(); itr != NET_array[i].cell_list.end(); itr++){
            if(block_check[(*itr)->get_current_block_num()])
                overlap_count++;
            else
                block_check[(*itr)->get_current_block_num()] = true;

        }
    }

    delete[] block_check;

    return overlap_count;
}