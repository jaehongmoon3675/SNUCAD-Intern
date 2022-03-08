#include <iostream>
#include <fstream>
#include <list>
#include <stack>
#include <algorithm>
#include <vector>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "WriteFile.h"

int gain_update_count = 0;



int read_output_part(Cell* &CELL_array, int C, Net* NET_array, int N){
    std::ifstream readFile;
    readFile.open("output.part");

    int block;
    std::string temp_cell_name;
    int* CELL_block_info = new int[C + 1];

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_name >> block;
            
            CELL_block_info[i] = block;
        }

        readFile.close();
    }
    else
        printf("No map file\n");

    int cutnet = 0;

    for(int i = 1; i <= N; i++){
        auto itr = NET_array[i].cell_list.begin();
        int current_block_num = CELL_block_info[(*itr)->get_cell_num()];
        itr++;

        for(; itr != NET_array[i].cell_list.end(); itr++){
            if(current_block_num != CELL_block_info[(*itr)->get_cell_num()]){
                cutnet++;
                break;
            }
        }
    }

    delete[] CELL_block_info;

    return cutnet;
}

int main(){
    int N, C; //the num of net and cell, respectively
    int P, W; //P: total pin num, W: total weight
    double r = 0.5; //balance factor
    int pass = 20; //how many pass we go through
    const int file_num = 0;
    int min_cutset_num;
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;
    const std::string file_name_arr[4] = {"aes_128", "ldpc", "jpeg", "initPlace"};
    std::string file_name = file_name_arr[file_num];


    P = read_hgr(N, C, NET_array, CELL_array, file_name);
    read_hgr_map(C, CELL_array, file_name);
    W = read_hgr_area(C, CELL_array, file_name); 
    
    int pmax = 0, smax = 0;
    double balance_low_bound = 0, balance_up_bound = 0;
    get_max(C, CELL_array, pmax, smax);
    balance_low_bound = r*W - smax * 10;
    balance_up_bound = r*W + smax * 10;

    printf("N: %d C: %d, P: %d, W: %d, R: %f\n", N, C, P, W, r);
    printf("pmax: %d, smax: %d\n", pmax, smax);

    Block A(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "A");
    Block B(pmax, balance_low_bound, balance_up_bound, C, N, W, r, "B");
    int A_cell_num = 0, B_cell_num = 0;
    int cutnet = read_output_part(CELL_array, C, NET_array, N);

    printf("cutnet: %d\n", cutnet);


    delete[] NET_array;
    delete[] CELL_array;

    return 0;
}