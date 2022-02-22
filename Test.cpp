#include <iostream>
#include <fstream>
#include <list>
#include <stack>
#include <algorithm>

#include "Cell.h"
#include "Net.h"
#include "ReadFile.h"
#include "FM_func.h"
#include "Block.h"
#include "WriteFile.h"

int gain_update_count = 0;

extern std::stack<Cell*> FreeCellList;

void read_output_part(Block &A, Block &B, const int C, Cell* &CELL_array, int &A_cell_num, int &B_cell_num){
    std::ifstream readFile;
    readFile.open("aes_128_1_500.part");

    //printf("read_hgr_map\n");

    int block;
    std::string temp_cell_name;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_name >> block;
            
            if(block == 1){
                CELL_array[i].set_current_block(&A);
                A.add_size(CELL_array[i].get_size());
                A_cell_num++;
            }
            else{
                CELL_array[i].set_current_block(&B);
                B.add_size(CELL_array[i].get_size());
                B_cell_num++;
            }
            FreeCellList.push(CELL_array + i);
        }

        readFile.close();
    }
    else
        printf("No map file\n");
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
    read_output_part(A, B, C, CELL_array, A_cell_num, B_cell_num);
    BlockReinitialization(C, A, B, CELL_array, NET_array, 0);

    
    printf("A size: %d, B size: %d\n", A.get_size(), B.get_size());


    printf("cutnet_num: %d\n", CountCutNet(A, NET_array, N));
    printf("A cell num: %d, B cell num: %d\n", A_cell_num, B_cell_num);


    delete[] NET_array;
    delete[] CELL_array;

    return 0;
}