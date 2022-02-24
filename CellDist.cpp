#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <list>
#include <stack>
#include <queue>

#include "Cell.h"
#include "Net.h"
#include "Block.h"
#include "CellDist.h"

extern std::stack<Cell*> FreeCellList;

bool CellDist::update(Cell* CELL_array, int C, Net* NET_array, int N, int _cutnet){
    if(_cutnet > cutnet)
        return false;
    
    /*
    if((_cutnet == cutnet) && (std::abs(_A_size - ideal_balance) >= std::abs(A_size - ideal_balance)))
        return false;
    */
   
    for(int i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() == A)
            distribution[i] = 1;
        else
            distribution[i] = 0;
    }

    A_size = A->get_size();
    B_size = B->get_size();
    A_count = A->get_cell_num(CELL_array, C);
    B_count = B->get_cell_num(CELL_array, C);
    A_uncut_count = A->get_uncut_count(NET_array, N);
    B_uncut_count = B->get_uncut_count(NET_array, N);
    cutnet = _cutnet;

    return true;
}

void CellDist::overWrite(Cell* CELL_array, int C, Net* NET_array, int N, int _cutnet){
    for(int i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() == A)
            distribution[i] = 1;
        else
            distribution[i] = 0;
    }

    A_size = A->get_size();
    B_size = B->get_size();
    A_count = A->get_cell_num(CELL_array, C);
    B_count = B->get_cell_num(CELL_array, C);
    A_uncut_count = A->get_uncut_count(NET_array, N);
    B_uncut_count = B->get_uncut_count(NET_array, N);
    cutnet = _cutnet;

    return;
}

//수정 필요
void CellDist::writeCellDist(Cell* CELL_array, int C, std::string _filename, int init_num, int pass, int bias) const{
    std::string filePath = _filename + "_" + std::to_string(init_num) + "_" + std::to_string(pass) + ".part";
    std::ofstream writeFile(filePath.data());

    if(writeFile.is_open()){
        for(int i = 1; i <= C; i++){
            writeFile << CELL_array[i].get_cell_name() << " ";

            if(distribution[i] == 1)
                writeFile << bias + 1 << std::endl;
            else
                writeFile << bias << std::endl;
        }

        writeFile.close();
    }
}

void CellDist::printCellDist() const{
    printf("Final min num of cutnet: %d\n", cutnet);
    printf("size of A: %d, size of B: %d\n", A_size, B_size);
    printf("num of cell in A: %d, num of cell in B: %d\n", A_count, B_count);
}

void LoadDistribution(CellDist &Distribution, Block &A, Block &B, Cell* CELL_array, int C){    
    Distribution.setBlockSize();
    for(int i = 1; i <= C; i++){
        CELL_array[i].set_current_block(Distribution.get_ith_cell_current_block(i));
    }
}