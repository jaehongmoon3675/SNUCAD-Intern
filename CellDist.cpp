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

bool CellDist::update(Cell* CELL_array, int C, int _A_size, int _B_size, int _A_count, int _B_count, int _cutnet){
    if(_cutnet > cutnet)
        return false;
    
    if((_cutnet == cutnet) && (std::abs(_A_size - ideal_balance) > std::abs(A_size - ideal_balance)))
        return false;
    
    for(int i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() == BlockA)
            distribution[i] = 1;
        else
            distribution[i] = 0;
    }

    A_size = _A_size;
    B_size = _B_size;
    A_count = _A_count;
    B_count = _B_count;
    cutnet = _cutnet;

    return true;
}

void CellDist::writeCellDist(Cell* CELL_array, int C) const{
    std::string filePath = "output.part";
    std::ofstream writeFile(filePath.data());

    if(writeFile.is_open()){
        for(int i = 1; i <= C; i++){
            writeFile << CELL_array[i].get_cell_name() << " ";

            if(CELL_array[i].get_current_block() == BlockA)
                writeFile << "1" << std::endl;
            else
                writeFile << "0" << std::endl;
        }

        writeFile.close();
    }
}

void CellDist::printCellDist() const{
    printf("cutnet: %d\n", cutnet);
    printf("size of A: %d, size of B: %d\n", A_size, B_size);
    printf("num of cell in A: %d, num of cell in B: %d\n", A_count, B_count);
}

void LoadDistribution(CellDist &Distribution, Cell* CELL_array, int C, int &cutnet){
    cutnet = Distribution.get_cutnet();
    Distribution.setBlockSize();
    for(int i = 1; i <= C; i++){
        CELL_array[i].set_current_block(Distribution.get_ith_cell_current_block(i));
    }
}