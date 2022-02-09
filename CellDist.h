#ifndef __CELLDIST_H__
#define __CELLDIST_H__

#include "Cell.h"
#include "Net.h"
#include "Block.h"

class CellDist{
public:
    CellDist(int C, int _ideal_balance, int _cutnet, int _A_size, int _B_size, int _A_count, int _B_count, Block* _BlockA, Block* _BlockB) 
        : ideal_balance(_ideal_balance), cutnet(_cutnet), A_size(_A_size), B_size(_B_size), A_count(_A_count), B_count(_B_count), BlockA(_BlockA), BlockB(_BlockB) {
        distribution = new int[C + 1];
    }
    bool update(Cell* CELL_array, int C, int _A_size, int _B_size, int _A_count, int _B_count, int _cutnet);
    void overWrite(Cell* CELL_array, int C, int _A_size, int _B_size, int _A_count, int _B_count, int _cutnet);
    void writeCellDist(Cell* CELL_array, int C) const;
    void printCellDist() const;
    int get_A_size() const { return A_size; }
    int get_B_size() const { return B_size; }
    int get_cutnet() const { return cutnet; }
    void setBlockSize(){
        BlockA->set_size(A_size);
        BlockB->set_size(B_size);
    }
    Block* get_ith_cell_current_block(int i){
        if(distribution[i] == 1)
            return BlockA;
        else
            return BlockB;
    }
    ~CellDist() {
        delete[] distribution;
    }
private:
    int* distribution;
    int cutnet;
    int ideal_balance;
    int A_size, B_size;
    int A_count, B_count;
    Block* BlockA, * BlockB;
};

void LoadDistribution(CellDist &Distribution, Cell* CELL_array, int C, int& cutnet);

#endif