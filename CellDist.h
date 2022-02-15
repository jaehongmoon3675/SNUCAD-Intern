#ifndef __CELLDIST_H__
#define __CELLDIST_H__

#include "Cell.h"
#include "Net.h"
#include "Block.h"

class CellDist{
public:
    CellDist(int _C, int _N, int _ideal_balance, int _cutnet, int _A_size, int _B_size, Block* _BlockA, Block* _BlockB) 
        : C(_C), N(_N), ideal_balance(_ideal_balance), cutnet(_cutnet), A_size(_A_size), B_size(_B_size), A_count(0), B_count(0), BlockA(_BlockA), BlockB(_BlockB) {
        distribution = new int[_C + 1];
    }
    CellDist(const CellDist& copy){
        C = copy.C;
        N = copy.N;
        cutnet = copy.cutnet;
        ideal_balance = copy.ideal_balance;
        A_size = copy.A_size;
        B_size = copy.B_size;
        A_count = copy.A_count;
        B_count = copy.B_count;
        BlockA = copy.BlockA;
        BlockB = copy.BlockB;
        
        distribution = new int[C + 1];

        for(int i = 1; i <= C; i++)
            distribution[i] = copy.distribution[i];
    }
    bool update(Cell* CELL_array, int C, int _A_size, int _B_size, int _cutnet);
    void overWrite(Cell* CELL_array, int C, int _A_size, int _B_size, int _cutnet);
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
    int C, N;
    int cutnet;
    int ideal_balance;
    int A_size, B_size;
    int A_count, B_count;
    Block* BlockA, * BlockB;
};

void LoadDistribution(CellDist &Distribution, Cell* CELL_array, int C, int& cutnet);

#endif