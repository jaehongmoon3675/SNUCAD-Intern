#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "Cell.h"
#include "Net.h"

class Block{
public:    
    Block(int init_pmax, double low_bound, double up_bound, int c, int n, int w, double r, std::string block_name)
        : PMAX(init_pmax), max_gain(-PMAX), lbound(low_bound), ubound(up_bound), C(c), N(n), W(w), R(r), size(0), name(block_name) {
        BUCKET = new Cell*[2 * PMAX + 1];
        
        for(int i = 0; i < 2 * PMAX + 1; i++)
            BUCKET[i] = nullptr;
        
        BUCKET += PMAX;

        Fdistribution = new int[N + 1];
        Ldistribution = new int[N + 1];
        gain = new int[C + 1];
    }
    void CalculateDistribution(Cell* CELL_array);
    void CellGainAdjustment(Block &T, Cell &c); //inner loop of implementation of the code prior to Proposition 2
    Cell* get_max_gain_cell() const;
    void remove_from_BUCKET(Cell* cell);
    double get_modified_balance_factor(Cell* cell){
        double modified_r = ((double)size + cell->get_size()) / W;

        return modified_r;
    }
    std::string get_block_name(){
        return name;
    }
    bool push_Cell(Cell* cell); //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환
    void print_Block();
    ~Block(){
        delete[] (BUCKET - PMAX);
        delete[] Fdistribution; //Free Distribution
        delete[] Ldistribution; //Locked Distribution
        delete[] gain;
    }
private:
    Cell** BUCKET;
    int* Fdistribution, * Ldistribution;
    int* gain;
    const int PMAX;
    int max_gain;
    double lbound, ubound;
    int size;
    std::string name;
    const int C, N, W;
    const double R;
};

//block의 사이즈도 여기서 계산해주어야 한다. BlockInitialization 실행 후 Reinitialization도 실행시켜주어야..
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C);

//implement how to choose the base cell, find base cell, remove it from block and push it into FreeCellList
Cell* ChooseBaseCell(Block &A, Block &B, int r); //r is a balance factor

//implementation of the code prior to Proposition 2
void BlockReinitialization(Block &A, Block &B, Cell* CELL_array);

#endif