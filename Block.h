#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "Cell.h"
#include "Net.h"

class Block{
public:
    int* Fdistribution, * Ldistribution;
    int* gain;

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
    void CellGainInitialization(Block &T, Cell &c); //inner loop of implementation of the code prior to Proposition 2
    Cell* get_max_gain_cell() const;
    void remove_from_BUCKET(Cell* cell);
    void adjust_maxgain();
    void add_size(int cell_size){ size += cell_size; }
    double get_modified_balance_factor(Cell* cell){
        double modified_r = ((double)size - cell->get_size()) / W;

        return modified_r;
    }
    std::string get_block_name(){
        return name;
    }
    void increase_cell_gain(Cell *cell);
    void decrease_cell_gain(Cell *cell);
    void increase_cell_gain_of_net(Net* net);
    void decrease_cell_gain_of_net(Net* net);
    Cell* find_cell_in_block(Net* net);
    bool push_Cell(Cell* cell); //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환
    void print_Block(Cell* CELL_array);
    void print_Block_short(Cell* CELL_array);
    void empty_BUCKET();
    int get_cell_num(Cell* CELL_array, int C);
    int ith_net_distribution(int i){
        return Fdistribution[i] + Ldistribution[i];
    }
    int get_size() { return size; }
    void set_size(int _size) { size = _size; }
    ~Block(){
        delete[] (BUCKET - PMAX);
        delete[] Fdistribution; //Free Distribution
        delete[] Ldistribution; //Locked Distribution
        delete[] gain;
    }
private:
    Cell** BUCKET;
    const int PMAX;
    int max_gain;
    double lbound, ubound;
    int size;
    std::string name;
    const int C, N, W;
    const double R;
};

//VERSION1
//block의 사이즈도 여기서 계산해주어야 한다. BlockInitialization 실행 후 Reinitialization도 실행시켜주어야..
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C);

//VERSION 2 Ver3
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N);
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, int ver3);

//implement how to choose the base cell, find base cell, remove it from block and push it into FreeCellList
Cell* ChooseBaseCell(Block &A, Block &B, double r); //r is a balance factor

//implementation of the code prior to Proposition 2
void BlockReinitialization(Block &A, Block &B, Cell* CELL_array);

void MoveCell(Block &F, Block &T, Cell* BaseCell);

#endif