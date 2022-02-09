#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "Cell.h"
#include "Net.h"

struct CellDist{
    int* distribution;
    int cutnet;
    int ideal_balance;
    int A_size, B_size;
    int A_count, B_count;
    Block* BlockA, * BlockB;
    CellDist(int C, int _ideal_balance, int _cutnet, int _A_size, int _B_size, int _A_count, int _B_count, Block* _BlockA, Block* _BlockB) 
        : ideal_balance(_ideal_balance), cutnet(_cutnet), A_size(_A_size), B_size(_B_size), A_count(_A_count), B_count(_B_count), BlockA(_BlockA), BlockB(_BlockB) {
        distribution = new int[C + 1];
    }
    bool update(Cell* CELL_array, int C, int _A_size, int _B_size, int _A_count, int _B_count, int _cutnet);
    void writeCellDist(Cell* CELL_array, int C);
    Block* get_ith_cell_current_block(int i){
        if(distribution[i] == 1)
            return BlockA;
        else
            return BlockB;
    }
    ~CellDist() {
        delete[] distribution;
    }
};

class Block{
public:
    int* Fdistribution, * Ldistribution;
    int* gain;

    Block(int init_pmax, double low_bound, double up_bound, int c, int n, int w, double r, std::string block_name)
        : PMAX(init_pmax), max_gain(-PMAX), lbound(low_bound), ubound(up_bound), cell_count(0), C(c), N(n), W(w), R(r), size(0), name(block_name) {
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
    void set_count_0(){ cell_count = 0; }
    int ith_net_distribution(int i){
        return Fdistribution[i] + Ldistribution[i];
    }
    void increase_cell_count() { cell_count++; }
    void decrease_cell_count() { cell_count--; }
    int get_size() { return size; }
    int get_cell_count() { return cell_count; }
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
    int cell_count;
    std::string name;
    const int C, N, W;
    const double R;
};

void LoadDistribution(CellDist &Distribution, Cell* CELL_array, int C);

//VERSION1
//block의 사이즈도 여기서 계산해주어야 한다. BlockInitialization 실행 후 Reinitialization도 실행시켜주어야..
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C);

//VERSION 2 Ver3
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N);

//implement how to choose the base cell, find base cell, remove it from block and push it into FreeCellList
Cell* ChooseBaseCell(Block &A, Block &B, double r); //r is a balance factor

//implementation of the code prior to Proposition 2
void BlockReinitialization(Block &A, Block &B, Cell* CELL_array);

void MoveCell(Block &F, Block &T, Cell* BaseCell);

#endif