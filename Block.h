#ifndef __BLOCK_H__
#define __BLOCK_H__

extern int ALPHA;

#include "Cell.h"
#include "Net.h"

class Block{
public:
    int* Fdistribution, * Ldistribution;
    int* gain;
    bool init_error;

    Block(int init_pmax, double low_bound, double up_bound, int c, int n, int w, double r, int _block_num_lb, int _block_num_ub, std::string block_name)
        : PMAX(20 * init_pmax * ALPHA), max_gain(-PMAX), lbound(low_bound), ubound(up_bound), C(c), N(n), W(w), R(r), block_num_lb(_block_num_lb), block_num_ub(_block_num_ub), size(0), name(block_name), init_error(true) {
        BUCKET = new Cell*[2 * PMAX + 1];
        
        for(int i = 0; i < 2 * PMAX + 1; i++)
            BUCKET[i] = nullptr;
        
        BUCKET += PMAX;

        Fdistribution = new int[N + 1];
        Ldistribution = new int[N + 1];
        gain = new int[C + 1];
    }
    Block(const Block&);
    bool operator==(const Block& compare) const;
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
    void increase_cell_gain(Cell *cell, int weight = 1);
    void decrease_cell_gain(Cell *cell, int weight = 1);
    void increase_cell_gain_of_net(Net* net);
    void decrease_cell_gain_of_net(Net* net);
    Cell* find_cell_in_block(Net* net);
    bool push_Cell_ub(Cell* cell); //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환
    bool push_Cell_r(Cell* cell); //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환
    void print_Block(Cell* CELL_array);
    void print_Block_short(Cell* CELL_array);
    void empty_BUCKET();
    int get_cell_num(const Cell* CELL_array, int C) const;
    int get_uncut_count(const Net* NET_array, int N) const;
    void set_current_block_of_net(Net* NET_array, int N);
    int ith_net_distribution(int i){
        return Fdistribution[i] + Ldistribution[i];
    }
    int get_W() const { return W; }
    int get_size() const { return size; }
    double get_ubound() const { return ubound; }
    void set_size(int _size) { size = _size; }
    bool get_balance() const {
        if(std::abs(size - R*W) < std::abs(ubound - R*W) / 2)
            return true;
        else
            return false;
    }
    bool check_block() const;
    bool bigger() const{
        if(size > R*W)
            return true;
        else
            return false;
    }
    void set_ubound(double _ubound) { ubound = _ubound; }
    bool is_this_fixed_cell_here(int test_block_num){
        if(block_num_lb <= test_block_num && test_block_num <= block_num_ub)
            return true;
        
        return false;
    }
    int get_block_num_lb() const { return block_num_lb; }
    int get_block_num_ub() const { return block_num_ub; }
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
    int block_num_ub, block_num_lb;
};

//VERSION1
//block의 사이즈도 여기서 계산해주어야 한다. BlockInitialization 실행 후 Reinitialization도 실행시켜주어야..
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C);
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C, int bin_num);
void BlockInitialization_r(Block &A, Block &B, Cell* CELL_array, int C);
void BlockInitialization_cell_bin(Block &A, Block &B, Cell* CELL_array, int C, int block_num, int bin_num = 0);

//VERSION 2
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N);
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, int bin_num);

//implement how to choose the base cell, find base cell, remove it from block and push it into FreeCellList
Cell* ChooseBaseCell_gain(Block &A, Block &B, double r); //r is a balance factor
Cell* ChooseBaseCell_balance(Block &A, Block &B, double r, bool destroy_balance); //r is a balance factor
Cell* ChooseBaseCell_balance(Block &A, Block &B, double r, bool destroy_balance, bool which_block);

//implementation of the code prior to Proposition 2
void BlockReinitialization(int C, Block &A, Block &B, Cell* CELL_array, Net* NET_array, int pass);

void MoveCell(Block &F, Block &T, Cell* BaseCell);

void StuckOut(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, bool bigger);

#endif