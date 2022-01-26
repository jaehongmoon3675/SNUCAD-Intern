#include <iostream>
#include <cmath>
#include <string>
#include <list>
#include <stack>

#include "Cell.h"
#include "Net.h"

std::stack<Cell*> FreeCellList;

/*
//Using CellNode to record the best distribution
struct CellNode{
    std::string* name;
    int* distribution;
    int length;
    int num_cutnet;
    int balance;
    CellNode(int C) : length(C), num_cutnet(0), balance(0) {
        name = new std::string[C + 1];
        distribution = new int[C + 1];
    }
    void operator=(const Cell* CELL_array){ //무엇이 올지는 좀 나중에 결정해야겠다...
        delete[] name;
        delete[] distribution;

        name = new std::string[copy.length];
        distribution = new std::string
    }
    ~CellNode() {
        delete[] name;
        delete[] distribution;
    }
};
*/

class Block{
public:    
    Block(int init_pmax, int low_bound, int up_bound, int c, int n, int w, int r)
        : pmax(init_pmax), lbound(low_bound), ubound(up_bound), C(c), N(n), W(w), R(r), size(0) {
        BUCKET = new Cell*[2 * pmax + 1];
        Fdistribution = new int[N + 1];
        Ldistribution = new int[N + 1];
        gain = new int[C + 1];
    }
    void CalculateDistribution(Cell* CELL_array){
        for(int i = 1; i <= N; i++){
            Fdistribution[i] = 0;
            Ldistribution[i] = 0;
        }

        for(int i = 1; i <= C; i++){
            if(this == CELL_array[i].get_current_block()){
                CELL_array[i].get_net_list(Fdistribution);
            }
        }
    }
    void CellGainAdjustment(Cell &c){ //inner loop of implementation of the code prior to Proposition 2
        int gain_adjust;
        if(this == c.get_current_block()){ //if this block is F (from block)
            for(auto i = c.net_list.begin(); i != c.net_list.end(); i++){
                if(Fdistribution[(*i)->get_net_num()] == 1)
                    gain_adjust++;
            }
        } 
        else{
            for(auto i = c.net_list.begin(); i != c.net_list.end(); i++){ //if this block is T (to block)
                if(Fdistribution[(*i)->get_net_num()] == 0)
                    gain_adjust--;
            }
        }

        gain[c.cell_num] = gain_adjust;
        c.BUCKETpre = nullptr;
        c.BUCKETnext = BUCKET[gain_adjust];
        BUCKET[gain_adjust]->BUCKETpre = &c;
        BUCKET[gain_adjust] = &c;

        if(gain_adjust >= pmax)
            pmax = gain_adjust;
    }
    Cell* get_max_gain_cell() const{
        int i = pmax;
        double modified_r;
        Cell* max_gain_cell = BUCKET[pmax];

        do{
            modified_r = ((double)size + max_gain_cell->get_size()) / W;

            if(lbound <= modified_r && modified_r <= ubound)
                break;
            
            max_gain_cell = max_gain_cell ->BUCKETnext;

            if(max_gain_cell == nullptr){
                i--;
                
                if(BUCKET[i] == nullptr)
                    break;
                else
                    max_gain_cell = BUCKET[i];
            }
        }while(BUCKET[i] != nullptr);

        return max_gain_cell;
    }
    double get_modified_balance_factor(Cell* cell){
        double modified_r = ((double)size + cell->get_size()) / W;

        return modified_r;
    }
    ~Block(){
        delete[] BUCKET;
        delete[] Fdistribution; //Free Distribution
        delete[] Ldistribution; //Locked Distribution
        delete[] gain;
    }
private:
    Cell** BUCKET;
    int* Fdistribution, * Ldistribution;
    int* gain;
    int pmax;
    int lbound, ubound;
    int size;
    const int C, N, W, R;
};

//block의 사이즈도 여기서 계산해주어야 한다.
void BlockInitialization(Block &A, Block &B, Cell* CELL_array){

}

//implement how to choose the base cell
Cell* ChooseBaseCell(Block &A, Block &B, int r){ //r is a balance factor
    Cell* max_gain_cellA, * max_gain_cellB;
    max_gain_cellA = A.get_max_gain_cell();
    max_gain_cellB = B.get_max_gain_cell();

    if(max_gain_cellA == nullptr && max_gain_cellB == nullptr)
        return nullptr;
    else if(max_gain_cellA == nullptr)
        return max_gain_cellB;
    else if(max_gain_cellB == nullptr)
        return max_gain_cellB;
    else{
        if(std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) < std::abs(B.get_modified_balance_factor(max_gain_cellB) - r))
            return max_gain_cellA;
        else
            return max_gain_cellB;
    }
}

//implementation of the code prior to Proposition 2
void BlockReinitialization(Block &A, Block &B, Cell* CELL_array){
    A.CalculateDistribution(CELL_array);
    B.CalculateDistribution(CELL_array);
    
    Cell *ctemp;

    while(!FreeCellList.empty()){
        ctemp = FreeCellList.top();
        FreeCellList.pop();

        A.CellGainAdjustment(*ctemp);
        B.CellGainAdjustment(*ctemp);
    }
}