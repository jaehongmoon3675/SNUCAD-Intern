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
        : PMAX(init_pmax), max_gain(-PMAX), lbound(low_bound), ubound(up_bound), C(c), N(n), W(w), R(r), size(0) {
        BUCKET = (new Cell*[2 * PMAX + 1]) + PMAX;
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

        if(gain_adjust >= max_gain)
            max_gain = gain_adjust;
    }
    Cell* get_max_gain_cell() const{
        int i = max_gain;
        double modified_r;
        Cell* max_gain_cell = BUCKET[max_gain];

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
    void remove_from_BUCKET(Cell* cell){
        if(this != cell->get_current_block())
            printf("error on remove_from_BUCKET() of Block Class\n");
        
        if(cell->BUCKETpre == nullptr)
            BUCKET[gain[cell->cell_num]] = cell->BUCKETnext;
        else
            cell->BUCKETpre->BUCKETnext = cell->BUCKETnext;
        
        if(cell->BUCKETnext != nullptr)
            cell->BUCKETnext->BUCKETpre = cell->BUCKETpre;
    }
    double get_modified_balance_factor(Cell* cell){
        double modified_r = ((double)size + cell->get_size()) / W;

        return modified_r;
    }
    bool push_Cell(Cell* cell){ //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환
        if(size + cell->get_size() <= ubound){
            cell->set_current_block(this);
            size += cell->get_size();

            return true;
        }
        else{
            return false;
        }
    }
    void print_Block(){
        printf("------------------------------------------\n");
        printf("max_gain: %d, lbound: %d , ubound: %d, size: %d\n", max_gain, lbound, ubound, size);
        
        printf("------------------------------------------\n");
        printf("Fdistribution & Ldistribution\n");
        for(int i = 1; i <= N; i++)
            printf("net %d: %d ", i, Fdistribution[i]);
        printf("\n");
        for(int i = 1; i <= N; i++)
            printf("net %d: %d ", i, Ldistribution[i]);
        printf("\n");

        printf("------------------------------------------\n");
        printf("gain of cell\n");
        for(int i = 1; i <= C; i++)
            printf("cell %d: %d", i, gain[i]);
        printf("\n");

        printf("------------------------------------------\n");
        printf("BUCKET");
        
        int head = max_gain;
        Cell *cell_head = nullptr;

        while(BUCKET[head] != nullptr){
            cell_head = BUCKET[head];

            printf("gain %d: ", head);

            while(cell_head != nullptr){
                printf("%d ", cell_head->cell_num);
                cell_head = cell_head->BUCKETnext;
            }

            printf("\n");
            head--;
        }
    }
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
    int lbound, ubound;
    int size;
    const int C, N, W, R;
};

//block의 사이즈도 여기서 계산해주어야 한다. BlockInitialization 실행 후 Reinitialization도 실행시켜주어야..
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C){
    int i;

    for(i = 1; i <= C; i++){
        if(!A.push_Cell(CELL_array + i));
            break;
        
        FreeCellList.push(CELL_array + i);
    }

    for(; i <= C; i++){
        if(!B.push_Cell(CELL_array + i));
            break;

        FreeCellList.push(CELL_array + i);
    }

    if(i == C)
        printf("Error on BlockInitialization");
}

//implement how to choose the base cell, find base cell, remove it from block and push it into FreeCellList
Cell* ChooseBaseCell(Block &A, Block &B, int r){ //r is a balance factor
    Cell* max_gain_cellA, * max_gain_cellB;
    max_gain_cellA = A.get_max_gain_cell();
    max_gain_cellB = B.get_max_gain_cell();

    if(max_gain_cellA == nullptr && max_gain_cellB == nullptr)
        return nullptr;
    else if(max_gain_cellA == nullptr){
        FreeCellList.push(max_gain_cellB);
        B.remove_from_BUCKET(max_gain_cellB);
        return max_gain_cellB;
    }
    else if(max_gain_cellB == nullptr){
        FreeCellList.push(max_gain_cellA);
        A.remove_from_BUCKET(max_gain_cellA);
        return max_gain_cellA;
    }
    else{
        if(std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) < std::abs(B.get_modified_balance_factor(max_gain_cellB) - r)){
            FreeCellList.push(max_gain_cellA);
            A.remove_from_BUCKET(max_gain_cellA);
            return max_gain_cellA;
        }
        else{
            FreeCellList.push(max_gain_cellB);
            B.remove_from_BUCKET(max_gain_cellB);
            return max_gain_cellB;
        }
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