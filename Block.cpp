#define NDEBUG

#include <iostream>
#include <ctime>
#include <cassert>
#include <fstream>
#include <cmath>
#include <string>
#include <list>
#include <stack>
#include <queue>


#include "Cell.h"
#include "Net.h"
#include "Block.h"

std::stack<Cell*> FreeCellList;

extern int gain_update_count;
double get_max_gain_cell_time = 0;


//Using CellNode to record the best distribution
Block::Block(const Block& copy)
    : PMAX(copy.PMAX), max_gain(copy.max_gain), lbound(copy.lbound), ubound(copy.ubound), C(copy.C), N(copy.N), W(copy.W), R(copy.R), size(copy.size), name(copy.name) {
    BUCKET = new Cell*[2 * PMAX + 1];
    
    BUCKET += PMAX;

    for(int i = -PMAX; i <= PMAX; i++)
        BUCKET[i] = copy.BUCKET[i];

    Fdistribution = new int[N + 1];
    Ldistribution = new int[N + 1];
    gain = new int[C + 1];

    for(int i = 1; i <= C; i++){
        Fdistribution[i] = copy.Fdistribution[i];
        Ldistribution[i] = copy.Ldistribution[i];
        gain[i] = copy.gain[i];
    }
}

bool Block::operator==(const Block& compare) const{
    if(max_gain != compare.max_gain){
        printf("not match max gain\n");
        printf("look!: %d, %d\n", max_gain, compare.max_gain);
        return false;
    }
    
    for(int i = -PMAX; i <= PMAX; i++)
        if(BUCKET[i] != compare.BUCKET[i]){
            printf("not match BUCKET\n");
            return false;
        }
    
    for(int i =1; i <= C; i++){
        if(Fdistribution[i] != compare.Fdistribution[i]){
            printf("not match Fdis\n");
            return false;
        }
        if(Ldistribution[i] != compare.Ldistribution[i]){
            printf("not match Ldis\n");
            return false;
        }
        
        if(gain[i] != compare.gain[i]){
            printf("not match gain\n");
            return false;
        }
    }

    return true;
}


void Block::CalculateDistribution(Cell* CELL_array){
    //printf("CalDist start\n");

    for(int i = 1; i <= N; i++){
        Fdistribution[i] = 0;
        Ldistribution[i] = 0;
    }

    for(int i = 1; i <= C; i++){
        if(this == CELL_array[i].get_current_block()){
            CELL_array[i].get_net_list(Fdistribution);
        }
    }

    //printf("CalDist end\n");
}

void Block::deactivate_large_net(Net* NET_array){
    for(int i = 1; i <= N; i++){
        if(!NET_array[i].get_activate())
            Ldistribution[i] = 1;
    }
}

void Block::CellGainInitialization(Block &T, Cell &c){ //inner loop of implementation of the code prior to Proposition 2
    //printf("CellGainAdjust start\n");
    
    int gain_adjust = 0;
    if(this == c.get_current_block()){ //if this block is F (from block)
        for(auto i = c.net_list.begin(); i != c.net_list.end(); i++){
            if(Fdistribution[(*i)->get_net_num()] == 1)
                gain_adjust++;
        }

        for(auto i = c.net_list.begin(); i != c.net_list.end(); i++){ //if this block is T (to block)
            if((T.Fdistribution[(*i)->get_net_num()] + T.Ldistribution[(*i)->get_net_num()]) == 0)
                gain_adjust--;
        }

        //printf("c.cell_num: %d, C: %d, length of array gain: %d\n", c.cell_num, C, sizeof(gain));
        gain[c.get_cell_num()] = gain_adjust;
        c.BUCKETpre = nullptr;
        c.BUCKETnext = BUCKET[gain_adjust];

        if(BUCKET[gain_adjust] != nullptr && c.get_cell_num() == BUCKET[gain_adjust]->get_cell_num()){
            printf("cellgaininitialization\n");
            printf("%d, %d\n", c.get_cell_num(), BUCKET[gain_adjust]->get_cell_num());
            assert(c.get_cell_num() != BUCKET[gain_adjust]->get_cell_num());
        }

        if(BUCKET[gain_adjust] != nullptr)
            BUCKET[gain_adjust]->BUCKETpre = &c;
        
        BUCKET[gain_adjust] = &c;

        if(gain_adjust >= max_gain)
            max_gain = gain_adjust;
    }
    //printf("CellGainAdjust end\n");
}

Cell* Block::get_max_gain_cell() const{
    //printf("start get_max_gain_cell()\n");

    if(max_gain < -PMAX)
        return nullptr;

    int i = max_gain;
    int modified_size;
    Cell* max_gain_cell = BUCKET[max_gain];

    if(max_gain_cell == nullptr)
        printf("gain adjustment error\n");

    do{
        modified_size = size - max_gain_cell->get_size();

        //printf("modified_size: %d, lbound: %f, ubound: %f\n", modified_size, lbound, ubound);

        if(lbound <= modified_size && modified_size <= ubound)
            break;
        
        max_gain_cell = max_gain_cell ->BUCKETnext;

        if(max_gain_cell == nullptr){
            //return nullptr;
            
            i--;
            
            if(i < -PMAX)
                return nullptr;

            if(BUCKET[i] == nullptr)
                break;
            else
                max_gain_cell = BUCKET[i];
            
        }
    }while(BUCKET[i] != nullptr);

    //printf("end get_max_gain_cell()\n");
    return max_gain_cell;
}

void Block::remove_from_BUCKET(Cell* cell){
    //printf("start remove_from_BUCKET\n");

    if(this != cell->get_current_block())
        printf("error on remove_from_BUCKET() of Block Class\n");

    if(cell->BUCKETpre == nullptr){
        //printf("no BUCKETpre\n");
        BUCKET[gain[cell->get_cell_num()]] = cell->BUCKETnext;
    }
    else{
        //printf("BUCKETpre exists\n");
        cell->BUCKETpre->BUCKETnext = cell->BUCKETnext;
    }

    if(cell->BUCKETnext != nullptr){
        //printf("BUCKETnext exists\n");
        cell->BUCKETnext->BUCKETpre = cell->BUCKETpre;
    }
    cell->BUCKETpre = nullptr;
    cell->BUCKETnext = nullptr;
    cell->locked = true;

    adjust_maxgain();

    //printf("end remove_from_BUCKET\n");
}

void Block::adjust_maxgain(){
    if(max_gain > PMAX)
        max_gain = PMAX;

    while(max_gain >= -PMAX && BUCKET[max_gain] == nullptr){
        max_gain--;
    }
}

bool Block::push_Cell_ub(Cell* cell){ //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환    
    //printf("push_Cell_ub\n");
    //printf("pushed cell num: %d ", cell->cell_num);
    //std::cout << name <<std::endl;
    
    //printf("%d %d\n", size, cell->get_size());
    //printf("%f", ubound);

    if(size + cell->get_size() <= ubound){
        //printf("push_Cell_ub: return true\n");
        cell->set_current_block(this);
        size += cell->get_size();

        return true;
    }
    else{
        return false;
    }
}

bool Block::push_Cell_r(Cell* cell){ //cell을 추가하였을 때 size가 ubound를 넘지 않으면 push하고 true를 반환, 아니면 false 반환    
    //printf("push_Cell_ub\n");
    //printf("pushed cell num: %d ", cell->cell_num);
    //std::cout << name <<std::endl;
    
    //printf("%d %d\n", size, cell->get_size());
    //printf("%f", ubound);

    if(size + cell->get_size() <= R * W){
        //printf("push_Cell_ub: return true\n");
        cell->set_current_block(this);
        size += cell->get_size();

        return true;
    }
    else{
        return false;
    }
}

void Block::increase_cell_gain(Cell* cell){
    //printf("start increase_cell_gain\n");
    //printf("update Cell %d gain!\n", cell->get_cell_num());
    gain_update_count++;

    if(cell->BUCKETpre == nullptr){
        //printf("no BUCKETpre\n");
        BUCKET[gain[cell->get_cell_num()]] = cell->BUCKETnext;
    }
    else{
        //printf("BUCKETpre exists\n");
        cell->BUCKETpre->BUCKETnext = cell->BUCKETnext;
    }

    if(cell->BUCKETnext != nullptr){
        //printf("BUCKETnext exists\n");
        cell->BUCKETnext->BUCKETpre = cell->BUCKETpre;
    }

    cell->BUCKETpre = nullptr;
    cell->BUCKETnext = nullptr;

    gain[cell->get_cell_num()]++; 
    cell->BUCKETnext = BUCKET[gain[cell->get_cell_num()]];
    
    if(BUCKET[gain[cell->get_cell_num()]] != nullptr)
        BUCKET[gain[cell->get_cell_num()]]->BUCKETpre = cell;

    BUCKET[gain[cell->get_cell_num()]] = cell;

    if(max_gain < PMAX)
        max_gain++;
    
    adjust_maxgain();
    //printf("end increase_cell_gain\n");
}

void Block::decrease_cell_gain(Cell* cell){
    //printf("start decrease_cell_gain\n");
    gain_update_count++;

    if(cell->BUCKETpre == nullptr){
        //printf("no BUCKETpre\n");
        BUCKET[gain[cell->get_cell_num()]] = cell->BUCKETnext;
    }
    else{
        //printf("BUCKETpre exists\n");
        cell->BUCKETpre->BUCKETnext = cell->BUCKETnext;
    }

    if(cell->BUCKETnext != nullptr){
        //printf("BUCKETnext exists\n");
        cell->BUCKETnext->BUCKETpre = cell->BUCKETpre;
    }
    cell->BUCKETpre = nullptr;
    cell->BUCKETnext = nullptr;
    gain[cell->get_cell_num()]--; 
    cell->BUCKETnext = BUCKET[gain[cell->get_cell_num()]];
    
    if(BUCKET[gain[cell->get_cell_num()]] != nullptr)
        BUCKET[gain[cell->get_cell_num()]]->BUCKETpre = cell;

    BUCKET[gain[cell->get_cell_num()]] = cell;
    adjust_maxgain();

    //printf("end decrease_cell_gain\n");
}

void Block::increase_cell_gain_of_net(Net* net){
    for(auto i = net->cell_list.begin(); i != net->cell_list.end(); i++)
        if((*i)->locked == false)
            increase_cell_gain(*i);
        else
            gain_update_count++;
    
}

void Block::decrease_cell_gain_of_net(Net* net){
    for(auto i = net->cell_list.begin(); i != net->cell_list.end(); i++)
        if((*i)->locked == false)
            decrease_cell_gain(*i);
        else
            gain_update_count++;
}

Cell* Block::find_cell_in_block(Net* net){
    for(auto i = net->cell_list.begin(); i != net->cell_list.end(); i++){
        if((*i)->get_current_block() == this)
            return (*i);
    }

    printf("NO Cell Error in Block::find_cell_in_block\n");
    return nullptr;
}

int Block::get_cell_num(Cell* Cell_array, int C){
    int cell_num = 0;

    for(int i = 1; i <= C; i++)
        if(Cell_array[i].get_current_block() == this)
            cell_num++;
        
    return cell_num;
}

void Block::print_Block(Cell* CELL_array){
    printf("------------------------------------------\n");
    printf("max_gain: %d, lbound: %f , ubound: %f, size: %d\n", max_gain, lbound, ubound, size);
    
    printf("------------------------------------------\n");
    printf("Fdistribution & Ldistribution\n");
    printf("Free: ");
    for(int i = 1; i <= N; i++)
        printf("net %d: %d ", i, Fdistribution[i]);
    printf("\n");
    printf("Lock: ");
    for(int i = 1; i <= N; i++)
        printf("net %d: %d ", i, Ldistribution[i]);
    printf("\n");

    printf("------------------------------------------\n");
    printf("gain of cell\n");
    for(int i = 1; i <= C; i++)
        if(CELL_array[i].get_current_block() == this)
            printf("cell %d: %d ", i, gain[i]);
    printf("\n");

    printf("------------------------------------------\n");
    printf("BUCKET\n");
    
    int head = max_gain;
    Cell *cell_head = nullptr;

    while(head >= -PMAX){
        cell_head = BUCKET[head];

        if(BUCKET[head] != nullptr)
            printf("gain %d: ", head);

        while(cell_head != nullptr){
            printf("%d ", cell_head->get_cell_num());
            cell_head = cell_head->BUCKETnext;
        }
        if(BUCKET[head] != nullptr)
            printf("\n");

        head--;
    }

    printf("\n");
}

void Block::print_Block_short(Cell* CELL_array){

    printf("------------------------------------------\n");
    printf("gain of cell\n");
    for(int i = 1; i <= C; i++)
        if(CELL_array[i].get_current_block() == this)
            printf("cell %d: %d ", i, gain[i]);
    printf("\n");

    printf("------------------------------------------\n");
    printf("BUCKET\n");
    
    int head = max_gain;
    Cell *cell_head = nullptr;

    while(head >= -PMAX){
        cell_head = BUCKET[head];

        if(BUCKET[head] != nullptr)
            printf("gain %d: ", head);

        while(cell_head != nullptr){
            printf("%d ", cell_head->get_cell_num());
            cell_head = cell_head->BUCKETnext;
        }
        if(BUCKET[head] != nullptr)
            printf("\n");

        head--;
    }

    printf("\n");
}

void Block::empty_BUCKET(){
    max_gain = -PMAX;
    for(int i = -PMAX; i <= PMAX; i++)
        BUCKET[i] = nullptr;
}

//VERSION 1 ver 1
//block의 사이즈도 여기서 계산해주어야 한다. BlockInitialization 실행 후 Reinitialization도 실행시켜주어야..
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, int C){
    int i;

    for(i = 1; i <= C; i++){
        if(!A.push_Cell_ub(CELL_array + i))
            break;
        //printf("push in A\n");
        //FreeCellList.push(CELL_array + i);
    }

    printf("%d th cell is on block B\n", i);

    for(; i <= C; i++){
        if(!B.push_Cell_ub(CELL_array + i))
            break;
        //printf("push in B\n");
        //FreeCellList.push(CELL_array + i);
    }
    //printf("end\n");

    if(i == C)
        printf("Error on BlockInitialization");
}

//ver2
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N){
    std::queue<Net *> net_queue;
    Net* temp_net = nullptr;
    bool check = true;
    //bool alter = true;
    //Block* current_block = &A;

    
    //Net_Num_Size* num_size_order = new Net_Num_Size[N + 1];
    bool* NET_check_array = new bool[N + 1];

    for(int i = 1; i <= N; i++){
        NET_check_array[i] = true;
        //num_size_order[i].num = NET_array[i].get_net_num();
        //num_size_order[i].size = NET_array[i].get_size();
    }

    //qsort(num_size_order + 1, N, sizeof(Net_Num_Size), net_compare);

    net_queue.push(NET_array + 1);
    NET_check_array[1] = false;

    while(check){
        if(net_queue.empty()){
            for(int i = 1; i <= N; i++){
                if(NET_check_array[i]){
                    net_queue.push(NET_array + i);
                    NET_check_array[i] = false;

                    break;
                }

            }
        }
        temp_net = net_queue.front();
        net_queue.pop();

        for(auto itr = temp_net->cell_list.begin(); itr != temp_net->cell_list.end(); itr++){
            if((*itr)->get_current_block() == &A)
                continue;

            if(A.push_Cell_ub(*itr)){
                //FreeCellList.push(*itr);
                for(auto jtr = (*itr)->net_list.begin(); jtr != (*itr)->net_list.end(); jtr++){
                    if(NET_check_array[(*jtr)->get_net_num()]){
                        net_queue.push(*jtr);
                        NET_check_array[(*jtr)->get_net_num()] = false;
                    }
                }
            }
            else{
                check = false;
                break; // A가 가득 찼으므로 break, 나머지는 B에 채우자                
            }
        }

        
        /*
        if(check){
            for(int i = 1; i <= N; i++)
                if(NET_check_array[i]){
                    net_queue.push(NET_array + i);
                    NET_check_array[i] = false;
                }
        } 
        */       
    }

    int i;

    for(i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() != &A){
            if(!B.push_Cell_ub(CELL_array + i))
                break;
            
            //FreeCellList.push(CELL_array + i);
        }
    }

    if(i <= C)
        printf("Initialization error!\n");

    delete[] NET_check_array;
}


//BlockInitialization ver3
void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, int ver3){
    //std::queue<Net *> net_queue;
    Net* temp_net = nullptr;
    bool check = true;
    bool alter = true; //true면 A, false면 B
    Block* current_block = &A;

    bool* NET_check_array = new bool[N + 1]; //Block의 위치(A, B)를 배당받지 않으면 true, 아니면 false

    for(int i = 1; i <= N; i++)
        NET_check_array[i] = true;

    check = false;

    for(int i = 1; i <= N; i++){
        for(auto itr = (NET_array + i)->cell_list.begin(); itr != (NET_array + i)->cell_list.end(); itr++){
            if((*itr)->get_current_block() == &A)
                continue;

            if(!A.push_Cell_ub(*itr)){
                check = true;
                break;
            }
        }

        if(check)
            break;
    }

    int i;

    for(i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() != &A){
            if(!B.push_Cell_ub(CELL_array + i))
                break;

        }
    }
    /*
    while(check){
        if(net_queue.empty())
            break;
        temp_net = net_queue.front();
        net_queue.pop();

        for(auto itr = temp_net->cell_list.begin(); itr != temp_net->cell_list.end(); itr++){
            if((*itr)->get_current_block() != nullptr)
                continue;

            if(current_block->push_Cell_ub(*itr)){
                //FreeCellList.push(*itr);
                for(auto jtr = (*itr)->net_list.begin(); jtr != (*itr)->net_list.end(); jtr++){
                    if(NET_check_array[(*jtr)->get_net_num()]){
                        net_queue.push(*jtr);
                        NET_check_array[(*jtr)->get_net_num()] = false;
                    }
                }
            }
            else{
                check = false;
                alter = !alter;
                if(alter)
                    current_block = &A;
                else   
                    current_block = &B;     
                break; // A가 가득 찼으므로 break, 나머지는 B에 채우자                
            }

            alter = !alter;
            if(alter)
                current_block = &A;
            else   
                current_block = &B;

        }

        
        if(check){
            for(int i = 1; i <= N; i++)
                if(NET_check_array[i]){
                    net_queue.push(NET_array + i);
                    NET_check_array[i] = false;
                }
        }
        
    }

    int i;

    for(i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() == nullptr){
            if(!current_block->push_Cell_ub(CELL_array + i))
                break;
            
            //FreeCellList.push(CELL_array + i);
        }
    }
    */


    if(i <= C)
        printf("Initialization error! i: %d, C: %d\n", i , C);

    delete[] NET_check_array;
}

struct Net_Num_Size{
    int num;
    int size;
};

int net_compare(const void* _a, const void* _b){
    Net_Num_Size* a = (Net_Num_Size*)_a;
    Net_Num_Size* b = (Net_Num_Size*)_b;
    if(a->size < b->size)
        return -1;
    else if(a->size > b->size)
        return 1;
    return 0;
}

void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, int ver4, int _ver4){
    std::queue<Net *> net_queue;
    Net* temp_net = nullptr;
    bool check = true;
    //int net_push_finish = 1;
    //bool alter = true;
    //Block* current_block = &A;

    
    Net_Num_Size* num_size_order = new Net_Num_Size[N + 1];
    bool* NET_check_array = new bool[N + 1];

    for(int i = 1; i <= N; i++){
        NET_check_array[i] = true;
        num_size_order[i].num = NET_array[i].get_net_num();
        num_size_order[i].size = NET_array[i].get_cell_count();
    }

    qsort(num_size_order + 1, N, sizeof(Net_Num_Size), net_compare);

    for(int i = 1; i <= N; i++){
        temp_net = NET_array + num_size_order[i].num;
        NET_check_array[num_size_order[i].num] = false;

        for(auto itr = temp_net->cell_list.begin(); itr != temp_net->cell_list.end(); itr++){
            if((*itr)->get_current_block() == &A)
                continue;

            if(!A.push_Cell_ub(*itr))
                break;
        }
    }

    int i;

    for(i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() != &A){
            if(!B.push_Cell_ub(CELL_array + i))
                break;
        }
    }

    if(i <= C)
        printf("Initialization error!, C: %d, i:%d\n", C, i);

    delete[] num_size_order;        
    delete[] NET_check_array;
}

void BlockInitialization(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, int ver5, int _ver5, int __ver5){
    std::queue<Net *> net_queueA;
    std::queue<Net *> net_queueB;
    Net* temp_netA = nullptr;
    Net* temp_netB = nullptr;
    bool checkA = true;
    bool checkB = true;
    //bool alter = true;
    //Block* current_block = &A;

    
    //Net_Num_Size* num_size_order = new Net_Num_Size[N + 1];
    bool* NET_check_array = new bool[N + 1];

    for(int i = 1; i <= N; i++){
        NET_check_array[i] = true;
        //num_size_order[i].num = NET_array[i].get_net_num();
        //num_size_order[i].size = NET_array[i].get_size();
    }

    //qsort(num_size_order + 1, N, sizeof(Net_Num_Size), net_compare);

    net_queueA.push(NET_array + 1);
    NET_check_array[1] = false;

    net_queueB.push(NET_array + N);
    NET_check_array[N] = false;

    while(checkA && checkB){
        if(net_queueA.empty()){
            for(int i = 1; i <= N; i++){
                if(NET_check_array[i]){
                    net_queueA.push(NET_array + i);
                    NET_check_array[i] = false;

                    break;
                }

            }
        }

        if(net_queueB.empty()){
            for(int i = 1; i <= N; i++){
                if(NET_check_array[i]){
                    net_queueB.push(NET_array + i);
                    NET_check_array[i] = false;

                    break;
                }

            }
        }

        if(net_queueA.empty()){
            break;
        }
        if(net_queueB.empty()){
            break;
        }

        temp_netA = net_queueA.front();
        temp_netB = net_queueB.front();
        net_queueA.pop();
        net_queueB.pop();

        for(auto itr = temp_netA->cell_list.begin(); itr != temp_netA->cell_list.end(); itr++){
            if((*itr)->get_current_block() != nullptr)
                continue;

            if(A.push_Cell_ub(*itr)){
                //FreeCellList.push(*itr);
                for(auto jtr = (*itr)->net_list.begin(); jtr != (*itr)->net_list.end(); jtr++){
                    if(NET_check_array[(*jtr)->get_net_num()]){
                        net_queueA.push(*jtr);
                        NET_check_array[(*jtr)->get_net_num()] = false;
                    }
                }
            }
            else{
                checkA = false;
                break; // A가 가득 찼으므로 break, 나머지는 B에 채우자                
            }
        }

        for(auto itr = temp_netB->cell_list.begin(); itr != temp_netB->cell_list.end(); itr++){
            if((*itr)->get_current_block() != nullptr)
                continue;

            if(B.push_Cell_ub(*itr)){
                //FreeCellList.push(*itr);
                for(auto jtr = (*itr)->net_list.begin(); jtr != (*itr)->net_list.end(); jtr++){
                    if(NET_check_array[(*jtr)->get_net_num()]){
                        net_queueB.push(*jtr);
                        NET_check_array[(*jtr)->get_net_num()] = false;
                    }
                }
            }
            else{
                checkB = false;
                break; // B가 가득 찼으므로 break, 나머지는 A에 채우자                
            }
        }
    }

    int i;

    if(!checkA){
        for(i = 1; i <= C; i++){
            if(CELL_array[i].get_current_block() == nullptr){
                if(!B.push_Cell_ub(CELL_array + i))
                    break;
                
                //FreeCellList.push(CELL_array + i);
            }
        }
    }
    else{
        for(i = 1; i <= C; i++){
            if(CELL_array[i].get_current_block() == nullptr){
                if(!A.push_Cell_ub(CELL_array + i))
                    break;
                
                //FreeCellList.push(CELL_array + i);
            }
        }
        
    }

    if(i <= C)
        printf("Initialization error!, i: %d, C: %d\n", i, C);

    delete[] NET_check_array;
}


//implementation of the code prior to Proposition 2
//first empty_BUCEKT()
void BlockReinitialization(int C, Block &A, Block &B, Cell* CELL_array, Net* NET_array, int pass){
    //printf("BlockReinit start\n");

    A.empty_BUCKET();
    B.empty_BUCKET();


    A.CalculateDistribution(CELL_array);
    B.CalculateDistribution(CELL_array);

    Cell *ctemp;
    int bias = (C / 16) * (3 * (pass / 16) % 16);
    bool reverse = ((pass / 47) % 2 == 0);
    //int bias = (C / 4) * ((pass / 8) % 4);
    //bool reverse = ((pass / 64) % 2 == 0);
    //int bias = (C / 4) * ((pass / 37) % 4);
    //bool reverse = ((pass / 64) % 2 == 0);
    //bool reverse = false;
    
    if(reverse){
        for(int i = bias + 1; i <= C; i++){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }

        for(int i = 1; i <= bias; i++){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }
    }
    else{
        for(int i = C; i >= bias + 1; i--){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }

        for(int i = bias; i >= 1; i--){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }
    }
    
    /*
    if(pass % 64 < 20){
        for(int i = 1; i <= C; i++){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }
    }
    else if(pass % 64 < 40){
        for(int i = C / 2; i <= C; i++){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }

        for(int i = C / 2 - 1; i >= 1; i--){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        }
    }
    else{
       for(int i = C; i >= 1; i--){
            ctemp = CELL_array + i;

            ctemp->locked = false;

            A.CellGainInitialization(B, *ctemp);
            B.CellGainInitialization(A, *ctemp);
        } 
    }
    */

    A.adjust_maxgain();
    B.adjust_maxgain();

    //printf("BlockReinit end\n");
}

//implement how to choose the base cell, find base cell, remove it from bucket and push it into FreeCellList
//아직 block에서 제거가 되지는 않았다.
Cell* ChooseBaseCell_gain(Block &A, Block &B, double r){ //r is a balance factor
    //printf("ChooseBaseCell start!\n");

    clock_t get_max_gain_cell_time_temp = clock();
    Cell* max_gain_cellA, * max_gain_cellB;
    max_gain_cellA = A.get_max_gain_cell();
    max_gain_cellB = B.get_max_gain_cell();

    get_max_gain_cell_time += (clock() - (double)get_max_gain_cell_time_temp);


    if(max_gain_cellA == nullptr && max_gain_cellB == nullptr)
        return nullptr;
    else if(max_gain_cellA == nullptr){
        //FreeCellList.push(max_gain_cellB);
        B.remove_from_BUCKET(max_gain_cellB);
        //printf("choose cell in B1\n");
        return max_gain_cellB;
    }
    else if(max_gain_cellB == nullptr){
        //FreeCellList.push(max_gain_cellA);
        A.remove_from_BUCKET(max_gain_cellA);
        //printf("choose cell in A1\n");
        return max_gain_cellA;
    }
    else{
        /*
        printf("max_gain_cell A is...\n");
        max_gain_cellA->print_Cell();

        printf("max_gain_cell B is...\n");
        max_gain_cellB->print_Cell();
        */

        if(A.gain[max_gain_cellA->get_cell_num()] > B.gain[max_gain_cellB->get_cell_num()]){
            //FreeCellList.push(max_gain_cellA);
            A.remove_from_BUCKET(max_gain_cellA);
            //printf("choose cell in A3\n");
            return max_gain_cellA;
        }
        else if(A.gain[max_gain_cellA->get_cell_num()] < B.gain[max_gain_cellB->get_cell_num()]){
            //FreeCellList.push(max_gain_cellB);
            B.remove_from_BUCKET(max_gain_cellB);
            //printf("choose cell in B3\n");
            return max_gain_cellB;
        }
            
    

        //std::cout << std::endl << "A: " << std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) << ", B: " << std::abs(1 - r - B.get_modified_balance_factor(max_gain_cellB)) <<std::endl;
        if(std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) < std::abs(1 - r - B.get_modified_balance_factor(max_gain_cellB))){
            //FreeCellList.push(max_gain_cellA);
            A.remove_from_BUCKET(max_gain_cellA);
            //printf("choose cell in A2\n");
            return max_gain_cellA;
        }
        else{
            //FreeCellList.push(max_gain_cellB);
            B.remove_from_BUCKET(max_gain_cellB);
            //printf("choose cell in B2\n");
            return max_gain_cellB;
        }
    }
}

Cell* ChooseBaseCell_balance(Block &A, Block &B, double r, bool destroy_balance){ //r is a balance factor
    //printf("ChooseBaseCell start!\n");

    Cell* max_gain_cellA, * max_gain_cellB;
    max_gain_cellA = A.get_max_gain_cell();
    max_gain_cellB = B.get_max_gain_cell();

    if(max_gain_cellA == nullptr && max_gain_cellB == nullptr)
        return nullptr;
    else if(max_gain_cellA == nullptr){
        //FreeCellList.push(max_gain_cellB);
        B.remove_from_BUCKET(max_gain_cellB);
        //printf("choose cell in B1\n");
        return max_gain_cellB;
    }
    else if(max_gain_cellB == nullptr){
        //FreeCellList.push(max_gain_cellA);
        A.remove_from_BUCKET(max_gain_cellA);
        //printf("choose cell in A1\n");
        return max_gain_cellA;
    }
    else{
        //std::cout << std::endl << "A: " << std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) << ", B: " << std::abs(1 - r - B.get_modified_balance_factor(max_gain_cellB)) <<std::endl;
        if(destroy_balance){
            if(std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) > std::abs(1 - r - B.get_modified_balance_factor(max_gain_cellB))){
                //push(max_gain_cellA);
                A.remove_from_BUCKET(max_gain_cellA);
                //printf("choose cell in A2\n");
                return max_gain_cellA;
            }
            else{
                //FreeCellList.push(max_gain_cellB);
                B.remove_from_BUCKET(max_gain_cellB);
                //printf("choose cell in B2\n");
                return max_gain_cellB;
            }
        }
        else{
            if(std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) < std::abs(1 - r - B.get_modified_balance_factor(max_gain_cellB))){
                //FreeCellList.push(max_gain_cellA);
                A.remove_from_BUCKET(max_gain_cellA);
                //printf("choose cell in A2\n");
                return max_gain_cellA;
            }
            else{
                //FreeCellList.push(max_gain_cellB);
                B.remove_from_BUCKET(max_gain_cellB);
                //printf("choose cell in B2\n");
                return max_gain_cellB;
            }
        }
    }
}

//BlockA가 크면 which_block은 true;
Cell* ChooseBaseCell_balance(Block &A, Block &B, double r, bool destroy_balance, bool which_block){ //r is a balance factor
    //printf("ChooseBaseCell start!\n");

    Cell* max_gain_cellA, * max_gain_cellB;
    max_gain_cellA = A.get_max_gain_cell();
    max_gain_cellB = B.get_max_gain_cell();

    if(max_gain_cellA == nullptr && max_gain_cellB == nullptr)
        return nullptr;
    else if(max_gain_cellA == nullptr){
        //FreeCellList.push(max_gain_cellB);
        B.remove_from_BUCKET(max_gain_cellB);
        //printf("choose cell in B1\n");
        return max_gain_cellB;
    }
    else if(max_gain_cellB == nullptr){
        //FreeCellList.push(max_gain_cellA);
        A.remove_from_BUCKET(max_gain_cellA);
        //printf("choose cell in A1\n");
        return max_gain_cellA;
    }
    else{
        //std::cout << std::endl << "A: " << std::abs(A.get_modified_balance_factor(max_gain_cellA) - r) << ", B: " << std::abs(1 - r - B.get_modified_balance_factor(max_gain_cellB)) <<std::endl;
        if(destroy_balance){
            if(which_block){
                //push(max_gain_cellA);
                A.remove_from_BUCKET(max_gain_cellA);
                //printf("choose cell in A2\n");
                return max_gain_cellA;
            }
            else{
                //FreeCellList.push(max_gain_cellB);
                B.remove_from_BUCKET(max_gain_cellB);
                //printf("choose cell in B2\n");
                return max_gain_cellB;
            }
        }
        else{
            if(which_block){
                //FreeCellList.push(max_gain_cellA);
                A.remove_from_BUCKET(max_gain_cellA);
                //printf("choose cell in A2\n");
                return max_gain_cellA;
            }
            else{
                //FreeCellList.push(max_gain_cellB);
                B.remove_from_BUCKET(max_gain_cellB);
                //printf("choose cell in B2\n");
                return max_gain_cellB;
            }
        }
    }
}

//실제 Move는 MoveCell에서 이루어짐
void MoveCell(Block &F, Block &T, Cell* BaseCell){
    //printf("start MoveCell\n");
    for(auto i = BaseCell->net_list.begin(); i != BaseCell->net_list.end(); i++){
        if(T.Ldistribution[(*i)->get_net_num()] == 0){
            if(T.Fdistribution[(*i)->get_net_num()] == 0){
                //printf("MoveCell 1\n");
                F.increase_cell_gain_of_net(*i);
            }
            else if(T.Fdistribution[(*i)->get_net_num()] == 1){
                //printf("MoveCell 2\n");
                T.decrease_cell_gain(T.find_cell_in_block(*i));
            }
        }

        F.Fdistribution[(*i)->get_net_num()]--;
        T.Ldistribution[(*i)->get_net_num()]++;
    }

    BaseCell->set_current_block(&T);
    F.add_size(-BaseCell->get_size());
    T.add_size(BaseCell->get_size());

    for(auto i = BaseCell->net_list.begin(); i != BaseCell->net_list.end(); i++){
        if(F.Ldistribution[(*i)->get_net_num()] == 0){
            if(F.Fdistribution[(*i)->get_net_num()] == 0){
                //printf("MoveCell 3\n");  
                T.decrease_cell_gain_of_net(*i);
            }
            else if(F.Fdistribution[(*i)->get_net_num()] == 1){
                //printf("MoveCell 4\n");
                F.increase_cell_gain(F.find_cell_in_block(*i));
            }
        }
    }

    //printf("end MoveCell\n");
}

//bigger가 true면 A가 더 크다.
void StuckOut(Block &A, Block &B, Cell* CELL_array, Net* NET_array, int C, int N, bool bigger){
    std::queue<Net *> net_queue;
    Net* temp_net = nullptr;
    bool check = true;
    //bool alter = true;
    Block* current_block = (bigger)? &B : &A;

    bool* NET_check_array = new bool[N + 1];

    for(int i = 1; i <= N; i++)
        NET_check_array[i] = true;

    for(int i = 1; i <= C; i++){
        if(CELL_array[i].get_current_block() == current_block){
            net_queue.push(CELL_array[i].net_list.front());
            NET_check_array[CELL_array[i].net_list.front()->get_net_num()] = false;
            break;
        }
    }

    while(check){
        temp_net = net_queue.front();
        net_queue.pop();

        for(auto itr = temp_net->cell_list.begin(); itr != temp_net->cell_list.end(); itr++){
            if((*itr)->get_current_block() == current_block)
                continue;

            if(current_block->push_Cell_ub(*itr)){
                for(auto jtr = (*itr)->net_list.begin(); jtr != (*itr)->net_list.end(); jtr++){
                    if(NET_check_array[(*jtr)->get_net_num()]){
                        net_queue.push(*jtr);
                        NET_check_array[(*jtr)->get_net_num()] = false;
                    }
                }
            }
            else{
                check = false;
                break; // A가 가득 찼으므로 break, 나머지는 B에 채우자                
            }
        }

        if(check){
            for(int i = 1; i <= N; i++)
                if(NET_check_array[i]){
                    net_queue.push(NET_array + i);
                    NET_check_array[i] = false;
                }
        }
    }

    delete[] NET_check_array;
}