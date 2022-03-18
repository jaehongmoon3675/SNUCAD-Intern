#ifndef __NET_H__
#define __NET_H__

class Cell;
class Block;

class Net{
public:
    std::list<Cell*> cell_list;

    Net() : net_num(0), size(0), cell_count(0), cut_state(false), activate(true), current_block(nullptr), cell_list(std::list<Cell*>()), weight(1), overlap(0) {}
    void push_cell(Cell* c){
        cell_list.push_back(c);
        size += c->get_size();
        cell_count++;
    }
    int get_net_num(){
        return net_num;
    }
    void print_Net(){
        printf("Net %d\t size: %d", net_num, size);

        printf("Net consists of ...\n");
        for(auto itr = cell_list.begin(); itr != cell_list.end(); itr++){
            printf("Cell %d ", (*itr)->get_cell_num());
        }
        printf("\n");
    }
    int get_size() const { return size; }
    int get_cell_count() const { return cell_count; }
    Block* get_current_block() const { return current_block; }
    void set_net_num(int num) { net_num = num; };
    void set_current_block(Block* _current_block) { current_block = _current_block; }
    bool get_activate() const { return activate; }
    int get_weight() const { return weight; }
    void set_weight(int _weight) { weight = _weight; }
    int adjust_overlap(bool** map, int ll_x, int ll_y);
    int count_overlap(bool** map, int ll_x, int ll_y, const int block_num) const;
    void adjust_weight(int max_overlap, int alpha_tune);
    int get_overlap() const { return overlap; }
    bool operator<(const Net& compare){
        if(this->size < compare.size)
            return true;
        
        return false;
    }
    ~Net(){

    }
private:
    int net_num;
    int size;
    int cell_count;
    bool cut_state;
    bool activate;
    int weight;
    int overlap;
    Block* current_block; //최종적으로 다음 FM으로 넘길때만 사용해야한다
};

void printNetInfo(Net* NET_array, int N);

#endif