#ifndef __CELL_H__
#define __CELL_H__

class Net;
class Block;

class Cell{
public:
    std::list<Net *> net_list;
    bool locked;
    Cell* BUCKETpre, * BUCKETnext;

    Cell() : size(1), pin(0), cell_num(0), current_block_num(0), BUCKETpre(nullptr), BUCKETnext(nullptr), CurrentBlock(nullptr), locked(false), net_list(std::list<Net*>()) {}
    void push_net(Net *n){
        net_list.push_back(n);
        pin++;
    }
    void set_cell_num(int num){
        cell_num = num;
        name = std::to_string(num);
    }
    int get_pin() const{
        return pin;
    }
    void set_current_block(Block* block){
        CurrentBlock = block;
    }
    Block* get_current_block() const{
        return CurrentBlock;
    }
    int get_cell_num() const{
        return cell_num;
    }
    int get_size() const{
        return size;
    }
    int get_current_block_num() const{ return current_block_num; }
    void set_size(int size){
        this->size = size;
    }
    void set_name(std::string name){
        this->name = name;
    }
    void set_current_block_num(int _current_block_num) { current_block_num = _current_block_num; }
    std::string get_cell_name() const { return name; }
    void print_Cell();
    void get_net_list(int* distribution); //add the distribution of net list to the array 'distribution'
    ~Cell(){
    }
private:
    int size;
    int pin;
    int cell_num;
    int current_block_num;
    Block* CurrentBlock;
    std::string name;
};

void printCellInfo(Cell* CELL_array, int C);

#endif