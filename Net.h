#ifndef __NET_H__
#define __NET_H__

class Cell;  //in ReadFile.cpp

class Net{
public:
    std::list<Cell*> cell_list;

    Net() : net_num(0), size(0), cut_state(false), cell_list(std::list<Cell*>()) {}
    void push_cell(Cell* c){
        cell_list.push_back(c);
        size += c->get_size();
    }
    int get_net_num(){
        return net_num;
    }
    void set_net_num(int num) { net_num = num; };
    ~Net(){

    }
private:
    int net_num;
    int size;
    bool cut_state;
};

#endif