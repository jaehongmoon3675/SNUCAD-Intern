#ifndef __NET_H__
#define __NET_H__

class Cell;
extern int NetCount;

class Net{
public:
    Net() : name(NetCount++), cut_state(false), cell_list(std::list<Cell*>()) {}
    void push_cell(Cell* c){
        cell_list.push_front(c);
    }
    int get_name(){
        return name;
    }
    ~Net(){

    }
private:
    int name;
    bool cut_state;
    std::list<Cell*> cell_list;
};

#endif