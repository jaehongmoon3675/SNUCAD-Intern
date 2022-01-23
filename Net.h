#ifndef __Net_h__
#define __Net_h__

class Net{
public:
    Net() : cut_state(false), cell_list(std::list<Cell*>()) {}
    void push_cell(Cell* c){
        cell_list.push_front(c);
    }
    ~Net(){

    }
private:
    bool cut_state;
    std::list<Cell*> cell_list;
};

#endif