#ifndef __Cell_h__
#define __Cell_h__

class Cell{
public:
    static int CellCount;

    Cell() : size(1), pin(0), gain(0), name(std::to_string(CellCount++)), locked(false), net_list(std::list<Net*>()) {}
    void push_net(Net *n){
        net_list.push_front(n);
        pin++;
    }
    void set_size(int size){
        this->size = size;
    }
    void set_name(std::string name){
        this->name = name;
    }
    ~Cell(){

    }
private:
    int size;
    int pin;
    int gain;
    bool locked;
    std::string name;

    std::list<Net *> net_list;
};

int Cell::CellCount = 1;

#endif