#include <iostream>
#include <list>

using namespace std;

class Block{

};

class Net{
public:
private:
    int gain;
    bool cut_state;
    list<Cell *> cell_array;
};

class Cell{
public:
private:
    int size;
    int pin;
    bool locked;

    list<Net *> net_array;
};

int main();