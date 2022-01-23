#include <iostream>

#include "ReadFile.h"
#include "Cell.h"
#include "Net.h"

#include <list>

class Block{

};


int main(){
    int N, C; //the num of net and cell, respectively
    Net* NET_array = nullptr;
    Cell* CELL_array = nullptr;

    read_hgr(N, C, NET_array, CELL_array);
    read_hgr_map(C, CELL_array);
    read_hgr_area(C, CELL_array);

}