#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "Block.h"

void Cell::print_Cell(){
    std::cout << "cell_num" << cell_num << ", name: " << name << ", size: " << size <<", pin_num: " << pin
                  << std::endl;
    
    std::cout << "Cell " << name << " is included in ";
    for(auto head =net_list.begin(); head != net_list.end(); head++){
        std::cout << "net " << (*head)->get_net_num() <<", ";
    }
    std::cout <<"Cell " << name << " is included in Block " << CurrentBlock->get_block_name() <<std::endl;
    std::cout << std::endl <<std::endl;
}

void Cell::get_net_list(int *distribution){
    for(auto head = net_list.begin(); head != net_list.end(); head++){
        distribution[(*head)->get_net_num()]++;
    }
}

void printCellInfo(Cell* CELL_array, int C){
    printf("Cell Info\n");
    for(int i = 1; i <= C; i++){
        printf("Cell num: %d, ", i);
        CELL_array[i].print_Cell();
    }
}