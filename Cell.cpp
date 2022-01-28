#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"

void Cell::print_Cell(){
    std::cout << "cell_num" << cell_num << "name: " << name << ", size: " << size <<", pin_num: " << pin
                  << std::endl;
    
    std::cout << "Cell " << name << " is included in ";
    for(auto head =net_list.begin(); head != net_list.end(); head++){
        std::cout << "net " << (*head)->get_net_num() <<", ";
    }
    std::cout << std::endl;
}

void Cell::get_net_list(int *distribution){
    for(auto head = net_list.begin(); head != net_list.end(); head++){
        distribution[(*head)->get_net_num()]++;
    }
}