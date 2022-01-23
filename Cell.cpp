#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"

void Cell::print_Cell(){
    std::cout << "name: " << name << ", size: " << size <<", pin_num: " << pin
                  << ", gain: " << gain << std::endl;
    
    std::cout << "Cell " << name << " is included in ";
    for(auto head =net_list.begin(); head != net_list.end(); head++){
        std::cout << "net " << (*head)->get_name() <<", ";
    }
    std::cout << std::endl;
}