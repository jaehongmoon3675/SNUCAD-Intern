#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "Block.h"

void write_output(Block &A, Cell* CEll_array, int C){
    std::string filePath = "output.part";

    std::ofstream writeFile(filePath.data());

    if(writeFile.is_open()){
        for(int i = 1; i <= C; i++){
            writeFile << CEll_array[i].get_cell_name() << " ";

            if(CEll_array[i].get_current_block() == &A)
                writeFile << "1" << std::endl;
            else
                writeFile << "0" << std::endl;
        }

        writeFile.close();
    }
}