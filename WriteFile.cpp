#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>

#include "Cell.h"
#include "Net.h"
#include "Block.h"

void write_output(Cell* CEll_array, int C, std::string _filename, int map_n, int map_m, int block_num, int init_num, int pass){
    //std::string filePath = _filename + "_" + std::to_string(block_num) + "_" + std::to_string(init_num) + "_" + std::to_string(pass) + ".part";

    std::string filePath = _filename + ".part";
    //std::string filePath = _filename + "_" + std::to_string(map_n) + "_" + std::to_string(map_m) + ".part";
    //std::string filePath = "output.part";
    
    std::ofstream writeFile(filePath.data());

    if(writeFile.is_open()){
        for(int i = 1; i <= C; i++){
            writeFile << CEll_array[i].get_cell_name() << " ";

            writeFile << CEll_array[i].get_current_block_num() << std::endl;
        }

        writeFile.close();
    }
}