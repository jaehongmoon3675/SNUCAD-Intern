#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "Net.h"
#include "Cell.h"

int CellCount = 0, NetCount = 0;

int read_hgr(int &N, int &C, Net* &NET_array, Cell* &CELL_array){
    int pin_num = 0;
    std::ifstream readFile;
    readFile.open("initPlace.def.hgr");

    //printf("read_hgr\n");

    if(readFile.is_open()){
        char c;
        int temp_cell;

        readFile >> N >> C;
        readFile.get(c);

        NET_array = new Net[N + 1];
        CELL_array = new Cell[C + 1];

        for(int i = 1; i <= N; i++){
            do{ 
                readFile >> temp_cell;
                pin_num++;;

                NET_array[i].push_cell(CELL_array + temp_cell);
                CELL_array[temp_cell].push_net(NET_array + i);

                if(readFile.eof())
                    break;
                readFile.get(c);
                if(readFile.eof())
                    break;
            }while(c != '\n');
        }

        readFile.close();
    }
    else
        printf("error for reading hgr file\n");

    return pin_num;
}

void read_hgr_map(const int C, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("initPlace.def.hgr.map");

    //printf("read_hgr_map\n");

    int temp_cell;
    std::string temp_cell_name;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell >> temp_cell_name;

            CELL_array[temp_cell].set_name(temp_cell_name);
        }

        readFile.close();
    }
    else
        printf("No map file\n");
}

void read_hgr_area(const int C, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("initPlace.def.area");

    //printf("read_hgr_area\n");

    int temp_cell_size;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_size;

            CELL_array[i].set_size(temp_cell_size);
        }

        readFile.close();
    }
    else{
        printf("No area file\n");
    }
}