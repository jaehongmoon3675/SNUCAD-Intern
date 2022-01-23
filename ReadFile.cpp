#include <iostream>
#include <fstream>
#include <string>

#include "Net.h"
#include "Cell.h"

void read_hgr(const int N, const int C, Net* &NET_array, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("initPlace.def.hgr");

    if(readFile.is_open()){
        char c;
        int temp_cell;

        scanf("%d %d", &N, &C);
        readFile.get(c);

        NET_array = new Net[N + 1];
        CELL_array = new Cell[C + 1];

        for(int i = 1; i <= N; i++){
            do{
                scanf("%d", &temp_cell);
                readFile.get(c);

                NET_array[i].push_cell(CELL_array + temp_cell);
                CELL_array[temp_cell].push_net(NET_array + i);
            }while(c != '\n');
        }

        readFile.close();
    }
    else
        printf("error for reading hgr file\n");
}

void read_hgr_map(const int C, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("initPlace.def.hgr.map");

    int temp_cell;
    std::string temp_cell_name;

    if(readFile.is_open()){
        for(int i = 1; i < C; i++){
            std::cin >> temp_cell >> temp_cell_name;

            CELL_array[temp_cell].set_name(temp_cell_name);
        }

        readFile.close();
    }
    else{
        printf("No map file\n");
    }
}

void read_hgr_area(const int C, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("initPlace.def.hgr.area");

    int temp_cell_size;

    if(readFile.is_open()){
        for(int i = 1; i < C; i++){
            scanf("%d", &temp_cell_size);

            CELL_array[i].set_size(temp_cell_size);
        }

        readFile.close();
    }
    else{
        printf("No area file\n");
    }
}