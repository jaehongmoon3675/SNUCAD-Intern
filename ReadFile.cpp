#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>

#include "Cell.h"
#include "Net.h"

int read_hgr(int &N, int &C, Net* &NET_array, Cell* &CELL_array, std::string _filename){
    int pin_num = 0;
    std::ifstream readFile;
    std::string filename = _filename + ".def.hgr";
    readFile.open(filename);

    //printf("read_hgr\n");

    if(readFile.is_open()){
        char c;
        int temp_cell;

        readFile >> N >> C;
        readFile.get(c);

        NET_array = new Net[N + 1];
        CELL_array = new Cell[C + 1];

        for(int i = 1; i <= C; i++)
            CELL_array[i].set_cell_num(i);

        for(int i = 1; i <= N; i++){
            std::vector<int> net_twice;
            do{ 
                
                readFile >> temp_cell;

                if(!CELL_array[temp_cell].net_list.empty() && CELL_array[temp_cell].net_list.back()->get_net_num() != i){
                    pin_num++;

                    NET_array[i].set_net_num(i);
                    NET_array[i].push_cell(CELL_array + temp_cell);
                    CELL_array[temp_cell].push_net(NET_array + i);
                }
                else if(CELL_array[temp_cell].net_list.empty()){
                    pin_num++;

                    NET_array[i].set_net_num(i);
                    NET_array[i].push_cell(CELL_array + temp_cell);
                    CELL_array[temp_cell].push_net(NET_array + i);
                }

                if(readFile.eof())
                    break;
                
                readFile.get(c);

                if(readFile.eof())
                    break;
            }while(readFile.peek() != '\n');
        }

        readFile.close();
    }
    else
        printf("error for reading hgr file\n");

    return pin_num;
}

void read_hgr_map(const int C, Cell* &CELL_array, std::string _filename){
    std::ifstream readFile;
    std::string filename = _filename + ".def.hgr.map";
    readFile.open(filename);

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

int read_hgr_area(const int C, Cell* &CELL_array, std::string _filename){
    int total_weight = 0;
    
    std::ifstream readFile;
    std::string filename = _filename + ".def.area";
    readFile.open(filename);

    //printf("read_hgr_area\n");

    int temp_cell_size;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_size;

            total_weight += temp_cell_size;
            CELL_array[i].set_size(temp_cell_size);
        }

        readFile.close();
    }
    else{
        printf("No area file\n");
        total_weight = C;
    }

    return total_weight;
}