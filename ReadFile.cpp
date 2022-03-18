#define NDEBUG
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>

#include "Cell.h"
#include "Net.h"
#include "Block.h"

int read_hgr(int &N, int &C, Net* &NET_array, Cell* &CELL_array, std::string _filename, int x_length, int y_length){
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

        NET_array = new Net[N + 1 + x_length * y_length];
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

        for(int i = N + 1; i <= N + x_length * y_length; i++){
            NET_array[i].set_net_num(i);
            NET_array[i].overlap_net = true;
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

void read_place(const int C, Cell* CELL_array, std::string _filename, int map_n, int map_m, std::vector<int> *BIN_array){
    std::ifstream readFile;
    std::string filename = _filename + ".def.place";
    readFile.open(filename);

    int ll_x, ll_y, ur_x, ur_y;

    //printf("read_hgr\n");

    if(readFile.is_open()){
        char c;
        int temp_cell;
        int cell_x, cell_y;
        int cell_bin;
        std::string cell_name;

        readFile >> ll_x >> ll_y >> ur_x >> ur_y;
        readFile.get(c);

        double bin_x_length = (double)(ur_x - ll_x) / map_m;
        double bin_y_length = (double)(ur_y - ll_y) / map_n;

        for(int i = 1; i <= C; i++){         
                readFile >> temp_cell >> cell_name >> cell_x >> cell_y;

                cell_bin = (int)((cell_y - ll_y) / bin_y_length) * map_m + (int)((cell_x - ll_x) / bin_x_length);

                assert(cell_bin < map_n * map_m);

                CELL_array[i].set_bin(cell_bin);
                BIN_array[cell_bin].push_back(i);

                if(readFile.eof())
                    break;
                
                readFile.get(c);

                if(readFile.eof())
                    break;
        }

        readFile.close();
    }
    else
        printf("error for reading place file\n");
}

void read_place(const int C, Cell* CELL_array, std::string _filename, int map_n, int map_m, std::vector<int> *BIN_array, int &_ll_x, int &_ll_y, int &_ur_x, int &_ur_y){
    std::ifstream readFile;
    std::string filename = _filename + ".def.scaled0.707.place";
    readFile.open(filename);

    double ll_x, ll_y, ur_x, ur_y;

    //printf("read_hgr\n");

    if(readFile.is_open()){
        char c;
        int temp_cell;
        double cell_x, cell_y;
        int cell_bin;
        std::string cell_name;

        readFile >> ll_x >> ll_y >> ur_x >> ur_y;
        _ll_x = ll_x;
        _ll_y = ll_y;
        _ur_x = ur_x;
        _ur_y = ur_y;

        if(CELL_array == nullptr)
            return;

        readFile.get(c);

        double bin_x_length = (double)(ur_x - ll_x) / map_m;
        double bin_y_length = (double)(ur_y - ll_y) / map_n;

        for(int i = 1; i <= C; i++){         
                readFile >> temp_cell >> cell_name >> cell_x >> cell_y;

                cell_bin = (int)((cell_y - ll_y) / bin_y_length) * map_m + (int)((cell_x - ll_x) / bin_x_length);

                assert(cell_bin < map_n * map_m);

                CELL_array[i].set_bin(cell_bin);
                CELL_array[i].set_ll(cell_x, cell_y);
                BIN_array[cell_bin].push_back(i);

                if(readFile.eof())
                    break;
                
                readFile.get(c);

                if(readFile.eof())
                    break;
        }
        /*
        printf("bin test\n");
        for(int i = 1; i <= C; i++){
            if(CELL_array[i].get_bin() == -1)
                printf("error1\n");
        }
        for(int i = 0; i < map_n * map_m; i++)
            printf("Bin %d size: %ld\n", i, BIN_array[i].size());
        */
        readFile.close();
    }
    else
        printf("error for reading place file\n");
}

void read_place(const int C, Cell* CELL_array, const int _N, const int N, Net* NET_array, std::string _filename, int map_n, int map_m, std::vector<int> *BIN_array, int &_ll_x, int &_ll_y, int &_ur_x, int &_ur_y, int overlap_x, int overlap_y){
    std::ifstream readFile;
    std::string filename = _filename + ".def.scaled0.707.place";
    readFile.open(filename);

    double ll_x, ll_y, ur_x, ur_y;
    int x_scale, y_scale;
    int tight = 10;

    //printf("read_hgr\n");

    if(readFile.is_open()){
        char c;
        int temp_cell;
        double cell_x, cell_y;
        int cell_bin;
        std::string cell_name;

        readFile >> ll_x >> ll_y >> ur_x >> ur_y;
        _ll_x = ll_x;
        _ll_y = ll_y;
        _ur_x = ur_x;
        _ur_y = ur_y;

        x_scale = (ur_x - ll_x) / (overlap_x + 1);
        y_scale = 1;

        if(CELL_array == nullptr)
            return;

        readFile.get(c);

        double bin_x_length = (double)(ur_x - ll_x) / map_m;
        double bin_y_length = (double)(ur_y - ll_y) / map_n;

        int overlap_criteria;

        for(int i = 1; i <= C; i++){         
                readFile >> temp_cell >> cell_name >> cell_x >> cell_y;

                cell_bin = (int)((cell_y - ll_y) / bin_y_length) * map_m + (int)((cell_x - ll_x) / bin_x_length);

                assert(cell_bin < map_n * map_m);

                overlap_criteria = (((cell_x + (double)x_scale / (tight + 5)) + (double)CELL_array[i].get_size() * 0.707 + 0.5) / x_scale);

                if(overlap_criteria > (((int)cell_x) / x_scale)){
                    NET_array[_N + (int)cell_y * overlap_x + overlap_criteria].push_cell(CELL_array + i);
                    CELL_array[i].push_net(NET_array + _N + (int)cell_y * overlap_x + overlap_criteria);
                }

                CELL_array[i].set_bin(cell_bin);
                CELL_array[i].set_ll(cell_x, cell_y);
                BIN_array[cell_bin].push_back(i);

                if(readFile.eof())
                    break;
                
                readFile.get(c);

                if(readFile.eof())
                    break;
        }
        /*
        printf("bin test\n");
        for(int i = 1; i <= C; i++){
            if(CELL_array[i].get_bin() == -1)
                printf("error1\n");
        }
        for(int i = 0; i < map_n * map_m; i++)
            printf("Bin %d size: %ld\n", i, BIN_array[i].size());
        */
        readFile.close();
    }
    else
        printf("error for reading place file\n");
}

void check_place(const int C, Cell* CELL_array, std::string _filename, int map_n, int map_m){
    std::ifstream readFile;
    std::string filename = _filename + ".def.place";
    readFile.open(filename);

    int ll_x, ll_y, ur_x, ur_y;

    //printf("read_hgr\n");

    if(readFile.is_open()){
        char c;
        int temp_cell;
        int cell_x, cell_y;
        int cell_bin;
        std::string cell_name;

        readFile >> ll_x >> ll_y >> ur_x >> ur_y;
        readFile.get(c);

        double bin_x_length = (double)(ur_x - ll_x) / map_m;
        double bin_y_length = (double)(ur_y - ll_y) / map_n;

        for(int i = 1; i <= C; i++){         
                readFile >> temp_cell >> cell_name >> cell_x >> cell_y;

                cell_bin = (int)((cell_y - ll_y) / bin_y_length) * map_m + (int)((cell_x - ll_x) / bin_x_length);

                assert(cell_bin < map_n * map_m);

                if(cell_bin != CELL_array[i].get_bin())
                    printf("bin error!\n");

                if(readFile.eof())
                    break;
                
                readFile.get(c);

                if(readFile.eof())
                    break;
        }

        readFile.close();
    }
    else
        printf("error for reading place file\n");
}


//cell의 순서가 .def.place와 완전히 동일하다고 가정
void read_partial_part(const int C, Cell* &CELL_array, std::string filename){
    std::ifstream readFile;
    readFile.open(filename +".partial.part");

    //printf("read_hgr_map\n");

    int block;
    std::string temp_cell_name;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_name >> block;
            
            while(CELL_array[i].get_cell_name() != temp_cell_name){
                i++;

                if(i > C)
                    break;
            }

            if(i > C)
                break;

            CELL_array[i].set_current_block_num(block);
            CELL_array[i].set_fixed(true);

            if(readFile.eof())
                break;
        }

        readFile.close();
    }
    else
        printf("No .partial.part file\n");
}

void check_partial_part(const int C, Cell* &CELL_array, std::string filename){
    std::ifstream readFile;
    readFile.open(filename +".partial.part");

    int block;
    std::string temp_cell_name;
    int count = 0;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_name >> block;
            
            while(CELL_array[i].get_cell_name() != temp_cell_name){
                i++;

                if(i > C)
                    break;
            }

            if(i > C)
                break;

            if(CELL_array[i].get_current_block_num() != block){
                printf("cell %d error! %d %d\n", i, CELL_array[i].get_current_block_num(), block);
                count++;
            }

            if(readFile.eof())
                break;
        }

        readFile.close();
    }
    else
        printf("No .partial.part file\n");

    printf("count: %d\n", count);
}

void read_output_part(Block &A, Block &B, const int C, Cell* &CELL_array){
    std::ifstream readFile;
    readFile.open("ldpc.part");

    //printf("read_hgr_map\n");

    int block;
    std::string temp_cell_name;

    if(readFile.is_open()){
        for(int i = 1; i <= C; i++){
            readFile >> temp_cell_name >> block;
            
            if(block == 1){
                CELL_array[i].set_current_block(&A);
                A.add_size(CELL_array[i].get_size());
            }
            else{
                CELL_array[i].set_current_block(&B);
                B.add_size(CELL_array[i].get_size());
            }
        }

        readFile.close();
    }
    else
        printf("No part file\n");
}