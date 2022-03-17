#ifndef __READFILE_H__
#define __READFILE_H__

//#include <vector>
#include "Cell.h"
#include "Net.h"
#include "Block.h"

//read_hgr returns P which is the sum of all cell's # of pin
int read_hgr(int &N, int &C, Net* &NET_array, Cell* &CELL_array, std::string _filename);
void read_hgr_map(const int C, Cell* &CELL_array, std::string _filename);
int read_hgr_area(const int C, Cell* &CELL_array, std::string _filename);
void read_place(const int C, Cell* CELL_array, std::string _filename, int map_n, int map_m, std::vector<int> *BIN_array);
void read_output_part(Block &A, Block &B, const int C, Cell* &CELL_array);
void read_partial_part(const int C, Cell* &CELL_array, std::string _filename);
bool check_partial_part(const int C, Cell* &CELL_array, std::string filename);
void check_place(const int C, Cell* CELL_array, std::string _filename, int map_n, int map_m);

#endif