#ifndef __FM_FUNC_H__
#define __FM_FUNC_H__

#include "Cell.h"
#include "Net.h"
#include "Block.h"
#include "CellDist.h"

//input C, Cell_array; output pmax, smax
void get_max(const int C, const Cell* CELL_array, int &pmax, int &smax);
int CountCutNet(Block &block, Net *NET_array, int N);
void Check(Block &A, Block &B, Net* NET_array, int N);
void CountCutNetAgain(Block &block, Net *NET_array, int N);
int get_max_cell_count(Net* NET_array, int N);

void FM_pass(int C, int N, double r, int pass_num, Cell* CELL_array, Net* NET_array, Block &A, Block &B, const bool stuck, CellDist& LocalMinDist, bool big_wave, bool alternate);
void FM_pass(int C, int N, double r, int pass_num, Cell* CELL_array, Net* NET_array, Block &A, Block &B, const bool stuck, CellDist& LocalMinDist, bool big_wave, bool alternate, int block_num);

//최초 호출시에는 current_block = nullptr
int FM(const int InitVer, const int pass, Cell* _CELL_array, Net* _NET_array, const int _C, const int _N, const int _P, const int _W, const int block_num, const Block* current_block, double skew, int bias, bool alternate, int bin_num = 0);
void bin_based_FM(const int InitVer, const int pass, Cell* _CELL_array, Net* _NET_array, const int _C, const int _N, const int _P, const int _W, const int block_num, double skew, int map_n, int map_m, std::vector<int> *BIN_array);
int calculate_degree(Cell* &CELL_array, int C, Net* NET_array, int N, int block_num, int &cutnet);

#endif