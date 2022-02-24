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

void FM_pass(int C, int N, double r, int pass_num, Cell* CELL_array, Net* NET_array, Block &A, Block &B, const bool stuck, CellDist& LocalMinDist, bool big_wave);
int FM(const int InitVer, const int pass, Cell* _CELL_array, Net* _NET_array, const int _C, const int _N, const int _P, const int _W, const int block_num, const Block* current_block, int bias);
int calculate_degree(Cell* &CELL_array, int C, Net* NET_array, int N, int block_num, int cutnet);

#endif