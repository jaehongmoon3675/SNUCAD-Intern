#ifndef __FM_FUNC_H__
#define __FM_FUNC_H__

#include "Cell.h"
#include "Net.h"
#include "Block.h"

//input C, Cell_array; output pmax, smax
void get_max(const int C, const Cell* CELL_array, int &pmax, int &smax);
int CountCutNet(Block &block, Net *NET_array, int N);
void Check(Block &A, Block &B, Net* NET_array, int N);
void CountCutNetAgain(Block &block, Net *NET_array, int N);
int get_max_cell_count(Net* NET_array, int N);

#endif