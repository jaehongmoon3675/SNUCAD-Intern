#ifndef __READFILE_H__
#define __READFILE_H__

#include "Cell.h"
#include "Net.h"

//read_hgr returns P which is the sum of all cell's # of pin
int read_hgr(int &N, int &C, Net* &NET_array, Cell* &CELL_array);
void read_hgr_map(const int C, Cell* &CELL_array);
int read_hgr_area(const int C, Cell* &CELL_array);

#endif