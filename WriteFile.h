#ifndef __WRITEFILE_H__
#define __WRITEFILE_H__

#include "Cell.h"
#include "Net.h"
#include "Block.h"

void write_output(Cell* CEll_array, int C, std::string _filename, int block_num, int init_num, int pass);

#endif