#include <iostream>
#include <list>

#include "Cell.h"
#include "Net.h"
#include "Block.h"

void get_max(const int C, const Cell* CELL_array, int &pmax, int &smax){
    pmax = -1;
    smax = -1;

    pmax = CELL_array[1].get_pin();
    smax = CELL_array[1].get_size();

    for(int i = 2; i <= C; i++){    
        if(pmax < CELL_array[i].get_pin())
            pmax = CELL_array[i].get_pin();

        if(smax < CELL_array[i].get_size())
            smax = CELL_array[i].get_size();
    }
}

int CountCutNet(Block &block, Net *NET_array, int N){
    int cutnet = 0;
    
    for(int i = 1; i <= N; i++){
        if((block.ith_net_distribution(i) < NET_array[i].cell_list.size()) && (block.ith_net_distribution(i) != 0))
            cutnet++;
        else if(block.ith_net_distribution(i) > NET_array[i].cell_list.size())
            printf("error\n");
    }

    return cutnet;
}

void CountCutNetAgain(Block &block, Net *NET_array, int N){
    
    for(int i = 1; i <= N; i++){
        if(block.ith_net_distribution(i) < NET_array[i].cell_list.size())
            printf("Net %d\n", i);
        else if(block.ith_net_distribution(i) > NET_array[i].cell_list.size())
            printf("error\n");
    }
}

void Check(Block &A, Block &B, Net* NET_array, int N){
    for(int i = 1; i <= N; i++)
        if(A.ith_net_distribution(i) + B.ith_net_distribution(i) != NET_array[i].cell_list.size())
            printf("%d th net error\n", i);
}