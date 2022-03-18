#include <iostream>
#include <list>
#include <vector>
#include <cmath>

#include "Cell.h"
#include "Net.h"

extern int ALPHA;

void printNetInfo(Net* NET_array, int N){
    printf("Net Info\n");
    for(int i = 1; i <= N; i++){
        NET_array[i].print_Net();
    }
}

int Net::adjust_overlap(bool** map, int ll_x, int ll_y){
    int cell_ll_x, cell_ll_y, cell_ur_x, cell_ur_y;

    for(auto itr = cell_list.begin(); itr != cell_list.end(); itr++){
        cell_ll_x = (*itr)->get_ll_x() - ll_x;
        cell_ll_y = (*itr)->get_ll_y() - ll_y;
        cell_ur_x = (*itr)->get_ll_x() - ll_x + (*itr)->get_size();
        cell_ur_y = (*itr)->get_ll_y() - ll_y + 1;

        for(int i = cell_ll_x; i <= cell_ur_x; i++){
            for(int j = cell_ll_y; j <= cell_ur_y; j++){
                if(map[i][j])
                    overlap++;
                else
                    map[i][j] = true;
            }
        }
    }

    for(auto itr = cell_list.begin(); itr != cell_list.end(); itr++){
        cell_ll_x = (*itr)->get_ll_x() - ll_x;
        cell_ll_y = (*itr)->get_ll_y() - ll_y;
        cell_ur_x = (*itr)->get_ll_x() - ll_x + (*itr)->get_size();
        cell_ur_y = (*itr)->get_ll_y() - ll_y + 1;

        for(int i = cell_ll_x; i <= cell_ur_x; i++){
            for(int j = cell_ll_y; j <= cell_ur_y; j++){
                map[i][j] = false;
            }
        }
    }

    return overlap;
}

void Net::adjust_weight(int max_overlap, int alpha_tune){
    if(cell_count == 0)
        weight = 0;
    else if(overlap_net){
        weight = ((std::sqrt(max_overlap + 1 - overlap) + max_overlap / 2) * (ALPHA - alpha_tune));
    }
    else
        weight = alpha_tune / 10 + 1;

    /*
    if(overlap == 0){
        weight = (max_overlap * alpha_tune) / ALPHA;
    }
    else
        weight = 1 + (max_overlap - overlap);
    */
}

int Net::count_overlap(bool** map, int ll_x, int ll_y, const int block_num) const {
    int cell_ll_x, cell_ll_y, cell_ur_x, cell_ur_y;
    int overlap_count = 0;

    for(auto itr = cell_list.begin(); itr != cell_list.end(); itr++){
        cell_ll_x = (*itr)->get_ll_x() - ll_x;
        cell_ll_y = (*itr)->get_ll_y() - ll_y;
        cell_ur_x = (*itr)->get_ll_x() - ll_x + (*itr)->get_size();
        cell_ur_y = (*itr)->get_ll_y() - ll_y + 1;

        if((*itr)->get_current_block_num() != block_num)
            continue;

        for(int i = cell_ll_x; i <= cell_ur_x; i++){
            for(int j = cell_ll_y; j <= cell_ur_y; j++){
                if(map[i][j])
                    overlap_count++;
                else
                    map[i][j] = true;
            }
        }
    }

    for(auto itr = cell_list.begin(); itr != cell_list.end(); itr++){
        cell_ll_x = (*itr)->get_ll_x() - ll_x;
        cell_ll_y = (*itr)->get_ll_y() - ll_y;
        cell_ur_x = (*itr)->get_ll_x() - ll_x + (*itr)->get_size();
        cell_ur_y = (*itr)->get_ll_y() - ll_y + 1;

        if((*itr)->get_current_block_num() != block_num)
            continue;

        for(int i = cell_ll_x; i <= cell_ur_x; i++){
            for(int j = cell_ll_y; j <= cell_ur_y; j++){
                map[i][j] = false;
            }
        }
    }

    return overlap_count;
}