#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "sharedblk.h"

int x_size = 0;
int y_size = 0;

int get_neighbors(unsigned int xpos, unsigned int ypos)
{
    unsigned int neighbor_count = MAX_NEIGHBORS;

    if (xpos == 0 || xpos == (x_size - 1))
        neighbor_count--;
    if (ypos == 0 || ypos == (y_size - 1))
        neighbor_count--;

    return neighbor_count;
}

bool have_south_neighbor (unsigned int ypos)
{
    if (ypos < y_size - 1)
        return true;
    else 
        return false;
}

bool have_north_neighbor (unsigned int ypos)
{
    if (ypos == 0)
        return false;
    else
        return true;
}

bool have_east_neighbor (unsigned int xpos)
{
    if (xpos < x_size - 1)
        return true;
    else
        return false;
}

bool have_west_neighbor (unsigned int xpos)
{
    if (xpos == 0)
        return false;
    else 
        return true;
}

void register_pe (pos_desc *mypos, int pe_number)
{
    mypos->xpos = pe_number % x_size;
    mypos->ypos = pe_number / x_size;
    mypos->neighbor_count = get_neighbors(mypos->xpos, mypos->ypos);
    memset(mypos->neighbors, -1, MAX_NEIGHBORS * sizeof(int));

    if (have_north_neighbor(mypos->ypos))
        mypos->neighbors[NORTH] = pe_number - x_size;
    if (have_south_neighbor(mypos->ypos))
        mypos->neighbors[SOUTH] = pe_number + x_size;
    if (have_east_neighbor(mypos->xpos))
        mypos->neighbors[EAST] = pe_number + 1;
    if (have_west_neighbor(mypos->xpos)) 
        mypos->neighbors[WEST] = pe_number - 1;

    return;
}

/* EOF */
