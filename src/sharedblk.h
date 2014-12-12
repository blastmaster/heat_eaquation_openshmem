#ifndef SHAREDBLK_H
#define SHAREDBLK_H

#include <stdbool.h>

#define MAX_NEIGHBORS 4

extern int x_size;
extern int y_size;

/*
 * These enumeration aliased the shared address indicies
 */
typedef enum eidx {NORTH, SOUTH, WEST, EAST} idx;

/*
 * position_description
 * shows how many neighbors exists.
 * Every pe got's his own position description.
 * the int *neighbors array contains the pe number of the according remote pe
 * can indicated with the eidx enum.
 */
typedef struct position_description
{
    double *v;
    double *u;
    unsigned int neighbor_count;
    unsigned int xpos; 
    unsigned int ypos;
    int neighbors[MAX_NEIGHBORS];
} pos_desc;

extern void register_pe (pos_desc *mypos, int penumber);
extern bool have_north_neighbor (unsigned int ypos);
extern bool have_south_neighbor (unsigned int ypos);
extern bool have_east_neighbor (unsigned int xpos);
extern bool have_west_neighbor (unsigned int xpos);

#endif /* SHAREDBLK_H */
