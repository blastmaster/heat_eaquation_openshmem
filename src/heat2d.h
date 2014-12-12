#ifndef HEAT2D_H
#define HEAT2D_H

#include "sharedblk.h"

/*
 * define macros
 */

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define BLKS(N,D) ((N)/(D))

#define COLUMS 2
#define ROWS 2

/*
 * declare functions
 */
extern void print_usage (void);
extern void allocate_all (pos_desc **mypos, double **a);
extern int *get_total_position (pos_desc *mp, int col, int row);
extern void initialize_matrix (pos_desc **mypos, double **a);
extern void my_free (pos_desc *mp, double *a);
extern void get_border_from_neighbor (pos_desc *mp, int idx, int remote_pe);
extern void print_square(pos_desc *mp);
extern void print_square_all(pos_desc *mp);
extern void print_alpha (double *a);

#endif /* HEAT2D_H */
