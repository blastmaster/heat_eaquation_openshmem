#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <mpp/shmem.h>
#include "heat2d.h"
#include "sharedblk.h"

#define CLOCK_PRECISION 1000000000L

long pSync[_SHMEM_BCAST_SYNC_SIZE];
double pWrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
double residual = 0.0;
double du = 0.0;
double *t = NULL;
uint32_t step = 0;
double n = 10.0;
double m = 10.0;
unsigned int colums;
unsigned int rows;

// iterators
int i = 0;
int j = 0;

int main (int argc, char *argv[])
{
  double ln = 2.0;
  // double lm = ln*m/n;
  double h = ln/n;
  double hi = 1.0/h;
  double hi2 = pow (hi,2);
  double alpha = 22;
  double dt = pow (h,2)/4/alpha;
  double eps = 1.0e-4;
  double *a = NULL;
  pos_desc *mp = NULL;

  // for measurements
  struct timespec start_time, end_time, start_step, end_step;
  //struct timespec res;

  // shmem variables
  int me, npes = 0;

  for (; i < _SHMEM_BCAST_SYNC_SIZE; i++) {
    pSync[i] = _SHMEM_SYNC_VALUE;
  }

  start_pes(0);
#ifdef OSH
  me = _my_pe();
  npes = _num_pes();
#else
  me = shmem_my_pe();
  npes = shmem_n_pes();
#endif

  //clock_getres(CLOCK_MONOTONIC, &res);

  if ( npes > 2 ) {
    colums = (int) sqrt ((double) npes);
    rows = colums;
  }
  else if ( npes == 2) {
    colums = COLUMS;
    rows = ROWS;
  }
  else {
    fprintf(stderr,"Error,Can not build block in these configuration.");
    return EXIT_FAILURE;
  }

  /*
   * eventually colum and row
   * distribution as commandline argument
   */

  if ( argc > 2 ) {
    n = atof (argv[1]);
    m = atof (argv[2]);
  }

  //setting global variables from header
  x_size = colums;
  y_size = rows;

  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);

  allocate_all(&mp, &a);
  register_pe (mp, me);
  initialize_matrix(&mp, &a);
  residual = eps + 1.0;
  shmem_barrier_all();

  clock_gettime(CLOCK_MONOTONIC, &start_time);
  while (residual > eps)
  {
      residual = 0.0;
      for (i = 1; i <= n_blks; i++) {
          for (j = 1; j <= m_blks; j++) {
              du =
                  (a[i * (m_blks+1) + (j-1)] * mp->u[i * (m_blks+2) + (j-1)] + a[i * (m_blks+1) + j] * mp->u[i * (m_blks+2) + (j+1)] +
                   a[(i-1) * (m_blks+1) + j] * mp->u[(i-1) * (m_blks+2) + j] + a[i * (m_blks+1) + j] * mp->u[(i+1) * (m_blks+2) + j] -
                   (a[i * (m_blks+1) + (j-1)] + a[i * (m_blks+1) + j] + a[(i-1) * (m_blks+1) + j] + a[i * (m_blks+1) + j]) *
                   mp->u[i * (m_blks + 2) + j]) * dt * hi2;

              mp->v[i * (m_blks + 2) + j] = mp->u[i * (m_blks + 2) + j] + du;
              du = MAX(du, -du);
              residual = MAX(residual, du);
          }
      }

      t = mp->u;
      mp->u = mp->v;
      mp->v = t;

      shmem_barrier_all();

      for ( i = 0; i < MAX_NEIGHBORS; i++ ) {
          if ( mp->neighbors[i] != -1 )
              get_border_from_neighbor(mp, i, mp->neighbors[i]);
      }

      step++;
      shmem_double_max_to_all(&residual, &residual, 1, 0, 0, npes, pWrk, pSync);
  }
  shmem_barrier_all();
  clock_gettime(CLOCK_MONOTONIC, &end_time);

  if (me == 0) {
      printf("%d %14.8lf \n", npes,
              ((double) end_time.tv_sec - (double) start_time.tv_sec) + ((double) end_time.tv_nsec - (double) start_time.tv_nsec) / CLOCK_PRECISION);
#ifdef STEP
      printf("steps %d \n", step);
#endif
  }

  my_free(mp, a);
  // sync before exit
  shmem_barrier_all();
#ifndef OSH
  shmem_finalize();
#endif
  return EXIT_SUCCESS;
}

/*
 * In case of it is a north neighbor the pe get the southside border
 * if it is a south neighbor get the northside border
 * if it is a west neighbor get the eastside border
 * if it is a east border get the westside border
 */

void get_border_from_neighbor (pos_desc *mp, int idx, int remote_pe)
{
  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);

  switch (idx) {
    case NORTH: shmem_double_get (&(mp->u[1]), &(mp->u[n_blks * (m_blks + 2)+1]), m_blks, remote_pe);
		break;
    case SOUTH: shmem_double_get (&(mp->u[(n_blks + 1) * (m_blks + 2) + 1]), &(mp->u[(m_blks + 2) + 1]), m_blks, remote_pe);
		break;
    case WEST: shmem_double_iget (&(mp->u[m_blks + 2]), &(mp->u[(2 * m_blks)+2]), m_blks + 2, m_blks + 2, n_blks, remote_pe);
	       break;
    case EAST: shmem_double_iget (&(mp->u[2*(m_blks + 1)+1]), &(mp->u[(m_blks + 2) + 1]), m_blks + 2, m_blks + 2, n_blks, remote_pe);
	       break;
  }

  return;
}

void print_square(pos_desc *mp)
{
  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);

  for (i = 0; i <= n_blks + 1; i++) {
    for (j = 0; j <= m_blks + 1; j++) {
      printf(" %lf ", mp->u[i * (m_blks + 2) + j]);
    }
    printf("\n\n");
  }

  return;
}

void print_square_all(pos_desc *mp)
{
#ifdef OSH
    int me = _my_pe();
    int npes = _num_pes();
#else
    int me = shmem_my_pe();
    int npes = shmem_n_pes();
#endif
  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);

  for ( i = 0; i < npes; i++ ) {
    if ( me == i ) {
      printf(" me = %d\n", me);
      for ( i = 0; i <= n_blks + 1; i++ ) {
          for ( j = 0; j <= m_blks + 1; j++ ) {
              printf(" %lf ", mp->u[i * (m_blks + 2) + j]);
          }
          printf("\n\n");
      }
    }
    shmem_barrier_all();
  }

  shmem_barrier_all();
  return;
}

/*
 * Should print an explanation of the
 * commandline arguments and there meaning.
 */
void print_usage (void)
{
  /**
   * its a trap!
   */
  return;
}

void allocate_all (pos_desc **mp, double **a)
{
  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);

  *mp = (pos_desc *) shmalloc (sizeof(pos_desc));
  assert (NULL != *mp);

  (*mp)->u = (double *) shmalloc ( (n_blks + 2) * (m_blks + 2) * sizeof(double));
  assert (NULL != (*mp)->u);

  (*mp)->v = (double *) shmalloc ( (n_blks + 2) * (m_blks + 2) * sizeof(double));
  assert (NULL != (*mp)->v);

  *a = (double *) shmalloc ((n_blks + 1) * (m_blks + 1) * sizeof(double));
  assert (NULL != *a);

  return;
}

int *get_total_position (pos_desc *mp, int col, int row)
{
  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);
  static int res[] = {0, 0};

  res[0] = mp->xpos * m_blks + row;
  res[1] = mp->ypos * n_blks + col;

  return res;
}

void initialize_matrix (pos_desc **mp, double **a)
{
  double u_min = 10.0;
  double u_max = 100.0;
  double alpha = 22;
  //double alpha_low = 0.035;
  int n_blks = BLKS(n, colums);
  int m_blks = BLKS(m, rows);
  int *total_pos = NULL;

  for ( i = 0; i <= n_blks + 1; i++ ) {
      for ( j = 0; j <= m_blks + 1; j++ ) {
          total_pos = get_total_position(*mp, i, j);
          if ( total_pos[1] < n / 2 ) {
              (*mp)->u[i * (m_blks + 2) + j] = u_max;
              (*mp)->v[i * (m_blks + 2) + j] = u_max;
          }
          else {
              (*mp)->u[i * (m_blks + 2) + j] = u_min;
              (*mp)->v[i * (m_blks + 2) + j] = u_min;
          }
      }
  }
  for (i = 0; i <= n_blks; i++) {
      for (j = 0; j <= m_blks; j++) {
          (*a)[i * (m_blks + 1) +j] = alpha;
      }
  }

  /*
  for (i = 0; i < 2*n_blks/4; i++) {
    for (j = 1*m_blks/2; j < 4*m_blks/7; j++) {
      (*a)[i * (m_blks + 1) + j] = alpha_low;
    }
  }
  */

  return;
}

void my_free (pos_desc *mp, double *a)
{
  shfree(mp->u);
  mp->u = NULL;
  shfree(mp->v);
  mp->v = NULL;
  shfree(mp);
  mp = NULL;

  shfree(a);
  a = NULL;

  return;
}
