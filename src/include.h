/*
 Copyright 2017 D-Wave Systems Inc.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

// ---------misc stuff used by habit
#ifndef QBSOLV_INCLUDE_H
#define QBSOLV_INCLUDE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#define  VERSION    "open source 2.4 (no Tabu search)"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define EPSILON 1.0e-7
#define BIGNEGFP -1.225E308
#define flt_equals(a, b) (fabs((a) - (b)) < EPSILON)

#define GETMEM(P, T, N) (P = (T*)malloc(sizeof(T) * N))
#define BADMALLOC {printf("\n  ------> Exit(%d) called by %s(%s.%d)\n\n", 9, __func__, __FILE__, __LINE__); exit(9); }
#define DL printf("-----> AT %s(%s.%d)\n",  __func__, __FILE__, __LINE__);
#define CPSECONDS ((double)(clock() - start_) / CLOCKS_PER_SEC)
#define DLT printf("%lf seconds ", CPSECONDS);
#define uint unsigned int


// ----------------------- STRUCTs -------------------------------------
struct nodeStr_ {
    int32_t n1, n2;
    double value;
};

// ------------------- Prototypes --------------------------------------

#ifdef __cplusplus
extern "C" {
#endif


int   main( int argc,  char *argv[]);
int   read_qubo(const char *inFileName, FILE *inFile);
void  write_qubo(double **val, int nMax, const char *filename);
void  solve( double **val,  int maxNodes, int nRepeats);
void  **malloc2D(uint rows, uint cols, uint size  );
void  fill_qubo(double **qubo, int maxNodes, struct nodeStr_ *nodes, int nNodes, struct nodeStr_ *couplers, int nCouplers);
void  print_qubo_format( void);
void  print_help( void);
void  print_solution_and_qubo(int8_t *solution, int maxNodes, double **qubo);
void  randomize_solution(int8_t *solution, int nbits);
void  shuffle_index(int *index, int n);
void  print_opts( );
void  print_output(int maxNodes, int8_t *Q, long numPartCalls, double energy, double seconds);
void  quick_sort_iterative_index(double arr[], int ind[], int n, int stack[]);
void  val_index_sort(int *index, double *val, int n);
void  index_sort(int *index, int n, short FWD);
bool dw_established();
void dw_init( );
void dw_solver( double **val, int maxNodes, int8_t *Q );
void dw_close();
void reduce(int *Icompress, double **qubo, uint sub_qubo_size, uint qubo_size, double **sub_qubo, int8_t *solution, int8_t *sub_solution);
#ifdef __cplusplus
}
#endif
#endif
