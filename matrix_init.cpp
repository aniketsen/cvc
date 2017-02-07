/****************************************************
 * matrix_init.cpp
 *
 * Mon Oct 17 08:10:58 CEST 2016
 *
 * PURPOSE:
 * DONE:
 * TODO:
 ****************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef HAVE_MPI
#  include <mpi.h>
#endif
#ifdef HAVE_OPENMP
#  include <omp.h>
#endif
#include "matrix_init.h"

#ifndef EXIT
#ifdef HAVE_MPI
#define EXIT(_i) { MPI_Abort(MPI_COMM_WORLD, (_i)); MPI_Finalize(); exit((_i)); }
#else
#define EXIT(_i) { exit(_i); }
#endif
#endif


namespace cvc {

/************************************************************************************
 * (de-)allocate 2-level buffer (n1 x n2 double matrix)
 ************************************************************************************/
int init_2level_buffer (double***buffer, unsigned int n1, unsigned int n2) {

  unsigned int i;
  size_t items = (size_t)n2;
  size_t bytes = items * n1 * sizeof(double);

  if(*buffer != NULL) {
    fprintf(stderr, "[init_2level_buffer] Error, buffer not NULL\n"); 
    return(1);
  }

  /* 1st, outer level */
  (*buffer) = (double**)malloc(n1 * sizeof(double*));
  if( *buffer == NULL ) {
    fprintf(stderr, "[init_2level_buffer] Error from malloc\n"); 
    EXIT(2);
  }

  /* 2nd, inner level */
  (*buffer)[0] = (double*)malloc(bytes);
  if( (*buffer)[0] == NULL ) {
    fprintf(stderr, "[init_2level_buffer] Error from malloc\n"); 
    EXIT(3);
  }
  for(i=1; i<n1; i++) (*buffer)[i] = (*buffer)[i-1] + items;
  memset( (*buffer)[0], 0, bytes);
  return(0);
}  /* end of init_2level_buffer */

int fini_2level_buffer (double***buffer) {
  if(*buffer != NULL) {
    if( (*buffer)[0] != NULL) { free( (*buffer)[0] ); }
    free( *buffer );
    *buffer = NULL;
  }
  return(0);
}  /* end of fini_2level_buffer */

/************************************************************************************
 * (de-)allocate 3-level buffer (n1 x n2 x n3 matrix)
 ************************************************************************************/
int init_3level_buffer (double****buffer, unsigned int n1, unsigned int n2, unsigned int n3) {

  unsigned int i, k, l;
  size_t items=0, bytes=0;

  if(*buffer != NULL) {
    fprintf(stderr, "[init_3level_buffer] Error, buffer not NULL\n"); 
    return(1);
  }

  /* 1st, outermost level */
  (*buffer) = (double***)malloc(n1 * sizeof(double**));
  if( *buffer == NULL ) {
    fprintf(stderr, "[init_3level_buffer] Error from malloc\n"); 
    EXIT(2);
  }

  /* 2nd, middle level */
  items = (size_t)n2;
  bytes = items * n1 * sizeof(double*);
  (*buffer)[0] = (double**)malloc(bytes);
  if( (*buffer)[0] == NULL ) {
    fprintf(stderr, "[init_3level_buffer] Error from malloc\n"); 
    EXIT(3);
  }
  for(i=1; i<n1; i++) (*buffer)[i] = (*buffer)[i-1] + items;

  /* 3rd, innermost level */
  items = (size_t)n3 * (size_t)n2 * (size_t)n1;
  bytes = items * sizeof(double);
  (*buffer)[0][0] = (double*)malloc(bytes);
  if( (*buffer)[0][0] == NULL ) {
    fprintf(stderr, "[init_3level_buffer] Error from malloc\n"); 
    EXIT(4);
  }
  l=0;
  for(i=0; i<n1; i++) {
  for(k=0; k<n2; k++) {
    (*buffer)[i][k] = (*buffer)[0][0] + l*n3;
    l++;
  }}
  memset( (*buffer)[0][0], 0, bytes);
  return(0);
}  /* end of init_3level_buffer */

int fini_3level_buffer (double****buffer) {
  if(*buffer != NULL) {
    if( (*buffer)[0] != NULL) { 
      if( (*buffer)[0][0] != NULL) { free( (*buffer)[0][0] ); }
      free( (*buffer)[0] ); 
    }
    free( *buffer );
    *buffer = NULL;
  }
  return(0);
}  /* end of fini_3level_buffer */


/************************************************************************************
 * (de-)allocate 4-level buffer (n1 x n2 x n3 x n4 matrix)
 ************************************************************************************/
int init_4level_buffer (double*****buffer, unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4) {

  unsigned int i, k, l, m;
  size_t items=0, bytes=0;

  if(*buffer != NULL) {
    fprintf(stderr, "[init_4level_buffer] Error, buffer not NULL\n"); 
    return(1);
  }

  /* 1st, outermost level */
  (*buffer) = (double****)malloc(n1 * sizeof(double***));
  if( *buffer == NULL ) {
    fprintf(stderr, "[init_4level_buffer] Error from malloc\n"); 
    EXIT(2);
  }

  /* 2nd level */
  items = (size_t)n2;
  bytes = items * n1 * sizeof(double**);
  (*buffer)[0] = (double***)malloc(bytes);
  if( (*buffer)[0] == NULL ) {
    fprintf(stderr, "[init_4level_buffer] Error from malloc\n"); 
    EXIT(3);
  }
  for(i=1; i<n1; i++) (*buffer)[i] = (*buffer)[i-1] + items;

  /* 3rd, penultimate level */
  items = (size_t)n3 * (size_t)n2 * (size_t)n1;
  bytes = items * sizeof(double*);
  (*buffer)[0][0] = (double**)malloc(bytes);
  if( (*buffer)[0][0] == NULL ) {
    fprintf(stderr, "[init_4level_buffer] Error from malloc\n"); 
    EXIT(4);
  }
  l=0;
  for(i=0; i<n1; i++) {
  for(k=0; k<n2; k++) {
    (*buffer)[i][k] = (*buffer)[0][0] + l*n3;
    l++;
  }}

  /* 4th, innermost level */
  items = (size_t)n3 * (size_t)n2 * (size_t)n1 * (size_t)n4;
  bytes = items * sizeof(double);
  (*buffer)[0][0][0] = (double*)malloc(bytes);
  if( (*buffer)[0][0][0] == NULL ) {
    fprintf(stderr, "[init_4level_buffer] Error from malloc\n"); 
    EXIT(5);
  }
  m=0;
  for(i=0; i<n1; i++) {
  for(k=0; k<n2; k++) {
  for(l=0; l<n3; l++) {
    (*buffer)[i][k][l] = (*buffer)[0][0][0] + m*n4;
    m++;
  }}}

  memset( (*buffer)[0][0][0], 0, bytes);
  return(0);
}  /* end of init_4level_buffer */

int fini_4level_buffer (double*****buffer) {
  if(*buffer != NULL) {
    if( (*buffer)[0] != NULL) { 
      if( (*buffer)[0][0] != NULL) { 
        if( (*buffer)[0][0][0] != NULL) { free( (*buffer)[0][0][0] );} 
        free( (*buffer)[0][0] ); 
      }
      free( (*buffer)[0] ); 
    }
    free( *buffer );
    *buffer = NULL;
  }
  return(0);
}  /* end of fini_4level_buffer */

}  /* end of namespace cvc */