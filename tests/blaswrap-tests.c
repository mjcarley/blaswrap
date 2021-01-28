/* This file is part of BLASWRAP, a set of C wrappers for BLAS routines
 *
 * Copyright (C) 2020 Michael Carley
 *
 * BLASWRAP is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  WBFMM is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with WBFMM.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#include <blaswrap.h>

gint random_matrix_d(gdouble *A, gint nr, gint nc) ;
gint matrix_transpose_d(gdouble *B, gdouble *A, gint nr, gint nc) ;

gint matrix_vector_mul_d(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble al, gdouble bt) ;
gint matrix_vector_mul_z(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble *al, gdouble *bt) ;
gint matrix_matrix_multiply_d(gdouble *A, gdouble *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gdouble al, gdouble bt,
			      gdouble *C, gint ldc) ;

gint matrix_multiply_test_d(gint nr, gint nc)

{
  gdouble B[2048], A[2048], x[128], y[128], z[128], al, bt, err ;
  gint i, incx, incy ;

  fprintf(stderr, "double-precision matrix-vector multiply\n") ;
  
  incx = 3 ; incy = 2 ;
  
  random_matrix_d(A, nr, nc) ;
  random_matrix_d(x, 1, nc*incx) ;
  random_matrix_d(y, nr*incy, 1) ;

  memcpy(z, y, nr*incy*sizeof(gdouble)) ;

  al = g_random_double() ; 
  bt = g_random_double() ; 

  matrix_vector_mul_d(A, nr, nc, x, incx, y, incy, al, bt) ;

  blaswrap_dgemv(FALSE, nr, nc, al, A, nc, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "plain error: %lg\n", err) ;

  /*set up for transpose multiplication*/
  memcpy(z, x, nc*incx*sizeof(gdouble)) ;

  blaswrap_dgemv(TRUE, nr, nc, al, A, nc, y, incy, bt, z, incx) ;

  matrix_transpose_d(B, A, nr, nc) ;

  matrix_vector_mul_d(B, nc, nr, y, incy, x, incx, al, bt) ;
  
  err = 0.0 ;
  
  for ( i = 0 ; i < nc ; i ++ ) {
    err = MAX(fabs(z[i*incx]-x[i*incx]), err) ;
  }

  fprintf(stderr, "transposed error: %lg\n", err) ;

  return 0 ;
}

gint matrix_matrix_multiply_test(gint m, gint n, gint k)

{
  gdouble B[2048], A[2048], C[2048], D[2048] ;
  gdouble x[128], y[128], z[128], al, bt, err ;
  gint i, j, lda, ldb, ldc ;

  fprintf(stderr, "double-precision matrix-matrix multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%dx%d]\n", m, k, k, n) ;
  
  lda = k + 3 ; ldb = n + 4 ; ldc = n + 5 ;
  
  random_matrix_d(A, m, lda) ;
  random_matrix_d(B, k, ldb) ;
  random_matrix_d(C, m, ldc) ;

  al = g_random_double() ;
  bt = g_random_double() ;
  
  memcpy(D, C, m*ldc*sizeof(gdouble)) ;

  matrix_matrix_multiply_d(A, B, m, n, k, lda, ldb, al, bt, D, ldc) ;

  blaswrap_dgemm(FALSE, FALSE, m, n, k, al, A, lda, B, ldb, bt, C, ldc) ;

  err = 0.0 ;

  for ( i = 0 ; i < m ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ ) {
      err = MAX(err, fabs(C[i*ldc+j]-D[i*ldc+j])) ;
    }
  }

  fprintf(stderr, "maximum error: %lg\n", err) ;
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gint nr, nc ;
  gint m, n, k ;

  m = 13 ; n = 4 ; k = 17 ;

  /* m = 2 ; n = 2 ; k = 2 ; */
  
  matrix_matrix_multiply_test(m, n, k) ;

  return 0 ;
  
  nr = 17 ; nc = 5 ;

  matrix_multiply_test_d(nr, nc) ;
  
  return 0 ;
}
