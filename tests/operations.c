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

gint random_matrix_d(gdouble *A, gint nr, gint nc, gboolean sym) ;
gint matrix_vector_mul_d(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble al, gdouble bt) ;
gint matrix_transpose_d(gdouble *B, gdouble *A, gint nr, gint nc) ;
gint matrix_matrix_multiply_d(gdouble *A, gdouble *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gdouble al, gdouble bt,
			      gdouble *C, gint ldc) ;
gint random_matrix_f(gfloat *A, gint nr, gint nc) ;
gint matrix_vector_mul_f(gfloat *A, gint nr, gint nc,
			 gfloat *x, gint incx,
			 gfloat *y, gint incy,
			 gfloat al, gfloat bt) ;
gint matrix_vector_mul_c(gfloat *A, gint nr, gint nc,
			 gfloat *x, gint incx,
			 gfloat *y, gint incy,
			 gfloat *al, gfloat *bt) ;
gint matrix_transpose_f(gfloat *B, gfloat *A, gint nr, gint nc) ;
gint matrix_matrix_multiply_f(gfloat *A, gfloat *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gfloat al, gfloat bt,
			      gfloat *C, gint ldc) ;
gint random_matrix_z(gdouble *A, gint nr, gint nc) ;
gint matrix_vector_mul_z(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble *al, gdouble *bt) ;
gint matrix_transpose_z(gdouble *B, gdouble *A, gint nr, gint nc) ;
gint matrix_matrix_multiply_z(gdouble *A, gdouble *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gdouble *al, gdouble *bt,
			      gdouble *C, gint ldc) ;


gint random_matrix_d(gdouble *A, gint nr, gint nc, gboolean sym)

{
  gint i, j ;

  if ( !sym ) {
    for ( i = 0 ; i < nr*nc ; i ++ ) A[i] = g_random_double() ;
    return 0 ;
  }

  if ( nr != nc ) {
    g_error("%s: symmetric matrix must be square (nr = %d, nc = %d)",
	    __FUNCTION__, nr, nc) ;
  }

  for ( i = 0 ; i < nc ; i ++ ) {
    for ( j = i ; j < nc ; j ++ ) {
      A[i*nc+j] = g_random_double() ;
      A[j*nc+i] = A[i*nc+j] ;
    }
  }
  
  return 0 ;
}

gint matrix_vector_mul_d(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble al, gdouble bt)

/*
  y := al*A*x + bt*y
*/

{
  gint i, j ;

  for ( i = 0 ; i < nr ; i ++ ) {
    y[i*incy] *= bt ;
    for ( j = 0 ; j < nc ; j ++ ) {
      y[i*incy] += al*A[i*nc+j]*x[j*incx] ;
    }
  }
  
  return 0 ;
}

gint matrix_transpose_d(gdouble *B, gdouble *A, gint nr, gint nc)

{
  gint i, j ;

  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      B[j*nr+i] = A[i*nc+j] ;
    }
  }
  
  return 0 ;
}

gint matrix_matrix_multiply_d(gdouble *A, gdouble *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gdouble al, gdouble bt,
			      gdouble *C, gint ldc)

/*
  A: m x k
  B: k x n
  C: m x n
*/
  
{
  gint ii, jj, kk ;

  for ( ii = 0 ; ii < m ; ii ++ ) {
    for ( jj = 0 ; jj < n ; jj ++ ) {
      C[ii*ldc + jj] *= bt ;
      for ( kk = 0 ; kk < k ; kk ++ ) {
	C[ii*ldc + jj] += al*A[ii*lda+kk]*B[kk*ldb+jj] ;
      }
    }
  }
  
  return 0 ;
}

gint random_matrix_f(gfloat *A, gint nr, gint nc)

{
  gint i ;

  for ( i = 0 ; i < nr*nc ; i ++ ) A[i] = g_random_double() ;
  
  return 0 ;
}

gint matrix_vector_mul_f(gfloat *A, gint nr, gint nc,
			 gfloat *x, gint incx,
			 gfloat *y, gint incy,
			 gfloat al, gfloat bt)

/*
  y := al*A*x + bt*y
*/

{
  gint i, j ;

  for ( i = 0 ; i < nr ; i ++ ) {
    y[i*incy] *= bt ;
    for ( j = 0 ; j < nc ; j ++ ) {
      y[i*incy] += al*A[i*nc+j]*x[j*incx] ;
    }
  }
  
  return 0 ;
}

gint matrix_vector_mul_c(gfloat *A, gint nr, gint nc,
			 gfloat *x, gint incx,
			 gfloat *y, gint incy,
			 gfloat *al, gfloat *bt)

/*
  y := al*A*x + bt*y
*/

{
  gint i, j ;
  gfloat ar, ai ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    ar = y[i*2*incy+0]*bt[0] - y[i*2*incy+1]*bt[1] ;
    ai = y[i*2*incy+1]*bt[0] + y[i*2*incy+0]*bt[1] ;
    y[i*2*incy+0] = ar ; y[i*2*incy+1] = ai ;
    for ( j = 0 ; j < nc ; j ++ ) {
      ar = A[i*2*nc+2*j+0]*x[j*2*incx+0] - A[i*2*nc+2*j+1]*x[j*2*incx+1] ;
      ai = A[i*2*nc+2*j+1]*x[j*2*incx+0] + A[i*2*nc+2*j+0]*x[j*2*incx+1] ;

      y[i*2*incy+0] += al[0]*ar - al[1]*ai ;
      y[i*2*incy+1] += al[1]*ar + al[0]*ai ;
    }
  }
  
  return 0 ;
}

gint matrix_transpose_f(gfloat *B, gfloat *A, gint nr, gint nc)

{
  gint i, j ;

  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      B[j*nr+i] = A[i*nc+j] ;
    }
  }
  
  return 0 ;
}

gint matrix_matrix_multiply_f(gfloat *A, gfloat *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gfloat al, gfloat bt,
			      gfloat *C, gint ldc)

{
  gint ii, jj, kk ;

  for ( ii = 0 ; ii < m ; ii ++ ) {
    for ( jj = 0 ; jj < n ; jj ++ ) {
      C[ii*ldc + jj] *= bt ;
      for ( kk = 0 ; kk < k ; kk ++ ) {
	C[ii*ldc + jj] += al*A[ii*lda+kk]*B[kk*ldb+jj] ;
      }
    }
  }
  
  return 0 ;
}

gint random_matrix_z(gdouble *A, gint nr, gint nc)

{
  gint i ;

  for ( i = 0 ; i < 2*nr*nc ; i ++ ) A[i] = g_random_double() ;
  
  return 0 ;
}

gint matrix_vector_mul_z(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble *al, gdouble *bt)

/*
  y := al*A*x + bt*y
*/

{
  gint i, j ;
  gdouble ar, ai ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    ar = y[i*2*incy+0]*bt[0] - y[i*2*incy+1]*bt[1] ;
    ai = y[i*2*incy+1]*bt[0] + y[i*2*incy+0]*bt[1] ;
    y[i*2*incy+0] = ar ; y[i*2*incy+1] = ai ;
    for ( j = 0 ; j < nc ; j ++ ) {
      ar = A[i*2*nc+2*j+0]*x[j*2*incx+0] - A[i*2*nc+2*j+1]*x[j*2*incx+1] ;
      ai = A[i*2*nc+2*j+1]*x[j*2*incx+0] + A[i*2*nc+2*j+0]*x[j*2*incx+1] ;

      y[i*2*incy+0] += al[0]*ar - al[1]*ai ;
      y[i*2*incy+1] += al[1]*ar + al[0]*ai ;
    }
  }
  
  return 0 ;
}

gint matrix_transpose_z(gdouble *B, gdouble *A, gint nr, gint nc)

{
  gint i, j ;

  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      B[2*(j*nr+i)+0] = A[2*(i*nc+j)+0] ;
      B[2*(j*nr+i)+1] = A[2*(i*nc+j)+1] ;
    }
  }
  
  return 0 ;
}

gint matrix_matrix_multiply_z(gdouble *A, gdouble *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gdouble *al, gdouble *bt,
			      gdouble *C, gint ldc)

{
  /* gint ii, jj, kk ; */

  /* for ( ii = 0 ; ii < m ; ii ++ ) { */
  /*   for ( jj = 0 ; jj < n ; jj ++ ) { */
  /*     C[ii*ldc + jj] *= bt ; */
  /*     for ( kk = 0 ; kk < k ; kk ++ ) { */
  /* 	C[ii*ldc + jj] += al*A[ii*lda+kk]*B[kk*ldb+jj] ; */
  /*     } */
  /*   } */
  /* } */
  
  return 0 ;
}
