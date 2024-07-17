/* This file is part of BLASWRAP, a set of C wrappers for BLAS routines
 *
 * Copyright (C) 2020, 2024 Michael Carley
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
gint random_matrix_z(gdouble *A, gint nr, gint nc, gboolean sym) ;
gint random_matrix_f(gfloat *A, gint nr, gint nc, gboolean sym) ;
gint matrix_transpose_d(gdouble *B, gdouble *A, gint nr, gint nc) ;
gint matrix_transpose_z(gdouble *B, gdouble *A, gint nr, gint nc) ;
gint matrix_transpose_f(gfloat *B, gfloat *A, gint nr, gint nc) ;

gint matrix_vector_mul_d(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble al, gdouble bt) ;
gint matrix_vector_mul_f(gfloat *A, gint nr, gint nc,
			 gfloat *x, gint incx,
			 gfloat *y, gint incy,
			 gfloat al, gfloat bt) ;
gint matrix_vector_mul_z(gdouble *A, gint nr, gint nc,
			 gdouble *x, gint incx,
			 gdouble *y, gint incy,
			 gdouble *al, gdouble *bt) ;

gint matrix_matrix_multiply_d(gdouble *A, gdouble *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gdouble al, gdouble bt,
			      gdouble *C, gint ldc) ;
gint matrix_matrix_multiply_f(gfloat *A, gfloat *B, gint m, gint n, gint k,
			      gint lda, gint ldb, gfloat al, gfloat bt,
			      gfloat *C, gint ldc) ;
gint matrix_matrix_multiply_test_d(gint m, gint n, gint k,
				   gint incx, gint incy, gint incz) ;
gint matrix_matrix_multiply_test_f(gint m, gint n, gint k,
				   gint incx, gint incy, gint incz) ;
gint matrix_multiply_test_d(gint nr, gint nc, gint incx, gint incy, gint incz) ;
gint matrix_multiply_test_f(gint nr, gint nc, gint incx, gint incy, gint incz) ;
gint matrix_multiply_test_z(gint nr, gint nc, gint incx, gint incy, gint incz) ;

gint matrix_symmetric_multiply_test_d(gint n,
				      gint incx, gint incy, gint incz) ;
gint matrix_symmetric_multiply_test_f(gint n,
				      gint incx, gint incy, gint incz) ;
gint matrix_symmetric_multiply_test_z(gint n,
				      gint incx, gint incy, gint incz) ;

gint matrix_multiply_test_d(gint nr, gint nc, gint incx, gint incy, gint incz)

{
  gdouble B[2048], A[2048], x[128], y[128], z[128], al, bt, err ;
  gint i ;

  fprintf(stderr, "double-precision real matrix-vector multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%d]\n", nr, nc, nc) ;
  
  random_matrix_d(A, nr, nc, FALSE) ;
  random_matrix_d(x, 1, nc*incx, FALSE) ;
  random_matrix_d(y, nr*incy, 1, FALSE) ;

  memcpy(z, y, nr*incy*sizeof(gdouble)) ;

  al = g_random_double() ; 
  bt = g_random_double() ; 

  matrix_vector_mul_d(A, nr, nc, x, incx, y, incy, al, bt) ;

  blaswrap_dgemv(FALSE, nr, nc, al, A, nc, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "plain error     : %lg\n", err) ;

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

gint matrix_multiply_test_z(gint nr, gint nc, gint incx, gint incy, gint incz)

{
  gdouble B[2048], A[2048], x[1024], y[1024], z[1024], al[2], bt[2], err ;
  gint i ;

  fprintf(stderr, "double-precision complex matrix-vector multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%d]\n", nr, nc, nc) ;
  
  random_matrix_z(A, nr, nc, FALSE) ;
  random_matrix_z(x, 1, nc*incx, FALSE) ;
  random_matrix_z(y, nr*incy, 1, FALSE) ;

  for ( i = 0 ; i < nr ; i ++ ) {
    z[2*i*incz+0] = y[2*i*incy+0] ;
    z[2*i*incz+1] = y[2*i*incy+1] ;
  }
    
  al[0] = g_random_double() ; al[1] = g_random_double() ; 
  bt[0] = g_random_double() ; bt[1] = g_random_double() ; 

  matrix_vector_mul_z(A, nr, nc, x, incx, y, incy, al, bt) ;

  blaswrap_zgemv(FALSE, nr, nc, al, A, nc, x, incx, bt, z, incz) ;

  err = 0.0 ;
  for ( i = 0 ; i < nr ; i ++ ) {
    err = MAX((z[i*2*incz+0]-y[i*2*incy+0])*(z[i*2*incz+0]-y[i*2*incy+0]) +
	      (z[i*2*incz+1]-y[i*2*incy+1])*(z[i*2*incz+1]-y[i*2*incy+1]),
	      err) ;
  }
  err = sqrt(err) ;
  
  fprintf(stderr, "plain error     : %lg\n", err) ;

  /*set up for transpose multiplication*/
  for ( i = 0 ; i < nc ; i ++ ) {
    z[2*i*incz+0] = x[2*i*incx+0] ;
    z[2*i*incz+1] = x[2*i*incx+1] ;
  }

  blaswrap_zgemv(TRUE, nr, nc, al, A, nc, y, incy, bt, z, incz) ;

  matrix_transpose_z(B, A, nr, nc) ;

  matrix_vector_mul_z(B, nc, nr, y, incy, x, incx, al, bt) ;
  
  err = 0.0 ;
  for ( i = 0 ; i < nc; i ++ ) {
    err = MAX((z[i*2*incz+0]-x[i*2*incx+0])*(z[i*2*incz+0]-x[i*2*incx+0]) +
	      (z[i*2*incz+1]-x[i*2*incx+1])*(z[i*2*incz+1]-x[i*2*incx+1]),
	      err) ;
  }
  err = sqrt(err) ;

  fprintf(stderr, "transposed error: %lg\n", err) ;
  
  return 0 ;
}

gint matrix_multiply_test_f(gint nr, gint nc, gint incx, gint incy, gint incz)

{
  gfloat B[2048], A[2048], x[128], y[128], z[128], al, bt, err ;
  gint i ;

  fprintf(stderr, "single-precision real matrix-vector multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%d]\n", nr, nc, nc) ;
  
  random_matrix_f(A, nr, nc, FALSE) ;
  random_matrix_f(x, 1, nc*incx, FALSE) ;
  random_matrix_f(y, nr*incy, 1, FALSE) ;

  memcpy(z, y, nr*incy*sizeof(gdouble)) ;

  al = g_random_double() ; 
  bt = g_random_double() ; 

  matrix_vector_mul_f(A, nr, nc, x, incx, y, incy, al, bt) ;

  blaswrap_sgemv(FALSE, nr, nc, al, A, nc, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "plain error     : %lg\n", err) ;

  /*set up for transpose multiplication*/
  memcpy(z, x, nc*incx*sizeof(gfloat)) ;

  blaswrap_sgemv(TRUE, nr, nc, al, A, nc, y, incy, bt, z, incx) ;

  matrix_transpose_f(B, A, nr, nc) ;

  matrix_vector_mul_f(B, nc, nr, y, incy, x, incx, al, bt) ;
  
  err = 0.0 ;
  
  for ( i = 0 ; i < nc ; i ++ ) {
    err = MAX(fabs(z[i*incx]-x[i*incx]), err) ;
  }

  fprintf(stderr, "transposed error: %lg\n", err) ;

  return 0 ;
}

gint matrix_symmetric_multiply_test_d(gint n,
				      gint incx, gint incy, gint incz)

{
  gdouble A[2048], x[128], y[128], z[128], al, bt, err ;
  gint i ;

  fprintf(stderr, "double-precision real symmetric matrix-vector multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%d]\n", n, n, n) ;
  
  random_matrix_d(A, n, n, TRUE) ;
  random_matrix_d(x, 1, n*incx, FALSE) ;
  random_matrix_d(y, n*incy, 1, FALSE) ;

  memcpy(z, y, n*incy*sizeof(gdouble)) ;

  al = g_random_double() ; 
  bt = g_random_double() ; 

  matrix_vector_mul_d(A, n, n, x, incx, y, incy, al, bt) ;

  blaswrap_dsymv(TRUE, n, al, A, n, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "upper error     : %lg\n", err) ;

  memcpy(z, y, n*incy*sizeof(gdouble)) ;
  matrix_vector_mul_d(A, n, n, x, incx, y, incy, al, bt) ;
  blaswrap_dsymv(TRUE, n, al, A, n, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "lower error     : %lg\n", err) ;

  return 0 ;
}

gint matrix_symmetric_multiply_test_f(gint n,
				      gint incx, gint incy, gint incz)

{
  gfloat A[2048], x[128], y[128], z[128], al, bt, err ;
  gint i ;

  fprintf(stderr, "single-precision real symmetric matrix-vector multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%d]\n", n, n, n) ;
  
  random_matrix_f(A, n, n, TRUE) ;
  random_matrix_f(x, 1, n*incx, FALSE) ;
  random_matrix_f(y, n*incy, 1, FALSE) ;

  memcpy(z, y, n*incy*sizeof(gfloat)) ;

  al = g_random_double() ; 
  bt = g_random_double() ; 

  matrix_vector_mul_f(A, n, n, x, incx, y, incy, al, bt) ;

  blaswrap_ssymv(TRUE, n, al, A, n, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "upper error     : %lg\n", err) ;

  memcpy(z, y, n*incy*sizeof(gfloat)) ;
  matrix_vector_mul_f(A, n, n, x, incx, y, incy, al, bt) ;
  blaswrap_ssymv(FALSE, n, al, A, n, x, incx, bt, z, incy) ;

  err = 0.0 ;
  
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX(fabs(z[i*incy]-y[i*incy]), err) ;
  }

  fprintf(stderr, "lower error     : %lg\n", err) ;

  return 0 ;
}

gint matrix_symmetric_multiply_test_z(gint n, gint incx, gint incy, gint incz)

{
  gdouble A[2048], x[1024], y[1024], z[1024], al[2], bt[2], err ;
  gint i ;

  fprintf(stderr,
	  "double-precision complex symmetric matrix-vector multiply\n") ;
  fprintf(stderr, "[%dx%d] x [%d]\n", n, n, n) ;
  
  random_matrix_z(A, n, n, TRUE) ;
  random_matrix_z(x, 1, n*incx, FALSE) ;
  random_matrix_z(y, n*incy, 1, FALSE) ;

  for ( i = 0 ; i < n ; i ++ ) {
    z[2*i*incz+0] = y[2*i*incy+0] ;
    z[2*i*incz+1] = y[2*i*incy+1] ;
  }
    
  al[0] = g_random_double() ; al[1] = g_random_double() ; 
  bt[0] = g_random_double() ; bt[1] = g_random_double() ; 

  matrix_vector_mul_z(A, n, n, x, incx, y, incy, al, bt) ;

  blaswrap_zsymv(TRUE, n, al, A, n, x, incx, bt, z, incz) ;

  err = 0.0 ;
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX((z[i*2*incz+0]-y[i*2*incy+0])*(z[i*2*incz+0]-y[i*2*incy+0]) +
	      (z[i*2*incz+1]-y[i*2*incy+1])*(z[i*2*incz+1]-y[i*2*incy+1]),
	      err) ;
  }
  err = sqrt(err) ;
  
  fprintf(stderr, "upper error     : %lg\n", err) ;

  for ( i = 0 ; i < n ; i ++ ) {
    z[2*i*incz+0] = y[2*i*incy+0] ;
    z[2*i*incz+1] = y[2*i*incy+1] ;
  }

  matrix_vector_mul_z(A, n, n, x, incx, y, incy, al, bt) ;

  blaswrap_zsymv(FALSE, n, al, A, n, x, incx, bt, z, incz) ;

  err = 0.0 ;
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX((z[i*2*incz+0]-y[i*2*incy+0])*(z[i*2*incz+0]-y[i*2*incy+0]) +
	      (z[i*2*incz+1]-y[i*2*incy+1])*(z[i*2*incz+1]-y[i*2*incy+1]),
	      err) ;
  }
  err = sqrt(err) ;
  
  fprintf(stderr, "lower error     : %lg\n", err) ;
  
  /* /\*set up for transpose multiplication*\/ */
  /* for ( i = 0 ; i < nc ; i ++ ) { */
  /*   z[2*i*incz+0] = x[2*i*incx+0] ; */
  /*   z[2*i*incz+1] = x[2*i*incx+1] ; */
  /* } */

  /* blaswrap_zgemv(TRUE, nr, nc, al, A, nc, y, incy, bt, z, incz) ; */

  /* matrix_transpose_z(B, A, nr, nc) ; */

  /* matrix_vector_mul_z(B, nc, nr, y, incy, x, incx, al, bt) ; */
  
  /* err = 0.0 ; */
  /* for ( i = 0 ; i < nc; i ++ ) { */
  /*   err = MAX((z[i*2*incz+0]-x[i*2*incx+0])*(z[i*2*incz+0]-x[i*2*incx+0]) + */
  /* 	      (z[i*2*incz+1]-x[i*2*incx+1])*(z[i*2*incz+1]-x[i*2*incx+1]), */
  /* 	      err) ; */
  /* } */
  /* err = sqrt(err) ; */

  /* fprintf(stderr, "transposed error: %lg\n", err) ; */
  
  return 0 ;
}

gint matrix_matrix_multiply_test_d(gint m, gint n, gint k,
				   gint incx, gint incy, gint incz)

{
  gdouble B[2048], A[2048], C[2048], D[2048], Bt[2048], At[2048] ;
  gdouble al, bt, err ;
  gint i, j, lda, ldb, ldc ;

  lda = k + incx - 1 ; ldb = n + incy - 1 ; ldc = n + incz - 1 ;
  
  fprintf(stderr, "double-precision real matrix-matrix multiply\n") ;
  fprintf(stderr, "[%dx%d(%d)] x [%dx%d(%d)]\n", m, k, lda, k, n, ldb) ;
  
  random_matrix_d(A, m, lda, FALSE) ;
  random_matrix_d(B, k, ldb, FALSE) ;
  random_matrix_d(C, m, ldc, FALSE) ;

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

  fprintf(stderr, "A x B error: %lg\n", err) ;

  /*A*B^T*/
  matrix_transpose_d(Bt, B, k, ldb) ;
  memcpy(D, C, m*ldc*sizeof(gdouble)) ;

  /*B:  k x n*/
  matrix_matrix_multiply_d(A, B, m, n, k, lda, ldb, al, bt, D, ldc) ;
  /*Bt: n x k*/
  blaswrap_dgemm(FALSE, TRUE, m, n, k, al, A, lda, Bt, k, bt, C, ldc) ;

  err = 0.0 ;
  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < k ; j ++ ) {
      err = MAX(err, fabs(C[i*ldc+j]-D[i*ldc+j])) ;
    }
  }

  fprintf(stderr, "A x B^T error: %lg\n", err) ;

  /*A^T*B*/
  matrix_transpose_d(At, A, m, lda) ;
  memcpy(D, C, m*ldc*sizeof(gdouble)) ;
  /*A:  m x k*/
  matrix_matrix_multiply_d(A, B, m, n, k, lda, ldb, al, bt, D, ldc) ;
  /*At: k x m*/
  blaswrap_dgemm(TRUE, FALSE, m, n, k, al, At, m, B, ldb, bt, C, ldc) ;

  err = 0.0 ;
  for ( i = 0 ; i < m ; i ++ ) {
    for ( j = 0 ; j < k ; j ++ ) {
      err = MAX(err, fabs(C[i*ldc+j]-D[i*ldc+j])) ;
    }
  }

  fprintf(stderr, "A^T x B error: %lg\n", err) ;
  
  return 0 ;
}

gint matrix_matrix_multiply_test_f(gint m, gint n, gint k,
				   gint incx, gint incy, gint incz)

{
  gfloat B[2048], A[2048], C[2048], D[2048], Bt[2048] ;
  gfloat al, bt, err ;
  gint i, j, lda, ldb, ldc ;

  lda = k + incx - 1 ; ldb = n + incy - 1 ; ldc = n + incz - 1 ;
  
  fprintf(stderr, "single-precision real matrix-matrix multiply\n") ;
  fprintf(stderr, "[%dx%d(%d)] x [%dx%d(%d)]\n", m, k, lda, k, n, ldb) ;
  
  random_matrix_f(A, m, lda, FALSE) ;
  random_matrix_f(B, k, ldb, FALSE) ;
  random_matrix_f(C, m, ldc, FALSE) ;

  al = g_random_double() ;
  bt = g_random_double() ;

  memcpy(D, C, m*ldc*sizeof(gfloat)) ;

  matrix_matrix_multiply_f(A, B, m, n, k, lda, ldb, al, bt, D, ldc) ;

  blaswrap_sgemm(FALSE, FALSE, m, n, k, al, A, lda, B, ldb, bt, C, ldc) ;

  err = 0.0 ;

  for ( i = 0 ; i < m ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ ) {
      err = MAX(err, fabs(C[i*ldc+j]-D[i*ldc+j])) ;
    }
  }

  fprintf(stderr, "A x B error: %lg\n", err) ;

  matrix_transpose_f(Bt, B, k, ldb) ;
  memcpy(D, C, m*ldc*sizeof(gfloat)) ;

  /*B:  k x n*/
  matrix_matrix_multiply_f(A, B, m, n, k, lda, ldb, al, bt, D, ldc) ;
  /*Bt: n x k*/
  blaswrap_sgemm(FALSE, TRUE, m, n, k, al, A, lda, Bt, k, bt, C, ldc) ;

  err = 0.0 ;
  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < k ; j ++ ) {
      err = MAX(err, fabs(C[i*ldc+j]-D[i*ldc+j])) ;
      /* fprintf(stderr, "%g %g\n", C[i*ldc+j], D[i*ldc+j]) ; */
    }
  }

  fprintf(stderr, "A x B^T error: %lg\n", err) ;
  
  return 0 ;
}

static gint triangular_solve_test_d(gint n)

{
  gdouble *A, *x, *y, *z, al, bt, err ;
  gint i, j, one = 1, info ;
  
  fprintf(stderr, "double-precision triangular solve\n") ;
  fprintf(stderr, "n = %d\n", n) ;

  A = (gdouble *)g_malloc(n*n*sizeof(gdouble)) ;
  x = (gdouble *)g_malloc(n  *sizeof(gdouble)) ;
  y = (gdouble *)g_malloc(n  *sizeof(gdouble)) ;
  z = (gdouble *)g_malloc(n  *sizeof(gdouble)) ;
  
  fprintf(stderr, "upper-triangular solve\n") ;
  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < i ; j ++ ) A[i*n+j] = 0.0 ;
    for ( j = i ; j < n ; j ++ ) A[i*n+j] = g_random_double() ;
    x[i] = y[i] = g_random_double() ;
  }
  /*upper triangular, not-transposed*/
  blaswrap_dtrtrs(TRUE, FALSE, FALSE, n, one, A, n, y, n, &info) ;

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemv(FALSE, n, n, al, A, n, y, one, bt, z, one) ;

  err = 0.0 ;
  for ( i = 0 ; i < n ; i ++ ) {
    err = MAX(err, fabs(x[i] - z[i])) ;
  }
  fprintf(stderr, "maximum error: %lg\n", err) ;
  
  return 0 ;
}

static gint parse_test(char *str)

{
  char *tests[] = {"matrix_vector", "matrix_matrix", "symmetric",
    "triangular", NULL} ;
  gint i ;

  for ( i = 0 ; tests[i] != NULL ; i ++ ) {
    if ( strcmp(tests[i], str) == 0 ) return i+1 ;
  }

  return 0 ;
}  

gint main(gint argc, char **argv)

{
  gint m, n, k, incx, incy, incz, test ;
  char ch, *progname ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  
  m = 13 ; n = 4 ; k = 17 ;
  incx = 1 ; incy = 1 ; incz = 1 ;
  test = -1 ;
  
  while ( (ch = getopt(argc, argv, "k:m:n:t:x:y:z:")) != EOF ) {
    switch (ch) {
    default:
      fprintf(stderr, "%s: option not recognized\n", progname) ;
      exit(1) ;
      break ;
    case 'k': k = atoi(optarg) ; break ;
    case 'm': m = atoi(optarg) ; break ;
    case 'n': n = atoi(optarg) ; break ;
    case 't': test = parse_test(optarg) ; break ;
    case 'x': incx = atoi(optarg) ; break ;
    case 'y': incy = atoi(optarg) ; break ;
    case 'z': incz = atoi(optarg) ; break ;
    }
  }

  if ( test == 0 ) {
    fprintf(stderr, "unrecognised test specified with -t option\n") ;
    return 1 ;
  }
  
  if ( test == 1 || test == -1 ) {
    matrix_multiply_test_d(m, n, incx, incy, incz) ;
    fprintf(stderr, "\n") ;
    matrix_multiply_test_f(m, n, incx, incy, incz) ;
    fprintf(stderr, "\n") ;
    matrix_multiply_test_z(m, n, incx, incy, incz) ;
    fprintf(stderr, "\n") ;

    if ( test > 0 ) return 0 ;
  }

  if ( test == 2 || test == -1 ) {
    matrix_matrix_multiply_test_d(m, n, k, incx, incy, incz) ;
    fprintf(stderr, "\n") ;
    matrix_matrix_multiply_test_f(m, n, k, incx, incy, incz) ;
    fprintf(stderr, "\n") ;

    if ( test > 0 ) return 0 ;
  }

  if ( test == 3 || test == -1 ) {
    matrix_symmetric_multiply_test_d(n, incx, incy, incz) ;
    fprintf(stderr, "\n") ;
    matrix_symmetric_multiply_test_f(n, incx, incy, incz) ;
    fprintf(stderr, "\n") ;
    matrix_symmetric_multiply_test_z(n, incx, incy, incz) ;
    fprintf(stderr, "\n") ;

    if ( test > 0 ) return 0 ;
  }
  
  if ( test == 4 || test == -1 ) {
    triangular_solve_test_d(n) ;

    if ( test > 0 ) return 0 ;
  }
  
  return 0 ;
}
