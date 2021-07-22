/* This file is part of BLASWRAP, a set of wrappers for using Fortran
 * BLAS in C programs.
 *
 * Copyright (C) 2020, 2021 Michael Carley
 *
 * BLASWRAP is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  BLASWRAP is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BLASWRAP.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef BLASWRAP_H_INCLUDED 
#define BLASWRAP_H_INCLUDED

#include <glib.h>

extern void dswap_(gint *N, gdouble *dx, gint *incx, gdouble *dy, gint *incy) ;
extern void sswap_(gint *N, gfloat *dx, gint *incx, gfloat *dy, gint *incy) ;

extern void dgemv_(gchar *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;
extern void sgemv_(gchar *trans, gint *m, gint *n, gfloat *alpha,
		   gfloat *A, gint *lda, gfloat *v, gint *incx,
		   gfloat *beta, gfloat *y, gint *incy) ;

extern void zgemv_(gchar *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;

extern void dgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;
extern void zgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;

extern void zgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;
extern void dgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;

extern void dgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;
extern void sgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gfloat *alpha,
		   gfloat *A, gint *lda,
		   gfloat *B, gint *ldb,
		   gfloat *beta, gfloat *C, gint *ldc) ;
extern void zgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;

extern void drotg_(gdouble *da, gdouble *db, gdouble *c, gdouble *s) ;
extern void srotg_(gfloat *da, gfloat *db, gfloat *c, gfloat *s) ;

extern void drot_(gint *N, gdouble *dx, gint *incx, gdouble *dy, gint *incy,
		  gdouble *C, gdouble *S) ;
extern void srot_(gint *N, gfloat *dx, gint *incx, gfloat *dy, gint *incy,
		  gfloat *C, gfloat *S) ;

extern gint dtrtrs_(gchar *uplo, gchar *trans, gchar *diag,
		    gint *n, gint *nrhs, gdouble *a, gint *lda, gdouble *b,
		    gint *ldb, gint *info) ;


extern void dscal_(gint *n, gdouble *da, gdouble *dx, gint *incx) ;
extern void sscal_(gint *n, gfloat  *da, gfloat  *dx, gint *incx) ;

extern gdouble dasum_ (gint *n, gdouble *x, gint *incx) ;
extern gdouble dzasum_(gint *n, gdouble *x, gint *incx) ;
extern gdouble dznrm2_(gint *n, gdouble *x, gint *incx) ;
extern gdouble dnrm2_ (gint *n, gdouble *x, gint *incx) ;
extern gfloat  snrm2_ (gint *n, gfloat *x, gint *incx) ;
extern gint    idamax_(gint *n, gdouble *x, gint *incx) ;
extern gint    izamax_(gint *n, gdouble *x, gint *incx) ;
extern gdouble ddot_  (gint *n, gdouble *x, gint *incx, 
		       gdouble *y, gint *incy) ;
extern gfloat  sdot_  (gint *n, gfloat *x, gint *incx, 
		       gfloat *y, gint *incy) ;
/* extern gsl_complex zdotu_ (gint *n, gdouble *x, gint *incx,  */
/* 			   gdouble *y, gint *incy) ; */
extern void    scopy_(gint *n, 
		      gfloat *x, gint *incx,
		      gfloat *y, gint *incy) ;
extern void    dcopy_(gint *n, 
		      gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    saxpy_(gint *n, 
		      gfloat *a, gfloat *x, gint *incx,
		      gfloat *y, gint *incy) ;
extern void    daxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    zaxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;

/*macros wrapping BLAS calls*/

/* copy x into y */
#define blaswrap_dcopy(_n,_x,_strx,_y,_stry)	\
  dcopy_(&(_n),(_x),&(_strx),(_y),&(_stry))
#define blaswrap_scopy(_n,_x,_strx,_y,_stry)	\
  scopy_(&(_n),(_x),&(_strx),(_y),&(_stry))

/* x := x*al */
#define blaswrap_dscal(_n,_al,_x,_strx) dscal_(&(_n),&(_al),(_x),&(_strx))
#define blaswrap_sscal(_n,_al,_x,_strx) sscal_(&(_n),&(_al),(_x),&(_strx))

/* swap x and y */
#define blaswrap_dswap(_n,_x,_strx,_y,_stry)	\
  dswap_(&(_n),(_x),&(_strx),(_y),&(_stry))
#define blaswrap_sswap(_n,_x,_strx,_y,_stry)	\
  sswap_(&(_n),(_x),&(_strx),(_y),&(_stry))

/*x'*x*/
#define blaswrap_dnrm2(_n,_x,_incx) dnrm2_(&(_n),(_x),&(_incx))
#define blaswrap_snrm2(_n,_x,_incx) snrm2_(&(_n),(_x),&(_incx))

/* sum x[i]*y[i] */
#define blaswrap_ddot(_n,_x,_strx,_y,_stry)	\
  ddot_(&(_n), (_x), &(_strx), (_y), &(_stry)) 
#define blaswrap_sdot(_n,_x,_strx,_y,_stry)	\
  sdot_(&(_n), (_x), &(_strx), (_y), &(_stry)) 

/* y := y + a*x */
#define blaswrap_saxpy(_n,_a,_x,_strx,_y,_stry)		\
  saxpy_(&(_n), &(_a), (_x), &(_strx), (_y), &(_stry))
#define blaswrap_daxpy(_n,_a,_x,_strx,_y,_stry)		\
  daxpy_(&(_n), &(_a), (_x), &(_strx), (_y), &(_stry))

/* y := al*A*x + bt*y */

/*trans nr nc al A lda x strx bt y stry*/
#define blaswrap_dgemv(_t,_m,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)	\
  do {									\
    if ( (_t) ) {							\
      dgemv_("N",&(_n),&(_m),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    } else {								\
      dgemv_("T",&(_n),&(_m),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    }									\
  } while (0)

#define blaswrap_sgemv(_t,_m,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)	\
  do {									\
    if ( (_t) ) {							\
      sgemv_("N",&(_n),&(_m),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    } else {								\
      sgemv_("T",&(_n),&(_m),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    }									\
  } while (0)

#define blaswrap_zgemv(_t,_m,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)	\
  do {									\
    if ( ((_t) == TRUE) ) {						\
      zgemv_("N",&(_n),&(_m),(_al),(_A),&(_lda),(_x),&(_incx),		\
	     (_bt),(_y),&(_incy)) ;					\
    } else {								\
      zgemv_("T",&(_n),&(_m),(_al),(_A),&(_lda),(_x),&(_incx),		\
	     (_bt),(_y),&(_incy)) ;					\
    }									\
  } while (0)

/* C := al*A*B + bt*C */
/* A [m x k] leading dimension lda*/
/* B [k x n] leading dimension ldb*/
/* C [m x n] leading dimension ldc*/
#define blaswrap_dgemm(_ta,_tb,_m,_n,_k,_al,_A,_lda,_B,_ldb,_bt,_C,_ldc) \
  do {									\
    if ( !(_ta) ) {							\
      if ( !(_tb) ) {							\
	dgemm_("N","N", &(_n),&(_m),&(_k),&(_al),(_B),&(_ldb),		\
	       (_A),&(_lda),&(_bt),(_C),&(_ldc)) ;			\
      } else {								\
	dgemm_("T","N", &(_n),&(_m),&(_k),&(_al),(_B),&(_ldb),		\
	       (_A),&(_lda),&(_bt),(_C),&(_ldc)) ;			\
      }									\
    } else {								\
      g_assert_not_reached() ;						\
    }									\
  } while (0)
  
#define blaswrap_sgemm(_ta,_tb,_m,_n,_k,_al,_A,_lda,_B,_ldb,_bt,_C,_ldc) \
  do {									\
    if ( !(_ta) ) {							\
      if ( !(_tb) ) {							\
	sgemm_("N","N", &(_n),&(_m),&(_k),&(_al),(_B),&(_ldb),		\
	       (_A),&(_lda),&(_bt),(_C),&(_ldc)) ;			\
      } else {								\
	sgemm_("T","N", &(_n),&(_m),&(_k),&(_al),(_B),&(_ldb),		\
	       (_A),&(_lda),&(_bt),(_C),&(_ldc)) ;			\
      }									\
    } else {								\
      g_assert_not_reached() ;						\
    }									\
  } while (0)

/*Givens rotations*/
#define blaswrap_drotg(_a,_b,_c,_s)					\
  do {									\
  gdouble _ta = _a, _tb = _b ;						\
  drotg_(&(_ta),&(_tb),(_c),(_s)) ;					\
  } while (0)
#define blaswrap_srotg(_a,_b,_c,_s)					\
  do {									\
  gfloat _ta = _a, _tb = _b ;						\
  srotg_(&(_ta),&(_tb),(_c),(_s)) ;					\
  } while (0)

/*triangular system solves*/
#define blaswrap_dtrtrs(_upper,_trans,_diag,_n,_nrhs,_a,_lda,_b,_ldb,_info) \
  do {									\
    gchar _upstr[1], _tstr[1], _dstr[1] ;				\
    _dstr[0] = ( _diag == TRUE ? 'U' : 'N') ;				\
    if ( (_upper) ) {							\
      if ( !(_trans) ) {						\
	_upstr[0] = 'L' ; _tstr[0] = 'T' ;				\
      }	else {								\
	g_assert_not_reached() ;					\
      }									\
    }									\
    if ( (!_upper) ) {							\
      if ( !(_trans) ) {						\
	g_assert_not_reached() ;					\
	_upstr[0] = 'L' ; _tstr[0] = 'N' ;				\
      } else {								\
	g_assert_not_reached() ;					\
      }									\
    }									\
    dtrtrs_(_upstr,_tstr,_dstr,&(_n),&(_nrhs),(_a),&(_lda),(_b),&(_ldb),(_info)) ;  \
  }  while (0) 

#endif /*BLASWRAP_H_INCLUDED*/
