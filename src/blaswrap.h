/* This file is part of BLASWRAP, a set of wrappers for using Fortran
 * BLAS in C programs.
 *
 * Copyright (C) 2020, 2021, 2024 Michael Carley
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

typedef struct {gdouble dat[2] ; } blaswrap_complex_double_t ;
typedef struct {gfloat dat[2] ; } blaswrap_complex_float_t ;

extern void dswap_(gint *N, gdouble *dx, gint *incx, gdouble *dy, gint *incy) ;
extern void sswap_(gint *N, gfloat *dx, gint *incx, gfloat *dy, gint *incy) ;

extern void dgemv_(char *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;
extern void sgemv_(char *trans, gint *m, gint *n, gfloat *alpha,
		   gfloat *A, gint *lda, gfloat *v, gint *incx,
		   gfloat *beta, gfloat *y, gint *incy) ;

extern void zgemv_(char *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;

extern void dsymv_(char *uplo, gint *n, gdouble *alpha, gdouble *A,
		   gint *lda, gdouble *x, gint *incx, gdouble *beta,
		   gdouble *y, gint *incy) ;
extern void ssymv_(char *uplo, gint *n, gfloat *alpha, gfloat *A,
		   gint *lda, gfloat *x, gint *incx, gfloat *beta,
		   gfloat *y, gint *incy) ;
extern void zsymv_(char *uplo, gint *n, gdouble *alpha, gdouble *A,
		   gint *lda, gdouble *x, gint *incx, gdouble *beta,
		   gdouble *y, gint *incy) ;

extern void dgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;
extern void zgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;

extern void zgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;
extern void dgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;

extern void dgbmv_(char *trans, gint *m, gint *n, gint *kl, gint *ku,
		   gdouble *alpha, gdouble *A, gint *lda, gdouble *x,
		   gint *incx, gdouble *beta, gdouble *y, gint *incy) ;

extern void dgemm_(char *transa, char *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;
extern void sgemm_(char *transa, char *transb,
		   gint *m, gint *n, gint *k, gfloat *alpha,
		   gfloat *A, gint *lda,
		   gfloat *B, gint *ldb,
		   gfloat *beta, gfloat *C, gint *ldc) ;
extern void zgemm_(char *transa, char *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;
extern void cgemm_(char *transa, char *transb,
		   gint *m, gint *n, gint *k, gfloat *alpha,
		   gfloat *A, gint *lda,
		   gfloat *B, gint *ldb,
		   gfloat *beta, gfloat *C, gint *ldc) ;

extern void drotg_(gdouble *da, gdouble *db, gdouble *c, gdouble *s) ;
extern void srotg_(gfloat *da, gfloat *db, gfloat *c, gfloat *s) ;
extern void zrotg_(gdouble *da, gdouble *db, gdouble *c, gdouble *s) ;

extern void dgeqp3_(gint *m, gint *n, gdouble *A, gint *lda, gint *jpvt,
		    gdouble *tau, gdouble *work, gint *lwork, gint *info) ;

extern void dlarf_(char *side, gint *m, gint *n, gdouble *v, gint *incv,
		   gdouble *tau, gdouble *C, gint *ldc, gdouble *work) ;

extern void drot_(gint *N, gdouble *dx, gint *incx, gdouble *dy, gint *incy,
		  gdouble *C, gdouble *S) ;
extern void srot_(gint *N, gfloat *dx, gint *incx, gfloat *dy, gint *incy,
		  gfloat *C, gfloat *S) ;
extern void dtrtrs_(char *uplo, char *trans, char *diag,
		    gint *n, gint *nrhs, gdouble *a, gint *lda, gdouble *b,
		    gint *ldb, gint *info) ;
extern void ztrtrs_(char *uplo, char *trans, char *diag,
		    gint *n, gint *nrhs, gdouble *a, gint *lda, gdouble *b,
		    gint *ldb, gint *info) ;
extern void dtrtri_(char *uplo, char *diag, gint *n, gdouble *a,
		    gint *lda, gint *info) ;

extern void dtpsv_(char *uplo, char *trans, char *diag, gint *n,
		   gdouble *ap, gdouble *x, gint *incx) ;

extern void dtptrs_(char *uplo, char *trans,char *diag, gint *n,
		    gint *nrhs, gdouble *ap, gdouble *b, gint *ldb,
		    gint *info) ;
extern void dgesvd_(char *jobu, char *jobvt, gint *m, gint *n,
		    gdouble *a, gint *lda, gdouble *s, gdouble *u,
		    gint *ldu, gdouble *vt,
		    gint *ldvt, gdouble *work, gint *lwork, gint *info) ;

extern void dgels_(char *trans, gint *m, gint *n, gint *nrhs,
		   gdouble *a, gint *lda, gdouble *b, gint *ldb,
		   gdouble *work, gint *lwork, gint *info) ;

/*LDL factorization*/
extern void dsptrf_(char *uplo, gint *N, gdouble *A, gint *ipiv, gint *info) ;

extern void dscal_(gint *n, gdouble *da, gdouble *dx, gint *incx) ;
extern void sscal_(gint *n, gfloat  *da, gfloat  *dx, gint *incx) ;
extern void zscal_(gint *n, gdouble *da, gdouble *dx, gint *incx) ;

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
extern blaswrap_complex_double_t zdotu_(gint *n, gdouble *x, gint *incx,
					gdouble *y, gint *incy) ;
extern blaswrap_complex_float_t cdotu_(gint *n, gfloat *x, gint *incx,
				       gfloat *y, gint *incy) ;
extern void    scopy_(gint *n, 
		      gfloat *x, gint *incx,
		      gfloat *y, gint *incy) ;
extern void    dcopy_(gint *n, 
		      gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    zcopy_(gint *n, 
		      gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    saxpy_(gint *n, 
		      gfloat *a, gfloat *x, gint *incx,
		      gfloat *y, gint *incy) ;
extern void    daxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    caxpy_(gint *n, 
		      gfloat *a, gfloat *x, gint *incx,
		      gfloat *y, gint *incy) ;
extern void    zaxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;

/*macros wrapping BLAS calls*/

/* copy x into y */
#define blaswrap_dcopy(_n,_x,_strx,_y,_stry)	\
  dcopy_(&(_n),(_x),&(_strx),(_y),&(_stry))
#define blaswrap_scopy(_n,_x,_strx,_y,_stry)	\
  scopy_(&(_n),(_x),&(_strx),(_y),&(_stry))
#define blaswrap_zcopy(_n,_x,_strx,_y,_stry)	\
  zcopy_(&(_n),(_x),&(_strx),(_y),&(_stry))

/* x := x*al */
#define blaswrap_dscal(_n,_al,_x,_strx) dscal_(&(_n),&(_al),(_x),&(_strx))
#define blaswrap_sscal(_n,_al,_x,_strx) sscal_(&(_n),&(_al),(_x),&(_strx))
#define blaswrap_zscal(_n,_al,_x,_strx) dscal_(&(_n), (_al),(_x),&(_strx))

/* swap x and y */
#define blaswrap_dswap(_n,_x,_strx,_y,_stry)	\
  dswap_(&(_n),(_x),&(_strx),(_y),&(_stry))
#define blaswrap_sswap(_n,_x,_strx,_y,_stry)	\
  sswap_(&(_n),(_x),&(_strx),(_y),&(_stry))

/*x'*x*/
#define blaswrap_dnrm2(_n,_x,_incx) dnrm2_(&(_n),(_x),&(_incx))
#define blaswrap_snrm2(_n,_x,_incx) snrm2_(&(_n),(_x),&(_incx))
#define blaswrap_dznrm2(_n,_x,_incx) dznrm2_(&(_n),(_x),&(_incx))

/* sum x[i]*y[i] */
#define blaswrap_ddot(_n,_x,_strx,_y,_stry)	\
  ddot_(&(_n), (_x), &(_strx), (_y), &(_stry)) 
#define blaswrap_sdot(_n,_x,_strx,_y,_stry)	\
  sdot_(&(_n), (_x), &(_strx), (_y), &(_stry)) 
#define blaswrap_cdotu(_out,_n,_x,_strx,_y,_stry)		\
  do {								\
    g_assert_not_reached() ; /*unchecked code*/			\
    blaswrap_complex_float_t _z ;				\
    _z = cdotu_(&(_n), (_x), &(_strx), (_y), &(_stry)) ;	\
    (_out)[0] = _z.dat[0] ; (_out)[1] = _z.dat[1] ;		\
  } while (0) 
#define blaswrap_zdotu(_out,_n,_x,_strx,_y,_stry)		\
  do {								\
    blaswrap_complex_double_t _z ;				\
    _z = zdotu_(&(_n), (_x), &(_strx), (_y), &(_stry)) ;	\
    (_out)[0] = _z.dat[0] ; (_out)[1] = _z.dat[1] ;		\
  } while (0) 

  
/* y := y + a*x */
#define blaswrap_saxpy(_n,_a,_x,_strx,_y,_stry)		\
  saxpy_(&(_n), &(_a), (_x), &(_strx), (_y), &(_stry))
#define blaswrap_daxpy(_n,_a,_x,_strx,_y,_stry)		\
  daxpy_(&(_n), &(_a), (_x), &(_strx), (_y), &(_stry))
#define blaswrap_caxpy(_n,_a,_x,_strx,_y,_stry)		\
  caxpy_(&(_n), (_a), (_x), &(_strx), (_y), &(_stry))
#define blaswrap_zaxpy(_n,_a,_x,_strx,_y,_stry)		\
  zaxpy_(&(_n), (_a), (_x), &(_strx), (_y), &(_stry))

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

#define blaswrap_dgbmv(_t,_m,_n,_kl,_ku,_al,_A,_lda,_x,_incx,_bt,_y,_incy) \
  do {									\
  if ((_t) == TRUE ) {							\
    g_assert_not_reached() ; /*untested code*/				\
  } else {								\
    dgbmv_("T", &(_m), &(_n), &(_kl), &(_ku), &(_al), (_A), &(_lda),	\
	   (_x), &(_incx), &(_bt), (_y), &(_incy)) ;			\
  }									\
  } while (0) 

/* y := al*A*x + bt*y, A symmetric */

/*upper n al A lda x strx bt y stry*/
#define blaswrap_dsymv(_u,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)		\
  do {									\
    if ( (_u) ) {							\
      dsymv_("U",&(_n),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    } else {								\
      dsymv_("L",&(_n),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    }									\
  } while (0)

#define blaswrap_ssymv(_u,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)		\
  do {									\
    if ( (_u) ) {							\
      ssymv_("U",&(_n),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    } else {								\
      ssymv_("L",&(_n),&(_al),(_A),&(_lda),(_x),&(_incx),		\
	     &(_bt),(_y),&(_incy)) ;					\
    }									\
  } while (0)

#define blaswrap_zsymv(_u,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)		\
  do {									\
    if ( (_u) ) {							\
      zsymv_("U",&(_n),(_al),(_A),&(_lda),(_x),&(_incx),		\
	     (_bt),(_y),&(_incy)) ;					\
    } else {								\
      zsymv_("L",&(_n),(_al),(_A),&(_lda),(_x),&(_incx),		\
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
      dgemm_("N","T", &(_n),&(_m),&(_k),&(_al),(_B),&(_ldb),		\
	   (_A),&(_lda),&(_bt),(_C),&(_ldc)) ;				\
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
#define blaswrap_zrotg(_a,_b,_c,_s)					\
  zrotg_((_a),(_b),(_c),(_s)) ;						

/*triangular system solves*/
#define blaswrap_dtrtrs(_upper,_trans,_diag,_n,_nrhs,_a,_lda,_b,_ldb,_info) \
  do {									\
    char _upstr[1], _tstr[1], _dstr[1] ;				\
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
#define blaswrap_ztrtrs(_upper,_trans,_diag,_n,_nrhs,_a,_lda,_b,_ldb,_info) \
  do {									\
    char _upstr[1], _tstr[1], _dstr[1] ;				\
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
    ztrtrs_(_upstr,_tstr,_dstr,&(_n),&(_nrhs),(_a),&(_lda),(_b),&(_ldb),(_info)) ; \
  }  while (0) 

#endif /*BLASWRAP_H_INCLUDED*/
