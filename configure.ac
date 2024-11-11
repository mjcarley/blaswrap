AC_INIT([blaswrap],[0.1])
AC_CONFIG_SRCDIR([src/blaswrap.h])
AC_CONFIG_MACRO_DIR([m4])

BLASWRAP_MAJOR_VERSION=0
BLASWRAP_MINOR_VERSION=1
BLASWRAP_MICRO_VERSION=1
BLASWRAP_INTERFACE_AGE=0
BLASWRAP_BINARY_AGE=0
BLASWRAP_VERSION=$BLASWRAP_MAJOR_VERSION.$BLASWRAP_MINOR_VERSION.$BLASWRAP_MICRO_VERSION

AC_SUBST(BLASWRAP_MAJOR_VERSION)
AC_SUBST(BLASWRAP_MINOR_VERSION)
AC_SUBST(BLASWRAP_MICRO_VERSION)
AC_SUBST(BLASWRAP_VERSION)

# libtool versioning
LT_RELEASE=$BLASWRAP_MAJOR_VERSION.$BLASWRAP_MINOR_VERSION
LT_CURRENT=`expr $BLASWRAP_MICRO_VERSION - $BLASWRAP_INTERFACE_AGE`
LT_REVISION=$BLASWRAP_INTERFACE_AGE
LT_AGE=`expr $BLASWRAP_BINARY_AGE - $BLASWRAP_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$BLASWRAP_VERSION
PACKAGE=blaswrap

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

##$PACKAGE, $VERSION, no-define)

# Specify a configuration file
AC_CONFIG_HEADERS([config.h])

AC_DEFINE_UNQUOTED(BLASWRAP_MAJOR_VERSION, $BLASWRAP_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(BLASWRAP_MINOR_VERSION, $BLASWRAP_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(BLASWRAP_MICRO_VERSION, $BLASWRAP_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(BLASWRAP_INTERFACE_AGE, $BLASWRAP_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(BLASWRAP_BINARY_AGE, $BLASWRAP_BINARY_AGE, [Binary age])

dnl Initialize libtool
LT_INIT

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_SEARCH_LIBS([strerror],[cposix])
AC_HEADER_STDC

PKG_CHECK_MODULES([GLIB], [glib-2.0])

AC_CHECK_LIB(m, cos)
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

ACX_BLAS([], [
  echo \#################################################################
  echo \#
  echo \# No BLAS located. If you do not have a BLAS installed on your
  echo \# system, you should install one. A suitable option is
  echo \# OpenBLAS, which is available from
  echo \#  
  echo \# https://www.openblas.net/
  echo \#
  echo \# If you do have a BLAS installed but it is not being located,
  echo \# you may need to specify the a link option with the LDFLAGS
  echo \# environment variable and/or use the --with-blas= option
  echo \# to specify the name of the BLAS library to be used.
  echo \#################################################################
  exit
])
ACX_LAPACK([
  echo LAPACK found at $LAPACK_LIBS
	], [
  echo \#################################################################
  echo \#
  echo \# LAPACK not located. If LAPACK is installed on your system but
  echo \# is not being located by the configure script, you may need to
  echo \# specify its location using the --with-lapack= option
  echo \#################################################################
  exit
])

CFLAGS="$CFLAGS $GLIB_CFLAGS"
LIBS="$GLIB_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_SUBST(CFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([
blaswrap.pc
Makefile
src/Makefile
src/blaswrap-config
tests/Makefile])

AC_OUTPUT

