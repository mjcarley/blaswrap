BLASWRAP is a (n incomplete) set of wrappers to ease the calling of
BLAS routines from C, mainly for use in my other numerical codes. It
defines a set of macros which rearrange function arguments so that
BLAS calls will work with C matrix ordering.

# Prerequisites

You will need to have a BLAS and LAPACK installed so that the
configuration script can find them.

# Installation

If you have downloaded the source from github or equivalent, you will
need the autotools suite. Generate the configure script with

`. autogen.sh`

To configure and install the code,

`./configure [OPTIONS]`

`make`

`make install`

For information on options to control configuration, including the
installation location and where to find any required libraries:

  `./configure --help`


