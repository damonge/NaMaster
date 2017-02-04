#ifndef _NM_DEFINE_
#define _NM_DEFINE_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#ifdef _WITH_NEEDLET
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#endif //_WITH_NEEDLET

#define COM_MAX(a,b)  (((a)>(b)) ? (a) : (b)) // maximum
#define COM_MIN(a,b)  (((a)<(b)) ? (a) : (b)) // minimum

#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

#ifdef _SPREC
typedef float flouble;
typedef float complex fcomplex;
#else //_SPREC
typedef double flouble;
typedef double complex fcomplex;
#endif //_SPREC

#endif //_NM_DEFINE_
