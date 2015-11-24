#ifndef _COMMON_FGRM
#define _COMMON_FGRM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <dam_utils.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

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

//Defined in master.c
void read_coupling_matrix(char *fname_in,int nbins_in,
			  gsl_matrix **coupling_matrix_b_out,
			  gsl_permutation **perm_out);
void compute_coupling_matrix(flouble *cl_mask,int n_lbin,
			     long nside_in,int lmax_in,int nbins_in,
			     gsl_matrix **coupling_matrix_b_out,
			     gsl_permutation **perm_out,char *write_matrix);
flouble *decouple_cl_l(flouble *cl_in,flouble *cl_noise_in,int nbins,int n_lbin,
		       gsl_matrix *coupling_matrix_b,gsl_permutation *perm);

//Defined in healpix_extra.c
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
long *he_query_strip(long nside,double theta1,double theta2,
		     long *npix_strip);
void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest);
double *he_generate_beam_window(int lmax,double fwhm_amin);
void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alms,double *window);
void he_anafast(flouble **maps,flouble **cls,int nside,int lmax,int nmaps,int pol);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed);

#endif //_COMMON_FGRM
