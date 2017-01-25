#ifndef _COMMON_
#define _COMMON_

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

typedef struct {
  int n_bands;
  int *nell_list;
  int **ell_list;
  flouble **w_list;
} BinSchm;

//Defined in common.c
int my_linecount(FILE *f);
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
BinSchm *create_bins(int nlb,int nside);
BinSchm *read_bins(char *fname,int nside);
void free_bins(BinSchm *bin);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t my_fread(void *ptr,size_t size,size_t count,FILE *stream);

//Defined in master.c
void read_coupling_matrix(char *fname_in,int nbins_in,
			  gsl_matrix **coupling_matrix_b_out,
			  gsl_permutation **perm_out,
			  int pol1,int pol2);
void compute_coupling_matrix(flouble *cl_mask,long nside_in,int lmax_in,BinSchm *bins,
			     gsl_matrix **coupling_matrix_b_out,
			     gsl_permutation **perm_out,
			     char *write_matrix,char *write_matrix_b,
			     int pol1,int pol2);
flouble **decouple_cl_l(flouble **cl_in,flouble **cl_noise_in,
			int n_cl,BinSchm *bins,
			gsl_matrix *coupling_matrix_b,gsl_permutation *perm);

//Defined in healpix_extra.c
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,int pol,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,int pol,flouble **maps,fcomplex **alms);
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
void he_get_file_params(char *fname,long *nside,int *nfields,int *isnest);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
long *he_query_strip(long nside,double theta1,double theta2,
		     long *npix_strip);
void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest);
double *he_generate_beam_window(int lmax,double fwhm_amin);
void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alm_in,fcomplex *alm_out,double *window);
void he_alm2cl(fcomplex **alms_1,fcomplex **alms_2,
	       int nmaps_1,int nmaps_2,
	       int pol_1,int pol_2,
	       flouble **cls,int lmax);
void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,int pol_1,int pol_2,
		flouble **cls,int nside,int lmax);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
//flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed);

#endif //_COMMON_
