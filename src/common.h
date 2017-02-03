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

//Defined in utils.c
int my_linecount(FILE *f);
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t my_fread(void *ptr,size_t size,size_t count,FILE *stream);

//Defined in field.c
typedef struct {
  long nside;
  long npix;
  int lmax;
  flouble *mask;
  int pol;
  int nmaps;
  flouble **maps;
  int ntemp;
  flouble ***temp;
  fcomplex ***a_temp;
  gsl_matrix *matrix_M;
} Field;
void field_free(Field *fl);
Field *field_alloc(long nside,flouble *mask,int pol,flouble **maps,int ntemp,flouble ***temp);
Field *field_read(char *fname_mask,char *fname_maps,char *fname_temp,int pol);

//Defined in bins.c
typedef struct {
  int n_bands;
  int *nell_list;
  int **ell_list;
  flouble **w_list;
} BinSchm;
BinSchm *bins_create(int nlb,int nside);
BinSchm *bins_read(char *fname,int nside);
void bins_free(BinSchm *bin);

//Defined in master.c
typedef struct {
  int lmax;
  int ncls;
  flouble *pcl_masks;
  flouble **coupling_matrix_unbinned;
  BinSchm *bin;
  gsl_matrix *coupling_matrix_binned;
  gsl_permutation *coupling_matrix_perm;
} MasterWorkspace;
MasterWorkspace *compute_coupling_matrix(Field *fl1,Field *fl2,BinSchm *bin); //
void write_master_workspace(MasterWorkspace *w,char *fname); //
MasterWorkspace *read_master_workspace(char *fname); //
void master_workspace_free(MasterWorkspace *w); //
void compute_deprojection_bias(Field *fl1,Field *fl2,flouble **cl_proposal,flouble **cl_bias);
void decouple_cl_l(MasterWorkspace *w,flouble **cl_in,flouble **cl_noise_in,flouble **cl_bias,flouble **cl_out);
MasterWorkspace *compute_power_spectra(Field *fl1,Field *fl2,BinSchm *bin,
				       flouble **cl_noise,flouble **cl_proposal,flouble **cl_out);

//Defined in healpix_extra.c
#define HE_NITER_DEFAULT 3
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,int pol,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,int pol,flouble **maps,fcomplex **alms,int niter);
void he_alm2cl(fcomplex **alms_1,fcomplex **alms_2,int pol_1,int pol_2,flouble **cls,int lmax);
void he_anafast(flouble **maps_1,flouble **maps_2,int pol_1,int pol_2,flouble **cls,int nside,int lmax);
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
void he_get_file_params(char *fname,long *nside,int *nfields,int *isnest);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
void he_query_strip(long nside,double theta1,double theta2,int *pixlist,long *npix_strip);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
void he_in_ring(int nside,int iz,flouble phi0,flouble dphi,int *listir,int *nir);
void he_query_disc(int nside,double cth0,double phi,flouble radius,int *listtot,int *nlist,int inclusive);
void he_udgrade(flouble *map_in,long nside_in,flouble *map_out,long nside_out,int nest);
double *he_generate_beam_window(int lmax,double fwhm_amin);
void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alm_in,fcomplex *alm_out,double *window);
void he_map_product(int nside,flouble *mp1,flouble *mp2,flouble *mp_out);
flouble he_map_dot(int nside,flouble *mp1,flouble *mp2);
//flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed);
#ifdef _WITH_NEEDLET
#define HE_NBAND_NX 512
#define HE_NORM_FT 2.2522836206907617
#define HE_NL_INTPREC 1E-6
#define HE_NT_NSIDE_MIN 32
typedef struct {
  double b;
  double inv_b;
  gsl_spline *b_spline;
  gsl_interp_accel *b_intacc;
  int niter;
  int nside0;
  int jmax_min;
  int nj;
  int *nside_arr;
  int *lmax_arr;
  flouble **b_arr;
} HE_nt_param;
void he_nt_end(HE_nt_param *par);
HE_nt_param *he_nt_init(flouble b_nt,int nside0,int niter);
void he_free_needlet(HE_nt_param *par,int pol,flouble ***nt);
flouble ***he_alloc_needlet(HE_nt_param *par,int pol);
void he_nt_get_window(HE_nt_param *par,int j,flouble *b);
fcomplex **he_needlet2map(HE_nt_param *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int input_TEB,int output_TEB);
fcomplex **he_map2needlet(HE_nt_param *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int input_TEB,int output_TEB);
#endif //_WITH_NEEDLET
#endif //_COMMON_
