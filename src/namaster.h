#ifndef _NAMASTER_
#define _NAMASTER_
#include "define.h"

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

#endif //_NAMASTER_
