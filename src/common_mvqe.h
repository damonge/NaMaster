#ifndef _COMMON_MVQE_
#define _COMMON_MVQE_

#include "common.h"

//Defined in cg.c
typedef struct {
  int n;
  flouble (*dotpr)(int *,flouble *,flouble *);
  void (*action_A)(flouble *,flouble *,void *);
  void (*action_E)(flouble *,flouble *,void *);
  flouble *xk;
  flouble *rk;
  flouble *zk;
  flouble *pk;
  flouble *vdum;
} ParamCG;
ParamCG *cg_param_new(int n);
void cg_param_free(ParamCG *cg);
void cg_solve(ParamCG *cg,flouble *b,flouble *x0,flouble *xf,flouble tol,int max_iter,void *params);

typedef struct {
  char prefixOut[256];

  char fnameMap[256];
  flouble *map;
  flouble *map_dum;
  char fnameMask[256];
  flouble *mask; 
  int nside;
  int npix;
  int npix_seen;
  int *ipix_seen;
  int lmax;
  fcomplex *alm_dum;

  char fnameBins[256];
  int *bins;
  int nbins;

  char fnameProposal[256];
  int lmax_proposal;
  flouble *cl_proposal;

  int lcut_precond;
  flouble *cl_precond;
} ParamMVQE;

//Defined in io_mvqe.h
void free_param_mvqe(ParamMVQE *par);
ParamMVQE *read_params(char *fname);
void invert_covar(flouble *map_in,flouble *map_out,ParamMVQE *par);

#endif //_COMMON_MVQE_

