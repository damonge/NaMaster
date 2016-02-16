#ifndef _COMMON_MVQE_
#define _COMMON_MVQE_

#include "common.h"

//Defined in cg.c
typedef struct {
  int n;
  flouble (*dotpr)(int ,flouble *,flouble *);
  void (*add_vscal)(int ,flouble *,flouble ,flouble *,flouble *);
  void (*action_A)(flouble *,flouble *,void *);
  void (*action_M)(flouble *,flouble *,void *);
  flouble *xk;
  flouble *rk;
  flouble *zk;
  flouble *pk;
  flouble *vdum;
} ParamCG;
ParamCG *cg_param_new(int n);
void cg_param_free(ParamCG *cg);
int cg_solve(ParamCG *cg,flouble *b,flouble *x0,flouble *xf,flouble tol,int max_iter,void *params);

//Defined in rng.c
#define RNG_NRAN 624
#define RNG_MRAN 397
#define RNG_MATRIX_A 0x9908b0df
#define RNG_UPPER_MASK 0x80000000UL
#define RNG_LOWER_MASK 0x7fffffffUL
typedef struct {
  unsigned long mt[RNG_NRAN];
  int mti;
  int calc_gauss;
  double u;
  double phi;
} Rng;
Rng *init_rng(unsigned long seed);
void end_rng(Rng *rng);
unsigned long rand_ulong(Rng *rng);
double rand_real01(Rng *rng);
double rand_gauss(Rng *rng);

typedef struct {
  char prefixOut[256];

  char fnameMap[256];
  flouble *map;
  char fnameMask[256];
  flouble *mask; 
  int nside;
  int npix;
  int npix_seen;
  int *ipix_seen;
  int lmax;
  flouble pixsize;
  flouble *map_sht;
  fcomplex *alm_sht;
  flouble *map_v;
  fcomplex *alm_v;
  fcomplex *alm_vcm1;
  fcomplex *alm_v2;

  char fnameBins[256];
  int *bins;
  int nbins;

  char fnameProposal[256];
  int lmax_proposal;
  flouble *cl_proposal;
  flouble sigma2_noise;

  int lcut_precond;
  flouble *cl_precond;
  flouble tol;
  int max_iter;
  flouble sigma2_precond;

  int seed;
  Rng **rngs;
} ParamMVQE;

//Defined in common_mvqe.c
void free_param_mvqe(ParamMVQE *par);
ParamMVQE *read_params(char *fname);
void invert_covar(flouble *map_in,flouble *map_out,ParamMVQE *par);
void fisher_avg(ParamMVQE *par,flouble *fisher);
//void create_unit_variance(ParamMVQE *par,flouble *v);

#endif //_COMMON_MVQE_

