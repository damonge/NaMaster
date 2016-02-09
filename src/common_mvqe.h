#ifndef _COMMON_MVQE_
#define _COMMON_MVQE_

#include "common.h"

typedef struct {
  char prefixOut[256];

  char fnameMap[256];
  flouble *map;
  char fnameMask[256];
  flouble *mask; 
  int nside;
  int npix;
  int npix_seen;
  int lmax;
  int *ipix_seen;

  char fnameBins[256];
  int *bins;
  int nbins;

  char fnameProposal[256];
  int lmax_proposal;
  flouble *cl_proposal;
} ParamMVQE;

//Defined in io_mvqe.h
void free_param_mvqe(ParamMVQE *par);
ParamMVQE *read_params(char *fname);

#endif //_COMMON_MVQE_

