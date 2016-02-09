#include "common_mvqe.h"

flouble dotprod(int n,flouble *v1,flouble *v2)
{
  int ii;
  double out=0;
  for(ii=0;ii<n;ii++)
    out+=v1[ii]*v2[ii];
  return (flouble)out;
}

void action_covar(flouble *map_in,flouble *map_out,ParamMVQE *par)
{
  int ii;

  memset(par->map_dum,0,par->npix*sizeof(flouble));
  for(ii=0;ii<par->npix_seen;ii++)
    par->map_dum[par->ipix_seen[ii]]=map_in[ii];

  he_map2alm(par->nside,par->lmax,1,&(par->map_dum),&(par->alm_dum));
  he_alter_alm(par->lmax,0,par->alm_dum,par->cl_proposal);
  he_alm2map(par->nside,par->lmax,1,&(par->map_dum),&(par->alm_dum));

  for(ii=0;ii<par->npix_seen;ii++)
    map_out[ii]=par->map_dum[par->ipix_seen[ii]];
}

void action_precond(flouble *map_in,flouble *map_out,ParamMVQE *par)
{
  int ii;

  memset(par->map_dum,0,par->npix*sizeof(flouble));
  for(ii=0;ii<par->npix_seen;ii++)
    par->map_dum[par->ipix_seen[ii]]=map_in[ii];

  he_map2alm(par->nside,par->lmax,1,&(par->map_dum),&(par->alm_dum));
  he_alter_alm(par->lmax,0,par->alm_dum,par->cl_precond);
  he_alm2map(par->nside,par->lmax,1,&(par->map_dum),&(par->alm_dum));

  for(ii=0;ii<par->npix_seen;ii++)
    map_out[ii]=par->map_dum[par->ipix_seen[ii]];
}

void callback_cg(flouble *mp,CGstate *cgs,ParamMVQE *par)
{
  return;
}
