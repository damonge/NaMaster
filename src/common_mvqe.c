#include "common_mvqe.h"

static void param_mvqe_print(ParamMVQE *par)
{
  printf(" * Read parameters :\n");
  printf("   - Input map : %s\n",par->fnameMap);
  printf("   - Input mask : %s\n",par->fnameMask);
  printf("   - Nside = %d\n",par->nside);
  printf("   - #seen pixels = %d/%d\n",par->npix_seen,par->npix);
  printf("   - Input bins : %s\n",par->fnameBins);
  printf("   - %d bins \n",par->nbins);
  printf("   - Input covariance : %s\n",par->fnameProposal);
  printf("   - lmax = %d\n",par->lmax);
  printf("   - lmax_proposal = %d\n",par->lmax_proposal);
  printf("   - lcut_precond = %d\n",par->lcut_precond);
  printf("   - Output prefix : %s\n",par->prefixOut);
  printf("\n");
}

static ParamMVQE *param_mvqe_new(void)
{
  ParamMVQE *par=my_malloc(sizeof(ParamMVQE));
  
  sprintf(par->prefixOut,"default"); 

  sprintf(par->fnameMap,"default"); 
  par->map=NULL;
  par->map_dum=NULL;
  sprintf(par->fnameMask,"default");
  par->mask=NULL;
  par->nside=-1;
  par->npix=-1;
  par->pixsize=-1.;
  par->npix_seen=-1;
  par->ipix_seen=NULL;
  par->lmax=-1;
  par->alm_dum=NULL;

  sprintf(par->fnameBins,"default"); 
  par->bins=NULL;
  par->nbins=-1;

  sprintf(par->fnameProposal,"default"); 
  par->lmax_proposal=-1;
  par->cl_proposal=NULL;

  par->lcut_precond=-1;
  par->cl_precond=NULL;

  return par;
}

void free_param_mvqe(ParamMVQE *par)
{
  if(par->nside>0) {
    free(par->map);
    free(par->mask);
    free(par->map_dum);
  }
  if(par->lmax>0)
    free(par->alm_dum);
  if(par->npix_seen>0)
    free(par->ipix_seen);
  if(par->nbins>0)
    free(par->bins);
  if(par->lmax_proposal>0)
    free(par->cl_proposal);
  if(par->lcut_precond>0)
    free(par->cl_precond);

  free(par);
}

ParamMVQE *read_params(char *fname)
{
  FILE *fi;
  char fname_in[256];
  int n_lin,ii;
  long nside_dum;
  ParamMVQE *par=param_mvqe_new();

  //Read parameters from file
  printf("*** Reading run parameters \n");
  fi=my_fopen(fname,"r");
  n_lin=my_linecount(fi); rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);

    if(!strcmp(s1,"input_map="))
      sprintf(par->fnameMap,"%s",s2);
    else if(!strcmp(s1,"input_mask="))
      sprintf(par->fnameMask,"%s",s2);
    else if(!strcmp(s1,"input_bins="))
      sprintf(par->fnameBins,"%s",s2);
    else if(!strcmp(s1,"input_proposal_cl="))
      sprintf(par->fnameProposal,"%s",s2);
    else if(!strcmp(s1,"output_prefix="))
      sprintf(par->prefixOut,"%s",s2);
    else if(!strcmp(s1,"lcut_preconditioner="))
      par->lcut_precond=atoi(s2);
   else
      fprintf(stderr,"MVQE: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  //Read map
  par->map=he_read_healpix_map(par->fnameMap,&nside_dum,0);
  par->nside=(int)(nside_dum);
  //Read mask
  par->mask=he_read_healpix_map(par->fnameMask,&nside_dum,0);
  if(nside_dum!=par->nside)
    report_error(1,"Mask and map resolutions don't match\n");
  par->npix=12*par->nside*par->nside;
  par->pixsize=4*M_PI/par->npix;
  //Get seen pixels
  par->npix_seen=0;
  for(ii=0;ii<par->npix;ii++) {
    if(par->mask[ii]>0.9)
      par->npix_seen++;
  }
  par->ipix_seen=my_malloc(par->npix_seen*sizeof(int));
  par->npix_seen=0;
  for(ii=0;ii<par->npix;ii++) {
    if(par->mask[ii]>0.9) {
      par->ipix_seen[par->npix_seen]=ii;
      par->npix_seen++;
    }
  }
  par->map_dum=my_malloc(par->npix*sizeof(flouble));

  //Read bins
  fi=my_fopen(par->fnameBins,"r");
  par->nbins=my_linecount(fi)-1;
  par->bins=my_malloc((par->nbins+1)*sizeof(int));
  for(ii=0;ii<=par->nbins;ii++) {
    int stat=fscanf(fi,"%d",&(par->bins[ii]));
    if(stat!=1)
      report_error(1,"Error reading bins file\n");
  }
  fclose(fi);

  //Read proposal power spectrum
  fi=my_fopen(par->fnameProposal,"r");
  par->lmax_proposal=my_linecount(fi)-1;
  par->cl_proposal=my_malloc((par->lmax_proposal+1)*sizeof(flouble));
  for(ii=0;ii<=par->lmax_proposal;ii++) {
    int ell;
#ifdef _SPREC
    int stat=fscanf(fi,"%d %f",ell,&(par->cl_proposal[ii]));
#else //_SPREC
    int stat=fscanf(fi,"%d %lf",ell,&(par->cl_proposal[ii]));
#endif //_SPREC
    if((stat!=2) || (ell!=ii))
      report_error(1,"Error reading proposal covariance\n");
    par->cl_proposal[ii]/=par->pixsize;
  }
  fclose(fi);

  //Choose lmax
  int lmax_1=3*par->nside-1;
  int lmax_2=par->lmax_proposal;
  int lmax_3=par->bins[par->nbins];
  par->lmax=lmax_1;
  if(lmax_2<par->lmax)
    par->lmax=lmax_2;
  if(lmax_3<par->lmax)
    par->lmax=lmax_3;
  par->lmax_proposal=par->lmax;
  ii=0;
  while(par->bins[ii+1]<=par->lmax)
    ii++;
  par->nbins=ii;
  par->alm_dum=my_malloc(he_nalms(par->lmax)*sizeof(fcomplex));

  //Prepare preconditioner power spectrum
  par->cl_precond=my_malloc(par->lmax*sizeof(flouble));
  for(ii=0;ii<par->lmax;ii++) {
    if(ii<par->lcut_precond)
      par->cl_precond[ii]=par->pixsize/par->cl_proposal[ii];
    else
      par->cl_precond[ii]=par->pixsize/par->cl_proposal[par->lcut_precond];
  }

  //Print parameters
  param_mvqe_print(par);

  return par;
}

static flouble dotprod(int n,flouble *v1,flouble *v2)
{
  int ii;
  double out=0;
  for(ii=0;ii<n;ii++)
    out+=v1[ii]*v2[ii];
  return (flouble)out;
}

static void action_covar(flouble *map_in,flouble *map_out,ParamMVQE *par)
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

static void action_precond(flouble *map_in,flouble *map_out,ParamMVQE *par)
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

//void callback_cg(flouble *mp,CGstate *cgs,ParamMVQE *par)
//{
//  return;
//}

void invert_covar(flouble *map_in,flouble *map_out,ParamMVQE *par)
{
  flouble *map_x0=my_malloc(par->npix_seen*sizeof(flouble));
  ParamCG *cg=cg_param_new(par->npix_seen);
  cg->action_A=&action_covar;
  cg->action_M=&action_precond;
  cg->dotpr=&dotprod;
  action_precond(map_in,map_x0,par);

  cg_solve(cg,map_in,map_x0,map_out,par->tol,par->max_iter,par);

  free(map_x0);
  cg_param_free(cg);
}
