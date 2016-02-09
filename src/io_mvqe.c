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
  printf("   - lmax_proposal = %d\n",par->lmax_proposal);
  printf("   - Output prefix : %s\n",par->prefixOut);
  printf("\n");
}

static ParamMVQE *param_mvqe_new(void)
{
  ParamMVQE *par=my_malloc(sizeof(ParamMVQE));
  
  sprintf(par->prefixOut,"default"); 

  sprintf(par->fnameMap,"default"); 
  par->map=NULL;
  sprintf(par->fnameMask,"default");
  par->mask=NULL;
  par->nside=-1;
  par->npix=-1;
  par->npix_seen=-1;
  par->ipix_seen=NULL;

  sprintf(par->fnameBins,"default"); 
  par->bins=NULL;
  par->nbins=-1;

  sprintf(par->fnameProposal,"default"); 
  par->lmax_proposal=-1;
  par->cl_proposal=NULL;

  return par;
}

void free_param_mvqe(ParamMVQE *par)
{
  if(par->nside>0) {
    free(par->map);
    free(par->mask);
  }
  if(par->npix_seen>0)
    free(par->ipix_seen);
  if(par->nbins>0)
    free(par->bins);
  if(par->lmax_proposal>0)
    free(par->cl_proposal);

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

  //Print parameters
  param_mvqe_print(par);

  return par;
}
