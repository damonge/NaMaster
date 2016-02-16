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
  printf("   - CG tolerance = %lE\n",par->tol);
  printf("   - CG max. iter. = %d\n",par->max_iter);
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
  par->pixsize=-1.;
  par->npix_seen=-1;
  par->ipix_seen=NULL;
  par->lmax=-1;

  par->map_sht=NULL;
  par->alm_sht=NULL;

  par->map_v=NULL;
  par->alm_v=NULL;
  par->alm_vcm1=NULL;
  par->alm_v2=NULL;

  sprintf(par->fnameBins,"default"); 
  par->bins=NULL;
  par->nbins=-1;

  sprintf(par->fnameProposal,"default"); 
  par->lmax_proposal=-1;
  par->cl_proposal=NULL;
  par->sigma2_noise=-1;

  par->lcut_precond=-1;
  par->cl_precond=NULL;
  par->sigma2_precond=-1;
  par->tol=-1;
  par->max_iter=0;

  par->seed=1234;
  par->rngs=NULL;

  return par;
}

void free_param_mvqe(ParamMVQE *par)
{
  int ii;
  if(par->nside>0) {
    free(par->map);
    free(par->mask);
    free(par->map_sht);
    free(par->map_v);
  }
  if(par->lmax>0) {
    free(par->alm_sht);
    free(par->alm_v);
    free(par->alm_vcm1);
    free(par->alm_v2);
  }
  if(par->npix_seen>0)
    free(par->ipix_seen);
  if(par->nbins>0)
    free(par->bins);
  if(par->lmax_proposal>0)
    free(par->cl_proposal);
  if(par->lcut_precond>0)
    free(par->cl_precond);

  for(ii=0;ii<omp_get_max_threads();ii++)
    end_rng(par->rngs[ii]);
  free(par->rngs);

  free(par);
}

ParamMVQE *read_params(char *fname)
{
  FILE *fi;
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
    else if(!strcmp(s1,"proposal_noise="))
      par->sigma2_noise=atof(s2);
    else if(!strcmp(s1,"output_prefix="))
      sprintf(par->prefixOut,"%s",s2);
    else if(!strcmp(s1,"lcut_preconditioner="))
      par->lcut_precond=atoi(s2);
    else if(!strcmp(s1,"cg_tolerance="))
      par->tol=atof(s2);
    else if(!strcmp(s1,"cg_max_iter="))
      par->max_iter=atoi(s2);
    else
      fprintf(stderr,"MVQE: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  //Initialize rng
  par->rngs=my_malloc(omp_get_max_threads()*sizeof(Rng *));
  for(ii=0;ii<omp_get_max_threads();ii++)
    par->rngs[ii]=init_rng(par->seed+ii);

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
  par->map_sht=my_malloc(par->npix*sizeof(flouble));
  par->map_v=my_malloc(par->npix*sizeof(flouble));

  //Read bins
  fi=my_fopen(par->fnameBins,"r");
  par->nbins=my_linecount(fi)-1; rewind(fi);
  par->bins=my_malloc((par->nbins+1)*sizeof(int));
  for(ii=0;ii<=par->nbins;ii++) {
    int stat=fscanf(fi,"%d",&(par->bins[ii]));
    if(stat!=1)
      report_error(1,"Error reading bins file\n");
  }
  fclose(fi);

  //Read proposal power spectrum
  fi=my_fopen(par->fnameProposal,"r");
  par->lmax_proposal=my_linecount(fi)-1; rewind(fi);
  par->cl_proposal=my_malloc((par->lmax_proposal+1)*sizeof(flouble));
  for(ii=0;ii<=par->lmax_proposal;ii++) {
    int ell;
#ifdef _SPREC
    int stat=fscanf(fi,"%d %f",&ell,&(par->cl_proposal[ii]));
#else //_SPREC
    int stat=fscanf(fi,"%d %lf",&ell,&(par->cl_proposal[ii]));
#endif //_SPREC
    if((stat!=2) || (ell!=ii))
      report_error(1,"Error reading proposal covariance\n");
    par->cl_proposal[ii]/=par->pixsize;
  }
  fclose(fi);
  par->sigma2_noise/=par->pixsize;

  //Choose lmax
  int lmax_1=3*par->nside-1;
  int lmax_2=3*par->nside-1;//par->lmax_proposal;
  int lmax_3=3*par->nside-1;//par->bins[par->nbins];
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
  par->alm_sht=my_malloc(he_nalms(par->lmax)*sizeof(fcomplex));
  par->alm_v=my_malloc(he_nalms(par->lmax)*sizeof(fcomplex));
  par->alm_vcm1=my_malloc(he_nalms(par->lmax)*sizeof(fcomplex));
  par->alm_v2=my_malloc(he_nalms(par->lmax)*sizeof(fcomplex));

  //Prepare preconditioner power spectrum
  par->cl_precond=my_malloc(par->lmax*sizeof(flouble));
  for(ii=0;ii<par->lmax;ii++) {
    if(ii<par->lcut_precond)
      par->cl_precond[ii]=1./par->cl_proposal[ii]-1./par->cl_proposal[par->lcut_precond];
    else
      par->cl_precond[ii]=0;
  }
  par->sigma2_precond=1./par->cl_proposal[par->lcut_precond];

  //Print parameters
  param_mvqe_print(par);

  return par;
}

static flouble dotprod(int n,flouble *v1,flouble *v2)
{
  double out=0;

#pragma omp parallel default(none)		\
  shared(n,v1,v2,out)
  {
    int ii;
    double out_thr=0;
    
#pragma omp for
    for(ii=0;ii<n;ii++) {
      out_thr+=v1[ii]*v2[ii];
    } //end omp for
    
#pragma omp atomic
    out+=out_thr;
  } //end omp parallel

  return (flouble)out;
}

static void add_vscal(int n,flouble *v1,flouble alpha,flouble *v2,flouble *vout)
{
#pragma omp parallel default(none)\
  shared(n,v1,alpha,v2,vout)
  {
    int ii;

#pragma omp for
    for(ii=0;ii<n;ii++) {
      vout[ii]=v1[ii]+alpha*v2[ii];
    } //end omp for
  } //end omp parallel
}

static void action_covar(flouble *map_in,flouble *map_out,void *params)
{
  int ii;
  ParamMVQE *par=(ParamMVQE *)params;

  memset(par->map_sht,0,par->npix*sizeof(flouble));

#pragma omp parallel for
  for(ii=0;ii<par->npix_seen;ii++) {
    par->map_sht[par->ipix_seen[ii]]=map_in[ii];
  }

  he_map2alm(par->nside,par->lmax,1,&(par->map_sht),&(par->alm_sht));
  he_alter_alm(par->lmax,0,par->alm_sht,par->alm_sht,par->cl_proposal);
  he_alm2map(par->nside,par->lmax,1,&(par->map_sht),&(par->alm_sht));

#pragma omp parallel for
  for(ii=0;ii<par->npix_seen;ii++) {
    map_out[ii]=par->map_sht[par->ipix_seen[ii]]+par->sigma2_noise*map_in[ii];
  }
}

static void action_precond(flouble *map_in,flouble *map_out,void *params)
{
  int ii;
  ParamMVQE *par=(ParamMVQE *)params;

  memset(par->map_sht,0,par->npix*sizeof(flouble));
#pragma omp parallel for
  for(ii=0;ii<par->npix_seen;ii++) {
    par->map_sht[par->ipix_seen[ii]]=map_in[ii];
  }

  he_map2alm(par->nside,par->lcut_precond,1,&(par->map_sht),&(par->alm_sht));
  he_alter_alm(par->lcut_precond,0,par->alm_sht,par->alm_sht,par->cl_precond);
  he_alm2map(par->nside,par->lcut_precond,1,&(par->map_sht),&(par->alm_sht));

#pragma omp parallel for
  for(ii=0;ii<par->npix_seen;ii++) {
    map_out[ii]=par->map_sht[par->ipix_seen[ii]]+par->sigma2_precond*map_in[ii];
  }
}

//void callback_cg(flouble *mp,CGstate *cgs,ParamMVQE *par)
//{
//  return;
//}

void invert_covar(flouble *map_in,flouble *map_out,ParamMVQE *par)
{
  int ii,converged;
  flouble *pix_in=my_malloc(par->npix_seen*sizeof(flouble));
  flouble *pix_out=my_malloc(par->npix_seen*sizeof(flouble));
  flouble *pix_x0=my_malloc(par->npix_seen*sizeof(flouble));
  ParamCG *cg=cg_param_new(par->npix_seen);
  cg->action_A=&action_covar;
  cg->action_M=&action_precond;
  cg->dotpr=&dotprod;
  cg->add_vscal=&add_vscal;

#pragma omp parallel for
  for(ii=0;ii<par->npix_seen;ii++) {
    pix_in[ii]=map_in[par->ipix_seen[ii]];
  }
  action_precond(pix_in,pix_x0,par);

  converged=cg_solve(cg,pix_in,pix_x0,pix_out,par->tol,par->max_iter,par);
  if(converged==0)
    report_error(0,"CG did not converge after %d iterations\n",par->max_iter);

  memset(map_out,0,par->npix*sizeof(flouble));
#pragma omp parallel for
  for(ii=0;ii<par->npix_seen;ii++) {
    map_out[par->ipix_seen[ii]]=pix_out[ii];
  }

  free(pix_x0);
  free(pix_in);
  free(pix_out);
  cg_param_free(cg); 
}

#define SQRT12 3.46410161514
static void create_unit_variance(ParamMVQE *par,flouble *v)
{
  memset(v,0,par->npix*sizeof(flouble));

#pragma omp parallel default(none)\
  shared(par,v) 
  {
    int ii;
    Rng *rng=par->rngs[omp_get_thread_num()];

#pragma omp for
    for(ii=0;ii<par->npix_seen;ii++) {
      //      v[par->ipix_seen[ii]]=rand_gauss(rng);
      //      v[par->ipix_seen[ii]]=SQRT12*(rand_real01(rng)-0.5);
      v[par->ipix_seen[ii]]=(flouble)(-1+2*((int)(rand_ulong(rng)%2)));
    }
  }
}

static void fisher_single(ParamMVQE *par,flouble *fisher)
{
  int ib1,ib2,ii;
  flouble *transfer=my_malloc((par->lmax+1)*sizeof(flouble));
  flouble *cl=my_malloc((par->lmax+1)*sizeof(flouble));

  //alm_v=SHT(~v)
  he_map2alm(par->nside,par->lmax,1,&(par->map_v),&(par->alm_v));
  //v=C^-1*v 
  invert_covar(par->map_v,par->map_v,par);
  //alm_vcm1=SHT(~(C^-1*v))
  he_map2alm(par->nside,par->lmax,1,&(par->map_v),&(par->alm_vcm1));
  for(ib1=0;ib1<par->nbins;ib1++) {
    memset(transfer,0,(par->lmax+1)*sizeof(flouble));
    for(ii=par->bins[ib1];ii<par->bins[ib1+1];ii++)
      transfer[ii]=1./par->pixsize;
    //alm_v2=transfer*SHT(~v)/pixsize
    he_alter_alm(par->lmax,0,par->alm_v,par->alm_v2,transfer);
    //v=SHT^-1(transfer*SHT(~v)/pixsize)
    he_alm2map(par->nside,par->lmax,1,&(par->map_v),&(par->alm_v2));
    //v=C^-1*SHT^-1(transfer*SHT(~v)/pixsize)
    invert_covar(par->map_v,par->map_v,par);
    //alm_v2=SHT(C^-1*SHT^-1(transfer*SHT(~v)/pixsize))
    he_map2alm(par->nside,par->lmax,1,&(par->map_v),&(par->alm_v2));
    he_alm2cl(&(par->alm_vcm1),&(par->alm_v2),1,1,0,0,&cl,par->lmax);
    for(ib2=ib1;ib2<par->nbins;ib2++) {
      flouble fish_here=0;
      for(ii=par->bins[ib2];ii<par->bins[ib2+1];ii++)
	fish_here+=(2*ii+1.)*cl[ii];
      fish_here/=(par->pixsize*par->pixsize);
      fisher[ib2+par->nbins*ib1]+=fish_here;
      if(ib2!=ib1)
	fisher[ib1+par->nbins*ib2]+=fish_here;
    }
  }

  free(transfer);
  free(cl);
}

#define NAVG 10
void fisher_avg(ParamMVQE *par,flouble *fisher)
{
  int isam;
  memset(fisher,0,par->nbins*par->nbins*sizeof(flouble));
  
  for(isam=0;isam<NAVG;isam++) {
    int jj;
    flouble mean_here=0;
    create_unit_variance(par,par->map_v);
    fisher_single(par,fisher);
    for(jj=0;jj<par->nbins*par->nbins;jj++)
      mean_here+=fisher[jj]/(isam+1);
    printf("Sample %d, mean %lE\n",isam,mean_here);
  }
 
  for(isam=0;isam<par->nbins*par->nbins;isam++)
    fisher[isam]/=NAVG;
}
