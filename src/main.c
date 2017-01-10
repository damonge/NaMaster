#include "common.h"
#include <unistd.h>

void run_master(char *fname_maps_1,char *fname_maps_2,
		char *fname_mask_1,char *fname_mask_2,
		int pol_1,int pol_2,
		char *fname_cl_noise,
		char *coupling_fname,
		char *coupling_b_fname,
		char *coupling_cov_fname,
		char *fname_out,
		char *fname_bins,
		int n_lbin)
{
  FILE *fi;
  long ip,nside_in,npix,nside_dum;
  int ii,nspec,lmax;
  flouble **cl_noise_bad_ub,**cl_maps_bad_ub,**dl_maps_good_b;
  flouble *mask_1,*mask_2,**maps_1,**maps_2;
  gsl_matrix *coupling_matrix_b;
  gsl_permutation *perm;
  int nmaps_1=1,nmaps_2=1;
  if(pol_1) nmaps_1=2;
  if(pol_2) nmaps_2=2;

  //Mask
  printf("Reading masks\n");
  mask_1=he_read_healpix_map(fname_mask_1,&nside_in,0);
  lmax=3*nside_in-1;
  npix=12*nside_in*nside_in;
  if(!strcmp(fname_mask_1,fname_mask_2))
    mask_2=mask_1;
  else {
    mask_2=he_read_healpix_map(fname_mask_2,&nside_dum,0);
    if(nside_dum!=nside_in)
      report_error(1,"Wrong nside %d\n",nside_dum);
  }

  //Binning
  BinSchm *bins;
  if(!strcmp(fname_bins,"none"))
    bins=create_bins(n_lbin,nside_in);
  else
    bins=read_bins(fname_bins,nside_in);

  //Maps
  printf("Reading maps\n");
  maps_1=my_malloc(nmaps_1*sizeof(flouble *));
  for(ii=0;ii<nmaps_1;ii++) {
    maps_1[ii]=he_read_healpix_map(fname_maps_1,&(nside_dum),ii);
    if(nside_dum!=nside_in)
      report_error(1,"Wrong nside %d\n",nside_dum);
    for(ip=0;ip<npix;ip++)
      maps_1[ii][ip]*=mask_1[ip];
  }
  if(!strcmp(fname_maps_1,fname_maps_2))
    maps_2=maps_1;
  else {
    maps_2=my_malloc(nmaps_2*sizeof(flouble *));
    for(ii=0;ii<nmaps_2;ii++) {
      maps_2[ii]=he_read_healpix_map(fname_maps_2,&(nside_dum),ii);
      if(nside_dum!=nside_in)
	report_error(1,"Wrong nside %d\n",nside_dum);
      for(ip=0;ip<npix;ip++)
	maps_2[ii][ip]*=mask_2[ip];
    }
  }

  //Allocate cl
  nspec=nmaps_1*nmaps_2;
  cl_noise_bad_ub=my_malloc(nspec*sizeof(flouble *));
  cl_maps_bad_ub=my_malloc(nspec*sizeof(flouble *));
  for(ii=0;ii<nspec;ii++) {
    cl_noise_bad_ub[ii]=my_calloc((lmax+1),sizeof(flouble));
    cl_maps_bad_ub[ii]=my_calloc((lmax+1),sizeof(flouble));
  }

  //Noise
  printf("Reading noise pseudo-cl\n");
  if(strcmp(fname_cl_noise,"none")) {
    fi=my_fopen(fname_cl_noise,"r");
    int nlin=my_linecount(fi); rewind(fi);
    if(nlin!=lmax+1)
      report_error(1,"Wrong number of multipoles for noise p.spec.\n");
    for(ii=0;ii<lmax+1;ii++) {
      int status,jj;
      flouble l;
      status=fscanf(fi,"%lf",&l);
      if(status!=1)
	report_error(1,"Error reading file %s\n",fname_cl_noise);
      for(jj=0;jj<nspec;jj++) {
	status=fscanf(fi,"%lf",&(cl_noise_bad_ub[jj][ii]));
	if(status!=1)
	  report_error(1,"Error reading file %s\n",fname_cl_noise);
      }
    }
    fclose(fi);
  }

  //Pseudo-cl
  printf("Computing data pseudo-cl\n");
  he_anafast(maps_1,maps_2,nmaps_1,nmaps_2,pol_1,pol_2,cl_maps_bad_ub,nside_in,lmax);

  //Coupling matrix
  if(access(coupling_b_fname,F_OK)!=-1) { //If file exists just read matrix
    printf("Reading coupling matrix\n");
    read_coupling_matrix(coupling_b_fname,bins->n_bands,
			 &coupling_matrix_b,&perm,pol_1,pol_2);
  }
  else { //Else, compute it
    flouble *cl_masks_bad_ub,*cl_mask_a_bad_ub,*cl_mask_b_bad_ub;
    cl_masks_bad_ub=my_malloc((lmax+1)*sizeof(flouble));
    printf("Computing mask pseudo-cl\n");
    he_anafast(&mask_1,&mask_2,1,1,0,0,&cl_masks_bad_ub,nside_in,lmax);
    if(strcmp(coupling_cov_fname,"none")) {
      flouble *mask11=my_malloc(npix*sizeof(flouble));
      flouble *mask12=my_malloc(npix*sizeof(flouble));
      flouble *mask22=my_malloc(npix*sizeof(flouble));
      for(ip=0;ip<npix;ip++) {
	mask11[ip]=mask_1[ip]*mask_1[ip];
	mask12[ip]=mask_1[ip]*mask_2[ip];
	mask22[ip]=mask_2[ip]*mask_2[ip];
      }
      printf("Computing mask pseudo-cl for covariance\n");
      he_anafast(&mask11,&mask22,1,1,0,0,&cl_mask_a_bad_ub,nside_in,lmax);
      he_anafast(&mask12,&mask12,1,1,0,0,&cl_mask_b_bad_ub,nside_in,lmax);
    }

    printf("Computing coupling matrix \n");
    compute_coupling_matrix(cl_masks_bad_ub,cl_mask_a_bad_ub,cl_mask_b_bad_ub,nside_in,lmax,bins,
			    &coupling_matrix_b,&perm,coupling_fname,coupling_b_fname,coupling_cov_fname,
			    pol_1,pol_2);
    free(cl_masks_bad_ub);
  }

  //Decouple cls
  printf("Decoupling cls\n");
  dl_maps_good_b=decouple_cl_l(cl_maps_bad_ub,cl_noise_bad_ub,nspec,
			       bins,coupling_matrix_b,perm);

  //Write output
  printf("Writing output\n");
  fi=my_fopen(fname_out,"w");
  for(ii=0;ii<bins->n_bands;ii++) {
    int jj;
    double l_here=0;
    for(jj=0;jj<bins->nell_list[ii];jj++)
      l_here+=bins->ell_list[ii][jj]*bins->w_list[ii][jj];
    fprintf(fi,"%.2lf ",l_here);
    for(jj=0;jj<nspec;jj++)
      fprintf(fi,"%lE ",dl_maps_good_b[jj][ii]);
    fprintf(fi,"\n");
  }
  fclose(fi);

  for(ii=0;ii<nspec;ii++) {
    free(cl_noise_bad_ub[ii]);
    free(cl_maps_bad_ub[ii]);
    free(dl_maps_good_b[ii]);
  }
  free(cl_noise_bad_ub);
  free(cl_maps_bad_ub);
  free(dl_maps_good_b);
  gsl_matrix_free(coupling_matrix_b);
  gsl_permutation_free(perm);
  for(ii=0;ii<nmaps_1;ii++)
    free(maps_1[ii]);
  free(maps_1);
  free(mask_1);
  if(strcmp(fname_maps_1,fname_maps_2)) {
    for(ii=0;ii<nmaps_2;ii++)
      free(maps_2[ii]);
    free(maps_2);
  }
  if(strcmp(fname_mask_1,fname_mask_2))
    free(mask_2);
}

int main(int argc,char **argv)
{
  int n_lbin=0,pol_1=0,pol_2=0;
  char fname_map_1[256]="none";
  char fname_map_2[256]="none";
  char fname_mask_1[256]="none";
  char fname_mask_2[256]="none";
  char fname_bins[256]="none";
  char fname_cl_noise[256]="none";
  char coupling_fname[256]="none";
  char coupling_cov_fname[256]="none";
  char coupling_b_fname[256]="none";
  char fname_out[256]="none";

  char **c;
  for(c=argv+1;*c;c++) {
    if(!strcmp(*c,"-map"))
      sprintf(fname_map_1,"%s",*++c);
    else if(!strcmp(*c,"-map_2"))
      sprintf(fname_map_2,"%s",*++c);
    else if(!strcmp(*c,"-mask"))
      sprintf(fname_mask_1,"%s",*++c);
    else if(!strcmp(*c,"-mask_2"))
      sprintf(fname_mask_2,"%s",*++c);
    else if(!strcmp(*c,"-pol"))
      pol_1=atoi(*++c);
    else if(!strcmp(*c,"-pol_2"))
      pol_2=atoi(*++c);
    else if(!strcmp(*c,"-cl_noise"))
      sprintf(fname_cl_noise,"%s",*++c);
    else if(!strcmp(*c,"-coupling"))
      sprintf(coupling_b_fname,"%s",*++c);
    else if(!strcmp(*c,"-coupling_unbinned"))
      sprintf(coupling_fname,"%s",*++c);
    else if(!strcmp(*c,"-coupling_covariance"))
      sprintf(coupling_cov_fname,"%s",*++c);
    else if(!strcmp(*c,"-out"))
      sprintf(fname_out,"%s",*++c);
    else if(!strcmp(*c,"-binning"))
      sprintf(fname_bins,"%s",*++c);
    else if(!strcmp(*c,"-nlb"))
      n_lbin=atoi(*++c);
    else if(!strcmp(*c,"-h")) {
      fprintf(stderr,"Usage: NaMaster -<opt-name> <option>\n");
      fprintf(stderr,"Options:\n");
      fprintf(stderr,"  -map      -> path to file containing map(s)\n");
      fprintf(stderr,"  -map_2    -> path to file containing 2nd map(s) (optional)\n");
      fprintf(stderr,"  -mask     -> path to file containing mask\n");
      fprintf(stderr,"  -mask_2   -> path to file containing mask for 2nd map(s) (optional)\n");
      fprintf(stderr,"  -pol      -> spin-0 (0) or spin-2 (1) input map(s)\n");
      fprintf(stderr,"  -pol_2    -> spin-0 (0) or spin-2 (1) 2nd input map(s)\n");
      fprintf(stderr,"  -cl_noise -> path to file containing noise Cl(s)\n");
      fprintf(stderr,"  -coupling -> path to file containing coupling matrix (optional)\n");
      fprintf(stderr,"  -coupling_unbinned -> path to file containing unbinned coupling matrix (optional)\n");
      fprintf(stderr,"               If non-existing, it will be computed and\n");
      fprintf(stderr,"               written there\n");
      fprintf(stderr,"  -out      -> output filename\n");
      fprintf(stderr,"  -binning  -> path to file containing binning scheme\n");
      fprintf(stderr,"  -nlb      -> number of ells per bin (used only if -binning isn't used)\n");
      fprintf(stderr,"  -h        -> this help\n\n");
      return 0;
    }
    else {
      fprintf(stderr,"Unknown option %s\n",*c);
      exit(1);
    }
  }

  if(!strcmp(fname_map_1,"none"))
    report_error(1,"Must provide map to correlate!\n");
  if(!strcmp(fname_mask_1,"none"))
    report_error(1,"Must provide mask\n");
  if(!strcmp(fname_out,"none"))
    report_error(1,"Must provide output filename\n");
  if(n_lbin<=0)
    report_error(1,"#ell per bin must be positive\n");

  if(!strcmp(fname_map_2,"none")) {
    sprintf(fname_map_2,"%s",fname_map_1);
    pol_2=pol_1;
  }
  if(!strcmp(fname_mask_2,"none"))
    sprintf(fname_mask_2,"%s",fname_mask_1);

  run_master(fname_map_1,fname_map_2,
	     fname_mask_1,fname_mask_2,
	     pol_1,pol_2,
	     fname_cl_noise,
	     coupling_fname,
	     coupling_b_fname,
	     coupling_cov_fname,
	     fname_out,fname_bins,n_lbin);

  return 0;
}
