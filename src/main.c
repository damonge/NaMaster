#include "common.h"

void run_master(Field *fl1,Field *fl2,
		char *fname_cl_noise,
		char *fname_cl_proposal,
		char *coupling_fname,
		char *coupling_b_fname,
		char *fname_out,
		char *fname_bins,
		int n_lbin)
{
  FILE *fi;
  int ii;
  flouble **cl_bias,**cl_noise_bad_ub,**cl_proposal,**cl_maps_bad_ub,**dl_maps_good_b;
  long nside_in=fl1->nside;
  int lmax=fl1->lmax;
  int nspec=fl1->nmaps*fl2->nmaps;
  gsl_matrix *coupling_matrix_b;
  gsl_permutation *perm;

  if(fl1->nside!=fl2->nside)
    report_error(1,"Can't correlate fields with different resolution\n");

  //Binning
  BinSchm *bins;
  if(!strcmp(fname_bins,"none"))
    bins=bins_create(n_lbin,nside_in);
  else
    bins=bins_read(fname_bins,nside_in);

  //Allocate cl
  cl_bias=my_malloc(nspec*sizeof(flouble *));
  cl_noise_bad_ub=my_malloc(nspec*sizeof(flouble *));
  cl_proposal=my_malloc(nspec*sizeof(flouble *));
  cl_maps_bad_ub=my_malloc(nspec*sizeof(flouble *));
  for(ii=0;ii<nspec;ii++) {
    cl_bias[ii]=my_calloc((lmax+1),sizeof(flouble));
    cl_noise_bad_ub[ii]=my_calloc((lmax+1),sizeof(flouble));
    cl_proposal[ii]=my_calloc((lmax+1),sizeof(flouble));
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

  //Read proposal cls
  if(strcmp(fname_cl_proposal,"none")) {
    fi=my_fopen(fname_cl_proposal,"r");
    int nlin=my_linecount(fi); rewind(fi);
    if(nlin!=lmax+1)
      report_error(1,"Wrong number of multipoles for noise p.spec.\n");
    for(ii=0;ii<lmax+1;ii++) {
      int status,jj;
      flouble l;
      status=fscanf(fi,"%lf",&l);
      if(status!=1)
	report_error(1,"Error reading file %s\n",fname_cl_proposal);
      for(jj=0;jj<nspec;jj++) {
	status=fscanf(fi,"%lf",&(cl_proposal[jj][ii]));
	if(status!=1)
	  report_error(1,"Error reading file %s\n",fname_cl_proposal);
      }
    }
    fclose(fi);
  }

  //Pseudo-cl
  printf("Computing data pseudo-cl\n");
  he_anafast(fl1->maps,fl2->maps,fl1->nmaps,fl2->nmaps,fl1->pol,fl1->pol,cl_maps_bad_ub,nside_in,lmax);

  //Compute deprojection bias
  printf("Computing deprojection bias\n");
  if((fl1->ntemp>0) || (fl2->ntemp>0))
    compute_deprojection_bias(fl1,fl2,cl_proposal,cl_bias);

  //Coupling matrix
  if(access(coupling_b_fname,F_OK)!=-1) { //If file exists just read matrix
    printf("Reading coupling matrix\n");
    read_coupling_matrix(coupling_b_fname,bins->n_bands,&coupling_matrix_b,&perm,fl1->nmaps*fl2->nmaps);
  }
  else { //Else, compute it
    printf("Computing coupling matrix \n");
    compute_coupling_matrix(fl1,fl2,bins,&coupling_matrix_b,&perm,coupling_fname,coupling_b_fname);
  }

  //Decouple cls
  printf("Decoupling cls\n");
  dl_maps_good_b=decouple_cl_l(cl_maps_bad_ub,cl_noise_bad_ub,cl_bias,
			       nspec,bins,coupling_matrix_b,perm);

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

  gsl_matrix_free(coupling_matrix_b);
  gsl_permutation_free(perm);
  bins_free(bins);

  for(ii=0;ii<nspec;ii++) {
    free(cl_bias[ii]);
    free(cl_noise_bad_ub[ii]);
    free(cl_proposal[ii]);
    free(cl_maps_bad_ub[ii]);
    free(dl_maps_good_b[ii]);
  }
  free(cl_bias);
  free(cl_proposal);
  free(cl_noise_bad_ub);
  free(cl_maps_bad_ub);
  free(dl_maps_good_b);
}

int main(int argc,char **argv)
{
  int n_lbin=1,pol_1=0,pol_2=0,is_auto=0;
  char fname_map_1[256]="none";
  char fname_map_2[256]="none";
  char fname_mask_1[256]="none";
  char fname_mask_2[256]="none";
  char fname_temp_1[256]="none";
  char fname_temp_2[256]="none";
  char fname_bins[256]="none";
  char fname_cl_noise[256]="none";
  char fname_cl_proposal[256]="none";
  char coupling_fname[256]="none";
  char coupling_b_fname[256]="none";
  char fname_out[256]="none";
  Field *fl1,*fl2;

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
    else if(!strcmp(*c,"-temp"))
      sprintf(fname_temp_1,"%s",*++c);
    else if(!strcmp(*c,"-temp_2"))
      sprintf(fname_temp_2,"%s",*++c);
    else if(!strcmp(*c,"-pol"))
      pol_1=atoi(*++c);
    else if(!strcmp(*c,"-pol_2"))
      pol_2=atoi(*++c);
    else if(!strcmp(*c,"-cl_noise"))
      sprintf(fname_cl_noise,"%s",*++c);
    else if(!strcmp(*c,"-cl_guess"))
      sprintf(fname_cl_proposal,"%s",*++c);
    else if(!strcmp(*c,"-coupling"))
      sprintf(coupling_b_fname,"%s",*++c);
    else if(!strcmp(*c,"-coupling_unbinned"))
      sprintf(coupling_fname,"%s",*++c);
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
      fprintf(stderr,"  -temp     -> path to file containing contaminant templates (optional)\n");
      fprintf(stderr,"  -temp_2   -> path to file containing contaminant templates\n");
      fprintf(stderr,"               for 2nd map(s) (optional)\n");
      fprintf(stderr,"  -pol      -> spin-0 (0) or spin-2 (1) input map(s)\n");
      fprintf(stderr,"  -pol_2    -> spin-0 (0) or spin-2 (1) 2nd input map(s)\n");
      fprintf(stderr,"  -cl_noise -> path to file containing noise Cl(s)\n");
      fprintf(stderr,"  -cl_guess -> path to file containing initial guess for the Cl(s)\n");
      fprintf(stderr,"  -coupling -> path to file containing coupling matrix (optional)\n");
      fprintf(stderr,"  -coupling_unbinned -> path to file containing unbinned coupling matrix (optional)\n");
      fprintf(stderr,"               If non-existing, it will be computed and written there\n");
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

  fl1=field_alloc(fname_mask_1,fname_map_1,fname_temp_1,pol_1);

  if(!strcmp(fname_map_2,"none")) {
    fl2=fl1;
    is_auto=1;
  }
  else {
    if(!strcmp(fname_mask_2,"none"))
      sprintf(fname_mask_2,"%s",fname_mask_1);
    if(!strcmp(fname_temp_2,"none"))
      sprintf(fname_temp_2,"%s",fname_temp_1);
    fl2=field_alloc(fname_mask_2,fname_map_2,fname_temp_2,pol_2);
  }

  run_master(fl1,fl2,
	     fname_cl_noise,
	     fname_cl_proposal,
	     coupling_fname,
	     coupling_b_fname,
	     fname_out,fname_bins,n_lbin);

  field_free(fl1);
  if(!is_auto)
    field_free(fl2);

  return 0;
}
