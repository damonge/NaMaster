#include "common.h"

void run_master(char *fname_maps_1,char *fname_maps_2,
		char *fname_mask_1,char *fname_mask_2,
		char *fname_temp_1,char *fname_temp_2,
		int pol_1,int pol_2,
		char *fname_cl_noise,
		char *fname_cl_proposal,
		char *coupling_fname,
		char *coupling_b_fname,
		char *fname_out,
		char *fname_bins,
		int n_lbin)
{
  FILE *fi;
  long ip,nside_in,npix,nside_dum;
  int ii,nspec,lmax;
  flouble **cl_noise_bad_ub,**cl_proposal,**cl_dum,**cl_maps_bad_ub,**dl_maps_good_b;
  flouble *mask_1,*mask_2,**maps_1,**maps_2,***temp_1,***temp_2;
  fcomplex ***a_temp_1,***a_temp_2;
  gsl_matrix *coupling_matrix_b;
  gsl_permutation *perm;
  int ntemp_1,ntemp_2;
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

  //Contaminants
  printf("Reading contaminant templates\n");
  if(strcmp(fname_temp_1,"none")) {
    int ncols,isnest,itemp;
    he_get_file_params(fname_temp_1,&nside_dum,&ncols,&isnest);
    if(nside_dum!=nside_in)
      report_error(1,"Wrong nside %d\n",nside_dum);
    if((ncols==0) || (ncols%nmaps_1!=0))
      report_error(1,"Not enough templates in file %s\n",fname_temp_1);
    ntemp_1=ncols/nmaps_1;
    temp_1=my_malloc(ntemp_1*sizeof(flouble **));
    a_temp_1=my_malloc(ntemp_1*sizeof(fcomplex **));
    for(itemp=0;itemp<ntemp_1;itemp++) {
      int imap;
      temp_1[itemp]=my_malloc(nmaps_1*sizeof(flouble *));
      a_temp_1[itemp]=my_malloc(nmaps_1*sizeof(fcomplex *));
      for(imap=0;imap<nmaps_1;imap++) {
	a_temp_1[itemp][imap]=my_malloc(he_nalms(lmax)*sizeof(fcomplex));
	temp_1[itemp][imap]=he_read_healpix_map(fname_temp_1,&(nside_dum),itemp*nmaps_1+imap);
	for(ip=0;ip<npix;ip++)
	  temp_1[itemp][imap][ip]*=mask_1[ip];
	he_map2alm(nside_in,lmax,1,pol_1,temp_1[itemp],a_temp_1[itemp]);
      }
    }
  }
  else
    ntemp_1=0;
  if(!strcmp(fname_temp_1,fname_temp_2)) {
    a_temp_2=a_temp_1;
    temp_2=temp_1;
    ntemp_2=ntemp_1;
  }
  else {
    if(strcmp(fname_temp_2,"none")) {
      int ncols,isnest,itemp;
      he_get_file_params(fname_temp_2,&nside_dum,&ncols,&isnest);
      if(nside_dum!=nside_in)
	report_error(1,"Wrong nside %d\n",nside_dum);
      if((ncols==0) || (ncols%nmaps_2!=0))
	report_error(1,"Not enough templates in file %s\n",fname_temp_2);
      ntemp_2=ncols/nmaps_2;
      temp_2=my_malloc(ntemp_2*sizeof(flouble **));
      a_temp_2=my_malloc(ntemp_2*sizeof(fcomplex **));
      for(itemp=0;itemp<ntemp_2;itemp++) {
	int imap;
	temp_2[itemp]=my_malloc(nmaps_2*sizeof(flouble *));
	a_temp_2[itemp]=my_malloc(nmaps_2*sizeof(fcomplex *));
	for(imap=0;imap<nmaps_2;imap++) {
	  a_temp_2[itemp][imap]=my_malloc(he_nalms(lmax)*sizeof(fcomplex));
	  temp_2[itemp][imap]=he_read_healpix_map(fname_temp_2,&(nside_dum),itemp*nmaps_2+imap);
	  for(ip=0;ip<npix;ip++)
	    temp_2[itemp][imap][ip]*=mask_2[ip];
	  he_map2alm(nside_in,lmax,1,pol_2,temp_2[itemp],a_temp_2[itemp]);
	}
      }
    }
    else
      ntemp_2=0;
  }
  printf("%d and %d templates maps 1 and 2\n",ntemp_1,ntemp_2);

  //Binning
  BinSchm *bins;
  if(!strcmp(fname_bins,"none"))
    bins=create_bins(n_lbin,nside_in);
  else
    bins=read_bins(fname_bins,nside_in);

  //Allocate cl
  nspec=nmaps_1*nmaps_2;
  cl_noise_bad_ub=my_malloc(nspec*sizeof(flouble *));
  cl_proposal=my_malloc(nspec*sizeof(flouble *));
  cl_dum=my_malloc(nspec*sizeof(flouble *));
  cl_maps_bad_ub=my_malloc(nspec*sizeof(flouble *));
  for(ii=0;ii<nspec;ii++) {
    cl_noise_bad_ub[ii]=my_calloc((lmax+1),sizeof(flouble));
    cl_proposal[ii]=my_calloc((lmax+1),sizeof(flouble));
    cl_dum[ii]=my_calloc((lmax+1),sizeof(flouble));
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

  //This could be OMPd TODO
  //Normalization matrices
  printf("Computing normalization matrices\n");
  flouble pixsize=4*M_PI/npix;
  gsl_matrix *matrix_M_1=NULL,*matrix_M_2;
  if(ntemp_1>0) {
    matrix_M_1=gsl_matrix_alloc(ntemp_1,ntemp_1);
    for(ii=0;ii<ntemp_1;ii++) {
      int jj;
      for(jj=ii;jj<ntemp_1;jj++) {
	int imap;
	double mat_element=0;
	for(imap=0;imap<nmaps_2;imap++) {
	  for(ip=0;ip<npix;ip++)
	    mat_element+=temp_1[ii][imap][ip]*temp_1[jj][imap][ip];
	}
	gsl_matrix_set(matrix_M_1,ii,jj,mat_element*pixsize);
	if(jj!=ii)
	  gsl_matrix_set(matrix_M_1,jj,ii,mat_element*pixsize);
      }
    }
    gsl_linalg_cholesky_decomp(matrix_M_1); //TODO: this won't necessarily be invertible
    gsl_linalg_cholesky_invert(matrix_M_1);
  }
  if(ntemp_2>0) { //TODO: this can be sped up by checking if M1=M2
    matrix_M_2=gsl_matrix_alloc(ntemp_2,ntemp_2);
    for(ii=0;ii<ntemp_2;ii++) {
      int jj;
      for(jj=ii;jj<ntemp_2;jj++) {
	int imap;
	double mat_element=0;
	for(imap=0;imap<nmaps_2;imap++) {
	  for(ip=0;ip<npix;ip++)
	    mat_element+=temp_2[ii][imap][ip]*temp_2[jj][imap][ip];
	}
	gsl_matrix_set(matrix_M_2,ii,jj,mat_element*pixsize);
	if(jj!=ii)
	  gsl_matrix_set(matrix_M_2,jj,ii,mat_element*pixsize);
      }
    }
    gsl_linalg_cholesky_decomp(matrix_M_2);
    gsl_linalg_cholesky_invert(matrix_M_2);
  }

  if((ntemp_1>0) || (ntemp_2>0)) {
    printf("Computing bias factors\n");
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

    //Allocate dummy maps and alms
    flouble **map_1_dum=my_malloc(nmaps_1*sizeof(flouble *));
    fcomplex **alm_1_dum=my_malloc(nmaps_1*sizeof(fcomplex *));
    for(ii=0;ii<nmaps_1;ii++) {
      map_1_dum[ii]=my_malloc(npix*sizeof(flouble));
      alm_1_dum[ii]=my_malloc(he_nalms(lmax)*sizeof(fcomplex));
    }
    flouble **map_2_dum=my_malloc(nmaps_2*sizeof(flouble *));
    fcomplex **alm_2_dum=my_malloc(nmaps_2*sizeof(fcomplex *));
    for(ii=0;ii<nmaps_2;ii++) {
      map_2_dum[ii]=my_malloc(npix*sizeof(flouble));
      alm_2_dum[ii]=my_malloc(he_nalms(lmax)*sizeof(fcomplex));
    }

    flouble **bias_F2=my_malloc(nspec*sizeof(flouble *));
    flouble **bias_F3=my_malloc(nspec*sizeof(flouble *));
    flouble **bias_F4=my_malloc(nspec*sizeof(flouble *));
    for(ii=0;ii<nspec;ii++) {
      bias_F2[ii]=my_calloc((lmax+1),sizeof(flouble));
      bias_F3[ii]=my_calloc((lmax+1),sizeof(flouble));
      bias_F4[ii]=my_calloc((lmax+1),sizeof(flouble));
    }

    if(ntemp_2>0) {
      printf("Computing F2\n");
      int iti;
      for(iti=0;iti<ntemp_2;iti++) {
	int itj;
	for(itj=0;itj<ntemp_2;itj++) {
	  int im1,im2;
	  double nij=gsl_matrix_get(matrix_M_2,iti,itj);
	  for(im2=0;im2<nmaps_2;im2++) {
	    for(ip=0;ip<npix;ip++)
	      map_2_dum[im2][ip]=temp_2[itj][im2][ip]*mask_2[ip]; //w*g^j
	    he_map2alm(nside_in,lmax,1,pol_2,map_2_dum,alm_2_dum); //SHT[w*g^j]
	  }
	  for(im1=0;im1<nmaps_1;im1++) {
	    for(im2=0;im2<nmaps_2;im2++)
	      he_alter_alm(lmax,-1.,alm_2_dum[im2],alm_1_dum[im1],cl_proposal[im1*nmaps_2+im2]); //C^ab*SHT[w*g^j]
	  }
	  he_alm2map(nside_in,lmax,1,pol_1,map_1_dum,alm_1_dum); //SHT^-1[C^ab*SHT[w*g^j]]
	  for(im1=0;im1<nmaps_1;im1++) {
	    for(ip=0;ip<npix;ip++)
	      map_1_dum[im1][ip]*=mask_1[ip]; //v*SHT^-1[C^ab*SHT[w*g^j]]
	  }
	  he_map2alm(nside_in,lmax,1,pol_1,map_1_dum,alm_1_dum); //SHT[v*SHT^-1[C^ab*SHT[w*g^j]]]
	  he_alm2cl(alm_1_dum,a_temp_2[iti],nmaps_1,nmaps_2,pol_1,pol_2,cl_dum,lmax);
	  //Sum_m(SHT[v*SHT^-1[C^ab*SHT[w*g^j]]]*g^i*)/(2l+1)
	  for(im1=0;im1<nspec;im1++) {
	    for(ip=0;ip<=lmax;ip++)
	      bias_F2[im1][ip]+=cl_dum[im1][ip]*nij;
	  }
	}
      }
    }

    if(ntemp_1>0) {
      printf("Computing F3\n");
      int iti;
      for(iti=0;iti<ntemp_1;iti++) {
	int itj;
	for(itj=0;itj<ntemp_1;itj++) {
	  int im1,im2;
	  double mij=gsl_matrix_get(matrix_M_1,iti,itj);
	  for(im1=0;im1<nmaps_1;im1++) {
	    for(ip=0;ip<npix;ip++)
	      map_1_dum[im1][ip]=temp_1[itj][im1][ip]*mask_1[ip]; //v*f^j
	    he_map2alm(nside_in,lmax,1,pol_1,map_1_dum,alm_1_dum); //SHT[v*f^j]
	  }
	  for(im2=0;im2<nmaps_2;im1++) {
	    for(im1=0;im1<nmaps_1;im2++)
	      he_alter_alm(lmax,-1.,alm_1_dum[im1],alm_2_dum[im2],cl_proposal[im1*nmaps_2+im2]); //C^abT*SHT[v*f^j]
	  }
	  he_alm2map(nside_in,lmax,1,pol_2,map_2_dum,alm_2_dum); //SHT^-1[C^abT*SHT[v*f^j]]
	  for(im2=0;im2<nmaps_2;im2++) {
	    for(ip=0;ip<npix;ip++)
	      map_2_dum[im1][ip]*=mask_2[ip]; //w*SHT^-1[C^abT*SHT[v*f^j]]
	  }
	  he_map2alm(nside_in,lmax,1,pol_2,map_2_dum,alm_2_dum); //SHT[w*SHT^-1[C^abT*SHT[v*f^j]]]
	  he_alm2cl(a_temp_1[iti],alm_2_dum,nmaps_1,nmaps_2,pol_1,pol_2,cl_dum,lmax);
	  //Sum_m(f^i*SHT[w*SHT^-1[C^abT*SHT[v*f^j]]]^*)/(2l+1)
	  for(im1=0;im1<nspec;im1++) {
	    for(ip=0;ip<=lmax;ip++)
	      bias_F3[im1][ip]+=cl_dum[im1][ip]*mij;
	  }
	}
      }
    }

    if((ntemp_1>0) && (ntemp_2>0)) {
      printf("Computing F4\n");
      int iti;
      for(iti=0;iti<ntemp_1;iti++) {
	int itj;
	for(itj=0;itj<ntemp_1;itj++) {
	  int itp;
	  double mij=gsl_matrix_get(matrix_M_1,iti,itj);
	  for(itp=0;itp<ntemp_2;itp++) {
	    int itq;
	    for(itq=0;itq<ntemp_2;itq++) {
	      double npq=gsl_matrix_get(matrix_M_2,itp,itq);
	      //HERE
	    }
	  }
	}
      }
    }
    
    for(ii=0;ii<nspec;ii++) {
      free(bias_F2[ii]);
      free(bias_F3[ii]);
      free(bias_F4[ii]);
    }
    free(bias_F2);
    free(bias_F3);
    free(bias_F4);
    for(ii=0;ii<nmaps_1;ii++) {
      free(map_1_dum[ii]);
      free(alm_1_dum[ii]);
    }
    free(map_1_dum);
    free(alm_1_dum);
    for(ii=0;ii<nmaps_2;ii++) {
      free(map_2_dum[ii]);
      free(alm_2_dum[ii]);
    }
    free(map_2_dum);
    free(alm_2_dum);
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
    flouble *cl_masks_bad_ub;
    cl_masks_bad_ub=my_malloc((lmax+1)*sizeof(flouble));
    printf("Computing mask pseudo-cl\n");
    he_anafast(&mask_1,&mask_2,1,1,0,0,&cl_masks_bad_ub,nside_in,lmax);

    printf("Computing coupling matrix \n");
    compute_coupling_matrix(cl_masks_bad_ub,nside_in,lmax,bins,
			    &coupling_matrix_b,&perm,coupling_fname,coupling_b_fname,
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
    free(cl_proposal[ii]);
    free(cl_dum[ii]);
    free(cl_maps_bad_ub[ii]);
    free(dl_maps_good_b[ii]);
  }
  free(cl_proposal);
  free(cl_dum);
  free(cl_noise_bad_ub);
  free(cl_maps_bad_ub);
  free(dl_maps_good_b);
  gsl_matrix_free(coupling_matrix_b);
  gsl_permutation_free(perm);
  if(ntemp_1>0) {
    for(ii=0;ii<ntemp_1;ii++) {
      int jj;
      for(jj=0;jj<nmaps_1;jj++) {
	free(temp_1[ii][jj]);
	free(a_temp_1[ii][jj]);
      }
      free(temp_1[ii]);
      free(a_temp_1[ii]);
    }
    free(temp_1);
    free(a_temp_1);
  }
  for(ii=0;ii<nmaps_1;ii++)
    free(maps_1[ii]);
  free(maps_1);
  free(mask_1);
  if(strcmp(fname_temp_1,fname_temp_2)) {
    if(ntemp_2>0) {
      for(ii=0;ii<ntemp_2;ii++) {
	int jj;
	for(jj=0;jj<nmaps_2;jj++) {
	  free(temp_2[ii][jj]);
	  free(a_temp_2[ii][jj]);
	}
	free(temp_2[ii]);
	free(a_temp_2[ii]);
      }
      free(temp_2);
      free(a_temp_2);
    }
  }
  if(strcmp(fname_maps_1,fname_maps_2)) {
    for(ii=0;ii<nmaps_2;ii++)
      free(maps_2[ii]);
    free(maps_2);
  }
  if(strcmp(fname_mask_1,fname_mask_2))
    free(mask_2);
  if(matrix_M_1!=NULL)
    gsl_matrix_free(matrix_M_1);
  if(matrix_M_2!=NULL)
    gsl_matrix_free(matrix_M_2);
}

int main(int argc,char **argv)
{
  int n_lbin=1,pol_1=0,pol_2=0;
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

  if(!strcmp(fname_map_2,"none")) {
    sprintf(fname_map_2,"%s",fname_map_1);
    pol_2=pol_1;
  }
  if(!strcmp(fname_mask_2,"none"))
    sprintf(fname_mask_2,"%s",fname_mask_1);

  run_master(fname_map_1,fname_map_2,
	     fname_mask_1,fname_mask_2,
	     fname_temp_1,fname_temp_2,
	     pol_1,pol_2,
	     fname_cl_noise,
	     fname_cl_proposal,
	     coupling_fname,
	     coupling_b_fname,
	     fname_out,fname_bins,n_lbin);

  return 0;
}
