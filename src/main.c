#include "common.h"
#include <unistd.h>

void run_master(char *fname_maps,
		char *fname_cl_noise,
		char *fname_mask,
		char *coupling_prefix,
		char *prefix_out,
		int nmaps,int n_lbin)
{
  long ip,nside_in,npix;
  int ii,nspec,lmax,nbins;
  flouble fsky;
  flouble **cl_noise_bad_ub,**cl_maps_bad_ub,**cl_masks_bad_ub,**dl_maps_good_b;
  flouble **masks,**maps;
  gsl_matrix **coupling_matrix_b;
  gsl_permutation **perm;
  char fname[256];

  //Mask
  printf("Reading masks\n");
  masks=(flouble **)dam_malloc(nmaps*sizeof(flouble *));
  for(ii=0;ii<nmaps;ii++) {
    long nside_dum;
    masks[ii]=he_read_healpix_map(fname_mask,&(nside_dum),ii);
    if(ii==0) {
      nside_in=nside_dum;
      lmax=3*nside_in-1;
      nbins=(lmax-1)/n_lbin;
      npix=12*nside_in*nside_in;
    }
    else {
      if(nside_dum!=nside_in)
	dam_report_error(1,"Wrong nside\n");
    }
    fsky=0;
    for(ip=0;ip<npix;ip++)
      fsky+=masks[ii][ip];
    fsky/=npix;
    printf("f_sky=%lE\n",fsky);
  }

  //Maps
  printf("Reading maps\n");
  maps=(flouble **)dam_malloc(nmaps*sizeof(flouble *));
  for(ii=0;ii<nmaps;ii++) {
    long nside_dum;
    maps[ii]=he_read_healpix_map(fname_maps,&(nside_dum),ii);
    if(nside_dum!=nside_in)
      dam_report_error(1,"Wrong nside\n");
    for(ip=0;ip<npix;ip++)
      maps[ii][ip]*=masks[ii][ip];
  }

  //Allocate cl
  nspec=nmaps*(nmaps+1)/2;
  cl_noise_bad_ub=(flouble **)dam_malloc(nspec*sizeof(flouble *));
  cl_maps_bad_ub=(flouble **)dam_malloc(nspec*sizeof(flouble *));
  cl_masks_bad_ub=(flouble **)dam_malloc(nspec*sizeof(flouble *));
  dl_maps_good_b=(flouble **)dam_malloc(nspec*sizeof(flouble *));
  coupling_matrix_b=(gsl_matrix **)dam_malloc(nspec*sizeof(gsl_matrix *));
  perm=(gsl_permutation **)dam_malloc(nspec*sizeof(gsl_permutation *));
  for(ii=0;ii<nspec;ii++) {
    cl_noise_bad_ub[ii]=(flouble *)dam_malloc((lmax+1)*sizeof(flouble));
    cl_maps_bad_ub[ii]=(flouble *)dam_malloc((lmax+1)*sizeof(flouble));
    cl_masks_bad_ub[ii]=(flouble *)dam_malloc((lmax+1)*sizeof(flouble));
  }

  //Noise
  printf("Reading noise pseudo-cl\n");
  FILE *fi=dam_fopen(fname_cl_noise,"r");
  int nlin=dam_linecount(fi); rewind(fi);
  if(nlin!=lmax+1)
    dam_report_error(1,"Wrong number of multipoles for noise p.spec.\n");
  for(ii=0;ii<lmax+1;ii++) {
    int status,jj;
    flouble l;
    status=fscanf(fi,"%lf",&l);
    if(status!=1)
      dam_report_error(1,"Error reading file %s\n",fname_cl_noise);
    for(jj=0;jj<nspec;jj++) {
      status=fscanf(fi,"%lf",&(cl_noise_bad_ub[jj][ii]));
      if(status!=1)
	dam_report_error(1,"Error reading file %s\n",fname_cl_noise);
    }
  }
  fclose(fi);

  //Pseudo-cl
  printf("Computing data pseudo-cl\n");
  he_anafast(maps,cl_maps_bad_ub,nside_in,lmax,nmaps,0);

  printf("Computing mask pseudo-cl\n");
  he_anafast(masks,cl_masks_bad_ub,nside_in,lmax,nmaps,0);

  //Coupling matrix
  int index=0;
  for(ii=0;ii<nmaps;ii++) {
    int jj;
    for(jj=ii;jj<nmaps;jj++) {
      char fname[256];
      sprintf(fname,"%s_%d%d.dat",coupling_prefix,ii,jj);
      if(access(fname,F_OK)!=-1) { //If file exists just read matrix
	printf("Reading coupling matrix\n");
	read_coupling_matrix(fname,nbins,
			     &(coupling_matrix_b[index]),&(perm[index]));
      }
      else { //Else, compute it
	printf("Computing coupling matrix %d %d %d\n",ii,jj,index);
	compute_coupling_matrix(cl_masks_bad_ub[index],n_lbin,nside_in,lmax,
				nbins,&(coupling_matrix_b[index]),&(perm[index]),
				fname);
      }
      index++;
    }
  }

  //Decouple cls
  printf("Decoupling cls\n");
  for(ii=0;ii<nspec;ii++) {
    dl_maps_good_b[ii]=decouple_cl_l(cl_maps_bad_ub[ii],cl_noise_bad_ub[ii],
				     nbins,n_lbin,coupling_matrix_b[ii],perm[ii]);
  }

  //Write output
  printf("Writing output\n");
  sprintf(fname,"%s_dl_decoupled.dat",prefix_out);
  fi=dam_fopen(fname,"w");
  for(ii=0;ii<nbins;ii++) {
    int jj;
    double l_here=2+ii*n_lbin+0.5*(n_lbin-1.);
    fprintf(fi,"%.2lf ",l_here);
    for(jj=0;jj<nspec;jj++)
      fprintf(fi,"%lE ",dl_maps_good_b[jj][ii]);
    fprintf(fi,"\n");
  }
  fclose(fi);

  for(ii=0;ii<nspec;ii++) {
    free(cl_noise_bad_ub[ii]);
    free(cl_maps_bad_ub[ii]);
    free(cl_masks_bad_ub[ii]);
    free(dl_maps_good_b[ii]);
    gsl_matrix_free(coupling_matrix_b[ii]);
    gsl_permutation_free(perm[ii]);
  }
  free(cl_noise_bad_ub);
  free(cl_maps_bad_ub);
  free(dl_maps_good_b);
  free(coupling_matrix_b);
  free(perm);
  for(ii=0;ii<nmaps;ii++) {
    free(maps[ii]);
    free(masks[ii]);
  }
  free(maps);
  free(masks);
}  

int main(int argc,char **argv)
{
  int nmaps,n_lbin;
  char fname_maps[256],fname_cl_noise[256],fname_mask[256],coupling_prefix[256],prefix_out[256];
  if(argc!=8) {
    fprintf(stderr,"Usage: namaster fname_maps fname_cl_noise fname_mask "
	    "coupling_prefix prefix_out nmaps n_lbin\n");
    exit(0);
  }
  sprintf(fname_maps,"%s",argv[1]);
  sprintf(fname_cl_noise,"%s",argv[2]);
  sprintf(fname_mask,"%s",argv[3]);
  sprintf(coupling_prefix,"%s",argv[4]);
  sprintf(prefix_out,"%s",argv[5]);
  nmaps=atoi(argv[6]);
  n_lbin=atoi(argv[7]);

  run_master(fname_maps,fname_cl_noise,fname_mask,coupling_prefix,prefix_out,nmaps,n_lbin);

  return 0;
}
