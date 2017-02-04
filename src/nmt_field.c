#include "utils.h"

void nmt_field_free(nmt_field *fl)
{
  int imap;
  for(imap=0;imap<fl->nmaps;imap++)
    free(fl->maps[imap]);
  free(fl->maps);
  free(fl->mask);
  if(fl->ntemp>0) {
    int itemp;
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      for(imap=0;imap<fl->nmaps;imap++)
	free(fl->temp[itemp][imap]);
      free(fl->temp[itemp]);
    }
    free(fl->temp);
    gsl_matrix_free(fl->matrix_M);
  }
  free(fl);
}

nmt_field *nmt_field_alloc(long nside,flouble *mask,int pol,flouble **maps,int ntemp,flouble ***temp)
{
  int ii,itemp,itemp2,imap;
  nmt_field *fl=my_malloc(sizeof(nmt_field));
  fl->nside=nside;
  fl->lmax=3*fl->nside-1;
  fl->npix=12*fl->nside*fl->nside;
  fl->pol=pol;
  if(pol) fl->nmaps=2;
  else fl->nmaps=1;
  fl->ntemp=ntemp;

  fl->mask=my_malloc(fl->npix*sizeof(flouble));
  memcpy(fl->mask,mask,fl->npix*sizeof(flouble));

  fl->maps=my_malloc(fl->nmaps*sizeof(flouble *));
  for(ii=0;ii<fl->nmaps;ii++) {
    fl->maps[ii]=my_malloc(fl->npix*sizeof(flouble));
    memcpy(fl->maps[ii],maps[ii],fl->npix*sizeof(flouble));
    he_map_product(fl->nside,fl->maps[ii],fl->mask,fl->maps[ii]);
  }

  if(fl->ntemp>0) {
    fl->temp=my_malloc(fl->ntemp*sizeof(flouble **));
    fl->a_temp=my_malloc(fl->ntemp*sizeof(fcomplex **));
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      fl->temp[itemp]=my_malloc(fl->nmaps*sizeof(flouble *));
      fl->a_temp[itemp]=my_malloc(fl->nmaps*sizeof(fcomplex *));
      for(imap=0;imap<fl->nmaps;imap++) {
	fl->temp[itemp][imap]=my_malloc(fl->npix*sizeof(flouble));
	fl->a_temp[itemp][imap]=my_malloc(he_nalms(fl->lmax)*sizeof(fcomplex));
	memcpy(fl->temp[itemp][imap],temp[itemp][imap],fl->npix*sizeof(flouble));
	he_map_product(fl->nside,fl->temp[itemp][imap],fl->mask,fl->temp[itemp][imap]); //Multiply by mask
      }
      he_map2alm(fl->nside,fl->lmax,1,fl->pol,fl->temp[itemp],fl->a_temp[itemp],HE_NITER_DEFAULT); //SHT
    }

    //Compute normalization matrix
    fl->matrix_M=gsl_matrix_alloc(fl->ntemp,fl->ntemp);
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      for(itemp2=itemp;itemp2<fl->ntemp;itemp2++) {
	flouble matrix_element=0;
	for(imap=0;imap<fl->nmaps;imap++)
	  matrix_element+=he_map_dot(fl->nside,fl->temp[itemp][imap],fl->temp[itemp2][imap]);
	gsl_matrix_set(fl->matrix_M,itemp,itemp2,matrix_element);
	if(itemp2!=itemp)
	  gsl_matrix_set(fl->matrix_M,itemp2,itemp,matrix_element);
      }
    }
    gsl_linalg_cholesky_decomp(fl->matrix_M); //TODO: this won't necessarily be invertible
    gsl_linalg_cholesky_invert(fl->matrix_M);

    //Deproject
    flouble *prods=my_calloc(fl->ntemp,sizeof(flouble));
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      for(imap=0;imap<fl->nmaps;imap++) 
	prods[itemp]+=he_map_dot(fl->nside,fl->temp[itemp][imap],fl->maps[imap]);
    }
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      flouble alpha=0;
      for(itemp2=0;itemp2<fl->ntemp;itemp2++) {
	double mij=gsl_matrix_get(fl->matrix_M,itemp,itemp2);
	alpha+=mij*prods[itemp2];
      }
#ifdef _DEBUG
      printf("alpha_%d = %lE\n",itemp,alpha);
#endif //_DEBUG
      for(imap=0;imap<fl->nmaps;imap++) {
	long ip;
	for(ip=0;ip<fl->npix;ip++)
	  fl->maps[imap][ip]-=alpha*fl->temp[itemp][imap][ip];
      }
    }
    free(prods);
  }

  return fl;
}

nmt_field *nmt_field_read(char *fname_mask,char *fname_maps,char *fname_temp,int pol)
{
  long nside,nside_dum;
  nmt_field *fl;
  flouble *mask;
  flouble **maps;
  flouble ***temp;
  int ii,ntemp,itemp,imap,nmaps=1;
  if(pol) nmaps=2;

  //Read mask and compute nside, lmax etc.
  mask=he_read_healpix_map(fname_mask,&nside,0);

  //Read data maps
  maps=my_malloc(nmaps*sizeof(flouble *));
  for(ii=0;ii<nmaps;ii++) {
    maps[ii]=he_read_healpix_map(fname_maps,&(nside_dum),ii);
    if(nside_dum!=nside)
      report_error(1,"Wrong nside %ld\n",nside_dum);
  }

  //Read templates and deproject
  if(strcmp(fname_temp,"none")) {
    int ncols,isnest;
    he_get_file_params(fname_temp,&nside_dum,&ncols,&isnest);
    if(nside_dum!=nside)
      report_error(1,"Wrong nside %ld\n",nside_dum);
    if((ncols==0) || (ncols%nmaps!=0))
      report_error(1,"Not enough templates in file %s\n",fname_temp);
    ntemp=ncols/nmaps;
    temp=my_malloc(ntemp*sizeof(flouble **));
    for(itemp=0;itemp<ntemp;itemp++) {
      temp[itemp]=my_malloc(nmaps*sizeof(flouble *));
      for(imap=0;imap<nmaps;imap++)
	temp[itemp][imap]=he_read_healpix_map(fname_temp,&(nside_dum),itemp*nmaps+imap);
    }
  }
  else {
    ntemp=0;
    temp=NULL;
  }

  fl=nmt_field_alloc(nside,mask,pol,maps,ntemp,temp);

  free(mask);
  for(imap=0;imap<nmaps;imap++)
    free(maps[imap]);
  free(maps);
  if(ntemp>0) {
    for(itemp=0;itemp<ntemp;itemp++) {
      for(imap=0;imap<nmaps;imap++)
	free(temp[itemp][imap]);
      free(temp[itemp]);
    }
    free(temp);
  }

  return fl;
}
