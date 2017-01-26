#include "common.h"

void field_free(Field *fl)
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

Field *field_alloc(char *fname_mask,char *fname_maps,char *fname_temp,int pol)
{
  int ii;
  long nside_dum;
  Field *fl=my_malloc(sizeof(Field));

  //Read mask and compute nside, lmax etc.
  fl->mask=he_read_healpix_map(fname_mask,&(fl->nside),0);
  fl->lmax=3*fl->nside-1;
  fl->npix=12*fl->nside*fl->nside;
  fl->pol=pol;
  if(pol) fl->nmaps=2;
  else fl->nmaps=1;

  //Read data maps
  fl->maps=my_malloc(fl->nmaps*sizeof(flouble *));
  for(ii=0;ii<fl->nmaps;ii++) {
    fl->maps[ii]=he_read_healpix_map(fname_maps,&(nside_dum),ii);
    if(nside_dum!=fl->nside)
      report_error(1,"Wrong nside %ld\n",nside_dum);
    he_map_product(fl->nside,fl->maps[ii],fl->mask,fl->maps[ii]);
  }

  //Read templates and deproject
  if(strcmp(fname_temp,"none")) {
    int ncols,isnest,itemp,itemp2,imap;
    he_get_file_params(fname_temp,&nside_dum,&ncols,&isnest);
    if(nside_dum!=fl->nside)
      report_error(1,"Wrong nside %ld\n",nside_dum);
    if((ncols==0) || (ncols%fl->nmaps!=0))
      report_error(1,"Not enough templates in file %s\n",fname_temp);
    fl->ntemp=ncols/fl->nmaps;
    fl->temp=my_malloc(fl->ntemp*sizeof(flouble **));
    fl->a_temp=my_malloc(fl->ntemp*sizeof(fcomplex **));
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      fl->temp[itemp]=my_malloc(fl->nmaps*sizeof(flouble *));
      fl->a_temp[itemp]=my_malloc(fl->nmaps*sizeof(fcomplex *));
      for(imap=0;imap<fl->nmaps;imap++) {
	fl->a_temp[itemp][imap]=my_malloc(he_nalms(fl->lmax)*sizeof(fcomplex));
	fl->temp[itemp][imap]=he_read_healpix_map(fname_temp,&(nside_dum),itemp*fl->nmaps+imap);
	he_map_product(fl->nside,fl->temp[itemp][imap],fl->mask,fl->temp[itemp][imap]); //Multiply by mask
	he_map2alm(fl->nside,fl->lmax,1,fl->pol,fl->temp[itemp],fl->a_temp[itemp]); //SHT
      }
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
      printf("alpha_%d = %lE\n",itemp2,alpha);
#endif //_DEBUG
      for(imap=0;imap<fl->nmaps;imap++) {
	long ip;
	for(ip=0;ip<fl->npix;ip++)
	  fl->maps[imap][ip]-=alpha*fl->temp[itemp][imap][ip];
      }
    }
    free(prods);
  }
  else
    fl->ntemp=0;
  return fl;
}
