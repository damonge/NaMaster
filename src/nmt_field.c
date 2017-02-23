#include "utils.h"

void nmt_field_free(nmt_field *fl)
{
  int imap;
  free(fl->beam);
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

nmt_field *nmt_field_alloc(long nside,flouble *mask,int pol,flouble **maps,
			   int ntemp,flouble ***temp,flouble *beam)
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

  fl->beam=my_malloc(3*fl->nside*sizeof(flouble));
  if(beam==NULL) {
    for(ii=0;ii<3*fl->nside;ii++)
      fl->beam[ii]=1.;
  }
  else
    memcpy(fl->beam,beam,3*fl->nside*sizeof(flouble));

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

nmt_field *nmt_field_read(char *fname_mask,char *fname_maps,char *fname_temp,char *fname_beam,int pol)
{
  long nside,nside_dum;
  nmt_field *fl;
  flouble *beam;
  flouble *mask;
  flouble **maps;
  flouble ***temp;
  int ii,ntemp,itemp,imap,nmaps=1;
  if(pol) nmaps=2;

  //Read mask and compute nside, lmax etc.
  mask=he_read_healpix_map(fname_mask,&nside,0);

  //Read beam
  if(!strcmp(fname_beam,"none"))
    beam=NULL;
  else {
    FILE *fi=my_fopen(fname_beam,"r");
    int nlines=my_linecount(fi); rewind(fi);
    if(nlines!=3*nside)
      report_error(1,"Beam file must have 3*nside rows and two columns\n");
    beam=my_malloc(3*nside*sizeof(flouble));
    for(ii=0;ii<3*nside;ii) {
      int l;
      double b;
      int stat=fscanf(fi,"%d %lf\n",&l,&b);
      if(stat!=2)
	report_error(1,"Error reading file %s, line %d\n",fname_beam,ii+1);
      if((l>3*nside-1) || (l<0))
	report_error(1,"Wrong multipole %d\n",l);
      beam[l]=b;
    }
  }

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

  fl=nmt_field_alloc(nside,mask,pol,maps,ntemp,temp,beam);

  if(beam!=NULL)
    free(beam);
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

static void apodize_mask_CX(long nside,flouble *mask_in,flouble *mask_out,flouble aposize,char *apotype)
{
  long ip,npix=he_nside2npix(nside);
  double aporad=aposize*M_PI/180;
  double x2_thr=1-cos(aporad);
  double inv_x2_thr=1./x2_thr;
  flouble *vec=my_malloc(3*npix*sizeof(flouble));
  flouble *cthv=my_malloc(npix*sizeof(flouble));
  flouble *phiv=my_malloc(npix*sizeof(flouble));
  int apotyp=0;
  if(!strcmp(apotype,"C1"))
    apotyp=0;
  else if(!strcmp(apotype,"C2"))
    apotyp=1;
  else
    report_error(1,"Unknown apodization type %s\n",apotype);

  if(mask_out!=mask_in)
    memcpy(mask_out,mask_in,npix*sizeof(flouble));

  //Get coords for each pixel
#pragma omp parallel default(none)		\
  shared(vec,npix,nside,cthv,phiv)
  {
    long ip;
#pragma omp for
    for(ip=0;ip<npix;ip++) {
      flouble *v=vec+3*ip;
      he_pix2vec_ring(nside,ip,v);
      cthv[ip]=v[2];
      if(v[2]==0)
	phiv[ip]=0;
      else {
	if(v[0]==0 && v[1]==0)
	  phiv[ip]=0;
	else
	  phiv[ip]=atan2(v[1],v[0]);
      }
    } //end omp for
  }//end omp parallel

#pragma omp parallel default(none)			\
  shared(vec,npix,x2_thr,inv_x2_thr,mask_in,mask_out)	\
  shared(nside,cthv,phiv,aporad,apotyp)
  {
    long ip;
    int lenlist0=(int)(2*npix*(1-cos(1.2*aporad)));
    int *listpix=my_malloc(lenlist0*sizeof(int));

#pragma omp for schedule(dynamic)
    for(ip=0;ip<npix;ip++) {
      if(mask_in[ip]>0) {
	int j;
	int lenlist_half=lenlist0/2;
	flouble *v0=vec+3*ip;
	flouble x2dist=1000;
	he_query_disc(nside,cthv[ip],phiv[ip],1.2*aporad,listpix,&lenlist_half,1);
	for(j=0;j<lenlist_half;j++) {
	  int ip2=listpix[j];
	  if(mask_in[ip2]<=0) {
	    flouble *v1=vec+3*ip2;
	    flouble x2=1-v0[0]*v1[0]-v0[1]*v1[1]-v0[2]*v1[2];
	    if(x2<x2dist) x2dist=x2;
	  }
	}
	if(x2dist<x2_thr) {
	  flouble f,xn;
	  if(x2dist<=0)
	    f=0;
	  else {
	    xn=sqrt(x2dist*inv_x2_thr);
	    if(apotyp==0)
	      f=xn-sin(xn*2*M_PI)/(2*M_PI);
	    else
	      f=0.5*(1-cos(xn*M_PI));
	  }
	  mask_out[ip]*=f;
	}
      }
    } //end omp for
    free(listpix);
  }//end omp parallel

  free(vec);
  free(cthv);
  free(phiv);
}

static void apodize_mask_smooth(long nside,flouble *mask_in,flouble *mask_out,flouble aposize)
{
  long ip,npix=he_nside2npix(nside);
  double aporad=aposize*M_PI/180;
  flouble *vec=my_malloc(3*npix*sizeof(flouble));
  flouble *cthv=my_malloc(npix*sizeof(flouble));
  flouble *phiv=my_malloc(npix*sizeof(flouble));
  flouble *mask_dum=my_malloc(npix*sizeof(flouble));
  fcomplex *alms_dum=my_malloc(he_nalms(3*nside-1)*sizeof(fcomplex));
  memcpy(mask_dum,mask_in,npix*sizeof(flouble));

  //Get coords for each pixel
#pragma omp parallel default(none)		\
  shared(vec,npix,nside,cthv,phiv)
  {
    long ip;
#pragma omp for
    for(ip=0;ip<npix;ip++) {
      flouble *v=vec+3*ip;
      he_pix2vec_ring(nside,ip,v);
      cthv[ip]=v[2];
      if(v[2]==0)
	phiv[ip]=0;
      else {
	if(v[0]==0 && v[1]==0)
	  phiv[ip]=0;
	else
	  phiv[ip]=atan2(v[1],v[0]);
      }
    } //end omp for
  }//end omp parallel

#pragma omp parallel default(none)	\
  shared(vec,npix,mask_in,mask_dum)	\
  shared(nside,cthv,phiv,aporad)
  {
    long ip;
    int lenlist0=(int)(2*npix*(1-cos(2.5*aporad)));
    int *listpix=my_malloc(lenlist0*sizeof(int));

#pragma omp for schedule(dynamic)
    for(ip=0;ip<npix;ip++) {
      if(mask_in[ip]<=0) {
	int j;
	int lenlist_half=lenlist0/2;
	he_query_disc(nside,cthv[ip],phiv[ip],2.5*aporad,listpix,&lenlist_half,1);
	for(j=0;j<lenlist_half;j++) {
	  int ip2=listpix[j];
	  mask_dum[ip2]=0;
	}
      }
    } //end omp for
    free(listpix);
  }//end omp parallel

  he_map2alm(nside,3*nside-1,1,0,&mask_dum,&alms_dum,3);
  he_alter_alm(3*nside-1,aporad*180*60*2.355/M_PI,alms_dum,alms_dum,NULL);
  he_alm2map(nside,3*nside-1,1,0,&mask_dum,&alms_dum);
  he_map_product(nside,mask_in,mask_dum,mask_out);

  free(vec);
  free(cthv);
  free(phiv);
  free(mask_dum);
  free(alms_dum);
}

void nmt_apodize_mask(long nside,flouble *mask_in,flouble *mask_out,flouble aposize,char *apotype)
{
  if((!strcmp(apotype,"C1")) || (!strcmp(apotype,"C2"))) {
    apodize_mask_CX(nside,mask_in,mask_out,aposize,apotype);
  }
  else if(!strcmp(apotype,"Smooth")) 
    apodize_mask_smooth(nside,mask_in,mask_out,aposize);
  else
    report_error(1,"Unknown apodization type %s. Allowed: \"Smooth\", \"C1\", \"C2\"");
}
