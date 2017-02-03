#include "common.h"

static flouble weigh_l(int l)
{
#ifdef _WEIGH_L2
  return (flouble)(l*(l+1))/(2*M_PI);
#else //_WEIGH_L2
  return 1.;
#endif //_WEIGH_L2
}

//Returns all non-zero wigner-3j symbols
// il2 (in) : l2
// il3 (in) : l3
// im2 (in) : m2
// im3 (in) : m3
// l1min_out (out) : min value for l1
// l1max_out (out) : max value for l1
// thrcof (out) : array with the values of the wigner-3j
// size (in) : size allocated for thrcof
static int drc3jj(int il2,int il3,int im2, int im3,int *l1min_out,
		  int *l1max_out,double *thrcof,int size)
{
  int sign1,sign2,nfin,im1,l1max,l1min,ii,lstep;
  int converging,nstep2,nfinp1,index,nlim;
  double newfac,c1,c2,sum1,sum2,a1,a2,a1s,a2s,dv,denom,c1old,oldfac,l1,l2,l3,m1,m2,m3;
  double x,x1,x2,x3,y,y1,y2,y3,sumfor,sumbac,sumuni,cnorm,thresh,ratio;
  double huge=sqrt(1.79E308/20.0);
  double srhuge=sqrt(huge);
  double tiny=1./huge;
  double srtiny=1./srhuge;

  im1=-im2-im3;
  l2=(double)il2; l3=(double)il3;
  m1=(double)im1; m2=(double)im2; m3=(double)im3;
  
  if((abs(il2+im2-il3+im3))%2==0)
    sign2=1;
  else
    sign2=-1;
  
  if((il2-abs(im2)<0)||(il3-abs(im3)<0))
    report_error(1,"Wrong arguments: %d %d %d %d\n",il2,il3,im2,im3);
  
  //l1 bounds
  l1max=il2+il3;
  l1min=COM_MAX((abs(il2-il3)),(abs(im1)));
  *l1max_out=l1max;
  *l1min_out=l1min;
  
  if(l1max-l1min<0) //Check for meaningful values
    report_error(0,"WTF?\n");
  
  if(l1max==l1min) { //If it's only one value:
    thrcof[0]=sign2/sqrt(l1min+l2+l3+1);
    return 0;
  }
  else {
    nfin=l1max-l1min+1;
    if(nfin>size) //Check there's enough space
      report_error(1,"Output array is too small %d\n",nfin);
    else {
      l1=l1min;
      newfac=0.;
      c1=0.;
      sum1=(l1+l1+1)*tiny;
      thrcof[0]=srtiny;
      
      lstep=0;
      converging=1;
      while((lstep<nfin-1)&&(converging)) { //Forward series
	lstep++;
	l1++; //order
	
	oldfac=newfac;
	a1=(l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+1);
	a2=(l1+m1)*(l1-m1);
	newfac=sqrt(a1*a2);
	
	if(l1>1) {
	  dv=-l2*(l2+1)*m1+l3*(l3+1)*m1+l1*(l1-1)*(m3-m2);
	  denom=(l1-1)*newfac;
	  if(lstep>1)
	    c1old=fabs(c1);
	  c1=-(l1+l1-1)*dv/denom;
	}
	else {
	  c1=-(l1+l1-1)*l1*(m3-m2)/newfac;
	}
	
	if(lstep<=1) {
	  x=srtiny*c1;
	  thrcof[1]=x;
	  sum1+=tiny*(l1+l1+1)*c1*c1;
	}
	else {
	  c2=-l1*oldfac/denom;
	  x=c1*thrcof[lstep-1]+c2*thrcof[lstep-2];
	  thrcof[lstep]=x;
	  sumfor=sum1;
	  sum1+=(l1+l1+1)*x*x;
	  if(lstep<nfin-1) {
	    if(fabs(x)>=srhuge) {
	      for(ii=0;ii<=lstep;ii++) {
		if(fabs(thrcof[ii])<srtiny)
		  thrcof[ii]=0;
		thrcof[ii]/=srhuge;
	      }
	      sum1/=huge;
	      sumfor/=huge;
	      x/=srhuge;
	    }
	    
	    if(c1old<=fabs(c1))
	      converging=0;
	  }
	}
      }
      
      if(nfin>2) {
	x1=x;
	x2=thrcof[lstep-1];
	x3=thrcof[lstep-2];
	nstep2=nfin-lstep-1+3;
	
	nfinp1=nfin+1;
	l1=l1max;
	thrcof[nfin-1]=srtiny;
	sum2=tiny*(l1+l1+1);
	
	l1+=2;
	lstep=0;
	while(lstep<nstep2-1) { //Backward series
	  lstep++;
	  l1--;
	  
	  oldfac=newfac;
	  a1s=(l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)*(-l1+l2+l3+2);
	  a2s=(l1+m1-1)*(l1-m1-1);
	  newfac=sqrt(a1s*a2s);
	  
	  dv=-l2*(l2+1)*m1+l3*(l3+1)*m1+l1*(l1-1)*(m3-m2);
	  denom=l1*newfac;
	  c1=-(l1+l1-1)*dv/denom;
	  if(lstep<=1) {
	    y=srtiny*c1;
	    thrcof[nfin-2]=y;
	    sumbac=sum2;
	    sum2+=tiny*(l1+l1-3)*c1*c1;
	  }
	  else {
	    c2=-(l1-1)*oldfac/denom;
	    y=c1*thrcof[nfin-lstep]+c2*thrcof[nfinp1-lstep]; //is the index ok??
	    if(lstep!=nstep2-1) {
	      thrcof[nfin-lstep-1]=y; //is the index ok??
	      sumbac=sum2;
	      sum2+=(l1+l1-3)*y*y;
	      if(fabs(y)>=srhuge) {
		for(ii=0;ii<=lstep;ii++) {
		  index=nfin-ii-1; //is the index ok??
		  if(fabs(thrcof[index])<srtiny)
		    thrcof[index]=0;
		  thrcof[index]=thrcof[index]/srhuge;
		}
		sum2/=huge;
		sumbac/=huge;
	      }
	    }
	  }
	}
	
	y3=y;
	y2=thrcof[nfin-lstep]; //is the index ok??
	y1=thrcof[nfinp1-lstep]; //is the index ok??
	
	ratio=(x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3);
	nlim=nfin-nstep2+1;
	
	if(fabs(ratio)<1) {
	  nlim++;
	  ratio=1./ratio;
	  for(ii=nlim-1;ii<nfin;ii++) //is the index ok??
	    thrcof[ii]*=ratio;
	  sumuni=ratio*ratio*sumbac+sumfor;
	}
	else {
	  for(ii=0;ii<nlim;ii++)
	    thrcof[ii]*=ratio;
	  sumuni=ratio*ratio*sumfor+sumbac;
	}
      }
      else
	sumuni=sum1;
      
      cnorm=1./sqrt(sumuni);
      if(thrcof[nfin-1]<0) sign1=-1;
      else sign1=1;
      
      if(sign1*sign2<=0)
	cnorm=-cnorm;
      if(fabs(cnorm)>=1) {
	for(ii=0;ii<nfin;ii++)
	  thrcof[ii]*=cnorm;
	return 0;
      }
      else {
	thresh=tiny/fabs(cnorm);
	for(ii=0;ii<nfin;ii++) {
	  if(fabs(thrcof[ii])<thresh)
	    thrcof[ii]=0;
	  thrcof[ii]*=cnorm;
	}
	return 0;
      }
    } //Size is good
  } //Doing for many l1s
  
  return 2;
}

static MasterWorkspace *master_workspace_new(int lmax,int ncls,BinSchm *bin)
{
  int ii;
  MasterWorkspace *w=my_malloc(sizeof(MasterWorkspace));
  w->lmax=lmax;
  w->ncls=ncls;

  w->pcl_masks=my_malloc((w->lmax+1)*sizeof(flouble));

  w->coupling_matrix_unbinned=my_malloc(w->ncls*(w->lmax+1)*sizeof(flouble *));
  for(ii=0;ii<w->ncls*(w->lmax+1);ii++)
    w->coupling_matrix_unbinned[ii]=my_calloc(w->ncls*(w->lmax+1),sizeof(flouble));

  w->bin=my_malloc(sizeof(BinSchm));
  w->bin->n_bands=bin->n_bands;
  w->bin->nell_list=my_malloc(w->bin->n_bands*sizeof(int));
  memcpy(w->bin->nell_list,bin->nell_list,w->bin->n_bands*sizeof(int));

  w->bin->ell_list=my_malloc(w->bin->n_bands*sizeof(int *));
  w->bin->w_list=my_malloc(w->bin->n_bands*sizeof(flouble *));
  for(ii=0;ii<w->bin->n_bands;ii++) {
    w->bin->ell_list[ii]=my_malloc(w->bin->nell_list[ii]*sizeof(int));
    w->bin->w_list[ii]=my_malloc(w->bin->nell_list[ii]*sizeof(flouble));
    memcpy(w->bin->ell_list[ii],bin->ell_list[ii],w->bin->nell_list[ii]*sizeof(int));
    memcpy(w->bin->w_list[ii],bin->w_list[ii],w->bin->nell_list[ii]*sizeof(flouble));
  }

  w->coupling_matrix_binned=gsl_matrix_alloc(w->ncls*w->bin->n_bands,w->ncls*w->bin->n_bands);
  w->coupling_matrix_perm=gsl_permutation_alloc(w->ncls*w->bin->n_bands);

  return w;
}

void master_workspace_free(MasterWorkspace *w)
{
  int ii;
  gsl_permutation_free(w->coupling_matrix_perm);
  gsl_matrix_free(w->coupling_matrix_binned);
  bins_free(w->bin);
  for(ii=0;ii<w->ncls*(w->lmax+1);ii++)
    free(w->coupling_matrix_unbinned[ii]);
  free(w->coupling_matrix_unbinned);
  free(w->pcl_masks);
  free(w);
}

MasterWorkspace *read_master_workspace(char *fname)
{
  int ii;
  MasterWorkspace *w=my_malloc(sizeof(MasterWorkspace));
  FILE *fi=my_fopen(fname,"rb");

  my_fread(&(w->lmax),sizeof(int),1,fi);
  my_fread(&(w->ncls),sizeof(int),1,fi);

  w->pcl_masks=my_malloc((w->lmax+1)*sizeof(flouble));
  my_fread(w->pcl_masks,sizeof(flouble),w->lmax+1,fi);

  w->coupling_matrix_unbinned=my_malloc(w->ncls*(w->lmax+1)*sizeof(flouble *));
  for(ii=0;ii<w->ncls*(w->lmax+1);ii++) {
    w->coupling_matrix_unbinned[ii]=my_malloc(w->ncls*(w->lmax+1)*sizeof(flouble));
    my_fread(w->coupling_matrix_unbinned[ii],sizeof(flouble),w->ncls*(w->lmax+1),fi);
  }

  w->bin=my_malloc(sizeof(BinSchm));
  my_fread(&(w->bin->n_bands),sizeof(int),1,fi);
  w->bin->nell_list=my_malloc(w->bin->n_bands*sizeof(int));
  w->bin->ell_list=my_malloc(w->bin->n_bands*sizeof(int *));
  w->bin->w_list=my_malloc(w->bin->n_bands*sizeof(flouble *));
  my_fread(w->bin->nell_list,sizeof(int),w->bin->n_bands,fi);
  for(ii=0;ii<w->bin->n_bands;ii++) {
    w->bin->ell_list[ii]=my_malloc(w->bin->nell_list[ii]*sizeof(int));
    w->bin->w_list[ii]=my_malloc(w->bin->nell_list[ii]*sizeof(flouble));
    my_fread(w->bin->ell_list[ii],sizeof(int),w->bin->nell_list[ii],fi);
    my_fread(w->bin->w_list[ii],sizeof(flouble),w->bin->nell_list[ii],fi);
  }

  w->coupling_matrix_binned=gsl_matrix_alloc(w->ncls*w->bin->n_bands,w->ncls*w->bin->n_bands);
  w->coupling_matrix_perm=gsl_permutation_alloc(w->ncls*w->bin->n_bands);
  gsl_matrix_fread(fi,w->coupling_matrix_binned);
  gsl_permutation_fread(fi,w->coupling_matrix_perm);
  
  fclose(fi);

  return w;
}

void write_master_workspace(MasterWorkspace *w,char *fname)
{
  int ii;
  FILE *fo=my_fopen(fname,"wb");

  my_fwrite(&(w->lmax),sizeof(int),1,fo);
  my_fwrite(&(w->ncls),sizeof(int),1,fo);
  my_fwrite(w->pcl_masks,sizeof(flouble),w->lmax+1,fo);
  for(ii=0;ii<w->ncls*(w->lmax+1);ii++)
    my_fwrite(w->coupling_matrix_unbinned[ii],sizeof(flouble),w->ncls*(w->lmax+1),fo);

  my_fwrite(&(w->bin->n_bands),sizeof(int),1,fo);
  my_fwrite(w->bin->nell_list,sizeof(int),w->bin->n_bands,fo);
  for(ii=0;ii<w->bin->n_bands;ii++) {
    my_fwrite(w->bin->ell_list[ii],sizeof(int),w->bin->nell_list[ii],fo);
    my_fwrite(w->bin->w_list[ii],sizeof(flouble),w->bin->nell_list[ii],fo);
  }

  gsl_matrix_fwrite(fo,w->coupling_matrix_binned);
  gsl_permutation_fwrite(fo,w->coupling_matrix_perm);
  
  fclose(fo);
}

static void bin_coupling_matrix(MasterWorkspace *w)
{
  int icl_a,icl_b,ib2,ib3,l2,l3,i2,i3,sig;
  
  for(icl_a=0;icl_a<w->ncls;icl_a++) {
    for(icl_b=0;icl_b<w->ncls;icl_b++) {
      for(ib2=0;ib2<w->bin->n_bands;ib2++) {
	for(ib3=0;ib3<w->bin->n_bands;ib3++) {
	  double coupling_b=0;
	  for(i2=0;i2<w->bin->nell_list[ib2];i2++) {
	    l2=w->bin->ell_list[ib2][i2];
	    for(i3=0;i3<w->bin->nell_list[ib3];i3++) {
	      l3=w->bin->ell_list[ib3][i3];
	      coupling_b+=w->coupling_matrix_unbinned[w->ncls*l2+icl_a][w->ncls*l3+icl_b]*
		w->bin->w_list[ib2][i2]*weigh_l(l2)/weigh_l(l3);
	    }
	  }
	  gsl_matrix_set(w->coupling_matrix_binned,w->ncls*ib2+icl_a,w->ncls*ib3+icl_b,coupling_b);
	}
      }
    }
  }

  gsl_linalg_LU_decomp(w->coupling_matrix_binned,w->coupling_matrix_perm,&sig);
}

//Computes binned coupling matrix
// fl1,fl2 (in) : fields we're correlating
// coupling_matrix_out (out) : unbinned coupling matrix
MasterWorkspace *compute_coupling_matrix(Field *fl1,Field *fl2,BinSchm *bin)
{
  int l2;
  MasterWorkspace *w;
  int n_cl=fl1->nmaps*fl2->nmaps;
  if(fl1->nside!=fl2->nside)
    report_error(1,"Can't correlate fields with different resolutions\n");
  w=master_workspace_new(fl1->lmax,n_cl,bin);
  he_anafast(&(fl1->mask),&(fl2->mask),0,0,&(w->pcl_masks),fl1->nside,fl1->lmax);
  for(l2=0;l2<=fl1->lmax;l2++)
    w->pcl_masks[l2]*=(2*l2+1.);

#pragma omp parallel default(none)				\
  shared(w)
  {
    int ll2,ll3;
    double *wigner_00=NULL,*wigner_22=NULL;
    int lstart=0;

    if((w->ncls==1) || (w->ncls==2))
      wigner_00=my_malloc(2*(w->lmax+1)*sizeof(double));
    if((w->ncls==2) || (w->ncls==4))
      wigner_22=my_malloc(2*(w->lmax+1)*sizeof(double));
   
    if(w->ncls>1)
      lstart=2;

#pragma omp for schedule(dynamic)
    for(ll2=lstart;ll2<=w->lmax;ll2++) {
      for(ll3=lstart;ll3<=w->lmax;ll3++) {
	int jj,l1,lmin_here,lmax_here;
	int lmin_here_00=0,lmax_here_00=2*(w->lmax+1)+1;
	int lmin_here_22=0,lmax_here_22=2*(w->lmax+1)+1;

	if((w->ncls==1) || (w->ncls==2))
	  drc3jj(ll2,ll3,0,0,&lmin_here_00,&lmax_here_00,wigner_00,2*(w->lmax+1));
	if((w->ncls==2) || (w->ncls==4))
	  drc3jj(ll2,ll3,2,-2,&lmin_here_22,&lmax_here_22,wigner_22,2*(w->lmax+1));

	lmin_here=COM_MAX(lmin_here_00,lmin_here_22);
	lmax_here=COM_MIN(lmax_here_00,lmax_here_22);
	  
	for(l1=lmin_here;l1<=lmax_here;l1++) {
	  if(l1<=w->lmax) {
	    double wfac;
	    int j00=l1-lmin_here_00;
	    int j22=l1-lmin_here_22;
	    if(w->ncls==1) {
	      wfac=w->pcl_masks[l1]*wigner_00[j00]*wigner_00[j00];
	      w->coupling_matrix_unbinned[1*ll2+0][1*ll3+0]+=wfac; //TT,TT
	    }
	    if(w->ncls==2) {
	      wfac=w->pcl_masks[l1]*wigner_00[j00]*wigner_22[j22];
	      w->coupling_matrix_unbinned[2*ll2+0][2*ll3+0]+=wfac; //TE,TE
	      w->coupling_matrix_unbinned[2*ll2+1][2*ll3+1]+=wfac; //TB,TB
	    }
	    if(w->ncls==4) {
	      int suml=l1+ll2+ll3;
	      wfac=w->pcl_masks[l1]*wigner_22[j22]*wigner_22[j22];
	      if(suml & 1) { //Odd sum
		w->coupling_matrix_unbinned[4*ll2+0][4*ll3+3]+=wfac; //EE,BB
		w->coupling_matrix_unbinned[4*ll2+1][4*ll3+2]-=wfac; //EB,BE
		w->coupling_matrix_unbinned[4*ll2+2][4*ll3+1]-=wfac; //BE,EB
		w->coupling_matrix_unbinned[4*ll2+3][4*ll3+0]+=wfac; //BB,EE
	      }
	      else {
		w->coupling_matrix_unbinned[4*ll2+0][4*ll3+0]+=wfac; //EE,EE
		w->coupling_matrix_unbinned[4*ll2+1][4*ll3+1]+=wfac; //EB,EB
		w->coupling_matrix_unbinned[4*ll2+2][4*ll3+2]+=wfac; //BE,BE
		w->coupling_matrix_unbinned[4*ll2+3][4*ll3+3]+=wfac; //BB,BB
	      }
	    }
	  }
	}
	for(jj=0;jj<w->ncls;jj++) {
	  int kk;
	  for(kk=0;kk<w->ncls;kk++)
	    w->coupling_matrix_unbinned[w->ncls*ll2+jj][w->ncls*ll3+kk]*=(2*ll3+1.)/(4*M_PI);
	}
      }
    } //end omp for
    if((w->ncls==1) || (w->ncls==2))
      free(wigner_00);
    if((w->ncls==2) || (w->ncls==4))
      free(wigner_22);
  } //end omp parallel

  bin_coupling_matrix(w);

  return w;
}

void compute_deprojection_bias(Field *fl1,Field *fl2,flouble **cl_proposal,flouble **cl_bias)
{
  int ii;
  flouble **cl_dum;
  long ip;
  int nspec=fl1->nmaps*fl2->nmaps;
  int lmax=fl1->lmax;

  cl_dum=my_malloc(nspec*sizeof(flouble *));
  for(ii=0;ii<nspec;ii++) {
    cl_dum[ii]=my_calloc((lmax+1),sizeof(flouble));
    for(ip=0;ip<=lmax;ip++)
      cl_bias[ii][ip]=0;
  }

  //TODO: some terms (e.g. C^ab*SHT[w*g^j]) could be precomputed
  //TODO: if fl1=fl2 F2=F3
  //Allocate dummy maps and alms
  flouble **map_1_dum=my_malloc(fl1->nmaps*sizeof(flouble *));
  fcomplex **alm_1_dum=my_malloc(fl1->nmaps*sizeof(fcomplex *));
  for(ii=0;ii<fl1->nmaps;ii++) {
    map_1_dum[ii]=my_malloc(fl1->npix*sizeof(flouble));
    alm_1_dum[ii]=my_malloc(he_nalms(fl1->lmax)*sizeof(fcomplex));
  }
  flouble **map_2_dum=my_malloc(fl2->nmaps*sizeof(flouble *));
  fcomplex **alm_2_dum=my_malloc(fl2->nmaps*sizeof(fcomplex *));
  for(ii=0;ii<fl2->nmaps;ii++) {
    map_2_dum[ii]=my_malloc(fl1->npix*sizeof(flouble));
    alm_2_dum[ii]=my_malloc(he_nalms(fl1->lmax)*sizeof(fcomplex));
  }

  if(fl2->ntemp>0) {
    printf("Computing F2\n");
    int iti;
    for(iti=0;iti<fl2->ntemp;iti++) {
      int itj;
      for(itj=0;itj<fl2->ntemp;itj++) {
	int im1,im2;
	double nij=gsl_matrix_get(fl2->matrix_M,iti,itj);
	//w*g^j
	for(im2=0;im2<fl2->nmaps;im2++)
	  he_map_product(fl2->nside,fl2->temp[itj][im2],fl2->mask,map_2_dum[im2]);
	//SHT[w*g^j]
	he_map2alm(fl2->nside,fl2->lmax,1,fl2->pol,map_2_dum,alm_2_dum,HE_NITER_DEFAULT);
	//C^ab*SHT[w*g^j]
	for(im1=0;im1<fl1->nmaps;im1++) {
	  for(im2=0;im2<fl2->nmaps;im2++)
	    he_alter_alm(lmax,-1.,alm_2_dum[im2],alm_1_dum[im1],cl_proposal[im1*fl2->nmaps+im2]);
	}
	//SHT^-1[C^ab*SHT[w*g^j]]
	he_alm2map(fl1->nside,fl1->lmax,1,fl1->pol,map_1_dum,alm_1_dum);
	//v*SHT^-1[C^ab*SHT[w*g^j]]
	for(im1=0;im1<fl1->nmaps;im1++)
	  he_map_product(fl1->nside,map_1_dum[im1],fl1->mask,map_1_dum[im1]);
	//SHT[v*SHT^-1[C^ab*SHT[w*g^j]]]
	he_map2alm(fl1->nside,fl1->lmax,1,fl1->pol,map_1_dum,alm_1_dum,HE_NITER_DEFAULT);
	//Sum_m(SHT[v*SHT^-1[C^ab*SHT[w*g^j]]]*g^i*)/(2l+1)
	he_alm2cl(alm_1_dum,fl2->a_temp[iti],fl1->pol,fl1->pol,cl_dum,lmax);
	for(im1=0;im1<nspec;im1++) {
	  for(ip=0;ip<=lmax;ip++)
	    cl_bias[im1][ip]-=cl_dum[im1][ip]*nij;
	}
      }
    }
  }

  if(fl1->ntemp>0) {
    printf("Computing F3\n");
    int iti;
    for(iti=0;iti<fl1->ntemp;iti++) {
      int itj;
      for(itj=0;itj<fl1->ntemp;itj++) {
	int im1,im2;
	double mij=gsl_matrix_get(fl1->matrix_M,iti,itj);
	//v*f^j
	for(im1=0;im1<fl1->nmaps;im1++)
	  he_map_product(fl1->nside,fl1->temp[itj][im1],fl1->mask,map_1_dum[im1]);
	//SHT[v*f^j]
	he_map2alm(fl1->nside,fl1->lmax,1,fl1->pol,map_1_dum,alm_1_dum,HE_NITER_DEFAULT);
	//C^abT*SHT[v*f^j]
	for(im2=0;im2<fl2->nmaps;im2++) {
	  for(im1=0;im1<fl1->nmaps;im1++)
	    he_alter_alm(lmax,-1.,alm_1_dum[im1],alm_2_dum[im2],cl_proposal[im1*fl2->nmaps+im2]);
	}
	//SHT^-1[C^abT*SHT[v*f^j]]
	he_alm2map(fl2->nside,fl2->lmax,1,fl2->pol,map_2_dum,alm_2_dum);
	//w*SHT^-1[C^abT*SHT[v*f^j]]
	for(im2=0;im2<fl2->nmaps;im2++)
	  he_map_product(fl2->nside,map_2_dum[im2],fl2->mask,map_2_dum[im2]);
	//SHT[w*SHT^-1[C^abT*SHT[v*f^j]]]
	he_map2alm(fl2->nside,fl2->lmax,1,fl2->pol,map_2_dum,alm_2_dum,HE_NITER_DEFAULT);
	//Sum_m(f^i*SHT[w*SHT^-1[C^abT*SHT[v*f^j]]]^*)/(2l+1)
	he_alm2cl(fl1->a_temp[iti],alm_2_dum,fl1->pol,fl2->pol,cl_dum,lmax);
	for(im1=0;im1<nspec;im1++) {
	  for(ip=0;ip<=lmax;ip++)
	    cl_bias[im1][ip]-=cl_dum[im1][ip]*mij;
	}
      }
    }
  }

  if((fl1->ntemp>0) && (fl2->ntemp>0)) {
    printf("Computing F4\n");
    int iti,itj,itp,itq,im1,im2;

    flouble *mat_prod=my_calloc(fl1->ntemp*fl2->ntemp,sizeof(flouble));
    for(itj=0;itj<fl1->ntemp;itj++) {
      for(itq=0;itq<fl2->ntemp;itq++) {
	//w*g^q
	for(im2=0;im2<fl2->nmaps;im2++)
	  he_map_product(fl2->nside,fl2->temp[itq][im2],fl2->mask,map_2_dum[im2]);
	//SHT[w*g^q]
	he_map2alm(fl2->nside,fl2->lmax,1,fl2->pol,map_2_dum,alm_2_dum,HE_NITER_DEFAULT);
	//C^ab*SHT[w*g^q]
	for(im1=0;im1<fl1->nmaps;im1++) {
	  for(im2=0;im2<fl2->nmaps;im2++)
	    he_alter_alm(lmax,-1.,alm_2_dum[im2],alm_1_dum[im1],cl_proposal[im1*fl2->nmaps+im2]);
	}
	//SHT^-1[C^ab*SHT[w*g^q]]
	he_alm2map(fl1->nside,fl1->lmax,1,fl1->pol,map_1_dum,alm_1_dum);
	for(im1=0;im1<fl1->nmaps;im1++) {
	  //v*SHT^-1[C^ab*SHT[w*g^q]]
	  he_map_product(fl1->nside,map_1_dum[im1],fl1->mask,map_1_dum[im1]);
	  //Int[f^jT*v*SHT^-1[C^ab*SHT[w*g^q]]]
	  mat_prod[itj*fl2->ntemp+itq]+=he_map_dot(fl1->nside,map_1_dum[im1],fl1->temp[itj][im1]);
	}
      }
    }

    for(iti=0;iti<fl1->ntemp;iti++) {
      for(itp=0;itp<fl2->ntemp;itp++) {
	//Sum_m(f^i*g^p*)/(2l+1)
	he_alm2cl(fl1->a_temp[iti],fl2->a_temp[itp],fl1->pol,fl2->pol,cl_dum,lmax);
	for(itj=0;itj<fl1->ntemp;itj++) {
	  double mij=gsl_matrix_get(fl1->matrix_M,iti,itj);
	  for(itq=0;itq<fl2->ntemp;itq++) {
	    double npq=gsl_matrix_get(fl2->matrix_M,itp,itq);
	    for(im1=0;im1<nspec;im1++) {
	      for(ip=0;ip<=lmax;ip++)
		cl_bias[im1][ip]+=cl_dum[im1][ip]*mat_prod[itj*fl2->ntemp+itq]*mij*npq;
	    }
	  }
	}
      }
    }

    free(mat_prod);
  }

  for(ii=0;ii<fl1->nmaps;ii++) {
    free(map_1_dum[ii]);
    free(alm_1_dum[ii]);
  }
  free(map_1_dum);
  free(alm_1_dum);
  for(ii=0;ii<fl2->nmaps;ii++) {
    free(map_2_dum[ii]);
    free(alm_2_dum[ii]);
  }
  free(map_2_dum);
  free(alm_2_dum);
  for(ii=0;ii<nspec;ii++)
    free(cl_dum[ii]);
  free(cl_dum);
}

void decouple_cl_l(MasterWorkspace *w,flouble **cl_in,flouble **cl_noise_in,flouble **cl_bias,flouble **cl_out)
{
  int icl,ib2,l2;
  gsl_vector *dl_map_bad_b=gsl_vector_alloc(w->ncls*w->bin->n_bands);
  gsl_vector *dl_map_good_b=gsl_vector_alloc(w->ncls*w->bin->n_bands);

  //Bin coupled power spectrum
  for(icl=0;icl<w->ncls;icl++) {
    for(ib2=0;ib2<w->bin->n_bands;ib2++) {
      int i2;
      double dl_b=0;
      for(i2=0;i2<w->bin->nell_list[ib2];i2++) {
	l2=w->bin->ell_list[ib2][i2];
	dl_b+=(cl_in[icl][l2]-cl_noise_in[icl][l2]-cl_bias[icl][l2])*weigh_l(l2)*w->bin->w_list[ib2][i2];
      }
      gsl_vector_set(dl_map_bad_b,w->ncls*ib2+icl,dl_b);
    }
  }

  printf("Solving for uncoupled Cls\n");
  gsl_linalg_LU_solve(w->coupling_matrix_binned,w->coupling_matrix_perm,dl_map_bad_b,dl_map_good_b);
  for(icl=0;icl<w->ncls;icl++) {
    for(ib2=0;ib2<w->bin->n_bands;ib2++)
      cl_out[icl][ib2]=gsl_vector_get(dl_map_good_b,w->ncls*ib2+icl);
  }

  gsl_vector_free(dl_map_bad_b);
  gsl_vector_free(dl_map_good_b);
}

MasterWorkspace *compute_power_spectra(Field *fl1,Field *fl2,BinSchm *bin,
				       flouble **cl_noise,flouble **cl_proposal,flouble **cl_out)
{
  int ii;
  flouble **cl_bias,**cl_data;
  MasterWorkspace *w=compute_coupling_matrix(fl1,fl2,bin);

  cl_bias=my_malloc(w->ncls*sizeof(flouble *));
  cl_data=my_malloc(w->ncls*sizeof(flouble *));
  for(ii=0;ii<w->ncls;ii++) {
    cl_bias[ii]=my_calloc((w->lmax+1),sizeof(flouble));
    cl_data[ii]=my_calloc((w->lmax+1),sizeof(flouble));
  }
  he_anafast(fl1->maps,fl2->maps,fl1->pol,fl2->pol,cl_data,fl1->nside,fl1->lmax);
  compute_deprojection_bias(fl1,fl2,cl_proposal,cl_bias);
  decouple_cl_l(w,cl_data,cl_noise,cl_bias,cl_out);
  for(ii=0;ii<w->ncls;ii++) {
    free(cl_bias[ii]);
    free(cl_data[ii]);
  }
  free(cl_bias);
  free(cl_data);

  return w;
}
