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

//Reads coupling matrix from file and computes its LU decomposition
void read_coupling_matrix(char *fname_in,int nbins_in,
			  gsl_matrix **coupling_matrix_b_out,
			  gsl_permutation **perm_out,
			  int pol1,int pol2)
{
  int sig,n_cl,stat;
  FILE *fi;
  gsl_permutation *perm;
  gsl_matrix *coupling_matrix_b;

  if(pol1) {
    if(pol2) n_cl=4;
    else n_cl=2;
  }
  else {
    if(pol2) n_cl=2;
    else n_cl=1;
  }
  perm=gsl_permutation_alloc(n_cl*nbins_in);
  coupling_matrix_b=gsl_matrix_alloc(n_cl*nbins_in,n_cl*nbins_in);

  fi=my_fopen(fname_in,"rb");
  stat=gsl_matrix_fread(fi,coupling_matrix_b);
  if(stat==GSL_EFAILED)
    report_error(1,"Error reading matrix from file %s\n",fname_in);
  fclose(fi);

  gsl_linalg_LU_decomp(coupling_matrix_b,perm,&sig);
  *coupling_matrix_b_out=coupling_matrix_b;
  *perm_out=perm;
}

//Computes binned coupling matrix
// fname_mask (in) : path to fits file containing the mask
// n_lbin (in) : number of multipoles per l-bin
// nside (out) : mask resolution
// lmax (out)  : maximum l used (3*nside-1)
// nbins (out) : number of l-bins
// coupling_matrix_b_out (out) : LU decomposition of the coupling matrix
// perm (out) : permutation used in the LU decomposition
// mask_out (out) : mask map
void compute_coupling_matrix(flouble *cl_mask,long nside_in,int lmax_in,BinSchm *bins,
			     gsl_matrix **coupling_matrix_b_out,
			     gsl_permutation **perm_out,
			     char *write_matrix,char *write_matrix_b,
			     int pol1,int pol2)
{
  int sig,ib2,ib3,l2,l3;
  double *cl_mask_bad_ub;
  double **coupling_matrix_ub;
  gsl_matrix *coupling_matrix_b;
  gsl_permutation *perm;
  int n_cl;

  if(pol1) {
    if(pol2) n_cl=4;
    else n_cl=2;
  }
  else {
    if(pol2) n_cl=2;
    else n_cl=1;
  }

  cl_mask_bad_ub=my_malloc((lmax_in+1)*sizeof(double));
  for(l2=0;l2<=lmax_in;l2++)
    cl_mask_bad_ub[l2]=cl_mask[l2]*(2*l2+1.);

  //Compute coupling matrix
  coupling_matrix_ub=my_malloc(n_cl*(lmax_in+1)*(sizeof(double *)));
  for(l2=0;l2<n_cl*(lmax_in+1);l2++)
    coupling_matrix_ub[l2]=my_calloc(n_cl*(lmax_in+1),(sizeof(double)));

  if((strcmp(write_matrix,"none")) && (access(write_matrix,F_OK)!=-1)) {
    int dum_n_cl,dum_lmax_in;
    FILE *fo=my_fopen(write_matrix,"rb");
    my_fread(&dum_n_cl,sizeof(int),1,fo);
    if(dum_n_cl!=n_cl)
      report_error(1,"Input covariance matrix has a wrong number of power spectra\n");
    my_fread(&dum_lmax_in,sizeof(int),1,fo);
    if(dum_lmax_in!=lmax_in)
      report_error(1,"Input covariance matrix has a wrong lmax\n");
    for(l2=0;l2<n_cl*(lmax_in+1);l2++)
      my_fread(coupling_matrix_ub[l2],sizeof(double),n_cl*(lmax_in+1),fo);
    fclose(fo);
  }
  else {
    printf("Computing unbinned coupling matrix\n");
#pragma omp parallel default(none)				\
  shared(lmax_in,cl_mask_bad_ub,coupling_matrix_ub,n_cl)
    {
      int ll2,ll3;
      double *wigner_00,*wigner_22;
      int lstart=0;
      
      if((n_cl==1) || (n_cl==2))
	wigner_00=my_malloc(2*(lmax_in+1)*sizeof(double));
      if((n_cl==2) || (n_cl==4))
	wigner_22=my_malloc(2*(lmax_in+1)*sizeof(double));
      
      if(n_cl>1)
	lstart=2;
      
#pragma omp for schedule(dynamic)
      for(ll2=lstart;ll2<=lmax_in;ll2++) {
	for(ll3=lstart;ll3<=lmax_in;ll3++) {
	  int jj,l1,lmin_here,lmax_here;
	  int lmin_here_00=0,lmax_here_00=2*(lmax_in+1)+1;
	  int lmin_here_22=0,lmax_here_22=2*(lmax_in+1)+1;
	  
	  if((n_cl==1) || (n_cl==2))
	    drc3jj(ll2,ll3,0,0,&lmin_here_00,&lmax_here_00,wigner_00,2*(lmax_in+1));
	  if((n_cl==2) || (n_cl==4))
	    drc3jj(ll2,ll3,2,-2,&lmin_here_22,&lmax_here_22,wigner_22,2*(lmax_in+1));
	  
	  lmin_here=COM_MAX(lmin_here_00,lmin_here_22);
	  lmax_here=COM_MIN(lmax_here_00,lmax_here_22);
	  
	  for(l1=lmin_here;l1<=lmax_here;l1++) {
	    if(l1<=lmax_in) {
	      double wfac;
	      int j00=l1-lmin_here_00;
	      int j22=l1-lmin_here_22;
	      if(n_cl==1) {
		wfac=cl_mask_bad_ub[l1]*wigner_00[j00]*wigner_00[j00];
		coupling_matrix_ub[1*ll2+0][1*ll3+0]+=wfac; //TT,TT
	      }
	      if(n_cl==2) {
		wfac=cl_mask_bad_ub[l1]*wigner_00[j00]*wigner_22[j22];
		coupling_matrix_ub[2*ll2+0][2*ll3+0]+=wfac; //TE,TE
		coupling_matrix_ub[2*ll2+1][2*ll3+1]+=wfac; //TB,TB
	      }
	      if(n_cl==4) {
		int suml=l1+ll2+ll3;
		wfac=cl_mask_bad_ub[l1]*wigner_22[j22]*wigner_22[j22];
		if(suml & 1) { //Odd sum
		  coupling_matrix_ub[4*ll2+0][4*ll3+3]+=wfac; //EE,BB
		  coupling_matrix_ub[4*ll2+1][4*ll3+2]-=wfac; //EB,BE
		  coupling_matrix_ub[4*ll2+2][4*ll3+1]-=wfac; //BE,EB
		  coupling_matrix_ub[4*ll2+3][4*ll3+0]+=wfac; //BB,EE
		}
		else {
		  coupling_matrix_ub[4*ll2+0][4*ll3+0]+=wfac; //EE,EE
		  coupling_matrix_ub[4*ll2+1][4*ll3+1]+=wfac; //EB,EB
		  coupling_matrix_ub[4*ll2+2][4*ll3+2]+=wfac; //BE,BE
		  coupling_matrix_ub[4*ll2+3][4*ll3+3]+=wfac; //BB,BB
		}
	      }
	    }
	  }
	  for(jj=0;jj<n_cl;jj++) {
	    int kk;
	    for(kk=0;kk<n_cl;kk++)
	      coupling_matrix_ub[n_cl*ll2+jj][n_cl*ll3+kk]*=(2*ll3+1.)/(4*M_PI);
	  }
	}
      } //end omp for
      if((n_cl==1) || (n_cl==2))
	free(wigner_00);
      if((n_cl==2) || (n_cl==4))
	free(wigner_22);
    } //end omp parallel
  }

  free(cl_mask_bad_ub);

  int nbins_in=bins->n_bands;
  coupling_matrix_b=gsl_matrix_alloc(n_cl*nbins_in,n_cl*nbins_in);
  int icl_a;
  for(icl_a=0;icl_a<n_cl;icl_a++) {
    int icl_b;
    for(icl_b=0;icl_b<n_cl;icl_b++) {
      for(ib2=0;ib2<nbins_in;ib2++) {
	for(ib3=0;ib3<nbins_in;ib3++) {
	  int i2;
	  double coupling_b=0;
	  for(i2=0;i2<bins->nell_list[ib2];i2++) {
	    int i3;
	    l2=bins->ell_list[ib2][i2];
	    for(i3=0;i3<bins->nell_list[ib3];i3++) {
	      l3=bins->ell_list[ib3][i3];
	      coupling_b+=coupling_matrix_ub[n_cl*l2+icl_a][n_cl*l3+icl_b]*bins->w_list[ib2][i2]*
		weigh_l(l2)/weigh_l(l3);
	    }
	  }
	  gsl_matrix_set(coupling_matrix_b,n_cl*ib2+icl_a,n_cl*ib3+icl_b,coupling_b);
	}
      }
    }
  }

  if((strcmp(write_matrix,"none")) && (access(write_matrix,F_OK)!=-1)) {
    FILE *fo=my_fopen(write_matrix,"wb");
    my_fwrite(&n_cl,sizeof(int),1,fo);
    my_fwrite(&lmax_in,sizeof(int),1,fo);
    for(l2=0;l2<n_cl*(lmax_in+1);l2++)
      my_fwrite(coupling_matrix_ub[l2],sizeof(double),n_cl*(lmax_in+1),fo);
    fclose(fo);
  }

  //Write matrix if required
  if(strcmp(write_matrix_b,"none")) {
    FILE *fo=my_fopen(write_matrix_b,"wb");
    gsl_matrix_fwrite(fo,coupling_matrix_b);
    fclose(fo);
  }

  for(l2=0;l2<n_cl*(lmax_in+1);l2++)
    free(coupling_matrix_ub[l2]);
  free(coupling_matrix_ub);

  if(flag_do_cov) {
    for(l2=0;l2<(lmax_in+1);l2++) {
      free(coupling_matrix_a_ub[l2]);
      free(coupling_matrix_b_ub[l2]);
    }
    free(coupling_matrix_a_ub);
    free(coupling_matrix_b_ub);
  }

  //Get LU decomposition
  perm=gsl_permutation_alloc(n_cl*nbins_in);
  gsl_linalg_LU_decomp(coupling_matrix_b,perm,&sig);
  *coupling_matrix_b_out=coupling_matrix_b;
  *perm_out=perm;
}

//Decouples unbinned cls (into binned dl)
// cl_in (in)  : input coupled cl
// nbins (in)  : number of l-bins
// n_lbin (in) : number of ls per bin
// coupling_matrix_b (in) : LU decomposition of the coupling matrix
// perm (int) : permutation used in the LU decomposition
// returns decoupled D_l
flouble **decouple_cl_l(flouble **cl_in,flouble **cl_noise_in,
			int n_cl,BinSchm *bins,
			gsl_matrix *coupling_matrix_b,gsl_permutation *perm)
{
  int icl,ib2,l2;
  int nbins=bins->n_bands;
  gsl_vector *dl_map_bad_b=gsl_vector_alloc(n_cl*nbins);
  gsl_vector *dl_map_good_b=gsl_vector_alloc(n_cl*nbins);
  flouble **dl_out=my_malloc(n_cl*sizeof(flouble *));

  for(icl=0;icl<n_cl;icl++)
    dl_out[icl]=my_malloc(nbins*sizeof(flouble));

  //Bin coupled power spectrum
  for(icl=0;icl<n_cl;icl++) {
    for(ib2=0;ib2<nbins;ib2++) {
      int i2;
      double dl_b=0;
      for(i2=0;i2<bins->nell_list[ib2];i2++) {
	l2=bins->ell_list[ib2][i2];
	dl_b+=(cl_in[icl][l2]-cl_noise_in[icl][l2])*weigh_l(l2)*bins->w_list[ib2][i2];
      }
      gsl_vector_set(dl_map_bad_b,n_cl*ib2+icl,dl_b);
    }
  }

  printf("Solving for uncoupled Cls\n");
  gsl_linalg_LU_solve(coupling_matrix_b,perm,
		      dl_map_bad_b,dl_map_good_b);
  for(icl=0;icl<n_cl;icl++) {
    for(ib2=0;ib2<nbins;ib2++)
      dl_out[icl][ib2]=gsl_vector_get(dl_map_good_b,n_cl*ib2+icl);
  }

  gsl_vector_free(dl_map_bad_b);
  gsl_vector_free(dl_map_good_b);

  return dl_out;
}
