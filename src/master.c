#include "common.h"

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
    dam_report_error(1,"Wrong arguments: %d %d %d %d\n",il2,il3,im2,im3);
  
  //l1 bounds
  l1max=il2+il3;
  l1min=DAM_MAX((abs(il2-il3)),(abs(im1)));
  *l1max_out=l1max;
  *l1min_out=l1min;
  
  if(l1max-l1min<0) //Check for meaningful values
    dam_report_error(0,"WTF?\n");
  
  if(l1max==l1min) { //If it's only one value:
    thrcof[0]=sign2/sqrt(l1min+l2+l3+1);
    return 0;
  }
  else {
    nfin=l1max-l1min+1;
    if(nfin>size) //Check there's enough space
      dam_report_error(1,"Output array is too small %d\n",nfin);
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
			  gsl_permutation **perm_out)
{
  int sig;
  FILE *fi=dam_fopen(fname_in,"rb");
  gsl_permutation *perm=gsl_permutation_alloc(nbins_in);
  gsl_matrix *coupling_matrix_b=gsl_matrix_alloc(nbins_in,nbins_in);
  int stat=gsl_matrix_fread(fi,coupling_matrix_b);
  if(stat==GSL_EFAILED)
    dam_report_error(1,"Error reading matrix from file %s\n",fname_in);
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
void compute_coupling_matrix(flouble *cl_mask,int n_lbin,
			     long nside_in,int lmax_in,int nbins_in,
			     gsl_matrix **coupling_matrix_b_out,
			     gsl_permutation **perm_out,char *write_matrix)
{
  int sig,ib2,ib3,l2,l3;
  double *cl_mask_bad_ub,**coupling_matrix_ub;
  gsl_matrix *coupling_matrix_b;
  gsl_permutation *perm;

  cl_mask_bad_ub=(double *)dam_malloc((lmax_in+1)*sizeof(double));
  for(l2=0;l2<=lmax_in;l2++)
    cl_mask_bad_ub[l2]=cl_mask[l2]*(2*l2+1.);

  //Compute coupling matrix
  coupling_matrix_ub=(double **)dam_malloc((lmax_in+1)*(sizeof(double *)));
  for(l2=0;l2<=lmax_in;l2++)
    coupling_matrix_ub[l2]=(double *)dam_malloc((lmax_in+1)*(sizeof(double)));

  printf("Computing unbinned coupling matrix\n");
#pragma omp parallel default(none)			\
  shared(lmax_in,cl_mask_bad_ub,coupling_matrix_ub)
  {
    int ll2,ll3;
    double *wigner_dum=(double *)dam_malloc(2*(lmax_in+1)*sizeof(double));

#pragma omp for schedule(dynamic)
    for(ll2=0;ll2<=lmax_in;ll2++) {
      for(ll3=0;ll3<=lmax_in;ll3++) {
	int jj,lmin_here,lmax_here,size_here;
	double coupling_here=0;
	drc3jj(ll2,ll3,0,0,&lmin_here,&lmax_here,wigner_dum,2*(lmax_in+1));
	size_here=lmax_here-lmin_here+1;
	for(jj=0;jj<size_here;jj++) {
	  int l1=lmin_here+jj;
	  if(l1<=lmax_in) {
	    double w=wigner_dum[jj];
	    coupling_here+=cl_mask_bad_ub[l1]*w*w;
	  }
	}
	coupling_matrix_ub[ll2][ll3]=(2*ll3+1.)*coupling_here/(4*M_PI);
      }
    } //end omp for
    free(wigner_dum);
  } //end omp parallel
  free(cl_mask_bad_ub);

  coupling_matrix_b=gsl_matrix_alloc(nbins_in,nbins_in);
  for(ib2=0;ib2<nbins_in;ib2++) {
    int l20=2+ib2*n_lbin;
    for(ib3=0;ib3<nbins_in;ib3++) {
      int l30=2+ib3*n_lbin;
      double coupling_b=0;
      for(l2=l20;l2<l20+n_lbin;l2++) {
	for(l3=l30;l3<l30+n_lbin;l3++) {
	  coupling_b+=coupling_matrix_ub[l2][l3]*
	    ((double)(l2*(l2+1))/(l3*(l3+1)));
	}
      }
      coupling_b/=n_lbin;
      gsl_matrix_set(coupling_matrix_b,ib2,ib3,coupling_b);
    }
  }
  for(l2=0;l2<=lmax_in;l2++)
    free(coupling_matrix_ub[l2]);
  free(coupling_matrix_ub);

  //Write matrix if required
  if(strcmp(write_matrix,"none")) {
    FILE *fo=dam_fopen(write_matrix,"wb");
    gsl_matrix_fwrite(fo,coupling_matrix_b);
    fclose(fo);
  }

  //Get LU decomposition
  perm=gsl_permutation_alloc(nbins_in);
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
flouble *decouple_cl_l(flouble *cl_in,flouble *cl_noise_in,int nbins,int n_lbin,
		       gsl_matrix *coupling_matrix_b,gsl_permutation *perm)
{
  int ib2,l2;
  gsl_vector *dl_map_bad_b=gsl_vector_alloc(nbins);
  gsl_vector *dl_map_good_b=gsl_vector_alloc(nbins);
  flouble *dl_out=(flouble *)dam_malloc(nbins*sizeof(flouble));

  //Bin coupled power spectrum
  for(ib2=0;ib2<nbins;ib2++) {
    int l20=2+ib2*n_lbin;
    double dl_b=0;
    for(l2=l20;l2<l20+n_lbin;l2++)
      dl_b+=(cl_in[l2]-cl_noise_in[l2])*((double)(l2*(l2+1))/(2*M_PI));;
    dl_b/=n_lbin;
    gsl_vector_set(dl_map_bad_b,ib2,dl_b);
  }

  printf("Solving for uncoupled Cls\n");
  gsl_linalg_LU_solve(coupling_matrix_b,perm,
		      dl_map_bad_b,dl_map_good_b);
  for(ib2=0;ib2<nbins;ib2++)
    dl_out[ib2]=gsl_vector_get(dl_map_good_b,ib2);

  gsl_vector_free(dl_map_bad_b);
  gsl_vector_free(dl_map_good_b);
  

  return dl_out;
}
