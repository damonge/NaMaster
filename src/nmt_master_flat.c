#include "utils.h"

void nmt_workspace_flat_free(nmt_workspace_flat *w)
{
  int ii;
  gsl_permutation_free(w->coupling_matrix_perm);
  gsl_matrix_free(w->coupling_matrix_binned);
  nmt_bins_flat_free(w->bin);
  nmt_flatsky_info_free(w->fs);
  for(ii=0;ii<w->ncls*w->nells;ii++)
    free(w->coupling_matrix_unbinned[ii]);
  free(w->coupling_matrix_unbinned);
  free(w->pcl_masks);
  free(w->l_arr);
  free(w->i_band);
  free(w->n_cells);
  free(w);
}

static nmt_workspace_flat *nmt_workspace_flat_new(int nl_rebin,int ncls,nmt_flatsky_info *fs,
						  nmt_binning_scheme_flat *bin)
{
  int ii,ib=0;
  nmt_workspace_flat *w=my_malloc(sizeof(nmt_workspace_flat));
  w->ncls=ncls;
  w->nl_rebin=nl_rebin;

  w->bin=nmt_bins_flat_create(bin->n_bands,bin->ell_0_list,bin->ell_f_list);
  w->lmax=w->bin->ell_f_list[w->bin->n_bands-1];

  w->fs=nmt_flatsky_info_alloc(fs->nx,fs->ny,fs->lx,fs->ly);

  if(w->lmax>w->fs->ell_min[w->fs->n_ell-1]+w->fs->dell)
    report_error(1,"You requested multipoles beyond the resolution of the input maps.\n");

  w->nells=0;
  for(ii=0;ii<w->fs->n_ell;ii++) {
    if(w->fs->ell_min[ii]<w->lmax)
      w->nells+=w->nl_rebin;
    else
      break;
  }

  w->pcl_masks=my_malloc(w->nells*sizeof(flouble));
  w->l_arr=my_malloc(w->nells*sizeof(flouble));
  w->i_band=my_malloc(w->nells*sizeof(int));
  w->n_cells=my_calloc(w->nells,sizeof(int));
  for(ii=0;ii<w->nells;ii++) {
    flouble l_here=ii*w->fs->dell/w->nl_rebin;
    ib=nmt_bins_flat_search_fast(w->bin,l_here,ib);
    w->l_arr[ii] =l_here;
    w->i_band[ii]=ib;
  }

  w->coupling_matrix_unbinned=my_malloc(w->ncls*w->nells*sizeof(flouble *));
  for(ii=0;ii<w->ncls*w->nells;ii++)
    w->coupling_matrix_unbinned[ii]=my_calloc(w->ncls*w->nells,sizeof(flouble));

  w->coupling_matrix_binned=gsl_matrix_alloc(w->ncls*w->bin->n_bands,w->ncls*w->bin->n_bands);
  w->coupling_matrix_perm=gsl_permutation_alloc(w->ncls*w->bin->n_bands);

  return w;
}

nmt_workspace_flat *nmt_workspace_flat_read(char *fname)
{
  int ii,nx,ny;
  flouble lx,ly;
  nmt_workspace_flat *w=my_malloc(sizeof(nmt_workspace_flat));
  FILE *fi=my_fopen(fname,"rb");

  my_fread(&(w->ncls),sizeof(int),1,fi);
  my_fread(&(w->nl_rebin),sizeof(int),1,fi);
  my_fread(&(w->nells),sizeof(int),1,fi);

  w->pcl_masks=my_malloc(w->nells*sizeof(flouble));
  my_fread(w->pcl_masks,sizeof(flouble),w->nells,fi);
  w->l_arr=my_malloc(w->nells*sizeof(flouble));
  my_fread(w->l_arr,sizeof(flouble),w->nells,fi);
  w->i_band=my_malloc(w->nells*sizeof(flouble));
  my_fread(w->i_band,sizeof(int),w->nells,fi);
  w->n_cells=my_malloc(w->nells*sizeof(flouble));
  my_fread(w->n_cells,sizeof(int),w->nells,fi);

  w->coupling_matrix_unbinned=my_malloc(w->ncls*w->nells*sizeof(flouble *));
  for(ii=0;ii<w->ncls*w->nells;ii++) {
    w->coupling_matrix_unbinned[ii]=my_malloc(w->ncls*w->nells*sizeof(flouble));
    my_fread(w->coupling_matrix_unbinned[ii],sizeof(flouble),w->ncls*w->nells,fi);
  }

  my_fread(&nx,sizeof(int),1,fi);
  my_fread(&ny,sizeof(int),1,fi);
  my_fread(&lx,sizeof(flouble),1,fi);
  my_fread(&ly,sizeof(flouble),1,fi);
  w->fs=nmt_flatsky_info_alloc(nx,ny,lx,ly);

  w->bin=my_malloc(sizeof(nmt_binning_scheme_flat));
  my_fread(&(w->bin->n_bands),sizeof(int),1,fi);
  w->bin->ell_0_list=my_malloc(w->bin->n_bands*sizeof(flouble));
  w->bin->ell_f_list=my_malloc(w->bin->n_bands*sizeof(flouble));
  my_fread(w->bin->ell_0_list,sizeof(flouble),w->bin->n_bands,fi);
  my_fread(w->bin->ell_f_list,sizeof(flouble),w->bin->n_bands,fi);

  w->coupling_matrix_binned=gsl_matrix_alloc(w->ncls*w->bin->n_bands,w->ncls*w->bin->n_bands);
  w->coupling_matrix_perm=gsl_permutation_alloc(w->ncls*w->bin->n_bands);
  gsl_matrix_fread(fi,w->coupling_matrix_binned);
  gsl_permutation_fread(fi,w->coupling_matrix_perm);
  
  fclose(fi);

  return w;
}

void nmt_workspace_flat_write(nmt_workspace_flat *w,char *fname)
{
  int ii;
  FILE *fo=my_fopen(fname,"wb");

  my_fwrite(&(w->ncls),sizeof(int),1,fo);
  my_fwrite(&(w->nl_rebin),sizeof(int),1,fo);
  my_fwrite(&(w->nells),sizeof(int),1,fo);

  my_fwrite(w->pcl_masks,sizeof(flouble),w->nells,fo);
  my_fwrite(w->l_arr,sizeof(flouble),w->nells,fo);
  my_fwrite(w->i_band,sizeof(int),w->nells,fo);
  my_fwrite(w->n_cells,sizeof(int),w->nells,fo);

  for(ii=0;ii<w->ncls*w->nells;ii++)
    my_fwrite(w->coupling_matrix_unbinned[ii],sizeof(flouble),w->ncls*w->nells,fo);

  my_fwrite(&(w->fs->nx),sizeof(int),1,fo);
  my_fwrite(&(w->fs->ny),sizeof(int),1,fo);
  my_fwrite(&(w->fs->lx),sizeof(flouble),1,fo);
  my_fwrite(&(w->fs->ly),sizeof(flouble),1,fo);

  my_fwrite(&(w->bin->n_bands),sizeof(int),1,fo);
  my_fwrite(w->bin->ell_0_list,sizeof(flouble),w->bin->n_bands,fo);
  my_fwrite(w->bin->ell_f_list,sizeof(flouble),w->bin->n_bands,fo);

  gsl_matrix_fwrite(fo,w->coupling_matrix_binned);
  gsl_permutation_fwrite(fo,w->coupling_matrix_perm);
  
  fclose(fo);
}

static void bin_coupling_matrix(nmt_workspace_flat *w)
{
  int icl_a,icl_b,il2,il3,ib2,ib3,sig;
  int *ncells_per_band=my_calloc(w->bin->n_bands,sizeof(int));
  flouble **coupling_b=my_malloc(w->ncls*w->bin->n_bands*sizeof(flouble *));
  for(ib2=0;ib2<w->ncls*w->bin->n_bands;ib2++)
    coupling_b[ib2]=my_calloc(w->ncls*w->bin->n_bands,sizeof(flouble *));

  for(il2=0;il2<w->nells;il2++) {
    ib2=w->i_band[il2];
    if(ib2>=0)
      ncells_per_band[ib2]+=w->n_cells[il2];
  }
  
  for(icl_a=0;icl_a<w->ncls;icl_a++) {
    for(icl_b=0;icl_b<w->ncls;icl_b++) {
      for(il2=0;il2<w->nells;il2++) {
	ib2=w->i_band[il2];
	if(ib2>=0) {
	  for(il3=0;il3<w->nells;il3++) {
	    ib3=w->i_band[il3];
	    if(ib3>=0) {
	      coupling_b[w->ncls*ib2+icl_a][w->ncls*ib3+icl_b]+=
		w->n_cells[il2]*w->coupling_matrix_unbinned[w->ncls*il2+icl_a][w->ncls*il3+icl_b];
	    }
	  }
	}
      }
      for(ib2=0;ib2<w->bin->n_bands;ib2++) {
	flouble norm=0;
	if(ncells_per_band[ib2]>0)
	  norm=1./((double)(ncells_per_band[ib2]));
	for(ib3=0;ib3<w->bin->n_bands;ib3++) {
	  flouble c=coupling_b[w->ncls*ib2+icl_a][w->ncls*ib3+icl_b]*norm;
	  gsl_matrix_set(w->coupling_matrix_binned,w->ncls*ib2+icl_a,w->ncls*ib3+icl_b,c);
	}
      }
    }
  }

  /*
  FILE *log;
  int di;
  log=fopen("log_b.log","w");
  for(di=0;di<w->ncls*w->bin->n_bands;di++) {
    int dj;
    for(dj=0;dj<w->ncls*w->bin->n_bands;dj++)
      fprintf(log,"%.3lf ",gsl_matrix_get(w->coupling_matrix_binned,di,dj));
    fprintf(log,"\n");
  }
  fclose(log);

  log=fopen("log_u.log","w");
  for(di=0;di<w->ncls*w->nells;di++) {
    int dj;
    for(dj=0;dj<w->ncls*w->nells;dj++)
      fprintf(log,"%.3lf ",w->coupling_matrix_unbinned[di][dj]);
    fprintf(log,"\n");
  }
  fclose(log);
  */

  gsl_linalg_LU_decomp(w->coupling_matrix_binned,w->coupling_matrix_perm,&sig);

  for(ib2=0;ib2<w->ncls*w->bin->n_bands;ib2++)
    free(coupling_b[ib2]);
  free(coupling_b);
  free(ncells_per_band);
}

typedef struct {
  flouble l,k;
  flouble *carr;
  flouble inv_dell;
  flouble qmax;
  int flag;
} intg_par;

#define INTG_FLAG_00 100
#define INTG_FLAG_02 101
#define INTG_FLAG_22 102
static double integ_coupling(double th,void *params)
{
  intg_par *p=(intg_par *)params;
  double cth=cos(th);
  double q=sqrt(p->k*p->k+p->l*p->l-2*p->k*p->l*cth);
  double cvw;

  if(q<0) cvw=p->carr[0];
  else if(q>=p->qmax) cvw=0;
  else cvw=p->carr[(int)(q*p->inv_dell)];

  if(p->flag==INTG_FLAG_00)
    return cvw;
  else {
    flouble c2th=2*cth*cth-1;
    if(p->flag==INTG_FLAG_02)
      return cvw*c2th;
    else if(p->flag==INTG_FLAG_22) {
      flouble c4th=2*c2th*c2th-1;
      return cvw*c4th;
    }
    else {
      report_error(1,"Unknown flag %d\n",p->flag);
      return 0;
    }
  }
}

static int check_flatsky_infos(nmt_flatsky_info *fs1,nmt_flatsky_info *fs2)
{
  if(fs1->nx!=fs2->nx) return 1;
  if(fs1->ny!=fs2->ny) return 1;
  if(fs1->lx!=fs2->lx) return 1;
  if(fs1->ly!=fs2->ly) return 1;
  return 0;
}

static nmt_workspace_flat *nmt_compute_coupling_matrix_flat_a(nmt_field_flat *fl1,nmt_field_flat *fl2,
							      nmt_binning_scheme_flat *bin,int nl_rebin)
{
  int i2;
  int n_cl=fl1->nmaps*fl2->nmaps;
  flouble *beam_prod;
  flouble *pcl_masks_raw=my_malloc(fl1->fs->n_ell*sizeof(flouble));
  nmt_workspace_flat *w=nmt_workspace_flat_new(nl_rebin,n_cl,fl1->fs,bin);
  flouble dell_here=w->fs->dell/w->nl_rebin;

  fs_anafast(fl1->fs,&(fl1->mask),&(fl2->mask),0,0,&(pcl_masks_raw));

  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  beam_prod=my_malloc(w->nells*sizeof(flouble));
  for(i2=0;i2<w->nells;i2++) {
    int isub=i2/w->nl_rebin;
    flouble l_here=w->l_arr[i2];
    flouble b1=nmt_k_function_eval(fl1->beam,l_here,intacc);
    flouble b2=nmt_k_function_eval(fl2->beam,l_here,intacc);
    beam_prod[i2]=b1*b2;
    w->n_cells[i2]=1;
    w->pcl_masks[i2]=dell_here*pcl_masks_raw[isub]*l_here/(4*M_PI*M_PI);
  }
  gsl_interp_accel_free(intacc);
  free(pcl_masks_raw);

#pragma omp parallel default(none)		\
  shared(w,beam_prod,dell_here)
  {
    int jj,ll1,ll2,ll3;

#pragma omp for
    for(ll1=0;ll1<w->nells;ll1++) {
      flouble l1=w->l_arr[ll1];
      for(ll2=0;ll2<w->nells;ll2++) {
	flouble l2=w->l_arr[ll2];
	for(ll3=0;ll3<w->nells;ll3++) {
	  flouble wfac;
	  flouble l3=w->l_arr[ll3];
	  if(l3<=fabs(l1-l2))
	    continue;
	  else if(l3>=(l1+l2))
	    break;
	  flouble j123,fac_spin=1.,prod=(l1+l2+l3)*(l1+l2-l3)*(l1-l2+l3)*(-l1+l2+l3);
	  if(prod<=0)
	    report_error(1,"This shouldn't have happened\n");
	  
	  j123=4./sqrt(prod);
	  if((w->ncls==2)||(w->ncls==4)) {
	    if((l1==0) || (l2==0))
	      fac_spin=0.;
	    else {
	      flouble l12=l1*l1,l22=l2*l2,l32=l3*l3;
	      flouble l14=l12*l12,l24=l22*l22,l34=l32*l32;
	      fac_spin=(l14+l24+l34-2*l12*l32-2*l22*l32)/(2*l12*l22);
	      if(w->ncls==4)
		fac_spin=fac_spin*fac_spin;
	    }
	  }

	  if(w->ncls==1) {
	    wfac=w->pcl_masks[ll3]*j123;
	    w->coupling_matrix_unbinned[1*ll1+0][1*ll2+0]+=wfac; //TT,TT
	  }
	  if(w->ncls==2) {
	    wfac=w->pcl_masks[ll3]*j123*fac_spin;
	    w->coupling_matrix_unbinned[2*ll1+0][2*ll2+0]+=wfac; //TE,TE
	    w->coupling_matrix_unbinned[2*ll1+1][2*ll2+1]+=wfac; //TB,TB
	  }
	  if(w->ncls==4) {
	    wfac=w->pcl_masks[ll3]*j123*fac_spin;
	    w->coupling_matrix_unbinned[4*ll1+0][4*ll2+0]+=wfac; //EE,EE
	    w->coupling_matrix_unbinned[4*ll1+1][4*ll2+1]+=wfac; //EE,EE
	    w->coupling_matrix_unbinned[4*ll1+2][4*ll2+2]+=wfac; //EE,EE
	    w->coupling_matrix_unbinned[4*ll1+3][4*ll2+3]+=wfac; //EE,EE
	    wfac=w->pcl_masks[ll3]*j123*(1-fac_spin);
	    w->coupling_matrix_unbinned[4*ll1+0][4*ll2+3]+=wfac; //EE,BB
	    w->coupling_matrix_unbinned[4*ll1+1][4*ll2+2]-=wfac; //EE,EE
	    w->coupling_matrix_unbinned[4*ll1+2][4*ll2+1]-=wfac; //EE,EE
	    w->coupling_matrix_unbinned[4*ll1+3][4*ll2+0]+=wfac; //EE,EE
	  }
	}
	for(jj=0;jj<w->ncls;jj++) {
	  int kk;
	  for(kk=0;kk<w->ncls;kk++)
	    w->coupling_matrix_unbinned[w->ncls*ll1+jj][w->ncls*ll2+kk]*=l2*beam_prod[ll2]*dell_here*0.5;
	}
      }
    } //end omp for
  } //end omp parallel

  bin_coupling_matrix(w);

  free(beam_prod);

  return w;
}

static nmt_workspace_flat *nmt_compute_coupling_matrix_flat_b(nmt_field_flat *fl1,nmt_field_flat *fl2,
							      nmt_binning_scheme_flat *bin,int nl_rebin)
{
  int ii,n_cl=fl1->nmaps*fl2->nmaps;
  nmt_workspace_flat *w=nmt_workspace_flat_new(nl_rebin,n_cl,fl1->fs,bin);

  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  flouble *beam_prod=my_malloc(w->nells*sizeof(flouble));
  for(ii=0;ii<w->nells;ii++) {
    flouble l_here=w->l_arr[ii];
    flouble b1=nmt_k_function_eval(fl1->beam,l_here,intacc);
    flouble b2=nmt_k_function_eval(fl2->beam,l_here,intacc);
    beam_prod[ii]=b1*b2;
    w->n_cells[ii]=1;
  }
  gsl_interp_accel_free(intacc);

  flouble *pcl_masks_raw=my_malloc(fl1->fs->n_ell*sizeof(flouble));
  fs_anafast(fl1->fs,&(fl1->mask),&(fl2->mask),0,0,&(pcl_masks_raw));

  gsl_set_error_handler_off();

#pragma omp parallel default(none)		\
  shared(w,beam_prod,pcl_masks_raw)
  {
    int ll1;
    gsl_function F;
    intg_par p;
    flouble dell_here=w->fs->dell/w->nl_rebin;
    gsl_integration_workspace *wi=gsl_integration_workspace_alloc(1000);
    p.qmax=w->fs->dell*w->fs->n_ell;
    p.carr=pcl_masks_raw;
    p.inv_dell=w->fs->i_dell;
    F.function=&integ_coupling;
    F.params=&p;

#pragma omp for
    for(ll1=0;ll1<w->nells;ll1++) {
      int ll2;
      flouble l1=w->l_arr[ll1];
      p.l=l1;
      for(ll2=0;ll2<w->nells;ll2++) {
	flouble l2=w->l_arr[ll2];
	p.k=l2;
	if(w->ncls==1) {
	  double integ,einteg;
	  p.flag=INTG_FLAG_00;
	  gsl_integration_qag(&F,0,M_PI,0,1E-4,1000,GSL_INTEG_GAUSS61,wi,&integ,&einteg);
	  w->coupling_matrix_unbinned[1*ll1+0][1*ll2+0]=integ; //TT,TT
	}
	if(w->ncls==2) {
	  double integ,einteg;
	  p.flag=INTG_FLAG_02;
	  gsl_integration_qag(&F,0,M_PI,0,1E-4,1000,GSL_INTEG_GAUSS61,wi,&integ,&einteg);
	  w->coupling_matrix_unbinned[2*ll1+0][2*ll2+0]=integ; //TE,TE
	  w->coupling_matrix_unbinned[2*ll1+1][2*ll2+1]=integ; //TB,TB
	}
	if(w->ncls==4) {
	  double integ_0,integ_4,einteg,integ_p,integ_m;
	  p.flag=INTG_FLAG_00;
	  gsl_integration_qag(&F,0,M_PI,0,1E-4,1000,GSL_INTEG_GAUSS61,wi,&integ_0,&einteg);
	  p.flag=INTG_FLAG_22;
	  gsl_integration_qag(&F,0,M_PI,0,1E-4,1000,GSL_INTEG_GAUSS61,wi,&integ_4,&einteg);
	  integ_p=0.5*(integ_0+integ_4);
	  integ_m=0.5*(integ_0-integ_4);
	  w->coupling_matrix_unbinned[4*ll1+0][4*ll2+0]= integ_p; //EE,EE
	  w->coupling_matrix_unbinned[4*ll1+1][4*ll2+1]= integ_p; //EB,EB
	  w->coupling_matrix_unbinned[4*ll1+2][4*ll2+2]= integ_p; //BE,BE
	  w->coupling_matrix_unbinned[4*ll1+3][4*ll2+3]= integ_p; //BB,BB
	  w->coupling_matrix_unbinned[4*ll1+0][4*ll2+3]= integ_m; //EE,BB
	  w->coupling_matrix_unbinned[4*ll1+1][4*ll2+2]=-integ_m; //EB,BE
	  w->coupling_matrix_unbinned[4*ll1+2][4*ll2+1]=-integ_m; //BE,EB
	  w->coupling_matrix_unbinned[4*ll1+3][4*ll2+0]= integ_m; //BB,EE
	}

	int jj,kk;
	flouble prefac=l2*beam_prod[ll2]*dell_here/(2*2*M_PI*M_PI);
	for(jj=0;jj<w->ncls;jj++) {
	  for(kk=0;kk<w->ncls;kk++)
	    w->coupling_matrix_unbinned[w->ncls*ll1+jj][w->ncls*ll2+kk]*=prefac;
	}
      }
    } //end omp for
    gsl_integration_workspace_free(wi);
  } //end omp parallel

  bin_coupling_matrix(w);

  free(beam_prod);
  free(pcl_masks_raw);

  return w;
}

static nmt_workspace_flat *nmt_compute_coupling_matrix_flat_c(nmt_field_flat *fl1,nmt_field_flat *fl2,
							      nmt_binning_scheme_flat *bin,int nl_rebin)
{
  int ii;
  nmt_workspace_flat *w=nmt_workspace_flat_new(1,fl1->nmaps*fl2->nmaps,fl1->fs,bin);
  nmt_flatsky_info *fs=fl1->fs;
  
  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  flouble *beam_prod=my_malloc(w->nells*sizeof(flouble));
  for(ii=0;ii<w->nells;ii++) {
    flouble l_here=w->l_arr[ii];
    flouble b1=nmt_k_function_eval(fl1->beam,l_here,intacc);
    flouble b2=nmt_k_function_eval(fl2->beam,l_here,intacc);
    beam_prod[ii]=b1*b2;
  }
  gsl_interp_accel_free(intacc);

  fcomplex *cmask1,*cmask2;
  flouble *maskprod,*cosarr,*sinarr;
  int *i_band;
  cmask1=dftw_malloc(fl1->fs->ny*(fs->nx/2+1)*sizeof(fcomplex));
  fs_map2alm(fl1->fs,1,0,&(fl1->mask),&cmask1);
  if(fl1==fl2)
    cmask2=cmask1;
  else {
    cmask2=dftw_malloc(fl2->fs->ny*(fs->nx/2+1)*sizeof(fcomplex));
    fs_map2alm(fl2->fs,1,0,&(fl2->mask),&cmask2);
  }
  maskprod=dftw_malloc(fl1->npix*sizeof(flouble));
  i_band=my_malloc(fl1->npix*sizeof(int));		    
  if(w->ncls>1) {
    cosarr=dftw_malloc(fl1->npix*sizeof(flouble));
    sinarr=dftw_malloc(fl1->npix*sizeof(flouble));
  }

#pragma omp parallel default(none)				\
  shared(fs,maskprod,cmask1,cmask2,w,i_band,cosarr,sinarr)
  {
    flouble dkx=2*M_PI/fs->lx;
    flouble dky=2*M_PI/fs->ly;
    int iy1,ix1;
    int *n_cells_thr=my_calloc(w->nells,sizeof(int));

#pragma omp for
    for(iy1=0;iy1<fs->ny;iy1++) {
      flouble ky;
      if(2*iy1<=fs->ny)
	ky=iy1*dky;
      else
	ky=-(fs->ny-iy1)*dky;
      for(ix1=0;ix1<fs->nx;ix1++) {
	flouble kx,kmod;
	int ix_here,index_here,index,ik;
	index=ix1+fs->nx*iy1;
	if(2*ix1<=fs->nx) {
	  kx=ix1*dkx;
	  ix_here=ix1;
	}
	else {
	  kx=-(fs->nx-ix1)*dkx;
	  ix_here=fs->nx-ix1;
	}
	index_here=ix_here+(fs->nx/2+1)*iy1;
	
	maskprod[index]=(creal(cmask1[index_here])*creal(cmask2[index_here])+
			 cimag(cmask1[index_here])*cimag(cmask2[index_here]));

	kmod=sqrt(kx*kx+ky*ky);
	ik=(int)(kmod*fs->i_dell);
	if(ik<w->nells) {
	  i_band[index]=ik;
	  n_cells_thr[ik]++;
	}
	else
	  i_band[index]=-1;

	if(w->ncls>1) {
	  flouble c,s;
	  if(kmod>0) {
	    c=kx/kmod;
	    s=ky/kmod;
	  }
	  else {
	    c=1.;
	    s=0.;
	  }
	  cosarr[index]=c*c-s*s;
	  sinarr[index]=2*s*c;
	}	
      }
    } //end omp for

#pragma omp critical
    {
      for(iy1=0;iy1<w->nells;iy1++)
	w->n_cells[iy1]+=n_cells_thr[iy1];
    } //end omp critical
    free(n_cells_thr);
  } //end omp parallel

#pragma omp parallel default(none)		\
  shared(fs,i_band,w,cosarr,sinarr,maskprod)
  {
    int iy1,ix1,ix2,iy2;
    flouble **coupling_matrix_thr=my_malloc(w->nells*w->ncls*sizeof(flouble *));
    for(iy1=0;iy1<w->nells*w->ncls;iy1++)
      coupling_matrix_thr[iy1]=my_calloc(w->nells*w->ncls,sizeof(flouble));

#pragma omp for
    for(iy1=0;iy1<fs->ny;iy1++) {
      for(ix1=0;ix1<fs->nx;ix1++) {
	int index1=ix1+fs->nx*iy1;
	int ik1=i_band[index1];
	if(ik1>=0) {
	  ik1*=w->ncls;
	  for(iy2=0;iy2<fs->ny;iy2++) {
	    for(ix2=0;ix2<fs->nx;ix2++) {
	      int index2=ix2+fs->nx*iy2;
	      int ik2=i_band[index2];
	      if(ik2>=0) {
		flouble cdiff=1,sdiff=1;
		int index;
		int iy=iy1-iy2;
		int ix=ix1-ix2;
		if(iy<0) iy+=fs->ny;
		if(ix<0) ix+=fs->nx;
		ik2*=w->ncls;
		index=ix+fs->nx*iy;

		if(w->ncls>1) {
		  cdiff=cosarr[index1]*cosarr[index2]+sinarr[index1]*sinarr[index2];
		  sdiff=sinarr[index1]*cosarr[index2]-cosarr[index1]*sinarr[index2];
		}

		if(w->ncls==1) {
		  coupling_matrix_thr[ik1+0][ik2+0]+=maskprod[index];
		}
		if(w->ncls==2) {
		  flouble fc=cdiff*maskprod[index];
		  flouble fs=sdiff*maskprod[index];
		  coupling_matrix_thr[ik1+0][ik2+0]+=fc; //CHECK SIGNS;
		  coupling_matrix_thr[ik1+0][ik2+1]-=fs;
		  coupling_matrix_thr[ik1+1][ik2+0]+=fs;
		  coupling_matrix_thr[ik1+1][ik2+1]+=fc;
		}
		if(w->ncls==4) {
		  flouble fcc=cdiff*cdiff*maskprod[index];
		  flouble fsc=sdiff*cdiff*maskprod[index];
		  flouble fss=sdiff*sdiff*maskprod[index];
		  coupling_matrix_thr[ik1+0][ik2+0]+=fcc; //CHECK SIGNS;
		  coupling_matrix_thr[ik1+0][ik2+1]-=fsc;
		  coupling_matrix_thr[ik1+0][ik2+2]-=fsc;
		  coupling_matrix_thr[ik1+0][ik2+3]+=fss;
		  coupling_matrix_thr[ik1+1][ik2+0]+=fsc;
		  coupling_matrix_thr[ik1+1][ik2+1]+=fcc;
		  coupling_matrix_thr[ik1+1][ik2+2]-=fss;
		  coupling_matrix_thr[ik1+1][ik2+3]-=fsc;
		  coupling_matrix_thr[ik1+2][ik2+0]+=fsc;
		  coupling_matrix_thr[ik1+2][ik2+1]-=fss;
		  coupling_matrix_thr[ik1+2][ik2+2]+=fcc;
		  coupling_matrix_thr[ik1+2][ik2+3]-=fsc;
		  coupling_matrix_thr[ik1+3][ik2+0]+=fss;
		  coupling_matrix_thr[ik1+3][ik2+1]+=fsc;
		  coupling_matrix_thr[ik1+3][ik2+2]+=fsc;
		  coupling_matrix_thr[ik1+3][ik2+3]+=fcc;
		}
	      }
	    }
	  }
	}
      }
    } //end omp for

#pragma omp critical
    {
      for(iy1=0;iy1<w->ncls*w->nells;iy1++) {
	for(iy2=0;iy2<w->ncls*w->nells;iy2++)
	  w->coupling_matrix_unbinned[iy1][iy2]+=coupling_matrix_thr[iy1][iy2];
      }
    } //end omp critical

    for(iy1=0;iy1<w->nells*w->ncls;iy1++)
      free(coupling_matrix_thr[iy1]);
    free(coupling_matrix_thr);
  } //end omp parallel

#pragma omp parallel default(none) \
  shared(w,beam_prod,fs)
  {
    int il1;
    flouble fac_norm=4*M_PI*M_PI/(fs->lx*fs->lx*fs->ly*fs->ly);

#pragma omp for
    for(il1=0;il1<w->nells;il1++) {
      int icl1;
      flouble norm;
      if(w->n_cells[il1]>0)
	norm=fac_norm/w->n_cells[il1];
      else
	norm=0;
      for(icl1=0;icl1<w->ncls;icl1++) {
	int il2;
	for(il2=0;il2<w->nells;il2++) {
	  int icl2;
	  for(icl2=0;icl2<w->ncls;icl2++)
	    w->coupling_matrix_unbinned[w->ncls*il1+icl1][w->ncls*il2+icl2]*=beam_prod[il2]*norm;
	}
      }
    } //end omp for
  } //end omp parallel

  bin_coupling_matrix(w);

  free(beam_prod);
  dftw_free(cmask1);
  if(fl1!=fl2)
    dftw_free(cmask2);
  dftw_free(maskprod);
  free(i_band);
  if(w->ncls>1) {
    dftw_free(cosarr);
    dftw_free(sinarr);
  }

  return w;
}

void nmt_compute_deprojection_bias_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
					int nl_prop,flouble *l_prop,flouble **cl_proposal,
					flouble **cl_bias)
{
  //Placeholder
  return;
}

void nmt_couple_cl_l_flat(nmt_workspace_flat *w,int nl,flouble *larr,flouble **cl_in,flouble **cl_out)
{
  int ii;
  flouble **cell_in=my_malloc(w->ncls*sizeof(flouble *));
  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  for(ii=0;ii<w->ncls;ii++) {
    int il;
    nmt_k_function *fcl=nmt_k_function_alloc(nl,larr,cl_in[ii],cl_in[ii][0],0.,0);
    cell_in[ii]=my_malloc(w->nells*sizeof(flouble));
    for(il=0;il<w->nells;il++)
      cell_in[ii][il]=nmt_k_function_eval(fcl,w->l_arr[il],intacc);
    nmt_k_function_free(fcl);
  }
  gsl_interp_accel_free(intacc);

  int icl1;
  for(icl1=0;icl1<w->ncls;icl1++) {
    for(ii=0;ii<w->fs->n_ell;ii++) {
      int ir;
      flouble clo=0;
      for(ir=0;ir<w->nl_rebin;ir++) {
	int icl2;
	int i1=ii*w->nl_rebin+ir;
	if(i1<w->nells) {
	  int ind1=i1*w->ncls+icl1;
	  for(icl2=0;icl2<w->ncls;icl2++) {
	    int i2;
	    for(i2=0;i2<w->nells;i2++) {
	      int ind2=i2*w->ncls+icl2;
	      flouble cli=cell_in[icl2][i2];
	      clo+=w->coupling_matrix_unbinned[ind1][ind2]*cli;
	    }
	  }
	}
      }
      cl_out[icl1][ii]=clo/w->nl_rebin;
    }
  }

  for(ii=0;ii<w->ncls;ii++)
    free(cell_in[ii]);
  free(cell_in);
}

void nmt_decouple_cl_l_flat(nmt_workspace_flat *w,flouble **cl_in,flouble **cl_noise_in,
			    flouble **cl_bias,flouble **cl_out)
{
  int icl,ib2;
  gsl_vector *dl_map_bad_b=gsl_vector_alloc(w->ncls*w->bin->n_bands);
  gsl_vector *dl_map_good_b=gsl_vector_alloc(w->ncls*w->bin->n_bands);

  //Bin coupled power spectrum
  for(icl=0;icl<w->ncls;icl++) {
    int il;
    flouble *sum=my_calloc(w->bin->n_bands,sizeof(flouble));
    int *n_per_band=my_calloc(w->bin->n_bands,sizeof(int));
    for(il=0;il<w->fs->n_ell;il++) {
      int ir;
      flouble dcl=cl_in[icl][il]-cl_noise_in[icl][il]-cl_bias[icl][il];
      for(ir=0;ir<w->nl_rebin;ir++) {
	int index=il*w->nl_rebin+ir;
	if(index<w->nells) {
	  int ind=w->i_band[il*w->nl_rebin+ir];
	  if(ind>=0) {
	    sum[ind]+=dcl;
	    n_per_band[ind]++;
	  }
	}
      }
    }
    for(il=0;il<w->bin->n_bands;il++) {
      if(n_per_band[il]>0)
	sum[il]/=n_per_band[il];
      gsl_vector_set(dl_map_bad_b,w->ncls*il+icl,sum[il]);
    }
    free(sum);
    free(n_per_band);
  }

  gsl_linalg_LU_solve(w->coupling_matrix_binned,w->coupling_matrix_perm,dl_map_bad_b,dl_map_good_b);
  for(icl=0;icl<w->ncls;icl++) {
    for(ib2=0;ib2<w->bin->n_bands;ib2++)
      cl_out[icl][ib2]=gsl_vector_get(dl_map_good_b,w->ncls*ib2+icl);
  }

  gsl_vector_free(dl_map_bad_b);
  gsl_vector_free(dl_map_good_b);
}

void nmt_compute_coupled_cell_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,flouble *larr,flouble **cl_out)
{
  int i;
  fs_alm2cl(fl1->fs,fl1->alms,fl2->alms,fl1->pol,fl2->pol,cl_out);
  for(i=0;i<fl1->fs->n_ell;i++)
    larr[i]=fl1->fs->ell_min[i]+(i+0.5)*fl1->fs->dell;
}

#define NMT_COUPLING_FLAT_LINT 201
#define NMT_COUPLING_FLAT_ANG  202
#define NMT_COUPLING_FLAT_FULL 203
nmt_workspace_flat *nmt_compute_coupling_matrix_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
						     nmt_binning_scheme_flat *bin,int nl_rebin,int method_flag)
{
  if(check_flatsky_infos(fl1->fs,fl2->fs))
    report_error(1,"Can only correlate fields defined on the same pixels!\n");

  if(method_flag==NMT_COUPLING_FLAT_LINT)
    return nmt_compute_coupling_matrix_flat_a(fl1,fl2,bin,nl_rebin);
  else if(method_flag==NMT_COUPLING_FLAT_ANG)
    return nmt_compute_coupling_matrix_flat_b(fl1,fl2,bin,nl_rebin);
  else if(method_flag==NMT_COUPLING_FLAT_FULL)
    return nmt_compute_coupling_matrix_flat_c(fl1,fl2,bin,nl_rebin);
  else {
    report_error(1,"Unknown method %d\n",method_flag);
    return NULL;
  }
}

nmt_workspace_flat *nmt_compute_power_spectra_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
						   nmt_binning_scheme_flat *bin,int nl_rebin,int method_flag,
						   nmt_workspace_flat *w0,flouble **cl_noise,
						   int nl_prop,flouble *l_prop,flouble **cl_prop,
						   flouble **cl_out)
{
  int ii;
  flouble **cl_bias,**cl_data,*larr;
  nmt_workspace_flat *w;

  if(w0==NULL)
    w=nmt_compute_coupling_matrix_flat(fl1,fl2,bin,nl_rebin,method_flag);
  else
    w=w0;

  larr=my_malloc(w->fs->n_ell*sizeof(flouble));
  cl_bias=my_malloc(w->ncls*sizeof(flouble *));
  cl_data=my_malloc(w->ncls*sizeof(flouble *));
  for(ii=0;ii<w->ncls;ii++) {
    cl_bias[ii]=my_calloc(w->fs->n_ell,sizeof(flouble));
    cl_data[ii]=my_calloc(w->fs->n_ell,sizeof(flouble));
  }
  nmt_compute_coupled_cell_flat(fl1,fl2,larr,cl_data);
  nmt_compute_deprojection_bias_flat(fl1,fl2,nl_prop,l_prop,cl_prop,cl_bias);
  nmt_decouple_cl_l_flat(w,cl_data,cl_noise,cl_bias,cl_out);
  for(ii=0;ii<w->ncls;ii++) {
    free(cl_bias[ii]);
    free(cl_data[ii]);
  }
  free(cl_bias);
  free(cl_data);

  return w;
}
