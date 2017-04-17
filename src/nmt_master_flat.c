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
  free(w);
}

static nmt_workspace_flat *nmt_workspace_flat_new(int nl_rebin,int ncls,nmt_flatsky_info *fs,
						  nmt_binning_scheme_flat *bin)
{
  int ii;
  nmt_workspace_flat *w=my_malloc(sizeof(nmt_workspace_flat));
  w->ncls=ncls;
  w->nl_rebin=nl_rebin;

  w->bin=nmt_bins_flat_create(bin->n_bands,bin->ell_0_list,bin->ell_f_list);

  w->fs=nmt_flatsky_info_alloc(fs->nx,fs->ny,fs->lx,fs->ly);

  w->nells=w->fs->n_ell*w->nl_rebin;
  w->pcl_masks=my_malloc(w->nells*sizeof(flouble));

  w->l_arr=my_malloc(w->nells*sizeof(flouble));
  w->i_band=my_malloc(w->nells*sizeof(int));
  for(ii=0;ii<w->fs->n_ell;ii++) {
    int jj,ib=0;
    for(jj=0;jj<w->nl_rebin;jj++) {
      flouble l_here=w->fs->ell_min[ii]+jj*w->fs->dell/w->nl_rebin;
      ib=nmt_bins_flat_search_fast(w->bin,l_here,ib);
      w->l_arr[ii*w->nl_rebin+jj] =l_here;
      w->i_band[ii*w->nl_rebin+jj]=ib;
    }
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
  flouble dell=w->fs->dell/w->nl_rebin;
  flouble **coupling_b=my_malloc(w->ncls*w->bin->n_bands*sizeof(flouble *));
  for(ib2=0;ib2<w->ncls*w->bin->n_bands;ib2++)
    coupling_b[ib2]=my_calloc(w->ncls*w->bin->n_bands,sizeof(flouble *));
  
  for(icl_a=0;icl_a<w->ncls;icl_a++) {
    for(icl_b=0;icl_b<w->ncls;icl_b++) {
      for(il2=0;il2<w->nells;il2++) {
	ib2=w->i_band[il2];
	if(ib2>=0) {
	  for(il3=0;il3<w->nells;il3++) {
	    ib3=w->i_band[il3];
	    if(ib3>=0)
	      coupling_b[w->ncls*ib2+icl_a][w->ncls*ib3+icl_b]+=w->coupling_matrix_unbinned[w->ncls*il2+icl_a][w->ncls*il3+icl_b];
	  }
	}
      }
      for(ib2=0;ib2<w->bin->n_bands;ib2++) {
	flouble norm=dell/(w->bin->ell_f_list[ib2]-w->bin->ell_0_list[ib2]);
	for(ib3=0;ib3<w->bin->n_bands;ib3++) {
	  flouble c=coupling_b[w->ncls*ib2+icl_a][w->ncls*ib3+icl_b]*norm;
	  gsl_matrix_set(w->coupling_matrix_binned,w->ncls*ib2+icl_a,w->ncls*ib3+icl_b,c);
	}
      }
    }
  }

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

  gsl_linalg_LU_decomp(w->coupling_matrix_binned,w->coupling_matrix_perm,&sig);

  for(ib2=0;ib2<w->ncls*w->bin->n_bands;ib2++)
    free(coupling_b[ib2]);
  free(coupling_b);
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

nmt_workspace_flat *nmt_compute_coupling_matrix_flat_b(nmt_field_flat *fl1,nmt_field_flat *fl2,
						       nmt_binning_scheme_flat *bin,int nl_rebin)
{
  int ii,n_cl=fl1->nmaps*fl2->nmaps;
  nmt_workspace_flat *w=nmt_workspace_flat_new(nl_rebin,n_cl,fl1->fs,bin);

  //TODO: check that fields are compatible
  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  flouble *beam_prod=my_malloc(w->nells*sizeof(flouble));
  for(ii=0;ii<w->fs->n_ell;ii++) {
    int ir;
    for(ir=0;ir<w->nl_rebin;ir++) {
      int ibin=ii*w->nl_rebin+ir;
      flouble l_here=w->l_arr[ibin];
      flouble b1=nmt_k_function_eval(fl1->beam,l_here,intacc);
      flouble b2=nmt_k_function_eval(fl2->beam,l_here,intacc);
      beam_prod[ibin]=b1*b2;
    }
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

nmt_workspace_flat *nmt_compute_coupling_matrix_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
						     nmt_binning_scheme_flat *bin,int nl_rebin)
{
  int l2;
  int n_cl=fl1->nmaps*fl2->nmaps;
  flouble *beam_prod;
  flouble *pcl_masks_raw=my_malloc(fl1->fs->n_ell*sizeof(flouble));
  nmt_workspace_flat *w=nmt_workspace_flat_new(nl_rebin,n_cl,fl1->fs,bin);
  flouble dell_here=w->fs->dell/w->nl_rebin;
  //TODO: check that fields are compatible

  fs_anafast(fl1->fs,&(fl1->mask),&(fl2->mask),0,0,&(pcl_masks_raw));

  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  beam_prod=my_malloc(w->nells*sizeof(flouble));
  for(l2=0;l2<w->fs->n_ell;l2++) {
    int ir;
    for(ir=0;ir<w->nl_rebin;ir++) {
      int ibin=l2*w->nl_rebin+ir;
      flouble l_here=w->l_arr[ibin];
      flouble b1=nmt_k_function_eval(fl1->beam,l_here,intacc);
      flouble b2=nmt_k_function_eval(fl2->beam,l_here,intacc);
      beam_prod[ibin]=b1*b2;
      w->pcl_masks[ibin]=dell_here*pcl_masks_raw[l2]*l_here/(4*M_PI*M_PI);
    }
  }
  gsl_interp_accel_free(intacc);

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
  free(pcl_masks_raw);

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
  flouble dell=w->fs->dell/w->nl_rebin;
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
      cl_out[icl1][ii]=clo*dell/w->nl_rebin;
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
  flouble dell=w->fs->dell/w->nl_rebin;
  gsl_vector *dl_map_bad_b=gsl_vector_alloc(w->ncls*w->bin->n_bands);
  gsl_vector *dl_map_good_b=gsl_vector_alloc(w->ncls*w->bin->n_bands);

  //Bin coupled power spectrum
  for(icl=0;icl<w->ncls;icl++) {
    int il;
    flouble *sum=my_calloc(w->bin->n_bands,sizeof(flouble));
    for(il=0;il<w->fs->n_ell;il++) {
      int ir;
      flouble dcl=cl_in[icl][il]-cl_noise_in[icl][il]-cl_bias[icl][il];
      for(ir=0;ir<w->nl_rebin;ir++) {
	int ind=w->i_band[il*w->nl_rebin+ir];
	if(ind>=0)
	  sum[ind]+=dcl;
      }
    }
    for(il=0;il<w->bin->n_bands;il++) {
      flouble norm=dell/(w->bin->ell_f_list[il]-w->bin->ell_0_list[il]);
      gsl_vector_set(dl_map_bad_b,w->ncls*il+icl,sum[il]*norm);
    }
    free(sum);
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

nmt_workspace_flat *nmt_compute_power_spectra_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
						   nmt_binning_scheme_flat *bin,int nl_rebin,
						   nmt_workspace_flat *w0,flouble **cl_noise,
						   int nl_prop,flouble *l_prop,flouble **cl_prop,
						   flouble **cl_out)
{
  int ii;
  flouble **cl_bias,**cl_data,*larr;
  nmt_workspace_flat *w;

  if(w0==NULL)
    w=nmt_compute_coupling_matrix_flat(fl1,fl2,bin,nl_rebin);
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
