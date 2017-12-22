#include "utils.h"

void nmt_workspace_flat_free(nmt_workspace_flat *w)
{
  int ii;
  gsl_permutation_free(w->coupling_matrix_perm);
  gsl_matrix_free(w->coupling_matrix_binned);
  for(ii=0;ii<w->ncls*w->bin->n_bands;ii++)
    free(w->coupling_matrix_unbinned[ii]);
  free(w->coupling_matrix_unbinned);
  free(w->n_cells);
  nmt_bins_flat_free(w->bin);
  nmt_flatsky_info_free(w->fs);
  free(w);
}

static nmt_workspace_flat *nmt_workspace_flat_new(int ncls,nmt_flatsky_info *fs,
						  nmt_binning_scheme_flat *bin,
						  flouble lmn_x,flouble lmx_x,flouble lmn_y,flouble lmx_y)
{
  int ii,ib=0;
  nmt_workspace_flat *w=my_malloc(sizeof(nmt_workspace_flat));
  w->ncls=ncls;

  w->ellcut_x[0]=lmn_x;
  w->ellcut_x[1]=lmx_x;
  w->ellcut_y[0]=lmn_y;
  w->ellcut_y[1]=lmx_y;

  w->bin=nmt_bins_flat_create(bin->n_bands,bin->ell_0_list,bin->ell_f_list);
  w->lmax=w->bin->ell_f_list[w->bin->n_bands-1];

  w->fs=nmt_flatsky_info_alloc(fs->nx,fs->ny,fs->lx,fs->ly);

  w->n_cells=my_calloc(w->bin->n_bands,sizeof(int));

  w->coupling_matrix_unbinned=my_malloc(w->ncls*w->bin->n_bands*sizeof(flouble *));
  for(ii=0;ii<w->ncls*w->bin->n_bands;ii++)
    w->coupling_matrix_unbinned[ii]=my_calloc(w->ncls*w->bin->n_bands,sizeof(flouble));

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

  my_fread(w->ellcut_x,sizeof(flouble),2,fi);
  my_fread(w->ellcut_y,sizeof(flouble),2,fi);

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
  w->lmax=w->bin->ell_f_list[w->bin->n_bands-1];

  w->n_cells=my_malloc(w->bin->n_bands*sizeof(int));
  my_fread(w->n_cells,sizeof(int),w->bin->n_bands,fi);

  w->coupling_matrix_unbinned=my_malloc(w->ncls*w->bin->n_bands*sizeof(flouble *));
  for(ii=0;ii<w->ncls*w->bin->n_bands;ii++) {
    w->coupling_matrix_unbinned[ii]=my_malloc(w->ncls*w->bin->n_bands*sizeof(flouble));
    my_fread(w->coupling_matrix_unbinned[ii],sizeof(flouble),w->ncls*w->bin->n_bands,fi);
  }

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
  my_fwrite(w->ellcut_x,sizeof(flouble),2,fo);
  my_fwrite(w->ellcut_y,sizeof(flouble),2,fo);

  my_fwrite(&(w->fs->nx),sizeof(int),1,fo);
  my_fwrite(&(w->fs->ny),sizeof(int),1,fo);
  my_fwrite(&(w->fs->lx),sizeof(flouble),1,fo);
  my_fwrite(&(w->fs->ly),sizeof(flouble),1,fo);

  my_fwrite(&(w->bin->n_bands),sizeof(int),1,fo);
  my_fwrite(w->bin->ell_0_list,sizeof(flouble),w->bin->n_bands,fo);
  my_fwrite(w->bin->ell_f_list,sizeof(flouble),w->bin->n_bands,fo);

  my_fwrite(w->n_cells,sizeof(int),w->bin->n_bands,fo);

  for(ii=0;ii<w->ncls*w->bin->n_bands;ii++)
    my_fwrite(w->coupling_matrix_unbinned[ii],sizeof(flouble),w->ncls*w->bin->n_bands,fo);

  gsl_matrix_fwrite(fo,w->coupling_matrix_binned);
  gsl_permutation_fwrite(fo,w->coupling_matrix_perm);
  
  fclose(fo);
}

static void bin_coupling_matrix(nmt_workspace_flat *w)
{
  int icl_a,icl_b,ib2,ib3,sig;
  
  for(icl_a=0;icl_a<w->ncls;icl_a++) {
    for(icl_b=0;icl_b<w->ncls;icl_b++) {
      for(ib2=0;ib2<w->bin->n_bands;ib2++) {
	for(ib3=0;ib3<w->bin->n_bands;ib3++) {
	  gsl_matrix_set(w->coupling_matrix_binned,w->ncls*ib2+icl_a,w->ncls*ib3+icl_b,
			 w->coupling_matrix_unbinned[w->ncls*ib2+icl_a][w->ncls*ib3+icl_b]);
	}
      }
    }
  }

  gsl_linalg_LU_decomp(w->coupling_matrix_binned,w->coupling_matrix_perm,&sig);
}

static int check_flatsky_infos(nmt_flatsky_info *fs1,nmt_flatsky_info *fs2)
{
  if(fs1->nx!=fs2->nx) return 1;
  if(fs1->ny!=fs2->ny) return 1;
  if(fs1->lx!=fs2->lx) return 1;
  if(fs1->ly!=fs2->ly) return 1;
  return 0;
}

static nmt_workspace_flat *nmt_compute_coupling_matrix_flat_c(nmt_field_flat *fl1,nmt_field_flat *fl2,
							      nmt_binning_scheme_flat *bin,
							      flouble lmn_x,flouble lmx_x,
							      flouble lmn_y,flouble lmx_y)
{
  int ii;
  nmt_workspace_flat *w=nmt_workspace_flat_new(fl1->nmaps*fl2->nmaps,fl1->fs,bin,
					       lmn_x,lmx_x,lmn_y,lmx_y);
  nmt_flatsky_info *fs=fl1->fs;
  
  fcomplex *cmask1,*cmask2;
  flouble *maskprod,*cosarr,*sinarr,*kmodarr;
  int *i_band;
  cmask1=dftw_malloc(fs->ny*(fs->nx/2+1)*sizeof(fcomplex));
  fs_map2alm(fl1->fs,1,0,&(fl1->mask),&cmask1);
  if(fl1==fl2)
    cmask2=cmask1;
  else {
    cmask2=dftw_malloc(fs->ny*(fs->nx/2+1)*sizeof(fcomplex));
    fs_map2alm(fl2->fs,1,0,&(fl2->mask),&cmask2);
  }
  maskprod=dftw_malloc(fl1->npix*sizeof(flouble));
  i_band=my_malloc(fl1->npix*sizeof(int));		    
  if(w->ncls>1) {
    kmodarr=dftw_malloc(fl1->npix*sizeof(flouble));
    cosarr=dftw_malloc(fl1->npix*sizeof(flouble));
    sinarr=dftw_malloc(fl1->npix*sizeof(flouble));
  }

  int *x_out_range,*y_out_range;
  x_out_range=my_calloc(fs->nx,sizeof(int));
  y_out_range=my_calloc(fs->ny,sizeof(int));
  for(ii=0;ii<fs->nx;ii++) {
    flouble k;
    if(2*ii<=fs->nx) k=ii*2*M_PI/fs->lx;
    else k=-(fs->nx-ii)*2*M_PI/fs->lx;
    if((k<=w->ellcut_x[1]) && (k>=w->ellcut_x[0]))
      x_out_range[ii]=1;
  }
  for(ii=0;ii<fs->ny;ii++) {
    flouble k;
    if(2*ii<=fs->ny) k=ii*2*M_PI/fs->ly;
    else k=-(fs->ny-ii)*2*M_PI/fs->ly;
    if((k<=w->ellcut_y[1]) && (k>=w->ellcut_y[0]))
      y_out_range[ii]=1;
  }

#pragma omp parallel default(none)					\
  shared(fs,maskprod,cmask1,cmask2,w,i_band,cosarr,sinarr,kmodarr)	\
  shared(x_out_range,y_out_range)
  {
    flouble dkx=2*M_PI/fs->lx;
    flouble dky=2*M_PI/fs->ly;
    int iy1,ix1;
    int *n_cells_thr=my_calloc(w->bin->n_bands,sizeof(int));

#pragma omp for
    for(iy1=0;iy1<fs->ny;iy1++) {
      flouble ky;
      int ik=0;
      if(2*iy1<=fs->ny)
	ky=iy1*dky;
      else
	ky=-(fs->ny-iy1)*dky;
      for(ix1=0;ix1<fs->nx;ix1++) {
	flouble kx,kmod;
	int ix_here,index_here,index;
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

	if(y_out_range[iy1] || x_out_range[ix1])
	  i_band[index]=-1;
	else {
	  kmod=sqrt(kx*kx+ky*ky);
	  ik=nmt_bins_flat_search_fast(w->bin,kmod,ik);
	  if(ik>=0) {
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
	    kmodarr[index]=kmod;
	  }	
	}
      }
    } //end omp for

#pragma omp critical
    {
      for(iy1=0;iy1<w->bin->n_bands;iy1++)
	w->n_cells[iy1]+=n_cells_thr[iy1];
    } //end omp critical
    free(n_cells_thr);
  } //end omp parallel
  free(x_out_range);
  free(y_out_range);

#pragma omp parallel default(none)				\
  shared(fs,i_band,w,cosarr,sinarr,kmodarr,maskprod,fl1,fl2)
  {
    int iy1,ix1,ix2,iy2;
    int pe1=fl1->pure_e,pe2=fl2->pure_e,pb1=fl1->pure_b,pb2=fl2->pure_b;
    int pure_any=pe1 || pb1 || pe2 || pb2;
    flouble **coupling_matrix_thr=my_malloc(w->bin->n_bands*w->ncls*sizeof(flouble *));
    for(iy1=0;iy1<w->bin->n_bands*w->ncls;iy1++)
      coupling_matrix_thr[iy1]=my_calloc(w->bin->n_bands*w->ncls,sizeof(flouble));

#pragma omp for
    for(iy1=0;iy1<fs->ny;iy1++) {
      for(ix1=0;ix1<fs->nx;ix1++) {
	int index1=ix1+fs->nx*iy1;
	int ik1=i_band[index1];
	if(ik1>=0) {
	  flouble inv_k1=0;
	  ik1*=w->ncls;
	  if((index1>0) && (w->ncls>1))
	    inv_k1=1./kmodarr[index1];
	  for(iy2=0;iy2<fs->ny;iy2++) {
	    for(ix2=0;ix2<fs->nx;ix2++) {
	      int index2=ix2+fs->nx*iy2;
	      int ik2=i_band[index2];
	      if(ik2>=0) {
		flouble cdiff=1,sdiff=0,kr=1,mp;
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
		  if((index1==0) && (index2==0))
		    kr=1;
		  else
		    kr=kmodarr[index2]*inv_k1;
		  kr*=kr;
		}
		mp=maskprod[index];

		if(w->ncls==1) {
		  coupling_matrix_thr[ik1+0][ik2+0]+=mp;
		}
		if(w->ncls==2) {
		  flouble fc[2],fs[2];
		  fc[0]=cdiff*mp;
		  fs[0]=sdiff*mp;
		  if(pure_any) {
		    fc[1]=kr*mp; fs[1]=0;
		  }
		  coupling_matrix_thr[ik1+0][ik2+0]+=fc[pe1+pe2]; //TE,TE
		  coupling_matrix_thr[ik1+0][ik2+1]-=fs[pe1+pe2]; //TE,TB
		  coupling_matrix_thr[ik1+1][ik2+0]+=fs[pb1+pb2]; //TB,TE
		  coupling_matrix_thr[ik1+1][ik2+1]+=fc[pb1+pb2]; //TB,TB
		}
		if(w->ncls==4) {
		  flouble fc[2],fs[2];
		  fc[0]=cdiff; fs[0]=sdiff;
		  if(pure_any) {
		    fc[1]=kr; fs[1]=0;
		  }
		  coupling_matrix_thr[ik1+0][ik2+0]+=fc[pe1]*fc[pe2]*mp; //EE,EE
		  coupling_matrix_thr[ik1+0][ik2+1]-=fc[pe1]*fs[pe2]*mp; //EE,EB
		  coupling_matrix_thr[ik1+0][ik2+2]-=fs[pe1]*fc[pe2]*mp; //EE,BE
		  coupling_matrix_thr[ik1+0][ik2+3]+=fs[pe1]*fs[pe2]*mp; //EE,BB
		  coupling_matrix_thr[ik1+1][ik2+0]+=fc[pe1]*fs[pb2]*mp; //EB,EE
		  coupling_matrix_thr[ik1+1][ik2+1]+=fc[pe1]*fc[pb2]*mp; //EB,EB
		  coupling_matrix_thr[ik1+1][ik2+2]-=fs[pe1]*fs[pb2]*mp; //EB,BE
		  coupling_matrix_thr[ik1+1][ik2+3]-=fs[pe1]*fc[pb2]*mp; //EB,BB
		  coupling_matrix_thr[ik1+2][ik2+0]+=fs[pb1]*fc[pe2]*mp; //BE,EE
		  coupling_matrix_thr[ik1+2][ik2+1]-=fs[pb1]*fs[pe2]*mp; //BE,EB
		  coupling_matrix_thr[ik1+2][ik2+2]+=fc[pb1]*fc[pe2]*mp; //BE,BE
		  coupling_matrix_thr[ik1+2][ik2+3]-=fc[pb1]*fs[pe2]*mp; //BE,BB
		  coupling_matrix_thr[ik1+3][ik2+0]+=fs[pb1]*fs[pb2]*mp; //BB,EE
		  coupling_matrix_thr[ik1+3][ik2+1]+=fs[pb1]*fc[pb2]*mp; //BB,EB
		  coupling_matrix_thr[ik1+3][ik2+2]+=fc[pb1]*fs[pb2]*mp; //BB,BE
		  coupling_matrix_thr[ik1+3][ik2+3]+=fc[pb1]*fc[pb2]*mp; //BB,BB
		}
	      }
	    }
	  }
	}
      }
    } //end omp for

#pragma omp critical
    {
      for(iy1=0;iy1<w->ncls*w->bin->n_bands;iy1++) {
	for(iy2=0;iy2<w->ncls*w->bin->n_bands;iy2++)
	  w->coupling_matrix_unbinned[iy1][iy2]+=coupling_matrix_thr[iy1][iy2];
      }
    } //end omp critical

    for(iy1=0;iy1<w->bin->n_bands*w->ncls;iy1++)
      free(coupling_matrix_thr[iy1]);
    free(coupling_matrix_thr);
  } //end omp parallel

#pragma omp parallel default(none) \
  shared(w,fs)
  {
    int il1;
    flouble fac_norm=4*M_PI*M_PI/(fs->lx*fs->lx*fs->ly*fs->ly);

#pragma omp for
    for(il1=0;il1<w->bin->n_bands;il1++) {
      int icl1;
      flouble norm;
      if(w->n_cells[il1]>0)
	norm=fac_norm/w->n_cells[il1];
      else
	norm=0;
      for(icl1=0;icl1<w->ncls;icl1++) {
	int il2;
	for(il2=0;il2<w->bin->n_bands;il2++) {
	  int icl2;
	  for(icl2=0;icl2<w->ncls;icl2++)
	    w->coupling_matrix_unbinned[w->ncls*il1+icl1][w->ncls*il2+icl2]*=norm;
	}
      }
    } //end omp for
  } //end omp parallel

  bin_coupling_matrix(w);

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
					nmt_binning_scheme_flat *bin,
					flouble lmn_x,flouble lmx_x,flouble lmn_y,flouble lmx_y,
					int nl_prop,flouble *l_prop,flouble **cl_proposal,
					flouble **cl_bias)
{
  //Placeholder
  int ii;
  long ip;
  int nspec=fl1->nmaps*fl2->nmaps;
  flouble **cl_dum=my_malloc(nspec*sizeof(flouble *));
  nmt_k_function **cl_proposal_f=my_malloc(nspec*sizeof(nmt_k_function *));
  for(ii=0;ii<nspec;ii++) {
    cl_dum[ii]=my_calloc(bin->n_bands,sizeof(flouble));
    cl_proposal_f[ii]=nmt_k_function_alloc(nl_prop,l_prop,cl_proposal[ii],cl_proposal[ii][0],0,0);
    for(ip=0;ip<bin->n_bands;ip++)
      cl_bias[ii][ip]=0;
  }

  //TODO: some terms (e.g. C^ab*SHT[w*g^j]) could be precomputed
  //TODO: if fl1=fl2 F2=F3
  //Allocate dummy maps and alms
  flouble **map_1_dum=my_malloc(fl1->nmaps*sizeof(flouble *));
  fcomplex **alm_1_dum=my_malloc(fl1->nmaps*sizeof(fcomplex *));
  for(ii=0;ii<fl1->nmaps;ii++) {
    map_1_dum[ii]=dftw_malloc(fl1->npix*sizeof(flouble));
    alm_1_dum[ii]=dftw_malloc(fl1->fs->ny*(fl1->fs->nx/2+1)*sizeof(fcomplex));
  }
  flouble **map_2_dum=my_malloc(fl2->nmaps*sizeof(flouble *));
  fcomplex **alm_2_dum=my_malloc(fl2->nmaps*sizeof(fcomplex *));
  for(ii=0;ii<fl2->nmaps;ii++) {
    map_2_dum[ii]=dftw_malloc(fl2->npix*sizeof(flouble));
    alm_2_dum[ii]=dftw_malloc(fl2->fs->ny*(fl2->fs->nx/2+1)*sizeof(fcomplex));
  }

  if(fl2->ntemp>0) {
    int iti;
    for(iti=0;iti<fl2->ntemp;iti++) {
      int itj;
      for(itj=0;itj<fl2->ntemp;itj++) {
	int im1,im2;
	double nij=gsl_matrix_get(fl2->matrix_M,iti,itj);
	//w*g^j
	for(im2=0;im2<fl2->nmaps;im2++)
	  fs_map_product(fl2->fs,fl2->temp[itj][im2],fl2->mask,map_2_dum[im2]);
	//DFT[w*g^j]
	fs_map2alm(fl2->fs,1,2*fl2->pol,map_2_dum,alm_2_dum);
	//C^ab*DFT[w*g^j]
	for(im1=0;im1<fl1->nmaps;im1++) {
	  fs_zero_alm(fl1->fs,alm_1_dum[im1]);
	  for(im2=0;im2<fl2->nmaps;im2++)
	    fs_alter_alm(fl2->fs,-1.,alm_2_dum[im2],alm_1_dum[im1],cl_proposal_f[im1*fl2->nmaps+im2],1);
	}
	//DFT^-1[C^ab*DFT[w*g^j]]
	fs_alm2map(fl1->fs,1,2*fl1->pol,map_1_dum,alm_1_dum);
	//v*DFT^-1[C^ab*DFT[w*g^j]]
	for(im1=0;im1<fl1->nmaps;im1++)
	  fs_map_product(fl1->fs,map_1_dum[im1],fl1->mask,map_1_dum[im1]);
	//DFT[v*DFT^-1[C^ab*DFT[w*g^j]]]
	fs_map2alm(fl1->fs,1,2*fl1->pol,map_1_dum,alm_1_dum);
	//Sum_m(DFT[v*DFT^-1[C^ab*DFT[w*g^j]]]*g^i*)/(2l+1)
	fs_alm2cl(fl1->fs,bin,alm_1_dum,fl2->a_temp[iti],fl1->pol,fl2->pol,cl_dum,lmn_x,lmx_x,lmn_y,lmx_y);
	for(im1=0;im1<nspec;im1++) {
	  for(ip=0;ip<bin->n_bands;ip++)
	    cl_bias[im1][ip]-=cl_dum[im1][ip]*nij;
	}
      }
    }
  }

  if(fl1->ntemp>0) {
    int iti;
    for(iti=0;iti<fl1->ntemp;iti++) {
      int itj;
      for(itj=0;itj<fl1->ntemp;itj++) {
	int im1,im2;
	double mij=gsl_matrix_get(fl1->matrix_M,iti,itj);
	//v*f^j
	for(im1=0;im1<fl1->nmaps;im1++)
	  fs_map_product(fl1->fs,fl1->temp[itj][im1],fl1->mask,map_1_dum[im1]);
	//DFT[v*f^j]
	fs_map2alm(fl1->fs,1,2*fl1->pol,map_1_dum,alm_1_dum);
	//C^abT*DFT[v*f^j]
	for(im2=0;im2<fl2->nmaps;im2++) {
	  fs_zero_alm(fl2->fs,alm_2_dum[im2]);
	  for(im1=0;im1<fl1->nmaps;im1++)
	    fs_alter_alm(fl1->fs,-1.,alm_1_dum[im1],alm_2_dum[im2],cl_proposal_f[im1*fl2->nmaps+im2],1);
	}
	//DFT^-1[C^abT*DFT[v*f^j]]
	fs_alm2map(fl2->fs,1,2*fl2->pol,map_2_dum,alm_2_dum);
	//w*DFT^-1[C^abT*DFT[v*f^j]]
	for(im2=0;im2<fl2->nmaps;im2++)
	  fs_map_product(fl2->fs,map_2_dum[im2],fl2->mask,map_2_dum[im2]);
	//DFT[w*DFT^-1[C^abT*DFT[v*f^j]]]
	fs_map2alm(fl2->fs,1,2*fl2->pol,map_2_dum,alm_2_dum);
	//Sum_m(f^i*DFT[w*DFT^-1[C^abT*DFT[v*f^j]]]^*)/(2l+1)
	fs_alm2cl(fl1->fs,bin,fl1->a_temp[iti],alm_2_dum,fl1->pol,fl2->pol,cl_dum,lmn_x,lmx_x,lmn_y,lmx_y);
	for(im1=0;im1<nspec;im1++) {
	  for(ip=0;ip<bin->n_bands;ip++)
	    cl_bias[im1][ip]-=cl_dum[im1][ip]*mij;
	}
      }
    }
  }

  if((fl1->ntemp>0) && (fl2->ntemp>0)) {
    int iti,itj,itp,itq,im1,im2;
    flouble *mat_prod=my_calloc(fl1->ntemp*fl2->ntemp,sizeof(flouble));
    for(itj=0;itj<fl1->ntemp;itj++) {
      for(itq=0;itq<fl2->ntemp;itq++) {
	//w*g^q
	for(im2=0;im2<fl2->nmaps;im2++)
	  fs_map_product(fl2->fs,fl2->temp[itq][im2],fl2->mask,map_2_dum[im2]);
	//DFT[w*g^q]
	fs_map2alm(fl2->fs,1,2*fl2->pol,map_2_dum,alm_2_dum);
	//C^ab*DFT[w*g^q]
	for(im1=0;im1<fl1->nmaps;im1++) {
	  fs_zero_alm(fl1->fs,alm_1_dum[im1]);
	  for(im2=0;im2<fl2->nmaps;im2++)
	    fs_alter_alm(fl2->fs,-1.,alm_2_dum[im2],alm_1_dum[im1],cl_proposal_f[im1*fl2->nmaps+im2],1);
	}
	//DFT^-1[C^ab*DFT[w*g^q]]
	fs_alm2map(fl1->fs,1,2*fl1->pol,map_1_dum,alm_1_dum);
	for(im1=0;im1<fl1->nmaps;im1++) {
	  //v*DFT^-1[C^ab*DFT[w*g^q]]
	  fs_map_product(fl1->fs,map_1_dum[im1],fl1->mask,map_1_dum[im1]);
	  //Int[f^jT*v*DFT^-1[C^ab*DFT[w*g^q]]]
	  mat_prod[itj*fl2->ntemp+itq]+=fs_map_dot(fl1->fs,map_1_dum[im1],fl1->temp[itj][im1]);
	}
      }
    }

    for(iti=0;iti<fl1->ntemp;iti++) {
      for(itp=0;itp<fl2->ntemp;itp++) {
	//Sum_m(f^i*g^p*)/(2l+1)
	fs_alm2cl(fl1->fs,bin,fl1->a_temp[iti],fl2->a_temp[itp],fl1->pol,fl2->pol,cl_dum,lmn_x,lmx_x,lmn_y,lmx_y);
	for(itj=0;itj<fl1->ntemp;itj++) {
	  double mij=gsl_matrix_get(fl1->matrix_M,iti,itj);
	  for(itq=0;itq<fl2->ntemp;itq++) {
	    double npq=gsl_matrix_get(fl2->matrix_M,itp,itq);
	    for(im1=0;im1<nspec;im1++) {
	      for(ip=0;ip<bin->n_bands;ip++)
		cl_bias[im1][ip]+=cl_dum[im1][ip]*mat_prod[itj*fl2->ntemp+itq]*mij*npq;
	    }
	  }
	}
      }
    }

    free(mat_prod);
  }

  for(ii=0;ii<fl1->nmaps;ii++) {
    dftw_free(map_1_dum[ii]);
    dftw_free(alm_1_dum[ii]);
  }
  free(map_1_dum);
  free(alm_1_dum);
  for(ii=0;ii<fl2->nmaps;ii++) {
    dftw_free(map_2_dum[ii]);
    dftw_free(alm_2_dum[ii]);
  }
  free(map_2_dum);
  free(alm_2_dum);
  for(ii=0;ii<nspec;ii++) {
    free(cl_dum[ii]);
    nmt_k_function_free(cl_proposal_f[ii]);
  }
  free(cl_proposal_f);
  free(cl_dum);

  return;
}

void nmt_couple_cl_l_flat(nmt_workspace_flat *w,int nl,flouble *larr,flouble **cl_in,flouble **cl_out)
{
  int ii;
  flouble **cell_in=my_malloc(w->ncls*sizeof(flouble *));
  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  for(ii=0;ii<w->ncls;ii++) {
    nmt_k_function *fcl=nmt_k_function_alloc(nl,larr,cl_in[ii],cl_in[ii][0],0.,0);
    cell_in[ii]=my_calloc(w->bin->n_bands,sizeof(flouble));

    int iy;
    flouble dkx=2*M_PI/w->fs->lx;
    flouble dky=2*M_PI/w->fs->ly;
    for(iy=0;iy<w->fs->ny;iy++) {
      flouble ky;
      int ik=0;
      if(2*iy<=w->fs->ny)
	ky=iy*dky;
      else
	ky=-(w->fs->ny-iy)*dky;
      if((ky>w->ellcut_y[1]) || (ky<w->ellcut_y[0])) {
	int ix;
	for(ix=0;ix<w->fs->nx;ix++) {
	  flouble kx;
	  if(2*ix<=w->fs->nx)
	    kx=ix*dkx;
	  else
	    kx=-(w->fs->nx-ix)*dkx;
	  if((kx>w->ellcut_x[1]) || (kx<w->ellcut_x[0])) {
	    double kmod=sqrt(kx*kx+ky*ky);
	    ik=nmt_bins_flat_search_fast(w->bin,kmod,ik);
	    if(ik>=0)
	      cell_in[ii][ik]+=nmt_k_function_eval(fcl,kmod,intacc);
	  }
	}
      }
    }

    for(iy=0;iy<w->bin->n_bands;iy++) {
      if(w->n_cells[iy]>0)
	cell_in[ii][iy]/=w->n_cells[iy];
      else
	cell_in[ii][iy]=0;
    }
    nmt_k_function_free(fcl);
  }
  gsl_interp_accel_free(intacc);

  int icl1;
  for(icl1=0;icl1<w->ncls;icl1++) {
    int i1;
    for(i1=0;i1<w->bin->n_bands;i1++) {
      int icl2;
      int ind1=i1*w->ncls+icl1;
      cl_out[icl1][i1]=0;
      for(icl2=0;icl2<w->ncls;icl2++) {
	int i2;
	for(i2=0;i2<w->bin->n_bands;i2++) {
	  int ind2=i2*w->ncls+icl2;
	  cl_out[icl1][i1]+=w->coupling_matrix_unbinned[ind1][ind2]*cell_in[icl2][i2];
	}
      }
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
    for(ib2=0;ib2<w->bin->n_bands;ib2++) {
      gsl_vector_set(dl_map_bad_b,w->ncls*ib2+icl,
		     cl_in[icl][ib2]-cl_noise_in[icl][ib2]-cl_bias[icl][ib2]);
    }
  }
  
  gsl_linalg_LU_solve(w->coupling_matrix_binned,w->coupling_matrix_perm,dl_map_bad_b,dl_map_good_b);
  for(icl=0;icl<w->ncls;icl++) {
    for(ib2=0;ib2<w->bin->n_bands;ib2++)
      cl_out[icl][ib2]=gsl_vector_get(dl_map_good_b,w->ncls*ib2+icl);
  }

  gsl_vector_free(dl_map_bad_b);
  gsl_vector_free(dl_map_good_b);
}

void nmt_compute_coupled_cell_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
				   nmt_binning_scheme_flat *bin,flouble **cl_out,
				   flouble lmn_x,flouble lmx_x,flouble lmn_y,flouble lmx_y)
{
  if(check_flatsky_infos(fl1->fs,fl2->fs))
    report_error(1,"Can only correlate fields defined on the same pixels!\n");
  fs_alm2cl(fl1->fs,bin,fl1->alms,fl2->alms,fl1->pol,fl2->pol,cl_out,lmn_x,lmx_x,lmn_y,lmx_y);
}

nmt_workspace_flat *nmt_compute_coupling_matrix_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
						     nmt_binning_scheme_flat *bin,
						     flouble lmn_x,flouble lmx_x,
						     flouble lmn_y,flouble lmx_y)
{
  if(check_flatsky_infos(fl1->fs,fl2->fs))
    report_error(1,"Can only correlate fields defined on the same pixels!\n");

  return nmt_compute_coupling_matrix_flat_c(fl1,fl2,bin,lmn_x,lmx_x,lmn_y,lmx_y);
}

nmt_workspace_flat *nmt_compute_power_spectra_flat(nmt_field_flat *fl1,nmt_field_flat *fl2,
						   nmt_binning_scheme_flat *bin,
						   flouble lmn_x,flouble lmx_x,flouble lmn_y,flouble lmx_y,
						   nmt_workspace_flat *w0,flouble **cl_noise,
						   int nl_prop,flouble *l_prop,flouble **cl_prop,
						   flouble **cl_out)
{
  int ii;
  flouble **cl_bias,**cl_data;
  nmt_workspace_flat *w;

  if(w0==NULL)
    w=nmt_compute_coupling_matrix_flat(fl1,fl2,bin,lmn_x,lmx_x,lmn_y,lmx_y);
  else
    w=w0;

  cl_bias=my_malloc(w->ncls*sizeof(flouble *));
  cl_data=my_malloc(w->ncls*sizeof(flouble *));
  for(ii=0;ii<w->ncls;ii++) {
    cl_bias[ii]=my_calloc(w->bin->n_bands,sizeof(flouble));
    cl_data[ii]=my_calloc(w->bin->n_bands,sizeof(flouble));
  }
  nmt_compute_coupled_cell_flat(fl1,fl2,bin,cl_data,lmn_x,lmx_x,lmn_y,lmx_y);
  nmt_compute_deprojection_bias_flat(fl1,fl2,bin,lmn_x,lmx_x,lmn_y,lmx_y,
				     nl_prop,l_prop,cl_prop,cl_bias);
  nmt_decouple_cl_l_flat(w,cl_data,cl_noise,cl_bias,cl_out);
  for(ii=0;ii<w->ncls;ii++) {
    free(cl_bias[ii]);
    free(cl_data[ii]);
  }
  free(cl_bias);
  free(cl_data);

  return w;
}
