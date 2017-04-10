#include "utils.h"

nmt_k_function *nmt_k_function_alloc(int nk,flouble *karr,flouble *farr,flouble y0,flouble yf,int is_const)
{
  nmt_k_function *f=my_malloc(sizeof(nmt_k_function));
  f->is_const=is_const;
  f->y0=y0;
  if(!is_const) {
    f->x0=karr[0];
    f->xf=karr[nk-1];
    f->yf=yf;
    f->spl=gsl_spline_alloc(gsl_interp_cspline,nk);
    gsl_spline_init(f->spl,karr,farr,nk);
  }
  return f;
}

void nmt_k_function_free(nmt_k_function *f)
{
  if(!(f->is_const))
    gsl_spline_free(f->spl);
  free(f);
}

flouble nmt_k_function_eval(nmt_k_function *f,flouble k,gsl_interp_accel *intacc)
{
  if((f->is_const) || (k<=f->x0))
      return f->y0;
  else if(k>=f->xf)
    return f->yf;
  else
    return gsl_spline_eval(f->spl,k,intacc);
}

#define N_DELL 2
nmt_flatsky_info *nmt_flatsky_info_alloc(int nx,int ny,flouble lx,flouble ly)
{
  nmt_flatsky_info *fs=my_malloc(sizeof(nmt_flatsky_info));
  fs->nx=nx;
  fs->ny=ny;
  fs->npix=nx*ny;
  fs->lx=lx;
  fs->ly=ly;
  fs->pixsize=lx*ly/(nx*ny);

  int ii;
  flouble dkx=2*M_PI/lx;
  flouble dky=2*M_PI/ly;
  flouble kmax_x=dkx*(nx/2);
  flouble kmax_y=dky*(ny/2);
  double dk=NMT_MAX(dkx,dky);
  double kmax=NMT_MAX(kmax_y,kmax_x);
  fs->dell=N_DELL*dk;
  fs->i_dell=1./fs->dell;
  fs->n_ell=0;
  while((fs->n_ell+1)*fs->dell<=kmax)
    fs->n_ell++;
  fs->ell_min=my_malloc(fs->n_ell*sizeof(flouble));
  fs->n_cells=my_calloc(fs->n_ell,sizeof(int));
  for(ii=0;ii<fs->n_ell;ii++)
    fs->ell_min[ii]=ii*fs->dell;

#pragma omp parallel default(none) \
  shared(fs,dkx,dky)
  {
    int iy;

#pragma omp for
    for(iy=0;iy<fs->ny;iy++) {
      int ix;
      flouble ky;
      if(2*iy<=fs->ny)
	ky=iy*dky;
      else
	ky=-(fs->ny-iy)*dky;
      for(ix=0;ix<=fs->nx/2;ix++) {
	flouble kx=ix*dkx;
	flouble kmod=sqrt(kx*kx+ky*ky);
	int ik=(int)(kmod*fs->i_dell);
	if(ik<fs->n_ell) {
#pragma omp atomic
	  fs->n_cells[ik]++;
	}
      }
    } //end omp for
  } //end omp parallel

  return fs;
}

void nmt_flatsky_info_free(nmt_flatsky_info *fs)
{
  free(fs->ell_min);
  free(fs->n_cells);
  free(fs);
}

void nmt_field_flat_free(nmt_field_flat *fl)
{
  int imap,itemp;
  nmt_k_function_free(fl->beam);

  nmt_flatsky_info_free(fl->fs);
  for(imap=0;imap<fl->nmaps;imap++) {
    dftw_free(fl->maps[imap]);
    dftw_free(fl->alms[imap]);
  }
  dftw_free(fl->mask);
  if(fl->ntemp>0) {
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      for(imap=0;imap<fl->nmaps;imap++)
	dftw_free(fl->temp[itemp][imap]);
    }
  }

  free(fl->alms);
  free(fl->maps);
  if(fl->ntemp>0) {
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      free(fl->temp[itemp]);
      free(fl->a_temp[itemp]);
    }
    free(fl->temp);
    free(fl->a_temp);
    gsl_matrix_free(fl->matrix_M);
  }
  free(fl);
}

static void nmt_purify_flat(nmt_field_flat *fl)
{
  return; //Placeholder
}

nmt_field_flat *nmt_field_flat_alloc(int nx,int ny,flouble lx,flouble ly,
				     flouble *mask,int pol,flouble **maps,int ntemp,flouble ***temp,
				     int nl_beam,flouble *l_beam,flouble *beam,
				     int pure_e,int pure_b)
{
  long ip;
  int ii,itemp,itemp2,imap;
  nmt_field_flat *fl=my_malloc(sizeof(nmt_field_flat));
  fl->fs=nmt_flatsky_info_alloc(nx,ny,lx,ly);
  fl->npix=nx*ny;
  fl->pol=pol;
  if(pol) fl->nmaps=2;
  else fl->nmaps=1;
  fl->ntemp=ntemp;

  fl->pure_e=0;
  fl->pure_b=0;
  if(pol) {
    if(pure_e)
      fl->pure_e=1;
    if(pure_b)
      fl->pure_b=1;
  }

  if(beam==NULL)
    fl->beam=nmt_k_function_alloc(-1,NULL,NULL,1.,1.,1);
  else
    fl->beam=nmt_k_function_alloc(nl_beam,l_beam,beam,beam[0],0.,0);

  fl->mask=dftw_malloc(fl->npix*sizeof(flouble));
  for(ip=0;ip<fl->npix;ip++)
    fl->mask[ip]=mask[ip];

  fl->maps=my_malloc(fl->nmaps*sizeof(flouble *));
  for(ii=0;ii<fl->nmaps;ii++) {
    fl->maps[ii]=dftw_malloc(fl->npix*sizeof(flouble));
    for(ip=0;ip<fl->npix;ip++)
      fl->maps[ii][ip]=maps[ii][ip];
  }

  if(fl->ntemp>0) {
    fl->temp=my_malloc(fl->ntemp*sizeof(flouble **));
    fl->a_temp=my_malloc(fl->ntemp*sizeof(fcomplex **));
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      fl->temp[itemp]=my_malloc(fl->nmaps*sizeof(flouble *));
      fl->a_temp[itemp]=my_malloc(fl->nmaps*sizeof(fcomplex *));
      for(imap=0;imap<fl->nmaps;imap++) {
	fl->temp[itemp][imap]=dftw_malloc(fl->npix*sizeof(flouble));
	fl->a_temp[itemp][imap]=dftw_malloc(fl->fs->ny*(fl->fs->nx/2+1)*sizeof(fcomplex));
	for(ip=0;ip<fl->npix;ip++)
	  fl->temp[itemp][imap][ip]=temp[itemp][imap][ip];
	fs_map_product(fl->fs,fl->temp[itemp][imap],fl->mask,fl->temp[itemp][imap]); //Multiply by mask
      }
      fs_map2alm(fl->fs,1,2*fl->pol,fl->temp[itemp],fl->a_temp[itemp]);
    }
    
    //Compute normalization matrix
    fl->matrix_M=gsl_matrix_alloc(fl->ntemp,fl->ntemp);
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      for(itemp2=itemp;itemp2<fl->ntemp;itemp2++) {
	flouble matrix_element=0;
	for(imap=0;imap<fl->nmaps;imap++)
	  matrix_element+=fs_map_dot(fl->fs,fl->temp[itemp][imap],fl->temp[itemp2][imap]);
	gsl_matrix_set(fl->matrix_M,itemp,itemp2,matrix_element);
	if(itemp2!=itemp)
	  gsl_matrix_set(fl->matrix_M,itemp2,itemp,matrix_element);
      }
    }
    gsl_linalg_cholesky_decomp(fl->matrix_M); //TODO: this won't necessarily be invertible
    gsl_linalg_cholesky_invert(fl->matrix_M);
  }

  fl->alms=my_malloc(fl->nmaps*sizeof(fcomplex *));
  for(ii=0;ii<fl->nmaps;ii++)
    fl->alms[ii]=dftw_malloc(fl->fs->ny*(fl->fs->nx/2+1)*sizeof(fcomplex));

  if(fl->pol && (fl->pure_e || fl->pure_b))
    nmt_purify_flat(fl);
  else {
    for(ii=0;ii<fl->nmaps;ii++)
      fs_map_product(fl->fs,fl->maps[ii],fl->mask,fl->maps[ii]);
    fs_map2alm(fl->fs,1,2*fl->pol,fl->maps,fl->alms);
  }

  if(fl->ntemp>0) {
    //Deproject
    flouble *prods=my_calloc(fl->ntemp,sizeof(flouble));
    for(itemp=0;itemp<fl->ntemp;itemp++) {
      for(imap=0;imap<fl->nmaps;imap++) 
	prods[itemp]+=fs_map_dot(fl->fs,fl->temp[itemp][imap],fl->maps[imap]);
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

flouble **nmt_synfast_flat(int nx,int ny,flouble lx,flouble ly,int nfields,int *spin_arr,
			   int nl_beam,flouble *l_beam,flouble **beam_fields,
			   int nl_cell,flouble *l_cell,flouble **cell_fields,
			   int seed)
{
  int ifield,imap;
  int nmaps=0,ncls=0;
  long npix=nx*ny;
  nmt_k_function **beam,**cell;
  flouble **maps;
  fcomplex **alms;
  nmt_flatsky_info *fs=nmt_flatsky_info_alloc(nx,ny,lx,ly);
  for(ifield=0;ifield<nfields;ifield++) {
    int nmp=1;
    if(spin_arr[ifield]) nmp=2;
    nmaps+=nmp;
  }

  imap=0;
  beam=my_malloc(nmaps*sizeof(nmt_k_function *));
  maps=my_malloc(nmaps*sizeof(flouble *));
  for(ifield=0;ifield<nfields;ifield++) {
    int imp,nmp=1;
    if(spin_arr[ifield]) nmp=2;
    for(imp=0;imp<nmp;imp++) {
      beam[imap+imp]=nmt_k_function_alloc(nl_beam,l_beam,beam_fields[ifield],beam_fields[ifield][0],0.,0);
      maps[imap+imp]=dftw_malloc(npix*sizeof(flouble));
    }
    imap+=nmp;
  }

  ncls=(nmaps*(nmaps+1))/2;
  cell=my_malloc(ncls*sizeof(nmt_k_function *));
  for(imap=0;imap<ncls;imap++)
    cell[imap]=nmt_k_function_alloc(nl_cell,l_cell,cell_fields[imap],cell_fields[imap][0],0.,0);

  alms=fs_synalm(nx,ny,lx,ly,nmaps,cell,beam,seed);

  for(imap=0;imap<nmaps;imap++)
    nmt_k_function_free(beam[imap]);
  free(beam);
  for(imap=0;imap<ncls;imap++)
    nmt_k_function_free(cell[imap]);
  free(cell);

  imap=0;
  for(ifield=0;ifield<nfields;ifield++) {
    int imp,nmp=1;
    if(spin_arr[ifield]) nmp=2;
    fs_alm2map(fs,1,spin_arr[ifield],&(maps[imap]),&(alms[imap]));
    for(imp=0;imp<nmp;imp++)
      dftw_free(alms[imap+imp]);
    imap+=nmp;
  }
  free(alms);
  nmt_flatsky_info_free(fs);

  return maps;
}
