#include "utils.h"

void *dftw_malloc(size_t n)
{
#ifdef _SPREC
  void *p=fftwf_malloc(n);
#else //_SPREC
  void *p=fftw_malloc(n);
#endif //_SPREC
  if(p==NULL)
    report_error(1,"Ran out of memory\n");
  return p;
}

void dftw_free(void *p)
{
#ifdef _SPREC
  fftwf_free(p);
#else //_SPREC
  fftw_free(p);
#endif //_SPREC
}

void fs_map_product(nmt_flatsky_info *fs,flouble *mp1,flouble *mp2,flouble *mp_out)
{
#pragma omp parallel default(none)		\
  shared(fs,mp1,mp2,mp_out)
  {
    long ip;

#pragma omp for
    for(ip=0;ip<fs->npix;ip++) {
      mp_out[ip]=mp1[ip]*mp2[ip];
    } //end omp for
  } //end omp parallel
}

flouble fs_map_dot(nmt_flatsky_info *fs,flouble *mp1,flouble *mp2)
{
  double sum=0;

#pragma omp parallel default(none)		\
  shared(mp1,mp2,sum,fs)
  {
    long ip;
    double sum_thr=0;
    
#pragma omp for
    for(ip=0;ip<fs->npix;ip++) {
      sum_thr+=mp1[ip]*mp2[ip];
    } //end omp for

#pragma omp critical
    {
      sum+=sum_thr;
    } //end omp critical
  } //end omp parallel

  return (flouble)(sum*fs->pixsize);
}  

static void qu2eb(nmt_flatsky_info *fs,int spin,fcomplex **alm)
{
#pragma omp parallel default(none) \
  shared(fs,spin,alm)
  {
    int iy;
    fcomplex sig=-pow(I,spin);
    flouble dkx=2*M_PI/fs->lx;
    flouble dky=2*M_PI/fs->ly;

#pragma omp for
    for(iy=0;iy<fs->ny;iy++) {
      int ix;
      flouble ky;
      if(2*iy<=fs->ny)
	ky=iy*dky;
      else
	ky=-(fs->ny-iy)*dky;
      for(ix=0;ix<=fs->nx/2;ix++) {
	flouble e,b,csphi,ssphi,cph,sph;
	int s=0;
	flouble kx=ix*dkx;
	long index=ix+(fs->nx/2+1)*iy;
	flouble kmod2=kx*kx+ky*ky;

	if(kmod2<=0) {
	  cph=1;
	  sph=0;
	}
	else {
	  flouble i_kmod=1./sqrt(kmod2);
	  cph=kx*i_kmod;
	  sph=ky*i_kmod;
	}

	csphi=1; ssphi=0;
	while(s<spin) {
	  flouble c2=csphi*cph-ssphi*sph;
	  flouble s2=ssphi*cph+csphi*sph;
	  csphi=c2;
	  ssphi=s2;
	  s++;
	}
	e=sig*(alm[0][index]*csphi-alm[1][index]*ssphi);
	b=sig*(alm[0][index]*ssphi+alm[1][index]*csphi);
	alm[0][index]=e;
	alm[1][index]=b;
      }
    } //end omp for
  } //end omp parallel
}

static void eb2qu(nmt_flatsky_info *fs,int spin,fcomplex **alm)
{
#pragma omp parallel default(none) \
  shared(fs,spin,alm)
  { 
    int iy;
    fcomplex sig=-pow(-I,spin);
    flouble dkx=2*M_PI/fs->lx;
    flouble dky=2*M_PI/fs->ly;

#pragma omp for
    for(iy=0;iy<fs->ny;iy++) {
      int ix;
      flouble ky;
      if(2*iy<=fs->ny)
	ky=iy*dky;
      else
	ky=-(fs->ny-iy)*dky;
      for(ix=0;ix<=fs->nx/2;ix++) {
	flouble q,u,csphi,ssphi,cph,sph;
	int s=0;
	flouble kx=ix*dkx;
	long index=ix+(fs->nx/2+1)*iy;
	flouble kmod2=kx*kx+ky*ky;
	
	if(kmod2<=0) {
	  cph=1;
	  sph=0;
	}
	else {
	  flouble i_kmod=1./sqrt(kmod2);
	  cph=kx*i_kmod;
	  sph=ky*i_kmod;
	}

	csphi=1; ssphi=0;
	while(s<spin) {
	  flouble c2=csphi*cph-ssphi*sph;
	  flouble s2=ssphi*cph+csphi*sph;
	  csphi=c2;
	  ssphi=s2;
	  s++;
	}

	q=sig*( alm[0][index]*csphi+alm[1][index]*ssphi);
	u=sig*(-alm[0][index]*ssphi+alm[1][index]*csphi);
	alm[0][index]=q;
	alm[1][index]=u;
      }
    } //end omp for
  } //end omp parallel
}

void fs_map2alm(nmt_flatsky_info *fs,int ntrans,int spin,flouble **map,fcomplex **alm)
{
  //TODO init threads??
#ifdef _SPREC
  fftwf_plan plan_ft;
#else //_SPREC
  fftw_plan plan_ft;
#endif //_SPREC
  int imap,nmaps=1;
  if(spin)
    nmaps=2;
  
  for(imap=0;imap<nmaps*ntrans;imap++) {
#ifdef _SPREC
    plan_ft=fftwf_plan_dft_r2c_2d(fs->ny,fs->nx,map[imap],alm[imap],FFTW_ESTIMATE);
    fftwf_execute(plan_ft);
    fftwf_destroy_plan(plan_ft);
#else //_SPREC
    plan_ft=fftw_plan_dft_r2c_2d(fs->ny,fs->nx,map[imap],alm[imap],FFTW_ESTIMATE);
    fftw_execute(plan_ft);
    fftw_destroy_plan(plan_ft);
#endif //_SPREC

#pragma omp parallel default(none) \
  shared(fs,alm,imap)
    {
      long ipix;
      flouble norm=fs->lx*fs->ly/(2*M_PI*fs->nx*fs->ny);
#pragma omp for
      for(ipix=0;ipix<fs->ny*(fs->nx/2+1);ipix++) {
	alm[imap][ipix]*=norm;
      } //end omp for
    } //end omp parallel
  }

  if(nmaps>1) { //Q,U -> E,B
    for(imap=0;imap<ntrans*nmaps;imap+=nmaps)
      qu2eb(fs,spin,&(alm[imap]));
  }
}

void fs_alm2map(nmt_flatsky_info *fs,int ntrans,int spin,flouble **map,fcomplex **alm)
{
  //TODO init threads??
#ifdef _SPREC
  fftwf_plan plan_ft;
#else //_SPREC
  fftw_plan plan_ft;
#endif //_SPREC
  int imap,nmaps=1;
  if(spin)
    nmaps=2;
  
  if(nmaps>1) { //E,B -> Q,U
    for(imap=0;imap<ntrans*nmaps;imap+=nmaps)
      eb2qu(fs,spin,&(alm[imap]));
  }

  for(imap=0;imap<nmaps*ntrans;imap++) {
#ifdef _SPREC
    plan_ft=fftwf_plan_dft_c2r_2d(fs->ny,fs->nx,alm[imap],map[imap],FFTW_ESTIMATE);
    fftwf_execute(plan_ft);
    fftwf_destroy_plan(plan_ft);
#else //_SPREC
    plan_ft=fftw_plan_dft_c2r_2d(fs->ny,fs->nx,alm[imap],map[imap],FFTW_ESTIMATE);
    fftw_execute(plan_ft);
    fftw_destroy_plan(plan_ft);
#endif //_SPREC

#pragma omp parallel default(none)		\
  shared(fs,map,imap)
    {
      long ipix;
      flouble norm=2*M_PI/(fs->lx*fs->ly);
#pragma omp for
      for(ipix=0;ipix<fs->npix;ipix++) {
	map[imap][ipix]*=norm;
      } //end omp for
    } //end omp parallel
  }

  if(nmaps>1) { //Q,U -> E,B
    for(imap=0;imap<ntrans*nmaps;imap+=nmaps)
      qu2eb(fs,spin,&(alm[imap]));
  }
}

void fs_alm2cl(nmt_flatsky_info *fs,fcomplex **alms_1,fcomplex **alms_2,
	       int pol_1,int pol_2,flouble **cls)
{
  int i1,il,nmaps_1=1,nmaps_2=1;
  if(pol_1) nmaps_1=2;
  if(pol_2) nmaps_2=2;

  for(i1=0;i1<nmaps_1;i1++) {
    int i2;
    fcomplex *alm1=alms_1[i1];
    for(i2=0;i2<nmaps_2;i2++) {
      int il;
      fcomplex *alm2=alms_2[i2];
      int index_cl=i2+nmaps_2*i1;
      for(il=0;il<fs->n_ell;il++)
	cls[index_cl][il]=0;

#pragma omp parallel default(none)		\
  shared(fs,alm1,alm2,index_cl,cls)
      {
	int iy;
	flouble dkx=2*M_PI/fs->lx;
	flouble dky=2*M_PI/fs->ly;

#pragma omp for
	for(iy=0;iy<fs->ny;iy++) {
	  int ix;
	  flouble ky;
	  if(2*iy<=fs->ny)
	    ky=iy*dky;
	  else
	    ky=-(fs->ny-iy)*dky;
	  for(ix=0;ix<=fs->nx/2;ix++) {
	    int s=0;
	    flouble kx=ix*dkx;
	    long index=ix+(fs->nx/2+1)*iy;
	    flouble kmod=sqrt(kx*kx+ky*ky);
	    int ik=(int)(kmod*fs->i_dell);
	    if(ik<fs->n_ell) {
#pragma omp atomic
	      cls[index_cl][ik]+=(creal(alm1[index])*creal(alm2[index])+cimag(alm1[index])*cimag(alm2[index]));
	    }
	  }
	} //end omp for
      } //end omp parallel

      for(il<0;il<fs->n_ell;il++) {
	if(fs->n_cells[il]<=0)
	  cls[index_cl][il]=0;
	else
	  cls[index_cl][il]/=fs->n_cells[il];
      }
    }
  }
}

void fs_anafast(nmt_flatsky_info *fs,flouble **maps_1,flouble **maps_2,
		int pol_1,int pol_2,flouble **cls)
{
  int i1;
  fcomplex **alms_1,**alms_2;
  int nmaps_1=1,nmaps_2=1;
  if(pol_1) nmaps_1=2;
  if(pol_2) nmaps_2=2;

  alms_1=my_malloc(nmaps_1*sizeof(fcomplex *));
  for(i1=0;i1<nmaps_1;i1++)
    alms_1[i1]=dftw_malloc(fs->ny*(fs->nx/2+1)*sizeof(fcomplex));
  fs_alm2map(fs,1,2*pol_1,maps_1,alms_1);

  if(maps_1==maps_2)
    alms_2=alms_1;
  else {
    alms_2=my_malloc(nmaps_2*sizeof(fcomplex *));
    for(i1=0;i1<nmaps_2;i1++)
      alms_2[i1]=dftw_malloc(fs->ny*(fs->nx/2+1)*sizeof(fcomplex));
    fs_alm2map(fs,1,2*pol_2,maps_2,alms_2);
  }

  fs_alm2cl(fs,alms_1,alms_2,pol_1,pol_2,cls);

  for(i1=0;i1<nmaps_1;i1++)
    dftw_free(alms_1[i1]);
  free(alms_1);
  if(maps_1!=maps_2) {
    for(i1=0;i1<nmaps_2;i1++)
      dftw_free(alms_2[i1]);
    free(alms_2);
  }
}
