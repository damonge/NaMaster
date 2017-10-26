#include <stdio.h>

#define CTEST_MAIN
#include "namaster.h"
#include "ctest.h"

#define NSIDE_TESTS 128

int *get_sequence(int n0,int nf)
{
  int i;
  int *seq=malloc((nf-n0)*sizeof(int));
  ASSERT_NOT_NULL(seq);
  for(i=0;i<(nf-n0);i++)
    seq[i]=n0+i;
  return seq;
}

CTEST(nmt,bins_constant) {
  nmt_binning_scheme *bin=nmt_bins_constant(4,2000);
  ASSERT_EQUAL(bin->n_bands,499);
  ASSERT_EQUAL(bin->ell_list[5][2],2+4*5+2);
  nmt_bins_free(bin);
}

CTEST(nmt,bins_var) {
  int i,j;
  nmt_binning_scheme *bin=nmt_bins_constant(4,2000);
  int *ells=get_sequence(2,1998);
  int *bpws=malloc(1996*sizeof(int));
  double *weights=malloc(1996*sizeof(double));
  for(i=0;i<1996;i++) {
    bpws[i]=i/4;
    weights[i]=0.25;
  }
  nmt_binning_scheme *bin2=nmt_bins_create(1996,bpws,ells,weights,2000);

  ASSERT_EQUAL(bin->n_bands,499);
  ASSERT_EQUAL(bin->n_bands,bin2->n_bands);
  ASSERT_EQUAL(bin->ell_list[5][2],2+4*5+2);
  ASSERT_EQUAL(bin2->ell_list[5][2],2+4*5+2);
  free(bpws);
  free(ells);
  free(weights);
  nmt_bins_free(bin);
  nmt_bins_free(bin2);
}

CTEST(nmt,field_t) {
  int i;
  int spins[1]={0};
  int *ell=get_sequence(0,3*NSIDE_TESTS);
  double *cl=malloc(3*NSIDE_TESTS*sizeof(double));
  double *bm=malloc(3*NSIDE_TESTS*sizeof(double));
  for(i=0;i<3*NSIDE_TESTS;i++) {
    cl[i]=pow((ell[i]+1.)/101.,-1.2);
    bm[i]=1.;
  }
  double **map=nmt_synfast_sph(NSIDE_TESTS,1,spins,3*NSIDE_TESTS-1,&cl,&bm,1234);
  double *mask=malloc(12*NSIDE_TESTS*NSIDE_TESTS*sizeof(double));
  for(i=0;i<12*NSIDE_TESTS*NSIDE_TESTS;i++)
    mask[i]=1.;
  nmt_field *fld=nmt_field_alloc_sph(NSIDE_TESTS,mask,0,map,0,NULL,NULL,0,0);
  ASSERT_EQUAL(fld->nside,NSIDE_TESTS);
  ASSERT_EQUAL(fld->npix,12*NSIDE_TESTS*NSIDE_TESTS);
  nmt_field_free(fld);
  free(mask);
  free(map[0]);
  free(map);
  free(ell);
  free(cl);
}

CTEST(nmt,field_r) {
  nmt_field *fld=nmt_field_read("test/mask.fits","test/maps.fits","none","none",0,0,0);
  ASSERT_EQUAL(fld->nside,256);
  nmt_field_free(fld);
}

int main(int argc,const char *argv[])
{
  int result=ctest_main(argc,argv);
  return result;
}
