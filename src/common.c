#include "common.h"

int my_linecount(FILE *f)
{
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr," Fatal error: %s",msg);
    exit(level);
  }
  else
    fprintf(stderr," Warning: %s",msg);
}

void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) report_error(1,"Out of memory\n");

  return outptr;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL)
    report_error(1,"Out of memory\n");

  return outptr;
}

FILE *my_fopen(const char *path,const char *mode)
{
  FILE *fout=fopen(path,mode);
  if(fout==NULL)
    report_error(1,"Couldn't open file %s\n",path);

  return fout;
}

size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error fwriting\n");

  return nmemb;
}

size_t my_fread(void *ptr,size_t size,size_t count,FILE *stream)
{
  if(fread(ptr,size,count,stream)!=count)
    report_error(1,"Error freading\n");

  return count;
}

void free_bins(BinSchm *bins)
{
  int ii;
  free(bins->nell_list);
  for(ii=0;ii<bins->n_bands;ii++) {
    free(bins->ell_list[ii]);
    free(bins->w_list[ii]);
  }
  free(bins->ell_list);
  free(bins->w_list);
}

BinSchm *create_bins(int nlb,int nside)
{
  int ii,lmax=3*nside-1;
  int nband_max=(lmax-1)/nlb;
  flouble w0=1./nlb;

  BinSchm *bins=my_malloc(sizeof(BinSchm));
  bins->n_bands=nband_max;
  bins->nell_list=my_calloc(nband_max,sizeof(int));
  bins->ell_list=my_malloc(nband_max*sizeof(int *));
  bins->w_list=my_malloc(nband_max*sizeof(flouble *));

  for(ii=0;ii<nband_max;ii++) {
    int jj;
    bins->ell_list[ii]=my_malloc(nlb*sizeof(int));
    bins->w_list[ii]=my_malloc(nlb*sizeof(int));
    for(jj=0;jj<nlb;jj++) {
      bins->ell_list[ii][jj]=2+ii*nlb+jj;
      bins->w_list[ii][jj]=w0;
    }
  }

  return bins;
}

BinSchm *read_bins(char *fname,int nside)
{
  FILE *fi=my_fopen(fname,"r");
  int ii,nlines=my_linecount(fi); rewind(fi);
  if(nlines!=3*nside)
    report_error(1,"Error reading binning table\n");

  int *band_number,*larr;
  flouble *warr;
  band_number=my_malloc(nlines*sizeof(int));
  larr=my_malloc(nlines*sizeof(int));
  warr=my_malloc(nlines*sizeof(flouble));
  
  int nband_max=0;
  for(ii=0;ii<nlines;ii++) {
    double w;
    int stat=fscanf(fi,"%d %d %lf",&(band_number[ii]),&(larr[ii]),&w);
    if(stat!=3)
      report_error(1,"Error reading %s, line %d\n",fname,ii+1);
    warr[ii]=w;
    if(band_number[ii]>nband_max)
      nband_max=band_number[ii];
  }
  fclose(fi);
  nband_max++;

  BinSchm *bins=my_malloc(sizeof(BinSchm));
  bins->n_bands=nband_max;
  bins->nell_list=my_calloc(nband_max,sizeof(int));
  bins->ell_list=my_malloc(nband_max*sizeof(int *));
  bins->w_list=my_malloc(nband_max*sizeof(flouble *));
  
  for(ii=0;ii<nlines;ii++) {
    if(band_number[ii]>=0)
      bins->nell_list[band_number[ii]]++;
  }

  for(ii=0;ii<nband_max;ii++) {
    bins->ell_list[ii]=my_malloc(bins->nell_list[ii]*sizeof(int));
    bins->w_list[ii]=my_malloc(bins->nell_list[ii]*sizeof(flouble));
  }

  for(ii=0;ii<nlines;ii++)
    bins->nell_list[band_number[ii]]=0;

  for(ii=0;ii<nlines;ii++) {
    int l=larr[ii];
    int b=band_number[ii];
    flouble w=warr[ii];

    if(b>=0) {
      bins->ell_list[b][bins->nell_list[b]]=l;
      bins->w_list[b][bins->nell_list[b]]=w;
      bins->nell_list[b]++;
    }
  }

  for(ii=0;ii<nband_max;ii++) {
    int jj;
    flouble norm=0;
    for(jj=0;jj<bins->nell_list[ii];jj++)
      norm+=bins->w_list[ii][jj];
    if(norm<=0)
      report_error(1,"Weights in band %d are wrong\n",ii);
    for(jj=0;jj<bins->nell_list[ii];jj++)
      bins->w_list[ii][jj]/=norm;
  }

  free(larr);
  free(band_number);
  free(warr);

  return bins;
}
