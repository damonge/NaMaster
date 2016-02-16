#include "common_mvqe.h"

int main(int argc,char **argv)
{
  char fnameInit[256];
  if(argc!=2) {
    fprintf(stderr,"Usage: ./mvqe <params_file>\n");
    exit(1);
  }
  sprintf(fnameInit,"%s",argv[1]);

  ParamMVQE *par=read_params(fnameInit);
  flouble *map_invcov=my_calloc(par->npix,sizeof(flouble));
  flouble *fisher=my_calloc(par->nbins*par->nbins,sizeof(flouble));
  invert_covar(par->map,map_invcov,par);
  fisher_avg(par,fisher);

  //  int ii;
  //  create_unit_variance(par,pix_data);
  //  flouble mean=0,rms=0;
  //  for(ii=0;ii<par->npix_seen;ii++) {
  //    mean+=pix_data[ii];
  //    rms+=pix_data[ii]*pix_data[ii];
  //  }
  //  mean/=par->npix_seen;
  //  rms=sqrt(rms/par->npix_seen-mean*mean);
  //  printf("%lE %lE\n",mean,rms);

  free_param_mvqe(par);
  free(map_invcov);
  free(fisher);

  return 0;
}
