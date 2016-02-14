#include "common_mvqe.h"

typedef struct {
  int n;
  flouble (*dotpr)(int *,flouble *,flouble *);
  void (*action_A)(flouble *,flouble *,void *);
  void (*action_E)(flouble *,flouble *,void *);
  flouble *xk;
  flouble *rk;
  flouble *zk;
  flouble *pk;
  flouble *vdum;
} ParamCG;

ParamCG *cg_param_new(int n)
{
  ParamCG *cg=my_malloc(sizeof(ParamCG));
  cg->n=n;
  cg->xk=my_malloc(n*sizeof(flouble));
  cg->rk=my_malloc(n*sizeof(flouble));
  cg->zk=my_malloc(n*sizeof(flouble));
  cg->pk=my_malloc(n*sizeof(flouble));
  cg->vdum=my_malloc(n*sizeof(flouble));
  
  return cg;
}

void cg_param_free(ParamCG *cg)
{
  if(cg->n>0) {
    free(cg->xk);
    free(cg->rk);
    free(cg->zk);
    free(cg->pk);
    free(cg->vdum);
  }

  free(cg);
}

void cg_solve(ParamCG *cg,flouble *b,flouble *x0,flouble *xf,flouble tol,int max_iter,void *params)
{
  int ii,restot,converged;
  flouble inv_bmod,tol2;

  tol2=tol*tol;
  inv_bmod=0;
  for(ii=0;ii<cg->n;ii++)
    inv_bmod+=b[ii]*b[ii];
  inv_bmod=1./inv_bmod;

  memcpy(cg->xk,x0,cg->n*sizeof(flouble)); //x0=x0
  cg->action_A(cg->xk,cg->rk,params);
  for(ii=0;ii<cg->n;ii++)
    cg->rk[ii]=b[ii]-cg->rk[ii]; //r0=b-A*x0
  cg->action_M(cg->rk,cg->zk,params); //z0=M*r0

  memcpy(cg->pk,cg->zk,cg->n*sizeof(flouble)); //p0=z0

  ii=0;
  converged=1;
  restot=2*tol2;
  while(1>0) {
    int jj;
    flouble rz,ak,bk;
    if(iter>max_iter) {
      converged=0;
      break;
    }      
    rz=cg->dotpr(cg->n,cg->rk,cg->zk);
    cg->action_A(cg->pk,cg->vdum,params);
    ak=rz/cg->dotpr(cg->n,cg->vdum,cg->pk); //a_k=rk*zk/(pk*A*pk)
    for(jj=0;jj<cg->n;jj++) {
      cg->xk[jj]+=ak*cg->pk[jj]; //x_{k+1}=x_k+a_k*p_k
      cg->rk[jj]-=ak*cg->vdum[jj]; //r_{k+1}=r_k-a_k*A*p_k
    }
    restot=cg->dotpr(cg->n,cg->rk,cg->rk)*inv_bmod;
    if(restot<tol2)
      break;
    cg->action_M(cg->rk,cg->zk,params); //z_{k+1}=M*r_{k+1}
    bk=cg->dotpr(cg->n,cg->rk,cg->zk)/rz; //b_k=r_{k+1}*z_{k+1}/(rk*zk)
    for(jj=0;jj<cg->n;jj++)
      cg->pk[jj]=cg->zk[jj]+bk*cg->pk[jj]; //p_{k+1}=z_{k+1}+b_k*p_k
    ii++;
  }

  memcpy(xf,cg->xk,cg->n*sizeof(flouble));
}
