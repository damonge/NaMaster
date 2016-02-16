#include "common_mvqe.h"

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

int cg_solve(ParamCG *cg,flouble *b,flouble *x0,flouble *xf,flouble tol,int max_iter,void *params)
{
  int ii,converged;
  flouble inv_bmod,tol2;

  tol2=tol*tol;
  inv_bmod=0;
  inv_bmod=cg->dotpr(cg->n,b,b);
  inv_bmod=1./inv_bmod;

  memcpy(cg->xk,x0,cg->n*sizeof(flouble)); //x0=x0
  cg->action_A(cg->xk,cg->rk,params); //r0=A*x0
  cg->add_vscal(cg->n,b,-1.,cg->rk,cg->rk); //r0=b-r0
  cg->action_M(cg->rk,cg->zk,params); //z0=M*r0

  memcpy(cg->pk,cg->zk,cg->n*sizeof(flouble)); //p0=z0

  ii=0;
  converged=1;
  while(1>0) {
    flouble rz,ak,bk;
    flouble restot;
    if(ii>max_iter) {
      converged=0;
      break;
    }      
    rz=cg->dotpr(cg->n,cg->rk,cg->zk);
    cg->action_A(cg->pk,cg->vdum,params);
    ak=rz/cg->dotpr(cg->n,cg->vdum,cg->pk); //a_k=rk*zk/(pk*A*pk)
    cg->add_vscal(cg->n,cg->xk,ak,cg->pk,cg->xk); //x_{k+1}=x_k+a_k*pk
    cg->add_vscal(cg->n,cg->rk,-ak,cg->vdum,cg->rk); //r_{k+1}=r_k-a_k*A*pk
    restot=cg->dotpr(cg->n,cg->rk,cg->rk)*inv_bmod;
    //    printf("%d %lE\n",ii,restot);
    if(restot<tol2)
      break;
    cg->action_M(cg->rk,cg->zk,params); //z_{k+1}=M*r_{k+1}
    bk=cg->dotpr(cg->n,cg->rk,cg->zk)/rz; //b_k=r_{k+1}*z_{k+1}/(rk*zk)
    cg->add_vscal(cg->n,cg->zk,bk,cg->pk,cg->pk); //p_{k+1}=z_{k+1}+b_k*p_k
    ii++;
  }

  //  printf(" CG converged after %d iterations\n",ii);
  memcpy(xf,cg->xk,cg->n*sizeof(flouble));

  return converged;
}
