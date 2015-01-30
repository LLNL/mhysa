#ifdef with_scalapack

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#ifndef serial
#include <mpi.h>
#endif
#include <tridiagLU.h>

extern void pddtsv_();

int tridiagScaLPK(double *a,double *b,double *c,double *x,
              int n,int ns,void *r,void *m)
{
  TridiagLU       *params = (TridiagLU*) r;
  int             rank,nproc,nglobal,nrhs,i,s,ia,ib,desca[9],descb[9],err,
                  lwork;
  double          *dl,*d,*du,*rhs,*work;
  struct timeval  start,end;

#ifdef serial
  rank  = 0;
  nproc = 1;
  nglobal=n;
#else
  MPI_Comm        *comm = (MPI_Comm*) m;
  int             ierr = 0;
  const int       nvar = 4;

  if (comm) {
    MPI_Comm_size(*comm,&nproc);
    MPI_Comm_rank(*comm,&rank);
  } else {
    rank  = 0;
    nproc = 1;
  }
  MPI_Allreduce(&n,&nglobal,1,MPI_INT,MPI_SUM,*comm);
#endif

  /* check */
  if (nglobal%n != 0) {
    if (!rank) {
      fprintf(stderr,"Error: The ScaLAPACK wrapper can only handle cases where the global ");
      fprintf(stderr,"size of system is an integer multiple of no. of processes.\n");
    }
    return(1);
  }

  if (!params) {
    fprintf(stderr,"Error in tridiagLU(): NULL pointer passed for parameters.\n");
    return(1);
  }

  /* start */
  gettimeofday(&start,NULL);

  nrhs = 1;
  ia = 1;
  ib = 1;

  lwork = (12*nproc+3*n) + ( (8*nproc) > (10*nproc+4*nrhs) ? (8*nproc) : (10*nproc+4*nrhs) );

  desca[0] = 501;
  desca[1] = params->blacs_ctxt;
  desca[2] = nglobal;
  desca[3] = n;
  desca[4] = 0;
  desca[5] = 0;
  desca[6] = 0;
  desca[7] = 0;
  desca[8] = 0;

  descb[0] = 502;
  descb[1] = params->blacs_ctxt;
  descb[2] = nglobal;
  descb[3] = n;
  descb[4] = 0;
  descb[5] = n;
  descb[6] = 0;
  descb[7] = 0;
  descb[8] = 0;

  dl    = (double*) calloc (n,sizeof(double));
  d     = (double*) calloc (n,sizeof(double));
  du    = (double*) calloc (n,sizeof(double));
  rhs   = (double*) calloc (n,sizeof(double));
  work  = (double*) calloc(lwork,sizeof(double));

  for (s=0; s<ns; s++) {

    for (i=0; i<n; i++) {
      dl[i] = a[i*ns+s];
      d [i] = b[i*ns+s];
      du[i] = c[i*ns+s];
      rhs[i]= x[i*ns+s];
    }

    /* call the ScaLAPACK function */
    pddtsv_(&nglobal,&nrhs,dl,d,du,&ia,desca,rhs,&ib,descb,work,&lwork,&err);
    if (err) return(err);
    
    for (i=0; i<n; i++) x[i*ns+s] = rhs[i];
  }

  /* end of stage 4 */
  gettimeofday(&end,NULL);

  free(dl);
  free(d);
  free(du);
  free(rhs);
  free(work);

  /* save runtimes if needed */
  params->stage1_time = 0.0;
  params->stage2_time = 0.0;
  params->stage3_time = 0.0;
  params->stage4_time = 0.0;
  long long walltime;
  walltime = ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
  params->total_time = (double) walltime / 1000000.0;
  return(0);
}

#endif
