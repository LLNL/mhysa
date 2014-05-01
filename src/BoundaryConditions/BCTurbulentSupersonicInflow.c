#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>

/* 
 * Turbulent Supersonic Inflow BC (specific to NavierStokes3D)
 * Used for turbulent supersonic inflow into the domain.
 *
 * A mean flow is specified and the flow fluctuations are 
 * added on this mean state.
 *
 * boundary->var is irrelevant, it acts on on the components
 * so no need to specify it for each component, just specify
 * it once with an arbitrary value for boundary->var 
*/

int BCTurbulentSupersonicInflowU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;

  double *inflow_data = boundary->UnsteadyDirichletData;
  int    *inflow_size = boundary->UnsteadyDirichletSize;

  if (ndims == 3) {

    /* create a fake physics object */
    double gamma;
    gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->on_this_proc) {
      /* the following bit is hardcoded for the inflow data
       * representing fluctuations in a domain of length 2pi */
      double  xt = boundary->FlowVelocity[dim] * waqt;
      int     N  = inflow_size[dim];
      double  L  = 2.0 * (4.0*atan(1.0));
      int     it = ((int) ((xt/L) * ((double)N))) % N;

      int bounds[ndims], indexb[ndims], indexi[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);
      int done = 0;
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        
        /* set the ghost point values - mean flow */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt      = boundary->FlowDensity;
        pressure_gpt = boundary->FlowPressure;
        uvel_gpt     = boundary->FlowVelocity[0];
        vvel_gpt     = boundary->FlowVelocity[1];
        wvel_gpt     = boundary->FlowVelocity[2];

        /* calculate the turbulent fluctuations */
        double duvel , dvvel , dwvel ;
        int index1[ndims]; _ArrayCopy1D_(indexb,index1,ndims);
        index1[dim] = it;
        int q; _ArrayIndex1D_(ndims,inflow_size,index1,0,q);
        duvel = inflow_data[q*nvars+1];
        dvvel = inflow_data[q*nvars+2];
        dwvel = inflow_data[q*nvars+3];

        /* add the turbulent fluctuations to the velocity field */
        uvel_gpt      += duvel;
        vvel_gpt      += dvvel;
        wvel_gpt      += dwvel;

        /* set the ghost point values */
        energy_gpt   = inv_gamma_m1*pressure_gpt
                       + 0.5 * rho_gpt 
                       * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);
        phi[nvars*p1+0] = rho_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt;
        phi[nvars*p1+3] = rho_gpt * wvel_gpt;
        phi[nvars*p1+4] = energy_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}

int BCTurbulentSupersonicInflowDU(void *b,void *m,int ndims,int nvars,int *size,int ghosts,double *phi,double *phi_ref,double waqt)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  int             v;

  if (boundary->on_this_proc) {
    int bounds[ndims], indexb[ndims];
    _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
    _ArraySetValue_(indexb,ndims,0); 
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,size  ,indexb,boundary->is,ghosts,p);
      for (v=0; v<nvars; v++) phi[nvars*p+v] = 0.0;
      _ArrayIncrementIndex_(ndims,bounds,indexb,done);
    }
  }
  return(0);
}

int BCReadTurbulentInflowData(void *b,void *m,int ndims,int nvars,int *DomainSize)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  char    *filename     = boundary->UnsteadyDirichletFilename;
  int     *inflow_size  = NULL;
  double  *inflow_data  = NULL;
  double  *buffer       = NULL;
  
  int dim = boundary->dim;
  int face= boundary->face;
  int d;

  if (!mpi->rank) {

    printf("Reading turbulent inflow boundary data from %s.\n",filename);

    FILE *in;
    int  ferr;

    /* calculate the number of processors that sit on unsteady boundary */
    int nproc = 1;
    for (d=0; d<ndims; d++) nproc *= mpi->iproc[d]; nproc /= mpi->iproc[dim];

    in = fopen(filename,"rb");
    if (!in) {
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): cannot open unsteady boundary data file %s.\n",filename);
      return(1);
    }
    int count = 0;
    while ((!feof(in)) && (count < nproc)) {
      int rank[ndims], size[ndims];
      ferr = fread(rank,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (1) in file reading, count %d.\n",count);
        return(1);
      }
      if (rank[dim] != (face > 0 ? 0 : mpi->iproc[dim]-1) ) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (2) in file reading, count %d.\n",count);
        return(1);
      }
      ferr = fread(size,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (3) in file reading, count %d.\n",count);
        return(1);
      }
      int flag = 1;
      for (d=0; d<ndims; d++) if ((d != dim) && (size[d] != DomainSize[d])) flag = 0;
      if (!flag) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (4) (dimension mismatch) in file reading, count %d.\n",count);
        return(1);
      }

      int data_size = nvars;
      for (d=0; d<ndims; d++) data_size *= size[d];
      buffer = (double*) calloc (data_size,sizeof(double));
      ferr = fread(buffer,sizeof(double),data_size,in);
      if (ferr != data_size) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (6) in file reading, count %d.\n",count);
        return(1);
      }

      int rank1D = MPIRank1D(ndims,mpi->iproc,rank);

      if (!rank1D) {

        int index[ndims];
        inflow_size = (int*) calloc (ndims, sizeof(int));
        _ArrayCopy1D_(size,inflow_size,ndims);
        inflow_data = (double*) calloc (data_size, sizeof(double));
        ArrayCopynD(ndims,buffer,inflow_data,size,0,0,index,nvars);

      } else {
#ifndef serial
        MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
        MPI_Isend(size,ndims,MPI_INT,rank1D,2152,mpi->world,&req[0]);
        MPI_Isend(buffer,data_size,MPI_DOUBLE,rank1D,2153,mpi->world,&req[1]);
        MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
#else
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): This is a serial run. Invalid (non-zero) rank read.\n");
#endif
      }

      free(buffer);
      count++;
    }

    if (count < nproc) {
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): missing data in unsteady boundary data file %s.\n",filename);
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): should contain data for %d processors, ", nproc);
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): but contains data for %d processors!\n", count);
      return(1);
    }

    fclose(in);

  } else {
#ifndef serial
    if (mpi->ip[dim] == (face > 0 ? 0 : mpi->iproc[dim]-1) ) {
      MPI_Request req = MPI_REQUEST_NULL;
      inflow_size = (int*) calloc (ndims,sizeof(int));
      MPI_Irecv(inflow_size,ndims,MPI_INT,0,2152,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
      int data_size = nvars;
      for (d=0; d<ndims; d++) data_size *= inflow_size[d];
      inflow_data = (double*) calloc (data_size,sizeof(double));
      MPI_Irecv(inflow_data,data_size,MPI_DOUBLE,0,2153,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
    }
#else
    fprintf(stderr,"Error in BCReadTurbulentInflowData(): Serial code should not be here!.\n");
#endif
  }
  
  boundary->UnsteadyDirichletSize = inflow_size;
  boundary->UnsteadyDirichletData = inflow_data;

  return(0);
}

