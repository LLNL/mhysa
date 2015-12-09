/*! @file MPIExchangeBoundariesnD.c
    @brief Exchange data and fill in ghost points for an n-dimensional array
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Exchange data across MPI ranks, and fill in ghost points for an n-dimensional array 
  (where \a n is the total number of spatial dimensions). If any of the physical boundaries
  are periodic, this function also exchanges data and fills in the ghost points for these 
  boundaries.

  The n-dimensional array must be stored in the memory as a single-index array, with the following order of mapping:
  + Number of variables (vector components)
  + Spatial dimension 0
  + Spatial dimension 1
  + ...
  + Spatial dimensions \a ndims-1

  For example, consider a 2D simulation (\a ndims = 2), of size \f$7 \times 3\f$, with \f$4\f$ vector components
  (\a nvars = 4). The following figure shows the layout (without the ghost points):
  @image html layout.png
  @image latex layout.eps width=0.9\textwidth

  The bold numbers in parentheses represent the 2D indices. The numbers below them are the indices of the array
  that correspond to that 2D location. Thus, elements 40,41,42, and 43 in the array are the 1st, 2nd, 3rd, and 
  4th vector components at location (1,3).

  If \f${\bf i}\left[{\rm ndims}\right]\f$ is an integer array representing an n-dimensional index 
  (for example, \f$\left(5,4\right)\f$ in 2D, \f$\left(3,5,2\right)\f$ in 3D), and the number of vector 
  components is \a nvars, then:
  + #_ArrayIndex1D_ computes the index \f$p\f$ in the array corresponding to \f${\bf i}\f$. In the above example, 
    \f${\bf i} = \left(1,3\right) \rightarrow p = 10\f$.
  + \a var[nvars*p+v] accesses the \a v-th component of the n-dimensional array \a var at location \f${\bf i}\f$.
    In the above example, to access the 3rd vector component at location \f$\left(1,3\right)\f$, we have \f$p=10\f$,
    so \a var [4*10+2] = \a var [42].
*/
int MPIExchangeBoundariesnD(
                              int     ndims,  /*!< Number of spatial dimensions */
                              int     nvars,  /*!< Number of variables (vector components) at each grid location */
                              int     *dim,   /*!< Integer array whose elements are the local size along each spatial dimension */
                              int     ghosts, /*!< Number of ghost points */
                              void    *m,     /*!< MPI object of type #MPIVariables */
                              double  *var    /*!< The array for which to exchange data and fill in ghost points */
                           )
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) m;
  int           d;
  
  int *ip     = mpi->ip;
  int *iproc  = mpi->iproc;
  int *bcflag = mpi->bcperiodic;

  int neighbor_rank[2*ndims], nip[ndims], index[ndims], bounds[ndims], offset[ndims];
  MPI_Request rcvreq[2*ndims], sndreq[2*ndims];
  for (d=0; d<2*ndims; d++) rcvreq[d] = sndreq[d] = MPI_REQUEST_NULL;

  /* each process has 2*ndims neighbors (except at non-periodic physical boundaries)  */
  /* calculate the rank of these neighbors (-1 -> none)                               */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(ip,nip,ndims); 
    if (ip[d] == 0) nip[d] = iproc[d]-1;
    else            nip[d]--;
    if ((ip[d] == 0) && (!bcflag[d])) neighbor_rank[2*d]   = -1;
    else                              neighbor_rank[2*d]   = MPIRank1D(ndims,iproc,nip);
    _ArrayCopy1D_(ip,nip,ndims); 
    if (ip[d] == (iproc[d]-1)) nip[d] = 0;
    else                       nip[d]++;
    if ((ip[d] == (iproc[d]-1)) && (!bcflag[d]))  neighbor_rank[2*d+1] = -1;
    else                                          neighbor_rank[2*d+1] = MPIRank1D(ndims,iproc,nip);
  }

  /* calculate dimensions of each of the send-receive regions */
  double *sendbuf = mpi->sendbuf;
  double *recvbuf = mpi->recvbuf;
  int    stride   = mpi->maxbuf;
  int    bufdim[ndims];
  for (d = 0; d < ndims; d++) {
    bufdim[d] = 1;
    int i;
    for (i = 0; i < ndims; i++) {
      if (i == d) bufdim[d] *= ghosts;
      else        bufdim[d] *= dim[i];
    }
  }

  /* post the receive requests */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Irecv(&recvbuf[2*d*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                mpi->world,&rcvreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Irecv(&recvbuf[(2*d+1)*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1631,
                mpi->world,&rcvreq[2*d+1]);
    }
  }

  /* count number of neighbors and copy data to send buffers */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      _ArraySetValue_(offset,ndims,0);
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
        int p2; _ArrayIndex1D_(ndims,bounds,index,0,p2);
        _ArrayCopy1D_((var+nvars*p1),(sendbuf+2*d*stride+nvars*p2),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
      }
    }
    if (neighbor_rank[2*d+1] != -1) {
      _ArraySetValue_(offset,ndims,0);offset[d] = dim[d]-ghosts;
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
        int p2; _ArrayIndex1D_(ndims,bounds,index,0,p2);
        _ArrayCopy1D_((var+nvars*p1),(sendbuf+(2*d+1)*stride+nvars*p2),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
      }
    }
  }

  /* send the data */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Isend(&sendbuf[2*d*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1631,
                mpi->world,&sndreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Isend(&sendbuf[(2*d+1)*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                mpi->world,&sndreq[2*d+1]);
    }
  }

  /* Wait till data is done received */
  MPI_Waitall(2*ndims,rcvreq,MPI_STATUS_IGNORE);

  /* copy received data to ghost points */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      _ArraySetValue_(offset,ndims,0); offset[d] = -ghosts;
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
        int p2; _ArrayIndex1D_(ndims,bounds,index,0,p2);
        _ArrayCopy1D_((recvbuf+2*d*stride+nvars*p2),(var+nvars*p1),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
      }
    }
    if (neighbor_rank[2*d+1] != -1) {
      _ArraySetValue_(offset,ndims,0); offset[d] = dim[d];
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
        int p2; _ArrayIndex1D_(ndims,bounds,index,0,p2);
        _ArrayCopy1D_((recvbuf+(2*d+1)*stride+nvars*p2),(var+nvars*p1),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
      }
    }
  }
  
  /* Wait till send requests are complete before freeing memory */
  MPI_Waitall(2*ndims,sndreq,MPI_STATUS_IGNORE);

#endif
  return(0);
}

