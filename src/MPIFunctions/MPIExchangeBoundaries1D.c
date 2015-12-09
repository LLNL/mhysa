/*! @file MPIExchangeBoundaries1D.c
    @brief Exchange data and fill ghost points for a 1D array
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Exchange the data across MPI ranks and fill in ghost points for a 1D array. In a multidimensional
  simulation, a 1D array is an array of data along one of the spatial dimensions, i.e. its an array
  storing a variable that varies in only one of the spatial dimension. For example, for a
  2D problem on a Cartesian grid (with spatial dimensions x and y), the array of x-coordinates is
  a 1D array along x, and the array of y-coordinates is a 1D array along y. Thus, the size of the 
  1D array is equal to the size of the domain along the spatial dimension corresponding to that array.

  Consider a two-dimensional problem, partitioned on 21 MPI ranks as follows:
  @image html mpi_ranks.png
  @image latex mpi_ranks.eps width=0.9\textwidth
  and consider rank 9. 
  
  If the argument \a dir is specified as 0, and thus we are dealing with a 1D array
  along dimension 0, then
  + Rank 9 will exchange data with ranks 8 and 10, and fill in its ghost points.

  If \a dir is specified as 1, and thus we are dealing with a 1D array along dimension
  1, then
  + Rank 9 will exchange data with ranks 2 and 16, and fill in its ghost points.
*/
int MPIExchangeBoundaries1D(
                              void    *m,     /*!< MPI object of type MPIVariables */
                              double  *x,     /*!< The 1D array for which to exchange data */
                              int     N,      /*!< Size of the array */
                              int     ghosts, /*!< Number of ghost points */
                              int     dir,    /*!< Spatial dimension corresponding to the 1D array */
                              int     ndims   /*!< Number of spatial dimensions in the simulation */
                           )
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) m;
  int           i;
  
  int *ip     = mpi->ip;
  int *iproc  = mpi->iproc;
  int non      = 0; /* number of neighbours */

  int neighbor_rank[2] = {-1,-1};
  int nip[ndims];

  /* each process has 2 neighbors (except at physical boundaries)       */
  /* calculate the rank of these neighbors (-1 -> none)                 */
  _ArrayCopy1D_(ip,nip,ndims);  nip[dir]--;
  if (ip[dir] == 0)             neighbor_rank[0] = -1;
  else                          neighbor_rank[0] = MPIRank1D(ndims,iproc,nip);
  _ArrayCopy1D_(ip,nip,ndims);  nip[dir]++;
  if (ip[dir] == (iproc[dir]-1))neighbor_rank[1] = -1;
  else                          neighbor_rank[1] = MPIRank1D(ndims,iproc,nip);

  /* Allocate send and receive buffers */
  double sendbuf[2][ghosts], recvbuf[2][ghosts];

  /* count number of neighbors and copy data to send buffers */
  non = 0;
  if (neighbor_rank[0] != -1) {
    non++;
    for (i = 0; i < ghosts; i++) sendbuf[0][i] = x[i+ghosts];
  }
  if (neighbor_rank[1] != -1) {
    non++;
    for (i = 0; i < ghosts; i++) sendbuf[1][i] = x[i+N];
  }
  MPI_Request requests[2*non];
  MPI_Status  statuses[2*non];

  /* exchange the data */
  int tick = 0;
  if (neighbor_rank[0]!= -1) {
    MPI_Irecv(recvbuf[0],ghosts,MPI_DOUBLE,neighbor_rank[0],1631,mpi->world,&requests[tick]);
    MPI_Isend(sendbuf[0],ghosts,MPI_DOUBLE,neighbor_rank[0],1631,mpi->world,&requests[tick+non]);
    tick++;
  }
  if (neighbor_rank[1] != -1) {
    MPI_Irecv(recvbuf[1],ghosts,MPI_DOUBLE,neighbor_rank[1],1631,mpi->world,&requests[tick]);
    MPI_Isend(sendbuf[1],ghosts,MPI_DOUBLE,neighbor_rank[1],1631,mpi->world,&requests[tick+non]);
    tick++;
  }

  /* Wait till data transfer is done */
  MPI_Waitall(2*non,requests,statuses);

  /* copy received data to ghost points */
  if (neighbor_rank[0] != -1) for (i = 0; i < ghosts; i++) x[i]          = recvbuf[0][i];
  if (neighbor_rank[1] != -1) for (i = 0; i < ghosts; i++) x[i+N+ghosts] = recvbuf[1][i];
  
#endif
  return(0);
}

