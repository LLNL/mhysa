/*! @file MPIGetFilename.c
    @brief Get filename indexed by MPI rank
    @author Debojyoti Ghosh
*/

#ifndef serial
#include <mpi.h>
#endif
#include <string.h>
#include <basic.h>
#include <mpivars.h>

static void GetStringFromInteger(int,char*,int);

/*! 
    Get a string representing a filename indexed by the MPI rank: 
    \a filename = \a root.index, where \a index is the string
    corresponding to the MPI rank.
*/
void MPIGetFilename(
                      char  *root,      /*!< filename root */
                      void  *c,         /*!< MPI communicator */
                      char  *filename   /*!< filename */
                   )
{
  char  tail[_MAX_STRING_SIZE_]="";
  int   rank;

#ifndef serial
  MPI_Comm  comm = *((MPI_Comm*)c);
  MPI_Comm_rank(comm,&rank);
#else
  rank = 0;
#endif

  GetStringFromInteger(rank,tail,4);
  strcpy(filename,"");
  strcat(filename,root);
  strcat(filename,"." );
  strcat(filename,tail);

  return;
}

/*!
  Get a string corresponding to an integer, i.e. 41 gives "00041" if
  \a width is 5, or "41" if \a width is 2, or "1" if \a width is 1.
*/
void GetStringFromInteger(
                            int   a,    /*!< the integer to convert to a string */
                            char  *A,   /*!< the string */
                            int   width /*!< desired width of the string */
                         )
{
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0'); 
    a /= 10;
    A[width-1-i] = digit;
  }
  return;
}
