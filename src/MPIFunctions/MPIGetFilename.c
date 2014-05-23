#ifndef serial
#include <mpi.h>
#endif
#include <string.h>
#include <basic.h>
#include <mpivars.h>

static void GetStringFromInteger(int,char*,int);

void MPIGetFilename(char *root,void *c,char *filename)
{
  char  tail[_MAX_STRING_SIZE_];
  int   rank;

#ifndef serial
  MPI_Comm  comm = *((MPI_Comm*)c);
  MPI_Comm_rank(comm,&rank);
#else
  rank = 0;
#endif

  GetStringFromInteger(rank,tail,4);
  strcat(filename,root);
  strcat(filename,"." );
  strcat(filename,tail);

  return;
}

void GetStringFromInteger(int a,char *A,int width)
{
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0'); 
    a /= 10;
    A[width-1-i] = digit;
  }
  return;
}
