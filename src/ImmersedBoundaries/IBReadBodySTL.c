/*! @file IBReadBodySTL.c
    @author Debojyoti Ghosh
    @brief Reads 3D body surface from STL file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <immersedboundaries.h>

/*! Function to read a 3D surface from a STL file. See 
    https://en.wikipedia.org/wiki/STL_(file_format)
    for details of the STL file format. 
    \b Notes: 
    + This function only reads ASCII STL files.
    + The body read from this file is distributed to
      all MPI ranks.
*/
int IBReadBodySTL(
                    Body3D **body,    /*!< 3D body to be allocated and read from file */
                    char   *filename, /*!< Filename */
                    void   *m,        /*!< MPI object of type #MPIVariables */
                    int    *stat      /*!< Status (0: success; non-0: failure) */
                 )
{
  MPIVariables  *mpi = (MPIVariables*) m;
  int n, ierr;
  *stat = 0;

  if ((*body) != NULL) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in IBReadBodySTL(): pointer to immersed body not NULL.\n");
    }
    *stat = -1;
    return(0);
  }

  /* Rank 0 reads in file */
  if (!mpi->rank) {

    FILE  *in;
    in = fopen(filename,"r");
    if (!in) {
      fprintf(stderr,"Error in IBReadBodySTL(): file %s not found or cannot be opened.\n",
              filename);
      *stat = 1;
    }
    else {

      int   nfacets = 0;
	    char  word[_MAX_STRING_SIZE_];
      /* Count number of facets in STL file */
      ierr = fscanf(in,"%s",word);
      if (!strcmp(word, "solid")) {
        nfacets = 0;
        while (strcmp(word, "endsolid")) {
          ierr = fscanf(in,"%s",word);
          if (!strcmp(word, "facet")) nfacets++;
        }
      }
      fclose(in);

      if (nfacets > 0) {
        (*body) = (Body3D*) calloc (1,sizeof(Body3D));
        (*body)->nfacets = nfacets;
        (*body)->surface = (Facet3D*) calloc (nfacets,sizeof(Facet3D));

        /* Read in STL data */
        in = fopen(filename,"r");
        ierr = fscanf(in,"%s",word);
        n = 0;
        while (strcmp(word, "endsolid")) {
          double t1, t2, t3;
          ierr = fscanf(in,"%s",word);
          if (!strcmp(word, "facet")) {
            ierr = fscanf(in,"%s",word);
            if (strcmp(word, "normal")) {
              fprintf(stderr,"Error in IBReadBodySTL(): Illegal keyword read in %s.\n",filename);
            } else {
              ierr = fscanf(in,"%lf",&t1);
              ierr = fscanf(in,"%lf",&t2);
              ierr = fscanf(in,"%lf",&t3);
              (*body)->surface[n].nx = t1; 
              (*body)->surface[n].ny = t2; 
              (*body)->surface[n].nz = t3;
            }
            ierr = fscanf(in,"%s",word);
            if (strcmp(word, "outer")) {
              fprintf(stderr,"Error in IBReadBodySTL(): Illegal keyword read in %s.\n",filename);
              *stat = 1;
            } else {
              ierr = fscanf(in,"%s",word);
              if (strcmp(word, "loop")) {
                fprintf(stderr,"Error in IBReadBodySTL(): Illegal keyword read in %s.\n",filename);
                *stat = 1;
              } else {
                ierr = fscanf(in,"%s",word);
                if (strcmp(word, "vertex")) {
                  fprintf(stderr,"Error in IBReadBodySTL(): Illegal keyword read in %s.\n",filename);
                  *stat = 1;
                } else {
                  ierr = fscanf(in,"%lf",&t1);
                  ierr = fscanf(in,"%lf",&t2);
                  ierr = fscanf(in,"%lf",&t3);
                  (*body)->surface[n].x1 = t1; 
                  (*body)->surface[n].y1 = t2; 
                  (*body)->surface[n].z1 = t3;
                }
                ierr = fscanf(in,"%s",word);
                if (strcmp(word, "vertex")) {
                  fprintf(stderr,"Error in IBReadBodySTL(): Illegal keyword read in %s.\n",filename);
                  *stat = 1;
                } else {
                  ierr = fscanf(in,"%lf",&t1);
                  ierr = fscanf(in,"%lf",&t2);
                  ierr = fscanf(in,"%lf",&t3);
                  (*body)->surface[n].x2 = t1; 
                  (*body)->surface[n].y2 = t2; 
                  (*body)->surface[n].z2 = t3;
                }
                ierr = fscanf(in,"%s",word);
                if (strcmp(word, "vertex")) {
                  fprintf(stderr,"Error in IBReadBodySTL(): Illegal keyword read in %s.\n",filename);
                  *stat = 1;
                } else {
                  ierr = fscanf(in,"%lf",&t1);
                  ierr = fscanf(in,"%lf",&t2);
                  ierr = fscanf(in,"%lf",&t3);
                  (*body)->surface[n].x3 = t1; 
                  (*body)->surface[n].y3 = t2; 
                  (*body)->surface[n].z3 = t3;
                }
              }
            }
            n++;
          }
        }
        fclose(in);
        /* done reading STL file */
        if (n != nfacets) {
          fprintf(stderr,"Error in IBReadBodySTL(): Inconsistency in number of facets read.\n");
          free((*body)->surface);
          free(*body);
          *stat = 1;
        }
      } else {
        fprintf(stderr,"Error in IBReadBodySTL(): nfacets = 0!!\n");
        *stat = 1;
      }

    }
  }
  MPIBroadcast_integer(stat,1,0,&mpi->world);

  if ((*stat)) return(0);

  /* Distribute the body to all MPI ranks */
  int nfacets;
  if (!mpi->rank) nfacets = (*body)->nfacets;
  MPIBroadcast_integer(&nfacets,1,0,&mpi->world);

  if (mpi->rank) {
    (*body) = (Body3D*) calloc (1,sizeof(Body3D));
    (*body)->nfacets = nfacets;
    (*body)->surface = (Facet3D*) calloc (nfacets,sizeof(Facet3D));
  }

  int bufdim = 12;
  double *buffer = (double*) calloc (bufdim*nfacets,sizeof(double));
  if (!mpi->rank) {
    for (n=0; n<nfacets; n++) {
      buffer[n*bufdim+ 0] = (*body)->surface[n].x1;
      buffer[n*bufdim+ 1] = (*body)->surface[n].x2;
      buffer[n*bufdim+ 2] = (*body)->surface[n].x3;
      buffer[n*bufdim+ 3] = (*body)->surface[n].y1;
      buffer[n*bufdim+ 4] = (*body)->surface[n].y2;
      buffer[n*bufdim+ 5] = (*body)->surface[n].y3;
      buffer[n*bufdim+ 6] = (*body)->surface[n].z1;
      buffer[n*bufdim+ 7] = (*body)->surface[n].z2;
      buffer[n*bufdim+ 8] = (*body)->surface[n].z3;
      buffer[n*bufdim+ 9] = (*body)->surface[n].nx;
      buffer[n*bufdim+10] = (*body)->surface[n].ny;
      buffer[n*bufdim+11] = (*body)->surface[n].nz;
    }
  }
  MPIBroadcast_double(buffer,(nfacets*bufdim),0,&mpi->world);
  if (mpi->rank) {
    for (n=0; n<nfacets; n++) {
      (*body)->surface[n].x1 = buffer[n*bufdim+ 0];
      (*body)->surface[n].x2 = buffer[n*bufdim+ 1];
      (*body)->surface[n].x3 = buffer[n*bufdim+ 2];
      (*body)->surface[n].y1 = buffer[n*bufdim+ 3];
      (*body)->surface[n].y2 = buffer[n*bufdim+ 4];
      (*body)->surface[n].y3 = buffer[n*bufdim+ 5];
      (*body)->surface[n].z1 = buffer[n*bufdim+ 6];
      (*body)->surface[n].z2 = buffer[n*bufdim+ 7];
      (*body)->surface[n].z3 = buffer[n*bufdim+ 8];
      (*body)->surface[n].nx = buffer[n*bufdim+ 9];
      (*body)->surface[n].ny = buffer[n*bufdim+10];
      (*body)->surface[n].nz = buffer[n*bufdim+11];
    }
  }
  free(buffer);

  return(0);
}
