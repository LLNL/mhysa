/*! @file tridiagLUInit.c
    @brief Initialize the tridiagLU solver
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string.h>
#include <tridiagLU.h>
#ifndef serial
#include <mpi.h>
#endif

/*!
Initialize the tridiagLU solver by setting the various parameters in
#TridiagLU to their default values. If the file \a lusolver.inp is 
available, read it and set the parameters.

The file \b lusolver.inp must be in the following format:\n

        begin
            <keyword>   <value>
            <keyword>   <value>
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

where the list of keywords are:\n
    Keyword name       | Type         | Variable                      | Default value
    ------------------ | ------------ | ----------------------------- | ------------------------
    evaluate_norm      | int          | #TridiagLU::evaluate_norm     | 1
    maxiter            | int          | #TridiagLU::maxiter           | 10
    atol               | double       | #TridiagLU::atol              | 1e-12
    rtol               | double       | #TridiagLU::rtol              | 1e-10
    verbose            | int          | #TridiagLU::verbose           | 0
    reducedsolvetype   | char[]       | #TridiagLU::reducedsolvetype  | #_TRIDIAG_JACOBI_

*/
int tridiagLUInit(
                    void *r,  /*!< Object of type TridiagLU */
                    void *c   /*!< MPI communicator */
                 ) 
{
  TridiagLU *t = (TridiagLU*) r;
  int       rank,ierr;
#ifdef serial
  rank  = 0;
#else
  MPI_Comm  *comm = (MPI_Comm*) c;
  if (!comm) rank = 0;
  else MPI_Comm_rank(*comm,&rank);
#endif

  /* default values */
  strcpy(t->reducedsolvetype,_TRIDIAG_JACOBI_);
  t->evaluate_norm = 1;
  t->maxiter       = 10;
  t->atol          = 1e-12;
  t->rtol          = 1e-10;
  t->verbose       = 0;

  /* read from file, if available */
  if (!rank) {
    FILE *in;
    in = fopen("lusolver.inp","r");
    if (!in) {
      printf("tridiagLUInit: File \"lusolver.inp\" not found. Using default values.\n");
    } else {
      char word[100];
      ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
      if (!strcmp(word, "begin")) {
	      while (strcmp(word, "end")) {
		      ierr = fscanf(in,"%s",word); if (ierr != 1) return(1);
     			if      (!strcmp(word, "evaluate_norm"   ))	ierr = fscanf(in,"%d" ,&t->evaluate_norm  );
     			else if (!strcmp(word, "maxiter"         ))	ierr = fscanf(in,"%d" ,&t->maxiter        );
   	  		else if (!strcmp(word, "atol"            ))	ierr = fscanf(in,"%lf",&t->atol           );
   		  	else if (!strcmp(word, "rtol"            ))	ierr = fscanf(in,"%lf",&t->rtol           );
   			  else if (!strcmp(word, "verbose"         ))	ierr = fscanf(in,"%d" ,&t->verbose        );
   			  else if (!strcmp(word, "reducedsolvetype"))	ierr = fscanf(in,"%s" ,t->reducedsolvetype);
          else if (strcmp(word,"end")) {
            char useless[100];
            ierr = fscanf(in,"%s",useless);
            printf("Warning: keyword %s in file \"lusolver.inp\" with value %s not recognized or extraneous. Ignoring.\n",word,useless);
          }
          if (ierr != 1) return(1);
        }
      } else {
   		  fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
        return(1);
      }
      fclose(in);
    }
  }

  /* broadcast to all processes */
#ifndef serial
  if (comm) {
    MPI_Bcast(t->reducedsolvetype,50,MPI_CHAR,0,*comm);
    MPI_Bcast(&t->evaluate_norm,1,MPI_INT,0,*comm);
    MPI_Bcast(&t->maxiter,1,MPI_INT,0,*comm);
    MPI_Bcast(&t->verbose,1,MPI_INT,0,*comm);
    MPI_Bcast(&t->atol,1,MPI_DOUBLE,0,*comm);
    MPI_Bcast(&t->rtol,1,MPI_DOUBLE,0,*comm);
  }
#endif

  return(0);
}

