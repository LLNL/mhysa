/*! @file ReadInputs.c
    @author Debojyoti Ghosh
    @brief Read the input parameters from \b solver.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <timeintegration.h>
#include <mpivars.h>
#include <hypar.h>

/*! Read the simulation inputs from the file \b solver.inp. 
    Rank 0 reads in the inputs and broadcasts them to all the
    processors.\n\n
    The format of \b solver.inp is as follows:\n

        begin
            <keyword>   <value>
            <keyword>   <value>
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords and their type are:\n
    Keyword name       | Type         | Variable                      | Default value
    ------------------ | ------------ | ----------------------------- | ------------------------------------------
    ndims              | int          | #HyPar::ndims                 | 1
    nvars              | int          | #HyPar::nvars                 | 1
    size               | int[ndims]   | #HyPar::dim_global            | must be specified
    iproc              | int[ndims]   | #MPIVariables::iproc          | must be specified
    ghost              | int          | #HyPar::ghosts                | 1
    n_iter             | int          | #HyPar::n_iter                | 0
    restart_iter       | int          | #HyPar::restart_iter          | 0
    time_scheme        | char[]       | #HyPar::time_scheme           | euler
    time_scheme_type   | char[]       | #HyPar::time_scheme_type      | none
    hyp_space_scheme   | char[]       | #HyPar::spatial_scheme_hyp    | 1
    hyp_flux_split     | char[]       | #HyPar::SplitHyperbolicFlux   | no
    hyp_interp_type    | char[]       | #HyPar::interp_type           | characteristic
    par_space_type     | char[]       | #HyPar::spatial_type_par      | nonconservative-1stage
    par_space_scheme   | char[]       | #HyPar::spatial_scheme_par    | 2
    dt                 | double       | #HyPar::dt                    | 0.0
    conservation_check | char[]       | #HyPar::ConservationCheck     | no
    screen_op_iter     | int          | #HyPar::screen_op_iter        | 1
    file_op_iter       | int          | #HyPar::file_op_iter          | 1000
    op_file_format     | char[]       | #HyPar::op_file_format        | text
    ip_file_type       | char[]       | #HyPar::ip_file_type          | ascii
    input_mode         | char[]       | #HyPar::input_mode            | serial
    output_mode        | char[]       | #HyPar::output_mode           | serial
    op_overwrite       | char[]       | #HyPar::op_overwrite          | no
    model              | char[]       | #HyPar::model                 | must be specified
    immersed_body      | char[]       | #HyPar::ib_filename           | "none"

    \b Notes:
    + "ndims" \b must be specified \b before "size" and "iproc".
    + if "input_mode" or "output_mode" are set to "parallel" or "mpi-io",
      the number of I/O ranks must be specified right after as an integer.
      For example:

          begin
              ...
              input_mode  parallel 4
              ...
          end

      This means that 4 MPI ranks will participate in file I/O (assuming
      total MPI ranks is more than 4) (see ReadArrayParallel(), 
      WriteArrayParallel(), ReadArrayMPI_IO() ).
      - The number of I/O ranks specified for "input_mode" and "output_mode"
        \b must \b be \b same. Otherwise, the value for the one specified last
        will be used.
      - The number of I/O ranks must be such that the total number of MPI ranks
        is an integer multiple. Otherwise, the code will use only 1 I/O rank.
    + If any of the keywords are not present, the default value is used, except
      the ones whose default values say "must be specified". Thus, keywords that
      are not required for a particular simulation may be left out of the 
      solver.inp input file. For example, 
      - a #Euler1D simulation does not need "par_space_type" or "par_space_scheme"
        because it does not have a parabolic term.
      - unless a conservation check is required, "conservation_check" can be left
        out and the code will not check for conservation.
      - "immersed_body" need not be specified if there are no immersed bodies present.
         \b NOTE: However, if it is specified, and a file of that filename does not
         exist, it will result in an error.
*/
int ReadInputs(
                void *s,  /*!< Solver object of type #HyPar */
                void *m   /*!< MPI object of type #MPIVariables */
              )
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ferr    = 0;

  /* set some default values for optional inputs */
  solver->ndims           = 1;
  solver->nvars           = 1;
  solver->ghosts          = 1;
  solver->dim_global      = NULL;
  solver->dim_local       = NULL;
  mpi->iproc              = NULL;
  mpi->N_IORanks          = 1;
  solver->dt              = -1;
  solver->cfl             = -1;
  solver->n_iter          = -1;
  solver->t_final         = -1;
  solver->restart_iter    = 0;
  solver->restart_time    = 0;
  solver->screen_op_iter  = 1;
  solver->file_op_iter    = 1000;
  solver->write_residual  = 0;
  solver->flag_ib         = 0;
  strcpy(solver->time_scheme        ,"euler"         );
  strcpy(solver->time_scheme_type   ," "             );
  strcpy(solver->spatial_scheme_hyp ,"1"             );
  strcpy(solver->spatial_type_par   ,_NC_1STAGE_     );
  strcpy(solver->spatial_scheme_par ,"2"             );
  strcpy(solver->interp_type        ,"characteristic");
  strcpy(solver->ip_file_type       ,"ascii"         );
  strcpy(solver->input_mode         ,"serial"        );
  strcpy(solver->output_mode        ,"serial"        );
  strcpy(solver->op_file_format     ,"text"          );
  strcpy(solver->op_overwrite       ,"no"            );
  strcpy(solver->model              ,"none"          );
  strcpy(solver->ConservationCheck  ,"no"            );
  strcpy(solver->SplitHyperbolicFlux,"no"            );
  strcpy(solver->ib_filename        ,"none"          );

  /* reading solver inputs */
  if (!mpi->rank) {
    FILE *in;
    printf("Reading solver inputs from file \"solver.inp\".\n");
    in = fopen("solver.inp","r");
    if (!in) {
      fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
      return(1);
    } else {
	    char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "ndims")) {
            ferr = fscanf(in,"%d",&solver->ndims); if (ferr != 1) return(1);
            solver->dim_global = (int*) calloc (solver->ndims,sizeof(int));
            mpi->iproc         = (int*) calloc (solver->ndims,sizeof(int));
          } else if (!strcmp(word, "nvars")) ferr = fscanf(in,"%d",&solver->nvars);
  			  else if   (!strcmp(word, "size")) {
            int i;
            if (!solver->dim_global) {
              fprintf(stderr,"Error in ReadInputs(): dim_global not allocated.\n");
              fprintf(stderr,"Please specify ndims before dimensions.\n"         );
              return(1);
            } else {
              for (i=0; i<solver->ndims; i++) ferr = fscanf(in,"%d",&solver->dim_global[i]);
              if (ferr != 1) return(1);
            }
          } else if (!strcmp(word, "iproc")) {
            int i;
            if (!mpi->iproc) {
              fprintf(stderr,"Error in ReadInputs(): iproc not allocated.\n");
              fprintf(stderr,"Please specify ndims before iproc.\n"         );
              return(1);
            } else {
              for (i=0; i<solver->ndims; i++) ferr = fscanf(in,"%d",&mpi->iproc[i]);
              if (ferr != 1) return(1);
            }
  			  } else if (!strcmp(word, "ghost"              ))	ferr = fscanf(in,"%d",&solver->ghosts             );
	    	  else if   (!strcmp(word, "n_iter"             ))  ferr = fscanf(in,"%d",&solver->n_iter             );
	    	  else if   (!strcmp(word, "t_final"            ))  ferr = fscanf(in,"%lf",&solver->t_final           );
	    	  else if   (!strcmp(word, "restart_iter"       ))  ferr = fscanf(in,"%d",&solver->restart_iter       );
	    	  else if   (!strcmp(word, "restart_time"       ))  ferr = fscanf(in,"%lf",&solver->restart_time      );
   			  else if   (!strcmp(word, "time_scheme"        ))  ferr = fscanf(in,"%s",solver->time_scheme         );
   		  	else if   (!strcmp(word, "time_scheme_type"   ))  ferr = fscanf(in,"%s",solver->time_scheme_type    );
   		  	else if   (!strcmp(word, "hyp_space_scheme"   ))  ferr = fscanf(in,"%s",solver->spatial_scheme_hyp  );
   		  	else if   (!strcmp(word, "hyp_flux_split"     ))  ferr = fscanf(in,"%s",solver->SplitHyperbolicFlux );
   		  	else if   (!strcmp(word, "hyp_interp_type"    ))  ferr = fscanf(in,"%s",solver->interp_type         );
   		  	else if   (!strcmp(word, "par_space_type"     ))  ferr = fscanf(in,"%s",solver->spatial_type_par    );
   		  	else if   (!strcmp(word, "par_space_scheme"   ))  ferr = fscanf(in,"%s",solver->spatial_scheme_par  );
   		  	else if   (!strcmp(word, "dt"                 ))  ferr = fscanf(in,"%lf",&solver->dt                );
   		  	else if   (!strcmp(word, "cfl"                ))  ferr = fscanf(in,"%lf",&solver->cfl               );
   		  	else if   (!strcmp(word, "conservation_check" ))  ferr = fscanf(in,"%s",solver->ConservationCheck   );
   		  	else if   (!strcmp(word, "screen_op_iter"     ))  ferr = fscanf(in,"%d",&solver->screen_op_iter     );
   		  	else if   (!strcmp(word, "file_op_iter"       ))  ferr = fscanf(in,"%d",&solver->file_op_iter       );
   		  	else if   (!strcmp(word, "op_file_format"     ))  ferr = fscanf(in,"%s",solver->op_file_format      );
   		  	else if   (!strcmp(word, "ip_file_type"       ))  ferr = fscanf(in,"%s",solver->ip_file_type        );
   		  	else if   (!strcmp(word, "input_mode"         ))  {
            ferr = fscanf(in,"%s",solver->input_mode);
            if (strcmp(solver->input_mode,"serial"))        ferr = fscanf(in,"%d",&mpi->N_IORanks);
   		  	} else if (!strcmp(word, "output_mode"        ))  {
            ferr = fscanf(in,"%s",solver->output_mode);
            if (strcmp(solver->output_mode,"serial"))       ferr = fscanf(in,"%d",&mpi->N_IORanks);
          } else if   (!strcmp(word, "op_overwrite"     ))  ferr = fscanf(in,"%s",solver->op_overwrite      );
   		  	else if   (!strcmp(word, "model"              ))  ferr = fscanf(in,"%s",solver->model             );
   		  	else if   (!strcmp(word, "immersed_body"      ))  ferr = fscanf(in,"%s",solver->ib_filename       );
          else if   ( strcmp(word, "end"                )) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            printf("Warning: keyword %s in file \"solver.inp\" with value %s not recognized or extraneous. Ignoring.\n",
                    word,useless);
          }
          if (ferr != 1) return(1);
        }
      } else {
   		  fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
        return(1);
      }
      fclose(in);

      if (solver->screen_op_iter <= 0)  solver->screen_op_iter = 1;
      if (solver->file_op_iter <= 0)    solver->file_op_iter   = solver->n_iter;

      if ((solver->ndims != 3) && (strcmp(solver->ib_filename,"none"))) {
        printf("Warning: immersed boundaries not implemented for ndims = %d. ",solver->ndims);
        printf("Ignoring input for \"immersed_body\" (%s).\n",solver->ib_filename);
        strcpy(solver->ib_filename,"none");
      }
      solver->flag_ib = strcmp(solver->ib_filename,"none");

      /* Print to screen the inputs read */
      int i;
	    printf("\tNo. of dimensions                          : %d\n",solver->ndims);
	    printf("\tNo. of variables                           : %d\n",solver->nvars);
	    printf("\tDomain size                                : ");
      for (i=0; i<solver->ndims; i++) printf ("%d ",solver->dim_global[i]);
      printf("\n");
#ifndef serial
	    printf("\tProcesses along each dimension             : ");
      for (i=0; i<solver->ndims; i++) printf ("%d ",mpi->iproc[i]);
      printf("\n");
#endif
	    printf("\tNo. of ghosts pts                          : %d\n"     ,solver->ghosts              );

	    if (solver->n_iter >= 0)  printf("\tNo. of iter.                               : %d\n",solver->n_iter);
	    if (solver->t_final>= 0)  printf("\tFinal time                                 : %f\n",solver->t_final);

	    printf("\tRestart iteration                          : %d\n"     ,solver->restart_iter        );
	    printf("\tRestart simulation time                    : %f\n"     ,solver->restart_time        );
#ifdef with_petsc
      if (solver->use_petscTS)
        printf("\tTime integration scheme                    : PETSc \n"                            );
      else {
        printf("\tTime integration scheme                    : %s ",solver->time_scheme             );
        if (strcmp(solver->time_scheme,_FORWARD_EULER_)) {
          printf("(%s)",solver->time_scheme_type                                                    );
        }
        printf("\n");
      }
#else
      printf("\tTime integration scheme                    : %s ",solver->time_scheme               );
      if (strcmp(solver->time_scheme,_FORWARD_EULER_)) {
        printf("(%s)",solver->time_scheme_type                                                      );
      }
      printf("\n");
#endif
      printf("\tSpatial discretization scheme (hyperbolic) : %s\n"     ,solver->spatial_scheme_hyp  );
      printf("\tSplit hyperbolic flux term?                : %s\n"     ,solver->SplitHyperbolicFlux );
      printf("\tInterpolation type for hyperbolic term     : %s\n"     ,solver->interp_type         );
      printf("\tSpatial discretization type   (parabolic ) : %s\n"     ,solver->spatial_type_par    );
      printf("\tSpatial discretization scheme (parabolic ) : %s\n"     ,solver->spatial_scheme_par  );

    	if (solver->dt  >= 0) printf("\tTime Step                                  : %E\n",solver->dt );
      if (solver->cfl >= 0) printf("\tCFL                                        : %E\n",solver->cfl);

    	printf("\tCheck for conservation                     : %s\n"     ,solver->ConservationCheck   );
      printf("\tScreen output iterations                   : %d\n"     ,solver->screen_op_iter      );
      printf("\tFile output iterations                     : %d\n"     ,solver->file_op_iter        );
      printf("\tInitial solution file type                 : %s\n"     ,solver->ip_file_type        );
      printf("\tInitial solution read mode                 : %s"       ,solver->input_mode          );
      if (strcmp(solver->input_mode,"serial"))    printf("\t[%d file IO rank(s)]\n",mpi->N_IORanks  );
      else                                        printf("\n");
      printf("\tSolution file write mode                   : %s"       ,solver->output_mode         );
      if (strcmp(solver->output_mode,"serial"))   printf("\t[%d file IO rank(s)]\n",mpi->N_IORanks  );
      else                                        printf("\n");
      printf("\tSolution file format                       : %s\n"     ,solver->op_file_format      );
      printf("\tOverwrite solution file                    : %s\n"     ,solver->op_overwrite        );
      printf("\tPhysical model                             : %s\n"     ,solver->model               );
      if (solver->flag_ib) {
        printf("\tImmersed Body                              : %s\n"     ,solver->ib_filename       );
      }
    }
    /* checks - restart only supported for binary output files */
    if ((solver->restart_iter != 0) && strcmp(solver->op_file_format,"binary")) {
      if (!mpi->rank) fprintf(stderr,"Error in ReadInputs(): Restart is supported only for binary output files.\n");
      return(1);
    }
  }

#ifndef serial
  /* Broadcast the input parameters */
  IERR MPIBroadcast_integer(&solver->ndims,1,0,&mpi->world); CHECKERR(ierr);
  if (mpi->rank) {
    solver->dim_global = (int*) calloc (solver->ndims,sizeof(int));
    mpi->iproc         = (int*) calloc (solver->ndims,sizeof(int));
  }
  IERR MPIBroadcast_integer(&solver->nvars              ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer( solver->dim_global         ,solver->ndims,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer( mpi->iproc                 ,solver->ndims,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&mpi->N_IORanks             ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&solver->ghosts             ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&solver->n_iter             ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&solver->restart_iter       ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&solver->screen_op_iter     ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&solver->file_op_iter       ,1            ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&solver->flag_ib            ,1            ,0,&mpi->world); CHECKERR(ierr);

  IERR MPIBroadcast_character(solver->time_scheme         ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->time_scheme_type    ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->spatial_scheme_hyp  ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->interp_type         ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->spatial_type_par    ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->spatial_scheme_par  ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->ConservationCheck   ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->SplitHyperbolicFlux ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->op_file_format      ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->ip_file_type        ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->input_mode          ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->output_mode         ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->op_overwrite        ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->model               ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(solver->ib_filename         ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);

  IERR MPIBroadcast_double(&solver->dt          ,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double(&solver->cfl         ,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double(&solver->t_final     ,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double(&solver->restart_time,1,0,&mpi->world); CHECKERR(ierr);
#endif

  return(0);
}
