#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <hypar.h>

int ReadInputs(void *s,void *m)
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ferr    = 0;

  /* set some default values for optional inputs */
  solver->ndims           = 1;
  solver->nvars           = 1;
  solver->dim_global      = NULL;
  solver->dim_local       = NULL;
  mpi->iproc              = NULL;
  mpi->N_IORanks          = 1;
  solver->dt              = 0.0;
  solver->n_iter          = 0;
  solver->restart_iter    = 0;
  solver->screen_op_iter  = 1;
  solver->file_op_iter    = 1000;
  solver->write_residual  = 0;
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
	    	  else if   (!strcmp(word, "restart_iter"       ))  ferr = fscanf(in,"%d",&solver->restart_iter       );
   			  else if   (!strcmp(word, "time_scheme"        ))  ferr = fscanf(in,"%s",solver->time_scheme         );
   		  	else if   (!strcmp(word, "time_scheme_type"   ))  ferr = fscanf(in,"%s",solver->time_scheme_type    );
   		  	else if   (!strcmp(word, "hyp_space_scheme"   ))  ferr = fscanf(in,"%s",solver->spatial_scheme_hyp  );
   		  	else if   (!strcmp(word, "hyp_flux_split"     ))  ferr = fscanf(in,"%s",solver->SplitHyperbolicFlux );
   		  	else if   (!strcmp(word, "hyp_interp_type"    ))  ferr = fscanf(in,"%s",solver->interp_type         );
   		  	else if   (!strcmp(word, "par_space_type"     ))  ferr = fscanf(in,"%s",solver->spatial_type_par    );
   		  	else if   (!strcmp(word, "par_space_scheme"   ))  ferr = fscanf(in,"%s",solver->spatial_scheme_par  );
   		  	else if   (!strcmp(word, "dt"                 ))  ferr = fscanf(in,"%lf",&solver->dt                );
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
	    printf("\tNo. of iter.                               : %d\n"     ,solver->n_iter              );
	    printf("\tRestart iteration                          : %d\n"     ,solver->restart_iter        );
#ifdef with_petsc
      if (solver->use_petscTS)
        printf("\tTime integration scheme                    : PETSc \n"                            );
      else
        printf("\tTime integration scheme                    : %s (%s)\n",
               solver->time_scheme,solver->time_scheme_type                                         );
#else
      printf("\tTime integration scheme                    : %s (%s)\n",
             solver->time_scheme,solver->time_scheme_type                                           );
#endif
      printf("\tSpatial discretization scheme (hyperbolic) : %s\n"     ,solver->spatial_scheme_hyp  );
      printf("\tSplit hyperbolic flux term?                : %s\n"     ,solver->SplitHyperbolicFlux );
      printf("\tInterpolation type for hyperbolic term     : %s\n"     ,solver->interp_type         );
      printf("\tSpatial discretization type   (parabolic ) : %s\n"     ,solver->spatial_type_par    );
      printf("\tSpatial discretization scheme (parabolic ) : %s\n"     ,solver->spatial_scheme_par  );
    	printf("\tTime Step                                  : %E\n"     ,solver->dt                  );
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

  IERR MPIBroadcast_double(&solver->dt,1,0,&mpi->world); CHECKERR(ierr);
#endif

  return(0);
}
