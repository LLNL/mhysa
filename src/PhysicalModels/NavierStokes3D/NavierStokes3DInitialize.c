/*! @file NavierStokes3DInitialize.c
    @author Debojyoti Ghosh
    @brief Initialization of the physics-related variables and function pointers for the 3D Navier-Stokes system
*/
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

double NavierStokes3DComputeCFL        (void*,void*,double,double);
int    NavierStokes3DFlux              (double*,double*,int,void*,double);
int    NavierStokes3DRoeAverage        (double*,double*,double*,void*);
int    NavierStokes3DParabolicFunction (double*,double*,void*,void*,double);
int    NavierStokes3DSource            (double*,double*,void*,void*,double);

int    NavierStokes3DLeftEigenvectors  (double*,double*,void*,int);
int    NavierStokes3DRightEigenvectors (double*,double*,void*,int);

int    NavierStokes3DUpwindRoe         (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes3DUpwindRF          (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes3DUpwindLLF         (double*,double*,double*,double*,double*,double*,int,void*,double);
int    NavierStokes3DUpwindRusanov     (double*,double*,double*,double*,double*,double*,int,void*,double);

int    NavierStokes3DIBAdiabatic  (void*,void*,double*,double);
int    NavierStokes3DIBIsothermal (void*,void*,double*,double);

int    NavierStokes3DPreStep           (double*,void*,void*,double);
int    NavierStokes3DIBForces          (void*,void*);

/*! Initialize the 3D Navier-Stokes (#NavierStokes3D) module:
    Sets the default parameters, read in and set physics-related parameters,
    and set the physics-related function pointers in #HyPar.

    This file reads the file "physics.inp" that must have the following format:

        begin
            <keyword>   <value>
            <keyword>   <value>
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords are:

    Keyword name       | Type                 | Variable                                                                | Default value
    ------------------ | -------------------- | ----------------------------------------------------------------------- | ------------------------
    nspecies           | int                  | #NavierStokes3D::n_species                                              | 1
    nvibeng            | int                  | #NavierStokes3D::n_vibeng                                               | 0
    gamma              | double               | #NavierStokes3D::gamma                                                  | 1.4
    Pr                 | double               | #NavierStokes3D::Pr                                                     | 0.72
    Re                 | double               | #NavierStokes3D::Re                                                     | -1  
    Minf               | double               | #NavierStokes3D::Minf                                                   | 1.0 
    upwinding          | char[]               | #NavierStokes3D::upw_choice                                             | "rusanov" (#_RUSANOV_)
    ib_wall_type       | char[]               | #NavierStokes3D::ib_wall_type                                           | "adiabatic" (#_IB_ADIABATIC_)

    + if "ib_wall_type" (#NavierStokes3D::ib_wall_type) is specified as "isothermal", 
      it should be followed by the wall temperature (##NavierStokes3D::T_ib_wall), i.e,

        begin
            ...
            ib_wall_type  isothermal 1.0
            ...
        end

    \b Note: "physics.inp" is \b optional; if absent, default values will be used.
*/
int NavierStokes3DInitialize(
                              void *s, /*!< Solver object of type #HyPar */
                              void *m  /*!< MPI object of type #MPIVariables */
                            )
{
  HyPar           *solver  = (HyPar*)          s;
  MPIVariables    *mpi     = (MPIVariables*)   m; 
  NavierStokes3D  *physics = (NavierStokes3D*) solver->physics;
  int             ferr     = 0;

  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in NavierStokes3DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->n_species  = 1;
  physics->n_vibeng   = 0;
  physics->gamma      = 1.4; 
  physics->Pr         = 0.72;
  physics->Re         = -1;
  physics->Minf       = 1.0;
  physics->C1         = 1.458e-6;
  physics->C2         = 110.4;
  physics->T_ib_wall = -DBL_MAX;
  strcpy(physics->upw_choice,_RUSANOV_);
  strcpy(physics->ib_write_surface_data,"yes");
  strcpy(physics->ib_wall_type,"adiabatic");

  /* reading physical model specific inputs */
  if (!mpi->rank) {
    FILE *in;
    printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word);                      if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word);                  if (ferr != 1) return(1);
          if (!strcmp(word, "gamma")) { 
            ferr = fscanf(in,"%lf",&physics->gamma);    if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upw_choice); if (ferr != 1) return(1);
          } else if (!strcmp(word, "nspecies")) { 
            ferr = fscanf(in,"%d",&physics->n_species); if (ferr != 1) return(1);
          } else if (!strcmp(word, "nvibeng")) { 
            ferr = fscanf(in,"%d",&physics->n_vibeng);  if (ferr != 1) return(1);
          } else if (!strcmp(word,"Pr")) {
            ferr = fscanf(in,"%lf",&physics->Pr);       if (ferr != 1) return(1);
          } else if (!strcmp(word,"Re")) {
            ferr = fscanf(in,"%lf",&physics->Re);       if (ferr != 1) return(1);
          } else if (!strcmp(word,"Minf")) {
            ferr = fscanf(in,"%lf",&physics->Minf);     if (ferr != 1) return(1);
          } else if (!strcmp(word,"ib_surface_data")) {
            ferr = fscanf(in,"%s",physics->ib_write_surface_data); if (ferr != 1) return(1);
          } else if (!strcmp(word,"ib_wall_type")) {
            ferr = fscanf(in,"%s",physics->ib_wall_type); 
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (!strcmp(physics->ib_wall_type,_IB_ISOTHERMAL_)) {
              ferr = fscanf(in,"%lf",&physics->T_ib_wall); 
              if (ferr != 1) {
                fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
                return 1;
              }
            }
            if (!solver->flag_ib) {
              printf("Warning: in NavierStokes3DInitialize().\n");
              printf("Warning: no immersed body present; specification of ib_wall_type unnecessary.\n");
            }
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",word,useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
	    } else {
    	  fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
        return(1);
	    }
    }
    fclose(in);
  }

  IERR MPIBroadcast_integer   (&physics->n_species            ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer   (&physics->n_vibeng             ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->upw_choice            ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->ib_write_surface_data ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->ib_wall_type          ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->gamma                ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Pr                   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Re                   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Minf                 ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->T_ib_wall            ,1                ,0,&mpi->world); CHECKERR(ierr);

  if (solver->nvars != (physics->n_species + physics->n_vibeng + 4)) {
    fprintf(  stderr,"Error in NavierStokes3DInitialize(): nvars has to be %d.\n",
              (physics->n_species + physics->n_vibeng + 4 ) );
    return(1);
  }
  physics->nvars = solver->nvars;

  /* if file output is disabled in HyPar, respect that */
  if (!strcmp(solver->op_file_format,"none")) {
    if (!strcmp(physics->ib_write_surface_data,"yes") && solver->flag_ib) {
      if (!mpi->rank) {
        printf("Warning from NavierStokes3DInitialize(): solver->op_file_format is set to \"none\", thus ");
        printf("setting physics->ib_write_surface_data to \"no\" (no solution files will be written).\n");
      }
    }
    strcpy(physics->ib_write_surface_data,"no");
  }

  /* Scaling Re by M_inf */
  physics->Re /= physics->Minf;

  /* check that solver has the correct choice of diffusion formulation, if viscous flow */
  if (strcmp(solver->spatial_type_par,_NC_2STAGE_) && (physics->Re > 0)) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in NavierStokes3DInitialize(): Parabolic term spatial discretization must be \"%s\"\n",_NC_2STAGE_);
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->PreStep    = NavierStokes3DPreStep;
  solver->ComputeCFL = NavierStokes3DComputeCFL;
  solver->FFunction  = NavierStokes3DFlux;
  solver->SFunction  = NavierStokes3DSource;
  solver->AveragingFunction = NavierStokes3DRoeAverage;
  if ((physics->n_species == 1) || (physics->n_vibeng == 1)) {
    solver->GetLeftEigenvectors   = NavierStokes3DLeftEigenvectors;
    solver->GetRightEigenvectors  = NavierStokes3DRightEigenvectors;
  }

  if (solver->flag_ib) {
    if (!strcmp(physics->ib_wall_type,_IB_ADIABATIC_)) {
      solver->IBFunction = NavierStokes3DIBAdiabatic;
    } else if (!strcmp(physics->ib_wall_type,_IB_ISOTHERMAL_)) {
      solver->IBFunction = NavierStokes3DIBIsothermal;
    } else {
      fprintf(stderr, "Error in NavierStokes3DInitialize()\n");
      fprintf(stderr, "  invalid value for IB wall type (%s).\n",
              physics->ib_wall_type );
    }
    if (!strcmp(physics->ib_write_surface_data,"yes")) {
      solver->PhysicsOutput     = NavierStokes3DIBForces;
    }
  }

  if      (!strcmp(physics->upw_choice,_ROE_)) {
    if ((physics->n_species > 1) || (physics->n_vibeng > 1)) {
      if (!mpi->rank) {
        fprintf(stderr,"Error in NavierStokes3DInitialize():\n");
        fprintf(stderr,"  %s is implemented for single-species only.\n",
                physics->upw_choice);
        fprintf(stderr,"  Please use %s for multispecies.\n",_RUSANOV_);
      }
      return 1;
    }
    solver->Upwind = NavierStokes3DUpwindRoe;
  } else if (!strcmp(physics->upw_choice,_RF_CHAR_)) {
    if ((physics->n_species > 1) || (physics->n_vibeng > 1)) {
      if (!mpi->rank) {
        fprintf(stderr,"Error in NavierStokes3DInitialize():\n");
        fprintf(stderr,"  %s is implemented for single-species only.\n",
                physics->upw_choice);
        fprintf(stderr,"  Please use %s for multispecies.\n",_RUSANOV_);
      }
      return 1;
    }
    solver->Upwind = NavierStokes3DUpwindRF;
  } else if (!strcmp(physics->upw_choice,_LLF_CHAR_)) {
    if ((physics->n_species > 1) || (physics->n_vibeng > 1)) {
      if (!mpi->rank) {
        fprintf(stderr,"Error in NavierStokes3DInitialize():\n");
        fprintf(stderr,"  %s is implemented for single-species only.\n",
                physics->upw_choice);
        fprintf(stderr,"  Please use %s for multispecies.\n",_RUSANOV_);
      }
      return 1;
    }
    solver->Upwind = NavierStokes3DUpwindLLF;
  } else if (!strcmp(physics->upw_choice,_RUSANOV_)) {
    solver->Upwind = NavierStokes3DUpwindRusanov;
  } else {
    if (!mpi->rank) {
      fprintf(stderr,"Error in NavierStokes3DInitialize(): %s is not a valid upwinding scheme. ",
              physics->upw_choice);
      fprintf(stderr,"Choices are %s, %s, %s, and %s.\n",_ROE_,_RF_CHAR_,_LLF_CHAR_,_RUSANOV_);
    }
    return(1);
  }

  /* finally, hijack the main solver's dissipation function pointer
   * to this model's own function, since it's difficult to express 
   * the dissipation terms in the general form                      */
  solver->ParabolicFunction = NavierStokes3DParabolicFunction;

  solver->nspecies = physics->n_species;
  return(0);
}
