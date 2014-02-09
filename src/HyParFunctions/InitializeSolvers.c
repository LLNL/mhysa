#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <hypar.h>
#include <tridiagLU.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <secondderivative.h>

/* Function declarations */
int WriteBinary                 (int,int,int*,double*,double*,char*,int*);
int WriteText                   (int,int,int*,double*,double*,char*,int*);
int WriteTecplot2D              (int,int,int*,double*,double*,char*,int*);
int WriteTecplot3D              (int,int,int*,double*,double*,char*,int*);
int ApplyBoundaryConditions     (void*,void*,double*);
int HyperbolicFunction          (double*,double*,void*,void*,double);
int ParabolicFunctionNC1Stage   (double*,double*,void*,void*,double);
int ParabolicFunctionCons1Stage (double*,double*,void*,void*,double);
int SourceFunction              (double*,double*,void*,void*,double);

int InitializeSolvers(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  if (!mpi->rank) printf("Initializing solvers.\n");

  solver->ApplyBoundaryConditions = ApplyBoundaryConditions;
  solver->HyperbolicFunction      = HyperbolicFunction;
  solver->SourceFunction          = SourceFunction;

  /* choose the type of parabolic discretization */
  if (!strcmp(solver->spatial_type_par,_NC_1STAGE_)) {
    solver->ParabolicFunction = ParabolicFunctionNC1Stage;
    if (!strcmp(solver->spatial_scheme_par,_SECOND_ORDER_)) {
      solver->SecondDerivativePar = SecondDerivativeSecondOrder; 
    } else {
      fprintf(stderr,"Error: %s is not a supported ",solver->spatial_scheme_par);
      fprintf(stderr,"spatial scheme of type %s for the parabolic terms.\n",
              solver->spatial_type_par);
    }
  } else if (!strcmp(solver->spatial_type_par,_CONS_1STAGE_)) {
    solver->ParabolicFunction = ParabolicFunctionCons1Stage;
    if (!strcmp(solver->spatial_scheme_par,_SECOND_ORDER_CENTRAL_)) {
      solver->InterpolateInterfacesPar = Interp2PrimSecondOrder; 
    } else {
      fprintf(stderr,"Error: %s is not a supported ",solver->spatial_scheme_par);
      fprintf(stderr,"spatial scheme of type %s for the parabolic terms.\n",
              solver->spatial_type_par);
    }
  } else {
    fprintf(stderr,"Error: %s is not a supported ",solver->spatial_type_par);
    fprintf(stderr,"spatial discretization type for the parabolic terms.\n");
    return(1);
  }

  /* Spatial interpolation for hyperbolic term */
  solver->interp    = NULL;
  solver->lusolver  = NULL;
  if (!strcmp(solver->spatial_scheme_hyp,_FIRST_ORDER_UPWIND_)) {
    /* First order upwind scheme */
    if (solver->nvars > 1) {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_))
        solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwindChar;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwind;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    } else {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) 
        solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwind;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwind;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    }
  } else if (!strcmp(solver->spatial_scheme_hyp,_THIRD_ORDER_MUSCL_)) {
    /* Third order MUSCL scheme */
    if (solver->nvars > 1) {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_))
        solver->InterpolateInterfacesHyp = Interp1PrimThirdOrderMUSCLChar;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimThirdOrderMUSCL;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    } else {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) 
        solver->InterpolateInterfacesHyp = Interp1PrimThirdOrderMUSCL;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimThirdOrderMUSCL;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    }
    solver->interp = (MUSCLParameters*) calloc(1,sizeof(MUSCLParameters));
    IERR MUSCLInitialize(solver->interp,mpi); CHECKERR(ierr);
  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_WENO_)) {
    /* Fifth order WENO scheme */
    if (solver->nvars > 1) {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENOChar;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENO;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    } else {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) 
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENO;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENO;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    }
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp); CHECKERR(ierr);
  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
    /* Fifth order CRWENO scheme */
    if (solver->nvars > 1) {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) {
        fprintf(stderr,"Error in InitializeSolvers(): Characteristic CRWENO not yet available.\n");
        return(1);
      } else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCRWENO;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    } else {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) 
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCRWENO;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCRWENO;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    }
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp); CHECKERR(ierr);
    solver->lusolver = (TridiagLU*) calloc (1,sizeof(TridiagLU));
    IERR tridiagLUInit(solver->lusolver,&mpi->world);CHECKERR(ierr);
  } else {
    fprintf(stderr,"Error: %s is a not a supported spatial interpolation scheme.\n",
            solver->spatial_scheme_hyp);
    return(1);
  }

  /* Time integration */
  if (!strcmp(solver->time_scheme,_FORWARD_EULER_)) { 
    solver->TimeIntegrate = TimeForwardEuler;
    solver->msti = NULL;
  } else if (!strcmp(solver->time_scheme,_RK_)) {
    solver->TimeIntegrate = TimeRK;
    solver->msti = (MSTIParameters*) calloc (1,sizeof(MSTIParameters));
    IERR TimeMSTIInitialize(solver->time_scheme,solver->time_scheme_type,
                              solver->msti); CHECKERR(ierr);
  } else {
    fprintf(stderr,"Error: %s is a not a supported time-integration scheme.\n",
            solver->time_scheme);
    return(1);
  }

  /* Solution output function */
  solver->WriteOutput = WriteText;
  if (!strcmp(solver->op_overwrite,"no")) strcpy(solver->op_filename,"op_00000");
  else                                    strcpy(solver->op_filename,"op");
  if (!strcmp(solver->op_file_format,"text")) {
    solver->WriteOutput = WriteText;
    strcat(solver->op_filename,".dat");
  } else if (!strcmp(solver->op_file_format,"tecplot2d")) {
    solver->WriteOutput = WriteTecplot2D;
    strcat(solver->op_filename,".dat");
  } else if (!strcmp(solver->op_file_format,"tecplot3d")) {
    solver->WriteOutput = WriteTecplot3D;
    strcat(solver->op_filename,".dat");
  } else if ((!strcmp(solver->op_file_format,"binary")) || (!strcmp(solver->op_file_format,"bin"))) {
    solver->WriteOutput = WriteBinary;
    strcat(solver->op_filename,".bin");
  } else if (!strcmp(solver->op_file_format,"none")) {
    solver->WriteOutput = NULL;
  } else {
    fprintf(stderr,"Error: %s is not a supported file format.\n",solver->op_file_format);
    return(1);
  }

  return(0);
}
