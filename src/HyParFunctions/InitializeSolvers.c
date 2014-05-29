#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <hypar.h>
#include <tridiagLU.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <firstderivative.h>
#include <secondderivative.h>

/* Function declarations */
int  WriteBinary                 (int,int,int*,double*,double*,char*,int*);
int  WriteText                   (int,int,int*,double*,double*,char*,int*);
int  WriteTecplot2D              (int,int,int*,double*,double*,char*,int*);
int  WriteTecplot3D              (int,int,int*,double*,double*,char*,int*);
int  ApplyBoundaryConditions     (void*,void*,double*,double*,int,double);
int  HyperbolicFunction          (double*,double*,void*,void*,double,int);
int  HyperbolicFunction1         (double*,double*,void*,void*,double,int);
int  HyperbolicFunction2         (double*,double*,void*,void*,double,int);
int  ParabolicFunctionNC1Stage   (double*,double*,void*,void*,double);
int  ParabolicFunctionNC2Stage   (double*,double*,void*,void*,double);
int  ParabolicFunctionCons1Stage (double*,double*,void*,void*,double);
int  SourceFunction              (double*,double*,void*,void*,double);
int  VolumeIntegral              (double*,double*,void*,void*);
int  BoundaryIntegral            (void*,void*);
void IncrementFilename           (char*);

int InitializeSolvers(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  if (!mpi->rank) printf("Initializing solvers.\n");

  solver->ApplyBoundaryConditions = ApplyBoundaryConditions;
  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    solver->HyperbolicFunction1     = HyperbolicFunction1;
    solver->HyperbolicFunction2     = HyperbolicFunction2;
    solver->HyperbolicFunction      = NULL;
  } else {
    solver->HyperbolicFunction      = HyperbolicFunction;
    solver->HyperbolicFunction1     = NULL;
    solver->HyperbolicFunction2     = NULL;
  }
  solver->SourceFunction            = SourceFunction;
  solver->VolumeIntegralFunction    = VolumeIntegral;
  solver->BoundaryIntegralFunction  = BoundaryIntegral;

  /* choose the type of parabolic discretization */
  if (!strcmp(solver->spatial_type_par,_NC_1STAGE_)) 
    solver->ParabolicFunction = ParabolicFunctionNC1Stage;
  else if (!strcmp(solver->spatial_type_par,_NC_2STAGE_))
    solver->ParabolicFunction = ParabolicFunctionNC2Stage;
  else if (!strcmp(solver->spatial_type_par,_CONS_1STAGE_))
    solver->ParabolicFunction = ParabolicFunctionCons1Stage;
  else {
    fprintf(stderr,"Error: %s is not a supported ",solver->spatial_type_par);
    fprintf(stderr,"spatial discretization type for the parabolic terms.\n");
    return(1);
  }
  if (!strcmp(solver->spatial_scheme_par,_SECOND_ORDER_CENTRAL_)) {
    solver->SecondDerivativePar      = SecondDerivativeSecondOrderCentral; 
    solver->FirstDerivativePar       = FirstDerivativeSecondOrderCentral; 
    solver->InterpolateInterfacesPar = Interp2PrimSecondOrder; 
  } else if (!strcmp(solver->spatial_scheme_par,_FOURTH_ORDER_CENTRAL_)) {
    solver->SecondDerivativePar      = SecondDerivativeFourthOrderCentral; 
    solver->FirstDerivativePar       = FirstDerivativeFourthOrderCentral; 
    solver->InterpolateInterfacesPar = NULL; /* not yet coded, setting to NULL so that the code crashes */
  } else {
    fprintf(stderr,"Error: %s is not a supported ",solver->spatial_scheme_par);
    fprintf(stderr,"spatial scheme of type %s for the parabolic terms.\n",
            solver->spatial_type_par);
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
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp,solver->interp_type); CHECKERR(ierr);
  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
    /* Fifth order CRWENO scheme */
    if (solver->nvars > 1) {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCRWENOChar;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
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
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp,solver->interp_type); CHECKERR(ierr);
    solver->lusolver = (TridiagLU*) calloc (1,sizeof(TridiagLU));
    IERR tridiagLUInit(solver->lusolver,&mpi->world);CHECKERR(ierr);
  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_)) {
    /* Fifth order HCWENO scheme */
    if (solver->nvars > 1) {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderHCWENOChar;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderHCWENO;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    } else {
      if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) 
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderHCWENO;
      else if (!strcmp(solver->interp_type,_COMPONENTS_))
        solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderHCWENO;
      else {
        fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
        fprintf(stderr,"supported interpolation type.\n");
        return(1);
      }
    }
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp,solver->interp_type); CHECKERR(ierr);
    solver->lusolver = (TridiagLU*) calloc (1,sizeof(TridiagLU));
    IERR tridiagLUInit(solver->lusolver,&mpi->world);CHECKERR(ierr);
  } else {
    fprintf(stderr,"Error: %s is a not a supported spatial interpolation scheme.\n",
            solver->spatial_scheme_hyp);
    return(1);
  }
  solver->SetInterpLimiterVar = InterpSetLimiterVar;

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
  solver->WriteOutput = NULL; /* default - no output */
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
  if ((!strcmp(solver->op_overwrite,"no")) && solver->restart_iter) {
    /* if it's a restart run, fast-forward the filename */
    int t;
    for (t=0; t<solver->restart_iter; t++) 
      if ((t+1)%solver->file_op_iter == 0) IncrementFilename(solver->op_filename);
  }

  return(0);
}
