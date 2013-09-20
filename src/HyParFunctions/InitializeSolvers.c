#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <secondderivative.h>

/* Function declarations */
int WriteText                   (int,int,int*,double*,double*,char*,int*);
int ApplyBoundaryConditions     (void*,void*,double*);
int HyperbolicFunction          (double*,double*,void*,void*,double);
int ParabolicFunctionNC1Stage   (double*,double*,void*,void*,double);
int ParabolicFunctionCons1Stage (double*,double*,void*,void*,double);
int SourceFunction              (double*,double*,void*,void*,double);

int InitializeSolvers(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0;

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
  if (!strcmp(solver->spatial_scheme_hyp,_FIRST_ORDER_UPWIND_)) {
    solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwind;
    solver->interp = NULL;
  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_WENO_)) {
    solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENO;
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    ierr = WENOInitialize(solver->interp,mpi); CHECKERR(ierr);
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
    ierr = TimeMSTIInitialize(solver->time_scheme,solver->time_scheme_type,
                              solver->msti); CHECKERR(ierr);
  } else {
    fprintf(stderr,"Error: %s is a not a supported time-integration scheme.\n",
            solver->time_scheme);
    return(1);
  }

  /* Solution output function */
  if (!strcmp(solver->op_file_format,"text")) {
    solver->WriteOutput = WriteText;
    if (!strcmp(solver->op_overwrite,"no")) strcpy(solver->op_filename,"op_00000.dat");
    else                                    strcpy(solver->op_filename,"op.dat");
  } else if (!strcmp(solver->op_file_format,"none")) {
    solver->WriteOutput = NULL;
  } else {
    fprintf(stderr,"Error: %s is not a supported file format.\n",solver->op_file_format);
    return(1);
  }

  return(0);
}
