/*! @file InitializeSolvers.c
    @author Debojyoti Ghosh
    @brief Initialize all solvers
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <io.h>
#include <hypar.h>
#include <tridiagLU.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <firstderivative.h>
#include <secondderivative.h>

/* Function declarations */
int  ApplyBoundaryConditions     (void*,void*,double*,double*,double);
int  ApplyIBConditions           (void*,void*,double*,double);
int  HyperbolicFunction          (double*,double*,void*,void*,double,int,
                                  int(*)(double*,double*,int,void*,double),
                                  int(*)(double*,double*,double*,double*,double*,
                                         double*,int,void*,double));
int  ParabolicFunctionNC1Stage   (double*,double*,void*,void*,double);
int  ParabolicFunctionNC2Stage   (double*,double*,void*,void*,double);
int  ParabolicFunctionNC1_5Stage (double*,double*,void*,void*,double);
int  ParabolicFunctionCons1Stage (double*,double*,void*,void*,double);
int  SourceFunction              (double*,double*,void*,void*,double);
int  VolumeIntegral              (double*,double*,void*,void*);
int  BoundaryIntegral            (void*,void*);
int  CalculateConservationError  (void*,void*);
void IncrementFilenameIndex      (char*,int);
int  NonLinearInterpolation      (double*,void*,void*,double,
                                  int(*)(double*,double*,int,void*,double));

/*! This function initializes all solvers-specific function pointers 
    depending on user input. The specific functions used for spatial
    discretization, time integration, and solution output are set here.
*/
int InitializeSolvers(
                        void *s, /*!< Solver object of type #HyPar */
                        void *m  /*!< MPI object of type #MPIVariables */
                     )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  if (!mpi->rank) printf("Initializing solvers.\n");

  solver->ApplyBoundaryConditions     = ApplyBoundaryConditions;
  solver->ApplyIBConditions           = ApplyIBConditions;
  solver->HyperbolicFunction          = HyperbolicFunction;
  solver->SourceFunction              = SourceFunction;
  solver->VolumeIntegralFunction      = VolumeIntegral;
  solver->BoundaryIntegralFunction    = BoundaryIntegral;
  solver->CalculateConservationError  = CalculateConservationError;
  solver->NonlinearInterp             = NonLinearInterpolation;

  /* choose the type of parabolic discretization */
  solver->ParabolicFunction         = NULL;
  solver->SecondDerivativePar       = NULL;
  solver->FirstDerivativePar        = NULL;
  solver->InterpolateInterfacesPar  = NULL;
  if (!strcmp(solver->spatial_type_par,_NC_1STAGE_)) {
    solver->ParabolicFunction = ParabolicFunctionNC1Stage;
    if (!strcmp(solver->spatial_scheme_par,_SECOND_ORDER_CENTRAL_)) {
      solver->SecondDerivativePar      = SecondDerivativeSecondOrderCentral; 
    } else if (!strcmp(solver->spatial_scheme_par,_FOURTH_ORDER_CENTRAL_)) {
      solver->SecondDerivativePar      = SecondDerivativeFourthOrderCentral; 
    } else {
      fprintf(stderr,"Error: %s is not a supported ",solver->spatial_scheme_par);
      fprintf(stderr,"spatial scheme of type %s for the parabolic terms.\n",
            solver->spatial_type_par);
    }
  } else if (!strcmp(solver->spatial_type_par,_NC_2STAGE_)) {
    solver->ParabolicFunction = ParabolicFunctionNC2Stage;
    if (!strcmp(solver->spatial_scheme_par,_SECOND_ORDER_CENTRAL_)) {
      solver->FirstDerivativePar       = FirstDerivativeFirstOrder; 
      /* why first order? see ParabolicFunctionNC2Stage.c. 2nd order central
         approximation to the 2nd derivative can be expressed as a conjugation
         of 1st order approximations to the 1st derivative (one forward and 
         one backward) -- this prevents odd-even decoupling */ 
    } else if (!strcmp(solver->spatial_scheme_par,_FOURTH_ORDER_CENTRAL_)) {
      solver->FirstDerivativePar       = FirstDerivativeFourthOrderCentral; 
      /* why 4th order? I could not derive the decomposition of the 
         4th order central approximation to the 2nd derivative! Some problems
         may show odd-even decoupling */ 
    } else {
      fprintf(stderr,"Error: %s is not a supported ",solver->spatial_scheme_par);
      fprintf(stderr,"spatial scheme of type %s for the parabolic terms.\n",
            solver->spatial_type_par);
    }
  } else if (!strcmp(solver->spatial_type_par,_NC_1_5STAGE_)) {
    solver->ParabolicFunction = ParabolicFunctionNC1_5Stage;
    if (!strcmp(solver->spatial_scheme_par,_SECOND_ORDER_CENTRAL_)) {
      solver->FirstDerivativePar       = FirstDerivativeSecondOrderCentral; 
      solver->SecondDerivativePar      = SecondDerivativeSecondOrderCentral; 
    } else if (!strcmp(solver->spatial_scheme_par,_FOURTH_ORDER_CENTRAL_)) {
      solver->FirstDerivativePar       = FirstDerivativeFourthOrderCentral; 
      solver->SecondDerivativePar      = SecondDerivativeFourthOrderCentral; 
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
  solver->interp                = NULL;
  solver->compact               = NULL;
  solver->lusolver              = NULL;
  solver->SetInterpLimiterVar   = NULL;
  solver->flag_nonlinearinterp  = 1;
  if (strcmp(solver->interp_type,_CHARACTERISTIC_) && strcmp(solver->interp_type,_COMPONENTS_)) {
    fprintf(stderr,"Error in InitializeSolvers(): %s is not a ",solver->interp_type);
    fprintf(stderr,"supported interpolation type.\n");
    return(1);
  }
  if (!strcmp(solver->spatial_scheme_hyp,_FIRST_ORDER_UPWIND_)) {

    /* First order upwind scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwindChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimFirstOrderUpwind;
    }

  } else if (!strcmp(solver->spatial_scheme_hyp,_SECOND_ORDER_CENTRAL_)) {

    /* Second order central scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimSecondOrderCentralChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimSecondOrderCentral;
    }

  } else if (!strcmp(solver->spatial_scheme_hyp,_THIRD_ORDER_MUSCL_)) {

    /* Third order MUSCL scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimThirdOrderMUSCLChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimThirdOrderMUSCL;
    }
    solver->interp = (MUSCLParameters*) calloc(1,sizeof(MUSCLParameters));
    IERR MUSCLInitialize(solver->interp,mpi); CHECKERR(ierr);

  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_UPWIND_)) {

    /* Fifth order upwind scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderUpwindChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderUpwind;
    }

  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_COMPACT_UPWIND_)) {

    /* Fifth order compact upwind scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCompactUpwindChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCompactUpwind;
    }
    solver->compact = (CompactScheme*) calloc(1,sizeof(CompactScheme));
    IERR CompactSchemeInitialize(solver,mpi,solver->interp_type);
    solver->lusolver = (TridiagLU*) calloc (1,sizeof(TridiagLU));
    IERR tridiagLUInit(solver->lusolver,&mpi->world);CHECKERR(ierr);

  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_WENO_)) {

    /* Fifth order WENO scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENOChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderWENO;
    }
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp,solver->interp_type); CHECKERR(ierr);
    solver->flag_nonlinearinterp = !(((WENOParameters*)solver->interp)->no_limiting);

  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {

    /* Fifth order CRWENO scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCRWENOChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderCRWENO;
    }
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp,solver->interp_type); CHECKERR(ierr);
    solver->flag_nonlinearinterp = !(((WENOParameters*)solver->interp)->no_limiting);
    solver->compact = (CompactScheme*) calloc(1,sizeof(CompactScheme));
    IERR CompactSchemeInitialize(solver,mpi,solver->interp_type);
    solver->lusolver = (TridiagLU*) calloc (1,sizeof(TridiagLU));
    IERR tridiagLUInit(solver->lusolver,&mpi->world);CHECKERR(ierr);

  } else if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_)) {

    /* Fifth order HCWENO scheme */
    if ((solver->nvars > 1) && (!strcmp(solver->interp_type,_CHARACTERISTIC_))) {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderHCWENOChar;
    } else {
      solver->InterpolateInterfacesHyp = Interp1PrimFifthOrderHCWENO;
    }
    solver->interp = (WENOParameters*) calloc(1,sizeof(WENOParameters));
    IERR WENOInitialize(solver,mpi,solver->spatial_scheme_hyp,solver->interp_type); CHECKERR(ierr);
    solver->flag_nonlinearinterp = !(((WENOParameters*)solver->interp)->no_limiting);
    solver->compact = (CompactScheme*) calloc(1,sizeof(CompactScheme));
    IERR CompactSchemeInitialize(solver,mpi,solver->interp_type);
    solver->lusolver = (TridiagLU*) calloc (1,sizeof(TridiagLU));
    IERR tridiagLUInit(solver->lusolver,&mpi->world);CHECKERR(ierr);

  } else {
    fprintf(stderr,"Error: %s is a not a supported spatial interpolation scheme.\n",
            solver->spatial_scheme_hyp);
    return(1);
  }

  /* Time integration */
  solver->time_integrator = NULL;
#ifdef with_petsc
  if (solver->use_petscTS) {
    /* dummy -- not used */
    solver->TimeIntegrate = TimeForwardEuler;
    solver->msti = NULL;
  } else {
    if (!strcmp(solver->time_scheme,_FORWARD_EULER_)) { 
      solver->TimeIntegrate = TimeForwardEuler;
      solver->msti = NULL;
    } else if (!strcmp(solver->time_scheme,_RK_)) {
      solver->TimeIntegrate = TimeRK;
      solver->msti = (ExplicitRKParameters*) calloc (1,sizeof(ExplicitRKParameters));
      IERR TimeExplicitRKInitialize(solver->time_scheme,solver->time_scheme_type,
                                    solver->msti,mpi); CHECKERR(ierr);
    } else if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
      solver->TimeIntegrate = TimeGLMGEE;
      solver->msti = (GLMGEEParameters*) calloc (1,sizeof(GLMGEEParameters));
      IERR TimeGLMGEEInitialize(solver->time_scheme,solver->time_scheme_type,
                                solver->msti,mpi); CHECKERR(ierr);
    } else {
      fprintf(stderr,"Error: %s is a not a supported time-integration scheme.\n",
              solver->time_scheme);
      return(1);
    }
  }
#else
  if (!strcmp(solver->time_scheme,_FORWARD_EULER_)) { 
    solver->TimeIntegrate = TimeForwardEuler;
    solver->msti = NULL;
  } else if (!strcmp(solver->time_scheme,_RK_)) {
    solver->TimeIntegrate = TimeRK;
    solver->msti = (ExplicitRKParameters*) calloc (1,sizeof(ExplicitRKParameters));
    IERR TimeExplicitRKInitialize(solver->time_scheme,solver->time_scheme_type,
                                  solver->msti,mpi); CHECKERR(ierr);
  } else if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
    solver->TimeIntegrate = TimeGLMGEE;
    solver->msti = (GLMGEEParameters*) calloc (1,sizeof(GLMGEEParameters));
    IERR TimeGLMGEEInitialize(solver->time_scheme,solver->time_scheme_type,
                              solver->msti,mpi); CHECKERR(ierr);
  } else {
    fprintf(stderr,"Error: %s is a not a supported time-integration scheme.\n",
            solver->time_scheme);
    return(1);
  }
#endif

  /* Solution output function */
  solver->WriteOutput    = NULL; /* default - no output */
  solver->filename_index = NULL;
  if (!strcmp(solver->output_mode,"serial")) {
    solver->index_length = 5;
    solver->filename_index = (char*) calloc (solver->index_length+1,sizeof(char));
    int i; for (i=0; i<solver->index_length; i++) solver->filename_index[i] = '0';
    solver->filename_index[solver->index_length] = (char) 0;
    if (!strcmp(solver->op_file_format,"text")) {
      solver->WriteOutput = WriteText;
      strcpy(solver->solnfilename_extn,".dat");
    } else if (!strcmp(solver->op_file_format,"tecplot2d")) {
      solver->WriteOutput = WriteTecplot2D;
      strcpy(solver->solnfilename_extn,".dat");
    } else if (!strcmp(solver->op_file_format,"tecplot3d")) {
      solver->WriteOutput = WriteTecplot3D;
      strcpy(solver->solnfilename_extn,".dat");
    } else if ((!strcmp(solver->op_file_format,"binary")) || (!strcmp(solver->op_file_format,"bin"))) {
      solver->WriteOutput = WriteBinary;
      strcpy(solver->solnfilename_extn,".bin");
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
        if ((t+1)%solver->file_op_iter == 0) IncrementFilenameIndex(solver->filename_index,solver->index_length);
    }
  } else if (!strcmp(solver->output_mode,"parallel")) {
    if (!strcmp(solver->op_file_format,"none")) solver->WriteOutput = NULL;
    else {
      /* only binary file writing supported in parallel mode */
      /* use post-processing scripts to convert              */
      solver->WriteOutput = WriteBinary;
      strcpy(solver->solnfilename_extn,".bin");
    }
  } else {
    fprintf(stderr,"Error: %s is not a supported output mode.\n",solver->output_mode);
    fprintf(stderr,"Should be \"serial\" or \"parallel\".    \n");
    return(1);
  }

  return(0);
}
