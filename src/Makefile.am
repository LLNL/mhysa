SUBDIRS = \
  ArrayFunctions \
  BoundaryConditions\
  FirstDerivative \
  HyParFunctions \
  InterpolationFunctions \
  PETScFunctions \
  PhysicalModels \
  MathFunctions \
  MPIFunctions \
  SecondDerivative \
  TimeIntegration \
  TridiagLU

bin_PROGRAMS = HyPar
HyPar_SOURCES = main.c
HyPar_LDADD = \
  HyParFunctions/libHyParFunctions.a \
  TimeIntegration/libTimeIntegration.a \
  PETScFunctions/libPETScFunctions.a \
  BoundaryConditions/libBoundaryConditions.a \
  PhysicalModels/LinearADR/libLinearADR.a \
  PhysicalModels/FPDoubleWell/libFPDoubleWell.a \
  PhysicalModels/FPPowerSystem/libFPPowerSystem.a \
  PhysicalModels/FPPowerSystem1Bus/libFPPowerSystem1Bus.a \
  PhysicalModels/FPPowerSystem3Bus/libFPPowerSystem3Bus.a \
  PhysicalModels/Euler1D/libEuler1D.a \
  PhysicalModels/Euler2D/libEuler2D.a \
  PhysicalModels/NavierStokes2D/libNavierStokes2D.a \
  PhysicalModels/NavierStokes3D/libNavierStokes3D.a \
  PhysicalModels/Numa2D/libNuma2D.a \
  PhysicalModels/Numa3D/libNuma3D.a \
  InterpolationFunctions/libInterpolationFunctions.a \
  FirstDerivative/libFirstDerivative.a \
  SecondDerivative/libSecondDerivative.a \
  TridiagLU/libTridiagLU.a \
  MPIFunctions/libMPIFunctions.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a
