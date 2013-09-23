SUBDIRS = \
  ArrayFunctions \
  BoundaryConditions\
  HyParFunctions \
  InterpolationFunctions \
  PETScFunctions \
  PhysicalModels \
  MathFunctions \
  MPIFunctions \
  SecondDerivative \
  TimeIntegration

bin_PROGRAMS = HyPar
HyPar_SOURCES = main.c
HyPar_LDADD = \
  HyParFunctions/libHyParFunctions.a \
  TimeIntegration/libTimeIntegration.a \
  PETScFunctions/libPETScFunctions.a \
  InterpolationFunctions/libInterpolationFunctions.a \
  SecondDerivative/libSecondDerivative.a \
  BoundaryConditions/libBoundaryConditions.a \
  MPIFunctions/libMPIFunctions.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a \
  PhysicalModels/LinearADR/libLinearADR.a \
  PhysicalModels/FPDoubleWell/libFPDoubleWell.a
