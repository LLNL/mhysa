SUBDIRS = \
  ArrayFunctions \
  BoundaryConditions\
  HyParFunctions \
  InterpolationFunctions \
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
  BoundaryConditions/libBoundaryConditions.a \
  PhysicalModels/LinearADR/libLinearADR.a \
  PhysicalModels/FPDoubleWell/libFPDoubleWell.a \
  PhysicalModels/FPPowerSystem/libFPPowerSystem.a \
  PhysicalModels/Euler1D/libEuler1D.a \
  InterpolationFunctions/libInterpolationFunctions.a \
  SecondDerivative/libSecondDerivative.a \
  MPIFunctions/libMPIFunctions.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a
