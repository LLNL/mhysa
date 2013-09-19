SUBDIRS = \
  ArrayFunctions \
  BoundaryConditions\
  HyParFunctions \
  InterpolationFunctions \
  LinearADR \
  MathFunctions \
  MPIFunctions \
  TimeIntegration

bin_PROGRAMS = HyPar
HyPar_SOURCES = main.c
HyPar_LDADD = \
  HyParFunctions/libHyParFunctions.a \
  TimeIntegration/libTimeIntegration.a \
  InterpolationFunctions/libInterpolationFunctions.a \
  LinearADR/libLinearADR.a \
  BoundaryConditions/libBoundaryConditions.a \
  MPIFunctions/libMPIFunctions.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a
