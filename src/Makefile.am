SUBDIRS = \
  AdvectionDiffusionReaction \
  ArrayFunctions \
  BoundaryConditions\
  HyParFunctions \
  MathFunctions \
  MPIFunctions \
  TimeIntegration

bin_PROGRAMS = HyPar
HyPar_SOURCES = main.c
HyPar_LDADD = \
  HyParFunctions/libHyParFunctions.a \
  TimeIntegration/libTimeIntegration.a \
  AdvectionDiffusionReaction/libAdvectionDiffusionReaction.a \
  BoundaryConditions/libBoundaryConditions.a \
  MPIFunctions/libMPIFunctions.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a
