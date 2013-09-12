SUBDIRS = \
  ArrayFunctions \
  BoundaryConditions\
  HyParFunctions \
  MathFunctions \
  MPIFunctions

bin_PROGRAMS = HyPar
HyPar_SOURCES = main.c
HyPar_LDADD = \
  HyParFunctions/libHyParFunctions.a \
  BoundaryConditions/libBoundaryConditions.a \
  MPIFunctions/libMPIFunctions.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a
