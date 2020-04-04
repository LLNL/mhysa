SUBDIRS = \
  ArrayFunctions \
  BandedMatrix \
  BoundaryConditions\
  FirstDerivative \
  CoreFunctions \
  ImmersedBoundaries \
  InterpolationFunctions \
  IOFunctions \
  PETScFunctions \
  PhysicalModels \
  MathFunctions \
  MPIFunctions \
  SecondDerivative \
  TimeIntegration \
  TridiagLU

bin_PROGRAMS = mhysa
mhysa_SOURCES = main.c
mhysa_LDADD = \
  CoreFunctions/libCoreFunctions.a \
  TimeIntegration/libTimeIntegration.a \
  PETScFunctions/libPETScFunctions.a \
  BoundaryConditions/libBoundaryConditions.a \
  PhysicalModels/Euler1D/libEuler1D.a \
  PhysicalModels/NavierStokes3D/libNavierStokes3D.a \
  ImmersedBoundaries/libImmersedBoundaries.a \
  InterpolationFunctions/libInterpolationFunctions.a \
  FirstDerivative/libFirstDerivative.a \
  SecondDerivative/libSecondDerivative.a \
  TridiagLU/libTridiagLU.a \
  IOFunctions/libIOFunctions.a \
  MPIFunctions/libMPIFunctions.a \
  BandedMatrix/libBandedMatrix.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a
