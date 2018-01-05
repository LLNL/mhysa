AC_INIT([mhysa],[1.0],[debojyoti.ghosh@gmail.com])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
echo -n AC_PACKAGE_VERSION>VERSION
AC_SUBST(PACKAGE_VERSION)
AC_MSG_NOTICE([mhysa AC_PACKAGE_VERSION])

AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT('$(top_srcdir)')

AC_ARG_WITH([mpi_dir],AC_HELP_STRING([--with-mpi-dir],[Specify path where MPI is installed.]))
AC_ARG_ENABLE([omp],AS_HELP_STRING([--enable-omp],[Enable OpenMP threads]))
AC_ARG_ENABLE([scalapack],AS_HELP_STRING([--enable-scalapack],[Enable ScaLAPACK]))
AC_ARG_WITH([blas_dir],AC_HELP_STRING([--with-blas-dir],[Specify path where BLAS libraries are installed.]))
AC_ARG_WITH([lapack_dir],AC_HELP_STRING([--with-lapack-dir],[Specify path where LAPACK libraries are installed.]))
AC_ARG_WITH([scalapack_dir],AC_HELP_STRING([--with-scalapack-dir],[Specify path where ScaLAPACK libraries are installed.]))

AC_PROG_CC

CFLAGS="$CFLAGS -I$includedir"
if test "$CC" = "gcc" ; then
  CFLAGS="$CFLAGS -std=c99"
  LIBS="$LIBS -lm"
  MPI_CC=mpicc
else 
  if test "$CC" = "xlc" ; then
    CFLAGS="$CFLAGS -qlanglvl=stdc99"
    LIBS="$LIBS"
    MPI_CC=mpixlc
  else 
    if test "$CC" = "xlc_r" ; then
      CFLAGS="$CFLAGS -qlanglvl=stdc99"
      LIBS="$LIBS"
      MPI_CC=mpixlc_r
    fi
  fi
fi

if test "x$enable_omp" = "xyes" ; then
  AC_MSG_NOTICE([Compiling with OpenMP support.])
  if test "$CC" = "gcc"  ; then
    CFLAGS="$CFLAGS -fopenmp"
    LIBS="$LIBS -lgomp"
  else
    if test "$CC" = "xlc" || test "$CC" = "xlc_r" ; then
      CFLAGS="$CFLAGS -qsmp=omp"
      LIBS="$LIBS"
    fi
  fi
  CFLAGS="$CFLAGS -Dwith_omp"
fi

if test "x$with_mpi_dir" != "x" ; then
  AC_MSG_CHECKING([for ${MPI_CC} in $with_mpi_dir/bin/])
  if test -e $with_mpi_dir/bin/${MPI_CC} ; then
    AC_MSG_RESULT([yes])
    MPICC=$with_mpi_dir/bin/${MPI_CC}
  else
    AC_MSG_RESULT([no])
  fi
else
  if test "x$MPI_DIR" != "x" ; then
    AC_MSG_NOTICE([Found environment variable MPI_DIR=${MPI_DIR}.])
    AC_MSG_CHECKING([for ${MPI_CC} in ${MPI_DIR}/bin/])
    if test -e ${MPI_DIR}/bin/${MPI_CC} ; then
      AC_MSG_RESULT([yes])
      MPICC=${MPI_DIR}/bin/${MPI_CC}
    else
      AC_MSG_RESULT([no])
    fi
  fi
fi

if test "x$MPICC" = "x" ; then
  AC_CHECK_PROG(MPICC,${MPI_CC},${MPI_CC})
fi

if test "x$MPICC" = "x" ; then
  AC_MSG_WARN([Cannot find ${MPI_CC}. Will compile with gcc. Use --with-mpi-dir to specify the location of ${MPI_CC}.])
  CFLAGS="$CFLAGS -Dserial"
else
  CC=$MPICC
fi

if test "x$enable_scalapack" = "xyes" ; then
  AC_MSG_NOTICE([Compiling with ScaLAPACK support.])
  CFLAGS="$CFLAGS -Dwith_scalapack"
  if test "x$with_scalapack_dir" != "x" ; then
    LDFLAGS="$LDFLAGS -Wl,-rpath,$with_scalapack_dir/lib -L$with_scalapack_dir/lib"
  fi
  if test "x$with_blas_dir" != "x" ; then
    LDFLAGS="$LDFLAGS -Wl,-rpath,$with_blas_dir/lib -L$with_blas_dir/lib"
  fi
  if test "x$with_lapack_dir" != "x" ; then
    LDFLAGS="$LDFLAGS -Wl,-rpath,$with_lapack_dir/lib -L$with_lapack_dir/lib"
  fi
  LIBS="$LIBS -lscalapack -lblas -llapack"
fi

if test "x$PETSC_DIR" = "x"; then
  AC_MSG_NOTICE([Environment variable PETSC_DIR not set. Compiling without PETSc interface.])
else
  if test [ ! -d "$PETSC_DIR" ]; then
    AC_MSG_NOTICE([Environment variable PETSC_DIR set, but PETSc not found there. Compiling without PETSc interface.])
  else
    if test "x$PETSC_ARCH" = "x"; then
      AC_MSG_NOTICE([Environment variable PETSC_ARCH not set. Compiling without PETSc interface.])
    else
      AC_MSG_NOTICE([Compiling with PETSc interface.])
cat <<MFILE > petscmake$$
PETSC_ARCH = ${PETSC_ARCH}
PETSC_DIR  = ${PETSC_DIR}
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
MFILE
      PETSC_INC="`make -f petscmake$$ getincludedirs`"
      PETSC_LIB="`make -f petscmake$$ getlinklibs`"
      rm petscmake$$
      CFLAGS="$CFLAGS -Dwith_petsc $PETSC_INC"
      CXXFLAGS=$CFLAGS
      LIBS="$LIBS $PETSC_LIB"
    fi
  fi
fi

AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/ArrayFunctions/Makefile
 src/BandedMatrix/Makefile
 src/BoundaryConditions/Makefile
 src/FirstDerivative/Makefile
 src/CoreFunctions/Makefile
 src/ImmersedBoundaries/Makefile
 src/InterpolationFunctions/Makefile
 src/IOFunctions/Makefile
 src/MathFunctions/Makefile
 src/MPIFunctions/Makefile
 src/PETScFunctions/Makefile
 src/PhysicalModels/Makefile
   src/PhysicalModels/Euler1D/Makefile
   src/PhysicalModels/NavierStokes3D/Makefile
 src/SecondDerivative/Makefile
 src/TimeIntegration/Makefile
 src/TridiagLU/Makefile
])

AC_OUTPUT
