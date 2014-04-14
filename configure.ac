AC_INIT([hypar],[0.1],[ghosh@mcs.anl.gov])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
echo -n AC_PACKAGE_VERSION>VERSION
AC_SUBST(PACKAGE_VERSION)
AC_MSG_NOTICE([hypar AC_PACKAGE_VERSION])

AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT('$(top_srcdir)')

AC_ARG_WITH([mpidir],AC_HELP_STRING([--with-mpidir],[Specify path where mpicc is installed.]))

if test "x$with_mpidir" != "x" ; then
  AC_MSG_CHECKING([for mpicc in $with_mpidir/bin/])
  if test -e $with_mpidir/bin/mpicc ; then
    AC_MSG_RESULT([yes])
    MPICC=$with_mpidir/bin/mpicc
  else
    AC_MSG_RESULT([no])
  fi
fi

if test "x$MPICC" = "x" ; then
  AC_CHECK_PROG(MPICC,mpicc,mpicc)
fi

AC_PROG_CC

if test "x$MPICC" = "x" ; then
  AC_MSG_WARN([Cannot find mpicc. Will compile with gcc. Use --with-mpidir to specify the location of mpicc.])
  CFLAGS="$CFLAGS -Dserial"
else
  CC=$MPICC
fi

CFLAGS="$CFLAGS -I/opt/local/include -I$includedir -std=c99"
LDFLAGS="$LDFLAGS -L/opt/local/lib"
LIBS="$LIBS -lm"

AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/ArrayFunctions/Makefile
 src/BoundaryConditions/Makefile
 src/FirstDerivative/Makefile
 src/HyParFunctions/Makefile
 src/InterpolationFunctions/Makefile
 src/MathFunctions/Makefile
 src/MPIFunctions/Makefile
 src/PhysicalModels/Makefile
   src/PhysicalModels/LinearADR/Makefile
   src/PhysicalModels/FPDoubleWell/Makefile
   src/PhysicalModels/FPPowerSystem/Makefile
   src/PhysicalModels/FPPowerSystem3Bus/Makefile
   src/PhysicalModels/Euler1D/Makefile
   src/PhysicalModels/Euler2D/Makefile
   src/PhysicalModels/NavierStokes2D/Makefile
   src/PhysicalModels/NavierStokes3D/Makefile
 src/SecondDerivative/Makefile
 src/TimeIntegration/Makefile
 src/TridiagLU/Makefile
])

AC_OUTPUT
