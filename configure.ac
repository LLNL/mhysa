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
if test "x$with_mpidir" != "x" ; then
  AC_MSG_CHECKING([for mpicxx in $with_mpidir/bin/])
  if test -e $with_mpidir/bin/mpicxx ; then
    AC_MSG_RESULT([yes])
    MPICXX=$with_mpidir/bin/mpicxx
  else
    AC_MSG_RESULT([no])
  fi
fi

if test "x$MPICC" = "x" ; then
  AC_CHECK_PROG(MPICC,mpicc,mpicc)
fi
if test "x$MPICXX" = "x" ; then
  AC_CHECK_PROG(MPICXX,mpicxx,mpicxx)
fi

AC_PROG_CC
AC_PROG_CXX

if test "x$MPICC" = "x" ; then
  AC_MSG_WARN([Cannot find mpicc. Will compile with gcc. Use --with-mpidir to specify the location of mpicc.])
  CFLAGS="$CFLAGS -Dserial"
else
  CC=$MPICC
fi
if test "x$MPICXX" = "x" ; then
  AC_MSG_WARN([Cannot find mpicxx. Will compile with g++. Use --with-mpidir to specify the location of mpicxx.])
else
  CXX=$MPICXX
fi

CFLAGS="$CFLAGS -I/opt/local/include -I$includedir"
CXXFLAGS=$CFLAGS
LDFLAGS="$LDFLAGS -L/opt/local/lib"
LIBS="$LIBS -lm"

AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/ArrayFunctions/Makefile
 src/BoundaryConditions/Makefile
 src/HyParFunctions/Makefile
 src/InterpolationFunctions/Makefile
 src/MathFunctions/Makefile
 src/MPIFunctions/Makefile
 src/PhysicalModels/Makefile
   src/PhysicalModels/LinearADR/Makefile
   src/PhysicalModels/FPDoubleWell/Makefile
   src/PhysicalModels/FPPowerSystem/Makefile
   src/PhysicalModels/Euler1D/Makefile
 src/SecondDerivative/Makefile
 src/TimeIntegration/Makefile
])

AC_OUTPUT
