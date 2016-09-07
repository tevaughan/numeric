#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([numeric], [1.2], [tevaughan@gmail.com])
AC_CONFIG_AUX_DIR([aux-scripts])
# AC_CONFIG_HEADERS reduces line length on standard output during build.
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([src/dimval.hpp])
AM_INIT_AUTOMAKE([-Wall foreign])
CPPFLAGS="$CPPFLAGS -finput-charset=UTF-8 -fexec-charset=UTF-8"
AC_LANG([C++])

AX_CODE_COVERAGE
if test -n "$CODE_COVERAGE_CXXFLAGS"; then
   CXXFLAGS="-Wall $CXXFLAGS"
else
   CXXFLAGS="-Wall -g -O2 $CXXFLAGS"
fi

# Checks for programs.

AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11

AM_PROG_AR
# LT_INIT must follow AM_PROG_AR.
LT_INIT

AC_PATH_PROG([PERL], [perl], [not-found])
if test "$PERL" = "not-found"; then
   echo "ERROR: perl not found"
   exit 1
fi

# Checks for libraries.

AC_CHECK_LIB([m], [pow])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([pow])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile])
AC_OUTPUT

