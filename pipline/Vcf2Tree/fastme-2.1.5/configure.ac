#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FastME],[2.1.5],[lefort@lirmm.fr])
AC_CONFIG_SRCDIR([src/interface_utilities.h])
AC_CONFIG_HEADERS([src/config.h])

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_CC

dnl for libraries.
AC_CHECK_LIB([m], [log])

dnl for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h float.h stdlib.h string.h])

dnl for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

dnl Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt strchr strstr])

dnl Add option for compiling with debugging informations
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Remove optimization options and add debug informations.])])
AS_IF([test "x$enable_debug" = "xyes"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O0 -g"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O2 -fomit-frame-pointer"])

dnl Add option for compiling with OpenMP
AC_ARG_ENABLE([OpenMP],
              [AS_HELP_STRING([--disable-OpenMP],
                              [Disable the use of OpenMP])],
              [],
              [AX_OPENMP])

dnl Add clang flags
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
#CFLAGS="$CFLAGS $OPENMP_CFLAGS -Weverything -Wno-padded"

LDFLAGS="$LDFLAGS $OPENMP_LDFLAG"

AC_CONFIG_FILES([Makefile
		 src/Makefile])
AC_OUTPUT

