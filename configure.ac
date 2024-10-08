dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -i
dnl
dnl You must have recent versions of Autoconf and Automake installed.

############################################################

AC_PREREQ(2.61)
AC_INIT([Specfem 2D], m4_normalize(m4_include([VERSION])), [see the wiki], [Specfem2D])
AC_CONFIG_SRCDIR([README.md])

AC_CONFIG_HEADERS([setup/config.h])

AC_CONFIG_MACRO_DIR([m4])
m4_include(m4/cit_backports.m4)

############################################################

AS_BOX([setting up compilation flags])
AC_CANONICAL_HOST

# 'configure' options

###
### precision
###

AC_ARG_ENABLE([double-precision],
    [AS_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test x"$want_double_precision" = xno; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])

###
### debugging
###

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [build with debugging options enabled @<:@default=no@:>@])],
    [want_debug="$enableval"],
    [want_debug=no])
AM_CONDITIONAL([COND_DEBUG], [test x"$want_debug" != xno])
export COND_DEBUG_FALSE
export COND_DEBUG_TRUE

###
### MPI
###

AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi],
        [build parallel version @<:@default=no@:>@])],
    [want_mpi="$withval"],
    [want_mpi=no])
AM_CONDITIONAL([COND_MPI], [test "$want_mpi" = yes])

###
### CUDA
###

AC_ARG_WITH([cuda],
    [AS_HELP_STRING([--with-cuda],
        [build CUDA GPU enabled version @<:@default=no@:>@])],
    [want_cuda="$withval"],
    [want_cuda=no])

# sets default: --with-cuda
if test x"$want_cuda" = xyes; then
    want_cuda=cuda4
fi
AM_CONDITIONAL([COND_CUDA], [test x"$want_cuda" != xno])
AM_CONDITIONAL([COND_CUDA4], [test x"$want_cuda" = xcuda4])
AM_CONDITIONAL([COND_CUDA5], [test x"$want_cuda" = xcuda5])
AM_CONDITIONAL([COND_CUDA6], [test x"$want_cuda" = xcuda6])
AM_CONDITIONAL([COND_CUDA7], [test x"$want_cuda" = xcuda7])
AM_CONDITIONAL([COND_CUDA8], [test x"$want_cuda" = xcuda8])
AM_CONDITIONAL([COND_CUDA9], [test x"$want_cuda" = xcuda9])
AM_CONDITIONAL([COND_CUDA10], [test x"$want_cuda" = xcuda10])
AM_CONDITIONAL([COND_CUDA11], [test x"$want_cuda" = xcuda11])
AM_CONDITIONAL([COND_CUDA12], [test x"$want_cuda" = xcuda12])

# cuda linking for cuda 5x and 6x and 7x and 8x and ..
AM_CONDITIONAL([COND_CUDA_PLUS],
  [test "$want_cuda" = cuda4 \
     -o "$want_cuda" = cuda5 \
     -o "$want_cuda" = cuda6 \
     -o "$want_cuda" = cuda7 \
     -o "$want_cuda" = cuda8 \
     -o "$want_cuda" = cuda9 \
     -o "$want_cuda" = cuda10 \
     -o "$want_cuda" = cuda11 \
     -o "$want_cuda" = cuda12 \]
)

###
### OpenMP
###

AC_ARG_WITH([openmp],
    [AS_HELP_STRING([--with-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [with_omp="$withval"],
    [with_omp=no])
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [enable_omp="$enableval"],
    [enable_omp=no])
AS_IF([test x"$enable_omp" != xno || test x"$with_omp" != xno],
    [want_omp=yes],[want_omp=no])
AM_CONDITIONAL([COND_OMP], [test x"$want_omp" != xno])

############################################################

# Checks for programs.

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi
# note: if not set, default FCFLAGS (and CFLAGS) will be set by autoconf/configure script (mostly -g -O2)
if test x"$FCFLAGS" = x && test x"$FLAGS_CHECK" != x ; then
    FCFLAGS=""
fi

AC_PROG_FC
export FC
export MPIFC
F77="$FC"
FFLAGS="$FCFLAGS"

AC_PROVIDE([AC_PROG_F77])
AC_SUBST([FCENV])
AC_SUBST(srcdir)

flags_guess="$SHELL $srcdir/flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags

AC_FC_WRAPPERS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

AC_FC_PP_DEFINE()
AC_SUBST([FC_DEFINE])
AC_FC_PP_SRCEXT(F90) dnl Because AC_FC_PP_DEFINE messes with things.

AC_PROG_CC

############################################################

#checks for Scotch
AS_BOX([SCOTCH])

export SCOTCH_DIR
export SCOTCH_LIBDIR
export SCOTCH_INCLUDEDIR
export USE_BUNDLED_SCOTCH

AC_ARG_WITH([scotch-dir],
AS_HELP_STRING([--with-scotch-dir=DIR],[define the root path to Scotch (e.g. /opt/scotch/)]),
[
 ac_scotch_dir="$withval";
])

AC_ARG_WITH([scotch-includedir],
AS_HELP_STRING([--with-scotch-includedir=DIR],[define the path to the Scotch headers (e.g. /opt/scotch/include)]),
[
 ac_scotch_include_dir="$withval";
])

AC_ARG_WITH([scotch-libdir],
AS_HELP_STRING([--with-scotch-libdir=DIR],[define the path to the Scotch libraries (e.g. /opt/scotch/lib)]),
[
 ac_scotch_lib_dir="$withval";
])

if test "${USE_BUNDLED_SCOTCH}" != "1"; then

   if test -z "${ac_scotch_lib_dir}"; then
      if test -n "${ac_scotch_dir}"; then
         ac_scotch_lib_dir="${ac_scotch_dir}/lib";
      else
        ac_scotch_lib_dir="/usr/lib";
      fi
   fi

  if test -z "${ac_scotch_include_dir}"; then
    if test -n "${ac_scotch_dir}"; then
      ac_scotch_include_dir="${ac_scotch_dir}/include";
    else
      ac_scotch_include_dir="/usr/include/scotch";
    fi
  fi

  scotch_lib=""
  ac_save_ldflags=${LDFLAGS}
  if test -n "${ac_scotch_lib_dir}"; then
    LDFLAGS="${LDFLAGS} -L${ac_scotch_lib_dir}";
  fi
  AC_CHECK_LIB(scotch,scotchfarchinit ,[scotch_lib="yes";],[scotch_lib="no";LDFLAGS=${ac_save_ldflags}], -lscotcherr)

  scotch_include=""
  ac_save_cppflags=${CPPFLAGS}
  # scotch only needed for fortran code
  #if test -n "${ac_scotch_include_dir}"; then
  #  CPPFLAGS="${CPPFLAGS} -I${ac_scotch_include_dir}";
  #fi

dnl This does not work because we are in Fortran mode, which does not
dnl handle headers.
dnl
dnl AC_CHECK_HEADER(scotchf.h,[scotch_include="yes";],[scotch_include="no";CPPFLAGS=${ac_save_cppflags}])

  scotch_usable="${scotch_lib}"

fi

AC_MSG_CHECKING([whether Scotch is usable])
if test "x${scotch_usable}" = "xyes"; then
 AC_DEFINE([HAVE_SCOTCH],[1],[defined if Scotch is installed])
 want_scotch=yes
 USE_BUNDLED_SCOTCH=0
 SCOTCH_DIR="${ac_scotch_dir}"
 SCOTCH_LIBDIR="${ac_scotch_lib_dir}"
 SCOTCH_INCLUDEDIR="${ac_scotch_include_dir}"
 AC_MSG_RESULT([yes])
else
 AC_DEFINE([HAVE_SCOTCH],[1],[defined if Scotch is installed])
 AC_MSG_RESULT([no, using bundled scotch instead])

 dnl autoconf 2.69 AC_PROG_LEX has no parameters
 dnl autoconf 2.70 AC_PROG_LEX prints a deprecation warning without params
 m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [-1], [dnl
   dnl autoconf < 2.70
   AC_PROG_LEX
 ], [
   dnl autoconf >= 2.70
   AC_PROG_LEX([noyywrap])
 ])
 if test -z "$LEX" || test "X$LEX" = "Xno"; then
   AC_MSG_ERROR([No suitable lex found])
 fi

 AC_PROG_YACC
 if test -z "$YACC" || test "X$YACC" = "Xno"; then
   AC_MSG_ERROR([No suitable yacc or bison found])
 fi

 ACX_PTHREAD(AC_MSG_RESULT([pthread found]), AC_MSG_ERROR([pthread not found]))

 # scotch only needed with mpi support
 if test "$want_mpi" = yes; then
  want_scotch=yes
  #daniel: scotch bundle
  # uses bundled scotch: current version pointed to by symbolic link scotch/
  USE_BUNDLED_SCOTCH=1

  SCOTCH_DIR="$srcdir/external_libs/scotch"
  SCOTCH_LIBDIR="${SCOTCH_DIR}/lib"
  SCOTCH_INCLUDEDIR="${SCOTCH_DIR}/include"
 else
  # no mpi support, scotch not needed
  want_scotch=no
  USE_BUNDLED_SCOTCH=0
 fi # want_mpi

fi
LDFLAGS=${ac_save_ldflags}

AM_CONDITIONAL([COND_SCOTCH], [test "$want_scotch" = yes])
AS_IF([test "$want_scotch" = yes],[
  AC_MSG_RESULT([SCOTCH is enabled])
],[
  AC_MSG_RESULT([SCOTCH is not enabled])
])


############################################################
# influential environment variables

AC_ARG_VAR(USE_BUNDLED_SCOTCH, [Set to 1 to always use the bundled Scotch library])
AC_ARG_VAR(SCOTCH_DIR, [Directory where Scotch is installed])
AC_ARG_VAR(SCOTCH_INCLUDEDIR, [Directory where Scotch headers are installed])
AC_ARG_VAR(SCOTCH_LIBDIR, [Directory where Scotch libraries are installed])

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])
AC_ARG_VAR(MPICC, [MPI C compiler command])

AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags])

# scratch disks
AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])

############################################################

# tests
AS_IF([test x"$MPIFC" = x],[
    MPIFC=mpif90
])

AS_IF([test x"$MPICC" = x],[
    MPICC=mpicc
])

AS_IF([test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x],[
    LOCAL_PATH_IS_ALSO_GLOBAL=true
])

# Checks for typedefs, structures, and compiler characteristics.
# check fortran modules flag
AS_BOX([module extensions])

AC_FC_MODULE_EXTENSION
AC_SUBST([FC_MODEXT])
AC_FC_MODULE_FLAG([], [])
AC_FC_MODULE_OUTPUT_FLAG([
	FCFLAGS_f90="$FC_MODOUT./obj $FC_MODINC./obj $FC_MODINC. $FCFLAGS_f90"
	FC_MODDIR=./obj
], [
	FC_MODDIR=.
])
AC_SUBST([FC_MODDIR])

# Checks for header files.

AC_LANG_PUSH(C)
AC_CHECK_HEADER(emmintrin.h,AC_DEFINE([HAVE_EMMINTRIN],[1],[Define if emmintrin.h]))
AC_CHECK_HEADER(xmmintrin.h,AC_DEFINE([HAVE_XMMINTRIN],[1],[Define if xmmintrin.h]))
# AIX doesn't have err.h so we need to check
AC_CHECK_HEADER(err.h,      AC_DEFINE([HAVE_ERR],      [1],[Define if err.h]))
AC_LANG_POP(C)

# Checks for libraries.

###
### MPI
###

AS_IF([test "$want_mpi" = yes], [
    AS_BOX([MPI])
    # checks MPI include directory
    CIT_MPI_INCDIR([$MPIFC])
])

###
### GPU
###

AS_IF([test x"$want_cuda" != xno], [
    AS_BOX([CUDA])
    CIT_CUDA_CONFIG
])

###
### OpenMP
###

AS_IF([test x"$want_omp" != xno], [
  AS_BOX([OpenMP])
  AC_MSG_NOTICE([OpenMP compilation is enabled])

  AC_ARG_VAR(OMP_FCFLAGS, [OpenMP Fortran compiler flags])
  AC_ARG_VAR(OMP_LIB, [Location of extra OpenMP libraries])

  # openmp checking
  AC_MSG_NOTICE([OpenMP flag $OMP_FCFLAGS])
  CIT_FC_OPENMP_MODULE([$FC], [$OMP_FCFLAGS])
])

##
## Git version info
##
AS_BOX([git version info])
CIG_PKG_GIT(SPECFEM2D)

AC_DEFINE_UNQUOTED([SPECFEM2D_VERSION], ["$PACKAGE_VERSION"], [Define SPECFEM2D version])
#AC_DEFINE_UNQUOTED([SPECFEM2D_DOI], ["10.5281/zenodo.xxxxx"], [Define SPECFEM2D doi])

AS_IF([test x"$SPECFEM2D_RELEASE_VERSION" = x0],[
  AC_MSG_NOTICE([building from git repository])
  GIT_PACKAGE_VERSION=$SPECFEM2D_GIT_REVISION
  GIT_COMMIT_VERSION=$SPECFEM2D_GIT_HASH
  GIT_DATE_VERSION=$SPECFEM2D_GIT_DATE
],[
  AC_MSG_NOTICE([not a git repository])
  GIT_PACKAGE_VERSION=$PACKAGE_VERSION
  GIT_COMMIT_VERSION="unknown"
  GIT_DATE_VERSION="unknown"
])
AC_SUBST([GIT_PACKAGE_VERSION])
AC_SUBST([GIT_COMMIT_VERSION])
AC_SUBST([GIT_DATE_VERSION])

############################################################

AS_BOX([setting up default simulation setup])

# Output results.
AC_CONFIG_FILES([
    Makefile
    setup/constants.h
    setup/constants_tomography.h
    setup/precision.h
    setup/config.fh
    setup/version.fh
])

# Scotch makefile
if test "$USE_BUNDLED_SCOTCH" = 1; then
  AC_CONFIG_FILES([
	  ${SCOTCH_DIR}/src/Makefile.inc
  ])
fi

AC_CONFIG_FILES([
    DATA/Par_file:DATA/Par_file
    DATA/SOURCE:DATA/SOURCE
])

AC_CONFIG_COMMANDS([bin], [AS_MKDIR_P(bin)])
AC_CONFIG_COMMANDS([obj], [AS_MKDIR_P(obj)])
AC_CONFIG_COMMANDS([DATA], [AS_MKDIR_P(DATA)])
AC_CONFIG_COMMANDS([OUTPUT_FILES], [AS_MKDIR_P(OUTPUT_FILES)])

AC_OUTPUT

dnl FYI, this is not a "result", but AS_ECHO is not available in 2.61, and it
dnl basically works like a portable echo.

AC_MSG_RESULT([])
AS_BOX([$PACKAGE_NAME $PACKAGE_VERSION])
AC_MSG_RESULT([])
AC_MSG_RESULT([$0 has completed and set up a default configuration to build.])
AC_MSG_RESULT([])
AC_MSG_RESULT([You may wish to modify the following files before running a simulation:])
AC_MSG_RESULT([AS_HELP_STRING([DATA/Par_file],
                              [Set parameters affecting the simulation.])])
AC_MSG_RESULT([AS_HELP_STRING([DATA/SOURCE],
                              [Set the source parameters before running the solver.])])
AC_MSG_RESULT([])

dnl end of configure.ac
