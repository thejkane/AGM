dnl Copyright 2013 The Trustees of Indiana University
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are permitted provided that the following conditions are met: 

dnl 1. Redistributions of source code must retain the above copyright notice, this
dnl    list of conditions and the following disclaimer. 
dnl 2. Redistributions in binary form must reproduce the above copyright notice,
dnl    this list of conditions and the following disclaimer in the documentation
dnl    and/or other materials provided with the distribution. 

dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
dnl ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
dnl WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
dnl DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
dnl ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
dnl (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
dnl LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
dnl ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
dnl (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl  Authors: Jeremiah Willcock
dnl           Andrew Lumsdaine
dnl 	      Marcin Zalewski

AC_PREREQ([2.69])
AC_INIT([AM++],[0.999],[ampp-users@crest.iu.edu],[amplusplus],[http://crest.iu.edu/projects/am++])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror])
AC_CONFIG_SRCDIR([am++/am++.hpp])
AC_CANONICAL_HOST

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
dnl AC_PROG_CXX
dnl AC_PROG_CC
AC_LANG([C++])
AX_PROG_CXX_MPI(,,[AC_MSG_FAILURE([Functional MPI compiler required.])])
# Does not work on crayc++
AX_CXXFLAGS_WARN_ALL
AM_PROG_AR

# Checks for libraries.
AX_BOOST_BASE([1.52],,[AC_MSG_FAILURE([BOOST libraries are required.])])
AX_BOOST_SYSTEM
AX_PTHREAD
AX_BOOST_THREAD
TEST_NBC(,,[AC_MSG_FAILURE([NBC libraries are required.])])

AC_SUBST([CPPFLAGS], ["$BOOST_CPPFLAGS $NBC_CPPFLAGS $CPPFLAGS"])
dnl $BOOST_CPPFLAGS 
AC_SUBST([LDFLAGS], ["$NBC_LDFLAGS $LDFLAGS"])
dnl $BOOST_LDFLAGS 
AC_SUBST([LIBS], ["$NBC_LIB $LIBS"])
dnl $BOOST_THREAD_LIB

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor sqrt])

# AM++ options
AMPP_THREADING
AMPP_SELFCHECK
AMPP_DEBUGGING

LT_INIT

# Until we can use other transports, we hardcode MPI
AC_DEFINE(TRANSPORT,[mpi],[Forces tests to use MPI transport. Currently, other transports are unsupported.])

AC_CONFIG_FILES([Makefile])

# Set this separately to work around bug in AX_CREATE_PKGCONFIG_INFO when $6 is non-empty
PACKAGE_LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
AX_CREATE_PKGCONFIG_INFO([libampp.pc], [], [$LIBS $BOOST_SYSTEM_LIB $BOOST_THREAD_LIB -lampp], [AM++ library], [$CPPFLAGS])

AC_OUTPUT

if test "x$HAVE_NBC_STUB" = "x1"; then
    AC_MSG_WARN([libnbc stub will only work with MPI-3 compatible MPI implementations. We do not check for MPI-3 compatibility explicitly, but a test libnbc program has been compiled successfully.])
fi
