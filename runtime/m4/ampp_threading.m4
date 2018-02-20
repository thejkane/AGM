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

AC_DEFUN([AMPP_THREADING],
[
AC_ARG_ENABLE([builtin-atomics],
  [AS_HELP_STRING([--enable-builtin-atomics],
    [enable to use <atomic> header from C++11. 
     default: disabled ])],
    [
    if test "$enableval" = "yes"; then
      AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
      AC_SUBST([AMPP_ATOMICS_CPPFLAGS], ["-DAMPLUSPLUS_BUILTIN_ATOMICS"])
      AC_SUBST([CPPFLAGS], ["$CPPFLAGS $AMPP_ATOMICS_CPPFLAGS"])
      AC_CHECK_HEADER([atomic],,[AC_MSG_FAILURE([<atomic> header not found.])],)
      HAVE_BUILTIN_ATOMICS=1
    elif test "$enableval" != "no"; then
      AC_MSG_ERROR(--enable-builtin-atomics does not accept an argument)
      HAVE_BUILTIN_ATOMICS=0
    fi
    AC_SUBST(HAVE_BUILTIN_ATOMICS)
    ],
    [
    HAVE_BUILTIN_ATOMICS=0
    AC_SUBST(HAVE_BUILTIN_ATOMICS)
    ])

AC_ARG_ENABLE([threading],
  [AS_HELP_STRING([--enable-threading@<:@=ARG@:>@],
    [disable for no threading,
     use thread serialized (ARG=serialized),
     or use thread multiple (ARG=multiple)
     default: @<:@ARG=serialized@:>@ ])],
    [
    if test "$enableval" = "no"; then
      AC_DEFINE(AMPLUSPLUS_SINGLE_THREADED,,[define to use single threaded am++])
      AC_DEFINE(BOOST_SP_DISABLE_THREADS,,[define to use non-atomic reference counts for Boost shared pointers])
    elif test "$enableval" = "serialized"; then
      AC_DEFINE(AMPLUSPLUS_USE_THREAD_SERIALIZED,,[define to use thread serialized in am++])
      if test $HAVE_BUILTIN_ATOMICS = "0"; then
        AC_CHECK_HEADER([boost/atomic.hpp],,[AC_MSG_FAILURE([<boost/atomic.hpp> header required when using am++ threading and builtin atomics are not enabled.])],)
      fi
      # Does not work with crayc++
      # AS_CASE([$host],
      #   [*86*-*-*],  [AC_CHECK_HEADER([xmmintrin.h],,[AC_MSG_FAILURE([<xmmintrin.h> header required when using am++ threading.])],)]
      # )
    elif test "$enableval" = "multiple"; then
      if test $HAVE_BUILTIN_ATOMICS = "0"; then
        AC_CHECK_HEADER([boost/atomic.hpp],,[AC_MSG_FAILURE([<boost/atomic.hpp> header required when using am++ threading and builtin atomics are not enabled.])],)
      fi
      # Does not work with crayc++
      # AS_CASE([$host],
      #   [*86*-*-*],  [AC_CHECK_HEADER([xmmintrin.h],,[AC_MSG_FAILURE([<xmmintrin.h> header required when using am++ threading.])],)]
      # )
    else
      AC_MSG_ERROR(--enable-threading must be single, serialized, or multiple)
    fi
    ],
    [
      if test $HAVE_BUILTIN_ATOMICS = "0"; then
        AC_CHECK_HEADER([boost/atomic.hpp],,[AC_MSG_FAILURE([<boost/atomic.hpp> header required when using am++ threading and builtin atomics are not enabled.])],)
      fi
      # Does not work with crayc++
      # AS_CASE([$host],
      #   [*86*-*-*],  [AC_CHECK_HEADER([xmmintrin.h],,[AC_MSG_FAILURE([<xmmintrin.h> header required when using am++ threading.])],)]
      # )
      AC_DEFINE(AMPLUSPLUS_USE_THREAD_SERIALIZED,,[define to use thread serialized in am++])
    ])
])
