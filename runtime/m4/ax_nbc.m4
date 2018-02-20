dnl Copyright 2013 The Trustees of Indiana University.
dnl
dnl Use, modification and distribution is subject to the Boost Software
dnl License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
dnl http://www.boost.org/LICENSE_1_0.txt)
dnl
dnl  Authors: Jeremiah Willcock
dnl           Andrew Lumsdaine
dnl 	      Marcin Zalewski

AC_DEFUN([TEST_NBC],
[
AC_ARG_WITH([nbc],
  [AS_HELP_STRING([--with-nbc@<:@=ARG@:>@],
    [use nbc library from a standard location (ARG=yes),
     from the specified location (ARG=<path>),
     or use NBC stub with MPI-3 (ARG=sutb)
     @<:@ARG=yes@:>@ ])],
    [
    if test "$withval" = "stub"; then
        want_nbc="stub"
    elif test "$withval" = "yes"; then
        want_nbc="yes"
        ac_nbc_path=""
    else
        want_nbc="yes"
        ac_nbc_path="$withval"
    fi
    ],
    [want_nbc="yes"])

AC_ARG_WITH([nbc-libdir],
        AS_HELP_STRING([--with-nbc-libdir=LIB_DIR],
        [Force given directory for nbc libraries. Note that this will override library path detection, so use this parameter only if default library detection fails and you know exactly where your nbc libraries are located.]),
        [
        if test -d "$withval"
        then
                ac_nbc_lib_path="$withval"
        else
                AC_MSG_ERROR(--with-nbc-libdir expected directory name)
        fi
        ],
        [ac_nbc_lib_path=""]
)

if test "x$want_nbc" = "xyes"; then

    dnl On 64-bit systems check for system libraries in both lib64 and lib.
    dnl The former is specified by FHS, but e.g. Debian does not adhere to
    dnl this (as it rises problems for generic multi-arch support).
    dnl The last entry in the list is chosen by default when no libraries
    dnl are found, e.g. when only header-only libraries are installed!
    libsubdirs="lib"
    ax_arch=`uname -m`
    case $ax_arch in
      x86_64|ppc64|s390x|sparc64|aarch64)
        libsubdirs="lib64 lib lib64"
        ;;
    esac

    dnl first we check the system location for nbc libraries
    if test "$ac_nbc_path" != ""; then
        NBC_CPPFLAGS="-I$ac_nbc_path/include"
        for ac_nbc_path_tmp in $libsubdirs; do
                if test -d "$ac_nbc_path"/"$ac_nbc_path_tmp" ; then
                        NBC_LDFLAGS="-L$ac_nbc_path/$ac_nbc_path_tmp"
                        break
                fi
        done
    elif test "$cross_compiling" != yes; then
        for ac_nbc_path_tmp in /usr /usr/local /opt /opt/local ; do
            if test -r "$ac_nbc_path_tmp/include/nbc.h"; then
                for libsubdir in $libsubdirs ; do
                    if ls "$ac_nbc_path_tmp/$libsubdir/libnbc_"* >/dev/null 2>&1 ; then break; fi
                done
                NBC_LDFLAGS="-L$ac_nbc_path_tmp/$libsubdir"
                NBC_CPPFLAGS="-I$ac_nbc_path_tmp/include"
                break;
            fi
        done
    fi

    dnl overwrite ld flags if we have required special directory with
    dnl --with-nbc-libdir parameter
    if test "$ac_nbc_lib_path" != ""; then
       NBC_LDFLAGS="-L$ac_nbc_lib_path"
    fi

else
    NBC_CPPFLAGS="-I$srcdir/libnbc-stub/include"
    HAVE_NBC_STUB=1
    AC_SUBST(HAVE_NBC_STUB)
fi

    CPPFLAGS_SAVED="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $NBC_CPPFLAGS"
    export CPPFLAGS

    LDFLAGS_SAVED="$LDFLAGS"
    LDFLAGS="$LDFLAGS $NBC_LDFLAGS"
    export LDFLAGS

    AC_REQUIRE([AX_PROG_CXX_MPI])
        AC_CACHE_CHECK(whether the NBC library is available,
                       ax_cv_nbc,
    [AC_LANG_PUSH(C++)
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <nbc.h>
    ]], [[
    NBC_Request nbcreq;
    MPI_Status mpistat;
    MPI_Comm mpicomm;
    MPI_Datatype mpidata;
    MPI_Op mpiop;

    NBC_Test(0, 0, 0);
    NBC_Ibarrier(mpicomm, 0);
    NBC_Iallreduce(0, 0, 0, mpidata, mpiop, mpicomm, 0);
    NBC_Wait(0, 0);
    ]])],[
    succeeded=yes
    found_system=yes
    ax_cv_nbc=$want_nbc
        ],[
	succeeded=no
	ax_cv_nbc=no
        ])
    AC_LANG_POP([C++])
    ])

    if test "$succeeded" != "yes" ; then
      AC_MSG_NOTICE([[We could not successuflly complie a program using nbclib.]])
      if test "$want_nbc" = "stub"; then
        AC_MSG_NOTICE([You have requested libnbc stub.  This requires MPI-3 compatible implementation.  Perhaps your implementations is not MPI-3 compatible.])
      fi
      # execute ACTION-IF-NOT-FOUND (if present):
      ifelse([$3], , :, [$3])
    else
    if test "x$ax_cv_nbc" = "xyes"; then
           NBCLIBDIR=`echo $NBC_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
                for libextension in `ls $NBCLIBDIR/libnbc*.{so,a}* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(nbc.*\)\.so.*$;\1;' -e 's;^lib\(nbc.*\)\.a*$;\1;'` ; do
                     ax_lib=${libextension}

                    AC_CHECK_LIB($ax_lib, exit,
                                 [NBC_LIB="-l$ax_lib"; AC_SUBST(NBC_LIB) link_nbc="yes"; break],
                                 [link_nbc="no"])
                done
                if test "x$link_nbc" != "xyes"; then
                for libextension in `ls $NBCLIBDIR/nbc*.{dll,a}* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(nbc.*\)\.dll.*$;\1;' -e 's;^\(nbc.*\)\.a*$;\1;'` ; do
                     ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                                 [NBC_LIB="-l$ax_lib"; AC_SUBST(NBC_LIB) link_nbc="yes"; break],
                                 [link_nbc="no"])
                done
		fi

		            if test "x$link_nbc" != "xyes"; then
                AC_MSG_ERROR(Could not link against $ax_lib !)
            fi
      fi

        AC_SUBST(NBC_CPPFLAGS)
        AC_SUBST(NBC_LDFLAGS)
        AC_DEFINE(HAVE_NBC,,[define if the Nbc library is available])
        # execute ACTION-IF-FOUND (if present):
        ifelse([$2], , :, [$2])
    fi



    CPPFLAGS="$CPPFLAGS_SAVED"
    LDFLAGS="$LDFLAGS_SAVED"



])
