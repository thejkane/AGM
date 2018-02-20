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

AC_DEFUN([AMPP_DEBUGGING],
[

AC_ARG_ENABLE([debugging],
  [AS_HELP_STRING([--enable-debugging],
    [enable am++ debugging,
     default: disabled])],
    [
    if test "$enableval" = "no"; then
      AC_DEFINE(NDEBUG,,[define to disable am++ debugging])
    elif test "$enableval" != "yes"; then
      AC_MSG_ERROR(--enable-debugging does not accept an argument)
    fi
    ],
    [
    AC_DEFINE(NDEBUG,,[define to disable am++ debugging])
    ])
])
