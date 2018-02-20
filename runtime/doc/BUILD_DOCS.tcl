#!/bin/sh
# \
exec tclsh "$0" "$@"

# Copyright (C) 2009 The Trustees of Indiana University.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Authors: Jeremiah Willcock, Andrew Lumsdaine

foreach input [glob *.rst] {
  set output [file join html "[file rootname $input].html"]
  puts "Processing $input -> $output"
  set processor [open "|rst2html-2.4.py --stylesheet=rst.css -gdt --link-stylesheet --traceback --trim-footnote-reference-space --footnote-references=superscript >$output" w]
  set inputfd [open $input r]
  set data [read $inputfd]
  close $inputfd
  foreach line [split $data \n] {
    if {[regexp {^\.\. image:: (http:.*)$} $line _ url]} {
      set tag $url
      regsub -all {.*/} $tag {} tag
      regsub -all {[^a-zA-Z0-9]} $tag _ tag
      set imageoutput [file join html "$tag.png"]
      puts "Getting image $url -> $imageoutput"
      exec wget -q -O $imageoutput $url
      puts $processor ".. image:: [file tail $imageoutput]"
    } else {
      puts $processor $line
    }
  }
  close $processor
}
