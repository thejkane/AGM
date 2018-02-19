.. Copyright (C) 2009-2013 The Trustees of Indiana University.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met: 

   1. Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer. 
   2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution. 

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

============
|Logo|
============

Parallel Boost Graph Library 2 (Active Messages Version)
********************************************************

Overview
--------

This is a beta version of a message-based implementation of the Parallel Boost Graph Library (Parallel
BGL), a generic parallel graph library built on and fashioned after
the Boost Graph Library (BGL). This version builds on Active Pebbles
model and its AM++ (http://www.crest.iu.edu/projects/am++/) implementation. If you have any questions,
comments, or require help with the Parallel BGL, please see the
Parallel BGL web page.

We strongly suggest building and executing the test suite before
embarking on writing your own programs.  Familiarity with the Boost
Graph Library
(http://www.boost.org/libs/graph/doc/table_of_contents.html) is
assumed.

The overview of concepts of PBGL2 can be found in the
cited papers.  The previous versions (on which this one builds), can
be found here:

* http://www.osl.iu.edu/research/pbgl/software/
* http://www.boost.org/libs/graph_parallel/

Get PBGL2
---------

PBGL2 can be downloaded here:

  `pbgl2-1.0.tar.gz`__

__ pbgl2-1.0.tar.gz

For installation, see the documentation included in the release.

Contact Information
-------------------

Contact ``pbgl2-users@crest.iu.edu`` with questions and bug reports.

Bibliography
------------

.. [Edmonds:PhD]
   Edmonds, Nick. 2013. Active Messages as a Spanning Model For Parallel 
   Graph Computation, available at http://www.cs.indiana.edu/~ngedmond/dissertation.pdf

.. [edmonds13:ics]
   Edmonds, Nick, Jeremiah Willcock, and Andrew Lumsdaine. (2013) 
   "Expressing Graph Algorithms Using Generalized Active Messages". In 
   (Eds.) *International Conference on Supercomputing*, Eugene, Oregon.

.. [edmonds10:hipc]
   Edmonds, Nick, Torsten Hoefler, and Andrew Lumsdaine. (2010) "A 
   Space-Efficient Parallel Algorithm for Computing Betweenness 
   Centrality in Distributed Memory". In (Eds.) *International 
   Conference on High Performance Computing*, Goa, India. 

.. [edmonds10:hipc-srs]
   Edmonds, Nicholas, et al. (2010) "Design of a Large-Scale 
   Hybrid-Parallel Graph Library". In (Eds.) *International Conference 
   on High Performance Computing, Student Research Symposium*, Goa, 
   India. 


----------------------------------------------------------------------------

Copyright (C) 2009-2013 The Trustees of Indiana University.

:Authors: 
          Nicholas Edmonds 
	  and Andrew Lumsdaine

.. |Logo| image:: pbgl-logo.png
            :align: middle
            :alt: PBGL2
            :target: http://crest.iu.edu/research/pbgl2
