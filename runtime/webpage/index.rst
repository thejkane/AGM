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

======
|Logo|
======

Overview
--------

AM++ [Willcock2010]_ is a user-level library for active messages based on the Active Pebbles [Willcock2011]_ programming model.  AM++ allows message handlers to be run in
an explicit loop that can be optimized and vectorized by the compiler and that can also be executed in parallel on multicore architectures. Runtime optimizations, such as message combining and
filtering, are also provided by the library, removing the need to implement that functionality at the application level.  

The overview of concepts of AM++ can be found in the cited papers.

Get AM++
--------

AM++ can be downloaded here:

  `am++-0.9999.tar.gz`__

__ am++-0.9999.tar.gz

For installation, see the documentation included in the release.  For online documentation, see here__.

__ doc/index.html

AM++ mailing list:

  ``ampp-users@crest.iu.edu``

----------------------------------------------------------------------------

Copyright (C) 2009-2013 The Trustees of Indiana University.

:Authors: 
          Jeremiah Willcock, 
	  Marcin Zalewski, 
	  and Andrew Lumsdaine

.. |Logo| image:: ampp-logo.png
            :align: middle
            :alt: AM++
            :target: http://crest.iu.edu/research/am++

.. _AM++\: a generalized active message framework: http://dl.acm.org/citation.cfm?id=1854323
.. _Active Pebbles\: Parallel Programming for Data-Driven Applications: http://dl.acm.org/citation.cfm?id=1995934
.. [Willcock2010] `AM++: a generalized active message framework`_ by J. Willcock, T. Hoefler, N. G. Edmonds, A. Lumsdaine. Proc. 19th Int. Conference on Parallel Architectures and Compilation Techniques (PACT`10).
.. [Willcock2011] `Active Pebbles: Parallel Programming for Data-Driven Applications`_ by J. Willcock, T. Hoefler, N. G. Edmonds, A. Lumsdaine.  Proc. Int. Conference on Supercomputing.
.. _Tutorial: tutorial.html
.. _Message type generators: generators.html
.. _Reductions: reductions.html
