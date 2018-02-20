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

==============================
|Logo| Message Type Generators
==============================

.. contents::

Overview
--------

Message types encapsulate orthogonal properties such as coalescing and caching.  To allow layering of these orthogonal properties, AM++ provides *message generators* that take as input another message type generator and tack on some additional features.  When all the desired properties of a message type are captured by layering of the generators, the final message type can then be extracted.  Not only does such design allow for orthogonal treatment of concerns, but it also allows for generic code to receive a generator as input and then, in a generic fashion, to add additional features to that generator in a generic fashion.

In general, generators have a compile time interface that generates the type of the resulting message type and a runtime interface for collecting the necessary inputs and for producing values of the final message type.

Generator Concept
-----------------

Notation
~~~~~~~~

``Gen``
  A type that models the Generator concept.

``ArgType``
  The argument type.

``HandlerType``
  A handler type that models one of the `handler objects`__.

__ handlers.html

``Owner``
  The owner map type that models the `Owner Map`__ concept.

__ oba.html#owner-map-concept

``Reduction``
  The reduction type that models the `reduction concept`__.

__ reductions.html#reduction-concept

``type``
  The result of message type computation.

Associated Type Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------+-----------------------+-----------------------------------+
|Message type computation|``Gen:: template       |The result type must model one of  |
|                        |call_result<ArgType,   |the `message type concepts`_.      |
|                        |HandlerType, Owner,    |                                   |
|                        |Reduction>::type``     |                                   |
+------------------------+-----------------------+-----------------------------------+
|Message type constructor|``type(const Gen& gen, |Constructor of the generated       |
|                        |transport, const       |message type.                      |
|                        |Owner&, const          |                                   |
|                        |Reduction&)``          |                                   |
+------------------------+-----------------------+-----------------------------------+



Simple Interface
~~~~~~~~~~~~~~~~

::

  template <typename ArgType, typename HandlerType, typename Gen, typename Owner, typename Reduction = no_reduction_t>
  struct message_type_generator_result;

This template allows a more direct access to the result of message type generation.  Instead having to first access ``call_result`` and then its member ``type``, the result can be extracted directly::

  typedef message_type_generator_result<arg_t, hadler_t, gen_t, owner_t, reduction_t>::type message_t;

Models
------

Built-in generators are defined in:
  
  <``am++/message_type_generators.hpp``>

Non-Coalesced Generator
~~~~~~~~~~~~~~~~~~~~~~~

::
  
  struct noncoalesced_generator;

This generator produces a message type that models the `Object-Addressable Message Type`_ concept.  In contrast to other generators, this generator does not require another generator as input.

This generator does not use the ``Reduction`` argument in the ``call_result`` type function, and the argument defaults to ``no_reduction_t``.

Simple Generator
~~~~~~~~~~~~~~~~

::

  template <typename CoalescingGen>
  struct simple_generator;

The simple generator takes an underlying `coalescing generator`_ and adds `object-based addressing`_.  This generator produces a message type that models the `Object-Addressable Message Type`_ concept.

.. rubric:: Constructor

::

  explicit simple_generator(const CoalescingGen& coalescing_gen);

Routing Generator
~~~~~~~~~~~~~~~~~

::

  template <typename CoalescingGen, typename Routing>
  struct routing_generator;

The routing generator takes an underlying `coalescing generator`_ and adds `object-based addressing`_ and routing_.

.. rubric:: Constructor

::

  explicit routing_generator(const CoalescingGen& c, const Routing& routing);


Cache Generator
~~~~~~~~~~~~~~~

::

  template <typename CoalescingGen, typename Routing>
  struct cache_generator

The cache generator takes an underlying `coalescing generator`_ and adds `object-based addressing`_, routing_, and caching_.  The type of caching is chosen appropriately, depending on the reduction type passed to ``call_result``.

.. rubric:: Constructor

::

  explicit cache_generator(const CoalescingGen& c, unsigned int lg_size, const Routing& routing = no_routing());

Per-Thread Cache Generator
~~~~~~~~~~~~~~~~~~~~~~~~~~

:: 

  template <typename CoalescingGen, typename Routing>
  struct per_thread_cache_generator;

Same as `Cache Generator`_ except that the cache is kept per thread, eliminating the need for synchronization.

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

.. _Object-Addressable Message Type: message_types.html#object-based-concept
.. _message type concepts: message_types.html#concepts
.. _coalescing generator: tbd
.. _object-based addressing: tbd
.. _routing: tbd
.. _caching: tbd
