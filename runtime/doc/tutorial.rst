.. Copyright (C) 2013 The Trustees of Indiana University.
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

===============
|Logo| Tutorial
===============

.. contents::

This tutorial gives a roundup of the AM++ interface.  The tutorial is based on the code in the `tutorial.cpp`_ file that is part of the AM++ distribution (use ``make tutorial`` in the distribution directory to compile it).  The goal is an overview of the interface, and the program does not do anything interesting.

First, one has to create an environment::

  int main(int argc, char** argv) {
  (void)argc; (void)argv;

  bool use_threads = true;
  environment env = mpi_environment(argc, argv, use_threads);

Environment takes the arguments ``argv`` and the flag ``use_threads`` that indicates whether threading is desired or not.  The environment sets up MPI as necessary, and if threading is requested it checks whether threading is supported by the MPI implementation.  Furthermore, it registers built-in datatypes (the programmer needs to register their own datatypes, see `MPI Datatypes`_.  When the environment object is destroyed, it cleans up the underlying MPI environment as necessary.  Next, we need to create a *transport*::

  transport trans = env.create_transport();

Transport contains all the information that AM++ needs to send and receive messages, to manage threading, and so on.  For example, transport can be used to obtain the current rank and the total number of ranks::

  transport::rank_type r = trans.rank(), s = trans.size();

Next, we will create some *message types* that we can then use to send messages.  The simplest message type is created directly by the transport::

  message_type<int> mt = trans.create_message_type<int>();
  mt->set_max_count(100);
  mt.set_handler(my_handler());

The message type ``mt`` can be used to send and receive values of type ``int``.  The ``max_count`` specifies the maximum number of messages to be received at one time (it cannot be 0).  Finally, ``set_handler`` is used to specify a function that will *handle* the incoming messages, and``my_handler`` is defined as follows::

  struct my_handler {
    // Raw handler
    void operator()(int source, const int* buf, int count) const {};
  };

When the handler is called by AM++, it receives the source rank of the messages ``source``, the data for the messages ``buf``, and the ``count`` of the received messages which is also the size of ``buf``.  The count is less or equal to the specified ``max_count``.  AM++ provides more complex message types.  For example, coalescing message type is created as follows::

    basic_coalesced_message_type<int, my_handler> mt2(basic_coalesced_message_type_gen(100), trans);
    mt2.set_handler(my_handler());

The basic coalesced message provides a simple coalescing layer.  The ``basic_coalesced_message_type_gen`` is a rather long-named way to specify the maximum count of coalesced messages.  In this case, it is probably an overkill, but it conforms to the general convention of AM++ where a message type takes a generator as the first argument, transport as the second argument, and additional, message type specific arguments after that.  Coalesced messages require extending ``my_handler`` with a new operator::

  void operator()(int source, int data) const {};

Coalesced message type puts many small messages into a one large message, and it handles ``max_count`` from the basic message type itself.  When it receives the large coalesced message, it invokes the handler function for each packed data (in this case an ``int``).  Because of that, the new handler has fewer arguments than a handler for a basic message type.

Both the raw and the coalesced message types require direct addressing (note that handlers receive source rank as one of the arguments).  AM++ provides *object-based addressing* where the destination rank of a message is inferred from the data in the message.  Object-based addressing (and many other features) can be added to messages using *message type generators*.  The simplest generator that adds object-based addressing is the aptly named ``simple_generator``::

  typedef simple_generator<basic_coalesced_message_type_gen> simple_gen;
  simple_gen::
    call_result<int, my_handler, owner_map_type, no_reduction_t>::type
    mt3(simple_gen(basic_coalesced_message_type_gen(100)),
        trans, owner_map, no_reduction);

A generator template, in general, takes some underlying generator type and extra arguments.  In the case of the ``simple_generator`` there is only one argument, ``basic_coalesced_message_type_gen``, that provides the underlying message type.  Every generator provides a member template ``call_result`` that takes the arguments necessary to create the message type and statically generates the resulting message type member called, of course, ``type``.  In the case of ``simple_generator``, ``call_result`` takes the data type sent in messages, ``int``, the handler, an *owner map*, and a *reduction*.  Owner map is an object that supports the `Boost Property Map`_ interface::

  struct owner_map_type {};

  template<typename T>
  amplusplus::transport::rank_type 
  get(const owner_map_type, const T) { return 0; } // something more complicated 

This owner map just demonstrates the interface, always returning rank 0 (not a good idea in general).  ``simple_generator`` also takes a reduction, but it is only to fit the general interface of generators.  Since it does not do anything with the reduction, we can pass in ``no_reduction`` (an object provided by AM++) without any effect on functionality.  Because object-based addressing infers the address from the message type itself, the handler interface for object-addressed messages only takes the data itself (``int`` in this case)::

  void operator()(int data) const {};  

Another layer that we can add to a message type is routing.  The idea behind routing is to use intermediate hops instead of sending messages directly to their destination.  One reason to do this is to decrease the number of buffers necessary to send messages.  With routing, each rank must keep buffers only for its intermediate neighbors in the routing scheme, not for all the other ranks as must be done without routing.  Routing can be layered on top of a coalescing message generator with ``routing_generator``::

  typedef routing_generator<basic_coalesced_message_type_gen, ring_routing> routing_gen;
  routing_gen::
    call_result<int, my_handler, owner_map_type, no_reduction_t>::type
    mt4(routing_gen(basic_coalesced_message_type_gen(100), ring_routing(trans.rank(), trans.size())), trans, owner_map, no_reduction);

In this example, ``routing_generator`` is instantiated with ``basic_coalesced_message_type_gen`` and ``ring_routing`` where each rank has two neighbors with all the ranks arranged into a ring.  Other examples of routing include ``rook_routing`` (think about a rook on a chessboard finding the shortest route between any two squares) and ``hypercube_routing`` (see Hypercube_).

Caching is another very useful feature in AM++.  A caching message type can be created with one of the caching generators, for example, with the ``per_thread_cache_generator``::

  typedef per_thread_cache_generator<amplusplus::counter_coalesced_message_type_gen, amplusplus::no_routing> cache_gen;
  cache_gen::call_result<std::pair<int, int>, my_handler, owner_map_type, idempotent_combination_t<max_value> >::type mt5(cache_gen(counter_coalesced_message_type_gen(1024), 20, no_routing(trans.rank(), trans.size())), trans, owner_map, idempotent_combination(max_value()));

``per_thread_cache_generator`` equips each thread with a cache of the previously sent messages.  In this particular case, the cache is based on an idempotent_ operation.  Because the operation is idempotent, the cache can be an optimized `write-through`_ cache: a message is sent if the reduction produced a value different than the one in cache and is ignored otherwise.  Note that the data in this message type is changed to ``std::pair<int, int>``.  The default instantiation of idempotent_combination extracts the first ``int`` in the pair as a key for the cache, and the second ``int`` as a value.  So, for example, two pairs ``(3,4)`` and ``(5,4)`` are unrelated as far as the combination is concerned, and both would be written through cache, assuming that the cache was empty.  Two values ``(3,4)`` and ``(3,2)`` are related because of the same key (``3`` for the first element of the pair), and if ``(3,4)`` was in cache, ``(3,2)`` would not get written through cache because :math:`2<4` (``max_value`` chooses the greater value).

The message types we created can now be used to actually send messages::

  shared_array<int> buf(new int[3]);
  {
    scoped_epoch e(trans);
    mt.message_being_built(0);
    mt.send(buf.get(), 3, dest, empty_deleter());
    mt2.send(1, dest);
    mt3.send(2);
    mt4.send(3);
    mt5.send(std::make_pair(1,1));
  } 

All AM++ messages must be sent within an *epoch*.  One way to create an epoch is with the ``scoped_epoch`` object which employs the RAII_ technique to acquire an epoch at the time of initialization and destroy an epoch at the end of the scope.  An epoch contains all the information about all sent messages, and it allows AM++ to perform tasks such as termination detection.  When the ``scoped_epoch`` object ``e`` is destroyed at the end of the scope, it invokes termination detection and it waits until no more messages are to be sent or delivered.  The important rule in AM++ is that no messages can be sent or received outside of an epoch.  

Note that non-coalesced message type needs a *deleter* in the send call.  The deleter is invoked by AM++ to handle potential deletion of the input buffer.  Such interface allows "fire and forget" sends with dynamically allocated buffers.  In our case, we have allocated the buffer with a ``shared_array``, so deletion is already taken care of, and we can use a "fake" ``empty_deleter``.

This concludes the tutorial on the basic interface of AM++.  For complete reference, see the complete documentation.

----------------------------------------------------------------------

Copyright (C) 2009-2013 The Trustees of Indiana University.

:Authors: 
          Jeremiah Willcock, 
	  Marcin Zalewski, 
	  and Andrew Lumsdaine

.. _tutorial.cpp: tutorial.cpp
.. _MPI Datatypes: mpi_datatypes.html
.. _Boost Property Map: http://www.boost.org/doc/libs/1_54_0/libs/property_map/doc/property_map.html
.. _Hypercube: http://en.wikipedia.org/wiki/Hypercube
.. _idempotent: http://en.wikipedia.org/wiki/Idempotence
.. _write-through: http://en.wikipedia.org/wiki/Cache_%28computing%29#Writing_policies
.. _RAII: http://en.wikipedia.org/wiki/Resource_Acquisition_Is_Initialization
.. |Logo| image:: ampp-logo.png
            :align: middle
            :alt: AM++
            :target: http://crest.iu.edu/research/am++
