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

=================
|Logo| Reductions
=================

.. contents::

Overview
--------

In AM++, reductions are used in caching of messages.  Reductions, in general, decide how to combine two values, one in cache and one that is being written.

Generators__ that produce caching message types must "know" how to deal with a given reduction.  AM++ provides a set of built-in reductions that support the most common caching scenarios.

__ generators.html

Reduction Concept
-----------------

Notation
~~~~~~~~

Reduction
  The reduction type

Valid Expressions
~~~~~~~~~~~~~~~~~

+--------------------+---------------------------------+--------------------+
|Check if a type is a|``is_valid_reduction<Reduction>``|This template must  |
|valid reduction     |                                 |be `Boost MPL true  |
|                    |                                 |or false`__.        |
+--------------------+---------------------------------+--------------------+

__ http://www.boost.org/doc/libs/1_55_0/libs/mpl/doc/refmanual/bool.html


Models
------

Void Reduction
~~~~~~~~~~~~~~

:: 

  typedef *implementation-defined* no_reduction_t;

This reduction is provided for the generators where no reduction is necessary.  If it is provided to a caching generator, then the resulting message type will not cache sent messages.

Duplicate Removal Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  template <typename Project, typename PolicyTag = make_default_policy>
  struct duplicate_removal_t;

Duplicate removal policy removes duplicate messages.  

``Project``
  The type of the projection that extracts a value from messages.  The values are compared in duplicate removal.  Note that the result of projection must support `boost::hash`__.

__ www.boost.org/libs/doc/html/boost/hash.html

``PolicyTag``
  The type describing the policy for comparison of values.

Two convenience functions are provided for producing values of ``duplicate_removal_t``.  The first function produces the default policy removal, and the second can be used with custom policies::

  template <typename Project> duplicate_removal_t<Project> duplicate_removal(const Project& p);
  template <typename Project, typename PolicyTag> duplicate_removal_t<Project, PolicyTag> duplicate_removal(const Project& p, const Policy& po);

The default policy compares the values extracted by the projection using the ``==`` operator.    

.. note::
  The default policy requires that the projection function's result be constructible from values of ``size_t``.

.. rubric:: Constructors

::

  explicit duplicate_removal_t(const Project& p);
  duplicate_removal_t(const Project& p, const Policy& po);

.. rubric:: Custom Policies

To provide a custom policy, a few steps are necessary:

1. A policy must be defined.  The policy must provide the following members:

``stored_type``
  The type of the values stored in the cache.

``stored_type stored_value(const Arg& d) const``
  A function which converts values of type ``Arg`` (the type of message date) into values of type ``stored_value``.

.. note:: 
  Step 3 describes how policies receive the ``Arg`` type.  

``bool match(const Arg& d, stored_type stored) const``
  A function matching stored values to values of type ``Arg``.

``size_t hash(stored_type v, size_t size) const``
  A function producing hashes for stored values.  Hashes produced must be smaller than ``size``.

``stored_type dummy_value(int h, size_t size) const``
  This function produces values of stored_type that do not hash to ``hash(h, size)``.


2. A tag must be defined.  The tag is a type that will be used to identify the policy.

3. The ``get_policy_type`` must be specialized for the new policy.  It is defined like this:

::

  template <typename Project, typename Policy, typename Arg>
  struct get_policy_type {
    typedef Policy type; 
    static type get_policy(const Project&, const Policy& p) {return p;}
  };

For the default policy, it is specialized as follows::

  struct get_policy_type<Project, make_default_policy, Arg> {
    typedef duplicate_policy_from_projection<Project, Arg> type; 
    static type get_policy(const Project& pr, const make_default_policy&) {return type(pr);}
  };

Combination Reduction
~~~~~~~~~~~~~~~~~~~~~

::

  template <typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
  struct combination_t;

Combination reduction combines values.  

``Combine``
  The operation used to combine values.

``OptIdentity``
  Optional identity type.  If identity type is specified, it is used to stop identity values from being sent.  If this is not desired, ``no_identity_t`` can be used to always send values.

``GetKey``
  The type of the function object used to extract keys from message data type.  By default, the first element of a pair is used as the key.   Note that the result of this function object must support `boost::hash`__.

__ www.boost.org/libs/doc/html/boost/hash.html

``GetValue``
  The type of the function object used to extract values from message data type.  By default, the second element of a pair is used as the value.

``MakeKeyval``
  The type of the function object used to combine keys and values into message data.  By default, the key and value are combined into a pair.

Three convenience functions are provided for making combination reductions::

  template <typename Combine>
  combination_t<Combine, no_identity_t, pair_first, pair_second, make_pair_t>
  combination(const Combine& c);

  template <typename Combine, typename Identity>
  combination_t<Combine, Identity, pair_first, pair_second, make_pair_t>
  combination(const Combine& c, const Identity& i);

  template <typename Combine, typename Identity, typename GetKey, typename GetValue, typename MakeKeyval>
  combination_t<Combine, Identity, GetKey, GetValue, MakeKeyval>
  combination(const Combine& c, const Identity& i, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval);

.. rubric:: Constructor

::

  combination_t(const Combine& combine, const OptIdentity& opt_identity, const GetKey& get_key, const GetValue& get_value, const MakeKeyval& make_keyval);

Idempotent Combination Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  template <typename Combine, typename OptIdentity = no_identity_t, typename GetKey = pair_first, typename GetValue = pair_second, typename MakeKeyval = make_pair_t>
  struct idempotent_combination_t;

Idempotent combination has an interface identical to combinations__.

__ #combination-reduction

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
