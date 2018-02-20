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

====================
|Logo| AM++ Concepts
====================

Architecture
============

.. figure:: AM++Architecture.png
   :alt: AM++ architecture
   :width: 400px
   :align: right

   **Figure 1:** AM++ architecture.

AM++ builds on an underlying communication layer.  The communication layer is encapsulated in a *transport*.  Multiple transports can be created.  Transports can also be cloned, creating transports with the same parameters.  Transports are designed to send messages produced by the higher-level components.  Transports provide MPI-style progress that is guaranteed only in calls to AM++.

Message types are registered with a transport.  Each message type is registered with the corresponding handler.  Message types have a static and a dynamic component: the C++ type of a message type is created statically, and the dynamic components are initialized dynamically, in constructors.  Such design allows AM++ to provide type safety and reduce casting.  Message types handle buffers automatically.  A sender-provided function is called when a buffer used in sending can be reused, and the handler buffers are created and managed by AM++.

Since AM++ handlers can send messages, AM++ provides *termination detection*.  Termination detection is handled in *epochs*, discussed in the next section.

AM++ provides message *coalescing*, *reductions* (cache), and *routing* built on the basic message types.  

Message coalescing may improve performance of small messages by putting many small messages into a larger buffer, that is then sent using the underlying message type.  AM++ automatically provides a handler that calls the user-provided handler on every small message in the large coalesced buffer.  Coalescing is "smart" in the sense that AM++ uses heuristics to decide when to send coalesced buffers, and it may send them when they contain as little as only one small message.

In some applications messages can be combined in one or another.  In AM++, this is achieved using reductions__.  There are 3 types of reductions supported by the AM++ generators__:

 - `Duplicate removal`__ should be used when sending the same message more then once does not change the result.
 - `Combination reduction`__ should be used when values can be combined before send, for example, using addition or multiplication.
 - `Idempotent combination reduction`__ should be used when repeatedly combining with the same value does not change the result.  Max and min functions are examples of idempotent reductions.

__ reductions.html
__ generators.html
__ reductions.html#duplicate-removal-reduction
__ reductions.html#combination-reduction
__ reductions.html#idempotent-combination-reduction

Finally, coalesced messages can be *routed*.  Routing is useful because coalesced messages must maintain a potentially large buffer for every destination.  With routing, the number of buffers can be reduced by adding intermediate hops.  Rook routing, for example, arranges ranks into a chessboard-like square.  To send a message to a rank in a different column than the origin rank, the message is first forwarded in the row of the chessboard, and then sent up in the column (think of a rook moving on the chessboard).  Routing can be combined with caching, providing an additional opportunity to apply reductions for messages-in-flight.

Summary of AM++ Components
==========================

.. _`Figure 2`:
.. figure:: pebble-model.png
  :alt: Summary of the AM++ model
  :width: 600px
  :align: right

  **Figure 2:** Summary of the AM++ model

`Figure 2`_ summarizes the AM++ components.  The example is based on a hypothetical insertions into a table, where each insertion results in an AM++ message that contains the cell number and the element to insert.  There are 4 nodes, P0, P1, P2, and P3.  Each node runs a multithreaded program (represented by the white cloud at node P0).  Local messages are handled immediately with an optional self-send check, and the ``table.insert(0, F)`` is handled immediately without any messages being sent.  The nodes are arranged into a hypercube, where each node can only send messages to its immediate neighbors (in this case).  Because of routing, AM++ can perform multi-source coalescing as done at node P1 with messages ``(6, D)`` and ``(6, A)`` reduced into the message ``(6, D A)``.  Node P3 sends two messages to the same table cell at node P1, and the messages can be reduced immediately, at node P3.  Node P2 sends two messages to node P3, to different table cells.  The messages are coalesced but not reduced (because of different destination cells).

Lifetime of various components.
===============================

.. figure:: AM++Runtime.png
   :alt: AM++ architecture
   :width: 600px
   :align: right

   **Figure 3:** The lifetime of various AM++ components

All messages in AM++ require a transport to be sent.  Furthermore, messages must be sent in *epochs*.  To send messages, an epoch must be started on all ranks involved in a transport.  Messages are sent and received when AM++ progress occurs, which is during the calls to the AM++ interface in the current implementation.  Since message handlers can spawn further messages, an epoch enters *termination detection* when end of the epoch is requested on a given node.  Additional messages can be sent and received during termination detection.  Epoch on any node will terminate only when all messages have been sent and received on all nodes.

Epochs cannot be nested.  All message types must be registered outside of epochs.  Furthermore, the lifetime of a transport must encompass all the epochs run in that transport.

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
