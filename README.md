# Abstract Graph Machine and Distributed, Shared-Memory Parallel Graph Kernels

### Overview
Abstract Graph Machine (AGM) models orderings in asynchronous parallel graph algorithms. The AGM model expresses a graph algorithms as a function (AKA processing function) and an ordering (strict weak ordering relation). This repository contains an implementation of the AGM model, set of graph kernels implemented using the AGM model. In addition to AGM graph kernels, repository also, contains few distributed, shared-memory parallel graph kernels that does not use the AGM model.

#### Authors of AGM : Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine

## Introduction

Figure~XX shows a very high-level overview of the Abstract Graph Machine. A graph kernel in AGM is expressed as a function and an ordering. Function takes a WorkItem as an input and produces zeor or more WorkItems. The definitiion of a WorkItem can be based on vertices or edges. More concretely an AGM for a particular graph kernel is defined with the following:

1. A definition of a Graph,
2. A definition of a WorkItem,
3. A definition of a set of states,
4. A definition of a processing function,
5. A set of initial WorkItems.
6. A definition of a strict weak ordering relation on WorkItem.

More detials of the abstract model can be found in [1], [2], [3].

### Extended Abstract Graph Machine (EAGM)
