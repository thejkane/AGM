# Abstract Graph Machine and Distributed, Shared-Memory Parallel Graph Kernels

### Overview
Abstract Graph Machine (AGM) models orderings in asynchronous parallel graph algorithms. The AGM model expresses a graph algorithms as a function (AKA processing function) and an ordering (strict weak ordering relation). This repository contains an implementation of the AGM model, set of graph kernels implemented using the AGM model. In addition to AGM graph kernels, repository also, contains few distributed, shared-memory parallel graph kernels that does not use the AGM model. For distributed communication, implementation uses MPI and a MPI based Active Messaging framework -- AM++. All implementations are in C++ and we make use of heavy template meta-programming, therefore, compilation times are quite high and execution include minimum overhead.

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

The extended AGM explores spatial and temporal ordering and derives less synchronous distributed, shared-memory parallel graph algorithms than AGM algorithms. For the EAGM orderings can be peformed at global memory level, node memory level, NUMA memory or at the thread memory level. A description of EAGM achieves spatial and temporal ordering is shown in Figure~YY.


### Graph Kernels Available in AGM/EAGM Model
The AGM/EAGM graph processing framework is implemented as part of Parallel Boost Graph Library, version 2 (PBGL2). Graph structure definitions are based on the graph structure definitions provided by the PBGL2 (Compressed Sparse Row and Adjacency List). Further, AGM/EAGM model uses 1D graph distributions.

1. Breadth First Search -- With a single processing function, we can achieve multiple algorithms by changing ordering. E.g., The chaotic BFS does not perform  any ordering, but Level synchronous BFS performs ordering by level.
2. Single Source Shortest Path -- Multiple algorithms can be derived by changing orderings. This includes Delta-Stepping, KLA, Bellman-Ford, Distributed-Control.
3. Connected Components
4. Maximal Independent Set (MIS) -- The FIX MIS.
5. Graph Coloring -- The original processing function is based on Jones-Plassman Graph Coloring.
6. k-Core Decomposition
7. PageRank
8. Triangle Counting (In Progress)

```cpp
typedef std::tuple<Vertex, Level> WorkItem;
```

### Other Graph Kernels (Non AGM/EAGM Graph Kernels)
In addition to AGM/EAGM graph kernels we have number of different graph kernels that does not use the AGM/EAGM model. They are as follows:

+ Triangle Counting (Authors : Thejaka Amila Kanewala, Andrew Lumsdaine)
  + Striped Triangle Counting
  + Blocked Triangle Counting
  + Traversal based Triangle Counting

+ Maximal Independent Set (Authors : Thejaka Amila Kanewala, Andrew Lumsdaine)
  + FIX
  + FIX-Bucket
  + FIX-PQ
  + Luby A
  + Luby B

+ Connected Components
  + Traversal Based Connected Components (Authors : Thejaka Amila Kanewala, Andrew Lumsdaine)
  + Priority Connected Components (Authors : Thejaka Amila Kanewala, Andrew Lumsdaine)
  + Delta Based Connected Components (Authors : Thejaka Amila Kanewala, Andrew Lumsdaine)
  + Shiloach-Vishkin Connected Components (Authors : Nick Edmonds, Andrew Lumsdaine)

+ Breadth First Search
  + Chaotic Breadth First Search (Authors : Thejaka Amila Kanewala, Andrew Lumsdaine)
  + Level-Synchronous Breadth First Search (Authors : Nick Edmonds, Andrew Lumsdaine)

+ Single Source Shortest-Paths
  
1. Level Synchronous Breadth First Search (Authors : Nick Edmonds, Andrew Lumsdaine)
2. Optimized Shiloach-Vishkin Connected Components (Authors : Nick Edmonds, Andrew Lumsdaine)
3. Delta-Stepping Shortes Paths (Authors : Nick Edmonds, Andrew Lumsdaine)
