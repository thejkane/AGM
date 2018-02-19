add_executable(breadth_first_search breadth_first_search.cpp)
target_link_libraries(breadth_first_search
  ${AM++_LIBRARIES} ${MPI_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )

set_target_properties(
  breadth_first_search
  PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}"
             COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

if (BUILD_TESTING)
  add_test(breadth_first_search
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}
    "${EXECUTABLE_OUTPUT_PATH}/\${CTEST_CONFIGURATION_TYPE}/breadth_first_search"
    ${MPIEXEC_POSTFLAGS} 
    "${PBGL_SOURCE_DIR}/libs/graph_parallel/example/weighted_graph.gr")
endif (BUILD_TESTING)

add_executable(shortest_paths shortest_paths.cpp)
target_link_libraries(shortest_paths
   ${AM++_LIBRARIES} ${MPI_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY})

set_target_properties(
  shortest_paths
  PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}"
             COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

if (BUILD_TESTING)
  add_test(shortest_paths
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}
    "${EXECUTABLE_OUTPUT_PATH}/\${CTEST_CONFIGURATION_TYPE}/shortest_paths"
    ${MPIEXEC_POSTFLAGS} 
    "${PBGL_SOURCE_DIR}/libs/graph_parallel/example/weighted_graph.gr")
endif (BUILD_TESTING)
