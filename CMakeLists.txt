##########################################################################
# CMake Build Rules for the Parallel Boost Graph Library                 #
##########################################################################
# Copyright (C) 2011, 2012 Nicholas Edmonds <ngedmond@cs.indiana.edu>    #
# Copyright (C) 2007, 2008 Douglas Gregor <doug.gregor@gmail.com>        #
# Copyright (C) 2007 Troy Straszheim                                     #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################
# Basic Usage:                                                           #
#                                                                        #
#   On Unix variants:                                                    #
#     ccmake PBGL_DIRECTORY                                              #
#                                                                        #
#     (c)onfigure options to your liking, then (g)enerate                #
#     makefiles. Use "make" to build, "make test" to test, "make         #
#     install" to install, and "make package" to build binary            #
#     packages.                                                          #
#                                                                        #
#   On Windows:                                                          #
#     run the CMake GNU, load the Parallel BGL directory, and generate   #
#     project files or makefiles for your environment.                   #
#                                                                        #
# For more information about CMake, see http://www.cmake.org             #
##########################################################################
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)
project(PBGL)

set(PBGL_VERSION 2.0)

# So we can find our custom modules
list(APPEND CMAKE_MODULE_PATH ${PBGL_SOURCE_DIR}/CMake)

# Put the libaries and binaries that get built into directories at the
# top of the build.
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

# Find Boost and the various Boost components we depend on
if(WIN32)
	set(Boost_USE_STATIC_LIBS ON)
endif(WIN32)
find_package(Boost 1.52.0 REQUIRED
  COMPONENTS thread system random)
include_directories(${Boost_INCLUDE_DIRS})

# Find libNBC
find_package(libNBC QUIET)
include_directories(${libNBC_INCLUDE_DIRS})

# Find libCDS
find_package(libCDS REQUIRED)
include_directories(${libCDS_INCLUDE_DIRS})

# Find AM++
find_package(AM++ REQUIRED)
include_directories(${AM++_INCLUDE_DIRS})

# Always include the Parallel BGL source directory at the beginning of
# the include list, since PBGL overrides some Boost header files.
include_directories(BEFORE ${PBGL_SOURCE_DIR})

# Testing support
include(CTest)
if(BUILD_TESTING)
  enable_testing()
  #add_subdirectory(libs/graph_parallel/test)
  add_subdirectory(libs/graph_parallel/drivers)
  add_subdirectory(libs/graph_parallel/drivers/tests)
  add_subdirectory(libs/graph_parallel/drivers/eagms/bfses)
  add_subdirectory(libs/graph_parallel/drivers/eagms/ssspes)
  add_subdirectory(libs/graph_parallel/drivers/eagms/cces)
  add_subdirectory(libs/graph_parallel/drivers/eagms/pageranks)
  add_subdirectory(libs/graph_parallel/drivers/eagms/kcores)
  add_subdirectory(libs/graph_parallel/drivers/eagms/mises)
  add_subdirectory(libs/graph_parallel/drivers/eagms/colorings)
endif(BUILD_TESTING)

# Example programs
set(BUILD_EXAMPLES TRUE CACHE BOOL "Whether to build Parallel BGL examples")
if (BUILD_EXAMPLES)
  add_subdirectory(libs/graph_parallel/example)
endif (BUILD_EXAMPLES)

# Documentation
find_program(RST2HTML
  NAMES rst2html.py rst2html
  DOC "Program that generates HTML from reStructuredText documents"
  )
set(RST2HTML_FLAGS "-t --no-xml-declaration"
    CACHE STRING "Extra flags to pass to RST2HTML")
mark_as_advanced(RST2HTML RST2HTML_FLAGS)

if (RST2HTML)
  option(BUILD_DOCUMENTATION "Whether to build documentation for the Parallel BGL" ON)
  if (BUILD_DOCUMENTATION)
    add_subdirectory(libs/graph_parallel/doc)
  endif (BUILD_DOCUMENTATION)
elseif (BUILD_DOCUMENTATION)
  message("error: rst2html is required to build Parallel BGL documentation")
endif ()

##########################################################################
# Installation                                                           #
##########################################################################
install(DIRECTORY boost 
        DESTINATION "include/pbgl-${PBGL_VERSION}"
        REGEX ".svn" EXCLUDE)

configure_file("${PBGL_SOURCE_DIR}/CMake/pbgl-config.cmake.in"
  "${PBGL_BINARY_DIR}/pbgl-config.cmake" 
  @ONLY)
install(FILES "${PBGL_BINARY_DIR}/pbgl-config.cmake"
  DESTINATION lib/pbgl-${PBGL_VERSION})
