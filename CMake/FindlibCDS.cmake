# - Try to find LIBCDS
# Once done this will define
#  LIBCDS_FOUND - System has LIBCDS
#  LIBCDS_INCLUDE_DIRS - The LIBCDS include directories
#  LIBCDS_LIBRARIES - The libraries needed to use LIBCDS
#  LIBCDS_DEFINITIONS - Compiler switches required for using LIBCDS

set(libCDS_DEFINITIONS "")

find_path(LIBCDS_INCLUDE_DIR init.hpp
	  HINTS ${LIBCDS_INCLUDEDIR} ${LIBCDS_ROOT} )

find_library(LIBCDS_LIBRARY NAMES cds libcds 
	     HINTS ${LIBCDS_LIBRARYDIR} ${LIBCDS_ROOT}/lib )

set(libCDS_LIBRARIES ${LIBCDS_LIBRARY} )
set(libCDS_INCLUDE_DIRS ${LIBCDS_INCLUDE_DIR}/ )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCDS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(libCDS  DEFAULT_MSG
                                  LIBCDS_LIBRARY LIBCDS_INCLUDE_DIR )

mark_as_advanced(LIBCDS_INCLUDE_DIR LIBCDS_LIBRARY )
