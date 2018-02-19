# - Try to find libNBC
# Once done this will define
#  LIBNBC_FOUND - System has libNBC
#  LIBNBC_INCLUDE_DIRS - The libNBC include directories
#  LIBNBC_LIBRARIES - The libraries needed to use libNBC
#  LIBNBC_DEFINITIONS - Compiler switches required for using libNBC

set(libNBC_DEFINITIONS "")

find_path(LIBNBC_INCLUDE_DIR nbc.h
	  HINTS ${LIBNBC_INCLUDEDIR} ${LIBNBC_ROOT}/include )

find_library(LIBNBC_LIBRARY NAMES nbc libnbc 
	     HINTS ${LIBNBC_LIBRARYDIR} ${LIBNBC_ROOT}/lib )

set(libNBC_LIBRARIES ${LIBNBC_LIBRARY} )
set(libNBC_INCLUDE_DIRS ${LIBNBC_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBNBC_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(libNBC  DEFAULT_MSG
                                  LIBNBC_LIBRARY LIBNBC_INCLUDE_DIR )

mark_as_advanced(LIBNBC_INCLUDE_DIR LIBNBC_LIBRARY )