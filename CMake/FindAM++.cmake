# - Try to find AM++
# Once done this will define
#  AM++_FOUND - System has AM++
#  AM++_INCLUDE_DIRS - The AM++ include directories
#  AM++_LIBRARIES - The libraries needed to use AM++
#  AM++_DEFINITIONS - Compiler switches required for using AM++

set(AM++_DEFINITIONS "")

find_path(AM++_INCLUDE_DIR am++.hpp
	  HINTS ${AM++_INCLUDEDIR} ${AMPP_ROOT}/am++ )

find_library(AM++_LIBRARY NAMES am++ libam++ 
	     HINTS ${AM++_LIBRARYDIR} ${AMPP_ROOT}/src )

set(AM++_LIBRARIES ${AM++_LIBRARY} )
set(AM++_INCLUDE_DIRS ${AM++_INCLUDE_DIR}/ )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set AM++_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(libam++  DEFAULT_MSG
                                  AM++_LIBRARY AM++_INCLUDE_DIR )

mark_as_advanced(AM++_INCLUDE_DIR AM++_LIBRARY )