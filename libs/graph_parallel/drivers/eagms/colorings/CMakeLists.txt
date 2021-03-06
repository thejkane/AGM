# This utility macro determines whether a particular string value
# occurs within a list of strings:
#
#  list_contains(result string_to_find arg1 arg2 arg3 ... argn)
# 
# This macro sets the variable named by result equal to TRUE if
# string_to_find is found anywhere in the following arguments.
macro(list_contains var value)
  set(${var})
  foreach (value2 ${ARGN})
    if (${value} STREQUAL ${value2})
      set(${var} TRUE)
    endif (${value} STREQUAL ${value2})
  endforeach (value2)
endmacro(list_contains)

# The PARSE_ARGUMENTS macro will take the arguments of another macro and
# define several variables. The first argument to PARSE_ARGUMENTS is a
# prefix to put on all variables it creates. The second argument is a
# list of names, and the third argument is a list of options. Both of
# these lists should be quoted. The rest of PARSE_ARGUMENTS are
# arguments from another macro to be parsed.
# 
#     PARSE_ARGUMENTS(prefix arg_names options arg1 arg2...) 
# 
# For each item in options, PARSE_ARGUMENTS will create a variable with
# that name, prefixed with prefix_. So, for example, if prefix is
# MY_MACRO and options is OPTION1;OPTION2, then PARSE_ARGUMENTS will
# create the variables MY_MACRO_OPTION1 and MY_MACRO_OPTION2. These
# variables will be set to true if the option exists in the command line
# or false otherwise.
# 
# For each item in arg_names, PARSE_ARGUMENTS will create a variable
# with that name, prefixed with prefix_. Each variable will be filled
# with the arguments that occur after the given arg_name is encountered
# up to the next arg_name or the end of the arguments. All options are
# removed from these lists. PARSE_ARGUMENTS also creates a
# prefix_DEFAULT_ARGS variable containing the list of all arguments up
# to the first arg_name encountered.
MACRO(PARSE_ARGUMENTS prefix arg_names option_names)
  SET(DEFAULT_ARGS)
  FOREACH(arg_name ${arg_names})
    SET(${prefix}_${arg_name})
  ENDFOREACH(arg_name)
  FOREACH(option ${option_names})
    SET(${prefix}_${option} FALSE)
  ENDFOREACH(option)

  SET(current_arg_name DEFAULT_ARGS)
  SET(current_arg_list)
  FOREACH(arg ${ARGN})
    LIST_CONTAINS(is_arg_name ${arg} ${arg_names})
    IF (is_arg_name)
      SET(${prefix}_${current_arg_name} ${current_arg_list})
      SET(current_arg_name ${arg})
      SET(current_arg_list)
    ELSE (is_arg_name)
      LIST_CONTAINS(is_option ${arg} ${option_names})
      IF (is_option)
      SET(${prefix}_${arg} TRUE)
      ELSE (is_option)
      SET(current_arg_list ${current_arg_list} ${arg})
      ENDIF (is_option)
    ENDIF (is_arg_name)
  ENDFOREACH(arg)
  SET(${prefix}_${current_arg_name} ${current_arg_list})
ENDMACRO(PARSE_ARGUMENTS)

macro(pbgl_add_test TESTNAME)
  parse_arguments(PBGL_ADD_TEST
    "ARGS;NUMPROCS"
    ""
    ${ARGN}
    )
  add_executable(${TESTNAME} ${TESTNAME}.cpp)
  target_link_libraries(${TESTNAME}
    ${AM++_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${MPI_LIBRARIES})
  set_target_properties(${TESTNAME}
    PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}"
    COMPILE_FLAGS "${MPI_CXX_INCLUDE_FLAGS}")
  

  if(NOT PBGL_ADD_TEST_NUMPROCS)
    set(PBGL_ADD_TEST_NUMPROCS ${MPIEXEC_MAX_NUMPROCS})
  endif (NOT PBGL_ADD_TEST_NUMPROCS)
  get_target_property(PBGL_ADD_TEST_LOCATION ${TESTNAME} LOCATION)
  string(REGEX REPLACE "\\$\\(.*\\)" "\${CTEST_CONFIGURATION_TYPE}"
        PBGL_ADD_TEST_LOCATION "${PBGL_ADD_TEST_LOCATION}")
  add_test(${TESTNAME}
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PBGL_ADD_TEST_NUMPROCS} 
    ${MPIEXEC_PREFLAGS}
    "${PBGL_ADD_TEST_LOCATION}"
    ${MPIEXEC_POSTFLAGS} ${PBGL_ADD_TEST_ARGS})
endmacro(pbgl_add_test)

macro(pbgl_add_target TESTNAME)
  parse_arguments(PBGL_ADD_TARGET
    "FILENAME;TARG"
    ""
    ${ARGN}
    )
  if(NOT PBGL_ADD_TARGET_FILENAME)
    add_executable(${TESTNAME} ${TESTNAME}.cpp)
  else()
    add_executable(${TESTNAME} ${PBGL_ADD_TARGET_FILENAME}.cpp)
  endif()
  if(PBGL_ADD_TARGET_TARG)
    ADD_DEPENDENCIES(${PBGL_ADD_TARGET_TARG} ${TESTNAME})
  endif()
  target_link_libraries(${TESTNAME}
    ${AM++_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${MPI_LIBRARIES})
  set_target_properties(${TESTNAME}
    PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}"
    COMPILE_FLAGS "${MPI_CXX_INCLUDE_FLAGS}")
endmacro(pbgl_add_target)

pbgl_add_target(eagm_gc_chaotic)

