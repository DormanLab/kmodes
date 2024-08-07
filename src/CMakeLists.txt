cmake_minimum_required (VERSION 3.5.0)

# set the project name
project(kmodes VERSION 1.0)

# specify the C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

# compile with debugging symbols
SET(GCC_COVERAGE_COMPILE_FLAGS "-g")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

# compile with optimization
#SET(GCC_COVERAGE_COMPILE_FLAGS "-O3")
#add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

add_definitions(-DMATHLIB_STANDALONE)

include(FindPackageHandleStandardArgs)

#aux_source_directory(. DIR_SRCS)
set(SRCS cluster.c cmdline.c error.c io.c io_kmodes.c
	kmodes.c matrix_exponential.c order.c r8slib.c
	run_kmodes.c sample.c timing_mach.c entropy.c hash.c
	rmath/RNG.c rmath/rgamma.c rmath/sunif.c rmath/rexp.c
	rmath/rnorm.c rmath/sexp.c rmath/snorm.c rmath/qnorm.c
	rmath/fmin2.c rmath/fmax2.c)


# following modified from RStudio's FindLibR.cmake to find R include and lib path 
# https://github.com/rstudio/rstudio/blob/master/cmake/modules/FindLibR.cmake 

#if(APPLE)
#
#   find_library(LIBR_LIBRARIES R)
#   
#   if(LIBR_LIBRARIES MATCHES ".*\\.framework")
#   set(LIBR_HOME "${LIBR_LIBRARIES}/Resources" CACHE PATH "R home directory")
#      set(LIBR_INCLUDE_DIRS "${LIBR_HOME}/include" CACHE PATH "R include directory")
#      set(LIBR_LIB_DIR "${LIBR_HOME}/lib" CACHE PATH "R lib directory")
#      set(LIBR_EXECUTABLE "${LIBR_HOME}/bin/R" CACHE PATH "R executable")
#   else()
#      get_filename_component(_LIBR_LIBRARIES "${LIBR_LIBRARIES}" REALPATH)
#      get_filename_component(_LIBR_LIBRARIES_DIR "${_LIBR_LIBRARIES}" PATH)
#      set(LIBR_EXECUTABLE "${_LIBR_LIBRARIES_DIR}/../bin/R")
#      execute_process(
#         COMMAND ${LIBR_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
#                   OUTPUT_VARIABLE LIBR_HOME
#      )
#      set(LIBR_HOME ${LIBR_HOME} CACHE PATH "R home directory")
#      set(LIBR_INCLUDE_DIRS "${LIBR_HOME}/include" CACHE PATH "R include directory")
#      set(LIBR_LIB_DIR "${LIBR_HOME}/lib" CACHE PATH "R lib directory")
#   endif()
#
#else()  # UNIX
#
## find executable
#      find_program(LIBR_EXECUTABLE R)
#      if(LIBR_EXECUTABLE-NOTFOUND)
#         message(STATUS "Unable to locate R executable")
#      endif()
#
## ask R for the home path
#      if(NOT LIBR_HOME)
#         execute_process(
#            COMMAND ${LIBR_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
#                      OUTPUT_VARIABLE LIBR_HOME
#         )
#         if(LIBR_HOME)
#           set(LIBR_HOME ${LIBR_HOME} CACHE PATH "R home directory")
#         endif()
#      endif()
#
# # ask R for the include dir
#      if(NOT LIBR_INCLUDE_DIRS)
#         execute_process(
#            COMMAND ${LIBR_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('include'))"
#            OUTPUT_VARIABLE LIBR_INCLUDE_DIRS
#         )
#         if(LIBR_INCLUDE_DIRS)
#           set(LIBR_INCLUDE_DIRS ${LIBR_INCLUDE_DIRS} CACHE PATH "R include directory")
#         endif()
#      endif()
#
#  # ask R for the lib dir
#      if(NOT LIBR_LIB_DIR)
#         execute_process(
#            COMMAND ${LIBR_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('lib'))"
#            OUTPUT_VARIABLE LIBR_LIB_DIR
#         )
#          if(LIBR_LIB_DIR)
#           set(LIBR_LIB_DIR ${LIBR_LIB_DIR} CACHE PATH "R lib directory")           
#         endif()
#      endif()
#
#
#   # sometimes lib and include are not under R.home()
#   get_filename_component(_LIBR_EXECUTABLE "${LIBR_EXECUTABLE}" PATH)
#   get_filename_component(_LIBR_EXECUTABLE_DIR "${_LIBR_EXECUTABLE}" PATH)
#   #message(STATUS "Found R upper directory: ${_LIBR_EXECUTABLE_DIR}")
#
#endif()
#
#message(STATUS "Found R: ${LIBR_EXECUTABLE}")
#message(STATUS "Found R home: ${LIBR_HOME}")
#message(STATUS "Found R include directory: ${LIBR_INCLUDE_DIRS}")
#message(STATUS "Found R lib directory: ${LIBR_LIB_DIR}")

include_directories(${LIBR_INCLUDE_DIRS})
link_directories(${LIBR_LIB_DIR})

# compile option 
add_compile_options(-D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500)
#add_compile_options(-Wall -Wextra -pedantic)
add_executable(run_kmodes ${SRCS})

# Rmath
#find_library(RMATH_LIBRARY NAMES Rmath  HINTS ${LIBR_LIB_DIR} ${_LIBR_EXECUTABLE_DIR}/lib ${_LIBR_EXECUTABLE_DIR}/lib64 .)
#find_package_handle_standard_args(LibRmath DEFAULT_MSG RMATH_LIBRARY)
#if(RMATH_LIBRARY)
#   target_link_libraries(run_kmodes ${RMATH_LIBRARY})
#else()
#   message(FATAL_ERROR "Failed to find libRmath !")
#endif()

# math
TARGET_LINK_LIBRARIES(run_kmodes -lm)
