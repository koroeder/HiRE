###########################################################
#
# CMake file for MD HiRE
# 
############################################################

cmake_minimum_required(VERSION 3.0)

# Silence warnings about PROJECT_VERSION not being set
cmake_policy(SET CMP0048 OLD)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(HIREMD)

set(PROJECT_VERSION "1.0")
set(PROJECT_CONTACT "kr366@cam.ac.uk")
get_filename_component(GIT_ROOT ${HIREMD_SOURCE_DIR}/../.. ABSOLUTE)
message(WARNING "Root directory is ${HIREMD_SOURCE_DIR}")
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(CMAKE_C_FLAGS_RELEASE "-O1 -DNDEBUG")
endif (CMAKE_C_COMPILER_ID MATCHES "Intel")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${GIT_ROOT}/CMakeModules)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Release Debug Debug_slow RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)

# Compiler blacklist. Feel free to add more as required.
if (${COMPILER_SWITCH} STREQUAL "gfortran")
    # Get the version number
    execute_process(COMMAND "${CMAKE_Fortran_COMPILER}" "--version" OUTPUT_VARIABLE GCC_VERSION)
    string(REGEX MATCH "[0-9]+(\\.[0-9]+)+" CMAKE_Fortran_COMPILER_VERSION ${GCC_VERSION})
    # Require at least 4.7.4 
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.7.4)
       message(WARNING "If you have run a 'module load', the configuration must be wiped and rerun.")
       message(FATAL_ERROR "GCC version must be at least 4.7.4 but ${CMAKE_Fortran_COMPILER_VERSION} is being used.")
    endif()
endif()
if (${COMPILER_SWITCH} STREQUAL "nag")
    # Get the version number
    execute_process(COMMAND "${CMAKE_Fortran_COMPILER}" "-V" ERROR_VARIABLE NAG_VERSION)
    string(REGEX MATCH "[0-9]+(\\.[0-9]+)*" NAG_MAJOR_VERSION ${NAG_VERSION})
    string(REGEX MATCH "[0-9][0-9][0-9]+" NAG_MINOR_VERSION ${NAG_VERSION})
    set(CMAKE_Fortran_COMPILER_VERSION "${NAG_MAJOR_VERSION}.${NAG_MINOR_VERSION}")
    # We need nagfor to be at least 6.0 to include features from the 2008 standard.
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.0.1028)
       message(WARNING "If you have run a 'module load', the configuration must be wiped and rerun.")
       message(FATAL_ERROR "nagfor version must be at least 6.0.1028 to include ERFC from the Fortran 2008 standard.")
    endif()
    # Ban 6.1.6109 as it has a problem with its preprocessor 
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 6.1.6109)
       message(WARNING "If you have run a 'module load', the configuration must be wiped and rerun.")
       message(FATAL_ERROR "nagfor version 6.1.6109 is disallowed due to a known preprocessor bug.")
    endif()
endif()
if (${COMPILER_SWITCH} STREQUAL "ifort")
    # Get the version number
    execute_process(COMMAND "${CMAKE_Fortran_COMPILER}" "--version" OUTPUT_VARIABLE IFORT_VERSION)
    string(REGEX MATCH "[0-9]+(\\.[0-9]+)+" CMAKE_Fortran_COMPILER_VERSION ${IFORT_VERSION})
    #  No current requirements 
endif()
if (${COMPILER_SWITCH} STREQUAL "pgi")
    # Get the version number
    execute_process(COMMAND "${CMAKE_Fortran_COMPILER}" "--version" OUTPUT_VARIABLE PGI_VERSION)
    string(REGEX MATCH "[0-9]+(\\.[0-9]+)+" CMAKE_Fortran_COMPILER_VERSION ${PGI_VERSION})
    # No current requirements 
endif()

message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
set(CMAKE_Fortran_MODULE_DIRECTORY ${HIREMD_BINARY_DIR}/modules)
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB DUMMY_HIRE ${GIT_ROOT}/HiRE_lib/hire_interface.F90)

set(ALL_DUMMIES ${DUMMY_HIRE})

# Glob all the sources
file(GLOB MD_SOURCES *.f90
                     *.f
		     *.F90
                     utils/*.f90)

file(GLOB NOT_MD_SRC myblas.f
                     mylapack.f)

list(REMOVE_ITEM MD_SOURCES ${NOT_MD_SRC} 
                            ${ALL_DUMMIES} )

# Collect all the extra source files into one place. It's just dummies at the moment.
set(ALL_EXTRA_SOURCES ${ALL_DUMMIES})

########################
# LIBRARIES START HERE #
########################

# Find BLAS and LAPACK libraries
find_package(MYBLAS)
find_package(MYLAPACK)

# Make an extras library for dummies, this is necessary because some of the
# files in gminlib use modules made by extralib
add_library(extralib ${ALL_EXTRA_SOURCES})
set_module_dir(extralib)
set_target_properties(extralib PROPERTIES LINKER_LANGUAGE "Fortran") 
set_target_properties(extralib PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};__SPARSE")


add_library(mdlib ${MD_SOURCES})
set_module_dir(mdlib)
set_module_depends(mdlib extralib)
set_target_properties(mdlib PROPERTIES LINKER_LANGUAGE "Fortran") 


# MPI flag 
option(WITH_MPI "Enable MPI compilation" OFF)
if(WITH_MPI)
  # append -DMPI to compile definitions
  add_definitions(-DMPI) 
  message("building with MPI: \n  ${CMAKE_Fortran_FLAGS}")   
endif(WITH_MPI)

option(ANALYSIS_VER "Enable printing for structural analysis" OFF)
if(ANALYSIS_VER)
        add_definitions(-DFOR_ANALYSIS)
endif(ANALYSIS_VER)

option(PROF_VER "Enable profiling" OFF)
if(PROF_VER)
	add_definitions(-DFOR_PROFILING)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
        SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif(PROF_VER)

add_definitions(-D_GIT_ROOT_='${GIT_ROOT}')

##########################
# EXECUTABLES START HERE #
##########################

add_subdirectory(${GIT_ROOT}/HiRE_lib hire_interface)



add_executable(HIREMD MD_HiRE.F90 ${DUMMY_HIRE})
set_module_dir(HIREMD)
set_module_depends(HIREMD mdlib hire_base)
set_target_properties(HIREMD PROPERTIES LINKER_LANGUAGE "Fortran") 
set_target_properties(HIREMD PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};__HIRE")
target_link_libraries(HIREMD mdlib hire_base ${MYLAPACK_LIBS})
