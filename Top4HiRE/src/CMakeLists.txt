cmake_minimum_required(VERSION 3.0)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(TOP4HIRE)

set(PROJECT_VERSION "1.0")
set(PROJECT_CONTACT "konstantin.roeder@kcl.ac.uk")

get_filename_component(GIT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

enable_language(Fortran)
set (CMAKE_Fortran_FLAGS "-ffixed-line-length-132 -ffree-line-length-0" CACHE TYPE STRING FORCE)
set (CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE TYPE STRING FORCE)
set (CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none" CACHE TYPE STRING FORCE)
set (FORTRAN_FREEFORM_FLAG "-ffree-form" CACHE TYPE STRING)

include_directories(${CMAKE_SOURCE_DIR})
# Glob all the sources
file(GLOB TOP_LIB_SOURCES ${GIT_ROOT}/src/*.f90)
file(GLOB NOT_SOURCES ${GIT_ROOT}/src/main.F90)
list(REMOVE_ITEM TOP_LIB_SOURCES ${NOT_SOURCES})

add_library(toplib ${TOP_LIB_SOURCES})

add_executable(Top4HiRE ${GIT_ROOT}/src/main.F90)
set_target_properties(Top4HiRE PROPERTIES LINKER_LANGUAGE "Fortran") 
target_link_libraries(Top4HiRE toplib)