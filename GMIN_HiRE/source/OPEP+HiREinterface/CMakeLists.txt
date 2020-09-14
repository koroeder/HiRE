###########################################################
#
# CMAKE FILE TO GENERATE libOPEP.a 
# 
############################################################

cmake_minimum_required(VERSION 3.0)

# Silence warnings about PROJECT_VERSION not being set
cmake_policy(SET CMP0048 OLD)

include_directories(${CMAKE_BINARY_DIR})
file(GLOB OPEP_SOURCES *.f90 *.f *.F *.F90)

# Remove the interface file
file(GLOB NOT_OPEP_SOURCES opep_interface.F90)
list(REMOVE_ITEM OPEP_SOURCES ${NOT_OPEP_SOURCES} )

# Create the library
add_library(OPEP+HiRE ${OPEP_SOURCES})
add_dependencies(OPEP+HiRE gminlib)
