cmake_minimum_required(VERSION 3.0)

project(hire_interface VERSION 0.0)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be run in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

# git root directory should be one directory above the source directory
get_filename_component(GIT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${GIT_ROOT}/CMakeModules)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  # Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)

# Add all the source files
file(GLOB HiRE_SOURCES ${GIT_ROOT}/HiRE_lib/*.F90 ${GIT_ROOT}/HiRE_lib/*.f90)
# Remove HiRE_interface.f90
file(GLOB HIREINTERFACE ${GIT_ROOT}/HiRE_lib/HiRE_interface.F90)

list(REMOVE_ITEM HiRE_SOURCES ${HIREINTERFACE})

# Build a base library
add_library(hire_base ${HiRE_SOURCES})
set_module_dir(hire_base)
    
# Build the interface library
add_library(hirelib ${HIREINTERFACE})
set_module_dir(hirelib)
set_module_depends(hirelib hire_base)
set_target_properties(hirelib PROPERTIES LINKER_LANGUAGE "Fortran")
set_target_properties(hirelib PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};__HIRE") 
target_link_libraries(hirelib hire_base)