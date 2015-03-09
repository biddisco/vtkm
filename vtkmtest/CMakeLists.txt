##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##
##  Copyright 2014 Sandia Corporation.
##  Copyright 2014 UT-Battelle, LLC.
##  Copyright 2014. Los Alamos National Security
##
##  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
##  the U.S. Government retains certain rights in this software.
##
##  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
##  Laboratory (LANL), the U.S. Government retains certain rights in
##  this software.
##============================================================================
cmake_minimum_required(VERSION 3.0.0)

#--------------------------------------------------
# Setup Output directories
#--------------------------------------------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Executables."
)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Libraries"
)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all static libraries."
)

#--------------------------------------------------
# load script for checking out projects from git
#--------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/CMake/GitExternal.cmake)

#--------------------------------------------------
# Windows link problems with GLUT
#--------------------------------------------------
#set(BUILD_SHARED_LIBS ON)

#--------------------------------------------------
# Find VTKm
#--------------------------------------------------
#find_package(VTKm)
if(WIN32)
  include_directories("D:/Code/vtkm")
  include_directories("D:/Build/vtkm")
elseif(APPLE)
  include_directories("/Users/biddisco/src/vtkm/")
  include_directories("/Users/biddisco/build/vtkm/")
endif()

#--------------------------------------------------
# Find OpenGL
#--------------------------------------------------
find_package(OpenGL)

#--------------------------------------------------
# Find Cuda
#--------------------------------------------------
find_package(Cuda)

#--------------------------------------------------
# Find Boost
#--------------------------------------------------
find_package(Boost 1.55 REQUIRED)

#--------------------------------------------------
# Find GLUT
#--------------------------------------------------
find_package(GLUT)
if(NOT GLUT_FOUND)
  option(GLUT_DOWNLOAD_AS_SUBPROJECT OFF) 
  if(GLUT_DOWNLOAD_AS_SUBPROJECT)

    # checkout master branch, but don't track it, verbose output
    Git_External(
        ${CMAKE_CURRENT_SOURCE_DIR}/freeglut
        https://github.com/dcnieho/FreeGLUT.git
        master
        NO_UPDATE
        VERBOSE
    )
    option(FREEGLUT_BUILD_DEMOS OFF "")
    option(FREEGLUT_BUILD_SHARED_LIBS OFF "")

    add_subdirectory(freeglut/freeglut/freeglut)
    link_directories(${CMAKE_CURRENT_BINARY_DIR}/freeglut/freeglut/freeglut/lib)
    add_definitions(-DFREEGLUT_STATIC)
    set(GLUT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/freeglut/freeglut/freeglut/include)
    set(GLUT_LIBRARIES freeglut_static)
  else()
    message("GLUT was not found, set GLUT_DOWNLOAD_AS_SUBPROJECT ON "
      "to download and build it as part of this project")
  endif()
ENDIF()

#--------------------------------------------------
# include dirs
#--------------------------------------------------
include_directories(
    ${Boost_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
)

add_executable(isosurface isosurface.cpp)
target_include_directories(isosurface PRIVATE ${GLUT_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR})
target_link_libraries(isosurface ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

cuda_add_executable(cuda_isosurface isosurface.cu)
target_link_libraries(cuda_isosurface ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
