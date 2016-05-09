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
##  Copyright 2014 Los Alamos National Security.
##
##  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
##  the U.S. Government retains certain rights in this software.
##
##  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
##  Laboratory (LANL), the U.S. Government retains certain rights in
##  this software.
##============================================================================

if (VTKm_HPX_initialize_complete)
  return()
endif (VTKm_HPX_initialize_complete)

vtkm_configure_device(Base)

if (VTKm_Base_FOUND)

  set(VTKm_HPX_FOUND ${VTKm_ENABLE_HPX})
  if (NOT VTKm_HPX_FOUND)
    message(STATUS "This build of VTK-m does not include HPX.")
  endif ()

  #---------------------------------------------------------------------------
  # Find HPX.
  #---------------------------------------------------------------------------
  if (VTKm_HPX_FOUND)
    find_package(HPX)
    if (NOT HPX_FOUND)
      message(STATUS "HPX not found")
      set(VTKm_HPX_FOUND)
    endif ()
  endif()

endif ()

#-----------------------------------------------------------------------------
# Set up the compiler flag optimizations
#-----------------------------------------------------------------------------
include(VTKmCompilerOptimizations)

#-----------------------------------------------------------------------------
# Set up all these dependent packages (if they were all found).
#-----------------------------------------------------------------------------
if (VTKm_HPX_FOUND)
  set(VTKm_INCLUDE_DIRS
    ${VTKm_INCLUDE_DIRS}
    ${HPX_INCLUDE_DIRS}
    )
  set(VTKm_LIBRARIES ${HPX_LIBRARIES})

  set(VTKm_HPX_initialize_complete TRUE)
endif()
