##=============================================================================
##
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##
##  Copyright 2016 Sandia Corporation.
##  Copyright 2016 UT-Battelle, LLC.
##  Copyright 2016 Los Alamos National Security.
##
##  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
##  the U.S. Government retains certain rights in this software.
##  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
##  Laboratory (LANL), the U.S. Government retains certain rights in
##  this software.
##
##=============================================================================

set(headers
  DeviceAdapterHPX.h
  internal/ArrayManagerExecutionHPX.h  
  internal/DeviceAdapterAlgorithmHPX.h 
  internal/DeviceAdapterTagHPX.h
)

add_definitions(-DVTKM_DEVICE_CONFIG_INCLUDE=\"vtkm/cont/hpx/vtkm_hpx.hpp\")

#-----------------------------------------------------------------------------
vtkm_declare_headers(${impl_headers} ${headers})

#-----------------------------------------------------------------------------
add_subdirectory(testing)
