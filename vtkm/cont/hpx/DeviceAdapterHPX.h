//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2016 Sandia Corporation.
//  Copyright 2016 UT-Battelle, LLC.
//  Copyright 2016 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//
//=============================================================================

#ifndef vtk_m_cont_DeviceAdapterHPX_h
#define vtk_m_cont_DeviceAdapterHPX_h

#include <vtkm/cont/hpx/internal/DeviceAdapterTagHPX.h>

#ifdef VTKM_ENABLE_HPX
# include <vtkm/cont/hpx/internal/ArrayManagerExecutionHPX.h>
# include <vtkm/cont/hpx/internal/DeviceAdapterAlgorithmHPX.h>
#endif

#endif //vtk_m_cont_DeviceAdapterHPX_h