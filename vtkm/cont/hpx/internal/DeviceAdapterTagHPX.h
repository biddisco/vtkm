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

#ifndef vtk_m_cont_internal_DeviceAdapterTagHPX_h
#define vtk_m_cont_internal_DeviceAdapterTagHPX_h

#include <vtkm/cont/internal/DeviceAdapterTag.h>

//We always create the tbb tag when included, but we only mark it as
//a valid tag when VTKM_ENABLE_HPX is true. This is for easier development
//of multi-backend systems
#ifdef VTKM_ENABLE_HPX
VTKM_VALID_DEVICE_ADAPTER(HPX, VTKM_DEVICE_ADAPTER_HPX);
#else
VTKM_INVALID_DEVICE_ADAPTER(HPX, VTKM_DEVICE_ADAPTER_HPX);
#endif

#endif //vtk_m_cont_internal_DeviceAdapterTagHPX_h
