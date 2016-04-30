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

set(benchmarks
	BenchmarkDeviceAdapter.cxx
  )

vtkm_clear_benchmarks()
vtkm_save_benchmarks(SERIAL ${benchmarks})
vtkm_benchmarks(VTKM_DEVICE_ADAPTER_SERIAL)

if (VTKm_ENABLE_CUDA)
  vtkm_clear_benchmarks()
  vtkm_save_benchmarks(CUDA ${benchmarks})
  vtkm_benchmarks(VTKM_DEVICE_ADAPTER_CUDA)
endif()
if (VTKm_ENABLE_TBB)
  vtkm_clear_benchmarks()
  vtkm_save_benchmarks(TBB ${benchmarks})
  vtkm_benchmarks(VTKM_DEVICE_ADAPTER_TBB)
endif()
if (VTKm_ENABLE_HPX)
  vtkm_clear_benchmarks()
  vtkm_save_benchmarks(HPX ${benchmarks})
  vtkm_benchmarks(VTKM_DEVICE_ADAPTER_HPX)
endif()

