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

#ifndef vtk_m_cont_internal_ArrayManagerExecutionHPX_h
#define vtk_m_cont_internal_ArrayManagerExecutionHPX_h

#include <vtkm/cont/internal/ArrayManagerExecution.h>
#include <vtkm/cont/internal/ArrayManagerExecutionShareWithControl.h>
#include <vtkm/cont/hpx/internal/DeviceAdapterTagHPX.h>

namespace vtkm {
namespace cont {
namespace internal {

template <typename T, class StorageTag>
class ArrayManagerExecution<T, StorageTag, vtkm::cont::DeviceAdapterTagHPX>
    : public vtkm::cont::internal::ArrayManagerExecutionShareWithControl
          <T, StorageTag>
{
public:
  typedef vtkm::cont::internal::ArrayManagerExecutionShareWithControl
      <T, StorageTag> Superclass;
  typedef typename Superclass::ValueType ValueType;
  typedef typename Superclass::PortalType PortalType;
  typedef typename Superclass::PortalConstType PortalConstType;
  typedef typename Superclass::StorageType StorageType;

  VTKM_CONT_EXPORT
  ArrayManagerExecution(StorageType *storage)
    : Superclass(storage) {  }

  VTKM_CONT_EXPORT
  PortalConstType PrepareForInput(bool updateData)
  {
    return this->Superclass::PrepareForInput(updateData);
  }

  VTKM_CONT_EXPORT
  PortalType PrepareForInPlace(bool updateData)
  {
    return this->Superclass::PrepareForInPlace(updateData);
  }

  VTKM_CONT_EXPORT
  PortalType PrepareForOutput(vtkm::Id numberOfValues)
  {
    return this->Superclass::PrepareForOutput(numberOfValues);
  }

};

}
}
} // namespace vtkm::cont::internal

#endif //vtk_m_cont_internal_ArrayManagerExecutionHPX_h
