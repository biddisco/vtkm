//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 Sandia Corporation.
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef vtkm_cont_internal_StorageError_h
#define vtkm_cont_internal_StorageError_h


namespace vtkm {
namespace cont {
namespace internal {

/// This is an invalid Storage. The point of this class is to include the
/// header file to make this invalid class the default Storage. From that
/// point, you have to specify an appropriate Storage or else get a compile
/// error.
///
struct StorageTagError
{
  // Not implemented.
};

}
}
} // namespace vtkm::cont::internal

#endif //vtkm_cont_internal_StorageError_h
