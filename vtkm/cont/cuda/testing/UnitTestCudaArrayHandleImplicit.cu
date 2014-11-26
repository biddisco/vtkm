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
//  Copyright 2014. Los Alamos National Security
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#define BOOST_SP_DISABLE_THREADS

#include <vtkm/cont/cuda/DeviceAdapterCuda.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleImplicit.h>
#include <vtkm/cont/DispatcherMapField.h>

#include <vtkm/exec/WorkletMapField.h>

#include <vtkm/cont/cuda/internal/testing/Testing.h>

namespace ut_implicit{

const vtkm::Id ARRAY_SIZE = 10;

struct Pass : public vtkm::exec::WorkletMapField
{
  typedef void ControlSignature(FieldIn, FieldOut);
  typedef _2 ExecutionSignature(_1);

  template<class ValueType>
  VTKM_EXEC_EXPORT
  ValueType operator()(const ValueType &inValue) const
  { return inValue; }

};

template<typename ValueType>
struct IndexSquared
{
  VTKM_EXEC_CONT_EXPORT
  ValueType operator()(vtkm::Id i) const
    { return ValueType(i*i); }
};


template< typename ValueType, typename FunctorType >
struct ImplicitTests
{
  typedef vtkm::cont::ArrayHandleImplicit<ValueType,FunctorType> ImplicitHandle;

  void operator()(const ValueType, FunctorType functor) const
  {
    ImplicitHandle implict =
            vtkm::cont::make_ArrayHandleImplicit<ValueType>(functor,ARRAY_SIZE);

    vtkm::cont::ArrayHandle< ValueType > result;
    vtkm::cont::DispatcherMapField< ut_implicit::Pass >().Invoke(
                                                     implict, result);

    //verify that the control portal works
    for(int i=0; i < ARRAY_SIZE; ++i)
      {
      const ValueType v = result.GetPortalConstControl().Get(i);
      const ValueType correct_value = IndexSquared<ValueType>()(i);
        DAX_TEST_ASSERT(v == correct_value,
                      "Implicit Handle with IndexSquared Failed");
      }
  }

};


template <typename T, typename F>
void RunImplicitTests(const T t, F f)
{
  ImplicitTests<T,F> tests;
  tests(t,f);
}

void TestArrayHandleImplicit()
{
  RunImplicitTests( vtkm::Id(0), IndexSquared<vtkm::Id>() );
  RunImplicitTests( vtkm::Scalar(0), IndexSquared<vtkm::Scalar>() );
  RunImplicitTests( vtkm::Vector3(vtkm::Scalar(0)), IndexSquared<vtkm::Vector3>() );
}



} // ut_implicit namespace

int UnitTestCudaArrayHandleImplicit(int, char *[])
{
  return vtkm::cont::cuda::internal::Testing::Run(
                                      ut_implicit::TestArrayHandleImplicit);
}
