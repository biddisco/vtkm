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

#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/cont/DispatcherMapField.h>

#include <vtkm/exec/WorkletMapField.h>

#include <vtkm/cont/cuda/internal/testing/Testing.h>

namespace ut_transform {

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
struct MySquare
{
  VTKM_EXEC_CONT_EXPORT
  MySquare() {}

  template<typename U>
  VTKM_EXEC_CONT_EXPORT
  ValueType operator()(U u) const
    { return ValueType(u*u); }
};


template< typename ValueType, typename FunctorType >
struct TransformTests
{
  typedef  vtkm::cont::ArrayHandleTransform< ValueType,
                                            vtkm::cont::ArrayHandle< vtkm::Id >,
                                            FunctorType > TransformHandle;

  typedef  vtkm::cont::ArrayHandleTransform< ValueType,
                                    vtkm::cont::ArrayHandleCounting< vtkm::Id >,
                                    FunctorType > CountingTransformHandle;

  void operator()(const ValueType, FunctorType functor) const
  {

    //test the make_ArrayHandleTransform method
    //test a transform handle with a counting handle as the values
    vtkm::cont::ArrayHandleCounting< vtkm::Id > counting =
            vtkm::cont::make_ArrayHandleCounting(vtkm::Id(0),ARRAY_SIZE);
    CountingTransformHandle countingTransformed =
            vtkm::cont::make_ArrayHandleTransform<ValueType>(counting,functor);

    {
    vtkm::cont::ArrayHandle< ValueType > result;
    vtkm::cont::DispatcherMapField< ut_transform::Pass >().Invoke(
                                                     countingTransformed,
                                                     result);

    DAX_TEST_ASSERT(ARRAY_SIZE == result.GetNumberOfValues(),
                    "result handle doesn't have the correct size");
    //verify that the control portal works
    for(int i=0; i < ARRAY_SIZE; ++i)
      {
      const ValueType v = result.GetPortalConstControl().Get(i);
      const ValueType correct_value = MySquare<ValueType>()(i);
      DAX_TEST_ASSERT(v == correct_value,
                      "Transform Handle with MySquare Failed");
      }
    }

    {
    //test a transform handle with a normal handle as the values
    //we are going to connect the two handles up, and than fill
    //the values and make the transform sees the new values in the handle
    vtkm::cont::ArrayHandle< vtkm::Id > input;
    TransformHandle thandle(input,functor);

    vtkm::cont::DispatcherMapField< ut_transform::Pass >().Invoke(
      vtkm::cont::make_ArrayHandleCounting(2,ARRAY_SIZE),
      input);

    vtkm::cont::ArrayHandle< ValueType > result;
    vtkm::cont::DispatcherMapField< ut_transform::Pass >().Invoke(
                                                     thandle,
                                                     result);

    //verify that the control portal works
    for(int i=0; i < ARRAY_SIZE; ++i)
      {
      const ValueType v = result.GetPortalConstControl().Get(i);
      const ValueType correct_value = MySquare<ValueType>()(i+2);
      DAX_TEST_ASSERT(v == correct_value,
                      "Transform Handle with MySquare Failed");
      }

    //now update the array handle values again, so that
    //we can verify the transform handle gets these updated values
    vtkm::cont::DispatcherMapField< ut_transform::Pass >().Invoke(
      vtkm::cont::make_ArrayHandleCounting(500,ARRAY_SIZE),
      input);

    //verify that the transform has the new values
    vtkm::cont::DispatcherMapField< ut_transform::Pass >().Invoke(
                                                     thandle,
                                                     result);
    for(int i=0; i < ARRAY_SIZE; ++i)
      {
      const ValueType v = result.GetPortalConstControl().Get(i);
      const ValueType correct_value = MySquare<ValueType>()(500+i);
      DAX_TEST_ASSERT(v == correct_value,
                      "Transform Handle with MySquare Failed");
      }
    }

  }

};


template <typename T, typename F>
void RunTransformTests(const T t, F f)
{
  TransformTests<T,F> tests;
  tests(t,f);
}

void TestArrayHandleTransform()
{
  RunTransformTests( vtkm::Id(0), MySquare<vtkm::Id>() );
  RunTransformTests( vtkm::Scalar(0), MySquare<vtkm::Scalar>() );
  // RunTransformTests( vtkm::Vector3(vtkm::Scalar(0)), MySquare<vtkm::Vector3>() );
}



} // ut_transform namespace

int UnitTestCudaArrayHandleTransform(int, char *[])
{
  return vtkm::cont::cuda::internal::Testing::Run(
                                     ut_transform::TestArrayHandleTransform);
}
