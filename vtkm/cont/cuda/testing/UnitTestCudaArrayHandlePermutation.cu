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
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/ArrayHandlePermutation.h>
#include <vtkm/cont/DispatcherMapField.h>
#include <vtkm/VectorTraits.h>

#include <vtkm/exec/WorkletMapField.h>

#include <vtkm/cont/cuda/internal/testing/Testing.h>

namespace ut_Permutation
{

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

template<typename T>
struct CountByThree
{
  VTKM_EXEC_CONT_EXPORT
  CountByThree(): Value() {}

  VTKM_EXEC_CONT_EXPORT
  explicit CountByThree(T t): Value(t) {}

  VTKM_EXEC_CONT_EXPORT
  CountByThree<T> operator+(vtkm::Id count) const
  { return CountByThree<T>(Value+(count*3)); }

  VTKM_CONT_EXPORT
  bool operator==(const CountByThree<T>& other) const
  { return Value == other.Value; }

  T Value;
};

}


namespace vtkm {

/// Implement Vector Traits for CountByThree so we can use it in the execution
/// environment
template<typename T>
struct VectorTraits< ut_Permutation::CountByThree<T> >
{
  typedef ut_Permutation::CountByThree<T> ComponentType;
  static const int NUM_COMPONENTS = 1;
  typedef VectorTraitsTagSingleComponent HasMultipleComponents;

  VTKM_EXEC_EXPORT
  static const ComponentType &GetComponent(const ComponentType &vector,
                                           int component) {
    return vector;
  }
  VTKM_EXEC_EXPORT
  static ComponentType &GetComponent(ComponentType &vector, int component) {
    return vector;
  }

  VTKM_EXEC_EXPORT static void SetComponent(ComponentType &vector,
                                           int component,
                                           ComponentType value) {
    vector = value;
  }

  VTKM_EXEC_CONT_EXPORT
  static vtkm::Tuple<ComponentType,NUM_COMPONENTS>
  ToTuple(const ComponentType &vector)
  {
    return vtkm::Tuple<T,1>(vector);
  }
};

}


namespace ut_Permutation {

template< typename ValueType>
struct TemplatedTests
{
  typedef vtkm::cont::ArrayHandleCounting<ValueType> CountingArrayHandleType;

  typedef vtkm::cont::ArrayHandlePermutation<
            vtkm::cont::ArrayHandle<vtkm::Id>, //key type
            CountingArrayHandleType > ArrayPermHandleType;

  typedef vtkm::cont::ArrayHandlePermutation<
            vtkm::cont::ArrayHandleCounting<vtkm::Id>, //key type
            CountingArrayHandleType > ArrayCountPermHandleType;

  void operator()( const ValueType startingValue )
  {
    vtkm::Id everyOtherBuffer[ARRAY_SIZE/2];
    vtkm::Id fullBuffer[ARRAY_SIZE];

    for (vtkm::Id index = 0; index < ARRAY_SIZE; index++)
      {
      everyOtherBuffer[index/2] = index; //1,3,5,7,9
      fullBuffer[index] = index;
      }

  {
  //verify the different constructors work
  CountingArrayHandleType counting(startingValue,ARRAY_SIZE);
  vtkm::cont::ArrayHandleCounting<vtkm::Id> keys(vtkm::Id(0),ARRAY_SIZE);

  ArrayCountPermHandleType permutation_constructor(keys,counting);

  ArrayCountPermHandleType make_permutation_constructor =
      vtkm::cont::make_ArrayHandlePermutation(
              vtkm::cont::make_ArrayHandleCounting(vtkm::Id(0),ARRAY_SIZE),
              vtkm::cont::make_ArrayHandleCounting(startingValue, ARRAY_SIZE));

  }

  //make a short permutation array, verify its length and values
  {
  vtkm::cont::ArrayHandle<vtkm::Id> keys =
      vtkm::cont::make_ArrayHandle(everyOtherBuffer, ARRAY_SIZE/2);
  CountingArrayHandleType values(startingValue,ARRAY_SIZE);

  ArrayPermHandleType permutation =
      vtkm::cont::make_ArrayHandlePermutation(keys,values);

  //verify the results by executing a worklet to copy the results
  vtkm::cont::ArrayHandle< ValueType > result;
  vtkm::cont::DispatcherMapField< ut_Permutation::Pass >().Invoke(
                                                     permutation, result);

  typename vtkm::cont::ArrayHandle< ValueType >::PortalConstControl portal =
                                      result.GetPortalConstControl();
  ValueType correct_value = startingValue;
  for(int i=0; i < ARRAY_SIZE/2; ++i)
    {
     //the permutation should have every other value
    correct_value = correct_value + 1;
    ValueType v = portal.Get(i);
    DAX_TEST_ASSERT(v == correct_value, "Count By Three permutation wrong");
    correct_value = correct_value + 1;
    }
  }

  //make a long permutation array, verify its length and values
  {
  vtkm::cont::ArrayHandle<vtkm::Id> keys =
      vtkm::cont::make_ArrayHandle(fullBuffer, ARRAY_SIZE);
  CountingArrayHandleType values(startingValue,ARRAY_SIZE);

  ArrayPermHandleType permutation =
      vtkm::cont::make_ArrayHandlePermutation(keys,values);

  //verify the results by executing a worklet to copy the results
  vtkm::cont::ArrayHandle< ValueType > result;
  vtkm::cont::DispatcherMapField< ut_Permutation::Pass >().Invoke(
                                                     permutation, result);

  typename vtkm::cont::ArrayHandle< ValueType >::PortalConstControl portal =
                                      result.GetPortalConstControl();
  ValueType correct_value = startingValue;
  for(int i=0; i < ARRAY_SIZE; ++i)
    {
    ValueType v = portal.Get(i);
    DAX_TEST_ASSERT(v == correct_value, "Full permutation wrong");
    correct_value = correct_value + 1;
    }
  }


  }
};

struct TestFunctor
{
  template <typename T>
  void operator()(const T t)
  {
    TemplatedTests<T> tests;
    tests(t);
  }
};

void TestArrayHandlePermutation()
{
  TestFunctor()( vtkm::Id(0) );
  TestFunctor()( vtkm::Scalar(0) );
  TestFunctor()( CountByThree<vtkm::Id>(12) );
  TestFunctor()( CountByThree<vtkm::Scalar>(1.2f) );
}



} // annonymous namespace

int UnitTestCudaArrayHandlePermutation(int, char *[])
{
  return vtkm::cont::cuda::internal::Testing::Run(
                                  ut_Permutation::TestArrayHandlePermutation);
}
