//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2015 Sandia Corporation.
//  Copyright 2015 UT-Battelle, LLC.
//  Copyright 2015 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#include <vtkm/cont/arg/TransportTagArrayInOut.h>

#include <vtkm/exec/FunctorBase.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DeviceAdapter.h>

#include <vtkm/cont/testing/Testing.h>

namespace {

static const vtkm::Id ARRAY_SIZE = 10;

template<typename PortalType>
struct TestKernel : public vtkm::exec::FunctorBase
{
  PortalType Portal;

  VTKM_EXEC_EXPORT
  void operator()(vtkm::Id index) const
  {
    typedef typename PortalType::ValueType ValueType;
    ValueType inValue = this->Portal.Get(index);
    this->Portal.Set(index, inValue + inValue);
  }
};

template<typename Device>
struct TryArrayInOutType
{
  template<typename T>
  void operator()(T) const
  {
    T array[ARRAY_SIZE];
    for (vtkm::Id index = 0; index < ARRAY_SIZE; index++)
    {
      array[index] = TestValue(index, T());
    }

    typedef vtkm::cont::ArrayHandle<T> ArrayHandleType;
    ArrayHandleType handle = vtkm::cont::make_ArrayHandle(array, ARRAY_SIZE);

    typedef typename ArrayHandleType::
        template ExecutionTypes<Device>::Portal PortalType;

    vtkm::cont::arg::Transport<
        vtkm::cont::arg::TransportTagArrayInOut, ArrayHandleType, Device>
        transport;

    TestKernel<PortalType> kernel;
    kernel.Portal = transport(handle, ARRAY_SIZE);

    vtkm::cont::DeviceAdapterAlgorithm<Device>::Schedule(kernel, ARRAY_SIZE);

    typename ArrayHandleType::PortalConstControl portal =
        handle.GetPortalConstControl();
    VTKM_TEST_ASSERT(portal.GetNumberOfValues() == ARRAY_SIZE,
                     "Portal has wrong number of values.");
    for (vtkm::Id index = 0; index < ARRAY_SIZE; index++)
    {
      T expectedValue = TestValue(index, T()) + TestValue(index, T());
      T retrievedValue = portal.Get(index);
      VTKM_TEST_ASSERT(test_equal(expectedValue, retrievedValue),
                       "Functor did not modify in place.");
    }
  }
};

template<typename Device>
void TryArrayInOutTransport(Device)
{
  vtkm::testing::Testing::TryTypes(TryArrayInOutType<Device>(),
                                   vtkm::TypeListTagCommon());
}

void TestArrayInOutTransport()
{
  TryArrayInOutTransport(VTKM_DEFAULT_DEVICE_ADAPTER_TAG());
}

} // anonymous namespace

int UnitTestTransportArrayInOut(int, char *[])
{
  return vtkm::cont::testing::Testing::Run(TestArrayInOutTransport);
}
