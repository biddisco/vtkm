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
#ifndef vtk_m_cont_testing_TestingArrayHandles_h
#define vtk_m_cont_testing_TestingArrayHandles_h

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/TypeTraits.h>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/cont/DeviceAdapterSerial.h>
#include <vtkm/cont/testing/Testing.h>

#include <boost/type_traits/is_same.hpp>
#include <algorithm>
#include <vector>

namespace vtkm {
namespace cont {
namespace testing {

namespace array_handle_testing
{
  template<class IteratorType, typename T>
  void CheckValues(IteratorType begin, IteratorType end, T)
  {

    vtkm::Id index = 0;
    for (IteratorType iter = begin; iter != end; iter++)
    {
      T expectedValue = TestValue(index, T());
      if (!test_equal(*iter, expectedValue))
      {
        std::stringstream message;
        message << "Got unexpected value in array." << std::endl
                << "Expected: " << expectedValue
                << ", Found: " << *iter << std::endl;
        VTKM_TEST_FAIL(message.str().c_str());
      }

      index++;
    }
  }

  template<typename T>
  void CheckArray(const vtkm::cont::ArrayHandle<T> &handle)
  {
    CheckPortal(handle.GetPortalConstControl());
  }

}

/// This class has a single static member, Run, that tests that all Fancy Array
/// Handles work with the given DeviceAdapter
///
template<class DeviceAdapterTag>
struct TestingArrayHandles
{

  struct PassThrough : public vtkm::worklet::WorkletMapField
  {
    typedef void ControlSignature(FieldIn<>, FieldOut<>);
    typedef _2 ExecutionSignature(_1);

    template<class ValueType>
    VTKM_EXEC_EXPORT
    ValueType operator()(const ValueType &inValue) const
    { return inValue; }

  };

  template<typename T, typename ExecutionPortalType>
  struct AssignTestValue : public vtkm::exec::FunctorBase
  {
    ExecutionPortalType Portal;
    VTKM_CONT_EXPORT
    AssignTestValue(ExecutionPortalType p): Portal(p) {}

    VTKM_EXEC_EXPORT
    void operator()(vtkm::Id index) const
    {
      this->Portal.Set(index, TestValue(index, T()) );
    }
  };

  template<typename T, typename ExecutionPortalType>
  struct InplaceFunctor : public vtkm::exec::FunctorBase
  {
    ExecutionPortalType Portal;
    VTKM_CONT_EXPORT
    InplaceFunctor(ExecutionPortalType p): Portal(p) {}

    VTKM_EXEC_EXPORT
    void operator()(vtkm::Id index) const
    {
      this->Portal.Set(index, T(this->Portal.Get(index) + T(1)));
    }
  };

private:
  static const vtkm::Id ARRAY_SIZE = 100;

  typedef vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapterTag> Algorithm;

  typedef vtkm::worklet::DispatcherMapField<PassThrough,
                                            DeviceAdapterTag> DispatcherPassThrough;
  struct VerifyEmptyArrays
  {
    template<typename T>
    VTKM_CONT_EXPORT void operator()(T) const
    {
      std::cout << "Try operations on empty arrays." << std::endl;
      // After each operation, reinitialize array in case something gets
      // allocated.
      vtkm::cont::ArrayHandle<T> arrayHandle = vtkm::cont::ArrayHandle<T>();
      VTKM_TEST_ASSERT(arrayHandle.GetNumberOfValues() == 0,
                       "Uninitialized array does not report zero values.");
      arrayHandle = vtkm::cont::ArrayHandle<T>();
      VTKM_TEST_ASSERT(
            arrayHandle.GetPortalConstControl().GetNumberOfValues() == 0,
            "Uninitialized array does not give portal with zero values.");
      arrayHandle = vtkm::cont::ArrayHandle<T>();
      arrayHandle.Shrink(0);
      arrayHandle = vtkm::cont::ArrayHandle<T>();
      arrayHandle.ReleaseResourcesExecution();
      arrayHandle = vtkm::cont::ArrayHandle<T>();
      arrayHandle.ReleaseResources();
      arrayHandle = vtkm::cont::ArrayHandle<T>();
      arrayHandle.PrepareForOutput(ARRAY_SIZE, DeviceAdapterTag());
    }
  };

  struct VerifyUserAllocatedHandle
  {
    template<typename T>
    VTKM_CONT_EXPORT void operator()(T) const
    {
      std::vector<T> buffer(ARRAY_SIZE);
	  for (vtkm::Id index = 0; index < ARRAY_SIZE; index++)
      {
        buffer[static_cast<std::size_t>(index)] = TestValue(index, T());
      }

      vtkm::cont::ArrayHandle<T> arrayHandle =
          vtkm::cont::make_ArrayHandle(buffer);

      VTKM_TEST_ASSERT(arrayHandle.GetNumberOfValues() == ARRAY_SIZE,
                       "ArrayHandle has wrong number of entries.");

      std::cout << "Check array with user provided memory." << std::endl;
      array_handle_testing::CheckArray(arrayHandle);

      std::cout << "Check out execution array behavior." << std::endl;
      { //as input
        typename vtkm::cont::ArrayHandle<T>::template
            ExecutionTypes<DeviceAdapterTag>::PortalConst
            executionPortal;
        executionPortal =
            arrayHandle.PrepareForInput(DeviceAdapterTag());

        //use a worklet to verify the input transfer worked properly
        vtkm::cont::ArrayHandle<T> result;
        DispatcherPassThrough().Invoke(arrayHandle, result);
        array_handle_testing::CheckArray(result);
      }

      std::cout << "Check out inplace." << std::endl;
      { //as inplace
        typename vtkm::cont::ArrayHandle<T>::template
            ExecutionTypes<DeviceAdapterTag>::Portal
            executionPortal;
        executionPortal =
            arrayHandle.PrepareForInPlace(DeviceAdapterTag());

        //use a worklet to verify the inplace transfer worked properly
        vtkm::cont::ArrayHandle<T> result;
        DispatcherPassThrough().Invoke(arrayHandle, result);
        array_handle_testing::CheckArray(result);
      }

      std::cout << "Check out output." << std::endl;
      { //as output with same length as user provided. This should work
        //as no new memory needs to be allocated
        typename vtkm::cont::ArrayHandle<T>::template
            ExecutionTypes<DeviceAdapterTag>::Portal
            executionPortal;
        executionPortal =
            arrayHandle.PrepareForOutput(ARRAY_SIZE,
                                        DeviceAdapterTag());

        //we can't verify output contents as those aren't fetched, we
        //can just make sure the allocation didn't throw an exception
      }

      std::cout << "Check CopyInto from control array" << std::endl;
      { //Release the execution resources so that data is only
        //in the control environment
        arrayHandle.ReleaseResourcesExecution();

        //Copy data from handle into iterator
        T array[ARRAY_SIZE];
        arrayHandle.CopyInto(array, DeviceAdapterTag());
        array_handle_testing::CheckValues(array, array+ARRAY_SIZE, T());
      }

      std::cout << "Check CopyInto from execution array" << std::endl;
      { //Copy the data to the execution environment
        vtkm::cont::ArrayHandle<T> result;
        DispatcherPassThrough().Invoke(arrayHandle, result);

        //Copy data from handle into iterator
        T array[ARRAY_SIZE];
        result.CopyInto(array, DeviceAdapterTag());
        array_handle_testing::CheckValues(array, array+ARRAY_SIZE, T());
      }

      if (!boost::is_same<DeviceAdapterTag, vtkm::cont::DeviceAdapterTagSerial>::value)
      {
        std::cout << "Check using different device adapter" << std::endl;
        //Copy the data to the execution environment
        vtkm::cont::ArrayHandle<T> result;
        DispatcherPassThrough().Invoke(arrayHandle, result);

        //CopyInto allows you to copy the data even
        //if you request it from a different device adapter
        T array[ARRAY_SIZE];
        result.CopyInto(array, vtkm::cont::DeviceAdapterTagSerial());
        array_handle_testing::CheckValues(array, array+ARRAY_SIZE, T());
      }

      { //as output with a length larger than the memory provided by the user
        //this should fail
        bool gotException = false;
        try
        {
          //you should not be able to allocate a size larger than the
          //user provided and get the results
          arrayHandle.PrepareForOutput(ARRAY_SIZE*2,DeviceAdapterTag());
          arrayHandle.GetPortalControl();
        }
        catch (vtkm::cont::Error &)
        {
          gotException = true;
        }
        VTKM_TEST_ASSERT(gotException,
                         "PrepareForOutput should fail when asked to "\
                         "re-allocate user provided memory.");
      }
    }
  };

  struct VerifyVTKMAllocatedHandle
  {
    template<typename T>
    VTKM_CONT_EXPORT void operator()(T) const
    {
      vtkm::cont::ArrayHandle<T> arrayHandle;

      VTKM_TEST_ASSERT(arrayHandle.GetNumberOfValues() == 0,
                       "ArrayHandle has wrong number of entries.");
      {
        typedef typename vtkm::cont::ArrayHandle<T>::template
            ExecutionTypes<DeviceAdapterTag>::Portal
              ExecutionPortalType;
          ExecutionPortalType executionPortal =
              arrayHandle.PrepareForOutput(ARRAY_SIZE*2,
                                           DeviceAdapterTag());

        //we drop down to manually scheduling so that we don't need
        //need to bring in array handle counting
        AssignTestValue<T, ExecutionPortalType> functor(executionPortal);
        Algorithm::Schedule(functor, ARRAY_SIZE*2);
      }

      VTKM_TEST_ASSERT(arrayHandle.GetNumberOfValues() == ARRAY_SIZE*2,
                       "Array not allocated correctly.");
      array_handle_testing::CheckArray(arrayHandle);

      std::cout << "Try shrinking the array." << std::endl;
      arrayHandle.Shrink(ARRAY_SIZE);
      VTKM_TEST_ASSERT(arrayHandle.GetNumberOfValues() == ARRAY_SIZE,
                       "Array size did not shrink correctly.");
      array_handle_testing::CheckArray(arrayHandle);

      std::cout << "Try reallocating array." << std::endl;
      arrayHandle.Allocate(ARRAY_SIZE*2);
      VTKM_TEST_ASSERT(arrayHandle.GetNumberOfValues() == ARRAY_SIZE*2,
                       "Array size did not allocate correctly.");
      // No point in checking values. This method can invalidate them.

      std::cout << "Try in place operation." << std::endl;
      {
        typedef typename vtkm::cont::ArrayHandle<T>::template
          ExecutionTypes<DeviceAdapterTag>::Portal
            ExecutionPortalType;
        ExecutionPortalType executionPortal =
            arrayHandle.PrepareForInPlace(DeviceAdapterTag());

        //in place can't be done through the dispatcher
        //instead we have to drop down to manually scheduling
        InplaceFunctor<T, ExecutionPortalType> functor(executionPortal);
        Algorithm::Schedule(functor,ARRAY_SIZE*2);
      }
      typename vtkm::cont::ArrayHandle<T>::PortalConstControl controlPortal =
          arrayHandle.GetPortalConstControl();
      for (vtkm::Id index = 0; index < ARRAY_SIZE; index++)
      {
        VTKM_TEST_ASSERT(test_equal(controlPortal.Get(index),
                                    TestValue(index, T()) + T(1)),
                         "Did not get result from in place operation.");
      }

      VTKM_TEST_ASSERT(arrayHandle == arrayHandle,
                       "Array handle does not equal itself.");
      VTKM_TEST_ASSERT(arrayHandle != vtkm::cont::ArrayHandle<T>(),
                       "Array handle equals different array.");
    }
  };

  struct TryArrayHandleType
    {
    void operator()() const
      {
      vtkm::testing::Testing::TryAllTypes(VerifyEmptyArrays());
      vtkm::testing::Testing::TryAllTypes(VerifyUserAllocatedHandle());
      vtkm::testing::Testing::TryAllTypes(VerifyVTKMAllocatedHandle());
      }
    };

public:
  static VTKM_CONT_EXPORT int Run()
  {
    return vtkm::cont::testing::Testing::Run(TryArrayHandleType());
  }
};

}
}
} // namespace vtkm::cont::testing

#endif //vtk_m_cont_testing_TestingArrayHandles_h
