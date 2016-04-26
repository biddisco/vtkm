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

#ifndef vtk_m_cont_internal_DeviceAdapterAlgorithmHPX_h
#define vtk_m_cont_internal_DeviceAdapterAlgorithmHPX_h

// include HPX headers before vtkm+boost to avoid problems with definitions
#include <vtkm/cont/hpx/internal/DeviceAdapterTagHPX.h>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/algorithms/inclusive_scan.hpp>
#include <hpx/parallel/algorithms/exclusive_scan.hpp>
#include <hpx/parallel/algorithms/sort.hpp>
#include <hpx/parallel/algorithms/reduce.hpp>
#include <hpx/parallel/algorithms/reduce_by_key.hpp>
#include <hpx/parallel/algorithms/copy.hpp>
//
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>
#include <vtkm/cont/ErrorExecution.h>
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>

#include <vtkm/exec/internal/ErrorMessageBuffer.h>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/utility/enable_if.hpp>


#include <algorithm>
#include <numeric>

namespace vtkm {
namespace cont {

    template <typename Op>
    struct initial_value
    {
        template <typename T>
        constexpr static T call() { return T(); }
    };

    template <typename T, typename F>
    struct initial_value<vtkm::cont::internal::WrappedBinaryOperator<T, F > >
    {
        constexpr static T call() { return T(); }
    };

    template <typename T>
    struct initial_value<vtkm::cont::internal::WrappedBinaryOperator<T, vtkm::internal::Multiply> >
    {
        constexpr static T call() { return T(1); }
    };

    template <typename T>
    struct initial_value<vtkm::cont::internal::WrappedBinaryOperator<T, std::multiplies<T> > >
    {
        constexpr static T call() { return T(1); }
    };

template<>
struct DeviceAdapterAlgorithm<vtkm::cont::DeviceAdapterTagHPX> :
    vtkm::cont::internal::DeviceAdapterAlgorithmGeneral<
        DeviceAdapterAlgorithm<vtkm::cont::DeviceAdapterTagHPX>,
        vtkm::cont::DeviceAdapterTagHPX>
{
private:
  typedef vtkm::cont::DeviceAdapterTagHPX Device;

public:

  //----------------------------------------------------------------------------
  template<typename T, class CIn, class COut>
  VTKM_CONT_EXPORT static T ScanInclusive(
      const vtkm::cont::ArrayHandle<T,CIn> &input,
      vtkm::cont::ArrayHandle<T,COut>& output)
  {
    typedef typename vtkm::cont::ArrayHandle<T,COut>
        ::template ExecutionTypes<Device>::Portal PortalOut;
    typedef typename vtkm::cont::ArrayHandle<T,CIn>
        ::template ExecutionTypes<Device>::PortalConst PortalIn;

    vtkm::Id numberOfValues = input.GetNumberOfValues();

    PortalIn inputPortal = input.PrepareForInput(Device());
    PortalOut outputPortal = output.PrepareForOutput(numberOfValues, Device());

    T result = T();
    if (numberOfValues <= 0) { return result; }

/*
    std::cout << "\nInput values " ;
    std::copy(
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal)-10,
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      std::ostream_iterator<T>(std::cout, ", ")
    );
*/

    hpx::parallel::inclusive_scan(hpx::parallel::par,
      vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      vtkm::cont::ArrayPortalToIteratorBegin(outputPortal));

/*
    std::cout << "\nOutput values " ;
    std::copy(
      vtkm::cont::ArrayPortalToIteratorEnd(outputPortal)-10,
      vtkm::cont::ArrayPortalToIteratorEnd(outputPortal),
      std::ostream_iterator<T>(std::cout, ", ")
    );
*/
    result =  outputPortal.Get(numberOfValues - 1);
    return result;
  }

  //----------------------------------------------------------------------------
  template<typename T, class CIn, class COut, class BinaryFunctor>
  VTKM_CONT_EXPORT static T ScanInclusive(
      const vtkm::cont::ArrayHandle<T,CIn> &input,
      vtkm::cont::ArrayHandle<T,COut> &output,
      BinaryFunctor binary_functor)
  {
      typedef typename vtkm::cont::ArrayHandle<T,COut>
          ::template ExecutionTypes<Device>::Portal PortalOut;
      typedef typename vtkm::cont::ArrayHandle<T,CIn>
          ::template ExecutionTypes<Device>::PortalConst PortalIn;

      vtkm::Id numberOfValues = input.GetNumberOfValues();

      PortalIn inputPortal = input.PrepareForInput(Device());
      PortalOut outputPortal = output.PrepareForOutput(numberOfValues, Device());

      T result = T();
      if (numberOfValues <= 0) { return result; }

      typedef internal::WrappedBinaryOperator<T, BinaryFunctor> wrapped_type;
      wrapped_type wrappedOp( binary_functor );
      hpx::parallel::inclusive_scan(hpx::parallel::par,
        vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
        vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
        vtkm::cont::ArrayPortalToIteratorBegin(outputPortal),
        initial_value<wrapped_type>().call(),
        wrappedOp);

      result =  outputPortal.Get(numberOfValues - 1);
      return result;
}

  //----------------------------------------------------------------------------
  template<typename T, class CIn, class COut>
  VTKM_CONT_EXPORT static T ScanExclusive(
      const vtkm::cont::ArrayHandle<T,CIn> &input,
      vtkm::cont::ArrayHandle<T,COut>& output)
  {
    typedef typename vtkm::cont::ArrayHandle<T,COut>
        ::template ExecutionTypes<Device>::Portal PortalOut;
    typedef typename vtkm::cont::ArrayHandle<T,CIn>
        ::template ExecutionTypes<Device>::PortalConst PortalIn;

    vtkm::Id numberOfValues = input.GetNumberOfValues();

    PortalIn inputPortal = input.PrepareForInput(Device());
    PortalOut outputPortal = output.PrepareForOutput(numberOfValues, Device());

    T result = T();
    if (numberOfValues <= 0) { return result; }

    // the calculation is 'in place' so get this value before it is overwritten
    T temp = inputPortal.Get(numberOfValues - 1);

    // vtkm::cont::ArrayPortalToIterators<PortalOut>::IteratorType fullValue =
    hpx::parallel::exclusive_scan(hpx::parallel::par,
      vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      vtkm::cont::ArrayPortalToIteratorBegin(outputPortal), T());

    result =  T(outputPortal.Get(numberOfValues - 1)) + temp;
    return result;
  }

  //----------------------------------------------------------------------------
  template<typename T, class CIn, class COut, class BinaryFunctor>
  VTKM_CONT_EXPORT static T ScanExclusive(
      const vtkm::cont::ArrayHandle<T,CIn> &input,
      vtkm::cont::ArrayHandle<T,COut> &output,
      BinaryFunctor binary_functor,
      const T& initialValue)
  {
    typedef typename vtkm::cont::ArrayHandle<T,COut>
        ::template ExecutionTypes<Device>::Portal PortalOut;
    typedef typename vtkm::cont::ArrayHandle<T,CIn>
        ::template ExecutionTypes<Device>::PortalConst PortalIn;

    vtkm::Id numberOfValues = input.GetNumberOfValues();

    PortalIn inputPortal = input.PrepareForInput(Device());
    PortalOut outputPortal = output.PrepareForOutput(numberOfValues, Device());

    T result = T();
    if (numberOfValues <= 0) { return result; }

    internal::WrappedBinaryOperator<T, BinaryFunctor> wrappedOp( binary_functor );

    // the calculation is 'in place' so get this value before it is overwritten
    T temp = inputPortal.Get(numberOfValues - 1);

    // vtkm::cont::ArrayPortalToIterators<PortalOut>::IteratorType fullValue =
    hpx::parallel::exclusive_scan(hpx::parallel::par,
      vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      vtkm::cont::ArrayPortalToIteratorBegin(outputPortal),
      initialValue,
      wrappedOp);

    result =  binary_functor(outputPortal.Get(numberOfValues - 1), temp);
    return result;
  }

  //----------------------------------------------------------------------------
  template<typename T, class SIn>
   VTKM_CONT_EXPORT static T Reduce(
       const vtkm::cont::ArrayHandle<T,SIn> &input,
       T initialValue)
   {
     const vtkm::Id numberOfValues = input.GetNumberOfValues();
     if (numberOfValues <= 0)
       {
       return initialValue;
       }

     auto inputPortal = input.PrepareForInput(Device());

     return hpx::parallel::reduce(
         hpx::parallel::par,
         vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
         vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
         std::forward<T>(initialValue));
   }

  //----------------------------------------------------------------------------
  template<typename T, class SIn, class BinaryFunctor>
   VTKM_CONT_EXPORT static T Reduce(
       const vtkm::cont::ArrayHandle<T,SIn> &input,
       T initialValue,
       BinaryFunctor binary_functor)
   {
     const vtkm::Id numberOfValues = input.GetNumberOfValues();
     if (numberOfValues <= 0)
       {
       return initialValue;
       }

     auto inputPortal = input.PrepareForInput(Device());

     return hpx::parallel::reduce(
         hpx::parallel::par,
         vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
         vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
         std::forward<T>(initialValue),
         std::forward<BinaryFunctor>(binary_functor));
   }

  //----------------------------------------------------------------------------
  template<typename T, typename U, class KIn, class VIn, class KOut, class VOut,
  class BinaryFunctor>
  VTKM_CONT_EXPORT static void ReduceByKey(
                                           const vtkm::cont::ArrayHandle<T,KIn> &keys,
                                           const vtkm::cont::ArrayHandle<U,VIn> &values,
                                           vtkm::cont::ArrayHandle<T,KOut> &keys_output,
                                           vtkm::cont::ArrayHandle<U,VOut> &values_output,
                                           BinaryFunctor binary_functor)
  {
    VTKM_ASSERT(keys.GetNumberOfValues() == values.GetNumberOfValues());

    const vtkm::Id numberOfKeys = keys.GetNumberOfValues();
    if (numberOfKeys <= 0)
    {
      return;
    }
    if (numberOfKeys <= 1)
    { //we only have a single key/value so that is our output
      Copy(keys, keys_output);
      Copy(values, values_output);
      return;
    }

    typedef const typename vtkm::cont::ArrayHandle<T,KIn>
        ::template ExecutionTypes<Device>::Portal PortalIn_k;
    typedef const typename vtkm::cont::ArrayHandle<U,VIn>
        ::template ExecutionTypes<Device>::Portal PortalIn_v;
    typedef typename vtkm::cont::ArrayHandle<T,KOut>
        ::template ExecutionTypes<Device>::PortalConst PortalOut_k;
    typedef typename vtkm::cont::ArrayHandle<U,VOut>
        ::template ExecutionTypes<Device>::PortalConst PortalOut_v;

    auto inputPortal_k = keys.PrepareForInput(Device());
    auto inputPortal_v = values.PrepareForInput(Device());
    auto outputPortal_k = keys_output.PrepareForOutput(numberOfKeys, Device());
    auto outputPortal_v = values_output.PrepareForOutput(numberOfKeys, Device());

    auto result =
        hpx::parallel::reduce_by_key(
            hpx::parallel::par,
            vtkm::cont::ArrayPortalToIteratorBegin(inputPortal_k),
            vtkm::cont::ArrayPortalToIteratorEnd(inputPortal_k),
            vtkm::cont::ArrayPortalToIteratorBegin(inputPortal_v),
            //
            vtkm::cont::ArrayPortalToIteratorBegin(outputPortal_k),
            vtkm::cont::ArrayPortalToIteratorBegin(outputPortal_v),
            std::equal_to<T>(),
            std::forward<BinaryFunctor>(binary_functor)
            );

    int reduced_size = std::distance(vtkm::cont::ArrayPortalToIteratorBegin(outputPortal_k), result.first);

    keys_output.Shrink( reduced_size );
    values_output.Shrink( reduced_size );
  }

  //--------------------------------------------------------------------------
  // Copy
  template<typename T, typename U, class CIn, class COut>
  VTKM_CONT_EXPORT static void Copy(const vtkm::cont::ArrayHandle<T, CIn> &input,
                                    vtkm::cont::ArrayHandle<U, COut> &output)
  {
    vtkm::Id arraySize = input.GetNumberOfValues();

    auto inputPortal = input.PrepareForInput(Device());
    auto oututPortal = output.PrepareForOutput(arraySize, Device());

    hpx::parallel::copy(
        hpx::parallel::par,
        vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
        vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
        vtkm::cont::ArrayPortalToIteratorBegin(oututPortal));
  }

private:
  //----------------------------------------------------------------------------
  // This runs in the execution environment.
  template<class FunctorType>
  class ScheduleKernel
  {
  public:
    ScheduleKernel(const FunctorType &functor)
      : Functor(functor) {  }

    //needed for when calling from schedule on a range
    VTKM_EXEC_EXPORT void operator()(vtkm::Id index) const
    {
      this->Functor(index);
    }

  private:
    const FunctorType Functor;
  };

public:
  //----------------------------------------------------------------------------
  template<class Functor>
  VTKM_CONT_EXPORT static void Schedule(Functor functor,
                                        vtkm::Id numInstances)
  {
    const vtkm::Id MESSAGE_SIZE = 1024;
    char errorString[MESSAGE_SIZE];
    errorString[0] = '\0';
    vtkm::exec::internal::ErrorMessageBuffer
        errorMessage(errorString, MESSAGE_SIZE);

    functor.SetErrorMessageBuffer(errorMessage);

    DeviceAdapterAlgorithm<Device>::ScheduleKernel<Functor> kernel(functor);

    hpx::parallel::for_each(
          hpx::parallel::par,
          ::boost::counting_iterator<vtkm::Id>(0),
          ::boost::counting_iterator<vtkm::Id>(numInstances),
          kernel);

    if (errorMessage.IsErrorRaised())
    {
      throw vtkm::cont::ErrorExecution(errorString);
    }
  }

  //----------------------------------------------------------------------------
  template<class FunctorType>
  VTKM_CONT_EXPORT
  static void Schedule(FunctorType functor, vtkm::Id3 rangeMax)
  {
    DeviceAdapterAlgorithm<Device>::Schedule(functor,
                                     rangeMax[0] * rangeMax[1] * rangeMax[2] );
  }

  //----------------------------------------------------------------------------
  template<typename T, class Storage>
  VTKM_CONT_EXPORT static void Sort(vtkm::cont::ArrayHandle<T,Storage>& values)
  {
    typedef typename vtkm::cont::ArrayHandle<T,Storage>
        ::template ExecutionTypes<Device>::Portal PortalType;

    //this is required to get sort to work with zip handles
    std::less< T > lessOp;
    Sort(values, lessOp );
  }

  //----------------------------------------------------------------------------
  template<typename T, class Storage, class BinaryCompare>
  VTKM_CONT_EXPORT static void Sort(vtkm::cont::ArrayHandle<T,Storage>& values,
      BinaryCompare binary_compare)
  {
      typedef typename vtkm::cont::ArrayHandle<T,Storage>::template
          ExecutionTypes<Device>::Portal PortalType;

      PortalType arrayPortal = values.PrepareForInPlace(Device());

      typedef vtkm::cont::ArrayPortalToIterators<PortalType> IteratorsType;
      IteratorsType iterators(arrayPortal);

      internal::WrappedBinaryOperator<bool,BinaryCompare> wrappedCompare(binary_compare);

      hpx::parallel::sort(hpx::parallel::par, iterators.GetBegin(), iterators.GetEnd(),
          wrappedCompare);
  }

  //----------------------------------------------------------------------------
  VTKM_CONT_EXPORT static void Synchronize()
  {
    // @TODO. not sure what's needed here
  }

};

}} // namespace vtkm::cont

namespace vtkm {
namespace cont {

template<typename T>
class DeviceAdapterAtomicArrayImplementation<T,vtkm::cont::DeviceAdapterTagHPX>
{
public:
  VTKM_CONT_EXPORT
  DeviceAdapterAtomicArrayImplementation(
               vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic> handle):
    Iterators( IteratorsType( handle.PrepareForInPlace(
                                      vtkm::cont::DeviceAdapterTagHPX())
                             ) )
  {
  }

  VTKM_EXEC_EXPORT
  T Add(vtkm::Id index, const T& value) const
  {
    T* lockedValue;
    lockedValue = (Iterators.GetBegin()+index);
    return vtkmAtomicAdd(lockedValue, value);
  }

  VTKM_EXEC_EXPORT
  T CompareAndSwap(vtkm::Id index, const T& newValue, const T& oldValue) const
  {
    T* lockedValue;
    lockedValue = (Iterators.GetBegin()+index);
    return vtkmCompareAndSwap(lockedValue, newValue, oldValue);
  }

private:
  typedef typename vtkm::cont::ArrayHandle<T,vtkm::cont::StorageTagBasic>
        ::template ExecutionTypes<DeviceAdapterTagHPX>::Portal PortalType;
  typedef vtkm::cont::ArrayPortalToIterators<PortalType> IteratorsType;

  IteratorsType Iterators;

  VTKM_EXEC_EXPORT
  vtkm::Int32 vtkmAtomicAdd(vtkm::Int32 *address, const vtkm::Int32 &value) const
  {
    return __sync_fetch_and_add(address,value);
  }

  VTKM_EXEC_EXPORT
  vtkm::Int64 vtkmAtomicAdd(vtkm::Int64 *address, const vtkm::Int64 &value) const
  {
    return __sync_fetch_and_add(address,value);
  }

  VTKM_EXEC_EXPORT
  vtkm::UInt32 vtkmAtomicAdd(vtkm::UInt32 *address, const vtkm::UInt32 &value) const
  {
    return __sync_fetch_and_add(address,value);
  }

  VTKM_EXEC_EXPORT
  vtkm::UInt64 vtkmAtomicAdd(vtkm::UInt64 *address, const vtkm::UInt64 &value) const
  {
    return __sync_fetch_and_add(address,value);
  }
/*
  VTKM_EXEC_EXPORT
  vtkm::UInt32 vtkmAtomicAdd(vtkm::Float32 *address, const vtkm::Float32 &value) const
  {
    return __sync_fetch_and_add(address,value);
  }
*/
  VTKM_EXEC_EXPORT
  vtkm::Int32 vtkmCompareAndSwap(vtkm::Int32 *address, const vtkm::Int32 &newValue, const vtkm::Int32 &oldValue) const
  {
    return __sync_val_compare_and_swap(address,oldValue, newValue);
  }

  VTKM_EXEC_EXPORT
  vtkm::Int64 vtkmCompareAndSwap(vtkm::Int64 *address,const vtkm::Int64 &newValue, const vtkm::Int64 &oldValue) const
  {
    return __sync_val_compare_and_swap(address,oldValue,newValue);
  }

};
}}


#endif //vtk_m_cont_internal_DeviceAdapterAlgorithmHPX_h
