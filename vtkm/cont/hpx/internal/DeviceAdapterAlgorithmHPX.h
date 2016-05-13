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

#define PREFIX_SCAN

// include HPX headers before vtkm+boost to avoid problems with definitions
#include <vtkm/cont/hpx/internal/DeviceAdapterTagHPX.h>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/algorithms/inclusive_scan.hpp>
#include <hpx/parallel/algorithms/exclusive_scan.hpp>
#include <hpx/parallel/algorithms/sort.hpp>
#include <hpx/parallel/algorithms/sort_by_key.hpp>
#include <hpx/parallel/algorithms/reduce.hpp>
#include <hpx/parallel/algorithms/reduce_by_key.hpp>
#include <hpx/parallel/algorithms/copy.hpp>
#ifdef PREFIX_SCAN
# include <hpx/parallel/algorithms/prefix_scan.hpp>
# include <hpx/parallel/algorithms/prefix_copy_if.hpp>
# include <hpx/parallel/algorithms/prefix_unique.hpp>
# define SCAN_I prefix_scan_inclusive
# define SCAN_E prefix_scan_exclusive
#else
# define SCAN_I inclusive_scan
# define SCAN_E exclusive_scan
#endif
#include <hpx/parallel/util/zip_iterator.hpp>
//
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>
#include <vtkm/cont/ErrorExecution.h>
#include <vtkm/cont/internal/DeviceAdapterAlgorithmGeneral.h>

#include <vtkm/exec/internal/ErrorMessageBuffer.h>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>

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
    vtkm::Id numberOfValues = input.GetNumberOfValues();

    T result = T();
    if (numberOfValues <= 0) { return result; }

    auto inputPortal = input.PrepareForInput(Device());
    auto outputPortal = output.PrepareForOutput(numberOfValues, Device());

    hpx::parallel::SCAN_I(hpx::parallel::par,
      vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      vtkm::cont::ArrayPortalToIteratorBegin(outputPortal));

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
      vtkm::Id numberOfValues = input.GetNumberOfValues();

      T result = T();
      if (numberOfValues <= 0) { return result; }

      auto inputPortal = input.PrepareForInput(Device());
      auto outputPortal = output.PrepareForOutput(numberOfValues, Device());

      typedef internal::WrappedBinaryOperator<T, BinaryFunctor> wrapped_type;
      wrapped_type wrappedOp( binary_functor );
      hpx::parallel::SCAN_I(hpx::parallel::par,
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
    vtkm::Id numberOfValues = input.GetNumberOfValues();

    T result = T();
    if (numberOfValues <= 0) { return result; }

    auto inputPortal = input.PrepareForInput(Device());
    auto outputPortal = output.PrepareForOutput(numberOfValues, Device());

    // the calculation is 'in place' so get this value before it is overwritten
    T temp = inputPortal.Get(numberOfValues - 1);

    // vtkm::cont::ArrayPortalToIterators<PortalOut>::IteratorType fullValue =
    hpx::parallel::SCAN_E(hpx::parallel::par,
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
    vtkm::Id numberOfValues = input.GetNumberOfValues();

    T result = T();
    if (numberOfValues <= 0) { return result; }

    auto inputPortal = input.PrepareForInput(Device());
    auto outputPortal = output.PrepareForOutput(numberOfValues, Device());

    internal::WrappedBinaryOperator<T, BinaryFunctor> wrappedOp( binary_functor );

    // the calculation is 'in place' so get this value before it is overwritten
    T temp = inputPortal.Get(numberOfValues - 1);

    hpx::parallel::SCAN_E(hpx::parallel::par,
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
      return Reduce(input,
          std::forward<T>(initialValue),
          std::plus<T>());
   }

  //----------------------------------------------------------------------------
  template<typename T, class SIn, class BinaryFunctor>
   VTKM_CONT_EXPORT static T Reduce(
       const vtkm::cont::ArrayHandle<T,SIn> &input,
       T initialValue,
       BinaryFunctor binary_functor)
   {
     const vtkm::Id numberOfValues = input.GetNumberOfValues();
     if (numberOfValues <= 0) { return initialValue; }

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
    if (numberOfKeys <= 0) { return; }
    if (numberOfKeys <= 1)
    { //we only have a single key/value so that is our output
      Copy(keys, keys_output);
      Copy(values, values_output);
      return;
    }

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

    std::size_t reduced_size = std::distance(vtkm::cont::ArrayPortalToIteratorBegin(outputPortal_k), result.first);

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

      try {
          hpx::parallel::for_each(
              hpx::parallel::par,
              ::boost::counting_iterator<vtkm::Id>(0),
              ::boost::counting_iterator<vtkm::Id>(numInstances),
              kernel);
      }
      catch(hpx::exception const& e) {
          throw vtkm::cont::ErrorExecution(e.what());
      }
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

  //--------------------------------------------------------------------------
  // Sort
  template<typename T, class Storage>
  VTKM_CONT_EXPORT static void Sort(vtkm::cont::ArrayHandle<T,Storage>& values)
  {
    Sort(values, std::less< T >() );
  }

  //----------------------------------------------------------------------------
  template<typename T, class Storage, class BinaryPredicate>
  VTKM_CONT_EXPORT static void Sort(vtkm::cont::ArrayHandle<T,Storage>& values,
      BinaryPredicate binary_compare)
  {
      auto arrayPortal = values.PrepareForInPlace(Device());

      internal::WrappedBinaryOperator<bool,BinaryPredicate>
          wrappedCompare(binary_compare);

      hpx::parallel::sort(hpx::parallel::par,
          vtkm::cont::ArrayPortalToIteratorBegin(arrayPortal),
          vtkm::cont::ArrayPortalToIteratorEnd(arrayPortal),
          wrappedCompare);
  }

  //--------------------------------------------------------------------------
  // SortByKey
  template<typename T, typename U, class StorageT,  class StorageU>
  VTKM_CONT_EXPORT static void SortByKey(
      vtkm::cont::ArrayHandle<T,StorageT> &keys,
      vtkm::cont::ArrayHandle<U,StorageU> &values)
  {
    SortByKey(keys, values, std::less<T>());
  }

  //----------------------------------------------------------------------------
  template<typename T, typename U,
           class StorageT, class StorageU,
           class BinaryPredicate>
  VTKM_CONT_EXPORT static void SortByKey(
      vtkm::cont::ArrayHandle<T,StorageT>& keys,
      vtkm::cont::ArrayHandle<U,StorageU>& values,
      BinaryPredicate binary_compare)
  {

      auto arrayPortal_k = keys.PrepareForInPlace(Device());
      auto arrayPortal_v = values.PrepareForInPlace(Device());

      internal::WrappedBinaryOperator<bool,BinaryPredicate>
          wrappedCompare(binary_compare);

      hpx::parallel::sort_by_key(hpx::parallel::par,
          vtkm::cont::ArrayPortalToIteratorBegin(arrayPortal_k),
          vtkm::cont::ArrayPortalToIteratorEnd(arrayPortal_k),
          vtkm::cont::ArrayPortalToIteratorBegin(arrayPortal_v),
          wrappedCompare);
  }

  //----------------------------------------------------------------------------
  VTKM_CONT_EXPORT static void Synchronize()
  {
    // @TODO. would like to add support for futures
  }

#ifdef PREFIX_SCAN
  //--------------------------------------------------------------------------
  // Stream Compact
  template<typename T, typename U, class CIn, class CStencil,
           class COut, class UnaryPredicate>
  VTKM_CONT_EXPORT static void StreamCompact(
      const vtkm::cont::ArrayHandle<T,CIn>& input,
      const vtkm::cont::ArrayHandle<U,CStencil>& stencil,
      vtkm::cont::ArrayHandle<T,COut>& output,
      UnaryPredicate unary_predicate)
  {
    VTKM_ASSERT(input.GetNumberOfValues() == stencil.GetNumberOfValues());
    vtkm::Id arraySize = stencil.GetNumberOfValues();

    // portal types for the input/output arrays
    typedef const typename vtkm::cont::ArrayHandle<T,CIn>
        ::template ExecutionTypes<Device>::PortalConst PortalIn;
    typedef const typename vtkm::cont::ArrayHandle<U,CStencil>
        ::template ExecutionTypes<Device>::PortalConst PortalStencil;
    typedef typename vtkm::cont::ArrayHandle<T,COut>
        ::template ExecutionTypes<Device>::Portal PortalOut;

    // input portal vars
    PortalIn inputPortal_i = input.PrepareForInput(Device());
    PortalStencil inputPortal_s = stencil.PrepareForInput(Device());

    // output portal var
    vtkm::cont::ArrayHandle<T,COut> temp_output;
    PortalOut tempPortal = temp_output.PrepareForOutput(arraySize, Device());

    typedef typename vtkm::cont::ArrayPortalToIterators<PortalIn>::IteratorType
        InIteratorType;
    typedef typename vtkm::cont::ArrayPortalToIterators<PortalStencil>::IteratorType
        StencilIteratorType;
    typedef typename vtkm::cont::ArrayPortalToIterators<PortalOut>::IteratorType
        OutIteratorType;
    //
    typedef typename hpx::util::zip_iterator<InIteratorType, StencilIteratorType> zip_type;
    typedef typename zip_type::reference zip_ref;
    typedef typename zip_type::value_type zip_value;

    InIteratorType iterators_i(inputPortal_i);
    StencilIteratorType iterators_s(inputPortal_s);

    auto end = hpx::parallel::prefix_copy_if_stencil(
        hpx::parallel::par,
        // begin
        vtkm::cont::ArrayPortalToIteratorBegin(inputPortal_i),
        vtkm::cont::ArrayPortalToIteratorEnd(inputPortal_i),
        vtkm::cont::ArrayPortalToIteratorBegin(inputPortal_s),
        vtkm::cont::ArrayPortalToIteratorBegin(tempPortal),
        unary_predicate
    );

    vtkm::Id finalSize = std::distance(vtkm::cont::ArrayPortalToIteratorBegin(tempPortal), end);
    temp_output.Shrink(finalSize);
    auto newPortal = temp_output.PrepareForInput(Device());
    auto outputPortal   = output.PrepareForOutput(finalSize, Device());
    //
    hpx::parallel::copy(
        hpx::parallel::par,
        vtkm::cont::ArrayPortalToIteratorBegin(newPortal),
        vtkm::cont::ArrayPortalToIteratorEnd(newPortal),
        vtkm::cont::ArrayPortalToIteratorBegin(outputPortal)
    );
  }

  template<typename T, typename U, class CIn, class CStencil, class COut>
  VTKM_CONT_EXPORT static void StreamCompact(
      /*const*/ vtkm::cont::ArrayHandle<T,CIn>& input,
      /*const*/ vtkm::cont::ArrayHandle<U,CStencil>& stencil,
      vtkm::cont::ArrayHandle<T,COut>& output)
  {
    ::vtkm::NotZeroInitialized unary_predicate;
    StreamCompact(input, stencil, output, unary_predicate);
  }

  template<typename T, class CStencil, class COut>
  VTKM_CONT_EXPORT static void StreamCompact(
      /*const*/ vtkm::cont::ArrayHandle<T,CStencil> &stencil,
      vtkm::cont::ArrayHandle<vtkm::Id,COut> &output)
  {
    vtkm::cont::ArrayHandleIndex input(stencil.GetNumberOfValues());
    StreamCompact(input, stencil, output);
  }

  //--------------------------------------------------------------------------
  // Unique
  template<typename T, class Storage>
  VTKM_CONT_EXPORT static void Unique(
      vtkm::cont::ArrayHandle<T,Storage> &values)
  {
    Unique(values, std::equal_to<T>());
  }

  template<typename T, class Storage, class BinaryCompare>
  VTKM_CONT_EXPORT static void Unique(
      vtkm::cont::ArrayHandle<T,Storage> &values,
      BinaryCompare binary_compare)
  {
      auto arrayPortal = values.PrepareForInPlace(Device());

      internal::WrappedBinaryOperator<bool,BinaryCompare>
          wrappedCompare(binary_compare);

      auto new_end = hpx::parallel::prefix_unique(hpx::parallel::par,
          vtkm::cont::ArrayPortalToIteratorBegin(arrayPortal),
          vtkm::cont::ArrayPortalToIteratorEnd(arrayPortal),
          wrappedCompare);

      std::size_t N = std::distance(
          vtkm::cont::ArrayPortalToIteratorBegin(arrayPortal), new_end);
      values.Shrink(N);
  }

#endif

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
