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
#ifndef vtk_m_cont_internal_DeviceAdapterAlgorithmHPX_h
#define vtk_m_cont_internal_DeviceAdapterAlgorithmHPX_h

// include HPX headers before vtkm+boost to avoid problems with definitions
#include <vtkm/cont/hpx/internal/DeviceAdapterTagHPX.h>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/include/parallel_scan.hpp> 
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

  /*
      std::cout << "\nInput values " ;
      std::copy(
        vtkm::cont::ArrayPortalToIteratorEnd(inputPortal)-10,
        vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
        std::ostream_iterator<T>(std::cout, ", ")
      );
  */
      internal::WrappedBinaryOperator<T, BinaryFunctor> wrappedOp( binary_functor );

      hpx::parallel::inclusive_scan(hpx::parallel::par,
        vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
        vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
        vtkm::cont::ArrayPortalToIteratorBegin(outputPortal),
        T(),
        wrappedOp);

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

    // vtkm::cont::ArrayPortalToIterators<PortalOut>::IteratorType fullValue = 
    hpx::parallel::exclusive_scan(hpx::parallel::par, 
      vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      vtkm::cont::ArrayPortalToIteratorBegin(outputPortal), T());

    result =  outputPortal.Get(numberOfValues - 1) + inputPortal.Get(numberOfValues - 1);
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
    // vtkm::cont::ArrayPortalToIterators<PortalOut>::IteratorType fullValue =
    hpx::parallel::exclusive_scan(hpx::parallel::par,
      vtkm::cont::ArrayPortalToIteratorBegin(inputPortal),
      vtkm::cont::ArrayPortalToIteratorEnd(inputPortal),
      vtkm::cont::ArrayPortalToIteratorBegin(outputPortal),
      initialValue,
      wrappedOp);

    result =  outputPortal.Get(numberOfValues - 1) + inputPortal.Get(numberOfValues - 1);
    return result;
  }

  //----------------------------------------------------------------------------
  // helper for ReduceByKey
  struct ReduceKeySeriesStates
  {
    bool fStart;    // START of a segment
    bool fEnd;      // END of a segment
    ReduceKeySeriesStates(bool start=false, bool end=false) : fStart(start), fEnd(end) {}
  };

  //----------------------------------------------------------------------------
  // helper for ReduceByKey
  template<typename InputPortalType>
  struct ReduceStencilGeneration : vtkm::exec::FunctorBase
  {
    typedef typename vtkm::cont::ArrayHandle< ReduceKeySeriesStates >::template ExecutionTypes<DeviceAdapterTagHPX>
    ::Portal KeyStatePortalType;

    InputPortalType Input;
    KeyStatePortalType KeyState;

    VTKM_CONT_EXPORT
    ReduceStencilGeneration(const InputPortalType &input,
                            const KeyStatePortalType &kstate)
    : Input(input),
    KeyState(kstate)
    {  }

    VTKM_EXEC_EXPORT
    void operator()(vtkm::Id centerIndex) const
    {
      typedef typename InputPortalType::ValueType ValueType;
      typedef typename KeyStatePortalType::ValueType KeyStateType;

      const vtkm::Id leftIndex = centerIndex - 1;
      const vtkm::Id rightIndex = centerIndex + 1;

      //we need to determine which of three states this
      //index is. It can be:
      // 1. Middle of a set of equivalent keys.
      // 2. Start of a set of equivalent keys.
      // 3. End of a set of equivalent keys.
      // 4. Both the start and end of a set of keys

      //we don't have to worry about an array of length 1, as
      //the calling code handles that use case

      if(centerIndex == 0)
      {
        //this means we are at the start of the array
        //means we are automatically START
        //just need to check if we are END
        const ValueType centerValue = this->Input.Get(centerIndex);
        const ValueType rightValue = this->Input.Get(rightIndex);
        const KeyStateType state = ReduceKeySeriesStates(true, rightValue != centerValue);
        this->KeyState.Set(centerIndex, state);
        std::cout << " Set state # " << true << (rightValue != centerValue) << std::endl;
      }
      else if(rightIndex == this->Input.GetNumberOfValues())
      {
        //this means we are at the end, so we are at least END
        //just need to check if we are START
        const ValueType centerValue = this->Input.Get(centerIndex);
        const ValueType leftValue = this->Input.Get(leftIndex);
        const KeyStateType state = ReduceKeySeriesStates(leftValue != centerValue, true);
        this->KeyState.Set(centerIndex, state);
        std::cout << " Set state # " << (leftValue != centerValue) << true << std::endl;
      }
      else
      {
        const ValueType centerValue = this->Input.Get(centerIndex);
        const bool leftMatches(this->Input.Get(leftIndex) == centerValue);
        const bool rightMatches(this->Input.Get(rightIndex) == centerValue);

        //assume it is the middle, and check for the other use-case
        KeyStateType state = ReduceKeySeriesStates(!leftMatches, !rightMatches);
        this->KeyState.Set(centerIndex, state);
        std::cout << " Set state # " << !leftMatches << !rightMatches << std::endl;
      }
    }
  };

  //----------------------------------------------------------------------------
  // helper for ReduceByKey
  template<typename BinaryFunctor>
  struct ReduceByKeyAdd
  {
    BinaryFunctor BinaryOperator;

    ReduceByKeyAdd(BinaryFunctor binary_functor):
    BinaryOperator( binary_functor )
    { }

    template<typename T>
    vtkm::Pair<T, ReduceKeySeriesStates> operator()(const vtkm::Pair<T, ReduceKeySeriesStates>& a,
                                                    const vtkm::Pair<T, ReduceKeySeriesStates>& b) const
    {
      typedef vtkm::Pair<T, ReduceKeySeriesStates> ReturnType;
      //need too handle how we are going to add two numbers together
      //based on the keyStates that they have

      // Make it work for parallel inclusive scan.  Will end up with all start bits = 1
      // the following logic should change if you use a different parallel scan algorithm.
      if (!b.second.fStart) {
        std::cout << "not second " << std::endl;
        // if b is not START, then it's safe to sum a & b.
        // Propagate a's start flag to b
        // so that later when b's START bit is set, it means there must exists a START between a and b
        return ReturnType(this->BinaryOperator(a.first , b.first),
                          ReduceKeySeriesStates(a.second.fStart, b.second.fEnd));
      }
      std::cout << "second " << std::endl;
      return b;
    }

  };

  //----------------------------------------------------------------------------
  // helper for ReduceByKey
  struct ReduceByKeyUnaryStencilOp
  {
    bool operator()(ReduceKeySeriesStates keySeriesState) const
    {
      return keySeriesState.fEnd;
    }
    
  };

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
    VTKM_ASSERT_CONT(keys.GetNumberOfValues() == values.GetNumberOfValues());
    const vtkm::Id numberOfKeys = keys.GetNumberOfValues();
    std::cout << "Using custom reduce by key algorithm " << std::endl;

    if(numberOfKeys <= 1)
    { //we only have a single key/value so that is our output
      Copy(keys, keys_output);
      Copy(values, values_output);
      return;
    }

    //we need to determine based on the keys what is the keystate for
    //each key. The states are start, middle, end of a series and the special
    //state start and end of a series
    vtkm::cont::ArrayHandle< ReduceKeySeriesStates > keystate;
    {
      typedef typename vtkm::cont::ArrayHandle<T,KIn>::template ExecutionTypes<DeviceAdapterTagHPX>
      ::PortalConst InputPortalType;

      typedef typename vtkm::cont::ArrayHandle< ReduceKeySeriesStates >::template ExecutionTypes<DeviceAdapterTagHPX>
      ::Portal KeyStatePortalType;

      InputPortalType inputPortal = keys.PrepareForInput(DeviceAdapterTagHPX());
      KeyStatePortalType keyStatePortal = keystate.PrepareForOutput(numberOfKeys,
                                                                    DeviceAdapterTagHPX());
      ReduceStencilGeneration<InputPortalType> kernel(inputPortal, keyStatePortal);
      Schedule(kernel, numberOfKeys);
    }

    //next step is we need to reduce the values for each key. This is done
    //by running an inclusive scan over the values array using the stencil.
    //
    // this inclusive scan will write out two values, the first being
    // the value summed currently, the second being 0 or 1, with 1 being used
    // when this is a value of a key we need to write ( END or START_AND_END)
    {
      typedef vtkm::cont::ArrayHandle<U,VIn> ValueInHandleType;
      typedef vtkm::cont::ArrayHandle<U,VOut> ValueOutHandleType;
      typedef vtkm::cont::ArrayHandle< ReduceKeySeriesStates> StencilHandleType;
      typedef vtkm::cont::ArrayHandleZip<ValueInHandleType,
      StencilHandleType> ZipInHandleType;
      typedef vtkm::cont::ArrayHandleZip<ValueOutHandleType,
      StencilHandleType> ZipOutHandleType;

      StencilHandleType stencil;
      ValueOutHandleType reducedValues;

      ZipInHandleType scanInput( values, keystate);
      ZipOutHandleType scanOutput( reducedValues, stencil);

      ScanInclusive(scanInput,
                    scanOutput,
                    ReduceByKeyAdd<BinaryFunctor>(binary_functor) );

      //at this point we are done with keystate, so free the memory
      keystate.ReleaseResources();

      // all we need know is an efficient way of doing the write back to the
      // reduced global memory. this is done by using StreamCompact with the
      // stencil and values we just created with the inclusive scan
      StreamCompact( reducedValues,
                                      stencil,
                                      values_output,
                                      ReduceByKeyUnaryStencilOp());
      
    } //release all temporary memory
    
    
    //find all the unique keys
    Copy(keys,keys_output);
    Unique(keys_output);
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

    std::cout << "Entering Sort algorithm" << std::endl;
    PortalType arrayPortal = values.PrepareForInPlace(Device());
    vtkm::cont::ArrayPortalToIterators<PortalType> iterators(arrayPortal);
    std::sort(iterators.GetBegin(), iterators.GetEnd());
  }

  //----------------------------------------------------------------------------
  template<typename T, class Storage, class BinaryFunctor>
  VTKM_CONT_EXPORT static void Sort(vtkm::cont::ArrayHandle<T,Storage>& values,
          BinaryFunctor binary_functor)
  {
    typedef typename vtkm::cont::ArrayHandle<T,Storage>
        ::template ExecutionTypes<Device>::Portal PortalType;
std::cout << "Entering Sort algorithm (functor)" << std::endl;
    internal::WrappedBinaryOperator<bool, BinaryFunctor> wrappedOp( binary_functor );

    PortalType arrayPortal = values.PrepareForInPlace(Device());
    vtkm::cont::ArrayPortalToIterators<PortalType> iterators(arrayPortal);
    std::sort(iterators.GetBegin(), iterators.GetEnd(), wrappedOp);
  }

  //----------------------------------------------------------------------------
  VTKM_CONT_EXPORT static void Synchronize()
  {
    // Nothing to do. This device is HPX and has no asynchronous operations.
  }

};

}
} // namespace vtkm::cont

#endif //vtk_m_cont_internal_DeviceAdapterAlgorithmHPX_h
