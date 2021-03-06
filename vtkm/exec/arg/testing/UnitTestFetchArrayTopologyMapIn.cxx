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

#include <vtkm/exec/arg/FetchTagArrayTopologyMapIn.h>

#include <vtkm/internal/FunctionInterface.h>
#include <vtkm/internal/Invocation.h>

#include <vtkm/testing/Testing.h>

namespace {

static const vtkm::Id ARRAY_SIZE = 10;

template<typename T>
struct TestPortal
{
  typedef T ValueType;

  VTKM_EXEC_CONT_EXPORT
  vtkm::Id GetNumberOfValues() const { return ARRAY_SIZE; }

  VTKM_EXEC_CONT_EXPORT
  ValueType Get(vtkm::Id index) const {
    VTKM_TEST_ASSERT(index >= 0, "Bad portal index.");
    VTKM_TEST_ASSERT(index < this->GetNumberOfValues(), "Bad portal index.");
    return TestValue(index, ValueType());
  }
};

struct NullParam {  };

template<vtkm::IdComponent InputDomainIndex,
         vtkm::IdComponent ParamIndex,
         typename T>
struct FetchArrayTopologyMapInTests
{

  template<typename Invocation>
  void TryInvocation(const Invocation &invocation) const
  {
    typedef vtkm::exec::arg::Fetch<
        vtkm::exec::arg::FetchTagArrayTopologyMapIn,
        vtkm::exec::arg::AspectTagDefault,
        Invocation,
        ParamIndex> FetchType;

    FetchType fetch;

    typename FetchType::ValueType value = fetch.Load(0, invocation);
    VTKM_TEST_ASSERT(value.GetNumberOfComponents() == 8,
                     "Topology fetch got wrong number of components.");

    VTKM_TEST_ASSERT(test_equal(value[0], TestValue(0, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[1], TestValue(1, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[2], TestValue(3, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[3], TestValue(2, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[4], TestValue(4, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[5], TestValue(5, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[6], TestValue(7, T())),
                     "Got invalid value from Load.");
    VTKM_TEST_ASSERT(test_equal(value[7], TestValue(6, T())),
                     "Got invalid value from Load.");
  }

  void operator()() const
  {
    std::cout << "Trying ArrayTopologyMapIn fetch on parameter " << ParamIndex
                 << " with type " << vtkm::testing::TypeName<T>::Name()
                 << std::endl;

    typedef vtkm::internal::FunctionInterface<
        void(NullParam,NullParam,NullParam,NullParam,NullParam)>
        BaseFunctionInterface;

    vtkm::internal::ConnectivityStructuredInternals<3> connectivityInternals;
    connectivityInternals.SetPointDimensions(vtkm::Id3(2,2,2));
    vtkm::exec::ConnectivityStructured<
        vtkm::TopologyElementTagPoint,vtkm::TopologyElementTagCell,3>
        connectivity(connectivityInternals);

    this->TryInvocation(vtkm::internal::make_Invocation<InputDomainIndex>(
                          BaseFunctionInterface()
                          .Replace<InputDomainIndex>(connectivity)
                          .template Replace<ParamIndex>(TestPortal<T>()),
                          NullParam(),
                          NullParam()));
  }

};

struct TryType
{
  template<typename T>
  void operator()(T) const
  {
    FetchArrayTopologyMapInTests<3,1,T>()();
    FetchArrayTopologyMapInTests<1,2,T>()();
    FetchArrayTopologyMapInTests<2,3,T>()();
    FetchArrayTopologyMapInTests<1,4,T>()();
    FetchArrayTopologyMapInTests<1,5,T>()();
  }
};

template<vtkm::IdComponent NumDimensions,
         vtkm::IdComponent ParamIndex,
         typename Invocation>
void TryStructuredPointCoordinatesInvocation(const Invocation &invocation)
{
  vtkm::exec::arg::Fetch<
      vtkm::exec::arg::FetchTagArrayTopologyMapIn,
      vtkm::exec::arg::AspectTagDefault,
      Invocation,
      ParamIndex> fetch;

  vtkm::Vec<vtkm::FloatDefault,3> origin =
      TestValue(0, vtkm::Vec<vtkm::FloatDefault,3>());
  vtkm::Vec<vtkm::FloatDefault,3> spacing =
      TestValue(1, vtkm::Vec<vtkm::FloatDefault,3>());

  vtkm::VecRectilinearPointCoordinates<NumDimensions> value =
      fetch.Load(0, invocation);
  VTKM_TEST_ASSERT(test_equal(value.GetOrigin(), origin), "Bad origin.");
  VTKM_TEST_ASSERT(test_equal(value.GetSpacing(), spacing), "Bad spacing.");

  origin[0] += spacing[0];
  value = fetch.Load(1, invocation);
  VTKM_TEST_ASSERT(test_equal(value.GetOrigin(), origin), "Bad origin.");
  VTKM_TEST_ASSERT(test_equal(value.GetSpacing(), spacing), "Bad spacing.");
}

template<vtkm::IdComponent NumDimensions>
void TryStructuredPointCoordinates(
    const vtkm::exec::ConnectivityStructured<
      vtkm::TopologyElementTagPoint,vtkm::TopologyElementTagCell,NumDimensions> &connectivity,
    const vtkm::internal::ArrayPortalUniformPointCoordinates &coordinates)
{
  typedef vtkm::internal::FunctionInterface<
      void(NullParam,NullParam,NullParam,NullParam,NullParam)>
      BaseFunctionInterface;

  // Try with topology in argument 1 and point coordinates in argument 2
  TryStructuredPointCoordinatesInvocation<NumDimensions,2>(
        vtkm::internal::make_Invocation<1>(
          BaseFunctionInterface()
          .Replace<1>(connectivity)
          .template Replace<2>(coordinates),
          NullParam(),
          NullParam())
        );
  // Try again with topology in argument 3 and point coordinates in argument 1
  TryStructuredPointCoordinatesInvocation<NumDimensions,1>(
        vtkm::internal::make_Invocation<3>(
          BaseFunctionInterface()
          .Replace<3>(connectivity)
          .template Replace<1>(coordinates),
          NullParam(),
          NullParam())
        );
}

void TryStructuredPointCoordinates()
{
  std::cout << "*** Fetching special case of uniform point coordinates. *****"
            << std::endl;

  vtkm::internal::ArrayPortalUniformPointCoordinates coordinates(
        vtkm::Id3(3,2,2),
        TestValue(0, vtkm::Vec<vtkm::FloatDefault,3>()),
        TestValue(1, vtkm::Vec<vtkm::FloatDefault,3>()));

  std::cout << "3D" << std::endl;
  vtkm::internal::ConnectivityStructuredInternals<3> connectivityInternals3d;
  connectivityInternals3d.SetPointDimensions(vtkm::Id3(3,2,2));
  vtkm::exec::ConnectivityStructured<
      vtkm::TopologyElementTagPoint, vtkm::TopologyElementTagCell, 3>
      connectivity3d(connectivityInternals3d);
  TryStructuredPointCoordinates(connectivity3d, coordinates);

  std::cout << "2D" << std::endl;
  vtkm::internal::ConnectivityStructuredInternals<2> connectivityInternals2d;
  connectivityInternals2d.SetPointDimensions(vtkm::Id2(3,2));
  vtkm::exec::ConnectivityStructured<
      vtkm::TopologyElementTagPoint, vtkm::TopologyElementTagCell, 2>
      connectivity2d(connectivityInternals2d);
  TryStructuredPointCoordinates(connectivity2d, coordinates);

  std::cout << "1D" << std::endl;
  vtkm::internal::ConnectivityStructuredInternals<1> connectivityInternals1d;
  connectivityInternals1d.SetPointDimensions(3);
  vtkm::exec::ConnectivityStructured<
      vtkm::TopologyElementTagPoint, vtkm::TopologyElementTagCell, 1>
      connectivity1d(connectivityInternals1d);
  TryStructuredPointCoordinates(connectivity1d, coordinates);
}

void TestArrayTopologyMapIn()
{
  vtkm::testing::Testing::TryTypes(TryType(),
                                   vtkm::TypeListTagCommon());

  TryStructuredPointCoordinates();
}

} // anonymous namespace

int UnitTestFetchArrayTopologyMapIn(int, char *[])
{
  return vtkm::testing::Testing::Run(TestArrayTopologyMapIn);
}
