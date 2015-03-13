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

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DynamicArrayHandle.h>

#include <vtkm/cont/testing/Testing.h>

void TestDataSet_Explicit1()
{
  std::cout << "1" << std::endl;
    std::vector<vtkm::cont::ArrayHandle<vtkm::FloatDefault, vtkm::cont::StorageTagBasic> > Fields;
    std::cout << "2" << std::endl;
    Fields.resize(5);
    std::cout << "3" << std::endl;

    const int nVerts = 5;
    vtkm::Float32 xVals[nVerts] = {10, 11, 12, 13, 14};
    vtkm::Float32 yVals[nVerts] = {20, 21, 22, 23, 24};

//    vtkm::cont::ArrayHandle<vtkm::FloatDefault,vtkm::cont::StorageTagBasic> field1;
//    vtkm::cont::ArrayHandle<vtkm::FloatDefault,vtkm::cont::StorageTagBasic> field2;
//    Fields[0] = field1;
//    Fields[1] = field2;

    vtkm::cont::ArrayHandle<vtkm::Float32> tmp1 = vtkm::cont::make_ArrayHandle(xVals, nVerts);
    vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::Copy(tmp1, Fields[0]);

    vtkm::cont::ArrayHandle<vtkm::Float32> tmp2 = vtkm::cont::make_ArrayHandle(yVals, nVerts);
    vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::Copy(tmp2, Fields[1]);

    for (int i=0; i<5; i++) {
      std::cout << "first val of field 1 = " << Fields[0].GetPortalControl().Get(i)<<std::endl;
      std::cout << "first val of field 2 = " << Fields[1].GetPortalControl().Get(i)<<std::endl;
    }
}

void TestDataSet_Explicit2()
{
    const int nVerts = 5;
    vtkm::Float32 xVals[nVerts] = {10, 11, 12, 13, 14};
    vtkm::Float32 yVals[nVerts] = {20, 21, 22, 23, 15};

    vtkm::cont::ArrayHandle<vtkm::FloatDefault,vtkm::cont::StorageTagBasic> field1;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault,vtkm::cont::StorageTagBasic> field2;

    vtkm::cont::ArrayHandle<vtkm::Float32> tmp1 = vtkm::cont::make_ArrayHandle(xVals, nVerts);
    vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::Copy(tmp1, field1);

    vtkm::cont::ArrayHandle<vtkm::Float32> tmp2 = vtkm::cont::make_ArrayHandle(yVals, nVerts);
    vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::Copy(tmp2, field2);

    std::cout << "first val of field 1 = " << field1.GetPortalControl().Get(0)<<std::endl;
    std::cout << "first val of field 2 = " << field2.GetPortalControl().Get(0)<<std::endl;
}

int main() {
  TestDataSet_Explicit1();
  TestDataSet_Explicit2();
  return 0;
}

