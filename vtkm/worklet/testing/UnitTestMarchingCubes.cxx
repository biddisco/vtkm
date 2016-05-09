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

#include <vtkm/Math.h>
#include <vtkm/cont/ArrayHandleUniformPointCoordinates.h>
#include <vtkm/cont/CellSetSingleType.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DynamicArrayHandle.h>
#include <vtkm/cont/testing/Testing.h>

#include <vtkm/worklet/MarchingCubes.h>
#include <vtkm/worklet/DispatcherMapField.h>

namespace {

class TangleField : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> vertexId, FieldOut<Scalar> v);
  typedef void ExecutionSignature(_1, _2);
  typedef _1 InputDomain;

  const vtkm::Id xdim, ydim, zdim;
  const vtkm::FloatDefault xmin, ymin, zmin, xmax, ymax, zmax;
  const vtkm::Id cellsPerLayer;

  VTKM_CONT_EXPORT
  TangleField(const vtkm::Id3 dims, const vtkm::FloatDefault mins[3], const vtkm::FloatDefault maxs[3]) : xdim(dims[0]), ydim(dims[1]), zdim(dims[2]),
              xmin(mins[0]), ymin(mins[1]), zmin(mins[2]), xmax(maxs[0]), ymax(maxs[1]), zmax(maxs[2]), cellsPerLayer((xdim) * (ydim)) { }

  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id &vertexId, vtkm::Float32 &v) const
  {
    const vtkm::Id x = vertexId % (xdim);
    const vtkm::Id y = (vertexId / (xdim)) % (ydim);
    const vtkm::Id z = vertexId / cellsPerLayer;

    const vtkm::FloatDefault fx = static_cast<vtkm::FloatDefault>(x) / static_cast<vtkm::FloatDefault>(xdim-1);
    const vtkm::FloatDefault fy = static_cast<vtkm::FloatDefault>(y) / static_cast<vtkm::FloatDefault>(xdim-1);
    const vtkm::FloatDefault fz = static_cast<vtkm::FloatDefault>(z) / static_cast<vtkm::FloatDefault>(xdim-1);

    const vtkm::Float32 xx = 3.0f*vtkm::Float32(xmin+(xmax-xmin)*(fx));
    const vtkm::Float32 yy = 3.0f*vtkm::Float32(ymin+(ymax-ymin)*(fy));
    const vtkm::Float32 zz = 3.0f*vtkm::Float32(zmin+(zmax-zmin)*(fz));

    v = (xx*xx*xx*xx - 5.0f*xx*xx + yy*yy*yy*yy - 5.0f*yy*yy + zz*zz*zz*zz - 5.0f*zz*zz + 11.8f) * 0.2f + 0.5f;
  }
};

vtkm::cont::DataSet MakeIsosurfaceTestDataSet(vtkm::Id3 dims)
{
  vtkm::cont::DataSet dataSet;

  const vtkm::Id3 vdims(dims[0] + 1, dims[1] + 1, dims[2] + 1);

  vtkm::FloatDefault mins[3] = {-1.0f, -1.0f, -1.0f};
  vtkm::FloatDefault maxs[3] = {1.0f, 1.0f, 1.0f};

  vtkm::cont::ArrayHandle<vtkm::Float32> fieldArray;
  vtkm::cont::ArrayHandleIndex vertexCountImplicitArray(vdims[0]*vdims[1]*vdims[2]);
  vtkm::worklet::DispatcherMapField<TangleField> tangleFieldDispatcher(TangleField(vdims, mins, maxs));
  tangleFieldDispatcher.Invoke(vertexCountImplicitArray, fieldArray);

  vtkm::Vec<vtkm::FloatDefault,3> origin(0.0f, 0.0f, 0.0f);
  vtkm::Vec<vtkm::FloatDefault,3> spacing(
        1.0f/static_cast<vtkm::FloatDefault>(dims[0]),
        1.0f/static_cast<vtkm::FloatDefault>(dims[2]),
        1.0f/static_cast<vtkm::FloatDefault>(dims[1]));

  vtkm::cont::ArrayHandleUniformPointCoordinates
      coordinates(vdims, origin, spacing);
  dataSet.AddCoordinateSystem(
          vtkm::cont::CoordinateSystem("coordinates", coordinates));

  dataSet.AddField(vtkm::cont::Field("nodevar", vtkm::cont::Field::ASSOC_POINTS, fieldArray));

  static const vtkm::IdComponent ndim = 3;
  vtkm::cont::CellSetStructured<ndim> cellSet("cells");
  cellSet.SetPointDimensions(vdims);
  dataSet.AddCellSet(cellSet);

  return dataSet;
}

class EuclideanNorm
{
public:
  VTKM_EXEC_CONT_EXPORT
  EuclideanNorm() : Reference(0.,0.,0.) {}
  VTKM_EXEC_CONT_EXPORT
  EuclideanNorm(vtkm::Vec<vtkm::Float32,3> reference):Reference(reference) {}

  VTKM_EXEC_CONT_EXPORT
  vtkm::Float32 operator()(vtkm::Vec<vtkm::Float32,3> v) const
  {
    vtkm::Vec<vtkm::Float32,3> d(v[0]-this->Reference[0],
                                 v[1]-this->Reference[1],
                                 v[2]-this->Reference[2]);
    return vtkm::Magnitude(d);
  }

private:
  vtkm::Vec<vtkm::Float32,3> Reference;
};

class CubeGridConnectivity
{
public:
  VTKM_EXEC_CONT_EXPORT
  CubeGridConnectivity() : Dimension(1),
                           DimSquared(1),
                           DimPlus1Squared(4) {}
  VTKM_EXEC_CONT_EXPORT
  CubeGridConnectivity(vtkm::Id dim) : Dimension(dim),
                                       DimSquared(dim*dim),
                                       DimPlus1Squared((dim+1)*(dim+1)) {}

  VTKM_EXEC_CONT_EXPORT
  vtkm::Id operator()(vtkm::Id vertex) const
  {
    typedef vtkm::CellShapeTagHexahedron HexTag;
    typedef vtkm::CellTraits<HexTag> HexTraits;

    vtkm::Id cellId = vertex/HexTraits::NUM_POINTS;
    vtkm::Id localId = vertex%HexTraits::NUM_POINTS;
    vtkm::Id globalId =
      (cellId + cellId/this->Dimension +
       (this->Dimension+1)*(cellId/(this->DimSquared)));

    switch (localId)
      {
      case 2: globalId += 1;
      case 3: globalId += this->Dimension;
      case 1: globalId += 1;
      case 0: break;
      case 6: globalId += 1;
      case 7: globalId += this->Dimension;
      case 5: globalId += 1;
      case 4: globalId += this->DimPlus1Squared; break;
      }

    return globalId;
  }

private:
  vtkm::Id Dimension;
  vtkm::Id DimSquared;
  vtkm::Id DimPlus1Squared;
};

class MakeRadiantDataSet
{
public:
  typedef vtkm::cont::ArrayHandleUniformPointCoordinates CoordinateArrayHandle;
  typedef vtkm::cont::ArrayHandleTransform<vtkm::Float32,
    vtkm::cont::ArrayHandleUniformPointCoordinates,
  EuclideanNorm> DataArrayHandle;
  typedef vtkm::cont::ArrayHandleTransform<vtkm::Id,
    vtkm::cont::ArrayHandleCounting<vtkm::Id>,
    CubeGridConnectivity> ConnectivityArrayHandle;

  typedef vtkm::cont::CellSetSingleType<
    vtkm::cont::ArrayHandleTransform<vtkm::Id,
      vtkm::cont::ArrayHandleCounting<vtkm::Id>,
      CubeGridConnectivity>::StorageTag> CellSet;

  vtkm::cont::DataSet Make3DRadiantDataSet(vtkm::IdComponent dim=5);
};

inline vtkm::cont::DataSet MakeRadiantDataSet::Make3DRadiantDataSet(vtkm::IdComponent dim)
{
  // create a cube from -.5 to .5 in x,y,z, consisting of <dim> cells on each
  // axis, with point values equal to the Euclidean distance from the origin.

  vtkm::cont::DataSet dataSet;

  typedef vtkm::CellShapeTagHexahedron HexTag;
  typedef vtkm::CellTraits<HexTag> HexTraits;

  typedef vtkm::Vec<vtkm::Float32,3> CoordType;

  const vtkm::IdComponent nCells = dim*dim*dim;

vtkm::Float32 spacing = vtkm::Float32(1./dim);
  CoordinateArrayHandle coordinates(vtkm::Id3(dim+1,dim+1,dim+1),
                                    CoordType(-.5,-.5,-.5),
                                    CoordType(spacing,spacing,spacing));

  DataArrayHandle distanceToOrigin(coordinates);
  DataArrayHandle distanceToOther(coordinates,
                                  EuclideanNorm(CoordType(1.,1.,1.)));

  ConnectivityArrayHandle connectivity(
    vtkm::cont::ArrayHandleCounting<vtkm::Id>(0,1,nCells*HexTraits::NUM_POINTS),
    CubeGridConnectivity(dim));

  dataSet.AddCoordinateSystem(
        vtkm::cont::CoordinateSystem("coordinates", coordinates));

  //Set point scalar
  dataSet.AddField(
    vtkm::cont::Field("distanceToOrigin", vtkm::cont::Field::ASSOC_POINTS,
                      vtkm::cont::DynamicArrayHandle(distanceToOrigin)));
  dataSet.AddField(
    vtkm::cont::Field("distanceToOther", vtkm::cont::Field::ASSOC_POINTS,
                      vtkm::cont::DynamicArrayHandle(distanceToOther)));

  CellSet cellSet(HexTag(), "cells");
  cellSet.Fill(connectivity);

  dataSet.AddCellSet(cellSet);

  return dataSet;
}

} // anonymous namespace

void TestMarchingCubesUniformGrid()
{
  std::cout << "Testing MarchingCubes filter on a uniform grid" << std::endl;

  vtkm::Id3 dims(4,4,4);
  vtkm::cont::DataSet dataSet = MakeIsosurfaceTestDataSet(dims);

  typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;
  vtkm::cont::CellSetStructured<3> cellSet;
  dataSet.GetCellSet().CopyTo(cellSet);
  vtkm::cont::ArrayHandle<vtkm::Float32> fieldArray;
  dataSet.GetField("nodevar").GetData().CopyTo(fieldArray);

  vtkm::worklet::MarchingCubes<vtkm::Float32,DeviceAdapter> isosurfaceFilter;

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3> > verticesArray;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3> > normalsArray;
  vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray;
  isosurfaceFilter.Run(0.5,
                       cellSet,
                       dataSet.GetCoordinateSystem(),
                       fieldArray,
                       verticesArray,
                       normalsArray);

  isosurfaceFilter.MapFieldOntoIsosurface(fieldArray,
                                          scalarsArray);

  std::cout << "vertices: ";
  vtkm::cont::printSummary_ArrayHandle(verticesArray, std::cout);
  std::cout << std::endl;
  std::cout << "normals: ";
  vtkm::cont::printSummary_ArrayHandle(normalsArray, std::cout);
  std::cout << std::endl;
  std::cout << "scalars: ";
  vtkm::cont::printSummary_ArrayHandle(scalarsArray, std::cout);
  std::cout << std::endl;

  VTKM_TEST_ASSERT(test_equal(verticesArray.GetNumberOfValues(), 480),
                   "Wrong result for Isosurface filter");
}

void TestMarchingCubesExplicit()
{
  std::cout << "Testing MarchingCubes filter on explicit data" << std::endl;

  typedef MakeRadiantDataSet DataSetGenerator;
  typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceTag;
  typedef vtkm::worklet::MarchingCubes<vtkm::Float32,DeviceTag> MarchingCubes;
  typedef vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3> > Vec3Handle;
  typedef vtkm::cont::ArrayHandle<vtkm::Float32> DataHandle;

  DataSetGenerator dataSetGenerator;

  vtkm::IdComponent Dimension = 10;
  vtkm::Float32 contourValue = vtkm::Float32(.45);

  vtkm::cont::DataSet dataSet =
    dataSetGenerator.Make3DRadiantDataSet(Dimension);

  DataSetGenerator::CellSet cellSet;
  dataSet.GetCellSet().CopyTo(cellSet);

  vtkm::cont::Field contourField = dataSet.GetField("distanceToOrigin");
  DataSetGenerator::DataArrayHandle contourArray;
  contourField.GetData().CopyTo(contourArray);
  Vec3Handle vertices;
  Vec3Handle normals;

  MarchingCubes marchingCubes;

  marchingCubes.Run(contourValue,
                    cellSet,
                    dataSet.GetCoordinateSystem(),
                    contourArray,
                    vertices,
                    normals);

  DataHandle scalars;

  vtkm::cont::Field projectedField = dataSet.GetField("distanceToOther");

  DataSetGenerator::DataArrayHandle projectedArray;
  projectedField.GetData().CopyTo(projectedArray);

  marchingCubes.MapFieldOntoIsosurface(projectedArray,
                                       scalars);

  std::cout << "vertices: ";
  vtkm::cont::printSummary_ArrayHandle(vertices, std::cout);
  std::cout << std::endl;
  std::cout << "normals: ";
  vtkm::cont::printSummary_ArrayHandle(normals, std::cout);
  std::cout << std::endl;
  std::cout << "scalars: ";
  vtkm::cont::printSummary_ArrayHandle(scalars, std::cout);
  std::cout << std::endl;

  VTKM_TEST_ASSERT(test_equal(vertices.GetNumberOfValues(), 2472),
                   "Wrong vertices result for MarchingCubes filter");
  VTKM_TEST_ASSERT(test_equal(normals.GetNumberOfValues(), 2472),
                   "Wrong normals result for MarchingCubes filter");
  VTKM_TEST_ASSERT(test_equal(scalars.GetNumberOfValues(), 2472),
                   "Wrong scalars result for MarchingCubes filter");
}

int UnitTestMarchingCubes(int, char *[])
{
  return vtkm::cont::testing::Testing::Run(TestMarchingCubesUniformGrid);
  return vtkm::cont::testing::Testing::Run(TestMarchingCubesExplicit);
}