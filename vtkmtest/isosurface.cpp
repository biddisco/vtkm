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

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#ifndef VTKM_DEVICE_ADAPTER
# define VTKM_DEVICE_ADAPTER VTKM_DEVICE_ADAPTER_SERIAL
#endif


#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/DynamicArrayHandle.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/Pair.h>

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

#include "isosurface.h"

#if defined (__APPLE__)
# include <GLUT/glut.h>
#else
# include <GL/glut.h>
#endif
#include "quaternion.h"

int WinId = 0;

template <typename FieldType, typename OutputType>
class IsosurfaceFilterUniformGrid;

/// Linear interpolation
template <typename T1, typename T2>
VTKM_EXEC_EXPORT
T1 lerp(T1 a, T1 b, T2 t)
{
  return a + t*(b-a);
}

/// Vector cross-product
template <typename T>
VTKM_EXEC_EXPORT
T cross(T a, T b)
{
  return vtkm::make_Vec(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

/// Vector normalization
template <typename T>
VTKM_EXEC_EXPORT
T normalize(T v)
{
  return ((1.0f / sqrt(vtkm::dot(v,v))) * v);
}

template< typename FieldType>
struct PortalTypes
{
private:
  typedef vtkm::cont::ArrayHandle<FieldType> HandleType;
  typedef typename HandleType::template ExecutionTypes<DeviceAdapter> ExecutionTypes;
public:
  typedef typename ExecutionTypes::Portal Portal;
  typedef typename ExecutionTypes::PortalConst PortalConst;
};



/// Global variables
///
Quaternion qrot;
bool render_enabled = true;
int lastx, lasty;
int mouse_state = 1;
IsosurfaceFilterUniformGrid<vtkm::Float32, vtkm::Float32>* isosurfaceFilter;


/// \brief Computes Marching Cubes case number for each cell, along with the number of vertices generated by that case
///
template <typename FieldType>
class ClassifyCell : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldOut<AllTypes> outputCaseInfo);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  typedef typename PortalTypes<FieldType>::PortalConst FieldPortalType;
  FieldPortalType pointData;

  typedef typename PortalTypes<vtkm::Id>::PortalConst TablePortalType;
  TablePortalType vertexTable;

  float isovalue;
  int xdim, ydim, zdim;
  int cellsPerLayer;
  int pointsPerLayer;

  VTKM_CONT_EXPORT
  ClassifyCell(FieldPortalType pointData, TablePortalType vertexTable, float isovalue, int xdim, int ydim, int zdim) :
         pointData(pointData),
         vertexTable(vertexTable),
         isovalue(isovalue),
         xdim(xdim), ydim(ydim), zdim(zdim),
         cellsPerLayer((xdim - 1) * (ydim - 1)),
         pointsPerLayer (xdim*ydim) {}

  VTKM_EXEC_EXPORT
  vtkm::Pair<vtkm::Id, vtkm::Id> operator()(const vtkm::Id &cellId) const
  {
    // Compute 3D indices of this cell
    const int x = cellId % (xdim - 1);
    const int y = (cellId / (xdim - 1)) % (ydim -1);
    const int z = cellId / cellsPerLayer;

    // Compute indices for the eight vertices of this cell
    const int i0 = x    + y*xdim + z * pointsPerLayer;
    const int i1 = i0   + 1;
    const int i2 = i0   + 1 + xdim;
    const int i3 = i0   + xdim;
    const int i4 = i0   + pointsPerLayer;
    const int i5 = i1   + pointsPerLayer;
    const int i6 = i2   + pointsPerLayer;
    const int i7 = i3   + pointsPerLayer;

    // Get the field values at these eight vertices
    const float f0 = this->pointData.Get(i0);
    const float f1 = this->pointData.Get(i1);
    const float f2 = this->pointData.Get(i2);
    const float f3 = this->pointData.Get(i3);
    const float f4 = this->pointData.Get(i4);
    const float f5 = this->pointData.Get(i5);
    const float f6 = this->pointData.Get(i6);
    const float f7 = this->pointData.Get(i7);

    // Compute the Marching Cubes case number for this cell
    unsigned int cubeindex = (f0 > isovalue);
    cubeindex += (f1 > isovalue)*2;
    cubeindex += (f2 > isovalue)*4;
    cubeindex += (f3 > isovalue)*8;
    cubeindex += (f4 > isovalue)*16;
    cubeindex += (f5 > isovalue)*32;
    cubeindex += (f6 > isovalue)*64;
    cubeindex += (f7 > isovalue)*128;

    // Return the Marching Cubes case number and the number of vertices this case generates
    return vtkm::make_Pair(cubeindex, this->vertexTable.Get(cubeindex));
  }
};


/// \brief Return whether the cell generates geometry or not
///
class IsValidCell : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<AllTypes> inputCaseInfo, FieldOut<IdType> outputIsValid);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  VTKM_CONT_EXPORT
  IsValidCell() { };

  VTKM_EXEC_EXPORT
  vtkm::Id operator()(const vtkm::Pair<vtkm::Id, vtkm::Id> &caseInfo) const
  {
    return caseInfo.second != 0;
  }
};


/// \brief Return a value of the specified field of pairs at the given index, similar to a permutation iterator
///
class Permute : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputIndex, FieldOut<IdType> outputVerticesEnum);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  int pairIndex;

  typedef PortalTypes< vtkm::Pair< vtkm::Id, vtkm::Id> >::PortalConst InputDomainPortalType;
  InputDomainPortalType inputDomain;

  VTKM_CONT_EXPORT
  Permute(int pairIndex, InputDomainPortalType inputDomain) : pairIndex(pairIndex), inputDomain(inputDomain) {}

  VTKM_EXEC_EXPORT
  vtkm::Id operator()(const vtkm::Id &index) const
  {
    if (pairIndex == 1) return this->inputDomain.Get(index).first;
    return this->inputDomain.Get(index).second;
  }
};


/// \brief Compute isosurface vertices, normals, and scalars
///
template <typename FieldType, typename OutputType>
class IsosurfaceFunctorUniformGrid : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldOut<IdType> successFlag);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  typedef typename PortalTypes< vtkm::Id >::PortalConst IdPortalType;
  IdPortalType validCellIndices, outputVerticesEnum, triangleTable;

  typedef typename PortalTypes< vtkm::Pair<vtkm::Id, vtkm::Id> >::PortalConst IdPairPortalType;
  IdPairPortalType caseInfo;

  typedef typename PortalTypes< FieldType >::PortalConst FieldPortalType;
  FieldPortalType field, source;

  typedef typename PortalTypes< OutputType >::Portal ScalarPortalType;
  ScalarPortalType scalars;

  typedef typename PortalTypes< vtkm::Vec<OutputType,3> >::Portal VertexPortalType;
  VertexPortalType vertices, normals;

  const int xdim, ydim, zdim, cellsPerLayer;
  const float isovalue, xmin, ymin, zmin, xmax, ymax, zmax;

  VTKM_CONT_EXPORT
  IsosurfaceFunctorUniformGrid(const float isovalue, const int dims[3], const float mins[3], const float maxs[3], IdPortalType validCellIndices, IdPortalType outputVerticesEnum, IdPairPortalType caseInfo,
                   FieldPortalType field, FieldPortalType source, IdPortalType triangleTable, VertexPortalType vertices, VertexPortalType normals, ScalarPortalType scalars) : isovalue(isovalue), xdim(dims[0]),
                   ydim(dims[1]), zdim(dims[2]), xmin(mins[0]), ymin(mins[1]), zmin(mins[2]), xmax(maxs[0]), ymax(maxs[1]), zmax(maxs[2]), validCellIndices(validCellIndices), outputVerticesEnum(outputVerticesEnum),
                   caseInfo(caseInfo), field(field), source(source), triangleTable(triangleTable), vertices(vertices), normals(normals), scalars(scalars), cellsPerLayer((xdim-1) * (ydim-1)) {}

  VTKM_EXEC_EXPORT
  vtkm::Id operator()(const vtkm::Id &index) const
  {
    // Get data for this cell
    const vtkm::Id cellId = validCellIndices.Get(index);
    const int outputVertId = this->outputVerticesEnum.Get(index); //cellId);
    const int cubeindex    = caseInfo.Get(cellId).first;
    const int numVertices  = caseInfo.Get(cellId).second;
    const int verticesForEdge[] = { 0, 1, 1, 2, 3, 2, 0, 3,
                              4, 5, 5, 6, 7, 6, 4, 7,
                              0, 4, 1, 5, 2, 6, 3, 7 };

    // Compute 3D indices of this cell
    int x = cellId % (xdim - 1);
    int y = (cellId / (xdim - 1)) % (ydim - 1);
    int z = cellId / cellsPerLayer;

    // Compute indices for the eight vertices of this cell
    int i[8];
    i[0] = x      + y*xdim + z * xdim * ydim;
    i[1] = i[0]   + 1;
    i[2] = i[0]   + 1 + xdim;
    i[3] = i[0]   + xdim;
    i[4] = i[0]   + xdim * ydim;
    i[5] = i[1]   + xdim * ydim;
    i[6] = i[2]   + xdim * ydim;
    i[7] = i[3]   + xdim * ydim;

    // Get the field values at these eight vertices
    float f[8];
    f[0] = this->field.Get(i[0]);
    f[1] = this->field.Get(i[1]);
    f[2] = this->field.Get(i[2]);
    f[3] = this->field.Get(i[3]);
    f[4] = this->field.Get(i[4]);
    f[5] = this->field.Get(i[5]);
    f[6] = this->field.Get(i[6]);
    f[7] = this->field.Get(i[7]);

    // Compute the coordinates of the uniform regular grid at each of the cell's eight vertices
    vtkm::Vec<FieldType, 3> p[8];
    p[0] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[1] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[2] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[3] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[4] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
    p[5] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
    p[6] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
    p[7] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));

    // Get the scalar source values at the eight vertices
    float s[8];
    s[0] = this->source.Get(i[0]);
    s[1] = this->source.Get(i[1]);
    s[2] = this->source.Get(i[2]);
    s[3] = this->source.Get(i[3]);
    s[4] = this->source.Get(i[4]);
    s[5] = this->source.Get(i[5]);
    s[6] = this->source.Get(i[6]);
    s[7] = this->source.Get(i[7]);

    // Interpolate for vertex positions and associated scalar values
    for (int v = 0; v < numVertices; v++)
    {
      const int edge = this->triangleTable.Get(cubeindex*16 + v);
      const int v0   = verticesForEdge[2*edge];
      const int v1   = verticesForEdge[2*edge + 1];
      const float t  = (isovalue - f[v0]) / (f[v1] - f[v0]);

      this->vertices.Set(outputVertId + v, lerp(p[v0], p[v1], t));
      this->scalars.Set(outputVertId + v, lerp(s[v0], s[v1], t));
    }

    // Generate normal vectors by cross product of triangle edges
    for (int v = 0; v < numVertices; v += 3)
    {
      vtkm::Vec<OutputType, 3> vertex0 = this->vertices.Get(outputVertId + v + 0);
      vtkm::Vec<OutputType, 3> vertex1 = this->vertices.Get(outputVertId + v + 1);
      vtkm::Vec<OutputType, 3> vertex2 = this->vertices.Get(outputVertId + v + 2);

      vtkm::Vec<OutputType, 3> curNorm = normalize(cross(vertex1-vertex0, vertex2-vertex0));
      this->normals.Set(outputVertId + v + 0, curNorm);
      this->normals.Set(outputVertId + v + 1, curNorm);
      this->normals.Set(outputVertId + v + 2, curNorm);
    }

    return 0;
  }
};


/// \brief Return the tangle field value for each vertex
///
class TangleField : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> vertexId, FieldOut<Scalar> tangleFieldValue);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  const int xdim, ydim, zdim, cellsPerLayer;
  const float xmin, ymin, zmin, xmax, ymax, zmax;

  VTKM_CONT_EXPORT
  TangleField(const int dims[3], const float mins[3], const float maxs[3]) : xdim(dims[0]), ydim(dims[1]), zdim(dims[2]),
              xmin(mins[0]), ymin(mins[1]), zmin(mins[2]), xmax(maxs[0]), ymax(maxs[1]), zmax(maxs[2]), cellsPerLayer((xdim) * (ydim)) { };

  VTKM_EXEC_EXPORT
  vtkm::Float32 operator()(const vtkm::Id &vertexId) const
  {
    // Compute 3D indices of this cell
    const int x = vertexId % (xdim);
    const int y = (vertexId / (xdim)) % (ydim);
    const int z = vertexId / cellsPerLayer;

    const float xx = 3.0*(xmin+(xmax-xmin)*(1.0*x/(xdim-1)));
    const float yy = 3.0*(ymin+(ymax-ymin)*(1.0*y/(xdim-1)));
    const float zz = 3.0*(zmin+(zmax-zmin)*(1.0*z/(xdim-1)));

    const float v = (xx*xx*xx*xx - 5.0f*xx*xx + yy*yy*yy*yy - 5.0f*yy*yy + zz*zz*zz*zz - 5.0f*zz*zz + 11.8f) * 0.2f + 0.5f;

    return v;
  }
};


/// \brief Compute isosurface
///
template <typename FieldType, typename OutputType>
class IsosurfaceFilterUniformGrid
{
public:
  IsosurfaceFilterUniformGrid() { };

  vtkm::cont::ArrayHandle<vtkm::Vec<OutputType,3> > verticesArray, normalsArray;

  template<typename T>
  void OutputArrayDebug(T outputArray)
  {
    typedef typename T::ValueType ValueType;
    typedef typename T::PortalConstControl PortalConstType;
    PortalConstType readPortal = outputArray.GetPortalConstControl();
    vtkm::cont::ArrayPortalToIterators<PortalConstType> iterators(readPortal);
    std::vector<ValueType> result(readPortal.GetNumberOfValues());
    std::copy(iterators.GetBegin(), iterators.GetEnd(), result.begin());
    std::copy(result.begin(), result.end(), std::ostream_iterator<FieldType>(std::cout, " "));  std::cout << std::endl;
  }

  void computeIsosurface(int dim, float isovalue, char* fileName = 0)
  {
    // Initialize parameters; define min and max in x, y, and z for a uniform structured grid
    int vdim = dim + 1;  int dim3 = dim*dim*dim;
    int vdims[3] = { vdim, vdim, vdim };
    float mins[3] = {-1.0f, -1.0f, -1.0f};
    float maxs[3] = {1.0f, 1.0f, 1.0f};

    vtkm::cont::ArrayHandle<FieldType> fieldArray;
    if (fileName != 0)
    {
      // Read the field from a file
      std::vector<FieldType> field;
      std::fstream in(fileName, std::ios::in);
      std::copy(std::istream_iterator<FieldType>(in), std::istream_iterator<FieldType>(), std::back_inserter(field));
      std::copy(field.begin(), field.end(), std::ostream_iterator<FieldType>(std::cout, " "));  std::cout << std::endl;
      fieldArray = vtkm::cont::make_ArrayHandle(&field[0], field.size());
    }
    else
    {
      // Generate tangle field
      vtkm::cont::ArrayHandleCounting<vtkm::Id> vertexCountImplicitArray(0, vdim*vdim*vdim);
      vtkm::worklet::DispatcherMapField<TangleField> tangleFieldDispatcher(TangleField(vdims, mins, maxs));
      tangleFieldDispatcher.Invoke(vertexCountImplicitArray, fieldArray);
    }

    // Set up the Marching Cubes tables
    vtkm::cont::ArrayHandle<vtkm::Id> vertexTableArray = vtkm::cont::make_ArrayHandle(numVerticesTable, 256);
    vtkm::cont::ArrayHandle<vtkm::Id> triangleTableArray = vtkm::cont::make_ArrayHandle(triTable, 256*16);

    // Call the ClassifyCell functor to compute the Marching Cubes case numbers for each cell, and the number of vertices to be generated
    vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(0, dim3);
    vtkm::cont::ArrayHandle<vtkm::Pair<vtkm::Id, vtkm::Id> > caseInfoArray;
    vtkm::worklet::DispatcherMapField<ClassifyCell<FieldType> > classifyCellDispatcher(ClassifyCell<FieldType>(fieldArray.PrepareForInput(DeviceAdapter()),
                                                                                       vertexTableArray.PrepareForInput(DeviceAdapter()),
                                                                                       isovalue, vdim, vdim, vdim));
    classifyCellDispatcher.Invoke(cellCountImplicitArray, caseInfoArray);

    // Determine which cells are "valid" (i.e., produce geometry), and perform an inclusive scan to get running total of the number of "valid" cells
    vtkm::cont::ArrayHandle<vtkm::Id> validCellEnumArray, validCellIndicesArray, validVerticesArray, outputVerticesEnumArray, successArray;
    vtkm::worklet::DispatcherMapField<IsValidCell> isValidCellDispatcher;
    isValidCellDispatcher.Invoke(caseInfoArray, validCellEnumArray);
    unsigned int numValidCells = vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::ScanInclusive(validCellEnumArray, validCellEnumArray);

    // Return if no cells generate geometry
    std::cout << "NumValidCells: " << numValidCells << std::endl;
    if (numValidCells == 0) { render_enabled = false;  return; }

    // Use UpperBounds, a "permutation", and an exclusive scan to compute the starting output vertex index for each cell
    vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numValidCells);
    vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::UpperBounds(validCellEnumArray, validCellCountImplicitArray, validCellIndicesArray);
    vtkm::worklet::DispatcherMapField<Permute> permuteDispatcher(Permute(2, caseInfoArray.PrepareForInput(DeviceAdapter())));
    permuteDispatcher.Invoke(validCellIndicesArray, validVerticesArray);
    int numTotalVertices = vtkm::cont::DeviceAdapterAlgorithm<VTKM_DEFAULT_DEVICE_ADAPTER_TAG>::ScanExclusive(validVerticesArray, outputVerticesEnumArray);

    // Call the IsosurfaceFunctor to actually compute all the output vertices, normals, and scalars
    vtkm::cont::ArrayHandle<OutputType> scalarsArray;
    vtkm::worklet::DispatcherMapField<IsosurfaceFunctorUniformGrid<FieldType, OutputType> > isosurfaceFunctorDispatcher(IsosurfaceFunctorUniformGrid<FieldType, OutputType>(isovalue, vdims, mins, maxs,
                                                                                                       validCellIndicesArray.PrepareForInput(DeviceAdapter()),
                                                                                                       outputVerticesEnumArray.PrepareForInput(DeviceAdapter()),
                                                                                                       caseInfoArray.PrepareForInput(DeviceAdapter()),
                                                                                                       fieldArray.PrepareForInput(DeviceAdapter()),
                                                                                                       fieldArray.PrepareForInput(DeviceAdapter()),
                                                                                                       triangleTableArray.PrepareForInput(DeviceAdapter()),
                                                                                                       verticesArray.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                                                                                                       normalsArray.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                                                                                                       scalarsArray.PrepareForOutput(numTotalVertices, DeviceAdapter())));
    isosurfaceFunctorDispatcher.Invoke(validCellCountImplicitArray, successArray);
  }
};


/// Initialize OpenGL parameters, including lighting
///
void initializeGL()
{
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);

  float white[] = { 0.8, 0.8, 0.8, 1.0 };
  float black[] = { 0.0, 0.0, 0.0, 1.0 };
  float lightPos[] = { 10.0, 10.0, 10.5, 1.0 };

  glLightfv(GL_LIGHT0, GL_AMBIENT, white);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
  glLightfv(GL_LIGHT0, GL_SPECULAR, black);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
}


/// Render the computed triangles
///
void displayCall()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective( 45.0f, 1.0f, 1.0f, 20.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

  glPushMatrix();
  float rotationMatrix[16];
  qrot.getRotMat(rotationMatrix);
  glMultMatrixf(rotationMatrix);

  glColor3f(0.1f, 0.1f, 0.6f);

  if (render_enabled)
  {
    glBegin(GL_TRIANGLES);
    for (unsigned int i=0; i<isosurfaceFilter->verticesArray.GetPortalConstControl().GetNumberOfValues(); i++)
    {
      vtkm::Vec<vtkm::Float32, 3> curNormal = isosurfaceFilter->normalsArray.GetPortalConstControl().Get(i);
      vtkm::Vec<vtkm::Float32, 3> curVertex = isosurfaceFilter->verticesArray.GetPortalConstControl().Get(i);
      glNormal3f(curNormal[0], curNormal[1], curNormal[2]);
      glVertex3f(curVertex[0], curVertex[1], curVertex[2]);
    }
    glEnd();
  }

  glPopMatrix();
  glutSwapBuffers();
}


/// Handle mouse button pushes
///
void mouseCall(int button, int state, int x, int y)
{
  if (button == 0) mouse_state = state;
  if ((button == 0) && (state == 0)) { lastx = x;  lasty = y; }
}


/// Handle mouse movements to rotate the object using quaternions
///
void mouseMove(int x, int y)
{
  int dx = x - lastx;
  int dy = y - lasty;

  if (mouse_state == 0)
  {
    Quaternion newRotX;
    newRotX.setEulerAngles(-0.2*dx*M_PI/180.0, 0.0, 0.0);
    qrot.mul(newRotX);

    Quaternion newRotY;
    newRotY.setEulerAngles(0.0, 0.0, -0.2*dy*M_PI/180.0);
    qrot.mul(newRotY);
  }
  lastx = x;
  lasty = y;

  glutPostRedisplay();
}

void keyboardCB( unsigned char key, int x, int y )
{
  switch ( key )
  {
    case 27: // Escape key
      glutDestroyWindow ( WinId );
      exit (0);
      break;
  }
  glutPostRedisplay();
}

/// Print a vector
///
std::string vec3String(const vtkm::Vec<vtkm::Float32,3>& data)
{
  std::ostringstream str;
  str << "(" << data[0] << ", " << data[1] << ", " << data[2] << ")";
  return str.str();
}


/// Compute an isosurface, and render it using GLUT
///
int main(int argc, char* argv[])
{
  // Abort if dimension and file name are not provided
  if (argc < 3)
  {
    std::cout << "Usage: isosurface {dimension} {isovalue} {optional-file-name} " << std::endl;
    return 0;
  }

  // Create the filter and compute the isosurface
  isosurfaceFilter = new IsosurfaceFilterUniformGrid<vtkm::Float32, vtkm::Float32>();
  if (argc == 3) isosurfaceFilter->computeIsosurface(atoi(argv[1]), atof(argv[2]));
  if (argc >= 4) isosurfaceFilter->computeIsosurface(atoi(argv[1]), atof(argv[2]), argv[3]);

  // Print the vertices
  // vtkm::cont::ArrayPortalToIterators<vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3> >::PortalConstControl> testIterators(isosurfaceFilter->verticesArray.GetPortalConstControl());
  // std::transform(testIterators.GetBegin(), testIterators.GetEnd(), std::ostream_iterator<std::string>(std::cout, " "), vec3String);
  // std::cout << std::endl;

  // Launch GLUT to render the resulting isosurface
  lastx = lasty = 0;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(1000, 1000);
  glutInitWindowPosition(300, 200);
  WinId = glutCreateWindow("VTK-m Isosurface");
  initializeGL();
  glutDisplayFunc(displayCall);
  glutMotionFunc(mouseMove);
  glutMouseFunc(mouseCall);
  glutKeyboardFunc(keyboardCB);
  glutMainLoop();

  return 0;
}


