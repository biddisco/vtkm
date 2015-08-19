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
#ifndef vtk_m_CellTraits_h
#define vtk_m_CellTraits_h

#include <vtkm/CellShape.h>

namespace vtkm {

/// \c vtkm::CellTraits::TopologyDimensionType is typedef to this with the
/// template parameter set to \c TOPOLOGICAL_DIMENSIONS. See \c
/// vtkm::CellTraits for more information.
///
template<vtkm::IdComponent dimension>
struct CellTopologicalDimensionsTag { };

/// \brief Information about a cell based on its tag.
///
/// The templated CellTraits struct provides the basic high level information
/// about cells (like the number of vertices in the cell or its
/// dimensionality).
///
template<class CellTag>
struct CellTraits
#ifdef VTKM_DOXYGEN_ONLY
{
  /// This defines the topological dimensions of the cell type. 3 for
  /// polyhedra, 2 for polygons, 1 for lines, 0 for points.
  ///
  const static vtkm::IdComponent TOPOLOGICAL_DIMENSIONS = 3;

  /// This tag is typedef'ed to
  /// vtkm::CellTopologicalDimensionsTag<TOPOLOGICAL_DIMENSIONS>. This provides
  /// a convenient way to overload a function based on topological dimensions
  /// (which is usually more efficient than conditionals).
  ///
  typedef vtkm::CellTopologicalDimensionsTag<TOPOLOGICAL_DIMENSIONS>
      TopologicalDimensionsTag;
};
#else // VTKM_DOXYGEN_ONLY
    ;
#endif // VTKM_DOXYGEN_ONLY

//-----------------------------------------------------------------------------

// Define traits for every cell type.

#define VTKM_DEFINE_CELL_TRAITS(name, dimensions) \
  template<> \
  struct CellTraits<vtkm::CellShapeTag ## name> { \
    const static vtkm::IdComponent TOPOLOGICAL_DIMENSIONS = dimensions; \
    typedef vtkm::CellTopologicalDimensionsTag<TOPOLOGICAL_DIMENSIONS> \
        TopologicalDimensionsTag; \
  }

VTKM_DEFINE_CELL_TRAITS(Empty, 0);
VTKM_DEFINE_CELL_TRAITS(Vertex, 0);
//VTKM_DEFINE_CELL_TRAITS(PolyVertex, 0);
VTKM_DEFINE_CELL_TRAITS(Line, 1);
//VTKM_DEFINE_CELL_TRAITS(PolyLine, 1);
VTKM_DEFINE_CELL_TRAITS(Triangle, 2);
//VTKM_DEFINE_CELL_TRAITS(TriangleStrip, 2);
VTKM_DEFINE_CELL_TRAITS(Polygon, 2);
VTKM_DEFINE_CELL_TRAITS(Pixel, 2);
VTKM_DEFINE_CELL_TRAITS(Quad, 2);
VTKM_DEFINE_CELL_TRAITS(Tetra, 3);
VTKM_DEFINE_CELL_TRAITS(Voxel, 3);
VTKM_DEFINE_CELL_TRAITS(Hexahedron, 3);
VTKM_DEFINE_CELL_TRAITS(Wedge, 3);
VTKM_DEFINE_CELL_TRAITS(Pyramid, 3);

#undef VTKM_DEFINE_CELL_TRAITS

} // namespace vtkm

#endif //vtk_m_CellTraits_h
