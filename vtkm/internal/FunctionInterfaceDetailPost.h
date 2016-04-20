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
// **** DO NOT EDIT THIS FILE!!! ****
// This file is automatically generated by FunctionInterfaceDetailPost.h.in

#ifndef vtk_m_internal_FunctionInterfaceDetailPost_h
#define vtk_m_internal_FunctionInterfaceDetailPost_h

#if !defined(vtk_m_internal_FunctionInterface_h) && !defined(VTKM_TEST_HEADER_BUILD)
#error FunctionInterfaceDetailPre.h must be included from FunctionInterface.h
#endif

#include <vtkm/internal/FunctionInterface.h>

#if VTKM_MAX_FUNCTION_PARAMETERS != 10
#error Mismatch of maximum parameters between FunctionInterfaceDatailPre.h.in and FunctionInterfaceDetailPost.h.in
#endif


namespace vtkm {
namespace internal {

namespace detail {

//============================================================================

template<typename Transform,
         typename R>
struct FunctionInterfaceStaticTransformType<R(), Transform> {
  typedef R(type)(
        );
};

template<typename Transform,
         typename R,
         typename P1>
struct FunctionInterfaceStaticTransformType<R(P1), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2>
struct FunctionInterfaceStaticTransformType<R(P1,P2), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4,P5), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type,
        typename Transform::template ReturnType<P5,5>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4,P5,P6), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type,
        typename Transform::template ReturnType<P5,5>::type,
        typename Transform::template ReturnType<P6,6>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4,P5,P6,P7), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type,
        typename Transform::template ReturnType<P5,5>::type,
        typename Transform::template ReturnType<P6,6>::type,
        typename Transform::template ReturnType<P7,7>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7,
         typename P8>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4,P5,P6,P7,P8), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type,
        typename Transform::template ReturnType<P5,5>::type,
        typename Transform::template ReturnType<P6,6>::type,
        typename Transform::template ReturnType<P7,7>::type,
        typename Transform::template ReturnType<P8,8>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7,
         typename P8,
         typename P9>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4,P5,P6,P7,P8,P9), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type,
        typename Transform::template ReturnType<P5,5>::type,
        typename Transform::template ReturnType<P6,6>::type,
        typename Transform::template ReturnType<P7,7>::type,
        typename Transform::template ReturnType<P8,8>::type,
        typename Transform::template ReturnType<P9,9>::type
        );
};

template<typename Transform,
         typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7,
         typename P8,
         typename P9,
         typename P10>
struct FunctionInterfaceStaticTransformType<R(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10), Transform> {
  typedef R(type)(
        typename Transform::template ReturnType<P1,1>::type,
        typename Transform::template ReturnType<P2,2>::type,
        typename Transform::template ReturnType<P3,3>::type,
        typename Transform::template ReturnType<P4,4>::type,
        typename Transform::template ReturnType<P5,5>::type,
        typename Transform::template ReturnType<P6,6>::type,
        typename Transform::template ReturnType<P7,7>::type,
        typename Transform::template ReturnType<P8,8>::type,
        typename Transform::template ReturnType<P9,9>::type,
        typename Transform::template ReturnType<P10,10>::type
        );
};


} // namespace detail

//============================================================================

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R()>
make_FunctionInterface(
      )
{
  FunctionInterface<R()> fi;
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1)>
make_FunctionInterface(
      const P1& p1
      )
{
  FunctionInterface<R(P1)> fi;
  fi.template SetParameter<1>(p1);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2
      )
{
  FunctionInterface<R(P1,P2)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3
      )
{
  FunctionInterface<R(P1,P2,P3)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4
      )
{
  FunctionInterface<R(P1,P2,P3,P4)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4,P5)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4,
      const P5& p5
      )
{
  FunctionInterface<R(P1,P2,P3,P4,P5)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  fi.template SetParameter<5>(p5);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4,P5,P6)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4,
      const P5& p5,
      const P6& p6
      )
{
  FunctionInterface<R(P1,P2,P3,P4,P5,P6)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  fi.template SetParameter<5>(p5);
  fi.template SetParameter<6>(p6);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4,
      const P5& p5,
      const P6& p6,
      const P7& p7
      )
{
  FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  fi.template SetParameter<5>(p5);
  fi.template SetParameter<6>(p6);
  fi.template SetParameter<7>(p7);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7,
         typename P8>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7,P8)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4,
      const P5& p5,
      const P6& p6,
      const P7& p7,
      const P8& p8
      )
{
  FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7,P8)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  fi.template SetParameter<5>(p5);
  fi.template SetParameter<6>(p6);
  fi.template SetParameter<7>(p7);
  fi.template SetParameter<8>(p8);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7,
         typename P8,
         typename P9>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7,P8,P9)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4,
      const P5& p5,
      const P6& p6,
      const P7& p7,
      const P8& p8,
      const P9& p9
      )
{
  FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7,P8,P9)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  fi.template SetParameter<5>(p5);
  fi.template SetParameter<6>(p6);
  fi.template SetParameter<7>(p7);
  fi.template SetParameter<8>(p8);
  fi.template SetParameter<9>(p9);
  return fi;
}

/// \brief Create a \c FunctionInterface
///
/// \c make_FunctionInterface is a function that takes a variable number of
/// arguments and returns a \c FunctionInterface object containing these
/// objects. Since the return type for the function signature is not specified,
/// you must always specify it as a template parameter
///
/// \code{.cpp}
/// vtkm::internal::FunctionInterface<void(int,double,char)> functionInterface =
///     vtkm::internal::make_FunctionInterface<void>(1, 2.5, 'a');
/// \endcode
///
VTKM_SUPPRESS_EXEC_WARNINGS
template<typename R,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5,
         typename P6,
         typename P7,
         typename P8,
         typename P9,
         typename P10>
VTKM_EXEC_CONT_EXPORT
FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10)>
make_FunctionInterface(
      const P1& p1,
      const P2& p2,
      const P3& p3,
      const P4& p4,
      const P5& p5,
      const P6& p6,
      const P7& p7,
      const P8& p8,
      const P9& p9,
      const P10& p10
      )
{
  FunctionInterface<R(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10)> fi;
  fi.template SetParameter<1>(p1);
  fi.template SetParameter<2>(p2);
  fi.template SetParameter<3>(p3);
  fi.template SetParameter<4>(p4);
  fi.template SetParameter<5>(p5);
  fi.template SetParameter<6>(p6);
  fi.template SetParameter<7>(p7);
  fi.template SetParameter<8>(p8);
  fi.template SetParameter<9>(p9);
  fi.template SetParameter<10>(p10);
  return fi;
}


}
} // namespace vtkm::internal

#endif //vtk_m_internal_FunctionInterfaceDetailPost_h
