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
#ifndef vtk_m_cont_arg_ControlSignatureTagBase_h
#define vtk_m_cont_arg_ControlSignatureTagBase_h

#include <vtkm/internal/ExportMacros.h>

VTKM_THIRDPARTY_PRE_INCLUDE
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_of.hpp>
VTKM_THIRDPARTY_POST_INCLUDE

namespace vtkm {
namespace cont {
namespace arg {

/// \brief The base class for all tags used in a \c ControlSignature.
///
/// If a new \c ControlSignature tag is created, it must be derived from this
/// class in some way. This helps identify \c ControlSignature tags in the \c
/// VTKM_IS_CONTROL_SIGNATURE_TAG macro and allows checking the validity of a
/// \c ControlSignature.
///
/// In addition to inheriting from this base class, a \c ControlSignature tag
/// must define the following three typedefs: \c TypeCheckTag, \c TransportTag
/// and \c FetchTag.
///
struct ControlSignatureTagBase {  };

namespace internal {

template<typename ControlSignatureTag>
struct ControlSignatureTagCheck
{
  static const bool Valid =
      boost::is_base_of<
          vtkm::cont::arg::ControlSignatureTagBase, ControlSignatureTag>::value;
};

} // namespace internal

/// Checks that the argument is a proper tag for an \c ControlSignature. This
/// is a handy concept check when modifying tags or dispatching to make sure
/// that a template argument is actually an \c ControlSignature tag. (You can
/// get weird errors elsewhere in the code when a mistake is made.)
///
#define VTKM_IS_CONTROL_SIGNATURE_TAG(tag) \
  BOOST_STATIC_ASSERT_MSG( \
    ::vtkm::cont::arg::internal::ControlSignatureTagCheck<tag>::Valid, \
    "Provided a type that is not a valid ControlSignature tag.")

}
}
} // namespace vtkm::cont::arg

#endif //vtk_m_cont_arg_ControlSignatureTagBase_h
