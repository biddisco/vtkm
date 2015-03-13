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
#ifndef vtk_hpx_h
#define vtk_hpx_h

// squash windows #defines before #including hpx
#define NOMINMAX

#ifdef _WIN32
 #define BOOST_PROGRAM_OPTIONS_DYN_LINK
#endif

// override int main() so that hpx is initialized on startup
#include <hpx/hpx_main.hpp>  

// squash windows #defines FTER #INCLUDING HPX
#undef GetMessage

#endif //vtk_hpx_h
