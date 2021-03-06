##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##
##  Copyright 2014 Sandia Corporation.
##  Copyright 2014 UT-Battelle, LLC.
##  Copyright 2014 Los Alamos National Security.
##
##  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
##  the U.S. Government retains certain rights in this software.
##
##  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
##  Laboratory (LANL), the U.S. Government retains certain rights in
##  this software.
##============================================================================

# This script is used to create the SystemInformation test. The test always
# passes. It just captures in its output the configuration of the system.
# This allows you to inspect the configuration of the system of a failed
# dashboard in case you don't have access to that dashboard.
#
# This script is called with a command like:
#
# cmake -D VTKm_BINARY_DIR=<top-of-build-tree> -D VTKm_SOURCE_DIR=<top-of-source-tree> -P <this-script>
#

set(FILES
  vtkm/internal/Configure.h
  CMakeCache.txt
  CMakeFiles/CMakeError.log
  )

function(print_file filename)
  set(full_filename "${VTKm_BINARY_DIR}/${filename}")
  message("

==============================================================================

Contents of \"${filename}\":
------------------------------------------------------------------------------")
  if(EXISTS "${full_filename}")
    file(READ ${full_filename} contents)
    message("${contents}")
  else()
    message("The file \"${full_filename}\" does not exist.")
  endif()
endfunction(print_file)


message("CTEST_FULL_OUTPUT (Avoid ctest truncation of output)")

execute_process(
  COMMAND git rev-parse -q HEAD
  WORKING_DIRECTORY "${VTKm_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA
  )

message("

==============================================================================

git SHA: ${GIT_SHA}")

foreach(filename ${FILES})
  print_file(${filename})
endforeach()
