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
#ifndef vtkm_exec_internal_ErrorMessageBuffer_h
#define vtkm_exec_internal_ErrorMessageBuffer_h

#include <vtkm/Types.h>

namespace vtkm {
namespace exec {
namespace internal {

/// Used to hold an error in the execution environment until the parallel
/// execution can complete. This is to be used in conjunction with a
/// DeviceAdapter's Schedule function to implement errors in execution
/// environments that cannot throw errors. This string should be global to all
/// threads. If the first entry in the string is '\0' (the C string
/// terminator), then we consider it as no error. Otherwise, the array contains
/// the string describing the error.
///
/// Before scheduling worklets, the global array should be cleared to have no
/// error. This can only be reliably done by the device adapter.
///
class ErrorMessageBuffer
{
public:
  VTKM_EXEC_CONT_EXPORT ErrorMessageBuffer()
    : MessageBuffer(), MessageBufferSize(0) {  }

  VTKM_EXEC_CONT_EXPORT
  ErrorMessageBuffer(char *messageBuffer, vtkm::Id bufferSize)
    : MessageBuffer(messageBuffer), MessageBufferSize(bufferSize) { }

  VTKM_EXEC_EXPORT void RaiseError(const char *message) const
  {
    // Only raise the error if one has not been raised yet. This check is not
    // guaranteed to work across threads. However, chances are that if two or
    // more threads simultaneously pass this test, they will be writing the
    // same error, which is fine. Even in the much less likely case that two
    // threads simultaneously write different error messages, the worst case is
    // that you get a mangled message. That's not good (and it's what we are
    // trying to avoid), but it's not critical.
    if (this->IsErrorRaised()) { return; }

    // Safely copy message into array.
    for (vtkm::Id index = 0; index < this->MessageBufferSize; index++)
      {
      this->MessageBuffer[index] = message[index];
      if (message[index] == '\0') { break; }
      }

    // Make sure message is null terminated.
    this->MessageBuffer[this->MessageBufferSize-1] = '\0';
  }

  VTKM_EXEC_CONT_EXPORT bool IsErrorRaised() const
  {
    if (this->MessageBufferSize > 0)
      {
      return (this->MessageBuffer[0] != '\0');
      }
    else
      {
      // If there is no buffer set, then always report an error.
      return true;
      }
  }

private:
  char *MessageBuffer;
  vtkm::Id MessageBufferSize;
};

}
}
} // namespace vtkm::exec::internal

#endif // vtkm_exec_internal_ErrorMessageBuffer_h