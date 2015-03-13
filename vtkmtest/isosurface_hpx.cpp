// tell vtkm that we want to use HPX backend
#define VTKM_DEVICE_ADAPTER VTKM_DEVICE_ADAPTER_HPX

// squash windows #defines before #including hpx
#define NOMINMAX

// override int main() so that hpx is initialized on startup
#include <hpx/hpx_main.hpp>  

// squash windows #defines FTER #INCLUDING HPX
#undef GetMessage

// if we have FreeGlut, change event handler
// #ifdef VTKM_USE_FREEGLUT
#undef VTKM_USE_FREEGLUT

#define HPX_TIMING

#include "isosurface.cpp"

