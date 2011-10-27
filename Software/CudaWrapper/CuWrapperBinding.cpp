#include "CuWrapperBinding.h"

// C Linkage
void ConstructCudaBinding_C()
{
#ifdef _USE_PAGELOCKED
	CuWrapper::CuInit(0, false);
#else
	// Don't use pitch if no page-locked memory is used
	CuWrapper::CuInit(0, true);
#endif
}

void KillCudaBinding_C()
{
	CuWrapper::CuCleanUp();
}

void Alloc3DPageLocked_C(CU_TYPE **ptr, int *xDim, int *yDim, int *zDim)
{
	CuWrapper::CuHostAllocPitch(ptr, *xDim, *yDim, *zDim);
	// Set the correct pitch so the calling method knows what the alignment is
	xDim[0] = CuWrapper::CuGetPitchSize<CU_TYPE>(*xDim);
}

void FreePageLocked_C(CU_TYPE **ptr)
{
	CuWrapper::CuHostFreePitch(*ptr);
}