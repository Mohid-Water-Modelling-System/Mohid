#pragma once
#include "stdafx.h"

// C Linkage
extern "C"
{
	void ConstructCudaBinding_C();
	void KillCudaBinding_C();
	void Alloc3DPageLocked_C(CU_TYPE **ptr, int *xDim, int *yDim, int *zDim);
	void FreePageLocked_C(CU_TYPE **ptr);
}

// Equivalent to size types in FORTRAN
struct T_Size3D
{
	int iLowerBound,
		iUpperBound,
		jLowerBound,
		jUpperBound,
		kLowerBound,
		kUpperBound;
};

struct T_Size2D
{
	int iLowerBound,
		iUpperBound,
		jLowerBound,
		jUpperBound;
};