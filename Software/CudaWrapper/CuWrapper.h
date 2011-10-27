#pragma once
#include "stdafx.h"
#include <vector>

const int BLOCK_DIM = 16;

class CuWrapper
{
private:
	static bool CuInitialized;
public:
	static cudaDeviceProp CuProperties;
	static dim3 CuBlockDim;
	static size_t CuPitchBytes;

	static void CuExe(cudaError_t cudaStatus, const char* errorMessage = "");
	static void CuCleanUp();

	static void CuInit(int device, bool ignorePitch = false);

	/// Get the pitch for a given dimension (in bytes)
	template <class T> inline static int CuGetPitchBytes(int minorDim)
	{
		return (int)ceil((float)(minorDim * sizeof(T)) / CuPitchBytes) * CuPitchBytes;
	}

	/// Allocate pagelocked padded host memory to replace cudaMallocPitch.
	/// cudaMemcpy2D cannot be used with streams, so use padded host matrix with cudaMemcpy.
	/// MinorDim will be padded (first dimension in FORTRAN, last dimension in C)
	template <class T> inline static void CuHostAllocPitch(T **ptr, int minorDim, int mediumDim, int majorDim)
	{
		int minorPitch = CuGetPitchBytes<T>(minorDim);
		CuExe(cudaHostAlloc((void**)ptr, minorPitch * mediumDim * majorDim, cudaHostAllocDefault));
		memset(*ptr, 0, minorPitch * mediumDim * majorDim);
	}

	template <class T> inline static void CuHostFreePitch(T *ptr)
	{
		CuExe(cudaFreeHost(ptr));
	}

	/// Get the pitch for a given dimension (in elements of size sizeof(T))
	template <class T> inline static int CuGetPitchSize(int minorDim)
	{
		return CuGetPitchBytes<T>(minorDim) / sizeof(T);
	}
};