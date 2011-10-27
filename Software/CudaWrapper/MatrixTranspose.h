#pragma once
#include "stdafx.h"

__global__ void DevMatrixTransposeItoJ(float *dst, float *src, 
	int iDim, int jDim, int kDim, int devIDim, int devJDim);

__global__ void DevMatrixTransposeItoK(float *dst, float *src, 
	int iDim, int jDim, int kDim, int devIDim, int devKDim);

__global__ void DevMatrixTransposeJtoK(float *dst, float *src, 
	int iDim, int jDim, int kDim, int devIDim);

__global__ void DevMatrixTransposeItoJNonCoalesced(float *dst, float *src, 
	int iDim, int jDim, int kDim, int devIDim, int devJDim);

__global__ void DevMatrixTransposeItoKNonCoalesced(float *dst, float *src, 
	int iDim, int jDim, int kDim, int devIDim, int devKDim);

class MatrixTranspose
{
private:
	MatrixTranspose();
public:
	~MatrixTranspose();

	static void TransposeItoJ(CU_TYPE *devDst, CU_TYPE *devSrc, int iDim, int jDim, int kDim, 
		size_t iPitchSize, size_t jPitchSize, dim3 gridDims, dim3 blockDims, cudaStream_t stream = (cudaStream_t)0);
	static void TransposeItoK(CU_TYPE *devDst, CU_TYPE *devSrc, int iDim, int jDim, int kDim, 
		size_t iPitchSize, size_t kPitchSize, dim3 gridDims, dim3 blockDims, cudaStream_t stream = (cudaStream_t)0);
	static void TransposeJtoK(CU_TYPE *devDst, CU_TYPE *devSrc, int iDim, int jDim, int kDim, 
		size_t iPitchSize, dim3 gridDims, dim3 blockDims, cudaStream_t stream = (cudaStream_t)0);
};