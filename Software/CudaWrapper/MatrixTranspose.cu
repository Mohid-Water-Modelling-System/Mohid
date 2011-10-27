#include "MatrixTranspose.h"

// Transpose a matrix from I to J
// Rules of thumb:
// - I should always be the minor dimension of the source matrix (column).
// - J should always be the medium dimension of the source matrix (row).
// - I will become the medium dimension of the source matrix (row).
// - J will become the minor dimension of the destination matrix (column).
__global__ void DevMatrixTransposeItoJ(CU_TYPE *dst, CU_TYPE *src, 
	int iDim, int jDim, int kDim, int devIDim, int devJDim)
{
	// BLOCK_DIM + 1 prevents bank conflicts when reading from shared memory
	__shared__ CU_TYPE tile[BLOCK_DIM][BLOCK_DIM + 1];
	// Use HI/LO to use less registers
	int tId = (threadIdx.x << 16) + threadIdx.y,
		offset = ((blockIdx.x * BLOCK_DIM) << 16) + blockIdx.y * BLOCK_DIM,
		srcOffset = HIWORD(tId) + HIWORD(offset) + (LOWORD(offset) + LOWORD(tId)) * devIDim,
		dstOffset = HIWORD(tId) + HIWORD(offset) * devJDim + LOWORD(offset) + LOWORD(tId) * devJDim,
		iMax = HIWORD(offset) + BLOCK_DIM > iDim ? iDim % BLOCK_DIM : BLOCK_DIM,
		jMax = LOWORD(offset) + BLOCK_DIM > jDim ? jDim % BLOCK_DIM : BLOCK_DIM;

	for(int k = 0; k < kDim; k++)
	{
		if(HIWORD(tId) < iMax && LOWORD(tId) < jMax)
		{
			// Src address = i + j * devIDim + k * devIDim * jDim
			tile[HIWORD(tId)][LOWORD(tId)] = src[srcOffset + k * devIDim * jDim];
		}
		__syncthreads();
		if(HIWORD(tId) < jMax && LOWORD(tId) < iMax)
		{
			dst[dstOffset + k * iDim * devJDim] = tile[LOWORD(tId)][HIWORD(tId)];
		}
	}
}

// Transpose a matrix from I to K
// Rules of thumb:
// - I should always be the minor dimension of the source matrix (column)
// - K should always be the major dimension of the source matrix (page)
// - I will become the major dimension of the source matrix (page)
// - K will become the minor dimension of the destination matrix (column)
__global__ void DevMatrixTransposeItoK(CU_TYPE *dst, CU_TYPE *src, 
	int iDim, int jDim, int kDim, int devIDim, int devKDim)
{
	// BLOCK_DIM + 1 prevents bank conflicts when reading from shared memory
	__shared__ CU_TYPE tile[BLOCK_DIM][BLOCK_DIM + 1];
	// Use HI/LO to use less registers
	int tId = (threadIdx.x << 16) + threadIdx.y,
		offset = ((blockIdx.x * BLOCK_DIM) << 16) + blockIdx.y * BLOCK_DIM,
		srcOffset = HIWORD(tId) + HIWORD(offset) + (LOWORD(offset) + LOWORD(tId)) * devIDim * jDim,
		dstOffset = HIWORD(tId) + HIWORD(offset) * devKDim * jDim + LOWORD(offset) + LOWORD(tId) * devKDim * jDim,
		iMax = HIWORD(offset) + BLOCK_DIM > iDim ? iDim % BLOCK_DIM : BLOCK_DIM,
		kMax = LOWORD(offset) + BLOCK_DIM > kDim ? kDim % BLOCK_DIM : BLOCK_DIM;
	
	for(int j = 0; j < jDim; j++)
	{
		if(HIWORD(tId) < iMax && LOWORD(tId) < kMax)
		{
			tile[HIWORD(tId)][LOWORD(tId)] = src[srcOffset + j * devIDim];
		}
		__syncthreads();
		if(HIWORD(tId) < kMax && LOWORD(tId) < iMax)
		{
			dst[dstOffset + j * devKDim] = tile[LOWORD(tId)][HIWORD(tId)];
		}
	}
}

// Transpose a matrix from J to K
// Rules of thumb:
// - I should always be the minor dimension of the source matrix (column)
// - J should always be the medium dimension of the source matrix (row)
// - K should always be the major dimension of the source matrix (page)
// - I will stay the minor dimension of the source matrix (column)
// - J will become the major dimension of the destination matrix (page)
// - K will become the medium dimension of the destination matrix (row)
__global__ void DevMatrixTransposeJtoK(CU_TYPE *dst, CU_TYPE *src, 
	int iDim, int jDim, int kDim, int devIDim)
{
	// No need for shared memory, since we do a coalesced read and a coalesced write (keep things simple for a change
	// No need to use HI/LO because we only use 4 registers
	int i = (blockIdx.x * BLOCK_DIM + threadIdx.x),
		j = (blockIdx.y * BLOCK_DIM + threadIdx.y),
		srcOffset = i + j * devIDim,
		dstOffset = i + j * devIDim * kDim;

	if(i < iDim && j < jDim)
	{
		for(int k = 0; k < kDim; k++)
		{
			// Coaleseced read to coalesced write
			dst[dstOffset + k * devIDim] = src[srcOffset + k * devIDim * jDim];
		}
	}
}

MatrixTranspose::MatrixTranspose()
{
}

MatrixTranspose::~MatrixTranspose()
{
}

void MatrixTranspose::TransposeItoJ(CU_TYPE *devDst, CU_TYPE *devSrc, int iDim, int jDim, int kDim, 
	size_t iPitchSize, size_t jPitchSize, dim3 gridDims, dim3 blockDims, cudaStream_t stream)
{
	DevMatrixTransposeItoK<<<gridDims, blockDims, 0, stream>>>(devDst, devSrc, iDim, jDim, kDim, iPitchSize, jPitchSize);
}

void MatrixTranspose::TransposeItoK(CU_TYPE *devDst, CU_TYPE *devSrc, int iDim, int jDim, int kDim, 
	size_t iPitchSize, size_t kPitchSize, dim3 gridDims, dim3 blockDims, cudaStream_t stream)
{
	DevMatrixTransposeItoK<<<gridDims, blockDims, 0, stream>>>(devDst, devSrc, iDim, jDim, kDim, iPitchSize, kPitchSize);
}

void MatrixTranspose::TransposeJtoK(CU_TYPE *devDst, CU_TYPE *devSrc, int iDim, int jDim, int kDim, 
	size_t iPitchSize, dim3 gridDims, dim3 blockDims, cudaStream_t stream)
{
	DevMatrixTransposeJtoK<<<gridDims, blockDims, 0, stream>>>(devDst, devSrc, iDim, jDim, kDim, iPitchSize);
}
