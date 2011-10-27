#pragma once
#include "stdafx.h"
#include <map>
#include <fstream>
#include <iomanip>

struct T_Size3D;
struct T_Size2D;

// C Binding
extern "C"
{
	void InitializeThomas_C(int *cudaObjID, T_Size3D *size);
	void KillThomas_C(int *cudaObjID);

	void SolveThomas_C(
		int *cudaObjID,
		int *xLBound, int *xUBound, int *yLBound, int *yUBound,
		int *zLBound, int *zUBound,
		CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
		CU_TYPE *TI, CU_TYPE *Res,
		int *dimensions);

	void SaveThomas_C(int* cudaObjId, CU_TYPE *results, int *dimension);
}

// Kernels
__global__ void DevThomasIK(
	int iLBound, int iUBound, int jLBound, int jUBound,
	int kLBound, int kUBound, size_t iPitch, size_t jDim,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F,
	CU_TYPE *W, CU_TYPE *G, 
	CU_TYPE *TI, CU_TYPE *Res);

__global__ void DevThomasJ(
	int iLBound, int iUBound, int jLBound, int jUBound,
	int kLBound, int kUBound, size_t iPitch, size_t jDim,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F,
	CU_TYPE *W, CU_TYPE *G, 
	CU_TYPE *TI, CU_TYPE *Res);

namespace NSThomas
{
	enum Dim
	{
		X = 0,
		Y = 1,
		Z = 2
	};
}

class Thomas
{
private:
	Thomas(int cudaID, int xDim, int yDim, int zDim);

	int CudaID;

	dim3 
		CuBlockDim,
		CuGridDimZY,
		CuGridDimXY,
		CuGridDimXZ,
		CuGridDimZX;

	dim3 Sizes;

	cudaStream_t 
		streamD, streamE, streamF,
		streamTI, streamRes;

	size_t 
		PitchX,
		PitchZ,
		PitchXSize,
		PitchZSize,
		MemXSize,
		MemZSize;

	// Device arrays (transpose to these arrays)
	CU_TYPE
		*DevW_XMinor, *DevG_XMinor,
		*DevD_XMinor, *DevE_XMinor, *DevF_XMinor,
		*DevTI_XMinor, *DevRes_XMinor;

	CU_TYPE
		*DevW_ZMinor, *DevG_ZMinor,
		*DevD_ZMinor, *DevE_ZMinor, *DevF_ZMinor,
		*DevTI_ZMinor, *DevRes_ZMinor;

	// Host arrays on device (copy from host to these arrays)
	CU_TYPE
		*DevHostD, *DevHostE, *DevHostF,
		*DevHostTI, *DevHostRes;

	// Testing only
	int OutputCounter[3];
	int FileCounter[3];

	void Initialize();
	void InitializeVariable(CU_TYPE **devPtr, size_t size);
	void CopyToDevice(CU_TYPE *devPtr, CU_TYPE *hostPtr, cudaStream_t stream = (cudaStream_t)0);
	void CopyToHost(CU_TYPE *hostPtr, CU_TYPE *devPtr, cudaStream_t stream = (cudaStream_t)0);

	static map<int, Thomas*> Instances;
public:
	~Thomas();
	
	static void CreateInstance(int cudaObjID, int xDim, int yDim, int zDi);
	static Thomas *GetInstance(int cudaObjID);
	static void KillInstance(int cudaObjID);

	void SolveThomasX(
		int xLBound, int xUBound, int yLBound, int yUBound,
		int zLBound, int zUBound,
		CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
		CU_TYPE *TI, CU_TYPE *Res);
	void SolveThomasY(
		int xLBound, int xUBound, int yLBound, int yUBound,
		int zLBound, int zUBound,
		CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
		CU_TYPE *TI, CU_TYPE *Res);
	void SolveThomasZ(
		int xLBound, int xUBound, int yLBound, int yUBound,
		int zLBound, int zUBound,
		CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
		CU_TYPE *TI, CU_TYPE *Res);

	// Testing only
	void SaveThomas(CU_TYPE *results, NSThomas::Dim dimension);
};