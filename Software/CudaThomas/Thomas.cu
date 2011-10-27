#include "Thomas.h"
#include "MatrixTranspose.h"

#pragma region C Binding
// Dimension 0 = X, 1 = Y, 2 = Z
void InitializeThomas_C(int *cudaObjID, T_Size3D *size)
{
	// NOTE: It is assumed that the lower bound is always 0, 
	// so the indices in the arrays always correspond correctly between FORTRAN and C++.
	// This actually deviates from the FORTRAN philosopy of 1-based arrays
	// When this method is called from FORTRAN, i = x, j = y, k = z
	int xDim = size->iUpperBound - size->iLowerBound + 1,
		yDim = size->jUpperBound - size->jLowerBound + 1,
		zDim = size->kUpperBound - size->kLowerBound + 1;
	Thomas::CreateInstance(*cudaObjID, xDim, yDim, zDim);
}

void KillThomas_C(int *cudaObjID)
{
	Thomas::KillInstance(*cudaObjID);
}

void SolveThomas_C(
	int *cudaObjID, 
	int *xLBound, int *xUBound, int *yLBound, int *yUBound,
	int *zLBound, int *zUBound,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
	CU_TYPE *TI, CU_TYPE *Res,
	int *dimension)
{
	Thomas *instance = Thomas::GetInstance(*cudaObjID);
	if(instance != NULL)
	{
		NSThomas::Dim dim = (NSThomas::Dim)*dimension;
		if(dim == NSThomas::X)
		{
			instance->SolveThomasX(*xLBound, *xUBound, *yLBound, *yUBound,
				*zLBound, *zUBound, D, E, F, TI, Res);
		}
		else if(dim == NSThomas::Y)
		{
			instance->SolveThomasY(*xLBound, *xUBound, *yLBound, *yUBound,
				*zLBound, *zUBound, D, E, F, TI, Res);
		}
		else
		{
			instance->SolveThomasZ(*xLBound, *xUBound, *yLBound, *yUBound,
				*zLBound, *zUBound, D, E, F, TI, Res);
		}
	}
}

void SaveThomas_C(int *cudaObjId, CU_TYPE *results, int *dim)
{
	Thomas::GetInstance(*cudaObjId)->SaveThomas(results, (NSThomas::Dim)*dim);
}
#pragma endregion

#pragma region Kernels
// Solve Thomas for dimension K where K is the major order and I is the minor order
// Solving for X or Z dimension
__global__ void DevThomasIK(
	int iLBound, int iUBound, int jLBound, int jUBound,
	int kLBound, int kUBound, size_t iPitchDim, size_t jDim,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F,
	CU_TYPE *W, CU_TYPE *G, CU_TYPE *TI,
	CU_TYPE *Res)
{
	int i = threadIdx.x + blockIdx.x * BLOCK_DIM,
		j = threadIdx.y + blockIdx.y * BLOCK_DIM;
	
	// Only threads that fall within the bounds should be used (inclusive bounds)
	if(i >= iLBound && i <= iUBound && j >= jLBound && j <= jUBound)
	{
		int ijDim = iPitchDim * jDim,
			index = i + j * iPitchDim + kLBound * ijDim;
		//1 + 1 * 3 + 1 * 3 * 3;
		CU_TYPE w = 0,
			wPrev = 0, 
			d = 0,
			e = E[index],
			edw = 0,
			g = 0,
			gPrev = 0,
			resPrev = 0;
		// Forward sweep (k = 0)
		// Start with first value
		wPrev = -F[index] / e;
		W[index] = wPrev;
		gPrev = TI[index] / e;
		// Need to set this in case kLBound == kUBound
		g = gPrev;
		G[index] = gPrev;

		// TODO: k could be replaced with index, using index = index; index <= ...; index += ijDim
		int k = 0;
		// Continue with second value, go up to last value
		for(k = kLBound + 1; k <= kUBound + 1; k++)
		{
			index += ijDim;

			// Calculate temp values
			d = D[index];
			edw = E[index] + d * wPrev;
			w = -F[index] / edw;
			g =  (TI[index] - d * gPrev) / edw;

			// Set actual values
			W[index] = w;
			G[index] = g;

			// Set temp values to prev values to use in next iteration
			wPrev = w;
			gPrev = g;
		}
		
		// Backward sweep
		// Start with last value
		// No need to recalculate first index, it's the same as from the last iteration of the previous loop
		// This also means that g is correct, no need to call G[index]
		resPrev = g;
		Res[index] = resPrev;
		// Continue with second last value, go down to first value
		for(k = kUBound; k >= kLBound; k--)
		{
			index -= ijDim;
			resPrev =  W[index] * resPrev + G[index];
			Res[index] = resPrev;
		}
	}
}

// Solve Thomas for dimension J where K is the major order and I is the minor order (requires I-K grid)
// Solving for Y dimension
__global__ void DevThomasJ(
	int iLBound, int iUBound, int jLBound, int jUBound,
	int kLBound, int kUBound, size_t iPitchDim, size_t jDim,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F,
	CU_TYPE *W, CU_TYPE *G, CU_TYPE *TI,
	CU_TYPE *Res)
{
	int i = threadIdx.x + blockIdx.x * BLOCK_DIM,
		k = threadIdx.y + blockIdx.y * BLOCK_DIM;
	
	// Only threads that fall within the bounds should be used (inclusive bounds)
	if(i >= iLBound && i <= iUBound && k >= kLBound && k <= kUBound)
	{
		int iPDim = iPitchDim,
			index = i + iPitchDim * jDim * k + jLBound * iPDim;
		
		CU_TYPE w = 0,
			wPrev = 0, 
			d = 0,
			e = E[index],
			edw = 0,
			g = 0,
			gPrev = 0,
			resPrev = 0;
		// Forward sweep (k = 0)
		// Start with first value
		wPrev = -F[index] / e;
		W[index] = wPrev;
		gPrev = TI[index] / e;
		// Need to set this in case kLBound == kUBound
		g = gPrev;
		G[index] = gPrev;

		int j = 0;
		// Continue with second value, go up to last value
		for(j = jLBound + 1; j <= jUBound; j++)
		{
			index += iPDim;

			// Calculate temp values
			d = D[index];
			edw = E[index] + d * wPrev;
			w = -F[index] / edw;
			g =  (TI[index] - d * gPrev) / edw;

			// Set actual values
			W[index] = w;
			G[index] = g;

			// Set temp values to prev values to use in next iteration
			wPrev = w;
			gPrev = g;
		}
		
		// Backward sweep
		// Start with last value
		// No need to recalculate first index, it's the same as from the last iteration of the previous loop
		// This also means that g is correct, no need to call G[index]
		resPrev = g;
		Res[index] = resPrev;
		// Continue with second last value, go down to first value
		for(j = jUBound - 1; j >= jLBound; j--)
		{
			index -= iPDim;
			resPrev =  W[index] * resPrev + G[index];
			Res[index] = resPrev;
		}
	}
}
#pragma endregion

map<int, Thomas*>Thomas::Instances;

Thomas::Thomas(int cudaID, int xDim, int yDim, int zDim)
{
	Sizes.x = xDim;
	Sizes.y = yDim;
	Sizes.z = zDim;
	this->CudaID = cudaID;

	Initialize();
}

Thomas::~Thomas()
{
	cudaFree(DevW_XMinor);
	cudaFree(DevG_XMinor);
	cudaFree(DevD_XMinor);
	cudaFree(DevE_XMinor);
	cudaFree(DevF_XMinor);
	cudaFree(DevTI_XMinor);
	cudaFree(DevRes_XMinor);

	cudaFree(DevW_ZMinor);
	cudaFree(DevG_ZMinor);
	cudaFree(DevD_ZMinor);
	cudaFree(DevE_ZMinor);
	cudaFree(DevF_ZMinor);
	cudaFree(DevTI_ZMinor);
	cudaFree(DevRes_ZMinor);

	cudaFree(DevHostD);
	cudaFree(DevHostE);
	cudaFree(DevHostF);
	cudaFree(DevHostTI);
	cudaFree(DevHostRes);

	cudaStreamDestroy(streamD);
	cudaStreamDestroy(streamE);
	cudaStreamDestroy(streamF);
	cudaStreamDestroy(streamTI);
	cudaStreamDestroy(streamRes);
}

void Thomas::CreateInstance(int cudaObjID, int xDim, int yDim, int zDim)
{
	// If no instance with this ID exists, create a new one
	if(Instances.count(cudaObjID) == 0)
	{
		// Create a new instance
		Instances[cudaObjID] = new Thomas(cudaObjID, xDim, yDim, zDim);
	}
}

Thomas *Thomas::GetInstance(int cudaObjID)
{
	if(Instances.count(cudaObjID) == 0)
	{
		return NULL;
	}
	return Instances[cudaObjID];
}

void Thomas::KillInstance(int cudaObjID)
{
	if(Instances[cudaObjID] != 0)
	{
		// Call the deconstructor
		delete Instances[cudaObjID];
		// Remove the (nullified) item from the map
		Instances.erase(cudaObjID);
	}
}

void Thomas::Initialize()
{
	// Source matrix is always in X,Y,Z ordering where X = minor, Y = medium, Z = major
	// If calculating in X direction, transpose to Z,Y,X so X becomes major and Z becomes minor
	// If calculating in Y direction, no need to transpose (use X-Z grid)
	// If calculating in Z direction, no need to transpose (use X-Y grid)
	// Rule of thumb: calculation direction cannot be minor

	PitchX = CuWrapper::CuGetPitchBytes<CU_TYPE>(Sizes.x);
	PitchZ = CuWrapper::CuGetPitchBytes<CU_TYPE>(Sizes.z);
	// PitchSizes equals element count of minor dimension (column)
	PitchXSize = PitchX / sizeof(CU_TYPE);
	PitchZSize = PitchZ / sizeof(CU_TYPE);

	MemXSize = PitchX * Sizes.y * Sizes.z;
	MemZSize = PitchZ * Sizes.y * Sizes.x;

	InitializeVariable(&DevW_XMinor, MemXSize);
	InitializeVariable(&DevG_XMinor, MemXSize);
	InitializeVariable(&DevD_XMinor, MemXSize);
	InitializeVariable(&DevE_XMinor, MemXSize);
	InitializeVariable(&DevF_XMinor, MemXSize);
	InitializeVariable(&DevTI_XMinor, MemXSize);
	InitializeVariable(&DevRes_XMinor, MemXSize);

	InitializeVariable(&DevW_ZMinor, MemZSize);
	InitializeVariable(&DevG_ZMinor, MemZSize);
	InitializeVariable(&DevD_ZMinor, MemZSize);
	InitializeVariable(&DevE_ZMinor, MemZSize);
	InitializeVariable(&DevF_ZMinor, MemZSize);
	InitializeVariable(&DevTI_ZMinor, MemZSize);
	InitializeVariable(&DevRes_ZMinor, MemZSize);

	// Allocate device memory in X,Y,Z ordering (X-minor). Host matrix is copied here before transposing (only for X dimension)
	InitializeVariable(&DevHostD, MemXSize);
	InitializeVariable(&DevHostE, MemXSize);
	InitializeVariable(&DevHostF, MemXSize);
	InitializeVariable(&DevHostTI, MemXSize);
	InitializeVariable(&DevHostRes, MemXSize);

#ifdef _USE_PAGELOCKED
	CuWrapper::CuExe(cudaStreamCreate(&streamD));
	CuWrapper::CuExe(cudaStreamCreate(&streamE));
	CuWrapper::CuExe(cudaStreamCreate(&streamF));
	CuWrapper::CuExe(cudaStreamCreate(&streamTI));
	CuWrapper::CuExe(cudaStreamCreate(&streamRes));
#else
	// Only use streams if _USE_PAGELOCKED is defined. Async copy cannot be used with page-able memory
	streamD = (cudaStream_t)0;
	streamE = (cudaStream_t)0;
	streamF = (cudaStream_t)0;
	streamTI = (cudaStream_t)0;
	streamRes = (cudaStream_t)0;
#endif _USE_PAGELOCKED

	CuBlockDim = CuWrapper::CuBlockDim;

	CuGridDimZY = dim3(
		ceil((float)Sizes.z / this->CuBlockDim.x),
		ceil((float)Sizes.y / this->CuBlockDim.y));
	CuGridDimXY = dim3(
		ceil((float)Sizes.x / this->CuBlockDim.x),
		ceil((float)Sizes.y / this->CuBlockDim.y));
	CuGridDimXZ = dim3(
		ceil((float)Sizes.x / this->CuBlockDim.x),
		ceil((float)Sizes.z / this->CuBlockDim.y));
	CuGridDimZX = dim3(
		ceil((float)Sizes.z / this->CuBlockDim.x),
		ceil((float)Sizes.x / this->CuBlockDim.y));

	for(int i = 0; i < 3; i++)
	{
		OutputCounter[i] = 0;
		FileCounter[i] = 0;
	}
}

void Thomas::InitializeVariable(CU_TYPE **devPtr, size_t size)
{
	CuWrapper::CuExe(cudaMalloc((void**)devPtr, size));
	CuWrapper::CuExe(cudaMemset(*devPtr, 0, size));
}

void Thomas::CopyToDevice(CU_TYPE *devPtr, CU_TYPE *hostPtr, cudaStream_t stream)
{
	// Matrix is always X minor when copied to device, so use MemXSize
	if(stream == NULL)
	{
		CuWrapper::CuExe(cudaMemcpy(devPtr, hostPtr, MemXSize, cudaMemcpyHostToDevice));
	}
	else
	{
		CuWrapper::CuExe(cudaMemcpyAsync(devPtr, hostPtr, MemXSize, cudaMemcpyHostToDevice, stream));
	}
}

void Thomas::CopyToHost(CU_TYPE *hostPtr, CU_TYPE *devPtr, cudaStream_t stream)
{
	// Matrix is always X minor when copied to host, so use MemXSize
	if(stream == NULL)
	{
		CuWrapper::CuExe(cudaMemcpy(hostPtr, devPtr, MemXSize, cudaMemcpyDeviceToHost));
	}
	else
	{
		CuWrapper::CuExe(cudaMemcpyAsync(hostPtr, devPtr, MemXSize, cudaMemcpyDeviceToHost, stream));
	}
}

void Thomas::SolveThomasX(
	int xLBound, int xUBound, int yLBound, int yUBound,
	int zLBound, int zUBound,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
	CU_TYPE *TI, CU_TYPE *Res)
{
	// Transpose X-minor to Z-minor
	CopyToDevice(DevHostD, D, streamD);
	MatrixTranspose::TransposeItoK(DevD_ZMinor, DevHostD, Sizes.x, Sizes.y, Sizes.z, 
		PitchXSize, PitchZSize, CuGridDimXZ, CuBlockDim, streamD);
		
	CopyToDevice(DevHostE, E, streamE);
	MatrixTranspose::TransposeItoK(DevE_ZMinor, DevHostE, Sizes.x, Sizes.y, Sizes.z, 
		PitchXSize, PitchZSize, CuGridDimXZ, CuBlockDim, streamE);

	CopyToDevice(DevHostF, F, streamF);
	MatrixTranspose::TransposeItoK(DevF_ZMinor, DevHostF, Sizes.x, Sizes.y, Sizes.z, 
		PitchXSize, PitchZSize, CuGridDimXZ, CuBlockDim, streamF);
		
	CopyToDevice(DevHostTI, TI, streamTI);
	MatrixTranspose::TransposeItoK(DevTI_ZMinor, DevHostTI, Sizes.x, Sizes.y, Sizes.z, 
		PitchXSize, PitchZSize, CuGridDimXZ, CuBlockDim, streamTI);

	CopyToDevice(DevHostRes, Res, streamRes);
	MatrixTranspose::TransposeItoK(DevRes_ZMinor, DevHostRes, Sizes.x, Sizes.y, Sizes.z, 
		PitchXSize, PitchZSize, CuGridDimXZ, CuBlockDim, streamRes);

	// Synchronize device, streams are used
	// (theoretically speaking not necessary because stream 0 always waits for all streams to finish)
	CuWrapper::CuExe(cudaDeviceSynchronize());

	// Solve thomas for X dimension
	DevThomasIK<<<CuGridDimZY, CuBlockDim>>>(
		zLBound, zUBound,
		yLBound, yUBound,
		xLBound, xUBound,
		PitchZSize, Sizes.y,
		DevD_ZMinor, DevE_ZMinor, DevF_ZMinor,
		DevW_ZMinor, DevG_ZMinor, 
		DevTI_ZMinor, DevRes_ZMinor);

	// Transpose result matrix back to X-minor
	MatrixTranspose::TransposeItoK(DevHostRes, DevRes_ZMinor, Sizes.z, Sizes.y, Sizes.x, 
		PitchZSize, PitchXSize, CuGridDimZX, CuBlockDim);
	
	// Copy result matrix back to host
	CopyToHost(Res, DevHostRes);
}

void Thomas::SolveThomasY(
	int xLBound, int xUBound, int yLBound, int yUBound,
	int zLBound, int zUBound,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
	CU_TYPE *TI, CU_TYPE *Res)
{
	// Copy input to device
	CopyToDevice(DevD_XMinor, D);
	CopyToDevice(DevE_XMinor, E);
	CopyToDevice(DevF_XMinor, F);
	CopyToDevice(DevTI_XMinor, TI);
	CopyToDevice(DevRes_XMinor, Res);

	// Execute Thomas for Y dimension
	DevThomasJ<<<CuGridDimXZ, CuBlockDim>>>(
		xLBound, xUBound,
		yLBound, yUBound,
		zLBound, zUBound,
		PitchXSize, Sizes.y,
		DevD_XMinor, DevE_XMinor, DevF_XMinor,
		DevW_XMinor, DevG_XMinor, 
		DevTI_XMinor, DevRes_XMinor);

	// Copy results back to host
	CopyToHost(Res, DevRes_XMinor);
}

void Thomas::SolveThomasZ(
	int xLBound, int xUBound, int yLBound, int yUBound,
	int zLBound, int zUBound,
	CU_TYPE *D, CU_TYPE *E, CU_TYPE *F, 
	CU_TYPE *TI, CU_TYPE *Res)
{
	// Copy input to device
	CopyToDevice(DevD_XMinor, D);
	CopyToDevice(DevE_XMinor, E);
	CopyToDevice(DevF_XMinor, F);
	CopyToDevice(DevTI_XMinor, TI);
	// TODO: find a way to avoid this copy. 
	// Could be done by storing the values outside the worksize into temporary arrays and write them back to Res
	CopyToDevice(DevRes_XMinor, Res);

	// Execute Thomas for Z dimension
	DevThomasIK<<<CuGridDimXY, CuBlockDim>>>(
		xLBound, xUBound,
		yLBound, yUBound,
		zLBound, zUBound,
		PitchXSize, Sizes.y,
		DevD_XMinor, DevE_XMinor, DevF_XMinor,
		DevW_XMinor, DevG_XMinor, 
		DevTI_XMinor, DevRes_XMinor);

	// Copy results back to host
	CopyToHost(Res, DevRes_XMinor);
}

void Thomas::SaveThomas(CU_TYPE *results, NSThomas::Dim dimension)
{
	int dim = (int)dimension;
	// Output every X time steps
	OutputCounter[dim] = OutputCounter[dim] % 25;
	if(OutputCounter[dim] == 0)
	{
		string dimStr;
		if(dim == 0)
			dimStr = "X";
		else if(dim == 1)
			dimStr = "Y";
		else
			dimStr = "Z";

		FileCounter[dim]++;
		cout << "Saving Thomas output " << FileCounter[dim] << " for model " << CudaID << ", dimension " << dimStr << endl << endl;

		stringstream fileName;
		fileName << "output_model" << CudaID << "_dim" << dimStr << "_" << FileCounter[dim] << ".txt";
		ofstream outputFile = ofstream(fileName.str());
		char buff[50];
		for(int z = 0; z < Sizes.z; z++)
		{
			outputFile << "Z = " << z << "---------------------------------------------" << endl;
			for(int y = 0; y < Sizes.y; y++)
			{
				int zy = z * PitchXSize * Sizes.y + y * PitchXSize;
				// One line for every row
				for(int x = 0; x < Sizes.x; x++)
				{
					CU_TYPE val = results[x + zy];
					// Sometimes FORTRAN has 0, C++ has -0
					if(val == -0)
					{
						val = 0;
					}
					else if(val <= -9900000000000000 || val >= 9900000000000000)
					{
						// Avoid having to use large padding before decimal point
						val = _Nan._Double;
					}
					// Print in 25 digits precision
					sprintf(buff, "%31.25f\t", val);
					outputFile << buff;
				}
				outputFile << endl;
			}
			outputFile << endl;
		}
		outputFile << endl;
		outputFile.close();
		OutputCounter[dim] = 0;
	}
	OutputCounter[dim]++;
}