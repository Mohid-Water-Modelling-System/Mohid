#include "CuWrapper.h"
#include "Generic.h"
#include <sstream>
#include <vector>

cudaDeviceProp CuWrapper::CuProperties;
dim3 CuWrapper::CuBlockDim(BLOCK_DIM, BLOCK_DIM, 1);
bool CuWrapper::CuInitialized = false;
size_t CuWrapper::CuPitchBytes = 0;

void CuWrapper::CuExe(cudaError_t cudaStatus, const char* errorMessage)
{
	if(cudaStatus != cudaSuccess)
	{
		cout << "A CUDA error occured: " << cudaGetErrorString(cudaStatus) << ".\n" << errorMessage << endl;
		CuCleanUp();
		exit(1);
	}
}

void CuWrapper::CuCleanUp()
{
	cudaDeviceSynchronize();
	// Reset the device for the current PROCESS (i.e. application). Does not reset the device for other applications.
	cudaDeviceReset();
	Generic::Print("Cleaned up CUDA.");
}

// WARNING: ignorePitch = true should only be used for testing!
void CuWrapper::CuInit(int device, bool ignorePitch)
{
	if(!CuInitialized)
	{
		Generic::Print("Initializing device");

		CuExe(cudaSetDevice(device), "cudaSetDevice failed, no device found");
		CuExe(cudaGetDeviceProperties(&CuProperties, 0), "Getting device properties failed");

		// Make sure the default block size does not exceed the maximum of the device
		if(CuBlockDim.x * CuBlockDim.y * CuBlockDim.z > CuProperties.maxThreadsPerBlock)
		{
			int dim = (int)floor(sqrt((float)CuProperties.maxThreadsPerBlock));
			dim = dim > 16 ? dim - dim % 16 : dim;
			CuBlockDim.x = dim;
			CuBlockDim.y = dim;
			CuBlockDim.z = 1;
		}

		// Determine the device pitch in bytes by allocating an integer array of dimension 1x1.
		// Use this to pad matrices on the host
		if(ignorePitch)
		{
			// WARNING: this will only make CuGetPitchSize work for CU_TYPE! IgnorePitch should only be used for testing
			CuPitchBytes = sizeof(CU_TYPE);
		}
		else
		{
			int *tmp;
			CuExe(cudaMallocPitch((void**)&tmp, &CuPitchBytes, sizeof(int), 1));
			CuExe(cudaFree(tmp));
		}

		CuInitialized = true;
	}
}