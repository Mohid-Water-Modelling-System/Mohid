// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#pragma once

// CU_TYPE = Floating point precision. Compute capability < 1.3 does not support double precision floating points, >= 1.3 does.
// Switch to single precision floating points for < 1.3
// All other CU_ constants are defined in CuWrapper.h
// Somehow declaring this in CuWrapper.h gives errors
#ifndef CU_TYPE_DEF
typedef double CU_TYPE;
#define CU_TYPE_DEF
#endif

// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <math.h>

using namespace std;

// CUDA Header Files
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Own Header Files. IMPORTANT: Keep in this order to prevent compilation errors
#include "Generic.h"
#include "CuWrapperBinding.h"
#include "CuWrapper.h"
#include "MatrixTranspose.h"
#include "Timer.h"