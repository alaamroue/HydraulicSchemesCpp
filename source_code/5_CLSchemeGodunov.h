/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include "definitions.h"

//Implementation of the 1st order accurate Godunov - type scheme for execution on the GPU.

// Structure definitions
typedef struct sFaceStructure {
	cl_double4	pN;
	cl_double4	pE;
	cl_double4	pS;
	cl_double4	pW;
} sFaceStructure;

#ifdef USE_FUNCTION_STUBS

// Function definitions
__kernel  REQD_WG_SIZE_FULL_TS
void gts_cacheDisabled(
	__constant	cl_double*,
	__global	cl_double const* restrict,
	__global	cl_double4*,
	__global	cl_double4*,
	__global    cl_double const* restrict
);

__kernel  REQD_WG_SIZE_FULL_TS
void gts_cacheEnabled(
	__constant	cl_double*,
	__global	cl_double const* restrict,
	__global	cl_double4*,
	__global	cl_double4*,
	__global    cl_double const* restrict
);

cl_uchar reconstructInterface(
	cl_double4,
	cl_double,
	cl_double4,
	cl_double,
	cl_double8*,
	cl_double8*,
	cl_uchar
);

#endif

