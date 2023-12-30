/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include "definitions.h"

//Calculate the timestep using a reduction procedure and increment the total model time.


#ifdef USE_FUNCTION_STUBS
 // Function definitions
__kernel  __attribute__((reqd_work_group_size(1, 1, 1)))
void tst_Advance_Normal(
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double4*,
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_uint*,
	__global	cl_uint*
);

__kernel  __attribute__((reqd_work_group_size(1, 1, 1)))
void tst_ResetCounters(
	__global	cl_double*,
	__global	cl_uint*,
	__global	cl_uint*
);

__kernel  __attribute__((reqd_work_group_size(1, 1, 1)))
void tst_UpdateTimestep(
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double*,
	__global	cl_double*
);

__kernel  REQD_WG_SIZE_LINE
void tst_Reduce(
	__global	cl_double4*,
	__global	cl_double const* restrict,
	__global	cl_double*
);

#endif
