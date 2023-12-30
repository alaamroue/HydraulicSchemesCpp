/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include "definitions.h"

//Implementation of the approximate HLLC Riemann solver for the GPU.

#ifdef USE_FUNCTION_STUBS

// Function definitions
cl_double4	riemannSolver(
	cl_uchar,
	cl_double8,
	cl_double8,
	bool
);

#endif