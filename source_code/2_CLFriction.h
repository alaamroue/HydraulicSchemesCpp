/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include "definitions.h"

//Calculate the timestep using a reduction procedure and increment the total model time.

void per_Friction(
	cl_double*,
	cl_double4*,
	cl_double*,
	cl_double*,
	cl_double* 
);

cl_double4 implicitFriction(
	cl_double4,
	cl_double,
	cl_double,
	cl_double
);