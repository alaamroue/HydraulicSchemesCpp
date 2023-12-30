/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
//#define DEBUG_OUTPUT
#define DEBUG_CELLX 2
#define DEBUG_CELLY 2


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <CL/opencl.h>
#include "GlobalHandlerClass.h"
#include "normalPlain.h"
#include "6_CLBoundaries.h"

//For Solver
#define QUITE_SMALL     1E-9
#define	VERY_SMALL      1E-10
#define DOMAIN_DIR_N	0
#define DOMAIN_DIR_E	1
#define DOMAIN_DIR_S	2
#define DOMAIN_DIR_W	3
#define	GRAVITY         9.81


//Compile time Definitions: for Cartesian Domain
#define	DOMAIN_CELLCOUNT 100
#define	DOMAIN_ROWS      10
#define	DOMAIN_COLS      10
#define	DOMAIN_DELTAX    1
#define	DOMAIN_DELTAY    1

//Dynamic Timesteps
#define TIMESTEP_EARLY_LIMIT			0.1
#define TIMESTEP_EARLY_LIMIT_DURATION	60.0
#define TIMESTEP_START_MINIMUM			1E-10
#define TIMESTEP_START_MINIMUM_DURATION	1.0
#define TIMESTEP_MINIMUM				1E-10
#define TIMESTEP_MAXIMUM				15.0

#define SCHEME_ENDTIME                  360000


//Boundries

// Hydrological timestep
// This should be low to capture velocities properly, but isn't
// always necessary
// TODO: Make configurable...
#define TIMESTEP_HYDROLOGICAL			1.0

// Boundary types
#define BOUNDARY_ATMOSPHERIC			0
#define BOUNDARY_FLOWCONDITIONS			1

// Boundary operating definitions
#define BOUNDARY_DEPTH_IGNORE			0
#define BOUNDARY_DEPTH_IS_FSL			1
#define BOUNDARY_DEPTH_IS_DEPTH			2
#define BOUNDARY_DEPTH_IS_CRITICAL		3

#define BOUNDARY_DISCHARGE_IGNORE		0
#define BOUNDARY_DISCHARGE_IS_DISCHARGE	1
#define BOUNDARY_DISCHARGE_IS_VELOCITY	2
#define BOUNDARY_DISCHARGE_IS_VOLUME	3

#define BOUNDARY_UNIFORM_RAIN_INTENSITY	0
#define BOUNDARY_UNIFORM_LOSS_RATE		1

#define BOUNDARY_GRIDDED_RAIN_INTENSITY 0
#define BOUNDARY_GRIDDED_RAIN_ACCUMUL	1
#define BOUNDARY_GRIDDED_MASS_FLUX		2


#define TIMESTEP_GROUPSIZE 12

#define COURANT_NUMBER 1
#define TIMESTEP_WORKERS 1


//Godnuov
#define GTS_DIM1 1
#define GTS_DIM2 1

cl_ulong	getCellID(cl_long lIdxX, cl_long lIdxY);
cl_ulong	getNeighbourByIndices(cl_long lIdxX, cl_long lIdxY, cl_uchar ucDirection);

cl_double4 riemannSolver(cl_uchar	ucDirection, cl_double8 pLeft, cl_double8 pRight, bool bDebug);

void gts_cacheDisabled(cl_double*,cl_double*,cl_double4*,cl_double4*,cl_double*, GlobalHandlerClass);
void solverFunctionPromaides(cl_double*, cl_double*, cl_double4*, cl_double4*, cl_double*, GlobalHandlerClass);
void rp(cl_double8, cl_double8);
cl_double8 d(cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double);

cl_uchar reconstructInterface(
	cl_double4		pStateLeft,
	cl_double		dBedLeft,
	cl_double4		pStateRight,
	cl_double		dBedRight,
	cl_double8* pOutputLeft,
	cl_double8* pOutputRight,
	cl_uchar		ucDirection
);

cl_double4 riemannSolver(cl_uchar	ucDirection, cl_double8 pLeft, cl_double8 pRight, bool bDebug);
void rp(cl_double8 d1, cl_double8 d2);
cl_double8 d(cl_double x1, cl_double x2, cl_double x3, cl_double x4, cl_double U, cl_double V, cl_double Zb, cl_double _);