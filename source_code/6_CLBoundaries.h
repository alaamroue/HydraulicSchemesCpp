/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include "definitions.h"


//Management functions for a domain boundaries.

//#ifdef USE_FUNCTION_STUBS

typedef struct sBdyCellConfiguration
{
	cl_ulong		TimeseriesEntries;
	cl_double		TimeseriesInterval;
	cl_double		TimeseriesLength;
	cl_ulong		RelationCount;
	cl_uint			DefinitionDepth;
	cl_uint			DefinitionDischarge;
} sBdyCellConfiguration;

typedef struct sBdyGriddedConfiguration
{
	cl_double		TimeseriesInterval;
	cl_double		GridResolution;
	cl_double		GridOffsetX;
	cl_double		GridOffsetY;
	cl_ulong		TimeseriesEntries;
	cl_ulong		Definition;
	cl_ulong		GridRows;
	cl_ulong		GridCols;
} sBdyGriddedConfiguration;

typedef struct sBdyUniformConfiguration
{
	cl_uint			TimeseriesEntries;
	cl_double		TimeseriesInterval;
	cl_double		TimeseriesLength;
	cl_uint			Definition;
} sBdyUniformConfiguration;

void bdy_Cell(
	sBdyCellConfiguration*,
	cl_ulong*,
	cl_double4*,
	cl_double*,
	cl_double*,
	cl_double*,
	cl_double4*,
	cl_double*,
	cl_double*
);

void bdy_Gridded(
	sBdyGriddedConfiguration*,
	cl_double*,
	cl_double*,
	cl_double*,
	cl_double*,
	cl_double4*,
	cl_double*,
	cl_double*
);

void bdy_Uniform(
	sBdyUniformConfiguration*,
	cl_double2*,
	cl_double*,
	cl_double*,
	cl_double*,
	cl_double4*,
	cl_double*,
	cl_double*,
	GlobalHandlerClass
);

//#endif