/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include "definitions.h"


cl_ulong	getNeighbourID(cl_ulong, cl_uchar);
cl_ulong	getNeighbourByIndices(cl_long, cl_long, cl_uchar);
cl_ulong	getCellID(cl_long, cl_long);
void		getCellIndices(cl_ulong, cl_long*, cl_long*);