/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "2_CLFriction.h"

//IMPLICIT FRICTION
// Reduce discharge values according to friction coefficients.


 //Point-implicit calculation of the friction effects
cl_double4 implicitFriction(
	cl_double4		pCellState,
	cl_double		dBedElevation,
	cl_double		dManningCoefficient,
	cl_double		dLclTimestep
)
{
	cl_double		dDepth, dQ;

	// Calculate depth and composite discharge
	dQ = sqrt(pCellState.z * pCellState.z + pCellState.w * pCellState.w);
	dDepth = pCellState.x - dBedElevation;

	// Low discharge or low depth means don't bother
	if (dDepth < VERY_SMALL || dQ < VERY_SMALL) return pCellState;

	// Coefficient of friction, etc. See Liang (2010)
	cl_double		dCf = (GRAVITY * dManningCoefficient * dManningCoefficient) / (pow((cl_double)dDepth, (cl_double)(1.0 / 3.0)));
	cl_double		dSfx = (-dCf / (dDepth * dDepth)) * pCellState.z * dQ;
	cl_double		dSfy = (-dCf / (dDepth * dDepth)) * pCellState.w * dQ;
	cl_double		dDx = 1.0 + dLclTimestep * (dCf / (dDepth * dDepth)) * (2 * (pCellState.z * pCellState.z) + (pCellState.w * pCellState.w)) / dQ;
	cl_double		dDy = 1.0 + dLclTimestep * (dCf / (dDepth * dDepth)) * ((pCellState.z * pCellState.z) + 2 * (pCellState.w * pCellState.w)) / dQ;
	cl_double		dFx = dSfx / dDx;
	cl_double		dFy = dSfy / dDy;

	// Friction can only stop flow, not reverse it
	if (pCellState.z >= 0.0)
	{
		if (dFx < -pCellState.z / dLclTimestep) dFx = -pCellState.z / dLclTimestep;
	}
	else {
		if (dFx > -pCellState.z / dLclTimestep) dFx = -pCellState.z / dLclTimestep;
	}
	if (pCellState.w >= 0.0)
	{
		if (dFy < -pCellState.w / dLclTimestep) dFy = -pCellState.w / dLclTimestep;
	}
	else {
		if (dFy > -pCellState.w / dLclTimestep) dFy = -pCellState.w / dLclTimestep;
	}

	// Update and commit data
	pCellState.z = pCellState.z + dLclTimestep * dFx;
	pCellState.w = pCellState.w + dLclTimestep * dFy;

	return pCellState;
}

/*
 *  Adjust the discharge with regard to friction
 */
void per_Friction(
	cl_double* dTimestep,
	cl_double4* pCellData,
	cl_double* dBedData,
	cl_double* dManningData,
	cl_double* dTime,			// TODO: Remove this, only required for temp rain
	GlobalHandlerClass ghc
)
{
	cl_double		dLclTimestep = *dTimestep;
	cl_long		lIdxX = ghc.get_global_id(0);
	cl_long		lIdxY = ghc.get_global_id(1);
	cl_ulong		ulIdx;

	cl_double4	pCellState;
	cl_double		dBedElevation;
	cl_double		dManningCoefficient;

	// Don't bother if we've gone beyond the domain bounds
	if (lIdxX >= DOMAIN_COLS - 1 || lIdxY >= DOMAIN_ROWS - 1 || lIdxX == 0 || lIdxY == 0)
		return;

	// Also don't bother if we've gone beyond the total simulation time
	if (dLclTimestep <= 0.0)
		return;

	ulIdx = getCellID(lIdxX, lIdxY);
	pCellState = pCellData[ulIdx];
	dBedElevation = dBedData[ulIdx];
	dManningCoefficient = dManningData[ulIdx];

	if (pCellState.x - dBedElevation < VERY_SMALL)
		return;

	pCellState = implicitFriction(
		pCellState,
		dBedElevation,
		dManningCoefficient,
		dLclTimestep
	);

	// TODO: Remvoe this
	// TEMPORARY ONLY
	// Introduce some rainfall to the domain at 10mm/hr
	//pCellState.x += 0.060/3600 * dLclTimestep;

	pCellData[ulIdx] = pCellState;
}
