/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "6_CLBoundaries.h"


//Management functions for a domain boundaries.

void bdy_Cell(
	sBdyCellConfiguration* pConfiguration,
	cl_ulong* pRelations,
	cl_double4* pTimeseries,
	cl_double* pTime,
	cl_double* pTimestep,
	cl_double* pTimeHydrological,
	cl_double4* pCellState,
	cl_double* pCellBed,
	cl_double* pCellManning,
	GlobalHandlerClass ghc

)
{
	cl_long				lRelationID = ghc.get_global_id(0);
	sBdyCellConfiguration pConfig = *pConfiguration;
	cl_double				dLocalTime = *pTime;
	cl_double				dLocalTimestep = *pTimestep;

	if (lRelationID >= pConfig.RelationCount || dLocalTime >= pConfig.TimeseriesLength || dLocalTimestep <= 0.0)
		return;

	cl_ulong				ulBaseTimestep = (cl_ulong)floor(dLocalTime / pConfig.TimeseriesInterval);
	cl_ulong				ulNextTimestep = ulBaseTimestep + 1;
	cl_ulong				ulCellID = pRelations[lRelationID];
	cl_double4			pCellData = pCellState[ulCellID];
	cl_double				dCellBed = pCellBed[ulCellID];
	cl_double4			pTSBase = pTimeseries[ulBaseTimestep];
	cl_double4			pTSNext = pTimeseries[ulNextTimestep];

	// Interpolate between timesteps
	cl_double4 pTSInterp = {
		pTSBase.s[0] + pTSNext.s[0] - pTSBase.s[0] * (fmod(dLocalTime, pConfig.TimeseriesInterval) / pConfig.TimeseriesInterval),
		pTSBase.s[1] + pTSNext.s[1] - pTSBase.s[1] * (fmod(dLocalTime, pConfig.TimeseriesInterval) / pConfig.TimeseriesInterval),
		pTSBase.s[2] + pTSNext.s[2] - pTSBase.s[2] * (fmod(dLocalTime, pConfig.TimeseriesInterval) / pConfig.TimeseriesInterval),
		pTSBase.s[3] + pTSNext.s[3] - pTSBase.s[3] * (fmod(dLocalTime, pConfig.TimeseriesInterval) / pConfig.TimeseriesInterval)
	};

	//cl_double4 pTSInterp = pTSBase + alaatemp * (fmod(dLocalTime, pConfig.TimeseriesInterval) / pConfig.TimeseriesInterval);

	// Apply depth/fsl
	if (pConfig.DefinitionDepth == BOUNDARY_DEPTH_IS_DEPTH)
	{
		#ifdef DEBUG_OUTPUT
		printf("Depth is fixed.\n");
		#endif
		pCellData.x = dCellBed + pTSInterp.y;			// Depth is fixed
	}
	else if (pConfig.DefinitionDepth == BOUNDARY_DEPTH_IS_FSL)
	{
		#ifdef DEBUG_OUTPUT
		printf("FSL is fixed.\n");
		#endif
		pCellData.x = fmax(dCellBed, pTSInterp.y);		// FSL is fixed
	}
	else
	{
		#ifdef DEBUG_OUTPUT
		printf("Depth and FSL are free.\n");
		#endif
		if (fabs(pTSInterp.z) > VERY_SMALL ||
			fabs(pTSInterp.w) > VERY_SMALL ||
			pConfig.DefinitionDischarge == BOUNDARY_DISCHARGE_IS_VOLUME)
		{
			// Calculate a suitable depth based
			cl_double dDepth = (fabs(pTSInterp.z) * dLocalTimestep) / DOMAIN_DELTAY + (fabs(pTSInterp.w) * dLocalTimestep) / DOMAIN_DELTAX;
			cl_double dNormalDepth = fmax(pow(pTSInterp.z, 2) / GRAVITY, pow(pTSInterp.w, 2) / GRAVITY);
			cl_double dCriticalDepth = fmax(pow(pow(pTSInterp.z, 2) / GRAVITY, 1.0 / 3.0), pow(pow(pTSInterp.w, 2) / GRAVITY, 1.0 / 3.0));

			// Not going to impose a direction if we're trying to represent
			// a surging discharge rate (e.g. manhole surge)
			if (pConfig.DefinitionDischarge == BOUNDARY_DISCHARGE_IS_VOLUME)
			{
				// In the case of volume boundaries, no scaling has taken place
				dNormalDepth = 0.0;
				dDepth = (fabs(pTSInterp.z) * dLocalTimestep) / (DOMAIN_DELTAX * DOMAIN_DELTAY);
				dCriticalDepth = 0.0;
				pTSInterp.z = 0.0;
				pTSInterp.w = 0.0;
			}

			pCellData.x = fmax(dCellBed + dCriticalDepth, pCellData.x + dDepth);

			#ifdef DEBUG_OUTPUT
			printf("Setting depth as %f.\n", dDepth);
			#endif
		}
	}

	if (pConfig.DefinitionDischarge == BOUNDARY_DISCHARGE_IS_DISCHARGE)
	{
		// Apply flow in X direction
		pCellData.z = pTSInterp.z;
	}
	else if (pConfig.DefinitionDischarge == BOUNDARY_DISCHARGE_IS_VELOCITY) {
		// Apply velocity in X direction
		pCellData.z = pTSInterp.z * (pCellData.x - dCellBed);
	}

	if (pConfig.DefinitionDischarge == BOUNDARY_DISCHARGE_IS_DISCHARGE)
	{
		// Apply flow in Y direction
		pCellData.w = pTSInterp.w;
	}
	else if (pConfig.DefinitionDischarge == BOUNDARY_DISCHARGE_IS_VELOCITY) {
		// Apply velocity in X direction
		pCellData.w = pTSInterp.w * (pCellData.x - dCellBed);
	}

	#ifdef DEBUG_OUTPUT
	printf("Final Cell Data:       { %f, %f, %f, %f }\n", pCellData.x, pCellData.y, pCellData.z, pCellData.w);
	#endif

	pCellState[ulCellID] = pCellData;
}

void bdy_Uniform(
	sBdyUniformConfiguration* pConfiguration,
	cl_double2* pTimeseries,
	cl_double* pTime,
	cl_double* pTimestep,
	cl_double* pTimeHydrological,
	cl_double4* pCellState,
	cl_double* pCellBed,
	cl_double* pCellManning,
	GlobalHandlerClass ghc
)
{
	// Which global series are we processing, and which cell
	// Global ID is X, Y cell, then Z for the series
	cl_long		lIdxX = ghc.get_global_id(0);
	cl_long		lIdxY = ghc.get_global_id(1);
	cl_ulong		ulIdx;

	// Don't bother if we've gone beyond the domain bounds
	if (lIdxX >= DOMAIN_COLS - 1 ||
		lIdxY >= DOMAIN_ROWS - 1 ||
		lIdxX <= 0 ||
		lIdxY <= 0)
		return;

	ulIdx = getCellID(lIdxX, lIdxY);

	// How far in to the simulation are we? And current cell data
	sBdyUniformConfiguration	pConfig = *pConfiguration;
	cl_double4				pCellData = pCellState[ulIdx];
	cl_double					dCellBedElev = pCellBed[ulIdx];
	cl_double					dLclTime = *pTime;
	cl_double					dLclRealTimestep = *pTimestep;
	cl_double					dLclTimestep = *pTimeHydrological;

	// Hydrological processes have their own timesteps
	if (dLclTimestep < TIMESTEP_HYDROLOGICAL || dLclRealTimestep <= 0.0)
		return;

	if (dLclTime >= pConfig.TimeseriesLength || pCellData.y <= -9999.0)
		return;

	// Calculate the right cell and stuff to be grabbing data from here...
	cl_ulong ulTimestep = (cl_ulong)floor(dLclTime / pConfig.TimeseriesInterval);
	cl_double2 dRecord = pTimeseries[ulTimestep];

	// Apply the value...
	if (pConfig.Definition == BOUNDARY_UNIFORM_RAIN_INTENSITY)
		pCellData.x += dRecord.y / 3600000.0 * dLclTimestep;

	if (pConfig.Definition == BOUNDARY_UNIFORM_LOSS_RATE)
		pCellData.x = std::max(dCellBedElev, pCellData.x - dRecord.y / 3600000.0 * dLclTimestep);

	// Return to global memory
	pCellState[ulIdx] = pCellData;
}

void bdy_Gridded(
	sBdyGriddedConfiguration* pConfiguration,
	cl_double* pTimeseries,
	cl_double* pTime,
	cl_double* pTimestep,
	cl_double* pTimeHydrological,
	cl_double4* pCellState,
	cl_double* pCellBed,
	cl_double* pCellManning,
	GlobalHandlerClass ghc
)
{
	// Which global series are we processing, and which cell
	// Global ID is X, Y cell, then Z for the series
	cl_long		lIdxX = ghc.get_global_id(0);
	cl_long		lIdxY = ghc.get_global_id(1);
	cl_ulong		ulIdx;

	// Don't bother if we've gone beyond the domain bounds
	if (lIdxX >= DOMAIN_COLS - 1 ||
		lIdxY >= DOMAIN_ROWS - 1 ||
		lIdxX <= 0 ||
		lIdxY <= 0)
		return;

	ulIdx = getCellID(lIdxX, lIdxY);

	// How far in to the simulation are we? And current cell data
	sBdyGriddedConfiguration	pConfig = *pConfiguration;
	cl_double4				pCellData = pCellState[ulIdx];
	cl_double					dCellBedElev = pCellBed[ulIdx];
	cl_double					dLclTime = *pTime;
	cl_double					dLclTimestep = *pTimeHydrological;

	// Cell disabled?
	if (pCellData.y <= -9999.0 || pCellData.x == -9999.0)
		return;

	// Hydrological processes have their own timesteps
	if (dLclTimestep < TIMESTEP_HYDROLOGICAL)
		return;

	// Calculate the right cell and stuff to be grabbing data from here...
	cl_ulong ulTimestep = (cl_ulong)floor(dLclTime / pConfig.TimeseriesInterval);
	if (ulTimestep >= pConfig.TimeseriesEntries) ulTimestep = pConfig.TimeseriesEntries;

	cl_double ulColumn = floor((((cl_double)lIdxX * (cl_double)DOMAIN_DELTAX) - pConfig.GridOffsetX) / pConfig.GridResolution);
	cl_double ulRow = floor((((cl_double)lIdxY * (cl_double)DOMAIN_DELTAY) - pConfig.GridOffsetY) / pConfig.GridResolution);
	cl_ulong ulBdyCell = (pConfig.GridRows * pConfig.GridCols) * ulTimestep +
		(pConfig.GridCols * (cl_ulong)ulRow) + (cl_ulong)ulColumn;
	cl_double dRate = pTimeseries[ulBdyCell];

	// Apply the value...
	if (pConfig.Definition == BOUNDARY_GRIDDED_RAIN_INTENSITY)
		pCellData.x += dRate / 3600000.0 * dLclTimestep;

	if (pConfig.Definition == BOUNDARY_GRIDDED_MASS_FLUX)
		pCellData.x += dRate / ((cl_double)DOMAIN_DELTAX * (cl_double)DOMAIN_DELTAY) * dLclTimestep;

	// Return to global memory
	pCellState[ulIdx] = pCellData;
}
