/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "4_CLDynamicTimestep.h"

//Calculate the timestep using a reduction procedure and increment the total model time.

/*
 *  Advance the total model time by the timestep specified
 */
//__kernel  __attribute__((reqd_work_group_size(1, 1, 1)))
void tst_Advance_Normal(
	cl_double * dTime,
	cl_double * dTimestep,
	cl_double * dTimeHydrological,
	cl_double * pReductionData,
	cl_double4 * pCellData,
	cl_double * dBedData,
	cl_double * dTimeSync,
	cl_double * dBatchTimesteps,
	cl_uint * uiBatchSuccessful,
	cl_uint * uiBatchSkipped
)
{
	cl_double	dLclTime = *dTime;
	cl_double	dLclTimestep = fmax(0.0, *dTimestep);
	cl_double	dLclTimeHydrological = *dTimeHydrological;
	cl_double	dLclSyncTime = *dTimeSync;
	cl_double dLclBatchTimesteps = *dBatchTimesteps;
	cl_uint uiLclBatchSuccessful = *uiBatchSuccessful;
	cl_uint uiLclBatchSkipped = *uiBatchSkipped;

	// Increment total time (only ever referenced in this kernel)
	dLclTime += dLclTimestep;
	dLclBatchTimesteps += dLclTimestep;

	if (dLclTimestep > 0.0)
	{
		uiLclBatchSuccessful++;
	}
	else {
		uiLclBatchSkipped++;
	}

	// Hydrological processes run with their own timestep which is larger
	if (dLclTimeHydrological > TIMESTEP_HYDROLOGICAL)
	{
		dLclTimeHydrological = dLclTimestep;
	}
	else {
		dLclTimeHydrological += dLclTimestep;
	}

	#ifdef TIMESTEP_DYNAMIC

	cl_double dCellSpeed, dMaxSpeed;
	cl_double dMinTime;

	dCellSpeed = 0.0;
	dMaxSpeed = 0.0;
	for (unsigned int i = 0; i < TIMESTEP_WORKERS; i++)
	{
		dCellSpeed = pReductionData[i];
		if (dCellSpeed > dMaxSpeed)
			dMaxSpeed = dCellSpeed;
	}

	// Convert velocity to a time (assumes domain deltaX=deltaY here)
	// Force progression at the start of a simulation.
	dMinTime = DOMAIN_DELTAX / dMaxSpeed;
	if (dLclTime < TIMESTEP_START_MINIMUM_DURATION && dMinTime < TIMESTEP_START_MINIMUM)
		dMinTime = TIMESTEP_START_MINIMUM;

	// Multiply by the Courant number
	dLclTimestep = COURANT_NUMBER * dMinTime;

	#endif
	#ifdef TIMESTEP_FIXED

	dLclTimestep = TIMESTEP_FIXED;

	#endif

	// Don't exceed the output interval
	// but also don't stop things at the start
	// Also don't exceed the synchronisation time
	/*
	if ( fmod( dLclTime, SCHEME_OUTPUTTIME ) < 1E-7 && dLclTime > 0.5 )
	{
		dLclTimestep = 0.0;
	} else {
		if ( ( dLclTime + dLclTimestep ) > ( trunc( dLclTime / SCHEME_OUTPUTTIME ) + 1 ) * SCHEME_OUTPUTTIME )
			dLclTimestep = ( ( trunc( dLclTime / SCHEME_OUTPUTTIME ) + 1 ) * SCHEME_OUTPUTTIME ) - dLclTime;
	}
	*/

	// Impose a minimum timestep
	if (dLclTimestep > 0.0 && dLclTimestep < TIMESTEP_MINIMUM)
		dLclTimestep = TIMESTEP_MINIMUM;

	// Don't exceed the sync time
	// A negative timestep suspends simulation but allows the value to be used
	// back on the host.
	if ((dLclTime + dLclTimestep) >= dLclSyncTime)
	{
		if (dLclSyncTime - dLclTime > VERY_SMALL)
			dLclTimestep = dLclSyncTime - dLclTime;
		if (dLclSyncTime - dLclTime <= VERY_SMALL)
			dLclTimestep = -dLclTimestep;
	}

	// Control the timestep initially to ensure it's not silly, because
	// boundary conditions may only just be kicking in (i.e. dry domain)
	if (dLclTime < TIMESTEP_EARLY_LIMIT_DURATION && dLclTimestep > TIMESTEP_EARLY_LIMIT)
		dLclTimestep = TIMESTEP_EARLY_LIMIT;

	// Don't exceed the total simulation time
	if ((dLclTime + dLclTimestep) > SCHEME_ENDTIME)
		dLclTimestep = SCHEME_ENDTIME - dLclTime;

	// A sensible maximum timestep
	if (dLclTimestep > TIMESTEP_MAXIMUM)
		dLclTimestep = TIMESTEP_MAXIMUM;

	// Commit to global memory
	*dTime = dLclTime;
	*dTimestep = dLclTimestep;
	*dTimeHydrological = dLclTimeHydrological;
	*dBatchTimesteps = dLclBatchTimesteps;
	*uiBatchSuccessful = uiLclBatchSuccessful;
	*uiBatchSkipped = uiLclBatchSkipped;
}

/*
 *  Advance the total model time by the timestep specified
 */
//__kernel  __attribute__((reqd_work_group_size(1, 1, 1)))
void tst_ResetCounters(
	cl_double* dBatchTimesteps,
	cl_uint* uiBatchSuccessful,
	cl_uint* uiBatchSkipped
)
{
	*uiBatchSuccessful = 0;
	*uiBatchSkipped = 0;
	*dBatchTimesteps = 0.0;
}

/*
 *  Reduce the timestep by calculating for each workgroup
 */
//__kernel  REQD_WG_SIZE_LINE
void tst_Reduce(
	cl_double4* pCellData,
	cl_double* dBedData,
	cl_double* pReductionData,
	GlobalHandlerClass ghc
)
{
	cl_double pScratchData[TIMESTEP_GROUPSIZE];

	// Get global ID for cell
	cl_uint		uiLocalID = ghc.get_local_id(0);
	cl_uint		uiLocalSize = ghc.get_local_size(0);

	cl_ulong	ulCellID = ghc.get_global_id(0);
	cl_double4	pCellState;
	cl_double	dBedElevation;
	cl_double	dCellSpeed, dDepth, dVelX, dVelY;
	cl_double	dMaxSpeed = 0.0;

	while (ulCellID < DOMAIN_CELLCOUNT)
	{
		// Calculate the velocity...
		pCellState = pCellData[ulCellID];
		dBedElevation = dBedData[ulCellID];

		dDepth = pCellState.x - dBedElevation;

		if (dDepth > QUITE_SMALL && pCellState.y > -9999.0)
		{
			#ifndef TIMESTEP_SIMPLIFIED

			dVelX = pCellState.z / dDepth;
			dVelY = pCellState.w / dDepth;
			if (dVelX < 0.0) dVelX = -dVelX;
			if (dVelY < 0.0) dVelY = -dVelY;

			dVelX += sqrt(GRAVITY * dDepth);
			dVelY += sqrt(GRAVITY * dDepth);

			#else

			dVelX = sqrt(GRAVITY * dDepth);
			dVelY = sqrt(GRAVITY * dDepth);

			#endif
			dCellSpeed = (dVelX < dVelY) ? dVelY : dVelX;

		}
		else {
			dCellSpeed = 0.0;
		}

		// Is this velocity higher, therefore a greater time constraint?
		if (dCellSpeed > dMaxSpeed)
			dMaxSpeed = dCellSpeed;

		// Move on to the next cell
		ulCellID += ghc.get_global_size(0);
	}

	// Commit to local memory
	pScratchData[uiLocalID] = dMaxSpeed;

	//Alaa
	// No progression until scratch memory is fully populated
	//barrier(CLK_LOCAL_MEM_FENCE);

	// 2nd stage of the reduction process
	// Funnelling style operation from the center
	for (int iOffset = uiLocalSize / 2;
		iOffset > 0;
		iOffset = iOffset / 2)
	{
		if (uiLocalID < iOffset)
		{
			cl_double	dComparison = pScratchData[uiLocalID + iOffset];
			cl_double	dMine = pScratchData[uiLocalID];
			pScratchData[uiLocalID] = (dMine < dComparison) ? dComparison : dMine;
		}

		//Alaa
		//barrier(CLK_LOCAL_MEM_FENCE);
	}

	// Only one workgroup to update the time
	if (uiLocalID == 0)
		pReductionData[ghc.get_group_id(0)] = pScratchData[0];
}

/*
 *  Update the timestep after a synchronisation or rollback
 *  Reduction will have been carried out again first.
 */
//__kernel  __attribute__((reqd_work_group_size(1, 1, 1)))
void tst_UpdateTimestep(
	cl_double* dTime,
	cl_double* dTimestep,
	cl_double* pReductionData,
	cl_double* dTimeSync,
	cl_double* dBatchTimesteps
)
{
	cl_double	dLclTime = *dTime;
	cl_double	dLclOriginalTimestep = fabs(*dTimestep);
	cl_double	dLclSyncTime = *dTimeSync;
	cl_double dLclBatchTimesteps = *dBatchTimesteps;
	cl_double dLclTimestep;

	//#ifdef TIMESTEP_DYNAMIC

	cl_double dCellSpeed, dMaxSpeed;
	cl_double dMinTime;

	dCellSpeed = 0.0;
	dMaxSpeed = 0.0;
	for (unsigned int i = 0; i < TIMESTEP_WORKERS; i++)
	{
		dCellSpeed = pReductionData[i];
		if (dCellSpeed > dMaxSpeed)
			dMaxSpeed = dCellSpeed;
	}

	// Convert velocity to a time (assumes domain deltaX=deltaY here)
	// Force progression at the start of a simulation.
	dMinTime = DOMAIN_DELTAX / dMaxSpeed;

	if (dLclTime < TIMESTEP_START_MINIMUM_DURATION && dMinTime < TIMESTEP_START_MINIMUM)
		dMinTime = TIMESTEP_START_MINIMUM;

	// Multiply by the Courant number
	dLclTimestep = COURANT_NUMBER * dMinTime;

	//#endif

	// We only adjust the timestep if it's SMALLER than our original
	// which already accounted for sync points etc.
	// TODO: Force time advancing if necessary instead of calling this...
	dLclTimestep = fmin(dLclTimestep, dLclOriginalTimestep);
	dLclBatchTimesteps = dLclBatchTimesteps - dLclOriginalTimestep + dLclTimestep;

	// Don't exceed the early limit
	if (dLclTime < TIMESTEP_EARLY_LIMIT_DURATION && dLclTimestep > TIMESTEP_EARLY_LIMIT)
		dLclTimestep = TIMESTEP_EARLY_LIMIT;

	// Don't exceed the sync time
	if ((dLclTime + dLclTimestep) >= dLclSyncTime)
		dLclTimestep = fmax((cl_double)0.0, dLclSyncTime - dLclTime);

	// A sensible maximum timestep
	if (dLclTimestep > TIMESTEP_MAXIMUM)
		dLclTimestep = TIMESTEP_MAXIMUM;

	// Commit to global memory
	*dTimestep = dLclTimestep;
	*dBatchTimesteps = dLclBatchTimesteps;
}
