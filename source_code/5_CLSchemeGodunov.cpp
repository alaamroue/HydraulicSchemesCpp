/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "5_CLSchemeGodunov.h"

//Implementation of the 1st order accurate Godunov-type scheme


 //Reconstruct the cell data in a non-negative way (depth positivity preserving)
cl_uchar reconstructInterface(
	cl_double4		pStateLeft,						// Left current state		Z, Zmax, Qx, Qy
	cl_double		dBedLeft,						// Left bed elevation
	cl_double4		pStateRight,					// Right current state
	cl_double		dBedRight,						// Right bed elevation
	cl_double8* pOutputLeft,						// Output data for LHS of Riemann
	cl_double8* pOutputRight,						// Output data for RHS of Riemann
	cl_uchar		ucDirection						// Direction under consideration
)
{
	cl_uchar		ucStop = 0;
	cl_double8		pReconstructionLeft, pReconstructionRight;
	cl_double		dDepthL = pStateLeft.x - dBedLeft;
	cl_double		dDepthR = pStateRight.x - dBedRight;

	// Initial values before reconstruction
	pReconstructionLeft = {
		pStateLeft.s[0],																		// Z	S0
			dDepthL,																			// H	S1
			pStateLeft.s[2],																	// Qx	S2
			pStateLeft.s[3],																	// Qy	S3
			(dDepthL < VERY_SMALL ? 0.0 : pStateLeft.s[2] / dDepthL),							// U	S4
			(dDepthL < VERY_SMALL ? 0.0 : pStateLeft.s[3] / dDepthL),							// V	S5
			dBedLeft,																			// Zb	S6
			0.0 };																			//		S7
	pReconstructionRight = {
		pStateRight.s[0],																		// Z	S0
			dDepthR,																			// H	S1
			pStateRight.s[2],																	// Qx	S2
			pStateRight.s[3],																	// Qy	S3
			(dDepthR < VERY_SMALL ? 0.0 : pStateRight.s[2] / dDepthR),							// U	S4
			(dDepthR < VERY_SMALL ? 0.0 : pStateRight.s[3] / dDepthR),							// V	S5
			dBedRight,																			// Zb	S6
			0.0};																				//		S7

	// Maximum bed elevation and vertical shift factor
	cl_double	dBedMaximum = (pReconstructionLeft.s[6] > pReconstructionRight.s[6] ? pReconstructionLeft.s[6] : pReconstructionRight.s[6]);
	cl_double	dShiftV = dBedMaximum - (ucDirection < DOMAIN_DIR_S ? pStateLeft : pStateRight).s[0];
	if (dShiftV < 0.0) dShiftV = 0.0;

	// Adjustment of depths and dependent elements
	pReconstructionLeft.s[1] = (pStateLeft.s[0] - dBedMaximum > 0.0 ? (pStateLeft.s[0] - dBedMaximum) : 0.0);
	pReconstructionLeft.s[0] = pReconstructionLeft.s[1] + dBedMaximum;
	pReconstructionLeft.s[2] = pReconstructionLeft.s[1] * pReconstructionLeft.s[4];
	pReconstructionLeft.s[3] = pReconstructionLeft.s[1] * pReconstructionLeft.s[5];

	pReconstructionRight.s[1] = (pStateRight.s[0] - dBedMaximum > 0.0 ? (pStateRight.s[0] - dBedMaximum) : 0.0);
	pReconstructionRight.s[0] = pReconstructionRight.s[1] + dBedMaximum;
	pReconstructionRight.s[2] = pReconstructionRight.s[1] * pReconstructionRight.s[4];
	pReconstructionRight.s[3] = pReconstructionRight.s[1] * pReconstructionRight.s[5];

	// Prevent draining from a dry cell
	// and the stopping conditions
	switch (ucDirection)
	{
	case DOMAIN_DIR_N:

		// NOTE: Do NOT include zero velocity in the check. Velocity must be non-negative for stopping conditions
		//		 to be required.
		if (pReconstructionLeft.s[1] <= VERY_SMALL && pStateLeft.w > 0.0) { ucStop++; }
		if (pReconstructionRight.s[1] <= VERY_SMALL && pReconstructionLeft.s[5] < 0.0) { ucStop++; pReconstructionLeft.s[5] = 0.0; }
		if (pReconstructionLeft.s[1] <= VERY_SMALL && pReconstructionRight.s[5] > 0.0) { ucStop++; pReconstructionRight.s[5] = 0.0; }

		break;
	case DOMAIN_DIR_S:

		if (pReconstructionRight.s[1] <= VERY_SMALL && pStateRight.w < 0.0) { ucStop++; }
		if (pReconstructionRight.s[1] <= VERY_SMALL && pReconstructionLeft.s[5] < 0.0) { ucStop++; pReconstructionLeft.s[5] = 0.0; }
		if (pReconstructionLeft.s[1] <= VERY_SMALL && pReconstructionRight.s[5] > 0.0) { ucStop++; pReconstructionRight.s[5] = 0.0; }

		break;
	case DOMAIN_DIR_E:

		if (pReconstructionLeft.s[1] <= VERY_SMALL && pStateLeft.z > 0.0) { ucStop++; }
		if (pReconstructionRight.s[1] <= VERY_SMALL && pReconstructionLeft.s[4] < 0.0) { ucStop++; pReconstructionLeft.s[4] = 0.0; }
		if (pReconstructionLeft.s[1] <= VERY_SMALL && pReconstructionRight.s[4] > 0.0) { ucStop++; pReconstructionRight.s[4] = 0.0; }

		break;
	case DOMAIN_DIR_W:

		if (pReconstructionRight.s[1] <= VERY_SMALL && pStateRight.z < 0.0) { ucStop++; }
		if (pReconstructionRight.s[1] <= VERY_SMALL && pReconstructionLeft.s[4] < 0.0) { ucStop++; pReconstructionLeft.s[4] = 0.0; }
		if (pReconstructionLeft.s[1] <= VERY_SMALL && pReconstructionRight.s[4] > 0.0) { ucStop++; pReconstructionRight.s[4] = 0.0; }

		break;
	}

	// Local modification of the bed level (and consequently, FSL to maintain depth)
	pReconstructionLeft.s[6] = dBedMaximum - dShiftV;
	pReconstructionRight.s[6] = dBedMaximum - dShiftV;
	pReconstructionLeft.s[0] -= dShiftV;
	pReconstructionRight.s[0] -= dShiftV;

	// Output vector: Z, H, Qx, Qy, U, V, Zb
	* pOutputLeft = pReconstructionLeft;
	*pOutputRight = pReconstructionRight;

	// Stop flow?
	return ucStop;
}

/*
 *  Calculate everything without using LDS caching
 */
//__kernel REQD_WG_SIZE_FULL_TS
void gts_cacheDisabled(
	cl_double* dTimestep,						// Timestep
	cl_double* dBedElevation,					// Bed elevation
	cl_double4* pCellStateSrc,					// Current cell state data
	cl_double4* pCellStateDst,					// Current cell state data
	cl_double* dManning,						// Manning values
	GlobalHandlerClass ghc
)
{

	//int gid = getCellID(5, 5);
	//printf("dBedElevation    :  %lf\n", dBedElevation[gid]);

	// Identify the cell we're reconstructing (no overlap)
	cl_long					lIdxX = ghc.get_global_id(0);
	cl_long					lIdxY = ghc.get_global_id(1);
	cl_ulong					ulIdx, ulIdxNeig;
	cl_uchar					ucDirection;

	ulIdx = getCellID(lIdxX, lIdxY);

	// Don't bother if we've gone beyond the domain bounds
	if (lIdxX >= DOMAIN_COLS - 1 ||
		lIdxY >= DOMAIN_ROWS - 1 ||
		lIdxX <= 0 ||
		lIdxY <= 0)
		//printf("Went Beyoooond {%i, %i} {%ld,%ld}\n",DOMAIN_COLS,DOMAIN_ROWS,lIdxX,lIdxY);
		return;

	cl_double		dLclTimestep = *dTimestep;
	cl_double		dManningCoef;
	cl_double		dCellBedElev, dNeigBedElevN, dNeigBedElevE, dNeigBedElevS, dNeigBedElevW;
	cl_double4	pCellData, pNeigDataN, pNeigDataE, pNeigDataS, pNeigDataW;					// Z, Zmax, Qx, Qy
	cl_double4	pSourceTerms, dDeltaValues;										// Z, Qx, Qy
	cl_double4	pFlux[4];																// Z, Qx, Qy
	cl_double8	pLeft, pRight;												// Z, H, Qx, Qy, U, V, Zb
	cl_uchar		ucStop = 0;
	cl_uchar		ucDryCount = 0;


	// Also don't bother if we've gone beyond the total simulation time
	if (dLclTimestep <= 0.0)
	{
		// TODO: Is there a way of avoiding this?!
		pCellStateDst[ulIdx] = pCellStateSrc[ulIdx];
		return;
	}

	// Load cell data
	dCellBedElev = dBedElevation[ulIdx];
	pCellData = pCellStateSrc[ulIdx];
	dManningCoef = dManning[ulIdx];

	// Cell disabled?
	if (pCellData.y <= -9999.0 || pCellData.x == -9999.0)
	{
		pCellStateDst[ulIdx] = pCellData;
		return;
	}

	ucDirection = DOMAIN_DIR_W;
	ulIdxNeig = getNeighbourByIndices(lIdxX, lIdxY, ucDirection);
	dNeigBedElevW = dBedElevation[ulIdxNeig];
	pNeigDataW = pCellStateSrc[ulIdxNeig];
	ucDirection = DOMAIN_DIR_S;
	ulIdxNeig = getNeighbourByIndices(lIdxX, lIdxY, ucDirection);
	dNeigBedElevS = dBedElevation[ulIdxNeig];
	pNeigDataS = pCellStateSrc[ulIdxNeig];
	ucDirection = DOMAIN_DIR_N;
	ulIdxNeig = getNeighbourByIndices(lIdxX, lIdxY, ucDirection);
	dNeigBedElevN = dBedElevation[ulIdxNeig];
	pNeigDataN = pCellStateSrc[ulIdxNeig];
	ucDirection = DOMAIN_DIR_E;
	ulIdxNeig = getNeighbourByIndices(lIdxX, lIdxY, ucDirection);
	dNeigBedElevE = dBedElevation[ulIdxNeig];
	pNeigDataE = pCellStateSrc[ulIdxNeig];

	#ifdef DEBUG_OUTPUT
	if (lIdxX == DEBUG_CELLX && lIdxY == DEBUG_CELLY)
	{
		//printf( "\n");
		//printf( "    ulIdx:  { %f )\n", ulIdx);
		//printf( "    ulIdxNeig:  { %f )\n", ulIdxNeig);
		//printf( "    dNeigBedElevW:  { %f )\n", dNeigBedElevW);
		//printf( "    pNeigDataW:  { %f )\n", pNeigDataW);

		//printf( "DOMAIN_CELLCOUNT:  { %f )\n", DOMAIN_CELLCOUNT);
		//printf( "DOMAIN_ROWS:  { %f )\n", DOMAIN_ROWS);
		//printf( "DOMAIN_COLS:  { %f )\n", DOMAIN_COLS);

		//printf( "GTS_DIM1:  { %f )\n", GTS_DIM1);
		//printf( "GTS_DIM2:  { %f )\n", GTS_DIM2);
		//printf( "get_group_id:  { %f )\n", get_group_id(0));


		// Print the value of dManning for the current work-item
		//printf("\nGlobalId: {%i,%i} \n", get_global_id(0),get_global_id(1));
		//printf("dTimestep    :  %lf\n", dTimestep[gid]);
		//printf("dBedElevation    :  %lf\n", dBedElevation[gid]);
		//printf("dManning     :  %lf \n",dManning[gid]);
		//printf("ulIdxNeig     :  %lld \n" , ulIdxNeig);

		//printf("pCellStateSrc = (%lf, %lf, %lf, %lf)\n", (pCellStateSrc)[ulIdxNeig].x, (pCellStateSrc)[ulIdxNeig].y, (pCellStateSrc)[ulIdxNeig].z, (pCellStateSrc)[ulIdxNeig]);
		//printf("pCellStateDst = (%lf, %lf, %lf, %lf)\n", (pCellStateDst)[ulIdxNeig], (pCellStateDst)[ulIdxNeig], (pCellStateDst)[ulIdxNeig], (pCellStateDst)[ulIdxNeig]);

		printf( "Current data:  { %f, %f, %f, %f )\n", pCellData.x, pCellData.y, pCellData.z, pCellData.w );
		printf( "Neighbour N:   { %f, %f, %f, %f )\n", pNeigDataN.x, dNeigBedElevN, pNeigDataN.z, pNeigDataN.w );
		printf( "Neighbour E:   { %f, %f, %f, %f )\n", pNeigDataE.x, dNeigBedElevE, pNeigDataE.z, pNeigDataE.w );
		printf( "Neighbour S:   { %f, %f, %f, %f )\n", pNeigDataS.x, dNeigBedElevS, pNeigDataS.z, pNeigDataS.w );
		printf( "Neighbour W:   { %f, %f, %f, %f )\n", pNeigDataW.x, dNeigBedElevW, pNeigDataW.z, pNeigDataW.w );

	}
	#endif

	if (pCellData.x - dCellBedElev < VERY_SMALL) ucDryCount++;
	if (pNeigDataN.x - dNeigBedElevN < VERY_SMALL) ucDryCount++;
	if (pNeigDataE.x - dNeigBedElevE < VERY_SMALL) ucDryCount++;
	if (pNeigDataS.x - dNeigBedElevS < VERY_SMALL) ucDryCount++;
	if (pNeigDataW.x - dNeigBedElevW < VERY_SMALL) ucDryCount++;

	// All neighbours are dry? Don't bother calculating
	if (ucDryCount >= 5) return;

	// Reconstruct interfaces
	// -> North
	ucStop += reconstructInterface(
		pCellData,							// Left cell data
		dCellBedElev,						// Left bed elevation
		pNeigDataN,							// Right cell data
		dNeigBedElevN,						// Right bed elevation
		&pLeft,								// Output for left
		&pRight,							// Output for right
		DOMAIN_DIR_N
	);
	pNeigDataN.x = pRight.s[0];
	dNeigBedElevN = pRight.s[6];
	#ifdef DEBUG_OUTPUT
	if (lIdxX == DEBUG_CELLX && lIdxY == DEBUG_CELLY)
	{
		printf( "Reconstruct NL:{ %f, %f, %f, %f )\n", pLeft.s[0], pLeft.s[6], pLeft.s[2], pLeft.s[3] );
		printf( "Reconstruct NR:{ %f, %f, %f, %f )\n", pRight.s[0], pRight.s[6], pRight.s[2], pRight.s[3] );
	}
	#endif
	pFlux[DOMAIN_DIR_N] = riemannSolver(DOMAIN_DIR_N, pLeft, pRight, false);

	// -> South
	ucStop += reconstructInterface(
		pNeigDataS,							// Left cell data
		dNeigBedElevS,						// Left bed elevation
		pCellData,							// Right cell data
		dCellBedElev,						// Right bed elevation
		&pLeft,								// Output for left
		&pRight,							// Output for right
		DOMAIN_DIR_S
	);
	pNeigDataS.x = pLeft.s[0];
	dNeigBedElevS = pLeft.s[6];
	pFlux[DOMAIN_DIR_S] = riemannSolver(DOMAIN_DIR_S, pLeft, pRight, false);

	// -> East
	ucStop += reconstructInterface(
		pCellData,							// Left cell data
		dCellBedElev,						// Left bed elevation
		pNeigDataE,							// Right cell data
		dNeigBedElevE,						// Right bed elevation
		&pLeft,								// Output for left
		&pRight,							// Output for right
		DOMAIN_DIR_E
	);
	pNeigDataE.x = pRight.s[0];
	dNeigBedElevE = pRight.s[6];
	pFlux[DOMAIN_DIR_E] = riemannSolver(DOMAIN_DIR_E, pLeft, pRight, false);

	// -> West
	ucStop += reconstructInterface(
		pNeigDataW,							// Left cell data
		dNeigBedElevW,						// Left bed elevation
		pCellData,							// Right cell data
		dCellBedElev,						// Right bed elevation
		&pLeft,								// Output for left
		&pRight,							// Output for right
		DOMAIN_DIR_W
	);
	pNeigDataW.x = pLeft.s[0];
	dNeigBedElevW = pLeft.s[6];
	pFlux[DOMAIN_DIR_W] = riemannSolver(DOMAIN_DIR_W, pLeft, pRight, false);

	// Source term vector
	// TODO: Somehow get these sorted too...
	pSourceTerms.x = 0.0;
	pSourceTerms.y = -1 * GRAVITY * ((pNeigDataE.x + pNeigDataW.x) / 2) * ((dNeigBedElevE - dNeigBedElevW) / DOMAIN_DELTAX);
	pSourceTerms.z = -1 * GRAVITY * ((pNeigDataN.x + pNeigDataS.x) / 2) * ((dNeigBedElevN - dNeigBedElevS) / DOMAIN_DELTAY);

	// Calculation of change values per timestep and spatial dimension
	dDeltaValues.x = (pFlux[1].x - pFlux[3].x) / DOMAIN_DELTAX +
		(pFlux[0].x - pFlux[2].x) / DOMAIN_DELTAY -
		pSourceTerms.x;
	dDeltaValues.z = (pFlux[1].y - pFlux[3].y) / DOMAIN_DELTAX +
		(pFlux[0].y - pFlux[2].y) / DOMAIN_DELTAY -
		pSourceTerms.y;
	dDeltaValues.w = (pFlux[1].z - pFlux[3].z) / DOMAIN_DELTAX +
		(pFlux[0].z - pFlux[2].z) / DOMAIN_DELTAY -
		pSourceTerms.z;

	// Round delta values to zero if small
	// TODO: Explore whether this can be rewritten as some form of clamp operation?
	if ((dDeltaValues.x > 0.0 && dDeltaValues.x < VERY_SMALL) ||
		(dDeltaValues.x < 0.0 && dDeltaValues.x > -VERY_SMALL))
		dDeltaValues.x = 0.0;
	if ((dDeltaValues.z > 0.0 && dDeltaValues.z < VERY_SMALL) ||
		(dDeltaValues.z < 0.0 && dDeltaValues.z > -VERY_SMALL))
		dDeltaValues.z = 0.0;
	if ((dDeltaValues.w > 0.0 && dDeltaValues.w < VERY_SMALL) ||
		(dDeltaValues.w < 0.0 && dDeltaValues.w > -VERY_SMALL))
		dDeltaValues.w = 0.0;

	// Stopping conditions
	if (ucStop > 0)
	{
		pCellData.z = 0.0;
		pCellData.w = 0.0;
	}

	// Update the flow state
	pCellData.x = pCellData.x - dLclTimestep * dDeltaValues.x;
	pCellData.z = pCellData.z - dLclTimestep * dDeltaValues.z;
	pCellData.w = pCellData.w - dLclTimestep * dDeltaValues.w;

	#ifdef FRICTION_ENABLED
	#ifdef FRICTION_IN_FLUX_KERNEL
	// Calculate the friction effects
	pCellData = implicitFriction(
		pCellData,
		dCellBedElev,
		dManningCoef,
		dLclTimestep
	);
	#endif
	#endif

	// New max FSL?
	if (pCellData.x > pCellData.y && pCellData.y > -9990.0)
		pCellData.y = pCellData.x;

	// Crazy low depths?
	if (pCellData.x - dCellBedElev < VERY_SMALL)
		pCellData.x = dCellBedElev;

	// Commit to global memory
	pCellStateDst[ulIdx] = pCellData;
}

void rp(cl_double8 d1, cl_double8 d2) {

	cl_double4 alaa = riemannSolver(
		DOMAIN_DIR_N,
		d1,
		d2,
		true
	);
	printf("Left: (%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n", d1.s0, d1.s1, d1.s2, d1.s3, d1.s4, d1.s5, d1.s6, d1.s7);
	printf("Right: (%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n", d2.s0, d2.s1, d2.s2, d2.s3, d2.s4, d2.s5, d2.s6, d2.s7);
	printf("Z, Qx, Qy: (%lf, %lf, %lf)\n", alaa.s0, alaa.s1, alaa.s2);
	std::cout << std::endl;
}

cl_double8 d(cl_double x1, cl_double x2,
	cl_double x3, cl_double x4,
	cl_double U, cl_double V, cl_double Zb, cl_double _) {

	return { x1, x2, x3, x4, U, V, Zb, _ };

}