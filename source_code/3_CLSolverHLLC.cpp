/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "3_CLSolverHLLC.h"

#define USE_ALTERNATE_CONSTRUCTS 1

// Calculate an approximate solution to the Riemann problem at the cell interface using the HLLC approach.

cl_double4 riemannSolver(
	cl_uchar	ucDirection,
	cl_double8		pLeft,
	cl_double8		pRight,
	bool		bDebug
)
{
	cl_uint2	uiDirectionVector;
	cl_double	FM_L, FM_R, F1_M, F2_M;
	cl_double	s_L, s_R, s_M, a_Avg, H_star, U_star, A_star;
	cl_double4	pFluxL, pFluxR, pFlux;
	cl_double2	dVel, dDis, dA;
	bool		bLeft, bRight, bMiddle_1, bMiddle_2;

	#ifdef USE_ALTERNATE_CONSTRUCTS

	cl_uint2 trail1 = { 0 , 1 };
	cl_uint2 trail2 = { 1 , 0 };

	uiDirectionVector = (ucDirection == DOMAIN_DIR_N || ucDirection == DOMAIN_DIR_S) ? trail1 : trail2;

	cl_double a1;
	cl_double a2;
	cl_double a3;
	cl_double a4;
	cl_double a5;
	cl_double a6;

	a1 = (pLeft.s[0] + pRight.s[0]) / 2;
	a2 = (pLeft.s[0] + pRight.s[0]) / 2;
	a3 = pLeft.s[6] * (pLeft.s[0] + pRight.s[0]);
	a4 = a2 - a3;
	a5 = (a1 * a2 - a3) * 0.5 * 9.81;
	a6 = 0.5 * GRAVITY * (
		
		((pLeft.s[0] + pRight.s[0]) / 2)     * (
													(pLeft.s[0] + pRight.s[0]) / 2) -
		pLeft.s[6]
													* (pLeft.s[0] + pRight.s[0])
		
		
		
		);

	// Are both sides dry? Simple solution if so...
	if (pLeft.s[1] < VERY_SMALL && pRight.s[1] < VERY_SMALL)
	{
		pFlux = {
				0.0,
				uiDirectionVector.s[0] * 0.5 * GRAVITY * (
				((pLeft.s[0] + pRight.s[0]) / 2) * ((pLeft.s[0] + pRight.s[0]) / 2) -
				pLeft.s[6] * (pLeft.s[0] + pRight.s[0])
				),
				uiDirectionVector.s[1] * 0.5 * GRAVITY * (
				((pLeft.s[0] + pRight.s[0]) / 2) * ((pLeft.s[0] + pRight.s[0]) / 2) -
				pLeft.s[6] * (pLeft.s[0] + pRight.s[0])
				),
				0.0
		};

		return pFlux;
	}
	#else
	uiDirectionVector = ((ucDirection == DOMAIN_DIR_N || ucDirection == DOMAIN_DIR_S) ? cl_uint2({ 0, 1 }) : cl_uint2({ 1, 0 }));

	// Are both sides dry? Simple solution if so...
	if (pLeft.s[1] < VERY_SMALL && pRight.s[1] < VERY_SMALL)
	{
		pFlux = {
			0.0,
			uiDirectionVector.s[0] * 0.5 * GRAVITY * (
				((pLeft.s[0] + pRight.s[0]) / 2) * ((pLeft.s[0] + pRight.s[0]) / 2) -
				pLeft.s[6] * (pLeft.s[0] + pRight.s[0])
				),
			uiDirectionVector.s[1] * 0.5 * GRAVITY * (
				((pLeft.s[0] + pRight.s[0]) / 2) * ((pLeft.s[0] + pRight.s[0]) / 2) -
				pLeft.s[6] * (pLeft.s[0] + pRight.s[0])
				),
			0.0
			};

		return pFlux;
	}
	#endif

	// Is one side dry?
	// -> Left
	pLeft.s[4] = (pLeft.s[1] < VERY_SMALL ? 0.0 : pLeft.s[2] / pLeft.s[1]);
	pLeft.s[5] = (pLeft.s[1] < VERY_SMALL ? 0.0 : pLeft.s[3] / pLeft.s[1]);

	// -> Right
	pRight.s[4] = (pRight.s[1] < VERY_SMALL ? 0.0 : pRight.s[2] / pRight.s[1]);
	pRight.s[5] = (pRight.s[1] < VERY_SMALL ? 0.0 : pRight.s[3] / pRight.s[1]);

	// Prerequisite calculations
	#ifdef USE_ALTERNATE_CONSTRUCTS
	dVel = {
			uiDirectionVector.s[0] * pLeft.s[4] + uiDirectionVector.s[1] * pLeft.s[5],				// Left
			uiDirectionVector.s[0] * pRight.s[4] + uiDirectionVector.s[1] * pRight.s[5]				// Right
	};
	dDis = {
			uiDirectionVector.s[0] * pLeft.s[2] + uiDirectionVector.s[1] * pLeft.s[3],				// Left
			uiDirectionVector.s[0] * pRight.s[2] + uiDirectionVector.s[1] * pRight.s[3]				// Right
	};
	dA = {
			sqrt(GRAVITY * pLeft.s[1]),														// Left
			sqrt(GRAVITY * pRight.s[1])														// Right
	};
	#else
	dVel = {
		uiDirectionVector.s[0] * pLeft.s[4] + uiDirectionVector.s[1] * pLeft.s[5],		// Left
		uiDirectionVector.s[0] * pRight.s[4] + uiDirectionVector.s[1] * pRight.s[5]		// Right
		};
	dDis = {
		uiDirectionVector.s[0] * pLeft.s[2] + uiDirectionVector.s[1] * pLeft.s[3],		// Left
		uiDirectionVector.s[0] * pRight.s[2] + uiDirectionVector.s[1] * pRight.s[3]		// Right
		};
	dA = {
		sqrt(GRAVITY * pLeft.s[1]),												// Left
		sqrt(GRAVITY * pRight.s[1])													// Right
		};
	#endif

	a_Avg = (dA.s[0] + dA.s[1]) / 2;
	H_star = ((a_Avg + (dVel.s[0] - dVel.s[1]) / 4) * (a_Avg + (dVel.s[0] - dVel.s[1]) / 4)) / GRAVITY;
	U_star = (dVel.s[0] + dVel.s[1]) / 2 + dA.s[0] - dA.s[1];
	A_star = sqrt(GRAVITY * H_star);

	// Calculate speed estimates
	if (pLeft.s[1] < VERY_SMALL)
	{
		s_L = dVel.s[1] - 2 * dA.s[1];
	}
	else {
		s_L = (((dVel.s[0] - dA.s[0]) > (U_star - A_star)) ? (U_star - A_star) : (dVel.s[0] - dA.s[0]));
	}
	if (pRight.s[1] < VERY_SMALL)
	{
		s_R = dVel.s[0] + 2 * dA.s[0];
	}
	else {
		s_R = (((dVel.s[1] + dA.s[1]) < (U_star + A_star)) ? (U_star + A_star) : (dVel.s[1] + dA.s[1]));
	}
	s_M = (s_L * pRight.s[1] * (dVel.s[1] - s_R) - s_R * pLeft.s[1] * (dVel.s[0] - s_L)) /
		(pRight.s[1] * (dVel.s[1] - s_R) - pLeft.s[1] * (dVel.s[0] - s_L));

	// Flux on left and right
	#ifdef USE_ALTERNATE_CONSTRUCTS
	pFluxL = {
						dDis.s[0],
						dVel.s[0] * pLeft.s[2] + uiDirectionVector.s[0] * 0.5 * GRAVITY * (pLeft.s[0] * pLeft.s[0] - 2 * pLeft.s[6] * pLeft.s[0]),
						dVel.s[0] * pLeft.s[3] + uiDirectionVector.s[1] * 0.5 * GRAVITY * (pLeft.s[0] * pLeft.s[0] - 2 * pLeft.s[6] * pLeft.s[0]),
						0.0
	};
	pFluxR = {
						dDis.s[1],
						dVel.s[1] * pRight.s[2] + uiDirectionVector.s[0] * 0.5 * GRAVITY * (pRight.s[0] * pRight.s[0] - 2 * pLeft.s[6] * pRight.s[0]),
						dVel.s[1] * pRight.s[3] + uiDirectionVector.s[1] * 0.5 * GRAVITY * (pRight.s[0] * pRight.s[0] - 2 * pLeft.s[6] * pRight.s[0]),
						0.0
	};
	#else
	pFluxL = {
		dDis.s[0],
		dVel.s[0] * pLeft.s[2] + uiDirectionVector.s[0] * 0.5 * GRAVITY * (pLeft.s[0] * pLeft.s[0] - 2 * pLeft.s[6] * pLeft.s[0]),
		dVel.s[0] * pLeft.s[3] + uiDirectionVector.s[1] * 0.5 * GRAVITY * (pLeft.s[0] * pLeft.s[0] - 2 * pLeft.s[6] * pLeft.s[0]),
		0.0
		};
	pFluxR = {
		dDis.s[1],
		dVel.s[1] * pRight.s[2] + uiDirectionVector.s[0] * 0.5 * GRAVITY * (pRight.s[0] * pRight.s[0] - 2 * pLeft.s[6] * pRight.s[0]),
		dVel.s[1] * pRight.s[3] + uiDirectionVector.s[1] * 0.5 * GRAVITY * (pRight.s[0] * pRight.s[0] - 2 * pLeft.s[6] * pRight.s[0]),
		0.0
		};
	#endif

	// Selection of the final result
	bLeft = s_L >= 0.0;
	bMiddle_1 = s_L < 0.0 && s_R >= 0.0 && s_M >= 0.0;
	bMiddle_2 = s_L < 0.0 && s_R >= 0.0 && !bMiddle_1;
	bRight = !bLeft && !bMiddle_1 && !bMiddle_2;

	if (bLeft)
	{
		#ifdef DEBUG_OUTPUT
		if (bDebug)
		{
			printf("(Dir %i) Using left fluxes\n", ucDirection);
		}
		#endif
		return pFluxL;
	}
	if (bRight)
	{
		#ifdef DEBUG_OUTPUT
		if (bDebug)
		{
			printf("(Dir %i) Using right fluxes\n", ucDirection);
		}
		#endif
		return pFluxR;
	}

	FM_L = uiDirectionVector.s[0] * pFluxL.y + uiDirectionVector.s[1] * pFluxL.z;
	FM_R = uiDirectionVector.s[0] * pFluxR.y + uiDirectionVector.s[1] * pFluxR.z;
	F1_M = (s_R * pFluxL.x - s_L * pFluxR.x + s_L * s_R * (pRight.s[0] - pLeft.s[0])) / (s_R - s_L);
	F2_M = (s_R * FM_L - s_L * FM_R + s_L * s_R * (dDis.s[1] - dDis.s[0])) / (s_R - s_L);

	#ifdef USE_ALTERNATE_CONSTRUCTS
	if (bMiddle_1)
	{
		pFlux = {
			F1_M,
			uiDirectionVector.s[0] * F2_M + uiDirectionVector.s[1] * F1_M * pLeft.s[4],
			uiDirectionVector.s[0] * F1_M * pLeft.s[5] + uiDirectionVector.s[1] * F2_M,
			0.0
		};
	}

	if (bMiddle_2)
	{
		pFlux = {
			F1_M,
			uiDirectionVector.s[0] * F2_M + uiDirectionVector.s[1] * F1_M * pRight.s[4],
			uiDirectionVector.s[0] * F1_M * pRight.s[5] + uiDirectionVector.s[1] * F2_M,
			0.0
		};
	}
	#else
	if (bMiddle_1)
	{
		pFlux = {
			F1_M,
			uiDirectionVector.s[0] * F2_M + uiDirectionVector.s[1] * F1_M * pLeft.s[4],
			uiDirectionVector.s[0] * F1_M * pLeft.s[5] + uiDirectionVector.s[1] * F2_M,
			0.0
			};
	}

	if (bMiddle_2)
	{
		pFlux = {
			F1_M,
			uiDirectionVector.s[0] * F2_M + uiDirectionVector.s[1] * F1_M * pRight.s[4],
			uiDirectionVector.s[0] * F1_M * pRight.s[5] + uiDirectionVector.s[1] * F2_M,
			0.0
			};
	}
	#endif

	return pFlux;
}

