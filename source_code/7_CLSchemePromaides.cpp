/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "7_CLSchemePromaides.h"


#define Cgg 9.8066
#define Cfacweir 2.95245

void solverFunctionPromaides(
	cl_double* dTimestep,						// Timestep
	cl_double* dBedElevation,					// Bed elevation
	cl_double4* pCellStateSrc,					// Current cell state data
	cl_double4* pCellStateDst,					// Current cell state data
	cl_double* dManning,						// Manning values
	GlobalHandlerClass ghc
)
{


	// Identify the cell we're reconstructing (no overlap)
	cl_long					lIdxX = ghc.get_global_id(0);
	cl_long					lIdxY = ghc.get_global_id(1);
	cl_ulong					ulIdx, ulIdxNeigN, ulIdxNeigE, ulIdxNeigS, ulIdxNeigW;

	ulIdx = getCellID(lIdxX, lIdxY);

	// Don't bother if we've gone beyond the domain bounds
	if (lIdxX >= DOMAIN_COLS - 1 ||
		lIdxY >= DOMAIN_ROWS - 1 ||
		lIdxX <= 0 ||
		lIdxY <= 0) {
		//printf("Went Beyoooond {%i, %i} {%ld,%ld}\n",DOMAIN_COLS,DOMAIN_ROWS,lIdxX,lIdxY);
		return;
	}


	cl_double		dLclTimestep = *dTimestep;
	cl_double		dManningCoef;
	cl_double		dCellBedElev, dNeigBedElevN, dNeigBedElevE, dNeigBedElevS, dNeigBedElevW;
	cl_double		opt_h, opt_hN, opt_hE, opt_hS, opt_hW;
	cl_double		opt_cN, opt_cE, opt_cS, opt_cW;
	cl_double		opt_z, opt_zNmax, opt_zEmax, opt_zSmax, opt_zWmax;
	cl_double		opt_s, opt_sN, opt_sE, opt_sS, opt_sW;
	cl_double4	pCellData, pNeigDataN, pNeigDataE, pNeigDataS, pNeigDataW;					// Z, Zmax, Qx, Qy


	cl_uchar		ucStop = 0;
	cl_uchar		ucDryCount = 0;
	cl_double		ds_dt_data = 0.0;
	cl_double		flow_depth = 0.0;
	cl_double		flow_depth_neigh = 0.0;
	cl_double		delta_h = 0.0;
	cl_double		abs_delta_h = 0.0;
	cl_double		ds_dt_buff = 0.0;
	cl_double		reduction_term = 0.0;
	cl_double		v_x = 0.0;
	cl_double		v_y = 0.0;
	cl_char			flowStates = 0;

	//TODO: Alaa 
	//This would be passed as a variable
	flowStates |= (true << 0);
	flowStates |= (false << 1);
	flowStates |= (false << 2);
	flowStates |= (false << 3);
	flowStates |= (false << 4);

	bool isFlowElement = (flowStates >> 0) & 1;
	bool noflow_x = (flowStates >> 1) & 1;
	bool noflow_y = (flowStates >> 2) & 1;
	bool opt_pol_x = (flowStates >> 3) & 1;
	bool opt_pol_y = (flowStates >> 4) & 1;

	//if (lIdxX==DOMAIN_COLS-1) {
	//	noflow_x = true;
	//}
	//if (lIdxY==DOMAIN_ROWS-1){
	//	noflow_y = true;
	//}

	// Also don't bother if we've gone beyond the total simulation time
	if (dLclTimestep <= 0.0)
	{
		//printf("90/n");
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
		printf("104/n");
		pCellStateDst[ulIdx] = pCellData;
		return;
	}



	ulIdxNeigN = getNeighbourByIndices(lIdxX, lIdxY, DOMAIN_DIR_N);
	dNeigBedElevN = dBedElevation[ulIdxNeigN];
	pNeigDataN = pCellStateSrc[ulIdxNeigN];

	ulIdxNeigE = getNeighbourByIndices(lIdxX, lIdxY, DOMAIN_DIR_E);
	dNeigBedElevE = dBedElevation[ulIdxNeigE];
	pNeigDataE = pCellStateSrc[ulIdxNeigE];

	ulIdxNeigS = getNeighbourByIndices(lIdxX, lIdxY, DOMAIN_DIR_S);
	dNeigBedElevS = dBedElevation[ulIdxNeigS];
	pNeigDataS = pCellStateSrc[ulIdxNeigS];

	ulIdxNeigW = getNeighbourByIndices(lIdxX, lIdxY, DOMAIN_DIR_W);
	dNeigBedElevW = dBedElevation[ulIdxNeigW];
	pNeigDataW = pCellStateSrc[ulIdxNeigW];


	//All neighbours are dry? Don't bother calculating
	if (pCellData.x - dCellBedElev < VERY_SMALL) ucDryCount++;
	if (pNeigDataN.x - dNeigBedElevN < VERY_SMALL) ucDryCount++;
	if (pNeigDataE.x - dNeigBedElevE < VERY_SMALL) ucDryCount++;
	if (pNeigDataW.x - dNeigBedElevW < VERY_SMALL) ucDryCount++;
	if (pNeigDataS.x - dNeigBedElevS < VERY_SMALL) ucDryCount++;
	if (ucDryCount == 5) {
		return;
	}
	//else{
	//	printf("116\n");
	//	
	//}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	//set data
	opt_z = dCellBedElev;
	if ((pCellData.x - opt_z) < VERY_SMALL) {
		opt_h = 0.0;
		opt_s = opt_z;
	}
	else {
		opt_h = pCellData.x - opt_z;
		opt_s = pCellData.x;
	}


	//Add Boundry and couple conditions
	//ds_dt_data = bound_cond_dsdt + coup_cond_dsdt;
	ds_dt_data = 10.0 / dLclTimestep * 0;

	//opt_h = pCellData.x  - dCellBedElev;
	opt_hN = pNeigDataN.x - dNeigBedElevN;
	opt_hE = pNeigDataE.x - dNeigBedElevE;
	opt_hS = pNeigDataS.x - dNeigBedElevS;
	opt_hW = pNeigDataW.x - dNeigBedElevW;

	opt_sN = pNeigDataN.x;
	opt_sE = pNeigDataE.x;
	opt_sS = pNeigDataS.x;
	opt_sW = pNeigDataW.x;


	opt_zNmax = dCellBedElev > dNeigBedElevN ? dCellBedElev : dNeigBedElevN;
	opt_zEmax = dCellBedElev > dNeigBedElevE ? dCellBedElev : dNeigBedElevE;
	opt_zSmax = dCellBedElev > dNeigBedElevS ? dCellBedElev : dNeigBedElevS;
	opt_zWmax = dCellBedElev > dNeigBedElevW ? dCellBedElev : dNeigBedElevW;

	opt_cN = 0.5 * (1 / dManning[ulIdx] + 1 / dManning[ulIdxNeigN]);
	opt_cE = 0.5 * (1 / dManning[ulIdx] + 1 / dManning[ulIdxNeigE]);
	opt_cS = 0.5 * (1 / dManning[ulIdx] + 1 / dManning[ulIdxNeigS]);
	opt_cW = 0.5 * (1 / dManning[ulIdx] + 1 / dManning[ulIdxNeigW]);

	//v_x = pCellData.z;
	//v_y = pCellData.w;

	if (!isFlowElement)
		return;

	//in x-direction
	if (!noflow_x) {
		if (!opt_pol_x) {
			//manning x
			if (opt_h > VERY_SMALL || opt_hE > VERY_SMALL) {

				if (flow_depth < 0.0) {
					flow_depth = 0.0;
				}
				flow_depth_neigh = opt_sE - opt_zEmax;
				if (flow_depth_neigh < 0.0) {
					flow_depth_neigh = 0.0;
				}

				flow_depth = std::max(flow_depth, flow_depth_neigh);

				ds_dt_buff = 0;
				if (flow_depth > VERY_SMALL) {
					//diffusive wave
					delta_h = opt_sE - opt_s;
					abs_delta_h = fabs(delta_h);

					if (abs_delta_h > VERY_SMALL) {
						ds_dt_buff = opt_cE * pow(flow_depth, (5.0 / 3.0));

						if (abs_delta_h <= 0.005078) {
							ds_dt_buff = ds_dt_buff * 0.10449968880528 * atan(159.877741951379 * delta_h); //0.0152
						}
						else {
							ds_dt_buff = ds_dt_buff * (delta_h / pow(abs_delta_h, 0.5));
						}

						//set the result
						ds_dt_data += ds_dt_buff;
						//printf("ds_dt_buff: %f \n" ,ds_dt_buff);
						v_x += -1 * ds_dt_buff * DOMAIN_DELTAX / flow_depth;
					}
				}
			}
		}
		else {
			printf("polini erorr");
			flow_depth = opt_s - opt_zEmax;
			flow_depth_neigh = opt_sE - opt_zEmax;

			//noFlow
			if ((flow_depth <= 0.0 && flow_depth_neigh <= 0.0) || (fabs(flow_depth - flow_depth_neigh) <= 0.0)) {
				ds_dt_buff = 0.0;
			}
			else { //Flow
				//flow out of this element without submerged weirflow reduction into the neihgbouring element
				if (flow_depth > 0.0 && flow_depth_neigh <= 0.0) {
					ds_dt_buff = -1.0 * Cfacweir * opt_cE * pow(flow_depth, (3.0 / 2.0));

					v_x = -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth;
				}
				//flow out of the neighbouring element without submerged weirflow reduction into this element
				else if (flow_depth <= 0.0 && flow_depth_neigh > 0.0) {

					ds_dt_buff = Cfacweir * opt_cE * pow(flow_depth_neigh, (3.0 / 2.0));
					v_x = -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth_neigh;
				}
				//submerged weirflow with reduction
				else if (flow_depth > 0.0 && flow_depth_neigh > 0.0) {
					//flow out of this element into the neihgbouring element
					if (flow_depth > flow_depth_neigh) {
						ds_dt_buff = Cfacweir * opt_cE * pow(flow_depth, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth_neigh / flow_depth);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = -1.0 * ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = -1.0 * ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}

						v_x = -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth;
					}
					//flow out of the neighbouring element into this element
					else {
						ds_dt_buff = Cfacweir * opt_cE * pow(flow_depth_neigh, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth / flow_depth_neigh);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}

						v_x = -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth_neigh;
					}
				}
				//set the result
				ds_dt_data = ds_dt_data + ds_dt_buff;
			}
		}
	}
	//in -x-direction
	if (!noflow_x) {
		if (!opt_pol_x) {
			//manning x
			if (opt_h > VERY_SMALL || opt_hW > VERY_SMALL) {

				flow_depth = opt_s - opt_zWmax;
				if (flow_depth < 0.0) {
					flow_depth = 0.0;
				}
				flow_depth_neigh = opt_sW - opt_zWmax;
				if (flow_depth_neigh < 0.0) {
					flow_depth_neigh = 0.0;
				}

				flow_depth = std::max(flow_depth, flow_depth_neigh);

				ds_dt_buff = 0;
				if (flow_depth > VERY_SMALL) {
					//diffusive wave
					delta_h = opt_sW - opt_s;
					abs_delta_h = fabs(delta_h);

					if (abs_delta_h > VERY_SMALL) {
						ds_dt_buff = opt_cW * pow(flow_depth, (5.0 / 3.0));

						if (abs_delta_h <= 0.005078) {
							ds_dt_buff = ds_dt_buff * 0.10449968880528 * atan(159.877741951379 * delta_h); //0.0152
						}
						else {
							ds_dt_buff = ds_dt_buff * (delta_h / pow(abs_delta_h, 0.5));
						}

						//set the result
						ds_dt_data += ds_dt_buff;
						v_x += ds_dt_buff * DOMAIN_DELTAX / flow_depth;
					}
				}
			}
		}
		else {
			printf("polini erorr");
			flow_depth = opt_s - opt_zWmax;
			flow_depth_neigh = opt_sW - opt_zWmax;

			//noFlow
			if ((flow_depth <= 0.0 && flow_depth_neigh <= 0.0) || (fabs(flow_depth - flow_depth_neigh) <= 0.0)) {
				ds_dt_buff = 0.0;
			}
			else { //Flow
				//flow out of this element without submerged weirflow reduction into the neihgbouring element
				if (flow_depth > 0.0 && flow_depth_neigh <= 0.0) {
					ds_dt_buff = -1.0 * Cfacweir * opt_cW * pow(flow_depth, (3.0 / 2.0));

					v_x += -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth;
				}
				//flow out of the neighbouring element without submerged weirflow reduction into this element
				else if (flow_depth <= 0.0 && flow_depth_neigh > 0.0) {

					ds_dt_buff = Cfacweir * opt_cW * pow(flow_depth_neigh, (3.0 / 2.0));
					v_x += -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth_neigh;
				}
				//submerged weirflow with reduction
				else if (flow_depth > 0.0 && flow_depth_neigh > 0.0) {
					//flow out of this element into the neihgbouring element
					if (flow_depth > flow_depth_neigh) {
						ds_dt_buff = Cfacweir * opt_cW * pow(flow_depth, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth_neigh / flow_depth);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = -1.0 * ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = -1.0 * ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}

						v_x += -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth;
					}
					//flow out of the neighbouring element into this element
					else {
						ds_dt_buff = Cfacweir * opt_cW * pow(flow_depth_neigh, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth / flow_depth_neigh);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}

						v_x += -1.0 * ds_dt_buff * DOMAIN_DELTAX / flow_depth_neigh;
					}
				}
				//set the result
				ds_dt_data = ds_dt_data + ds_dt_buff;
			}
		}
	}

	//in y-direction
	if (!noflow_y) {
		if (!opt_pol_y) {

			if (opt_h > VERY_SMALL || opt_hN > VERY_SMALL) {
				//calculate the mid of the flow depth
				flow_depth = opt_s - opt_zNmax;
				if (flow_depth < 0.0) {
					flow_depth = 0.0;
				}
				flow_depth_neigh = opt_sN - opt_zNmax;
				if (flow_depth_neigh < 0.0) {
					flow_depth_neigh = 0.0;
				}
				//flow_depth=(flow_depth+flow_depth_neigh)*0.5;
				flow_depth = std::max(flow_depth, flow_depth_neigh);

				ds_dt_buff = 0;
				if (flow_depth > VERY_SMALL) {
					//diffusive wave
					delta_h = opt_sN - opt_s;
					abs_delta_h = fabs(delta_h);

					if (abs_delta_h > VERY_SMALL) {
						ds_dt_buff = opt_cN * pow(flow_depth, (5.0 / 3.0));

						if (abs_delta_h <= 0.005078) {
							ds_dt_buff *= 0.10449968880528 * atan(159.877741951379 * delta_h); //0.0152
						}
						else {
							ds_dt_buff *= delta_h / pow(abs_delta_h, 0.5);
						}

						//set the result
						ds_dt_data += ds_dt_buff;
						v_y += -1 * ds_dt_buff * DOMAIN_DELTAY / flow_depth;
					}
				}
			}
		}
		else {
			printf("polini erorr");
			flow_depth = opt_s - opt_zNmax;
			flow_depth_neigh = opt_sN - opt_zNmax;

			//noFlow
			if ((flow_depth <= 0.0 && flow_depth_neigh <= 0.0) || (fabs(flow_depth - flow_depth_neigh) <= 0.0)) {
				ds_dt_buff = 0.0;
			}
			//flow
			else {
				//flow out of this element without submerged weirflow reduction into the neihgbouring element
				if (flow_depth > 0.0 && flow_depth_neigh <= 0.0) {
					ds_dt_buff = -1.0 * Cfacweir * opt_cN * pow(flow_depth, (3.0 / 2.0));
					v_y = -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth;
				}
				//flow out of the neighbouring element without submerged weirflow reduction into this element
				else if (flow_depth <= 0.0 && flow_depth_neigh > 0.0) {
					ds_dt_buff = Cfacweir * opt_cN * pow(flow_depth_neigh, (3.0 / 2.0));
					v_y = -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth_neigh;
				}
				//submerged weirflow with reduction
				else if (flow_depth > 0.0 && flow_depth_neigh > 0.0) {
					//flow out of this element into the neihgbouring element
					if (flow_depth > flow_depth_neigh) {
						ds_dt_buff = Cfacweir * opt_cN * pow(flow_depth, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth_neigh / flow_depth);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = -1.0 * ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = -1.0 * ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}
						v_y = -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth;
					}
					//flow out of the neighbouring element into this element
					else {
						ds_dt_buff = Cfacweir * opt_cN * pow(flow_depth_neigh, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth / flow_depth_neigh);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}
						v_y = -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth_neigh;
					}
				}

				//set result data
				ds_dt_data = ds_dt_data + ds_dt_buff;
			}
		}
	}

	//in -y-direction
	if (!noflow_y) {
		if (!opt_pol_y) {

			if (opt_h > VERY_SMALL || opt_hS > VERY_SMALL) {
				//calculate the mid of the flow depth
				flow_depth = opt_s - opt_zSmax;
				if (flow_depth < 0.0) {
					flow_depth = 0.0;
				}
				flow_depth_neigh = opt_sS - opt_zSmax;
				if (flow_depth_neigh < 0.0) {
					flow_depth_neigh = 0.0;
				}
				//flow_depth=(flow_depth+flow_depth_neigh)*0.5;
				flow_depth = std::max(flow_depth, flow_depth_neigh);

				ds_dt_buff = 0;
				if (flow_depth > VERY_SMALL) {
					//diffusive wave
					delta_h = opt_sS - opt_s;
					abs_delta_h = fabs(delta_h);

					if (abs_delta_h > VERY_SMALL) {
						ds_dt_buff = opt_cS * pow(flow_depth, (5.0 / 3.0));

						if (abs_delta_h <= 0.005078) {
							ds_dt_buff *= 0.10449968880528 * atan(159.877741951379 * delta_h); //0.0152
						}
						else {
							ds_dt_buff *= delta_h / pow(abs_delta_h, 0.5);
						}

						//set the result
						ds_dt_data += ds_dt_buff;
						v_y += ds_dt_buff * DOMAIN_DELTAY / flow_depth;
					}
				}
			}
		}
		else {
			printf("polini erorr");
			flow_depth = opt_s - opt_zSmax;
			flow_depth_neigh = opt_sS - opt_zSmax;

			//noFlow
			if ((flow_depth <= 0.0 && flow_depth_neigh <= 0.0) || (fabs(flow_depth - flow_depth_neigh) <= 0.0)) {
				ds_dt_buff = 0.0;
			}
			//flow
			else {
				//flow out of this element without submerged weirflow reduction into the neihgbouring element
				if (flow_depth > 0.0 && flow_depth_neigh <= 0.0) {
					ds_dt_buff = -1.0 * Cfacweir * opt_cS * pow(flow_depth, (3.0 / 2.0));
					v_y += -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth;
				}
				//flow out of the neighbouring element without submerged weirflow reduction into this element
				else if (flow_depth <= 0.0 && flow_depth_neigh > 0.0) {
					ds_dt_buff = Cfacweir * opt_cS * pow(flow_depth_neigh, (3.0 / 2.0));
					v_y += -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth_neigh;
				}
				//submerged weirflow with reduction
				else if (flow_depth > 0.0 && flow_depth_neigh > 0.0) {
					//flow out of this element into the neihgbouring element
					if (flow_depth > flow_depth_neigh) {
						ds_dt_buff = Cfacweir * opt_cS * pow(flow_depth, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth_neigh / flow_depth);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = -1.0 * ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = -1.0 * ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}
						v_y += -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth;
					}
					//flow out of the neighbouring element into this element
					else {
						ds_dt_buff = Cfacweir * opt_cS * pow(flow_depth_neigh, (3.0 / 2.0));
						//reduction of the discharge (submerged weirflow)
						reduction_term = (1.0 - flow_depth / flow_depth_neigh);
						//replace the ^(1/3) by a fitted arctan-function; at the boundary they have the same values
						if (reduction_term <= 0.000463529) {
							ds_dt_buff = ds_dt_buff * 0.057965266895 * atan(8984.365582471040 * reduction_term);
						}
						else {
							ds_dt_buff = ds_dt_buff * pow(reduction_term, (1.0 / 3.0));
						}
						v_y += -1.0 * ds_dt_buff * DOMAIN_DELTAY / flow_depth_neigh;
					}
				}

				//set result data
				ds_dt_data = ds_dt_data + ds_dt_buff;
			}
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////



	// Update the flow state
	pCellData.x = pCellData.x + dLclTimestep * ds_dt_data;
	//pCellData.z = v_x;
	//pCellData.w = v_y;

	//pNeigDataN.x		= pNeigDataN.x	- dLclTimestep * ds_dt_dataN;
	//pNeigDataE.x		= pNeigDataE.x	- dLclTimestep * ds_dt_dataE;

	// New max FSL?
	if (pCellData.x > pCellData.y && pCellData.y > -9990.0)
		pCellData.y = pCellData.x;

	// Crazy low depths?
	if (pCellData.x - dCellBedElev < VERY_SMALL)
		pCellData.x = dCellBedElev;

	if (lIdxX == 3 && lIdxY == 3)
	{
		//printf("pCellData.x: %f pCellData.y: %f pCellData.z: %f pCellData.w: %f   \n", pCellData.x, pCellData.y, pCellData.z, pCellData.w);
	}

	// Commit to global memory
	pCellStateDst[ulIdx] = pCellData;

}