/*
 * ------------------------------------------------------------------
 *  Author: Alaa Mroue - 2023
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "main.h"
#include <iostream>

using namespace std;

int main() {

	// Initializations
	unsigned long long iterationToPerform = 100;
	unsigned long long nextBatchIterations = iterationToPerform;
	cl_double dTimestep = 0.0001;
	cl_double pTimeHydrological = 0;
	cl_double pTime = 0;

	/* Testing riemannSolver
	cl_double8 left = d(158.227507, 0.000007, 0.000000, 0.000000, 0.000000, 0.000000, 158.227500, 0.0);
	cl_double8 right = d(158.227500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 158.227500, 0.0);
	
	rp(left, right);
	cout << endl;
	*/

	// Define a uniform grid with mountain-like terrain
	normalPlain np = normalPlain(10, 10);
	normalPlain np2 = normalPlain(10, 10);
	np.SetBedElevationMountain();

	// Define Boundary Conditions
	sBdyUniformConfiguration pConfiguration;
	pConfiguration.TimeseriesEntries = 2;
	pConfiguration.TimeseriesInterval = 1000;
	pConfiguration.TimeseriesLength = 1000000.00;
	pConfiguration.Definition = 0;
	
	// Define Time series 
	cl_double2* pTimeseries = new cl_double2[2];
	pTimeseries[0] = { 0,11.5 };
	pTimeseries[1] = { 360000,11.5 };

	// Define water levels
	cl_double* dBedElevation = new cl_double[100];
	cl_double4* pCellStateSrc = new cl_double4[100];
	cl_double4* pCellStateDst = new cl_double4[100];
	cl_double* dManning = new cl_double[100];
	for (int i = 0; i < 100; i++){ 
		dBedElevation[i] = np.getBedElevation(i);
		pCellStateSrc[i] = { np.getBedElevation(i) + 0.1,0,0,0 };
		pCellStateDst[i] = pCellStateSrc[i];
		dManning[i] = 100;
	
	}

	np.outputShape();
	
	//Main Program Loops

	while(iterationToPerform > 0 ){

		//Apply Rain
		GlobalHandlerClass ghc1 = GlobalHandlerClass(0, 0);
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				bdy_Uniform(&pConfiguration,pTimeseries, &pTime,&dTimestep,&pTimeHydrological,pCellStateSrc,dBedElevation,dManning,ghc1);
				ghc1.globalintY++;
			}
			ghc1.globalintY = 0;
			ghc1.globalintX++;
		}


		//Apply Scheme
		GlobalHandlerClass ghc2 = GlobalHandlerClass(0, 0);
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				gts_cacheDisabled(&dTimestep, dBedElevation, pCellStateSrc, pCellStateDst, dManning, ghc2);
				//solverFunctionPromaides(&dTimestep, dBedElevation, pCellStateSrc, pCellStateDst, dManning, ghc2);
				ghc2.globalintY++;
			}
			ghc2.globalintY = 0;
			ghc2.globalintX++;
		}

		//Advance Time
		pTime += dTimestep;
		pTimeHydrological += dTimestep;

		//Set Results
		pCellStateSrc = pCellStateDst;

		//Output progress
		if (iterationToPerform % 876 ==0)
			cout << "\rIteration Left:" << iterationToPerform << "      Time Spent: " << pTime << " s";
		iterationToPerform--;

		//Output Results
		if (iterationToPerform == 0) {
			np2.setBedElevation(pCellStateDst);
			cout << "\rAfter " << nextBatchIterations << " Iterations, Spent: " << pTime << " s                          " << endl;
			np2.outputShape();
			cout << "How many Iterations to perform?: ";
			cin >> nextBatchIterations;
			iterationToPerform = nextBatchIterations;
			cout << endl;
		}
	}

	return 0;
}

/*

Input of gts_cacheDisabled:
dTimestep			: Timestamp between Calling, higher=faster but possibility of no convergence
dBedElevation		: Elevation of the bed
pCellStateSrc		: Z, Zmax, Qx, Qy
pCellStateDst		: Z, Zmax, Qx, Qy
dManning			: Manning

Input of Riemann Solver
	ucDirection: Direction (4 choices)
	pLeft, pRight:
				0. Z   : Water level?
				1. H   : WaterDepth
				2. Qx  : Momentum in x direction
				3. Qy  : Momentum in y direction
				4. U   : is calculated: Velocity Components
				5. V   : is calculated: Velocity Components
				6. Zb  : Bed Elevation
				7. ??  : Not Used

Output of Reimann Solver
	Z : flux
	Qx: Momentum in x direction
	Qy: Momentum in y direction


bdy_Uniform:
sBdyUniformConfiguration:   {Total Number of entries, entry2-entry1, 0=Rain / 1=Loss}
pTimeseries					{0, 11.5}  {time, uniform_amount}
pTime;						Current Simulation Time	
pTimestep;					Timestep used (Changes based on complexity)
pTimeHydrological;			Current Simulation Time % 1.0


*/