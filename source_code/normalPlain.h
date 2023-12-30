/*
 * ------------------------------------------------------------------
 *  Author: Alaa Mroue - 2023
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#pragma once
#include <iostream>
#include <iomanip>
#include <CL/cl.h>


class normalPlain {
public:
	int sizex;
	int sizey;
	int size;
	float** bedElevation;

	normalPlain(int, int);
	int getSize();
	double getBedElevation(int index);
	float getBedElevation(int, int);
	void setBedElevation(int, int, float);
	void setBedElevation(cl_double4* src);
	void SetBedElevationMountain();
	void outputShape();
};