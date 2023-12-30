/*
 * ------------------------------------------------------------------
 *  Author: Alaa Mroue - 2023
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "normalPlain.h"

normalPlain::normalPlain (int sizex,int sizey) {
	this->sizex = sizex;
	this->sizey = sizey;
	this->size = sizex * sizey;

	// Allocate memory for the data array
	this->bedElevation = new float* [sizex];
	for (int i = 0; i < sizex; i++) {
		this->bedElevation[i] = new float[sizey];
	}

}

int normalPlain::getSize() {
	return this->size;
}

double normalPlain::getBedElevation(int index) {
	return this->bedElevation[(int) floor(index/10)][index % 10];
}

float normalPlain::getBedElevation(int indexX, int indexY) {
	return this->bedElevation[indexX][indexY];
}

void normalPlain::setBedElevation(int indexX, int indexY, float value) {
	this->bedElevation[indexX][indexY] = value;
}

void normalPlain::setBedElevation(cl_double4* src) {
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			this->setBedElevation(i, j, src[i * 10 + j].s[0]);
		}
	}
}

void normalPlain::SetBedElevationMountain() {
	int const SIZE = this->getSize()/10;
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			this->setBedElevation(i, j, 0);
			this->setBedElevation(i, j, 0);
		}
	}
	this->setBedElevation(6, 6, 0.16);
	this->setBedElevation(6, 7, 0.16);
	this->setBedElevation(7, 6, 0.16);
	this->setBedElevation(7, 7, 0.16);
}

void normalPlain::outputShape() {
	int size = this->getSize()/10;
	double value;
	std::cout << std::fixed;
	std::cout << std::setprecision(2);
	std::cout << std::endl;

	for (int i = size-1; i > -1 ; i--) {
		for (int j = 0; j < size; j++) {
			value = this->getBedElevation(i, j);
			if (value-2 > 100+i*10+j) {
				std::cout <<  value << " ";
			}
			else if (value > 100 + i * 10 + j) {
				std::cout << value << " ";
			
			}
			else if (value < 100 + i * 10 + j) {
				std::cout << value << " ";
			}
			else {
				std::cout << value << " ";
			
			}

		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}