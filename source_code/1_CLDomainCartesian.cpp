/*
 * ------------------------------------------------------------------
 *  Code adapted by Alaa Mroue from HiPIMS by Luke S. Smith and Qiuhua Liang
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "1_CLDomainCartesian.h"

//Management functions for a Cartesian domain.

 /*
  *  Fetch the ID for a cell using its X and Y indices
  */
cl_ulong	getCellID(cl_long lIdxX, cl_long lIdxY)
{
	cl_long	lCols = DOMAIN_COLS;
	return (lIdxY * lCols) + lIdxX;
}

/*
 *  Fetch the X and Y indices for a cell using its ID
 */
void	getCellIndices(cl_ulong ulID, cl_long* lIdxX, cl_long* lIdxY)
{
	*lIdxX = ulID % DOMAIN_COLS;
	*lIdxY = (ulID - *lIdxX) / DOMAIN_COLS;
}

/*
 *  Fetch the ID for a neighbouring cell in the domain
 */
cl_ulong	getNeighbourID(cl_ulong ulCellID, cl_uchar ucDirection)
{
	cl_long lIdxX = 0;
	cl_long lIdxY = 0;
	getCellIndices(ulCellID, &lIdxX, &lIdxY);

	switch (ucDirection)
	{
	case DOMAIN_DIR_N:
		++lIdxY;
		break;
	case DOMAIN_DIR_E:
		++lIdxX;
		break;
	case DOMAIN_DIR_S:
		--lIdxY;
		break;
	case DOMAIN_DIR_W:
		--lIdxX;
		break;
	}

	return getCellID(lIdxX, lIdxY);
}

/*
 *  Fetch the ID for a neighbouring cell in the domain
 */
cl_ulong	getNeighbourByIndices(cl_long lIdxX, cl_long lIdxY, cl_uchar ucDirection)
{
	switch (ucDirection)
	{
	case DOMAIN_DIR_N:
		++lIdxY;
		break;
	case DOMAIN_DIR_E:
		++lIdxX;
		break;
	case DOMAIN_DIR_S:
		--lIdxY;
		break;
	case DOMAIN_DIR_W:
		--lIdxX;
		break;
	}

	return getCellID(lIdxX, lIdxY);
}
