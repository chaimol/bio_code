//
//  Copyright 2002-2007 Rick Desper, Olivier Gascuel
//  Copyright 2007-2014 Olivier Gascuel, Stephane Guindon, Vincent Lefort
//
//  This file is part of FastME.
//
//  FastME is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FastME is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastME.  If not, see <http://www.gnu.org/licenses/>
//


#include "random.h"

static unsigned long mt[N];	/* the array for the state vector */
static int mti=N+1;			/* mti==N+1 means mt[N] is not initialized */

/*********************************************************/

/* Initializing the array with a seed */
void sgenrand (unsigned long seed)
{
	int i;

	for (i=0; i<N; i++)
	{
		mt[i] = seed & 0xffff0000;
		seed = 69069 * seed + 1;
		mt[i] |= (seed & 0xffff0000) >> 16;
		seed = 69069 * seed + 1;
	}

	mti = N;

	return;
}

/*********************************************************/

/* Generating reals
 * unsigned long for integer generation */
double uniformGenerator ()
{
	int kk;
	unsigned long y;
	double temp;
	static unsigned long mag01[2]={0x0, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N)	/* generate N words at one time */
	{
		if (mti == N+1)		/* if sgenrand() has not been called, */
			sgenrand(4357);	/* a default initial seed is used */

		for (kk=0; kk<N-M; kk++)
		{
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (; kk<N-1; kk++)
		{
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
		}

		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

		mti = 0;
	}

	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

	temp = ( (double)y * 2.3283064365386963e-10 );	/* reals: [0,1)-interval */

	return (temp);	/* return y; for integer generation */
}

/*********************************************************/

/* random integer between 0 and maxValue-1, uniform */
int getIntRandom (int maxValue)
{
	double temp;
	int returnValue;

	temp = maxValue * uniformGenerator ();
	returnValue = (int) temp;

	return (returnValue);
}

/*********************************************************/

int getWeightedIntRandom (double *weights)
{
	double temp, partialSum;
	int returnValue = 0;

	partialSum = weights[0];
	temp = uniformGenerator();
	while (temp < partialSum)
	{
		returnValue++;
		partialSum += weights[returnValue];
	}
	returnValue--;

	return (returnValue);
}

/*********************************************************/

double getMatrixMean (double **mat, double **P, int height, int width)
{
	int i,j;
	double mu = 0.0;

	for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			mu += mat[i][j]*P[i][j];

	return (mu);
}

/*********************************************************/

/* Returns exponentially distributed, positive, random deviate of unit
 * mean, using ran(idum) as the source of uniform deviates */
double StandardExponential ()
{
	double dum;

	do
		dum = uniformGenerator ();
	while ( fabs (dum) < DBL_EPSILON );

	return (-log(dum));
}

/*********************************************************/

/* Returns  scale *exponential + lag, where exponential is a random
 * sample from an exponential distribution */
double transformedExponential (double scale, double lag)
{
	double value;

	value = StandardExponential () * scale + lag;

	return (value);
}

/*********************************************************/

/* Returns a random sample from a standard normal distribution */
double StandardGaussian ()
{
	static int iset = 0;
	static double gset;
	double fac, rsq,v1,v2;

	if (0 == iset)	/* no random number ready to go */
	{
		do	/* pick two random numbers in the square from -1 to + 1 */
		{
			v1 = 2.0 * uniformGenerator () - 1.0;
			v2 = 2.0 * uniformGenerator () - 1.0;
			rsq = (v1*v1) + (v2 * v2);	/* check if they are in the unit circle */
		} while (rsq >= 1.0 || rsq == 0.0);

		fac = sqrt (-2.0 * log (rsq) / rsq);
		gset = v1 * fac;
		iset = 1;

		return (v2 * fac);
	}
	else
	{
		iset = 0;
		return (gset);
	}
}

/*********************************************************/

/* Returns a random sample from a normal distribution of
 * specified mean and standard deviation */
double transformedGaussian (double mean, double stddev)
{
	double stdGaussian, returnValue;

	stdGaussian = StandardGaussian ();
	returnValue = stdGaussian * stddev + mean;

	return (returnValue);
}

/*********************************************************/

boolean coinToss (double p)
{
	double x;

	x = uniformGenerator ();

	if (x < p)
		return (TRUE);

	else
		return (FALSE);
}

/*********************************************************/

void subsetSelect (double p, int seqlength, int *selected)
{
	int i;

	for (i=0; i<seqlength; i++)
		selected[i] = coinToss (p) ? 1:0;

	return;
}

/*********************************************************/

void bootstrapSelect (int seqlength, int *selected)
{
	int i,j;

	for (i=0; i<seqlength; i++)
	{
		j = getIntRandom (seqlength);
		(selected[j])++;
	}

	return;
}
