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


#ifndef RANDOM_H_
#define RANDOM_H_

#include "p_utils.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1 -1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df		/* constant vector a */
#define UPPER_MASK 0x80000000	/* most significant w-r bits */
#define LOWER_MASK 0x7fffffff	/* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >18)

void sgenrand (unsigned long seed);
double uniformGenerator (void);
int getIntRandom (int maxValue);
int getWeightedIntRandom (double *weights);
double getMatrixMean (double **mat, double **P, int height, int width);
double StandardExponential (void);
double transformedExponential (double scale, double lag);
double StandardGaussian (void);
double transformedGaussian (double mean, double stddev);
boolean coinToss (double p);
void subsetSelect (double p, int seqlength, int *selected);
void bootstrapSelect (int seqlength, int *selected);

#endif /*RANDOM_H_*/

