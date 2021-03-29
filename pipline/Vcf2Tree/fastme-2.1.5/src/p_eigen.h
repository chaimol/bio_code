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

#ifndef P_EIGEN_H
#define P_EIGEN_H

#include "p_utils.h"

#define BASE		 2	// base of floating point arithmetic
#define DIGITS		40	// no. of digits to the base BASE in the fraction
#define MAXITER		30	// max2. no. of iterations to converge

#define pos(i,j,n)	((i)*(n)+(j))
#define csize(a)	(fabs(a.re)+fabs(a.im))


int Eigen (int job, double *A, int n, double *rr, double *ri,
	double *vr, double *vi, double *w);
void balance (double *mat, int n, int *low, int *hi, double *scale);
void unbalance (int n, double *vr, double *vi, int low, int hi, double *scale);
int realeig (int job, double *mat, int n,int low, int hi, double *valr,
	double *vali, double *vr, double *vi);
void elemhess (int job, double *mat, int n, int low, int hi,
	double *vr, double *vi, int *work);
int ludcmp (double **a, int n, double *d);

/* complex functions */

typedef struct
{
	double re, im;
} complex;

complex compl (double re, double im);
complex _conj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
complex cfactor (complex x, double a);
int cxtoy (complex *x, complex *y, int n);
int cmatby (complex *a, complex *b, complex *c, int n, int m, int k);
int cmatout (FILE * fout, complex *x, int n, int m);
int cmatinv (complex *x, int n, int m, double *space);


#endif /*P_EIGEN_H_*/

