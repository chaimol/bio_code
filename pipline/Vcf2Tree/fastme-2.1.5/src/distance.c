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


#include "distance.h"


/*********************************************************/

/* Returns the number of c1 -> c2 changes in s->t */
int countStateChanges (char *s, char *t, int length, char c1, char c2,
	int *filter)
{
	int i;
	int matches = 0;

	for (i=0; i<length; i++)
		if ((c1 == s[i]) && (c2 == t[i]))
			matches += filter[i];

	return (matches);
}

/*********************************************************/

int *copyFilter (int *filter, int l)
{
	int i;
	int *ret;
	
	ret = (int *) mCalloc (l, sizeof (int));
	
	for (i=0; i<l; i++)
		ret[i] = filter[i];
	
	return ret;
}

/*********************************************************/

void ijFilter (int *filter, char *s1, char *s2, int itype, int seqlength)
{
	int i;

	for (i=0; i<seqlength; i++)
	{
		if (!isBoostrap)
			filter[i] = 1;
		
		if (PROTEIN == itype)
		{
			if ((NULL == strchr (PROTEIN_ALPHABET, s1[i])) ||
				(NULL == strchr (PROTEIN_ALPHABET, s2[i])))
			{
				if (isBoostrap)
				{
					if (filter[i] > 0)
						filter[i] = filter[i] -1;
				}
				else
				{
					filter[i] = 0;
					if (verbose > 2)
						Debug ( (char*)"Removing site %d.", i);
				}

			}
			else if (('*' == s1[i]) || ('?' == s1[i]) || ('-' == s1[i]) ||
					('*' == s2[i]) || ('?' == s2[i]) || ('-' == s2[i]))
			{
				if (isBoostrap)
				{
					if (filter[i] > 0)
						filter[i] = filter[i] -1;
				}
				else
				{
					filter[i] = 0;
					if (verbose > 2)
						Debug ( (char*)"Removing site %d.", i);
				}

			}
		}
		else
		{
			if ((NULL == strchr (DNA_ALPHABET, s1[i])) ||
				(NULL == strchr (DNA_ALPHABET, s2[i])))
			{
				if (isBoostrap)
				{
					if (filter[i] > 0)
						filter[i] = filter[i] -1;
				}
				else
				{
					filter[i] = 0;
					if (verbose > 2)
						Debug ( (char*)"Removing site %d.", i);
				}

			}
			else if (('*' == s1[i]) || ('?' == s1[i]) || ('-' == s1[i]) ||
					('*' == s2[i]) || ('?' == s2[i]) || ('-' == s2[i]))
			{
				if (isBoostrap)
				{
					if (filter[i] > 0)
						filter[i] = filter[i] -1;
				}
				else
				{
					filter[i] = 0;
					if (verbose > 2)
						Debug ( (char*)"Removing site %d.", i);
				}
			}
		}
	}

	return;
}

/*********************************************************/

int seqCharMatches (char *s, int length, char c, int *filter)
{
	int i;
	int matches = 0;

	for (i=0; i<length; i++)
		if (c == s[i])
			matches += filter[i];

	return (matches);
}

/*********************************************************/

/* called when calculating stationary probabilities */
int matrixCharMatches (char **s, int numSeqs, int length, char c, int *filter)
{
	int i;
	int matches = 0;

	for (i=0; i<numSeqs; i++)
		matches += seqCharMatches (s[i], length, c, filter);

	return (matches);
}

/*********************************************************/

double *calcStationaryProbsGlobal (char **s, int numSeqs, int length,
	int *filter, int numSelected, int alphabetSize, const char *alphabet)
{
	int i;
	double *p;

	p = (double *) mCalloc (alphabetSize, sizeof (double));
	for (i=0; i<alphabetSize; i++)
	{
		p[i] = (double) (matrixCharMatches (s, numSeqs, length,
			alphabet[i], filter)) / (numSeqs * numSelected);
	}

	return (p);
}

/*********************************************************/

void calcTransitionProbs (double **P, char *s1, char *s2, int length,
	int *filter, int numSelected, int alphabetSize, const char *alphabet)
{
	int i, j;

	for(i=0; i<alphabetSize; i++)
		for (j=0; j<alphabetSize; j++)
			P[i][j] = (double) ( countStateChanges (s1, s2, length,
				alphabet[i], alphabet[j], filter) ) / (double) numSelected;

	// note: P is NOT a probability matrix for this program
	// sum_j P[i][j] = Pi[i]

	return;
}

/* should check formulas for following two functions
 * transition rate is probability of witnessing a transition change from
 * sequence i to sequence j. This rate should be equal to the sum of the
 * four probabilities below */

/*********************************************************/

double calcTransitionRate (double **P)
{
	double a;

	a = P[ADENINE][GUANINE] + P[GUANINE][ADENINE]
		+ P[CYTOSINE][THYMINE] + P[THYMINE][CYTOSINE];

	return (a);
}

/* The transversion rate is the probability of seeing a non-transition
 * change from sequence i to sequence j.
 * Essentially is the probability of a change minus the transitition probability */

/*********************************************************/

double calcTransversionRate (double **P)
{
	double b=0.0;

	b += P[ADENINE][CYTOSINE];
	b += P[ADENINE][THYMINE];
	b += P[GUANINE][CYTOSINE];
	b += P[GUANINE][THYMINE];
	b += P[CYTOSINE][ADENINE];
	b += P[CYTOSINE][GUANINE];
	b += P[THYMINE][GUANINE];
	b += P[THYMINE][ADENINE];

	return (b);
}

/*********************************************************/

double calcRYSYM (double b, boolean use_gamma, float gamma)	/* RY symetric distance */
{
	double returnValue, loc;

	if ( fabs (b) - DBL_EPSILON < 0 )
		return (0.0);

	loc = 1 - ( 2 * b );

	if (0 >= loc)
		return (DNA_DIST_MAX);

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = gamma * (0.5 * (pow (loc, -1.0 / gamma) - 1.0)) ;
	}
	else
	{
		returnValue = -0.5 * (log (loc));
	}

	return (returnValue);
}

/*********************************************************/

double calcJC69 (double b, boolean use_gamma, float gamma)	/* Jukes-Cantor distance */
{
	double returnValue, loc;

	if ( fabs (b) - DBL_EPSILON < 0 )
		return (0.0);

	loc = 1 - ( 4 * b / 3 );

	if (0 >= loc)
		return (DNA_DIST_MAX);

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = gamma * (0.75 * (pow (loc, -1.0 / gamma) - 1.0)) ;
	}
	else
	{
		returnValue = -0.75 * (log (loc));
	}

	return (returnValue);
}

/*********************************************************/

double calcF81 (double loc, double b, boolean use_gamma, float gamma)	/* Felsenstein 1981 distance */
{
	double returnValue;

	if ( fabs (b) - DBL_EPSILON < 0 )
		return (0.0);

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = gamma * loc * (pow ( (1.0 - (b / loc)), -1.0 / gamma) - 1) ;
	}
	else
	{
		returnValue = -1.0 * loc * (log (1.0 - (b / loc)));
	}

	return (returnValue);
}

/*********************************************************/

double calcK2P (double a, double b, boolean use_gamma, float gamma)	/* Kimura 2-parameter distance */
{
	double returnValue, loc1, loc2;

	if ( (fabs (a) - DBL_EPSILON < 0) && (fabs (b) - DBL_EPSILON < 0) )
		return (0.0);

	loc1 = 1.0 - 2.0*a - b;
	loc2 = 1.0 - 2.0*b;

	if ((0 >= loc1) || (0 >= loc2))
	{
		return (DNA_DIST_MAX);
	}

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = gamma * (0.5 * (pow (loc1, -1.0 / gamma)) +
					0.25 * (pow (loc2, -1.0 / gamma)) - 0.75);
	}
	else
	{
		returnValue = -0.5 * (log (loc1)) - 0.25 * (log (loc2));
	}

	return (returnValue);
}

/*********************************************************/

void calcF84AuxProbs (double *Pi, double *Mloc, double *Nloc, double *Ploc)
{
	double PAG, PCT, PR, PY;

	PAG = Pi[ADENINE]  * Pi[GUANINE];
	PCT = Pi[CYTOSINE] * Pi[THYMINE];
	PR  = Pi[ADENINE]  + Pi[GUANINE];
	PY  = Pi[CYTOSINE] + Pi[THYMINE];

	*Mloc = PAG/PR + PCT/PY;
	*Nloc = PAG + PCT;
	*Ploc = PR * PY;

	return;
}

/* Pi stationary frequencies P transition probabilities */

/*********************************************************/

double calcF84 (double a, double b, boolean use_gamma, float gamma, double Mloc, double Nloc, double Ploc)
{
	double returnValue;
	double loc1, loc2;

	if ( (fabs (a) - DBL_EPSILON < 0) && (fabs (b) - DBL_EPSILON < 0) )
		return (0.0);

	loc1 = 1.0 - a/(2.0*Mloc) - b*(Mloc - Nloc)/(2.0*Mloc*Ploc);
	loc2 = 1.0 - b/(2.0*Ploc);

	if ((0 >= loc1) || (0 >= loc2))
	{
		return (DNA_DIST_MAX);
	}

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = 2.0 * gamma * (Mloc * pow (loc1, -1.0 / gamma) +
					(Nloc + Ploc - Mloc) * pow (loc2, -1.0 / gamma) -
					Nloc - Ploc);
	}
	else
	{
		returnValue = -2.0 * Mloc * log (loc1) -
					 2.0 * (Nloc + Ploc - Mloc) * log (loc2);
	}

	return (returnValue);
}

/*********************************************************/

void calcTNAuxProbs (double *Pi, double *PAPG, double *PCPT, double *PR,
	double *PY)
{
	*PR   = Pi[ADENINE]  + Pi[GUANINE];
	*PY   = Pi[CYTOSINE] + Pi[THYMINE];
	*PAPG = Pi[ADENINE]  * Pi[GUANINE];
	*PCPT = Pi[CYTOSINE] * Pi[THYMINE];

	return;
}

/*********************************************************/

double calcTN93 (double aR, double aY, double b, double PR, double PY,
	double PAPG, double PCPT, boolean use_gamma, float gamma)
{
	double loc1, loc2, loc3;
	double returnValue;

	if ( (fabs (aR) - DBL_EPSILON < 0) && (fabs (aY) - DBL_EPSILON < 0) && (fabs (b) - DBL_EPSILON < 0) )
		return (0.0);

	loc1 = 1.0 - ( (PR * aR) / (2.0 * PAPG) ) - (b / (2.0 * PR));
	loc2 = 1.0 - ( (PY * aY) / (2.0 * PCPT) ) - (b / (2.0 * PY));
	loc3 = 1.0 - (b / (2.0 * PR * PY));

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = (2.0 * gamma * PAPG / PR * pow (loc1, -1.0 / gamma) ) +
					(2.0 * gamma * PCPT / PY * pow (loc2, -1.0 / gamma) ) +
					(2.0 * gamma * ( (PR * PY) - (PAPG * PY / PR) - (PCPT * PR / PY) ) * pow (loc3, -1.0 / gamma) ) -
					(2.0 * gamma * (PAPG + PCPT + (PR * PY)) );
	}
	else
	{
		returnValue = (-2.0 * PAPG / PR * log (loc1) ) -
					(2.0 * PCPT / PY * log (loc2) ) -
					(2.0 * ( (PR * PY) - (PAPG * PY / PR) - (PCPT * PR / PY) ) * log (loc3) );
	}

//DEBUG
//Debug ("fTsR %lf, fTsY %lf, fTr %lf, t1 %lf, t2 %lf, t3 %lf, d %lf", aR, aY, b, loc1, loc2, loc3, returnValue);

	return (returnValue);
}

/* transversion only distance */

/*********************************************************/

double XX (double PR, double PY, double b, boolean use_gamma, float gamma)
{
	double returnValue, Z;

	Z = 1 - PR * PR - PY * PY;

	if (b/Z >= 1)
	{
		return (DNA_DIST_MAX);
	}

	if (use_gamma)
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		returnValue = gamma * Z * (pow (1 - b / Z, -1 / gamma) -1);
	}
	else
	{
		returnValue = -Z * log (1- b / Z);
	}

	return (returnValue);
}

/* called only with positive integers */

/*********************************************************/

int factorial (int n)
{
	if (1 == n)
		return (n);

	else
		return (n * factorial (n-1));
}

/* index refers to index of the permutation in the (I think) Gray code
 * length is number of permtuations considered
 * index runs from 1 to factorial(length) */

/*********************************************************/

int *nextPerm (int *p, int index, int size, int length)
{
	int temp;
	if (0 ==  index % factorial (size))
	{
		temp = p[length-1];
		p[length-1] = p[length-1-size];
		p[length-1-size] = temp;
		return (p);
	}
	else
		return (nextPerm (p, index, size-1, length));
}

/*********************************************************/

double permDiagProduct (double **P, int *p, int d)
{
	int i;
	double prod = 1.0;

	for (i=0; i<d; i++)
		prod = prod * P[i][p[i]];

	return (prod);
}

/*********************************************************/

double det (double **P, int d)
{
	int *p;
	int signum = 1;
	int i, numPerms;
	double det = 0;

	p = initPerm (d);
	numPerms = factorial (d);

	for (i=0; i<numPerms; i++)
	{
		det += signum * permDiagProduct (P, p, d);
		p = nextPerm (p, i+1, d-1, d);
		signum = -1 * signum;
	}

	free (p);

	return(det);
}

/*********************************************************/

double logdet (double **P, double *Pi1, double *Pi2)
{
	int i;
	double returnValue, detP;

	detP = det (P, DNA_ALPHABET_SIZE);

	if (0 >= detP)
	{
		return (DNA_DIST_MAX);
	}

	returnValue = -0.5 * log (detP);

	for (i=0; i<DNA_ALPHABET_SIZE; i++)
	{
		if ((0 >= Pi1[i]) || (0 >= Pi2[i]))
			Exit ( (char*)"Logdet value of Pi1[i] is %f, of Pi2[i] is %f, i is %d.",
				Pi1[i], Pi2[i], i);

		returnValue += (log(Pi1[i]) + log(Pi2[i])) / 8;
	}

	return (returnValue);
}

/*********************************************************/

int support (int *v, int length)
{
	int i;
	int count = 0;

	for (i=0; i<length; i++)
		if (v[i])
			count++;

	return (count);
}

/*********************************************************/

double HammingDistance (char *v1, char *v2, int *filter, int length, int numSelected)
{
	int i;
	int d=0;

	for (i=0; i<length; i++)
		if (v1[i] != v2[i])
			d += filter[i];

	return ((double) d / numSelected);
}

/*********************************************************/

double protDiff (double *P)
{
	int i;
	double sum = 0.0;

	for (i=0; i<PROTEIN_ALPHABET_SIZE; i++)
		sum += P[i] * P[i];

	return (1 - sum);
}

/*********************************************************/

double protFormula (double b, boolean use_gamma, float gamma, double pdiff)
{
	double d, y;

	y = 1 - b / pdiff;

	if (y <= 0)
		return (PROT_DIST_MAX);

	if (! use_gamma)
		d = -1.0 * pdiff * log (y);

	else
	{
		gamma = (gamma < DBL_EPSILON) ? DBL_EPSILON : gamma;
		d = gamma * pdiff * (-1 + pow (y, -1 / gamma));
	}

	return (d);
}

/*********************************************************/

int aaIndex (char s, char *alphabet, int d)
{
	int i;

	for (i=0; i<d; i++)
		if (s == alphabet[i])
			return (i);

	Exit ( (char*)"Looking for character %c in protein string.", s);

	//return 0;
}

/*********************************************************/

double simScore (char *s1, char *s2, double **scoreMatrix, int seqlength,
	char *alphabet, int alphabetSize)
{
	int i;
	double sum = 0.0;

	for (i=0; i<seqlength; i++)
		sum += scoreMatrix [aaIndex (s1[i], alphabet, alphabetSize)]
						[aaIndex (s2[i], alphabet, alphabetSize)];

	return (sum);
}

/*********************************************************/

double expectedProtSimScore (double *P, double **scoreMatrix, int alphabetSize)
{
	int i, j;
	double sum = 0;

	for (i=0; i<alphabetSize; i++)
		for (j=0; j<alphabetSize; j++)
			sum += P[i] * P[j] * scoreMatrix[i][j];

	return (sum);
}

/*********************************************************/

double scoreDistij (int i, int j, char *si, char *sj, int seqLength,
	double simExp, double *simDiags, double **scoreMatrix, char *alphabet,
	int alphabetSize)
{
	double simij, simUpper, ratio;

	simij = simScore (si, sj, scoreMatrix, seqLength, alphabet, alphabetSize);
	simUpper = 0.5 * (simDiags[i] + simDiags[j]);
	ratio = (simij - simExp) / (simUpper - simExp);

	if (ratio <= 0 )
		return (DNA_DIST_MAX);
	else
		return (-100.0 * log (ratio));
}

/*********************************************************/

void scoreDist (double *P, char **data, int numSpecies, int seqLength,
	double **scoreMatrix, double **D, char *alphabet, int alphabetSize)
{
	int i, j;
	double simExp;
	double *simDiags;

	simDiags = (double *) mCalloc (numSpecies, sizeof (double));
	memset (simDiags, 0, (unsigned long) numSpecies * sizeof(double) );
	simExp = seqLength * expectedProtSimScore (P, scoreMatrix, alphabetSize);
	
#ifdef _OPENMP
	#pragma omp parallel
	{
		#pragma omp for
#endif

	for (i=0; i<numSpecies; i++)
	{
		D[i] = (double *) mCalloc (numSpecies, sizeof (double));
		D[i][i] = 0.0;
		simDiags[i] = simScore (data[i], data[i], scoreMatrix, seqLength,
				alphabet, alphabetSize);
	}
	
#ifdef _OPENMP
		#pragma omp for private (i, j)
#endif

	for (i=0; i<numSpecies-1; i++)
	{
		for (j=i+1; j<numSpecies; j++)
			D[i][j] = D[j][i] = scoreDistij (i, j, data[i], data[j], seqLength,
					simExp, simDiags, scoreMatrix, alphabet, alphabetSize);
	}

#ifdef _OPENMP
	}
#endif

	free (simDiags);

	return;
}

/*********************************************************/

void gapCheckFilter (int *filter, int itype, int seqlength, int numSeqs, char **data)
{
	int i, j;

	if (PROTEIN == itype)
	{
		for (i=0; i<seqlength; i++)
			for (j=0; j<numSeqs; j++)
			{
				if (NULL == strchr (PROTEIN_ALPHABET, data[j][i]))
				{
					filter[i] = 0;
					if (verbose > 2 && !isBoostrap)
						Debug ( (char*)"Removing site %d.", i);

					break;
				}
				else if (('*' == data[j][i]) || ('?' == data[j][i]) ||
						('-' == data[j][i]))
				{
					filter[i] = 0;
					if (verbose > 2 && !isBoostrap)
						Debug ( (char*)"Removing site %d.",i);

					break;
				}
			}
	}
	else
	{
		for (i=0; i<seqlength; i++)
			for (j=0; j<numSeqs; j++)
			{
				if (NULL == strchr (DNA_ALPHABET, data[j][i]))
				{
					filter[i] = 0;
					if (verbose > 2 && !isBoostrap)
						Debug ( (char*)"Removing site %d.", i);

					break;
				}
				else if (('*' == data[j][i]) || ('?' == data[j][i]) ||
						('-' == data[j][i]))
				{
					filter[i] = 0;
					if (verbose > 2 && !isBoostrap)
						Debug ( (char*)"Removing site %d.", i);

					break;
				}
			}
	}

	return;
}

/*********************************************************/

int gapCheckProportion (char **data, int numSeqs, int numSites, int itype,
	int *filter, FILE *fpO, boolean gapCheck)
{
	int numSelected, numGaps, i;
	int *gapFilterWarn;
	
	gapFilterWarn = initOneArray (numSites);
	gapCheckFilter (gapFilterWarn, itype, numSites, numSeqs, data);
	numSelected = support (gapFilterWarn, numSites);
	numGaps = numSites - numSelected;

	if (!isBoostrap && numSelected < numSites / 2)
	{
		float p = roundf ((float)numGaps / (float)numSites * 100);
		fprintf (fpO, "\t%d%% of sites contain gaps.\n\n", (int) p);
		Warning ( (char*)"%d%% of sites contain gaps.", (int) p);	
	}

	if (gapCheck)
	{
		for (i=0; i<numSites; i++)
			filter[i] = gapFilterWarn[i];
	}
	else
		numSelected = numSites;

	free (gapFilterWarn);

	return numSelected;
}

/*********************************************************/

double **makeDistMatrix (char **data, int numSeqs, int numSites,
	boolean use_gamma, float gamma, int model, int itype, int *filter,
	boolean gapCheck, FILE *fpO, boolean outputMatrix)
{
	int numSelected;
	double **D;

	D = initDoubleMatrix (numSeqs);

	// Here we create 'filter' for global computation of equilibrium frequencies
	// The frequencies are counted on the entire alignment
	// 'filter' may contain '0' values for gap containing sites if 'gapCheck' is TRUE
	numSelected = gapCheckProportion (data, numSeqs, numSites, itype,
		filter, fpO, gapCheck);

	switch (model)
	{
		case F81LIKE:
			computePoisson (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case PDIST:
			computePDIST (data, numSeqs, numSites, numSelected, itype,
				filter, D, gapCheck, outputMatrix);
			break;

		case RY:
			computeRY (data, numSeqs, numSites, numSelected,
				use_gamma, gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case RYSYM:
			computeRYSYM (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case JC69:
			computeJC69 (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case F81:
			computeF81 (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case F84:
			computeF84 (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case K2P:
			computeK2P (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;

		case TN93:
			computeTN93 (data, numSeqs, numSites, numSelected, use_gamma,
				gamma, itype, filter, D, gapCheck, outputMatrix);
			break;
			
		case LOGDET:
			computeLOGDET (data, numSeqs, numSites, numSelected, itype,
				filter, D, gapCheck, outputMatrix);
			break;
			
		default:
			Exit ( (char*)"Please specify model for sequence data.");
	}

	if (warnCheckMaxDist (D, numSeqs) && !isBoostrap)
		Warning ( (char*)"Give up this dataset because at least one distance exceeds %.2f.", DNA_DIST_MAX);

	return (D);
}

/*********************************************************/

boolean warnCheckMaxDist (double **D, int numSeqs)
{
	int i, j;
	boolean ret = FALSE;

	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
			// PROT_DIST_MAX is used here when the distance computation is impossible
			// The output value will be 'NA'
			if (D[i][j] >= DNA_DIST_MAX && D[i][j] < PROT_DIST_MAX)
			{
				ret = TRUE;
				D[i][j] = D[j][i] = DNA_DIST_MAX;
			}
		}
	}

	return ret;
}

/*********************************************************/

void computeRYSYM (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	double b;
	double **PStateChanges;
	int i, j, numS;
	int *f;
	boolean abort = FALSE;
				
#ifdef _OPENMP
	#pragma omp parallel for private (i, j, b, f, numS, PStateChanges) if (!isBoostrap)
#endif
	
	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
					
					PStateChanges = initDoubleMatrix (DNA_ALPHABET_SIZE);

					calcTransitionProbs (PStateChanges, data[i], data[j],
						numSites, f, numS, DNA_ALPHABET_SIZE, DNA_ALPHABET);
				
					b = calcTransversionRate (PStateChanges);

					D[j][i] = D[i][j] = calcRYSYM (b, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					freeMatrix (PStateChanges, DNA_ALPHABET_SIZE);
					free (f);
				}
			} // end if abort
		}
	} //end for - end parallel region
	
	if (abort)
		Exit ( (char*)"Unable to compute all distances");
		
	return;
}

/*********************************************************/

void computeK2P (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	double a, b;
	double **PStateChanges;
	int i, j, numS;
	int *f;
	boolean abort = FALSE;
				
#ifdef _OPENMP
	#pragma omp parallel for private (i, j, a, b, f, numS, PStateChanges) if (!isBoostrap)
#endif
	
	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
				
					PStateChanges = initDoubleMatrix (DNA_ALPHABET_SIZE);
					calcTransitionProbs (PStateChanges, data[i], data[j],
						numSites, f, numS, DNA_ALPHABET_SIZE, DNA_ALPHABET);
				
					a = calcTransitionRate (PStateChanges);
					b = calcTransversionRate (PStateChanges);
					D[j][i] = D[i][j] = calcK2P (a, b, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					freeMatrix (PStateChanges, DNA_ALPHABET_SIZE);
					free (f);
				}
			} // end if abort
		}
	} //end for - end parallel region
	
	if (abort)
		Exit ( (char*)"Unable to compute all distances");
	
	return;
}

/*********************************************************/

void computeJC69 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	double b;
	int i, j, numS;
	int *f;
	boolean abort = FALSE;

#ifdef _OPENMP
	#pragma omp parallel for private (i, j, b, f, numS) if (!isBoostrap)
#endif

	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
					
					b = HammingDistance (data[i], data[j], f, numSites, numS);
				
					D[j][i] = D[i][j] = calcJC69 (b, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					free (f);
				}
			} // end if abort
		}
	} //end for - end parallel region

	if (abort)
		Exit ( (char*)"Unable to compute all distances");
	
	return;
}

/*********************************************************/

void computePoisson (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	double loc, b;
	int i, j, numS;
	int *f;
	double *PStationary;
	boolean abort = FALSE;

	// Count the equilibrium frequencies on the entire proteic alignment
	PStationary = calcStationaryProbsGlobal (data, numSeqs, numSites,
		filter, numSelected, 20, "ACDEFGHIKLMNPQRSTVWY");
	loc = 1.0;
	for (i=0; i<20; i++)
		loc -= PStationary[i] * PStationary[i];

#ifdef _OPENMP
	#pragma omp parallel if (!isBoostrap)
	{
#endif

#ifdef _OPENMP
	#pragma omp for private (i, j, b, f, numS)
#endif
	
	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck) {	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
				
					b = HammingDistance (data[i], data[j], f, numSites, numS);
					
					D[j][i] = D[i][j] = calcF81 (loc, b, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					free (f);
				}
			} // end if abort
		}
	} //end for
#ifdef _OPENMP
	} //end parallel region
#endif

	free (PStationary);
	if (abort)
		Exit ( (char*)"Unable to compute all distances");

	return;
}

/*********************************************************/

void computePDIST (char **data, int numSeqs, int numSites, int numSelected,
	int itype, int *filter, double **D, boolean gapCheck, boolean outputMatrix)
{
	double b;
	int i, j, numS;
	int *f;
	boolean abort = FALSE;

#ifdef _OPENMP
	#pragma omp parallel for private (i, j, b, f, numS) if (!isBoostrap)
#endif
	
	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
				
					b = HammingDistance (data[i], data[j], f, numSites, numS);
				
					D[i][j] = D[j][i] = b;

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					free (f);
				}
			} // end if abort
		}
	} //end for - end parallel region

	if (abort)
		Exit ( (char*)"Unable to compute all distances");
		
	return;
}

/*********************************************************/

void computeF81 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	double loc, b;
	int i, j, numS;
	int *f;
	double *PStationary;
	boolean abort = FALSE;
	
#ifdef _OPENMP
	#pragma omp parallel if (!isBoostrap)
	{
#endif

	// Count equilibrium frequencies on entire the alignment
	// 'filter' was computed regarding to the gaps on the alignment
	// if 'gapCheck' is TRUE, 'filter' may contain '0' values else it is full of '1'
	PStationary = calcStationaryProbsGlobal (data, numSeqs, numSites,
		filter, numSelected, DNA_ALPHABET_SIZE, DNA_ALPHABET);
	
	loc = 1 - (PStationary[ADENINE] * PStationary[ADENINE]) -
				(PStationary[GUANINE] * PStationary[GUANINE]) -
				(PStationary[CYTOSINE] * PStationary[CYTOSINE]) -
				(PStationary[THYMINE] * PStationary[THYMINE]);
	
#ifdef _OPENMP
	#pragma omp for private (i, j, b, f, numS)
#endif

	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort)
			{
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;

					b = HammingDistance (data[i], data[j], f, numSites, numS);
				
					D[j][i] = D[i][j] = calcF81 (loc, b, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					free (f);
				}
			} // end if abort
		}
	} //end for

#ifdef _OPENMP
	} //end parallel region
#endif

	free (PStationary);
	if (abort)
		Exit ( (char*)"Unable to compute all distances");

	return;
}

/*********************************************************/

void computeF84 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	int i, j, numS;
	int *f;
	double a, b;
	double Mloc, Nloc, Ploc;
	double *PStationary;
	double **PStateChanges;
	boolean abort = FALSE;

#ifdef _OPENMP
	#pragma omp parallel if (!isBoostrap)
	{
#endif

	// Count equilibrium frequencies on the alignment
	// 'filter' was computed regarding to the gaps on the alignment
	// if 'gapCheck' is TRUE, 'filter' may contain '0' values else it is full of '1'
	PStationary = calcStationaryProbsGlobal (data, numSeqs, numSites,
			filter, numSelected, DNA_ALPHABET_SIZE, DNA_ALPHABET);
	calcF84AuxProbs (PStationary, &Mloc, &Nloc, &Ploc);
	
#ifdef _OPENMP
	#pragma omp for private (i, j, a, b, f, numS, PStateChanges)
#endif

	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
					
					PStateChanges = initDoubleMatrix (DNA_ALPHABET_SIZE);
					calcTransitionProbs (PStateChanges, data[i], data[j],
						numSites, f, numS, DNA_ALPHABET_SIZE, DNA_ALPHABET);
				
					a = calcTransitionRate (PStateChanges);
					b = calcTransversionRate (PStateChanges);
					D[i][j] = D[j][i] = calcF84 (a, b, use_gamma, gamma,
						Mloc, Nloc, Ploc);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					freeMatrix (PStateChanges, DNA_ALPHABET_SIZE);
					free (f);
				}
			} // end if abort
		}
	} //end for

#ifdef _OPENMP
	} //end parallel region
#endif

	free (PStationary);
	if (abort)
		Exit ( (char*)"Unable to compute all distances");

	return;
}

/*********************************************************/

void computeTN93 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	int i, j, numS;
	int *f;
	double a1, a2, b;
	double PAPG, PCPT, PR, PY;
	double *PStationary;
	double **PStateChanges;
	boolean abort = FALSE;
	
#ifdef _OPENMP
	#pragma omp parallel if (!isBoostrap)
	{
#endif

	// Count equilibrium frequencies on the alignment
	// 'filter' was computed regarding to the gaps on the alignment
	// if 'gapCheck' is TRUE, 'filter' may contain '0' values else it is full of '1'
	PStationary = calcStationaryProbsGlobal (data, numSeqs, numSites,
		filter, numSelected, DNA_ALPHABET_SIZE, DNA_ALPHABET);
	calcTNAuxProbs (PStationary, &PAPG, &PCPT, &PR, &PY);

#ifdef _OPENMP
	#pragma omp for private (i, j, a1, a2, b, f, numS, PStateChanges)
#endif

	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
					
					PStateChanges = initDoubleMatrix (DNA_ALPHABET_SIZE);
					calcTransitionProbs (PStateChanges, data[i], data[j],
						numSites, f, numS, DNA_ALPHABET_SIZE, DNA_ALPHABET);
					
					a1 = PStateChanges[ADENINE][GUANINE] +
						PStateChanges[GUANINE][ADENINE];	//purine transition rate
					a2 = PStateChanges[CYTOSINE][THYMINE] +
						PStateChanges[THYMINE][CYTOSINE];	//pyrimidine transition rate
					b = calcTransversionRate (PStateChanges);
					D[i][j] = D[j][i] = calcTN93 (a1, a2, b, PR, PY, PAPG,
						PCPT, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					freeMatrix (PStateChanges, DNA_ALPHABET_SIZE);
					free (f);
				}
			} // end if abort
		}
	} //end for

#ifdef _OPENMP
	} //end parallel region
#endif

	free (PStationary);
	if (abort)
		Exit ( (char*)"Unable to compute all distances");

	return;
}

/*********************************************************/

void computeRY (char **data, int numSeqs, int numSites,
	int numSelected, boolean use_gamma, float gamma, int itype,
	int *filter, double **D, boolean gapCheck, boolean outputMatrix)
{
	int i, j, numS;
	int *f;
	double b;
	double PR, PY;
	double *PStationary;
	double **PStateChanges;
	boolean abort = FALSE;

#ifdef _OPENMP
	#pragma omp parallel if (!isBoostrap)
	{
#endif

	PStationary = calcStationaryProbsGlobal (data, numSeqs, numSites,
		filter, numSelected, DNA_ALPHABET_SIZE, DNA_ALPHABET);
	PR = PStationary[ADENINE] + PStationary[GUANINE];
	PY = PStationary[CYTOSINE] + PStationary[THYMINE];

#ifdef _OPENMP
	#pragma omp for private (i, j, b, f, numS, PStateChanges)
#endif
	
	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
						
					PStateChanges = initDoubleMatrix (DNA_ALPHABET_SIZE);
					calcTransitionProbs (PStateChanges, data[i], data[j],
						numSites, f, numS, DNA_ALPHABET_SIZE, DNA_ALPHABET);
				
					b = calcTransversionRate (PStateChanges);
					D[i][j] = D[j][i] = XX (PR, PY, b, use_gamma, gamma);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					freeMatrix (PStateChanges, DNA_ALPHABET_SIZE);
					free (f);
				}
			} // end if abort
		}
	} //end for

#ifdef _OPENMP
	} //end parallel region
#endif

	free (PStationary);
	if (abort)
		Exit ( (char*)"Unable to compute all distances");

	return;
}

/*********************************************************/

void computeLOGDET (char **data, int numSeqs, int numSites,
	int numSelected, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix)
{
	int i, j, numS;
	int *f;
	double **PStateChanges, **Pi2;
	boolean abort = FALSE;

#ifdef _OPENMP
	#pragma omp parallel if (!isBoostrap)
	{
#endif

	Pi2 = (double **) mCalloc (numSeqs, sizeof (double *));

#ifdef _OPENMP
	#pragma omp for private (i)
#endif

	for (i=0; i<numSeqs; i++)
		Pi2[i] = calcStationaryProbsGlobal (data + i, 1, numSites,
			filter, numSelected, DNA_ALPHABET_SIZE, DNA_ALPHABET);

#ifdef _OPENMP
	#pragma omp for private (i, j, f, numS, PStateChanges)
#endif
	
	for (i=0; i<numSeqs-1; i++)
	{
		for (j=i; j<numSeqs; j++)
		{
#ifdef _OPENMP
			#pragma omp flush (abort)
#endif
			if (!abort) {
				if (i == j)
					D[i][j] = 0.0;
				else
				{
					f = copyFilter (filter, numSites);
					if (!gapCheck)	// Pairwise deletion of gaps
						// => compute a new 'filter' for each pairwise distance computation
					{
						ijFilter (f, data[i], data[j], itype, numSites);
						numS = support (f, numSites);
					}
					else
						numS = numSelected;
				
					PStateChanges = initDoubleMatrix (DNA_ALPHABET_SIZE);
					calcTransitionProbs (PStateChanges, data[i], data[j],
						numSites, f, numS, DNA_ALPHABET_SIZE, DNA_ALPHABET);

					D[i][j] = D[j][i] = logdet (PStateChanges, Pi2[i], Pi2[j]);

					if (numS == 0) {
						if (outputMatrix)
							D[j][i] = D[i][j] = PROT_DIST_MAX +1;
						else
						{
							abort = TRUE;
#ifdef _OPENMP
							#pragma omp flush (abort)
#endif
						}
					}
					
					freeMatrix (PStateChanges, DNA_ALPHABET_SIZE);
					free (f);
				}
			} // end if abort
		}
	} //end for

#ifdef _OPENMP
	} //end parallel region
#endif

	freeMatrix (Pi2, numSeqs);
	if (abort)
		Exit ( (char*)"Unable to compute all distances");
	
	return;
}

/*********************************************************/

void symmetrizeDoubleMatrix (double **X, int n)
{
	int i, j;

	for (i=0; i<n-1; i++)
		for (j=i+1; j<n; j++)
			X[i][j] = X[j][i] = 0.5 * (X[i][j] + X[j][i]);

	return;
}

