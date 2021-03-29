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


#ifndef DISTANCE_H_
#define DISTANCE_H_


#include "random.h"
#include "inputs.h"


int countStateChanges (char *s, char *t, int length, char c1, char c2,
	int *filter);

int *copyFilter (int *filter, int l);

void ijFilter (int *filter, char *s1, char *s2, int itype, int seqlength);

int seqCharMatches (char *s, int length, char c, int *filter);

int matrixCharMatches (char **s, int numSeqs, int length, char c,
	int *filter);

double *calcStationaryProbsGlobal (char **s, int numSeqs, int length,
	int *filter, int numSelected, int alphabetSize, const char *alphabet);

void calcTransitionProbs (double **P, char *s1, char *s2, int length,
	int *filter, int numSelected, int alphabetSize, const char *alphabet);

double calcTransitionRate (double **P);

double calcTransversionRate (double **P);

double calcRYSYM (double b, boolean use_gamma, float gamma);

double calcJC69 (double b, boolean use_gamma, float gamma);

double calcF81 (double loc, double b, boolean use_gamma, float gamma);

double calcK2P (double a, double b, boolean use_gamma, float gamma);

void calcF84AuxProbs (double *Pi, double *Mloc, double *Nloc,
	double *Ploc);

double calcF84 (double a, double b, boolean use_gamma, float gamma,
	double Mloc, double Nloc, double Ploc);

void calcTNAuxProbs (double *Pi, double *PAPG, double *PCPT, double *PR,
	double *PY);

double calcTN93 (double aR, double aY, double b, double PR, double PY,
	double PAPG, double PCPT, boolean use_gamma, float gamma);

double XX (double PR, double PY, double b, boolean use_gamma,
	float gamma);

int factorial (int n);

int *nextPerm (int *p, int index, int size, int length);

double permDiagProduct (double **P, int *p, int d);

double det (double **P, int d);

double logdet (double **P, double *Pi1, double *Pi2);

int support (int *v, int length);

double HammingDistance (char *v1, char *v2, int *filter, int length,
	int numSelected);

double protDiff (double *P);

double protFormula (double b, boolean use_gamma, float gamma,
	double pdiff);

int aaIndex (char s, char *alphabet, int d);

double simScore (char *s1, char *s2, double **scoreMatrix, int seqlength,
	char *alphabet, int alphabetSize);

double expectedProtSimScore (double *P, double **scoreMatrix,
	int alphabetSize);

double scoreDistij (int i, int j, char *si, char *sj, int seqLength,
	double simExp, double *simDiags, double **scoreMatrix, char *alphabet,
	int alphabetSize);

void scoreDist (double *P, char **data, int numSpecies, int seqLength,
	double **scoreMatrix, double **D, char *alphabet, int alphabetSize);

void gapCheckFilter (int *filter, int itype, int seqlength, int numSeqs,
	char **data);
	
int gapCheckProportion (char **data, int numSeqs, int numSites, int itype,
	int *filter, FILE *fpO, boolean gapCheck);

double **makeDistMatrix (char **data, int numSeqs, int numSites,
	boolean use_gamma, float gamma, int model, int itype, int *filter,
	boolean gapCheck, FILE *fpO, boolean outputMatrix);

boolean warnCheckMaxDist (double **D, int numSeqs);

void computeRY (char **data, int numSeqs, int numSites,
	int numSelected, boolean use_gamma, float gamma, int itype,
	int *filter, double **D, boolean gapCheck, boolean outputMatrix);

void computeRYSYM (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computeK2P (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computeJC69 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computePoisson (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computePDIST (char **data, int numSeqs, int numSites, int numSelected,
	int itype, int *filter, double **D, boolean gapCheck, boolean outputMatrix);

void computeF81 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computeF84 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computeTN93 (char **data, int numSeqs, int numSites, int numSelected,
	boolean use_gamma, float gamma, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void computeLOGDET (char **data, int numSeqs, int numSites,
	int numSelected, int itype, int *filter, double **D,
	boolean gapCheck, boolean outputMatrix);

void symmetrizeDoubleMatrix (double **X, int n);


#endif /*DISTANCE_H_*/

