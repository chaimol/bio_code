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


#ifndef NNI_H_
#define NNI_H_

#include "gme.h"
#include "heap.h"

double **buildAveragesTable (tree *T, double **D);
double wf2 (double lambda, double D_AD, double D_BC, double D_AC,
	double D_BD, double D_AB, double D_CD);
int NNIEdgeTest (edge *e, tree *T, double **A, double *weight);
void NNIupdateAverages (double **A, edge *e, edge *par, edge *skew,
	edge *swap, edge *fixed, tree *T);
void NNItopSwitch (tree *T, edge *e, int direction, double **A);
void NNIRetestEdge (int *p, int *q, edge *e,tree *T, double **avgDistArray,
	double *weights, int *location, int *possibleSwaps);
void NNI (tree *T, double **avgDistArray, int *count, FILE *statfile);


#endif /*NNI_H_*/

