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


#ifndef BNNI_H_
#define BNNI_H_

#include "bme.h"
#include "heap.h"

void bNNIRetestEdge (int *p, int *q, edge *e,tree *T, double **avgDistArray,
	double *weights, int *location, int *possibleSwaps);
void bNNItopSwitch (edge *e, int direction, double **A);
void bNNI (tree *T, double **avgDistArray, int *count, FILE *statfile);
void updateSubTreeAfterNNI (double **A, node *v, edge *rootEdge, node *closer,
	node *further, double dcoeff, int direction);
void bNNIupdateAverages (double **A, node *v, edge *par, edge *skew,
	edge *swap, edge *fixed);
double wf5 (double D_AD, double D_BC, double D_AC, double D_BD,
	double D_AB, double D_CD);
int bNNIEdgeTest (edge *e, tree *T, double **A, double *weight);
void limitedFillTableUp (edge *e, edge *f, double **A, edge *trigger);

#endif /*BNNI_H_*/

