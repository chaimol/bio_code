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


#ifndef SPR_H_
#define SPR_H_

#include "bme.h"
#include "traverse.h"
#include "newick.h"
#include "inputs.h"

#ifndef SPR_EPSILON
#define SPR_EPSILON .0000001
#endif

void zero3DMatrix (double ***X, int h, int l, int w);
void findTableMin (int *imin, int *jmin, int *kmin, int n, double ***X, double *min);
void SPR (tree *T, double **D, double **A, int *count, FILE *statfile);
void assignSPRWeights (node *vtest, double **A, double ***swapWeights);
void assignDownWeightsUp (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignDownWeightsSkew (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignDownWeightsDown (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignUpWeights (edge *etest, node *vtest, node *va, edge *back, node *cprev,
	double oldD_AB, double coeff, double **A, double ***swapWeights);
void pruneSubtree (edge *p, edge *u, edge *d);
void SPRsplitEdge (edge *e, edge *p, edge *d);
void SPRDownShift (node *v, edge *e);
void SPRUpShift (node *vmove, edge *esplit);
void SPRTopShift (node *vmove, edge *esplit, int UpOrDown);

#endif /*SPR_H_*/
