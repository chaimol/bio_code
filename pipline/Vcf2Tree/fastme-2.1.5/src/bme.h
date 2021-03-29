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


#ifndef BME_H_
#define BME_H_

#include "traverse.h"

void BalWFext (edge *e, double **A);
void BalWFint (edge *e, double **A);
void assignBMEWeights (tree *T, double **A);
void BMEcalcDownAverage (node *v, edge *e, double **D, double **A);
void BMEcalcUpAverage (tree *T, node *v, edge *e, double **D, double **A);
void BMEcalcNewvAverages (tree *T, node *v, double **D, double **A);
void updatePair (double **A, edge *nearEdge, edge *farEdge, node *v,
	node *root, double dcoeff, int direction);
void updateSubTree (double **A, edge *nearEdge, node *v, node *root,
	node *newNode, double dcoeff, int direction);
void BMEupdateAveragesMatrix (double **A, edge *e, node *v,node *newNode);
double wf3 (double D_AB, double D_AC, double D_kB, double D_kC);
void BMEtestEdge (edge *e, node *v, double **A);
void BMEsplitEdge (tree *T, node *v, edge *e, double **A);
tree *BMEaddSpecies (tree *T,node *v, double **D, double **A);
void calcUpAverages (double **D, double **A, edge *e, edge *g);
void makeBMEAveragesTable (tree *T, double **D, double **A);

#endif /*BME_H_*/

