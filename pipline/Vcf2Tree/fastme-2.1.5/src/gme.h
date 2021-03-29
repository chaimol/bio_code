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


#ifndef GME_H_
#define GME_H_

#include "traverse.h"

void fillTableUp (edge *e, edge *f, double **A, double **D, tree *T);
void OLSext (edge *e, double **A);
double wf (double lambda, double D_LR, double D_LU, double D_LD,
	double D_RU, double D_RD, double D_DU);
void OLSint (edge *e, double **A);
void assignOLSWeights (tree *T, double **A);
void makeOLSAveragesTable (tree *T, double **D, double **A);
void GMEcalcDownAverage (node *v, edge *e, double **D, double **A);
void GMEcalcUpAverage (node *v, edge *e, double **D, double **A);
void GMEcalcNewvAverages (tree *T, node *v, double **D, double **A);
double wf4 (double lambda, double lambda2, double D_AB, double D_AC,
	double D_BC, double D_Av, double D_Bv, double D_Cv);
void testEdge (edge *e, node *v, double **A);
tree *GMEaddSpecies (tree *T,node *v, double **D, double **A);
void GMEupdateAveragesMatrix (double **A, edge *e, node *v, node *newNode);
void GMEsplitEdge (tree *T, node *v, edge *e, double **A);
void updateSubTreeAverages (double **A, edge *e, node *v, int direction);
void assignBottomsize (edge *e);
void assignTopsize (edge *e, int numLeaves);
void assignAllSizeFields (tree *T);

#endif /*GME_H_*/

