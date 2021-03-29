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


#ifndef INPUTS_H_
#define INPUTS_H_

#include "graph.h"
#include "heap.h"

void compareSets (tree *T, set *S);
void freeCharMatrix (char **D, int size);
void freeMatrix (double **D, int size);
boolean isTrgMatrix (FILE *ifile);
double **loadM (FILE *ifile, int *size, set *S);
double **loadMatrix (FILE *ifile, int *size, set *S);
void loadTrgMatrix (FILE *ifile, int size, set *S, double **table);
void loadRectMatrix (FILE *ifile, int size, set *S, double **table);
void partitionSizes (tree *T);
double **loadScoreMatrix (int *d, FILE *ifile, char *alphabet);

#endif /*INPUTS_H_*/

