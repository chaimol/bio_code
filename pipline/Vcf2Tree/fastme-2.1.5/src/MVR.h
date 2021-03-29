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


#ifndef MVR_H_
#define MVR_H_

#include "BIONJ.h"

tree *unj (double **D, set *species, int n, const char *format);
int SymmetrizeMVR (double **delta, int n);
double Finish_branch_length_MVR (int i, int j, int k, double **delta, int n);
void FinishStrMVR (double **delta, int n, POINTERS *trees, char *StrTree, const char *format);
void Branch_lengthMVR (int a, int b, double *la, double *lb, double **delta, int n);
double Reduction10MVR (int a, double la, int b, double lb, int i, double lamda, double **delta);
double Reduction11MVR (int x, int y, int i, double **delta);
double LamdaMVR (int x, int y, int i, double **delta);
double WeightMVR (int x, int y, int i, double **delta, double u);
double mu (int a, int b, double **delta, int n);

#endif /*MVR_H_*/

