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


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         BIONJ program                                     ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;
;                         GERAD - Montreal- Canada                          ;
;                         olivierg@crt.umontreal.ca                         ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                           ;
;                         UNIX version, written in C                        ;
;                         by Hoa Sien Cuong (Univ. Montreal)                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


#ifndef BIONJ_H_
#define BIONJ_H_

#include "newick.h"

typedef struct word
{
	char name[MAX_NAME_LENGTH];
	struct word *suiv;
} WORD;

typedef struct pointers
{
	WORD *head;
	WORD *tail;
} POINTERS;


tree *bionj (double **D, set *species, int n, boolean isNJ, const char *format);
void Initialize (double **D, set *species, double **delta, POINTERS *trees, int n);
void Print_outputChar (int i, POINTERS *trees, char *output);
boolean Symmetrize (double **delta, int n);
void Concatenate (char chain1[MAX_NAME_LENGTH], int ind, POINTERS *trees, int post);
double Distance (int i, int j, double **delta);
double Variance (int i, int j, double **delta);
int Emptied (int i, double **delta);
double Sum_S (int i, double **delta);
void Compute_sums_Sx (double **delta, int n);
void Best_pair (double **delta, int r, int *a, int *b, int n);
double Finish_branch_length (int i, int j, int k, double **delta);
void FinishStr (double **delta, int n, POINTERS *trees, char *StrTree, const char *format);
double Agglomerative_criterion (int i, int j, double **delta, int r);
double Branch_length (int a, int b, double **delta, int r);
double Reduction4 (int a, double la, int b, double lb, int i, double lamda, double **delta);
double Reduction10 (int a, int b, int i, double lamda, double vab, double **delta);
double Lamda (int a, int b, double vab, double **delta, int n, int r);

#endif /*BIONJ_H_*/

