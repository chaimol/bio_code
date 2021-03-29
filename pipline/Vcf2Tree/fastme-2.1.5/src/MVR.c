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
;                         MVR implementation of UNJ                         ;
;                                                                           ;
;                         Vincent Lefort                                    ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         vincent.lefort@lirmm.fr                           ;
;                                                                           ;
;                         code is heavily borrowed from BIONJ :             ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

#include "MVR.h"

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Main program                                      ;
;                                                                           ;
; input       :                                                             ;
;              double **D     : the distance matrix                         ;
;              set *species   : the correponding taxa                       ;
;              int n          : the number of species                       ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

tree *unj (double **D, set *species, int n, const char *format)
{
	POINTERS *trees;			/* list of subtrees            */
	tree *ret;					/* the returned tree           */
	char *chain1;				/* stringized branch-length    */
	char *str;					/* the string containing final tree */
	int *a, *b;					/* pair to be agglomerated     */
	double **delta;				/* delta matrix                */
	double *la;					/* first taxon branch-length   */
	double *lb;					/* second taxon branch-length  */
	double lamda;
	int i;
	int ok;
	int r;						/* number of subtrees          */
	int x, y;

	a = (int*) mCalloc (1, sizeof(int));
	b = (int*) mCalloc (1, sizeof(int));
	la = (double*) mCalloc (1, sizeof(double));
	lb = (double*) mCalloc (1, sizeof(double));
	chain1 = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
	str = (char *) mCalloc (MAX_INPUT_SIZE, sizeof(char));

	/*      Create the delta matrix     */
	delta = (double **) mCalloc (n + 1, sizeof (double*));
	for (i=1; i<= n; i++)
		delta[i] = (double *) mCalloc (n + 1, sizeof (double));

	trees = (POINTERS *) mCalloc (n + 1, sizeof (POINTERS));

	/*   initialise and symmetrize the running delta matrix    */
	r=n;
	*a=0;
	*b=0;
	Initialize (D, species, delta, trees, n);

	ok = SymmetrizeMVR (delta, n);
	if (!ok)
		Warning ( (char*)"The matrix is not symmetric.");

	while (r > 3)
	{
		Compute_sums_Sx(delta, n);			/* compute the sum Sx       */
		Best_pair(delta, r, a, b, n);		/* find the best pair by    */
											/* minimizing (1)           */

		Branch_lengthMVR (*a, *b, la, lb, delta, n);
		for (i=1; i <= n; i++)
		{
			if (!Emptied (i,delta) && (i != *a) && (i != *b))
			{
				lamda = LamdaMVR (*a, *b, i, delta);	/* compute lambda using (10)*/
				if (*a > i)
				{
					x=*a;
					y=i;
				}
				else
				{
					x=i;
					y=*a;				/* apply reduction formulae */
				}						/* 10 and 11 to delta       */

				delta[x][y] = Reduction10MVR (*a, *la, *b, *lb, i, lamda, delta);
				delta[y][x] = Reduction11MVR (*a, *b, i, delta);
			}
		}

		strncpy (chain1, "", 1);				/* agglomerate the subtrees */
		strncat (chain1, "(",1);				/* a and b together with the*/
		Concatenate (chain1, *a, trees, 0);		/* branch-lengths according */
		strncpy (chain1, "", 1);				/* to the NEWWICK format    */
		strncat (chain1, ":", 1);
		//snprintf (chain1+strlen(chain1), MAX_NAME_LENGTH, "%f", *la);
		snprintf (chain1+strlen(chain1), MAX_NAME_LENGTH, format, *la);

		strncat (chain1, ",", 1);
		Concatenate (chain1, *a, trees, 1);
		trees[*a].tail->suiv = trees[*b].head;
		trees[*a].tail = trees[*b].tail;
		strncpy (chain1, "", 1);
		strncat( chain1, ":", 1);
		//snprintf (chain1+strlen(chain1), MAX_NAME_LENGTH, "%f", *lb);
		snprintf (chain1+strlen(chain1), MAX_NAME_LENGTH, format, *lb);

		strncat (chain1, ")", 1);
		Concatenate (chain1, *a, trees, 1);
		delta[*b][0] = 1.0;						/* make the b line empty     */
		trees[*b].head = NULL;
		trees[*b].tail = NULL;
		r = r-1;
	}

	FinishStrMVR (delta, n, trees, str, format);	/* compute the branch-lengths*/
													/* of the last three subtrees*/
	ret = readNewickString (str);
	ret = detrifurcate (ret);

	for (i=1; i<= n; i++)
		free (delta[i]);

	free (delta);
	free (trees);
	free (str);
	free (chain1);
	free (a);
	free (b);
	free (la);
	free (lb);

	return (ret);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                             Utilities                                     ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/



/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Symmetrize  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function checks the delta matrix is symmetric;         ;
;               if not the matrix is made symmetric.                        ;
;                                                                           ;
; input       :                                                             ;
;              double **delta : delta matrix                                ;
;              int n          : number of taxa                              ;
;                                                                           ;
; return value:                                                             ;
;              int symmetric  : indicate if the matrix has been made        ;
;                               symmetric or not                            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

int SymmetrizeMVR (double **delta, int n)
{
	int lig;			/* matrix line        */
	int col;			/* matrix column      */
	double value;		/* symmetrized value  */
	int symmetric;

	symmetric = 1;

	for (lig=1; lig<=n; lig++)
	{
		for (col=1; col<lig; col++)
		{
			if ( (delta[lig][col] - delta[col][lig] > DBL_EPSILON) || (delta[col][lig] - delta[lig][col] > DBL_EPSILON) )
			{
				value = (delta[lig][col] + delta[col][lig]) / 2;
				delta[lig][col] = value;
				symmetric = 0;
			}
			delta[col][lig] = 1;
		}
	}

	return (symmetric);
}


/*;;;;;;;;;;;;;;;;;;;;;;Finish_branch_length;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description :  Compute the length of the branch attached                 ;
;                 to the subtree i, during the final step                   ;
;                                                                           ;
;  input       :                                                            ;
;                int i          : position of subtree i                     ;
;                int j          : position of subtree j                     ;
;                int k          : position of subtree k                     ;
;                double **delta :                                           ;
;                                                                           ;
;  return value:                                                            ;
;                double length  : The length of the branch                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

double Finish_branch_length_MVR (int i, int j, int k, double **delta, int n)
{
	double length, Wi, Dxu, u;

	u = mu (i, j, delta, n);
	Wi = WeightMVR (i, j, k, delta, u);
	Dxu = Wi * (Distance (i, j, delta) + Distance (i, k, delta) -
		Distance (j, k, delta));
	length = Dxu / (2 * Wi);

	return (length);
}

/*********************************************************/

void FinishStrMVR (double **delta, int n, POINTERS *trees, char *StrTree, const char *format)
{
	int l = 1;
	int i = 0;
	double length;
	char *tmp;
	WORD *bidon;
	WORD *ele;
	int last[3];

	while (l <= n)
	{								// find the last subtree
		if (!Emptied (l, delta))
		{
			last[i] = l;
			i++;
		}
		l++;
	}

	tmp = (char*) mCalloc (MAX_NAME_LENGTH, sizeof(char));
	StrTree[0] = '(';

	length = Finish_branch_length_MVR (last[0], last[1], last[2], delta, n);
	Print_outputChar (last[0], trees, StrTree);
	//snprintf (tmp, 16, "%.6f,", length);
	snprintf (tmp, DECIMAL_DIG+2, format, length);
	if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-2)
	{
		strncat (StrTree, ":", 1);
		strncat (StrTree, tmp, strlen (tmp));
		strncat (StrTree, ",", 1);
	}

	length = Finish_branch_length_MVR (last[1], last[0], last[2], delta, n);
	Print_outputChar (last[1], trees, StrTree);
	snprintf (tmp, INPUT_SIZE, format, length);
	if (n>2)
		strncat (tmp, ",", 1);
/*		snprintf (tmp, 16, "%.6f,", length);
	else
		snprintf (tmp, 16, "%.6f", length);
*/
	if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-1)
	{
		strncat (StrTree, ":", 1);
		strncat (StrTree, tmp, strlen (tmp));
	}

	if (n>2) 
	{
		length = Finish_branch_length_MVR (last[2], last[1], last[0], delta, n);
		Print_outputChar (last[2], trees, StrTree);
		//snprintf (tmp, 16, "%.6f", length);
		snprintf (tmp, INPUT_SIZE, format, length);
		if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-1)
		{
			strncat (StrTree, ":", 1);
			strncat (StrTree, tmp, strlen (tmp));
			//strncat (StrTree, ",", 1);
		}
	}

	if (strlen (StrTree) < MAX_INPUT_SIZE-3)
		strncat (StrTree, ");\n", 3);

	if (n>3)
		n=3;
		
	for(i=0; i<n; i++)
	{
		bidon = trees[last[i]].head;
		ele = bidon;
		while (bidon != NULL)
		{
			ele = ele->suiv;
			free (bidon);
			bidon = ele;
		}
	}

	free (tmp);

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                          Formulae                                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void Branch_lengthMVR (int a, int b, double *la, double *lb, double **delta, int n)
{
	int i;
	double Dxu, Dyu, Wi, u;

	u = mu (a, b, delta, n);
	Dxu = Dyu = 0;

	for(i=1; i<=n ; i++)
	{
		if (!Emptied (i, delta) && (i != a) && (i != b))
		{
			Wi = WeightMVR (a, b, i, delta, u);
			Dxu = Dxu + (Wi * (Distance (a, b, delta) + Distance (a, i, delta) -
				Distance (b, i, delta)));
			Dyu = Dyu + (Wi * (Distance (a, b, delta) + Distance (b, i, delta) -
				Distance (a, i, delta)));
		}
	}

	*la = Dxu;
	*lb = Dyu;

	return;
}

/*********************************************************/

double Reduction10MVR (int a, double la, int b, double lb, int i, double lamda, double **delta)
{
	double Dui;

	Dui = (lamda * Distance (a, i, delta)) + ((1-lamda) * (Distance (b, i, delta))) -
		(lamda * la) - ((1-lamda) * lb);

	return (Dui);
}

/*********************************************************/

double Reduction11MVR (int x, int y, int i, double **delta)
{
	double Vui, varsum;

	varsum = Variance (x, i, delta) + Variance (y, i, delta);

	if (varsum < DBL_EPSILON)
		varsum = DBL_EPSILON;

	Vui = (Variance (x, i, delta) * Variance (y, i, delta)) / varsum;

	return (Vui);
}

/*********************************************************/

double LamdaMVR (int x, int y, int i, double **delta)
{
	double lamda, varsum;										/* Formula (10) */

	varsum = Variance (x, i, delta) + Variance (y, i, delta);

	if (varsum < DBL_EPSILON)
		varsum = DBL_EPSILON;

	lamda = Variance (y, i, delta) / varsum;

	return (lamda);
}

/*********************************************************/

double WeightMVR (int x, int y, int i, double **delta, double u)
{
	double Wi, varsum;											/* Formula (14) */

	varsum = Variance (x, i, delta) + Variance (y, i, delta);

	if (varsum < DBL_EPSILON)
		varsum = DBL_EPSILON;

	Wi = u / varsum;

	return (Wi);
}

/*********************************************************/

double mu (int a, int b, double **delta, int n)
{
	double u, Sv, varsum;										/* Formula (14) */
	int i;

	Sv = 0;
	for(i=1; i<=n ; i++)
	{
		if (!Emptied (i, delta) && (i != a) && (i != b))
		{
			varsum = Variance (a, i, delta) + Variance (b, i, delta);

			if (varsum < DBL_EPSILON)
				varsum = DBL_EPSILON;

			Sv += 1 / varsum;
		}
	}

	if (Sv < DBL_EPSILON)
		Sv = DBL_EPSILON;

	u = 0.5 * (1 / Sv);

	return u;
}

