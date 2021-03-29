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

#include "BIONJ.h"

/*;;;;;;;;;;;  INPUT, OUTPUT, INITIALIZATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                                                                           ;
;              The delta matrix is read from the input-file.                ;
;              It is recommended to put it and the executable in            ;
;              a special directory. The input-file and output-file          ;
;              can be given as arguments to the executable by               ;
;              typing them after the executable (Bionj input-file           ;
;              output-file) or by typing them when asked by the             ;
;              program. The input-file has to be formated according         ;
;              the PHYLIP standard. The output file is formated             ;
;              according to the NEWICK standard.                            ;
;                                                                           ;
;              The lower-half of the delta matrix is occupied by            ;
;              dissimilarities. The upper-half of the matrix is             ;
;              occupied by variances. The first column                      ;
;              is initialized as 0; during the algorithm some               ;
;              indices are no more used, and the corresponding              ;
;              positions in the first column are set to 1.                  ;
;                                                                           ;
;              This delta matix is made symmetrical using the rule:         ;
;              Dij = Dji <- (Dij + Dji)/2. The diagonal is set to 0;        ;
;              during the further steps of the algorithm, it is used        ;
;              to store the sums Sx.                                        ;
;                                                                           ;
;              A second array, trees, is used to store taxon names.         ;
;              During the further steps of the algoritm, some               ;
;              positions in this array are emptied while the others         ;
;              are used to store subtrees.                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Main program                                      ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

tree *bionj (double **D, set *species, int n, boolean isNJ, const char *format)
{
	POINTERS *trees;		/* list of subtrees            */
	tree *ret;				/* the returned tree           */
	char *chain1;			/* stringized branch-length    */
	int *a, *b;				/* pair to be agglomerated     */
	double **delta;			/* delta matrix                */
	double la;				/* first taxon branch-length   */
	double lb;				/* second taxon branch-length  */
	double vab;				/* variance of Dab             */
	int r;					/* number of subtrees          */
	char *str;				/* the string containing the final tree */
	double lamda = 0.5;
	int i, x, y;
	boolean ok;

	a = (int*) mCalloc (1, sizeof (int));
	b = (int*) mCalloc (1, sizeof (int));
	chain1 = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
	str = (char *) mCalloc (MAX_INPUT_SIZE, sizeof (char));

	/* Create the delta matrix */
	delta = (double **) mCalloc (n+1, sizeof (double*));
	for (i=1; i<= n; i++)
		delta[i] = (double *) mCalloc(n+1, sizeof (double));

	trees = (POINTERS *) mCalloc (n+1, sizeof (POINTERS));

	/* Initialise and symmetrize the running delta matrix */
	r = n;
	*a = 0;
	*b = 0;
	Initialize (D, species, delta, trees, n);
	ok = Symmetrize (delta, n);
	if (!ok)
		Warning ( (char*)"The matrix is not symmetric.");

	while (r > 3)										/* until r=3                */
	{
		Compute_sums_Sx (delta, n);						/* compute the sum Sx       */
		Best_pair (delta, r, a, b, n);					/* find the best pair by    */
		vab = Variance (*a, *b, delta);					/* minimizing (1)           */
		la = Branch_length (*a, *b, delta, r);			/* compute branch-lengths   */
		lb = Branch_length (*b, *a, delta, r);			/* using formula (2)        */
		if (!isNJ)
			lamda = Lamda (*a, *b, vab, delta, n, r);	/* compute lambda* using (9)*/

		for (i=1; i <=n; i++)
		{
			if (!Emptied (i, delta) && (i != *a) && (i != *b))
			{
				if (*a > i)
				{
					x=*a;
					y=i;
				}
				else
				{
					x=i;
					y=*a;					/* apply reduction formulae */
				}							/* 4 and 10 to delta        */
				delta[x][y] = Reduction4 (*a, la, *b, lb, i, lamda, delta);
				delta[y][x] = Reduction10 (*a, *b, i, lamda, vab, delta);
			}
		}

		strncpy (chain1, "", 1);				/* agglomerate the subtrees  */
		strncat (chain1, "(", 1);				/* a and b together with the */
		Concatenate (chain1, *a, trees, 0);		/* branch-lengths according  */
		strncpy (chain1, "", 1);				/* to the NEWWICK format     */
		strncat (chain1, ":", 1);
		//snprintf (chain1 + strlen (chain1), 16, "%.6f", la);
		snprintf (chain1 + strlen (chain1), MAX_NAME_LENGTH, format, la);

		strncat (chain1, ",", 1);
		Concatenate (chain1, *a, trees, 1);
		trees[*a].tail->suiv = trees[*b].head;
		trees[*a].tail = trees[*b].tail;
		strncpy (chain1, "", 1);
		strncat (chain1, ":", 1);
		//snprintf (chain1 + strlen (chain1), 16, "%.6f", lb);
		snprintf (chain1 + strlen (chain1), MAX_NAME_LENGTH, format, lb);

		strncat (chain1, ")", 1);
		Concatenate (chain1, *a, trees, 1);
		delta[*b][0]=1.0;					/* make the b line empty */
		trees[*b].head = NULL;
		trees[*b].tail = NULL;
		r = r-1;
	}

	FinishStr (delta, n, trees, str, format);	/* compute the branch-lengths */
												/* of the last three subtrees */	

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

	return (ret);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;; Initialize ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function reads an input distance matrix and list of    ;
;               taxa and return the corresponding BIONJ delta matrix        ;
                and trees (list of taxa) in the right data structure.       ;
;                                                                           ;
; input       :                                                             ;
;              double **D      : input distance matrix                      ;
;              set *species    : the corresponding taxa                     ;
;              double **delta  : output delta matrix                        ;
;              POINTERS *trees : list of taxa                               ;
;              int n           : number of species                          ;
;                                                                           ;
; return value:                                                             ;
;              double **delta : delta matrix                                ;
;              char *trees    : list of taxa                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void Initialize (double **D, set *species, double **delta, POINTERS *trees, int n)
{
	int lig;			/* matrix line       */
	int col;			/* matrix column     */
	WORD *name;
	set *reader;

	// read taxon name from species
	lig = 0;
	reader = species;

	do
	{
		lig++;
		name = (WORD *) mCalloc (1, sizeof(WORD));
		strncpy (name->name, reader->firstNode->label, MAX_NAME_LENGTH);
		name->suiv = NULL;
		trees[lig].head = name;
		trees[lig].tail = name;
		reader = reader->secondNode;
	} while (NULL != reader);
	reader = NULL;

	// copy Distance matrix into delta
	for (lig=0; lig<n; lig++)
		for (col=0; col<n; col++)
			delta[lig+1][col+1] = D[lig][col];

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Print_outputChar ;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function prints out the tree in the output string.     ;
;                                                                           ;
; input       :                                                             ;
;              int i           : indicate the subtree i to be printed.      ;
;              POINTERS *trees : pointer to the subtrees.                   ;
:              char *output    : pointer to the output string.              ;
;                                                                           ;
; return value: The phylogenetic tree in the output string.                 ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void Print_outputChar (int i, POINTERS *trees, char *output)
{
	WORD *parcour;

	parcour = trees[i].head;
	while (parcour != NULL && (strlen (output) + strlen (parcour->name) < MAX_INPUT_SIZE))
	{
		output = strncat (output, parcour->name, strlen (parcour->name));
		parcour = parcour->suiv;
	}

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                             Utilities                                     ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/



/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Symmetrize ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifies if the delta matrix is symmetric;    ;
;               if not the matrix is made symmetric.                        ;
;                                                                           ;
; input       :                                                             ;
;              double **delta    : delta matrix                             ;
;              int n             : number of taxa                           ;
;                                                                           ;
; return value:                                                             ;
;              boolean symmetric : indicate if the matrix has been made     ;
;                                  symmetric or not                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

boolean Symmetrize (double **delta, int n)
{
	int lig;			/* matrix line        */
	int col;			/* matrix column      */
	double value;		/* symmetrized value  */
	boolean symmetric = TRUE;

	for (lig=1; lig <= n; lig++)
	{
		for (col=1; col < lig; col++)
		{
			if ( (delta[lig][col] - delta[col][lig] > DBL_EPSILON) || (delta[col][lig] - delta[lig][col] > DBL_EPSILON) )
			{
				value = (delta[lig][col] + delta[col][lig]) / 2;
				delta[lig][col] = value;
				delta[col][lig] = value;
				symmetric = FALSE;
			}
		}
	}

	return (symmetric);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Concatenate ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function concatenates a string to another.             ;
;                                                                           ;
; input       :                                                             ;
;      char *chain1    : the string to be concatenated.                     ;
;      int ind         : indicate the subtree to which concatenate the      ;
;                        string                                             ;
;      POINTERS *trees : pointer to subtrees.                               ;
;      int post        : position to which concatenate (front (0) or        ;
;                        end (1))                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void Concatenate (char chain1[MAX_NAME_LENGTH], int ind, POINTERS *trees, int post)
{
	WORD *bran;

	bran = (WORD *) mCalloc (1, sizeof (WORD));
	strncpy (bran->name, chain1, MAX_NAME_LENGTH);
	bran->suiv = NULL;

	if (post == 0)
	{
		bran->suiv = trees[ind].head;
		trees[ind].head = bran;
	}
	else
	{
		trees[ind].tail->suiv = bran;
		trees[ind].tail = trees[ind].tail->suiv;
	}

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Distance ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve and return de distance between taxa  ;
;               i and j from the delta matrix.                              ;
;                                                                           ;
; input       :                                                             ;
;               int i           : taxon i                                   ;
;               int j           : taxon j                                   ;
;               double **delta  : the delta matrix                          ;
;                                                                           ;
; return value:                                                             ;
;               double distance : dissimilarity between the two taxa        ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

double Distance (int i, int j, double **delta)
{
	if (i > j)
		return (delta[i][j]);

	else
		return (delta[j][i]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Variance ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve and return the variance of the       ;
;               distance between i and j, from the delta matrix.            ;
;                                                                           ;
; input       :                                                             ;
;               int i           : taxon i                                   ;
;               int j           : taxon j                                   ;
;               double **delta  : the delta matrix                          ;
;                                                                           ;
; return value:                                                             ;
;               double distance : the variance of Dij                       ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

double Variance (int i, int j, double **delta)
{
	if (i > j)
		return (delta[j][i]);

	else
		return (delta[i][j]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Emptied ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function check if a line is emptied or not.            ;
;                                                                           ;
; input       :                                                             ;
;               int i          : subtree (or line) i                        ;
;               double **delta : the delta matrix                           ;
;                                                                           ;
; return value:                                                             ;
;               0              : if not emptied.                            ;
;               1              : if emptied.                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

int Emptied (int i, double **delta)
{
	return ((int) delta[i][0]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Sum_S ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function retrieves the sum Sx from the diagonal       ;
;                of the delta matrix.                                       ;
;                                                                           ;
;  input       :                                                            ;
;               int i          : subtree i                                  ;
;               double **delta : the delta matrix                           ;
;                                                                           ;
;  return value:                                                            ;
;                double delta[i][i] : sum Si                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

double Sum_S (int i, double **delta)
{
	return (delta[i][i]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Compute_sums_Sx ;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function computes the sums Sx and store them in the    ;
;               diagonal of the delta matrix.                               ;
;                                                                           ;
; input       :                                                             ;
;     	         double **delta : the delta matrix.                         ;
;     	         int n          : the number of taxa                        ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void Compute_sums_Sx (double **delta, int n)
{
	int i, j;
	double sum = 0;

	for (i=1; i <= n; i++)
	{
		if (!Emptied (i, delta))
		{
			sum = 0;
			for (j=1; j<=n; j++)
			{
				if (i!=j && !Emptied (j, delta))	/* compute the sum Si */
					sum += Distance (i, j, delta);
			}
		}
		delta[i][i] = sum;		/* store the sum Si in */
	}							/* delta's diagonal    */

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Best_pair ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function finds the best pair to be agglomerated by    ;
;                minimizing the agglomerative criterion (1).                ;
;                                                                           ;
;  input       :                                                            ;
;                double **delta : the delta matrix                          ;
;                int r          : number of subtrees                        ;
;                int *a         : contain the first taxon of the pair       ;
;                int *b         : contain the second taxon of the pair      ;
;                int n          : number of taxa                            ;
;                                                                           ;
;  return value:                                                            ;
;                int *a         : the first taxon of the pair               ;
;                int *b         : the second taxon of the pair              ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void Best_pair (double **delta, int r, int *a, int *b, int n)
{
	double Qxy;			/* value of the criterion calculated */
	int x, y;			/* the pair which is tested          */
	double Qmin;		/* current minimun of the criterion  */

	Qmin=1.0e300;

	for (x=1; x<=n; x++)
	{
		if (!Emptied (x, delta))
		{
			for (y=1; y<x; y++)
			{
				if (!Emptied (y, delta))
				{
					Qxy = Agglomerative_criterion (x, y, delta, r);
					if (Qxy < Qmin - DBL_EPSILON)
					{
						Qmin = Qxy;
						*a = x;
						*b = y;
					}
				}
			}
		}
	}

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; Finish_branch_length ;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description :  Compute the length of the branch attached                 ;
;                 to the subtree i, during the final step                   ;
;                                                                           ;
;  input       :                                                            ;
;                int i          : position of subtree i                     ;
;                int j          : position of subtree j                     ;
;                int k          : position of subtree k                     ;
;                double **delta : distance matrix                           ;
;                                                                           ;
;  return value:                                                            ;
;                double length  : The length of the branch                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

double Finish_branch_length (int i, int j, int k, double **delta)
{
	double length;

	length = 0.5 * (Distance (i, j, delta) + Distance (i, k, delta) -
			Distance (j, k, delta));

	return (length);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;; FinishStr ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function compute the length of the lasts three        ;
;                subtrees and write the tree in the output string.          ;
;                                                                           ;
;  input       :                                                            ;
;                double **delta  : the delta matrix                         ;
;                int n           : the number of taxa                       ;
;                POINTERS *trees : list of subtrees                         ;
;                char *StrTree   : the output string                        ;
;                                                                           ;
;  return value:                                                            ;
;                char *StrTree   : the output string                        ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*********************************************************/

void FinishStr (double **delta, int n, POINTERS *trees, char *StrTree, const char *format)
{
	int l = 1;
	int i = 0;
	double length;
	char *tmp;
	WORD *bidon;
	WORD *ele;
	int last[3];				/* the last three subtrees     */

	while(l <= n)
	{							/* find the last tree subtree  */
		if (!Emptied (l, delta))
		{
			last[i] = l;
			i++;
		}
		l++;
	}

	tmp = (char*) mCalloc (16, sizeof (char));
	StrTree[0]='(';
	length = Finish_branch_length (last[0], last[1], last[2], delta);
	Print_outputChar (last[0], trees, StrTree);
	//snprintf (tmp, 16, "%.6f,", length);
	snprintf (tmp, INPUT_SIZE, format, length);
	if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-2)
	{
		strncat (StrTree, ":", 1);
		strncat (StrTree, tmp, strlen (tmp));
		strncat (StrTree, ",", 1);
	}

	length = Finish_branch_length (last[1], last[0], last[2], delta);
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
		length = Finish_branch_length (last[2], last[1], last[0], delta);
		Print_outputChar (last[2], trees, StrTree);
		//snprintf (tmp, 16, "%.6f", length);
		snprintf (tmp, INPUT_SIZE, format, length);
		if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-2)
		{
			strncat (StrTree, ":", 1);
			strncat (StrTree, tmp, strlen (tmp));
		}
	}

	if (strlen (StrTree) < MAX_INPUT_SIZE - 3)
		strncat (StrTree, ");\n", 3);

	if (n>3)
		n=3;
	
	for(i=0; i<n; i++)
	{
		bidon = trees[last[i]].head;
		ele = bidon;
		while (NULL != bidon)
		{
			ele = ele->suiv;
			free (bidon);
			bidon = ele;
		}
	}

	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                          Formulae                                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


/*********************************************************/

double Agglomerative_criterion (int i, int j, double **delta, int r)
{
	double Qij;											/* Formula (1) */

	Qij = (r - 2) * Distance (i, j, delta) - Sum_S (i, delta) - Sum_S (j, delta);

	return (Qij);
}

/*********************************************************/

double Branch_length (int a, int b, double **delta, int r)
{
	double length;										/* Formula (2) */

	length = (Distance (a, b, delta) + (Sum_S (a, delta) - Sum_S (b, delta)) / (r-2)) / 2;

	return (length);
}

/*********************************************************/

double Reduction4 (int a, double la, int b, double lb, int i, double lamda, double **delta)
{
	double Dui;											/* Formula (4) */

	Dui = lamda * (Distance (a, i, delta) - la) +
		(1-lamda) * (Distance (b, i, delta) - lb);
	
	return (Dui);
}

/*********************************************************/

double Reduction10 (int a, int b, int i, double lamda, double vab, double **delta)
{
	double Vci;											/* Formula (10) */

	Vci = lamda * Variance (a, i, delta) +
		(1-lamda) * Variance (b, i, delta) -
		lamda * (1-lamda) * vab;

	return (Vci);
}

/*********************************************************/

double Lamda (int a, int b, double vab, double **delta, int n, int r)
{
	int i;
	double lamda = 0.0;

	if (vab == 0.0)
		lamda = 0.5;

	else
	{
		for (i=1; i<=n; i++)
		{
			if (a != i && b != i && !Emptied (i, delta))
				lamda = lamda + (Variance (b, i, delta) - Variance (a, i, delta));
		}
		lamda = 0.5 + lamda / (2 * (r - 2) * vab);
	}													/* Formula (9) and the   */
														/* constraint that lamda */
	if (lamda > 1.0)									/* belongs to [0,1]      */
    	lamda = 1.0;

	if (lamda < 0.0)
		lamda = 0.0;

	return (lamda);
}

