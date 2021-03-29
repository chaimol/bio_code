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


#include "inputs.h"

/*********************************************************/

void compareSets (tree *T, set *S)
{
	edge *e;
	node *v,*w;
	set *X;

	e = depthFirstTraverse (T, NULL);
	while (NULL != e)
	{
		v = e->head;
		for (X=S; NULL!=X; X=X->secondNode)
		{
			w = X->firstNode;
			if (0 == strcmp (v->label, w->label))
			{
				v->index2 = w->index2;
				w->index2 = -1;
				break;
			}
		}
		e = depthFirstTraverse (T, e);
	}

	v = T->root;
	for (X=S; NULL!=X; X=X->secondNode)
	{
		w = X->firstNode;
		if (0 == strcmp (v->label, w->label))
		{
			v->index2 = w->index2;
			w->index2 = -1;
			break;
		}
	}
	if (-1 == v->index2)
		Exit ( (char*)"Leaf (1) '%s' in tree not in distance matrix.", v->label);

	e = depthFirstTraverse (T, NULL);
	while (NULL != e)
	{
		v = e->head;
		if ((leaf(v)) && (-1 == v->index2))
			Exit ( (char*)"Leaf (2) '%s' in tree not in distance matrix.", v->label);

		e = depthFirstTraverse (T, e);
	}

	for (X=S; NULL!=X; X=X->secondNode)
		if (X->firstNode->index2 > -1)
			Exit ( (char*)"Node '%s' in matrix but not a leaf in tree.", X->firstNode->label);

  return;
}

/*********************************************************/

void freeCharMatrix (char **D, int size)
{
	int i;

	for(i=0; i<size; i++)
	{
		if (NULL != D[i])
			free (D[i]);
	}

	if (NULL != D)
		free (D);

	return;
}

/*********************************************************/

void freeMatrix (double **D, int size)
{
	int i=0;

	if (NULL != D)
	{
		for (i=0; i<size; i++)
			if (NULL != D[i])
				free(D[i]);
		free (D);
	}

	return;
}

/*********************************************************/

boolean isTrgMatrix (FILE *ifile)
{
	char nextString[MAX_NAME_LENGTH];
	fpos_t pos;
	boolean ret = false;
	int c;

	// record the FILE stream position
	if (fgetpos (ifile, &pos) != 0)
		Exit ( (char*)"Cannot access matrix file.");

	// read first label
	if (!(fscanf (ifile, "%s", nextString)))
		Exit ( (char*)"Cannot read matrix file.");

	// read first matrix value
	if (!(fscanf (ifile, "%s", nextString)))
		Exit ( (char*)"Cannot read matrix values.");

	// read end of line
	c = fgetc (ifile);
	// This is the Linux end of line
	if ('\n' == c)
		ret = true;
	// This is the Windows end of line
	else if ('\r'  == c)
	{
		if ('\n' == fgetc (ifile))
			ret = true;
	}

	// rewind the FILE stream position to the previously recorded
	if (fsetpos (ifile, &pos) != 0)
		Exit ( (char*)"Cannot move through matrix file.");

	return ret;
}

/*********************************************************/

double **loadM (FILE *ifile, int *size, set *S)
{
	double **D, val;
	int i, j, zsize;
	char *line = NULL;
	char *tok;
	char **tokens;
	node *v = NULL;
	
	// Read first line to get matrix size
	char fline[64];
	if (!(fscanf (ifile, "%s", fline)))
		Exit ( (char*)"Cannot load input matrix.");
	// Read until end of line
	char c;
	do
	{
		c = (char)fgetc (ifile);
	}
	while (c != '\n' && c != EOF);
	if (!isNumeric (fline))
		Exit ( (char*)"Expecting an integer value on first line of input matrix (number of taxa).");
		
	zsize = atoi(fline);

	if ((zsize < 0) || (zsize > MAXSIZE))
		Exit ( (char*)"Number of taxa is out of bounds.");
		
	D = (double **) mCalloc (zsize, sizeof (double *));
	for (i=0; i<zsize; i++)
	{
		D[i] = (double *) mCalloc (zsize, sizeof (double));
		memset (D[i], 0, (unsigned long) zsize * sizeof(double));
	}

	// Read matrix
	for (i=0; i<zsize; i++)
	{
		// Read line, remove trailing EOL, replace tabulation by blank space
		// line length = MAX_NAME_LENGTH + 16 blank spaces + ( N taxa * ([DECIMAL_DIG] digits + 2 blank spaces) )
		if (NULL != (line = str_replace ( str_replace ( str_replace ( getLine (ifile, line, (MAX_NAME_LENGTH + 16 + (zsize * (DECIMAL_DIG+2)))), '\n', " " ), '\r', " " ), '\t', " " )))
		{
			// Split line on blank spaces
			tokens = str_split (line, ' ');
			if (tokens)
			{
				// Read each token
				for (j=0; j<=zsize; j++)
				{
					if (NULL != tokens[j])
					{
						tok = tokens[j];
						// First token is the sequence name
						if (j==0)
						{
							if (strlen (tok) > MAX_NAME_LENGTH)
								Exit ( (char*)"Taxa name length is limited to %d char.", MAX_NAME_LENGTH);
							
							if (checkLabelExist (S, tok))
								Exit ( (char*)"Duplicated taxon name: '%s'.", tok);

							v = makeNode (tok, -1);
							v->index2 = i;
							S = addToSet (v, S);
						}
						// Other tokens are distance values
						else
						{
							if (strlen (tok) > DECIMAL_DIG)
								Exit ( (char*)"Distance precision must not exceed %s digits.", DECIMAL_DIG);
							
							if (! isNumeric (tok))
								Exit ( (char*)"Invalid distance matrix : numerical value expected for taxon '%s' instead of '%s'.",
									v->label, tok);
						
							val = atof (tok);
						
							if (0 > val)
								Exit ( (char*)"Distance matrix expected : input of %s off diagonal is inappropriate.", tok);

							D[j-1][i] = val;
						}
						free (tokens[j]);
					}
					else
						Exit ( (char*)"Invalid matrix format.");
				}
			}
			free(tokens);
		}
		else
			Exit ( (char*)"Cannot read matrix line %d.", i);
	}
	
	if (NULL != line)
		free (line);
	
	*size = zsize;
	
	return D;
}

/*********************************************************/

double **loadMatrix (FILE *ifile, int *size, set *S)
{
	char nextString[MAX_NAME_LENGTH];
	double **table;
	int zsize;

	if (!(fscanf (ifile, "%s", nextString)))
		Exit ( (char*)"Cannot load input matrix.");

	if (!isNumeric (nextString))
		Exit ( (char*)"Expecting an integer value on first line of input matrix (number of taxa).");

	zsize = atoi(nextString);

	if ((zsize < 0) || (zsize > MAXSIZE))
		Exit ( (char*)"Number of taxa is out of bounds.");

	table = (double **) mCalloc (zsize, sizeof (double *));

	// is input matrix triangular or square ?
	if (isTrgMatrix (ifile))
		loadTrgMatrix (ifile, zsize, S, table);

	else
		loadRectMatrix (ifile, zsize, S, table);

	*size = zsize;

	return (table);
}

/*********************************************************/

void loadTrgMatrix (FILE *ifile, int size, set *S, double **table)
{
	char nextString[MAX_NAME_LENGTH];
	node *v;
	double val;
	int i, j;

	for (i=0; i<size; i++)
	{
		j = 0;
		table[i] = (double *) mCalloc (size, sizeof (double));
		if (!(fscanf (ifile, "%s", nextString)))
			Exit ( (char*)"Cannot load label %d.", i);

		v = makeNode (nextString, -1);
		v->index2 = i;
		S = addToSet (v, S);

		while (j < i+1)
		{
			if (!(fscanf (ifile, "%s", nextString)))
				Exit ( (char*)"Cannot load (%d,%d)-entry.", i, j);

			if ((nextString[0] < '0' || nextString[0] > '9') && nextString[0] != '.')
				Exit ( (char*)"Invalid distance matrix : numerical value expected for taxon '%s' instead of '%s'.",
					v->label, nextString);

			val = atof (nextString);
			if ((i != j) && (0 > val))
				Exit ( (char*)"Distance matrix expected : input of %s off diagonal is inappropriate.", nextString);

			table[j][i] = val;
			table[i][j++] = val;
		}
	}

	return;
}

/*********************************************************/

void loadRectMatrix (FILE *ifile, int size, set *S, double **table)
{
	char nextString[MAX_NAME_LENGTH];
	node *v;
	double val;
	int i, j;

	for (i=0; i<size; i++)
	{
		j = 0;
		table[i] = (double *) mCalloc (size, sizeof (double));
		if (!(fscanf (ifile, "%s", nextString)))
			Exit ( (char*)"Cannot load label %d.", i);

		v = makeNode (nextString, -1);
		v->index2 = i;
		S = addToSet (v, S);

		while (j < size)
		{
			if (!(fscanf (ifile, "%s", nextString)))
				Exit ( (char*)"Cannot load (%d,%d)-entry.", i, j);

			if ((nextString[0] < '0' || nextString[0] > '9') && nextString[0] != '.')
				Exit ( (char*)"Invalid distance matrix : numerical value expected for taxon '%s' instead of '%s'.",
					v->label, nextString);

			val = atof (nextString);
			if ((i != j) && (0 > val))
				Exit ( (char*)"Distance matrix expected : input of %s off diagonal is inappropriate.", nextString);

			table[i][j++] = val;
		}
	}
	return;
}

/*********************************************************/

void partitionSizes (tree *T)
{
	edge *e;
	e = depthFirstTraverse (T, NULL);

	while (NULL != e)
	{
		if (leaf (e->head))
			e->bottomsize = 1;
		else
			e->bottomsize = e->head->leftEdge->bottomsize +
						e->head->rightEdge->bottomsize;
		e->topsize = (T->size + 2) / 2 - e->bottomsize;
		e = depthFirstTraverse (T, e);
	}

	return;
}

/*********************************************************/

double **loadScoreMatrix (int *d, FILE *ifile, char *alphabet)
{
	char nextString[MAX_NAME_LENGTH];
	double **table;
	int i, j;

	if (!(fscanf (ifile, "%s", nextString)))
		Exit ( (char*)"Cannot load size of score matrix.");

	*d = atoi (nextString);
	if ((*d < 1) || (*d > PROTEIN_ALPHABET_SIZE))
		Exit ( (char*)"Invalid size of score matrix.");

	table = initDoubleMatrix (*d);
	/* *d is actual size of input matrix, which may be less than
	 * PROTEIN_ALPHABET_SIZE, since the ALPHABET contains special characters */

	for (i=0; i<*d; i++)
	{
		j = 0;
		if (!(fscanf (ifile, "%s", nextString)))
			Exit ( (char*)"Cannot load label %d.", i);

		alphabet[i] = nextString[0];
		while (j < *d)
		{
			if (!(fscanf (ifile, "%s", nextString)))
				Exit ( (char*)"Cannot load (%d,%d)-entry.", i, j);

			table[i][j++] = atof (nextString);
		}	/* j loop */
	}	/* i loop */

	if (*d < PROTEIN_ALPHABET_SIZE)
		alphabet[*d] = '\0';

	return (table);
}

