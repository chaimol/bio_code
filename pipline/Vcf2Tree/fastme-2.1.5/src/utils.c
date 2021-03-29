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


#include "utils.h"


/*********************************************************/

void constantToStr (int c, char *str)
{
	switch (c)
	{
		case TaxAddBAL:
			strncpy (str, "TaxAdd_BalME", MAX_NAME_LENGTH);
			break;
		case TaxAddOLS:
			strncpy (str, "TaxAdd_OLSME", MAX_NAME_LENGTH);
			break;
		case BALNNI:
			strncpy (str, "NNI_BalME", MAX_NAME_LENGTH);
			break;
		case OLSNNI:
			strncpy (str, "NNI_OLSME", MAX_NAME_LENGTH);
			break;
		case NJ:
			strncpy (str, "NJ", MAX_NAME_LENGTH);
			break;
		case UNJ:
			strncpy (str, "UNJ", MAX_NAME_LENGTH);
			break;
		case BIONJ:
			strncpy (str, "BIONJ", MAX_NAME_LENGTH);
			break;
		case BrBAL:
			strncpy (str, "BalLS", MAX_NAME_LENGTH);
			break;
		case BrOLS:
			strncpy (str, "OLS", MAX_NAME_LENGTH);
			break;
		case USER:
			strncpy (str, "User", MAX_NAME_LENGTH);
			break;
		case NONE:
			strncpy (str, "None", MAX_NAME_LENGTH);
			break;
			
		case PDIST:
			strncpy (str, "p-distance", MAX_NAME_LENGTH);
			break;
		case RY:
			strncpy (str, "RY", MAX_NAME_LENGTH);
			break;
		case RYSYM:
			strncpy (str, "RY symetric", MAX_NAME_LENGTH);
			break;
		case F81:
			strncpy (str, "F81", MAX_NAME_LENGTH);
			break;
		case F84:
			strncpy (str, "F84", MAX_NAME_LENGTH);
			break;
		case TN93:
			strncpy (str, "TN93", MAX_NAME_LENGTH);
			break;
		case K2P:
			strncpy (str, "K2P", MAX_NAME_LENGTH);
			break;
		case JC69:
			strncpy (str, "JC69", MAX_NAME_LENGTH);
			break;
		case LOGDET:
			strncpy (str, "LogDet", MAX_NAME_LENGTH);
			break;
		case PROTEIN:
			strncpy (str, "Protein", MAX_NAME_LENGTH);
			break;
		case MATRIX:
			strncpy (str, "Matrix", MAX_NAME_LENGTH);
			break;
		case DNA:
			strncpy (str, "DNA", MAX_NAME_LENGTH);
			break;
//		case SCOREDIST:
//			strncpy (str, "Scoring matrix", MAX_NAME_LENGTH);
//			break;
		case WAG:
			strncpy (str, "WAG", MAX_NAME_LENGTH);
			break;
		case DAYHOFF:
			strncpy (str, "Dayhoff", MAX_NAME_LENGTH);
			break;
		case JTT:
			strncpy (str, "JTT", MAX_NAME_LENGTH);
			break;
		case BLOSUM62:
			strncpy (str, "BLOSUM62", MAX_NAME_LENGTH);
			break;
		case MTREV:
			strncpy (str, "MtREV", MAX_NAME_LENGTH);
			break;
		case RTREV:
			strncpy (str, "RtREV", MAX_NAME_LENGTH);
			break;
		case CPREV:
			strncpy (str, "CpREV", MAX_NAME_LENGTH);
			break;
		case DCMUT:
			strncpy (str, "DCMut", MAX_NAME_LENGTH);
			break;
		case VT:
			strncpy (str, "VT", MAX_NAME_LENGTH);
			break;
		case LG:
			strncpy (str, "LG", MAX_NAME_LENGTH);
			break;
		case F81LIKE:
			strncpy (str, "F81-like", MAX_NAME_LENGTH);
			break;
		case HIVB:
			strncpy (str, "HIVb", MAX_NAME_LENGTH);
			break;
		case HIVW:
			strncpy (str, "HIVw", MAX_NAME_LENGTH);
			break;
		case FLU:
			strncpy (str, "FLU", MAX_NAME_LENGTH);
			break;
		default:
			break;
	}
	return;
}

/*********************************************************/

int *initZeroArray (int l)
{
	int *x;
	int i;

	x = (int *) mCalloc (l, sizeof (int));

	for (i=0; i<l; i++)
		x[i] = 0;

	return (x);
}

/*********************************************************/

int *initOneArray (int l)
{
	int *x;
	int i;

	x = (int *) mCalloc (l, sizeof (int));

	for (i=0; i<l; i++)
		x[i] = 1;

	return (x);
}

/*********************************************************/

double **initDoubleMatrix (int d)
{
	int i,j;
	double **A;

	A = (double **) mCalloc (d, sizeof (double *));

	for (i=0; i<d; i++)
	{
		A[i] = (double *) mCalloc (d, sizeof (double));
		for (j=0; j<d; j++)
			A[i][j] = 0.0;
	}

	return (A);
}

/*********************************************************/

void fillZeroMatrix (double ***A, int d)
{
	int i,j;
	
	for (i=0; i<d; i++)
	{
		for (j=0; j<d; j++)
			(*A)[i][j] = 0.0;
	}
	
	return;
}

/*********************************************************/

boolean whiteSpace (char c)
{
	if ((' ' == c) || ('\t' == c) || ('\n' == c))
		return (TRUE);

	else
		return (FALSE);
}

/*********************************************************/

void Uppercase (char *str)
{
	char ch;
	int i;

	for (i=0; i< (int) strlen (str); i++)
	{
		ch = (char) toupper (str[i]);
		str[i] = ch;
	}

	return;
}

/*********************************************************/

void Exit (char *message, ...)
{
	va_list ptr;

	fprintf (stderr, "\n . Error: ");
	va_start (ptr, message);
	vfprintf (stderr, message, ptr);
	va_end (ptr);
	fprintf (stderr, "\n");

	fflush (NULL);

	exit (EXIT_FAILURE);

	//return;
}

/*********************************************************/

void Warning (char *message, ...)
{
	va_list ptr;

	printf ("\n . Warning: ");
	va_start (ptr, message);
	vprintf (message, ptr);
	va_end (ptr);
	printf ("\n");

	fflush (NULL);

	return;
}

/*********************************************************/

void Message (char *message, ...)
{
	va_list ptr;

	printf ("\n . ");
	va_start (ptr, message);
	vprintf (message, ptr);
	va_end (ptr);
	printf ("\n");

	fflush (NULL);

	return;
}

/*********************************************************/

void Debug (char *message, ...)
{
	va_list ptr;

	printf ("\n ... ");
	va_start (ptr, message);
	vprintf (message, ptr);
	va_end (ptr);

	fflush (NULL);

	return;
}

/*********************************************************/

void *mCalloc (int nb, size_t size)
{
	void *allocated = NULL;

#ifdef _OPENMP
	#pragma omp critical (memAlloc)
	{
#endif

	if ((allocated = calloc ((size_t) nb, (size_t) size)) == NULL)
		Exit ( (char*)"Low memory! nb %d size %d", nb, size);

#ifdef _OPENMP
	}
#endif

	return allocated;
}

/*********************************************************/

boolean isNumeric (const char *p)
{
	if (*p)
	{
		char c;
		while ((c=*p++))
		{
			if (!isdigit(c) && '.' != c && ',' != c)
				return false;
		}
		return true;
	}

	return false;
}

/*********************************************************/

char *getLine (FILE *file, char *line, const int len)
{
	line = (char*) mCalloc (len, sizeof(char));
    
    memset (line, '\0', (size_t) len);
    
    if (! fgets (line, len, file))
    {
		Exit ( (char*)"Cannot read line.");
	}

    return line;
}

/*********************************************************/

char *str_replace (const char *s, char ch, const char *repl)
{
	int count = 0;
	const char *t;
	for (t=s; *t; t++)
		count += (*t == ch);

	size_t rlen = strlen(repl);
	char *res = (char*) mCalloc ( (int) strlen(s) + ( ( (int)rlen - 1) * count) + 1, sizeof(char));
	char *ptr = res;
	for(t=s; *t; t++)
	{
		if(*t == ch)
		{
			memcpy(ptr, repl, rlen);
			ptr += rlen;
		}
		else
		{
			*ptr++ = *t;
		}
	}
	*ptr = 0;
	
	return res;
}

/*********************************************************/

int countFields (char *str, const char c)
{
	int i, count, len;
	count = 0;
	
	len = (int) strlen(str);
	
	// increment 'count' each time a 'c' (separator) is encountered if previous char was not 'c'
	for (i=1; i<len; i++)
	{
		if (str[i] == c && str[i-1] != c)
			count++;
	}
	
	// decrement 'count' if last char of the string is 'c'
	// (i.e. cannot be considered as a separator)
	if (str[len-1] == c)
		count--;
	
	// the number of fields equals the number of separators +1
	count++;
	
	return count;
}

/*********************************************************/

char **str_split (char *str, const char delim)
{
    char** result = NULL;
    int i, count;
    char d[2];
    d[0] = delim;
    d[1] = '\0';

    /* Count how many elements will be extracted. */
    count = countFields (str, delim);

    result = (char**) mCalloc (count, sizeof(char*));
    if (result)
    {
        i = 0;
        char* token = strtok (str, d);

        while (token && i < count)
        {
            //assert(i < count);
            result[i] = (char*) mCalloc ( (int) strlen(token) +1, sizeof(char));
            strncpy (result[i++], token, strlen(token));
            token = strtok (NULL, d);
        }
    }

    return result;
}


