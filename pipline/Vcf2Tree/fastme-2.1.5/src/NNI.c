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


#include "NNI.h"

/* NNI functions for unweighted OLS topological switches */

/*********************************************************/

double **buildAveragesTable (tree *T, double **D)
{
	int i, j;
	double **A;

	A = (double **) mCalloc (T->size, sizeof(double *));

	for (i = 0; i < T->size; i++)
	{
		A[i] = (double *) mCalloc (T->size, sizeof(double));
		for (j=0; j<T->size; j++)
			A[i][j] = 0.0;
	}

	makeOLSAveragesTable (T, D, A);

	return (A);
}

/*********************************************************/

double wf2 (double lambda, double D_AD, double D_BC, double D_AC,
		double D_BD, double D_AB, double D_CD)
{
	double weight;

	weight = 0.5 * (lambda * (D_AC + D_BD) + (1 - lambda) * (D_AD + D_BC) + (D_AB + D_CD));

	return (weight);
}

/*********************************************************/

int NNIEdgeTest (edge *e, tree *T, double **A, double *weight)
{
	int a, b, c, d;
	edge *f;
	double *lambda;
	double D_LR, D_LU, D_LD, D_RD, D_RU, D_DU;
	double w1, w2, w0;

	if ((leaf(e->tail)) || (leaf(e->head)))
		return (NONE);

	lambda = (double *) mCalloc (3, sizeof(double));
	a = e->tail->parentEdge->topsize;
	f = siblingEdge (e);
	b = f->bottomsize;
	c = e->head->leftEdge->bottomsize;
	d = e->head->rightEdge->bottomsize;

	lambda[0] = ((double) b*c + a*d) / ((a + b) * (c+d));
	lambda[1] = ((double) b*c + a*d) / ((a + c) * (b+d));
	lambda[2] = ((double) c*d + a*b) / ((a + d) * (b+c));

	D_LR = A[e->head->leftEdge->head->index][e->head->rightEdge->head->index];
	D_LU = A[e->head->leftEdge->head->index][e->tail->index];
	D_LD = A[e->head->leftEdge->head->index][f->head->index];
	D_RU = A[e->head->rightEdge->head->index][e->tail->index];
	D_RD = A[e->head->rightEdge->head->index][f->head->index];
	D_DU = A[e->tail->index][f->head->index];

	w0 = wf2 (lambda[0], D_RU, D_LD, D_LU, D_RD, D_DU, D_LR);
	w1 = wf2 (lambda[1], D_RU, D_LD, D_DU, D_LR, D_LU, D_RD);
	w2 = wf2 (lambda[2], D_DU, D_LR, D_LU, D_RD, D_RU, D_LD);

	free(lambda);

	if (w0 <= w1)
	{
		if (w0 <= w2)	// w0 <= w1,w2
		{
			*weight = 0.0;
			return (NONE);
		}
		else			// w2 < w0 <= w1
		{
			*weight = w2 - w0;
			if (verbose > 2 && !isBoostrap)
			{
				Debug ( (char*)"Possible swap across '%s'. Weight dropping by %f.", e->label, w0 - w2);
				Debug ( (char*)"New tree length should be %f.", T->weight + w2 - w0);
			}
			return (RIGHT);
		}
	}
	else if (w2 <= w1)	// w2 <= w1 < w0
	{
		*weight = w2 - w0;
		if (verbose > 2 && !isBoostrap)
		{
			Debug ( (char*)"Possible swap across '%s'. Weight dropping by %f.", e->label, w0 - w2);
			Debug ( (char*)"New tree length should be %f.", T->weight + w2 - w0);
		}
		return (RIGHT);
	}
	else				// w1 < w2, w0
	{
		*weight = w1 - w0;
		if (verbose > 2 && !isBoostrap)
		{
			Debug ( (char*)"Possible swap across '%s'. Weight dropping by %f.", e->label, w0 - w1);
			Debug ( (char*)"New tree length should be %f.", T->weight + w1 - w0);
		}
		return (LEFT);
	}
}

/*********************************************************/

void NNIupdateAverages (double **A, edge *e, edge *par, edge *skew,
	edge *swap, edge *fixed, tree *T)
{
	node *v;
	edge *elooper;
	v = e->head;

	/* first, v */

	A[e->head->index][e->head->index] = (swap->bottomsize *
			((skew->bottomsize*A[skew->head->index][swap->head->index] +
			fixed->bottomsize*A[fixed->head->index][swap->head->index]) /
			e->bottomsize) + par->topsize *
			((skew->bottomsize*A[skew->head->index][par->head->index] +
			fixed->bottomsize*A[fixed->head->index][par->head->index]) /
			e->bottomsize)) / e->topsize;

	elooper = findBottomLeft (e);

	/* next we loop over all the edges which are below e */

	while (e != elooper)
	{
		A[e->head->index][elooper->head->index] =
			A[elooper->head->index][v->index] =
			(swap->bottomsize*A[elooper->head->index][swap->head->index] +
			par->topsize*A[elooper->head->index][par->head->index]) / e->topsize;

		elooper = depthFirstTraverse (T, elooper);
	}
	elooper = findBottomLeft (swap);

	/* next we loop over all the edges below and including swap */

	while (swap != elooper)
	{
		A[e->head->index][elooper->head->index] =
			A[elooper->head->index][e->head->index] =
			(skew->bottomsize * A[elooper->head->index][skew->head->index] +
			fixed->bottomsize*A[elooper->head->index][fixed->head->index]) / e->bottomsize;

		elooper = depthFirstTraverse (T, elooper);
	}

	/* now elooper = skew */

	A[e->head->index][elooper->head->index] =
		A[elooper->head->index][e->head->index] =
		(skew->bottomsize * A[elooper->head->index][skew->head->index] +
		fixed->bottomsize* A[elooper->head->index][fixed->head->index]) / e->bottomsize;

	/* finally, we loop over all the edges in the tree on the far side of parEdge */

	elooper = T->root->leftEdge;
	while ((elooper != swap) && (elooper != e))	// start a top-first traversal
	{
		A[e->head->index][elooper->head->index] =
			A[elooper->head->index][e->head->index] =
			(skew->bottomsize * A[elooper->head->index][skew->head->index] +
			fixed->bottomsize* A[elooper->head->index][fixed->head->index]) / e->bottomsize;

		elooper = topFirstTraverse (T, elooper);
	}

	/* At this point, elooper = par.
	 * We finish the top-first traversal, excluding the subtree below par */

	elooper = moveUpRight (par);

	while (NULL != elooper)
	{
		A[e->head->index][elooper->head->index] =
			A[elooper->head->index][e->head->index] =
			(skew->bottomsize * A[elooper->head->index][skew->head->index] +
			fixed->bottomsize* A[elooper->head->index][fixed->head->index]) / e->bottomsize;

		elooper = topFirstTraverse (T, elooper);
	}

	return;
}

/*********************************************************/

void NNItopSwitch (tree *T, edge *e, int direction, double **A)
{
	edge *par, *fixed;
	edge *skew, *swap;

	if (verbose > 2 && !isBoostrap)
	{
		if (LEFT == direction)
			Debug ( (char*)"Performing branch swap across edge '%s' with left subtree.", e->label);

		else
			Debug ( (char*)"Performing branch swap across edge '%s' with right subtree.", e->label);
	}

	if (LEFT == direction)
		swap = e->head->leftEdge;

	else
		swap = e->head->rightEdge;

	skew = siblingEdge (e);
	fixed = siblingEdge (swap);
	par = e->tail->parentEdge;

	/* Perform topological switch by changing f from (u,b) to (v,b)
	 * and g from (v,c) to (u,c), necessitatates also changing parent fields */

	swap->tail = e->tail;
	skew->tail = e->head;

	if (LEFT == direction)
		e->head->leftEdge = skew;

	else
		e->head->rightEdge = skew;

	if (skew == e->tail->rightEdge)
		e->tail->rightEdge = swap;

	else
		e->tail->leftEdge = swap;

	/* Both topsize and bottomsize change for e, but nowhere else */

	e->topsize = par->topsize + swap->bottomsize;
	e->bottomsize = fixed->bottomsize + skew->bottomsize;
	NNIupdateAverages (A, e, par, skew, swap, fixed,T);

	return;
}

/*********************************************************/

void NNIRetestEdge (int *p, int *q, edge *e,tree *T, double **avgDistArray,
	double *weights, int *location, int *possibleSwaps)
{
	int tloc;

	tloc = location[e->head->index+1];
	location[e->head->index+1] = NNIEdgeTest (e, T, avgDistArray, weights + e->head->index+1);

	if (NONE == location[e->head->index+1])
	{
		if (NONE != tloc)
			popHeap (p, q, weights, (*possibleSwaps)--, q[e->head->index+1]);
	}
	else
	{
		if (NONE == tloc)
			pushHeap (p, q, weights, (*possibleSwaps)++, q[e->head->index+1]);

		else
			reHeapElement (p, q, weights, *possibleSwaps, q[e->head->index+1]);
	}

	return;
}

/*********************************************************/

void NNI (tree *T, double **avgDistArray, int *count, FILE *statfile)
{
	edge *e, *centerEdge;
	edge **edgeArray;
	int *location;
	int *p,*q;
	int i;
	int possibleSwaps;
	double *weights;

	p = initPerm (T->size+1);
	q = initPerm (T->size+1);
	edgeArray = (edge **) mCalloc ((T->size+1), sizeof (edge *));
	weights = (double *) mCalloc ((T->size+1), sizeof (double));
	location = (int *) mCalloc ((T->size+1), sizeof (int));

	for (i=0; i<T->size+1; i++)
	{
		weights[i] = 0.0;
		location[i] = NONE;
	}

	assignOLSWeights (T, avgDistArray);
	weighTree (T);

	if (!isBoostrap)
	{
		fprintf (statfile, "\tBefore NNI:     tree length is %f.\n", T->weight);
		if (verbose > 2)
			Debug ( (char*)"Before NNI: tree length is %f.", T->weight);
		else if (verbose > 1)
			Message ( (char*)". Before NNI: tree length is %f.", T->weight);
	}

	e = findBottomLeft (T->root->leftEdge);

	while (NULL != e)
	{
		edgeArray[e->head->index+1] = e;
		location[e->head->index+1] = NNIEdgeTest (e, T, avgDistArray,
										weights + e->head->index + 1);
		e = depthFirstTraverse (T,e);
	}

	possibleSwaps = makeThreshHeap (p, q, weights, T->size+1, 0.0);
	permInverse (p, q, T->size+1);

	/* We put the negative values of weights into a heap, indexed by p
	 * with the minimum value pointed to by p[1]
	 * p[i] is index (in edgeArray) of edge with i-th position in the
	 * heap, q[j] is the position of edge j in the heap */

	while (weights[p[1]] < -DBL_EPSILON)
	{
		centerEdge = edgeArray[p[1]];
		(*count)++;
		T->weight = T->weight + weights[p[1]];
		if (!isBoostrap)
		{
			fprintf (statfile, "\tNNI  %5d: new tree length is %f.\n", *count, T->weight);
			if (verbose > 2)
				Debug ( (char*)"NNI %5d: new tree length is %f.", *count, T->weight);
			else if (verbose > 1)
				Message ( (char*)". NNI %5d: new tree length is %f.", *count, T->weight);
		}

		NNItopSwitch (T, edgeArray[p[1]], location[p[1]], avgDistArray);
		location[p[1]] = NONE;
		weights[p[1]] = 0.0;	// after the NNI, this edge is in optimal configuration
		popHeap (p, q, weights, possibleSwaps--, 1);

		// but we must retest the other four edges

		e = centerEdge->head->leftEdge;
		NNIRetestEdge (p, q, e, T, avgDistArray, weights, location, &possibleSwaps);
		e = centerEdge->head->rightEdge;
		NNIRetestEdge (p, q, e, T, avgDistArray, weights, location, &possibleSwaps);
		e = siblingEdge (centerEdge);
		NNIRetestEdge (p, q, e, T, avgDistArray, weights, location, &possibleSwaps);
		e = centerEdge->tail->parentEdge;
		NNIRetestEdge (p, q, e, T, avgDistArray, weights, location, &possibleSwaps);
	}

	free (p);
	free (q);
	free (location);
	free (edgeArray);
	free (weights);

	return;
}

