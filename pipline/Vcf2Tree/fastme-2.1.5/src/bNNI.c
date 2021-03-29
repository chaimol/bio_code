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


#include "bNNI.h"


/*********************************************************/

void bNNIRetestEdge (int *p, int *q, edge *e,tree *T, double **avgDistArray,
	double *weights, int *location, int *possibleSwaps)
{
	int tloc;

	tloc = location[e->head->index+1];
	location[e->head->index+1] = bNNIEdgeTest (e, T, avgDistArray, weights + e->head->index+1);

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

void bNNItopSwitch (edge *e, int direction, double **A)
{
	edge *down, *swap, *fixed;
	node *u, *v;

	if (verbose > 2 && !isBoostrap)
	{
		if (LEFT == direction)
			Debug ( (char*)"Performing branch swap across edge '%s' with left subtree.", e->label);
		else
			Debug ( (char*)"Performing branch swap across edge '%s' with right subtree.", e->label);
	}

	down = siblingEdge (e);
	u = e->tail;
	v = e->head;

	if (LEFT == direction)
	{
		swap = e->head->leftEdge;
		fixed = e->head->rightEdge;
		v->leftEdge = down;
	}
	else
	{
		swap = e->head->rightEdge;
		fixed = e->head->leftEdge;
		v->rightEdge = down;
	}

	swap->tail = u;
	down->tail = v;

	if (e->tail->leftEdge == e)
		u->rightEdge = swap;

	else
		u->leftEdge = swap;

	bNNIupdateAverages (A, v, e->tail->parentEdge, down, swap, fixed);

	return;
}

/*********************************************************/

void bNNI (tree *T, double **avgDistArray, int *count, FILE *statfile)
{
	edge *e;
	edge **edgeArray;
	int *p, *location, *q;
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

	assignBMEWeights (T, avgDistArray);
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
		location[e->head->index+1] = bNNIEdgeTest (e, T, avgDistArray,
			weights + e->head->index + 1);
		e = depthFirstTraverse (T,e);
	}

	possibleSwaps = makeThreshHeap (p, q, weights, T->size+1,0.0);
	permInverse (p, q, T->size+1);

	/* We put the negative values of weights into a heap, indexed by p
	 * with the minimum value pointed to by p[1]
	 * p[i] is index (in edgeArray) of edge with i-th position in the
	 * heap, q[j] is the position of edge j in the heap */

	while (weights[p[1]] < -DBL_EPSILON)
	{
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

		bNNItopSwitch (edgeArray[p[1]], location[p[1]], avgDistArray);
		location[p[1]] = NONE;
		weights[p[1]] = 0.0;	//after the bNNI, this edge is in optimal configuration
		popHeap (p, q, weights, possibleSwaps--, 1);

		/* but we must retest the other edges of T */
		/* CHANGE 2/28/2003 expanding retesting to _all_ edges of T */

		e = depthFirstTraverse (T, NULL);
		while (NULL != e)
		{
			bNNIRetestEdge (p, q, e, T, avgDistArray, weights, location, &possibleSwaps);
			e = depthFirstTraverse (T, e);
		}
	}

	free (p);
	free (q);
	free (location);
	free (edgeArray);
	free (weights);
	assignBMEWeights (T, avgDistArray);

	return;
}

/* This function is the meat of the average distance matrix recalculation.
 * Idea is: we are looking at the subtree rooted at rootEdge. The subtree
 * rooted at closer is closer to rootEdge after the NNI, while the subtree
 * rooted at further is further to rootEdge after the NNI. direction tells
 * the direction of the NNI with respect to rootEdge */

/*********************************************************/

void updateSubTreeAfterNNI (double **A, node *v, edge *rootEdge,
	node *closer, node *further, double dcoeff, int direction)
{
	edge *sib;

	switch (direction)
	{
		case UP :	/* rootEdge is below the center edge of the NNI
					 * recursive calls to subtrees, if necessary */
			if (NULL != rootEdge->head->leftEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->leftEdge,
					closer, further, 0.5 * dcoeff, UP);

			if (NULL != rootEdge->head->rightEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->rightEdge,
					closer, further, 0.5 * dcoeff, UP);

			updatePair (A, rootEdge, rootEdge, closer, further, dcoeff, UP);
			sib = siblingEdge (v->parentEdge);
			A[rootEdge->head->index][v->index] =
			A[v->index][rootEdge->head->index] =
				0.5 * A[rootEdge->head->index][sib->head->index] +
				0.5 * A[rootEdge->head->index][v->parentEdge->tail->index];
			break;

		case DOWN :	/* rootEdge is above the center edge of the NNI */
			sib = siblingEdge (rootEdge);
			if (NULL != sib)
				updateSubTreeAfterNNI (A, v, sib, closer, further,
					0.5 * dcoeff, SKEW);

			if (NULL != rootEdge->tail->parentEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->tail->parentEdge,
					closer, further, 0.5 * dcoeff, DOWN);

			updatePair (A, rootEdge, rootEdge, closer, further, dcoeff, DOWN);
			A[rootEdge->head->index][v->index] =
			A[v->index][rootEdge->head->index] =
				0.5 * A[rootEdge->head->index][v->leftEdge->head->index] +
				0.5 * A[rootEdge->head->index][v->rightEdge->head->index];
			break;

		case SKEW :	/* rootEdge is in subtree skew to v */
			if (NULL != rootEdge->head->leftEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->leftEdge,
					closer, further, 0.5 * dcoeff, SKEW);

			if (NULL != rootEdge->head->rightEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->rightEdge,
					closer, further, 0.5 * dcoeff, SKEW);

			updatePair (A, rootEdge, rootEdge, closer, further, dcoeff, UP);
			A[rootEdge->head->index][v->index] =
			A[v->index][rootEdge->head->index] =
				0.5 * A[rootEdge->head->index][v->leftEdge->head->index] +
				0.5 * A[rootEdge->head->index][v->rightEdge->head->index];
			break;
	}

	return;
}

/* swapping across edge whose head is v */

/*********************************************************/

void bNNIupdateAverages (double **A, node *v, edge *par, edge *skew,
	edge *swap, edge *fixed)
{
	A[v->index][v->index] = 0.25 *
		(A[fixed->head->index][par->head->index] +
		A[fixed->head->index][swap->head->index] +
		A[skew->head->index][par->head->index] +
		A[skew->head->index][swap->head->index]);

	updateSubTreeAfterNNI (A, v, fixed, skew->head, swap->head, 0.25, UP);
	updateSubTreeAfterNNI (A, v, par, swap->head, skew->head, 0.25, DOWN);
	updateSubTreeAfterNNI (A, v, skew, fixed->head, par->head, 0.25, UP);
	updateSubTreeAfterNNI (A, v, swap, par->head, fixed->head, 0.25, SKEW);

	return;
}

/*********************************************************/

double wf5 (double D_AD, double D_BC, double D_AC, double D_BD,
	double D_AB, double D_CD)
{
	double weight;

	weight = 0.25 * (D_AC + D_BD + D_AD + D_BC) + 0.5 * (D_AB + D_CD);

	return (weight);
}

/*********************************************************/

int bNNIEdgeTest (edge *e, tree *T, double **A, double *weight)
{
	edge *f;
	double D_LR, D_LU, D_LD, D_RD, D_RU, D_DU;
	double w1, w2, w0;

	if ((leaf(e->tail)) || (leaf(e->head)))
		return (NONE);

	f = siblingEdge (e);

	D_LR = A[e->head->leftEdge->head->index][e->head->rightEdge->head->index];
	D_LU = A[e->head->leftEdge->head->index][e->tail->index];
	D_LD = A[e->head->leftEdge->head->index][f->head->index];
	D_RU = A[e->head->rightEdge->head->index][e->tail->index];
	D_RD = A[e->head->rightEdge->head->index][f->head->index];
	D_DU = A[e->tail->index][f->head->index];

	w0 = wf5 (D_RU, D_LD, D_LU, D_RD, D_DU, D_LR);	// weight of current config
	w1 = wf5 (D_RU, D_LD, D_DU, D_LR, D_LU, D_RD);	// weight with L<->D switch
	w2 = wf5 (D_DU, D_LR, D_LU, D_RD, D_RU, D_LD);	// weight with R<->D switch

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
		return(RIGHT);
	}
	else				// w1 < w2, w0
	{
		*weight = w1 - w0;
		if (verbose > 2 && !isBoostrap)
		{
			Debug ( (char*)"Possible swap across '%s'. Weight dropping by %f.", e->label, w0 - w1);
			Debug ( (char*)"New tree length should be %f.", T->weight + w1 - w0);
		}
		return(LEFT);
	}
}

/* limitedFillTableUp fills all the entries in D associated with
 * e->head, f->head and those edges g->head above e->head, working
 * recursively and stopping when trigger is reached */

/*********************************************************/

void limitedFillTableUp (edge *e, edge *f, double **A, edge *trigger)
{
	edge *g, *h;
	g = f->tail->parentEdge;
	if (f != trigger)
		limitedFillTableUp (e, g, A, trigger);

	h = siblingEdge (f);
	A[e->head->index][f->head->index] =
	A[f->head->index][e->head->index] =
	0.5 * (A[e->head->index][g->head->index] + A[e->head->index][h->head->index]);

	return;
}

