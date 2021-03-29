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


#include "gme.h"

/* fillTableUp fills all the entries in D associated with
 * e->head,f->head and those edges g->head above e->head */

/*********************************************************/

void fillTableUp (edge *e, edge *f, double **A, double **D, tree *T)
{
	edge *g,*h;

	if (T->root == f->tail)
	{
		if (leaf (e->head))
			A[e->head->index][f->head->index] =
			A[f->head->index][e->head->index] =
			D[e->head->index2][f->tail->index2];

		else
		{
			g = e->head->leftEdge;
			h = e->head->rightEdge;
			A[e->head->index][f->head->index] =
			A[f->head->index][e->head->index] =
				(g->bottomsize*A[f->head->index][g->head->index] +
				h->bottomsize*A[f->head->index][h->head->index]) / e->bottomsize;
		}
	}
	else
	{
		g = f->tail->parentEdge;
		fillTableUp (e, g, A, D, T);	// recursive call
		h = siblingEdge (f);
		A[e->head->index][f->head->index] =
		A[f->head->index][e->head->index] =
			(g->topsize*A[e->head->index][g->head->index] +
			h->bottomsize*A[e->head->index][h->head->index]) / f->topsize;
	}

	return;
}

/* OLSint and OLSext use the average distance array to calculate weights
 * instead of using the edge average weight fields */

/*********************************************************/

void OLSext (edge *e, double **A)
{
	edge *f, *g;

	if (leaf (e->head))
	{
		f = siblingEdge (e);
		e->distance = 0.5 * (A[e->head->index][e->tail->index] +
				A[e->head->index][f->head->index] -
				A[f->head->index][e->tail->index]);
	}
	else
	{
		f = e->head->leftEdge;
		g = e->head->rightEdge;
		e->distance = 0.5 * (A[e->head->index][f->head->index] +
				A[e->head->index][g->head->index] -
				A[f->head->index][g->head->index]);
	}

	return;
}

/*********************************************************/

double wf (double lambda, double D_LR, double D_LU, double D_LD,
	double D_RU, double D_RD, double D_DU)
{
	double weight;

	weight = 0.5 * (lambda * (D_LU  + D_RD) +
		(1-lambda) * (D_LD + D_RU) -
		(D_LR + D_DU));

	return (weight);
}

/*********************************************************/

void OLSint (edge *e, double **A)
{
	double lambda;
	edge *left, *right, *sib;

	left = e->head->leftEdge;
	right = e->head->rightEdge;
	sib = siblingEdge (e);

	lambda = ((double) sib->bottomsize * left->bottomsize +
		right->bottomsize * e->tail->parentEdge->topsize) /
		(e->bottomsize * e->topsize);

	e->distance = wf (lambda, A[left->head->index][right->head->index],
			A[left->head->index][e->tail->index],
			A[left->head->index][sib->head->index],
			A[right->head->index][e->tail->index],
			A[right->head->index][sib->head->index],
			A[sib->head->index][e->tail->index]);

	return;
}

/*********************************************************/

void assignOLSWeights (tree *T, double **A)
{
	edge *e;

	e = depthFirstTraverse (T, NULL);

	while (NULL != e)
	{
		if ((leaf (e->head)) || (leaf(e->tail)))
			OLSext (e, A);

		else
			OLSint (e, A);

		e = depthFirstTraverse (T, e);
	}

	return;
}

/* makes table of average distances from scratch */

/*********************************************************/

void makeOLSAveragesTable (tree *T, double **D, double **A)
{
	edge *e, *f, *g, *h, *exclude;

	e = f = NULL;
	e = depthFirstTraverse (T, e);

	while (NULL != e)
	{
		f = e;
		exclude = e->tail->parentEdge;

		/* We want to calculate A[e->head][f->head] for all edges
		 * except those edges which are ancestral to e.
		 * For those edges, we will calculate A[e->head][f->head]
		 * to have a different meaning, later */

		if (leaf (e->head))
		{
			while (NULL != f)
			{
				if (exclude != f)
				{
					if (leaf (f->head))
						A[e->head->index][f->head->index] =
						A[f->head->index][e->head->index] =
						D[e->head->index2][f->head->index2];

					else
					{
						g = f->head->leftEdge;
						h = f->head->rightEdge;
						A[e->head->index][f->head->index] =
						A[f->head->index][e->head->index] =
							(g->bottomsize * A[e->head->index][g->head->index] +
							h->bottomsize * A[e->head->index][h->head->index]) /
							f->bottomsize;
					}
				}
				else	/* exclude == f */
					exclude = exclude->tail->parentEdge;

				f = depthFirstTraverse (T, f);
			}
		}
		else	/* e->head is not a leaf, so we go recursively
				 * to values calculated for the nodes below e */
		{
			while (NULL !=f )
			{
				if (exclude != f)
				{
					g = e->head->leftEdge;
					h = e->head->rightEdge;
					A[e->head->index][f->head->index] =
					A[f->head->index][e->head->index] =
						(g->bottomsize * A[f->head->index][g->head->index] +
						h->bottomsize * A[f->head->index][h->head->index]) /
						e->bottomsize;
				}
				else
					exclude = exclude->tail->parentEdge;

				f = depthFirstTraverse (T, f);
			}
		}

		/* Now we move to fill up the rest of the table :
	 	* We want A[e->head->index][f->head->index] for those cases
	 	* where e is an ancestor of f, or vice versa.
	 	* We'll do this by choosing e via a depth first-search,
	 	* and the recursing for f up the path to the root */

		f = e->tail->parentEdge;
		if (NULL != f)
			fillTableUp (e, f, A, D, T);

		e = depthFirstTraverse (T, e);
	}

	return;

	/* We are indexing this table by vertex indices,
	 * but really the underlying object is the edge set.
	 * Thus, the array is one element too big in each direction,
	 * but we'll ignore the entries involving the root,
	 * and instead refer to each edge by the head of that edge.
	 * The head of the root points to the edge ancestral to the rest of the tree,
	 * so we'll keep track of up distances by pointing to that head */
}

/*********************************************************/

void GMEcalcDownAverage (node *v, edge *e, double **D, double **A)
{
	edge *left, *right;

	if (leaf(e->head))
		A[e->head->index][v->index] = D[v->index2][e->head->index2];

	else
	{
		left = e->head->leftEdge;
		right = e->head->rightEdge;
		A[e->head->index][v->index] =
			(left->bottomsize * A[left->head->index][v->index] +
			right->bottomsize * A[right->head->index][v->index]) /
			e->bottomsize;
	}

	return;
}

/*********************************************************/

void GMEcalcUpAverage (node *v, edge *e, double **D, double **A)
{
	edge *up, *down;

	if (NULL == e->tail->parentEdge)
		A[v->index][e->head->index] =  D[v->index2][e->tail->index2];

	else
	{
		up = e->tail->parentEdge;
		down = siblingEdge (e);
		A[v->index][e->head->index] =
			(up->topsize * A[v->index][up->head->index] +
			down->bottomsize * A[down->head->index][v->index]) /
			e->topsize;
	}

	return;
}

/* This function calculates average distance D_Xv for each X which is
 * a set of leaves of an induced subtree of T */

/*********************************************************/

void GMEcalcNewvAverages (tree *T, node *v, double **D, double **A)
{
	/* loop over edges
	 * depth-first search */
	edge *e;
	e = NULL;

	e = depthFirstTraverse (T, e);	/* the downward averages need to be
									 * calculated from bottom to top */
	while (NULL != e)
	{
		GMEcalcDownAverage (v, e, D, A);
		e = depthFirstTraverse (T, e);
	}

	e = topFirstTraverse (T, e);	/* the upward averages need to be
									 * calculated from top to bottom */
	while (NULL != e)
	{
		GMEcalcUpAverage (v, e, D, A);
		e = topFirstTraverse (T, e);
	}

	return;
}

/*********************************************************/

double wf4 (double lambda, double lambda2, double D_AB, double D_AC,
	double D_BC, double D_Av, double D_Bv, double D_Cv)
{
	return (((1-lambda) * (D_AC + D_Bv) + (lambda2-1) * (D_AB + D_Cv) +
		(lambda - lambda2) * (D_BC + D_Av)));
}

/* testEdge cacluates what the OLS weight would be if v were inserted
 * into T along e.
 * Compare against known values for inserting along f = e->parentEdge
 * edges are tested by a top-first, left-first scheme.
 * We presume all distances are fixed to the correct weight for
 * e->parentEdge, if e is a left-oriented edge */

/*********************************************************/

void testEdge (edge *e, node *v, double **A)
{
  double lambda, lambda2;
  edge *par, *sib;

  sib = siblingEdge (e);
  par = e->tail->parentEdge;

  /* C is set above e->tail
   * B is set below e
   * A is set below sib
   * following the nomenclature of Desper & Gascuel */

   lambda =  (((double) (sib->bottomsize + e->bottomsize*par->topsize)) /
   			((1 + par->topsize) * (par->bottomsize)));
	lambda2 = (((double) (sib->bottomsize + e->bottomsize*par->topsize)) /
			((1 + e->bottomsize) * (e->topsize)));
	e->totalweight = par->totalweight +
		wf4 (lambda, lambda2, A[e->head->index][sib->head->index],
			A[sib->head->index][e->tail->index],
			A[e->head->index][e->tail->index],
			A[sib->head->index][v->index],A[e->head->index][v->index],
			A[v->index][e->tail->index]);

	return;
}

/* The key function of the program addSpecies inserts the node v to the
 * tree T. It uses testEdge to see what the weight would be if v split
 * a particular edge. Weights are assigned by OLS formula */

/*********************************************************/

tree *GMEaddSpecies (tree *T, node *v, double **D, double **A)
{
	tree *T_e;
	edge *e;				// loop variable
	edge *e_min;			//	points to best edge seen thus far
	double w_min = 0.0;		// used to keep track of tree weights

	if (verbose > 2 && !isBoostrap)
		Debug ( (char*)"Adding %s.", v->label);

	/* initialize variables as necessary */

	/* CASE 1: T is empty, v is the first node */
	if (NULL == T)	//create a tree with v as only vertex, no edges
	{
		T_e = newTree();
		T_e->root = v;
		/* note that we are rooting T arbitrarily at a leaf.
		 * T->root is not the phylogenetic root */
		v->index = 0;
		T_e->size = 1;
		return (T_e);
	}

	/* CASE 2: T is a single-vertex tree */
	if (1 == T->size)
	{
		v->index = 1;
		e = makeEdge ( (char*)"", T->root, v, 0.0);
		snprintf (e->label,2,"E1");
		e->topsize = 1;
		e->bottomsize = 1;
		A[v->index][v->index] = D[v->index2][T->root->index2];
		T->root->leftEdge = v->parentEdge = e;
		T->size = 2;
		return (T);
	}

	/* CASE 3: T has at least two nodes and an edge.
	 * Insert new node by breaking one of the edges */
	v->index = T->size;
	GMEcalcNewvAverages (T, v, D, A);
	/* calcNewvAverges will assign values to all the edge averages of T
	 * which include the node v. Will do so using pre-existing averages
	 * in T and information from A,D*/

	e_min = T->root->leftEdge;
	e = e_min->head->leftEdge;
	while (NULL != e)
	{
		testEdge (e, v, A);
		/* testEdge tests weight of tree if loop variable e is the edge
		 * split, places this weight in e->totalweight field */
		if (e->totalweight < w_min)
		{
			e_min = e;
			w_min = e->totalweight;
		}
		e = topFirstTraverse (T, e);
	}
	// e_min now points at the edge we want to split

	GMEsplitEdge (T, v, e_min, A);

	return (T);
}

/* The ME version of updateAveragesMatrix does not update the entire
 * matrix A, but updates A[v->index][w->index] whenever this represents
 * an average of 1-distant or 2-distant subtrees */

/*********************************************************/

void GMEupdateAveragesMatrix (double **A, edge *e, node *v, node *newNode)
{
	edge *sib, *par, *left, *right;

	sib = siblingEdge (e);
	left = e->head->leftEdge;
	right = e->head->rightEdge;
	par = e->tail->parentEdge;

	/* We need to update the matrix A so all 1-distant, 2-distant,
	 * and 3-distant averages are correct */

	/* First, initialize the newNode entries */
	/* 1-distant */
	A[newNode->index][newNode->index] =
		(e->bottomsize * A[e->head->index][e->head->index] +
		A[v->index][e->head->index]) / (e->bottomsize + 1);

	/* 1-distant for v */
	A[v->index][v->index] =
		(e->bottomsize * A[e->head->index][v->index] +
		e->topsize * A[v->index][e->head->index]) /
		(e->bottomsize + e->topsize);

	/* 2-distant for v,newNode */
	A[v->index][newNode->index] =
	A[newNode->index][v->index] =
	A[v->index][e->head->index];

	/* second 2-distant for newNode */
	A[newNode->index][e->tail->index] =
	A[e->tail->index][newNode->index] =
		(e->bottomsize * A[e->head->index][e->tail->index]+
		A[v->index][e->tail->index]) / (e->bottomsize + 1);

	/* third 2-distant for newNode */
	A[newNode->index][e->head->index] =
	A[e->head->index][newNode->index] =
	A[e->head->index][e->head->index];

	if (NULL != sib)
	{
		/* fourth and last 2-distant for newNode */
		A[newNode->index][sib->head->index] =
		A[sib->head->index][newNode->index] =
			(e->bottomsize * A[sib->head->index][e->head->index] +
			A[sib->head->index][v->index]) / (e->bottomsize + 1);
		updateSubTreeAverages (A, sib, v, SKEW);	/* updates sib and below */
	}

	if (NULL != par)
	{
		if (e->tail->leftEdge == e)
			updateSubTreeAverages (A, par, v, LEFT);	/* updates par and above */

		else
			updateSubTreeAverages (A, par, v, RIGHT);
	}

	if (NULL != left)
		updateSubTreeAverages (A, left, v, UP);		/* updates left and below */

	if (NULL != right)
		updateSubTreeAverages (A, right, v, UP);	/* updates right and below */

	/* 1-dist for e->head */
	A[e->head->index][e->head->index] =
		(e->topsize * A[e->head->index][e->head->index] +
		A[e->head->index][v->index]) / (e->topsize+1);

	/* 2-dist for e->head (v, newNode, left, right)
	 * taken care of elsewhere
	 * 3-dist with e->head either taken care of elsewhere (below)
	 * or unchanged (sib, e->tail) */

	/* Symmetrize the matrix (at least for distant-2 subtrees) */
	A[v->index][e->head->index] = A[e->head->index][v->index];

	/* and distant-3 subtrees */
	A[e->tail->index][v->index] = A[v->index][e->tail->index];

	if (NULL != left)
		A[v->index][left->head->index] = A[left->head->index][v->index];

	if (NULL != right)
		A[v->index][right->head->index] = A[right->head->index][v->index];

	if (NULL != sib)
		A[v->index][sib->head->index] = A[sib->head->index][v->index];

	return;
}

/*********************************************************/

void GMEsplitEdge (tree *T, node *v, edge *e, double **A)
{
	char nodelabel[MAX_NAME_LENGTH];
	char edgelabel[MAX_NAME_LENGTH];
	edge *newPendantEdge;
	edge *newInternalEdge;
	node *newNode;

	snprintf(nodelabel,1," ");
	newNode = makeNode (nodelabel, T->size + 1);

	snprintf (edgelabel, MAX_NAME_LENGTH, "E%d", T->size);
	newPendantEdge = makeEdge (edgelabel, newNode, v, 0.0);

	snprintf (edgelabel, MAX_NAME_LENGTH, "E%d", T->size+1);
	newInternalEdge = makeEdge (edgelabel, newNode, e->head, 0.0);

	if (verbose > 2 && !isBoostrap)
	{
		if ((NULL == e->tail->label) || (strlen (e->tail->label) == 0))
		{
			if ((NULL == e->head->label) || (strlen (e->head->label) == 0))
				Debug ( (char*)"Inserting node '%s' on edge '%s' between internal nodes.",
					v->label, e->label);
			else
				Debug ( (char*)"Inserting node '%s' on edge '%s' between node '%s' and an internal node.",
					v->label, e->label, e->head->label);
		}
		else
		{
			if ((NULL == e->head->label) || (strlen (e->head->label) == 0))
				Debug ( (char*)"Inserting node '%s' on edge '%s' between node '%s' and an internal node.",
					v->label, e->label, e->tail->label);
			else
				Debug ( (char*)"Inserting node '%s' on edge '%s' between nodes '%s' and '%s'",
					v->label, e->label, e->tail->label, e->head->label);
		}
	}

	/* update the matrix of average distances
	 * also updates the bottomsize, topsize fields*/

	GMEupdateAveragesMatrix (A, e, v, newNode);

	newNode->parentEdge = e;
	e->head->parentEdge = newInternalEdge;
	v->parentEdge = newPendantEdge;
	e->head = newNode;

	T->size = T->size + 2;

	if (e->tail->leftEdge == e)
	{
		newNode->leftEdge = newInternalEdge;
		newNode->rightEdge = newPendantEdge;
	}
	else
	{
		newNode->leftEdge = newInternalEdge;
		newNode->rightEdge = newPendantEdge;
	}

	// assign proper topsize, bottomsize values to the two new Edges
	newPendantEdge->bottomsize = 1;
	newPendantEdge->topsize = e->bottomsize + e->topsize;

	newInternalEdge->bottomsize = e->bottomsize;
	newInternalEdge->topsize = e->topsize;	// off by one, but we adjust that below

	// and increment these fields for all other edges
	updateSizes (newInternalEdge, UP);
	updateSizes (e, DOWN);

	return;
}

/* The monster function of this program */

/*********************************************************/

void updateSubTreeAverages (double **A, edge *e, node *v, int direction)
{
	edge *sib, *left, *right, *par;

	left = e->head->leftEdge;
	right = e->head->rightEdge;
	sib = siblingEdge (e);
	par = e->tail->parentEdge;

	switch (direction)
	{
		/* Want to preserve correctness of all
		 * 1-distant, 2-distant, and 3-distant averages
		 *
		 * 1-distant updated at edge splitting the two trees
		 * 2-distant updated :
		 * 	(left->head,right->head) and (head,tail) updated at	a given edge.
		 *  Note, NOT updating (head,sib->head)!
		 * 	(That would lead to multiple updating)
		 * 3-distant updated at edge in center of quartet */

		case UP :	/* point of insertion is above e */
			/* 1-distant average of nodes below e to nodes above e */
			A[e->head->index][e->head->index] =
				(e->topsize * A[e->head->index][e->head->index] +
				A[e->head->index][v->index]) / (e->topsize + 1);
			/* 2-distant average of nodes below e to nodes above parent of e */
			A[e->head->index][par->head->index] =
			A[par->head->index][e->head->index] =
				(par->topsize * A[par->head->index][e->head->index] +
				A[e->head->index][v->index]) / (par->topsize + 1);

			/* must do both 3-distant averages involving par */
			if (NULL != left)
			{
				updateSubTreeAverages (A, left, v, UP);	/* and recursive call */
				/* 3-distant average */
				A[par->head->index][left->head->index] =
				A[left->head->index][par->head->index] =
					(par->topsize * A[par->head->index][left->head->index] +
					A[left->head->index][v->index]) / (par->topsize + 1);
			}

			if (NULL != right)
			{
				updateSubTreeAverages (A, right, v, UP);
				A[par->head->index][right->head->index] =
				A[right->head->index][par->head->index] =
					(par->topsize * A[par->head->index][right->head->index] +
					A[right->head->index][v->index]) / (par->topsize + 1);
			}
			break;

		case SKEW :	/* point of insertion is skew to e */
			/* 1-distant average of nodes below e to nodes above e */
			A[e->head->index][e->head->index] =
				(e->topsize * A[e->head->index][e->head->index] +
				A[e->head->index][v->index]) / (e->topsize + 1);

			/* no 2-distant averages to update in this case
			 * updating 3-distant averages involving sib */
			if (NULL != left)
			{
				updateSubTreeAverages (A, left, v, UP);
				A[sib->head->index][left->head->index] =
				A[left->head->index][sib->head->index] =
					(sib->bottomsize * A[sib->head->index][left->head->index] +
					A[left->head->index][v->index]) / (sib->bottomsize + 1);
			}

			if (NULL != right)
			{
				updateSubTreeAverages (A, right, v, UP);
				A[sib->head->index][right->head->index] =
				A[right->head->index][sib->head->index] =
					(sib->bottomsize * A[par->head->index][right->head->index] +
					A[right->head->index][v->index]) / (sib->bottomsize + 1);
			}
			break;

		case LEFT :	/* point of insertion is below the edge left */
			/* 1-distant average */
			A[e->head->index][e->head->index] =
				(e->bottomsize * A[e->head->index][e->head->index] +
				A[v->index][e->head->index]) / (e->bottomsize + 1);
			/* 2-distant averages */
			A[e->head->index][e->tail->index] =
			A[e->tail->index][e->head->index] =
				(e->bottomsize * A[e->head->index][e->tail->index] +
				A[v->index][e->tail->index]) / (e->bottomsize + 1);
			A[left->head->index][right->head->index] =
			A[right->head->index][left->head->index] =
				(left->bottomsize * A[right->head->index][left->head->index] +
				A[right->head->index][v->index]) / (left->bottomsize + 1);
			/* 3-distant avereages involving left */
			if (NULL != sib)
			{
				updateSubTreeAverages (A, sib, v, SKEW);
				A[left->head->index][sib->head->index] =
				A[sib->head->index][left->head->index] =
					(left->bottomsize * A[left->head->index][sib->head->index] +
					A[sib->head->index][v->index]) / (left->bottomsize + 1);
			}

			if (NULL != par)
			{
				 if (e->tail->leftEdge == e)
				 	updateSubTreeAverages (A, par, v, LEFT);

				else
					updateSubTreeAverages (A, par, v, RIGHT);

				A[left->head->index][par->head->index] =
				A[par->head->index][left->head->index] =
					(left->bottomsize * A[left->head->index][par->head->index] +
					A[v->index][par->head->index]) / (left->bottomsize + 1);
			}
			break;

		case RIGHT :	/* point of insertion is below the edge right */
			/* 1-distant average */
			A[e->head->index][e->head->index] =
				(e->bottomsize * A[e->head->index][e->head->index] +
				A[v->index][e->head->index]) / (e->bottomsize + 1);
			/* 2-distant averages */
			A[e->head->index][e->tail->index] =
			A[e->tail->index][e->head->index] =
				(e->bottomsize * A[e->head->index][e->tail->index] +
				A[v->index][e->tail->index]) / (e->bottomsize + 1);
			A[left->head->index][right->head->index] =
			A[right->head->index][left->head->index] =
				(right->bottomsize * A[right->head->index][left->head->index] +
				A[left->head->index][v->index]) / (right->bottomsize + 1);
			/* 3-distant avereages involving right */
			if (NULL != sib)
			{
				updateSubTreeAverages (A, sib, v, SKEW);
				A[right->head->index][sib->head->index] =
				A[sib->head->index][right->head->index] =
					(right->bottomsize * A[right->head->index][sib->head->index] +
					A[sib->head->index][v->index]) / (right->bottomsize + 1);
			}

			if (NULL != par)
			{
				if (e->tail->leftEdge == e)
					updateSubTreeAverages (A, par, v, LEFT);

				else
					updateSubTreeAverages (A, par, v, RIGHT);

				A[right->head->index][par->head->index] =
				A[par->head->index][right->head->index] =
					(right->bottomsize * A[right->head->index][par->head->index] +
					A[v->index][par->head->index]) / (right->bottomsize + 1);
			}
			break;
	}

	return;
}

/*********************************************************/

void assignBottomsize (edge *e)
{
	if (leaf(e->head))
		e->bottomsize = 1;

	else
	{
		assignBottomsize (e->head->leftEdge);
		assignBottomsize (e->head->rightEdge);
		e->bottomsize = e->head->leftEdge->bottomsize + e->head->rightEdge->bottomsize;
	}

	return;
}

/*********************************************************/

void assignTopsize (edge *e, int numLeaves)
{
	if (NULL != e)
	{
		e->topsize = numLeaves - e->bottomsize;
		assignTopsize (e->head->leftEdge,numLeaves);
		assignTopsize (e->head->rightEdge,numLeaves);
	}

	return;
}

/*********************************************************/

void assignAllSizeFields (tree *T)
{
	assignBottomsize (T->root->leftEdge);
	assignTopsize (T->root->leftEdge, T->size / 2 + 1);

	return;
}

