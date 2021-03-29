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


#ifndef GRAPH_H_
#define GRAPH_H_

#include "utils.h"

typedef struct node
{
	char label[MAX_NAME_LENGTH];
	struct edge *parentEdge;
	struct edge *leftEdge;
	struct edge *middleEdge;
	struct edge *rightEdge;
	int index;
	int index2;
} node;

typedef struct edge
{
	char label[MAX_NAME_LENGTH];
	struct node *tail;		/* for edge (u,v), u is the tail, v is the head */
	struct node *head;
	int bottomsize;			/* number of nodes below edge */
	int topsize;			/* number of nodes above edge */
	double distance;
	double totalweight;
} edge;

typedef struct tree
{
	struct node *root;
	int size;
	double weight;
} tree;

typedef struct set
{
	struct node *firstNode;
	struct set *secondNode;
} set;

/* from traverse.c */
edge *depthFirstTraverse (tree *T, edge *e);


boolean leaf (node *v);
set *addToSet (node *v, set *X);
set *copySet (set *X);
node *makeNode (char *label, int index);
edge *makeEdge (char *label, node *tail, node *head, double weight);
tree *newTree (void);
void freeSubTree (edge *e);
void freeTree (tree *T);
void freeSet (set *S);
void freeNode (node *n);
node *copyNode (node *v);
edge *copyEdge (edge *e);
tree *detrifurcate (tree *T);
edge *siblingEdge (edge *e);
void updateSizes (edge *e, int direction);
node *copySubtree (node *v);
tree *copyTree (tree *T);
void weighTree (tree *T);
edge *findEdge (tree *T, edge *e);
node *indexedNode (tree *T, int i);
edge *indexedEdge (tree *T, int i);
boolean checkLabelExist (set *S, char *label);

#endif /*GRAPH_H_*/

