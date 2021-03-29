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


#ifndef NEWICK_H_
#define NEWICK_H_

#include "graph.h"


typedef struct Unode
{
	double  value;
	char*   name;
	int     length;
	struct Unode** child;
	struct Unode* father;
} Unode;


node *decodeNewickSubtree (char *treeString, tree *T, int *uCount,
	int *nodeCount, int *edgeCount);
	
tree *readNewickString (char *str);

tree *loadNewickTree (FILE *ifile, int numLeaves);

void NewickPrintSubtree (tree *T, edge *e, FILE *ofile, const char *format);

void NewickPrintBinaryTree (tree *T, FILE *ofile, const char *format);

void NewickPrintTrinaryTree (tree *T, FILE *ofile, const char *format);

void NewickPrintTree (tree *T, FILE *ofile, int precision);

void NewickPrintSubtreeStr (tree *T, edge *e, char *str, const char *format);

void NewickPrintBinaryTreeStr (tree *T, char *str, const char *format);

void NewickPrintTrinaryTreeStr (tree *T, char *str, const char *format);

void NewickPrintTreeStr (tree *T, char *str, int precision);

boolean isNwkRootedTree (char *str);

Unode *UnewNode (void);

void UfreeNodes (Unode* tree);

Unode *UreadNewick (char* line, int length);

void UprintTree (char *f, Unode *tree);

void UprintsubTree (char *f, Unode *tree);

#endif /*NEWICK_H_*/
