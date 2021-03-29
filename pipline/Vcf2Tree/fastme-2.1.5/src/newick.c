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


#include "newick.h"

//int nodeCount;
//int edgeCount;

/* decodeNewickSubtree is used to turn a string of the form
 * "(v1:d1,v2:d2,(subtree) v3:d3....vk:dk) subroot:d," into a subtree
 * rooted at subrooted, with corresponding subtrees and leaves at v1
 * through vk. It is called recursively on subtrees */

/*********************************************************/

node *decodeNewickSubtree (char *treeString, tree *T, int *uCount,
	int *nodeCount, int *edgeCount)
{
	node *thisNode = NULL;
	node *centerNode;
	double thisWeight = 0;
	edge *thisEdge;
	char stringWeight[MAX_NAME_LENGTH];
	int state;
	int i = 0;
	int j;
	int left,right;
	int parcount;

	left = right = 0;
	parcount = 0;
	state = ReadOpenParenthesis;

	if ('(' == treeString[0])
		parcount++;
	
	centerNode = makeNode ( (char*)"", *nodeCount);
	*nodeCount = *nodeCount +1;
	
	T->size++;

	while (parcount > 0)
	{
		while (whiteSpace (treeString[i]))
			i++;

		switch (state)
		{
			case (ReadOpenParenthesis) :
				if ('(' != treeString[0])
					Exit ( (char*)"Invalid Newick string.");

				i++;
				state = ReadSubTree;
				break;

			case (ReadSubTree) :
				if ('(' == treeString[i])	/* if treeString[i] is a left parenthesis,
											 * we scan down the string until we find its partner.
											 * The relative number of '('s vs. ')'s is counted
											 * by the variable parcount */
				{
					left = i++;
					parcount++;
					while (parcount > 1)
					{
						while (('(' != treeString[i]) && (')' != treeString[i]))
							i++;	/* skip over characters which are not parentheses */

						if ('(' == treeString[i])
							parcount++;

						else if (')' == treeString[i])
							parcount--;

						i++;
					}	/* end while */

					right = i;		/* at this point, the subtree string goes from
									 *  treeString[left] to treeString[right - 1] */

					thisNode = decodeNewickSubtree (treeString + left, T, uCount, nodeCount, edgeCount);
							  /*note that this step will put thisNode in T*/

					i = right;	/* having created the node for the subtree, we move
								 * to assigning the label for the new node.
								 * treeString[right] will be the start of this label */
				}	/* end if ('(' == treeString[i]) condition */
				else
				{
					thisNode = makeNode ( (char*)" ", *nodeCount);
					*nodeCount = *nodeCount +1;
					
					T->size++;
				}

				state = ReadLabel;
				break;

			case (ReadLabel) :
				left = i;					/* recall "left" is the left marker for the substring, "right" the right */
				if (':' == treeString[i])	/* this means an internal node? */
				{
					snprintf (thisNode->label, 1, " ");
					right = i;
				}
				else
				{
					while ((':' != treeString[i]) && (',' != treeString[i]) && (')' != treeString[i]))
						i++;

					right = i;
					j = 0;

					for (i = left; i < right; i++)
						if (! (whiteSpace (treeString[i])))
							thisNode->label[j++] = treeString[i];

					thisNode->label[j] = '\0';
				}

				if (':' == treeString[right])
					state = ReadWeight;
				else
				{
					state = AddEdge;
					thisWeight = 0.0;
				}
				i = right + 1;
				break;

			case (ReadWeight) :
				left = i;
				while ((',' != treeString[i]) && (')' != treeString[i]))
					i++;

				right = i;
				j = 0;
				for (i = left; i < right; i++)
					stringWeight[j++] = treeString[i];

				stringWeight[j] = '\0';
				thisWeight = atof (stringWeight);
				state = AddEdge;
				break;

			case (AddEdge) :
				thisEdge = makeEdge ( (char*)" ", centerNode, thisNode, thisWeight);
				thisNode->parentEdge = thisEdge;
				if (NULL == centerNode->leftEdge)
					centerNode->leftEdge = thisEdge;

				else if (NULL == centerNode->rightEdge)
					centerNode->rightEdge = thisEdge;

				else if (NULL == centerNode->middleEdge)
					centerNode->middleEdge = thisEdge;

				else
					Exit ( (char*)"Node %s has too many (>3) children.\n", centerNode->label);

				snprintf (thisEdge->label, MAX_NAME_LENGTH, "E%d", *edgeCount);
				*edgeCount = *edgeCount +1;
				i = right + 1;
				if (',' == treeString[right])
					state = ReadSubTree;
				else
					parcount--;

				break;
		}
	}

	return (centerNode);
}

/*********************************************************/

tree *readNewickString (char *str)
{
	tree *T;
	node *centerNode;
	int i = 0;
	int j = 0;
	int inputLength, uCount, parCount, nodeCount, edgeCount;
	char rootLabel[MAX_NAME_LENGTH];

	uCount = parCount = nodeCount = edgeCount = 0;

	T = newTree();

	if ('(' != str[0])
		Exit ( (char*)"Generated tree does not start with '('.");

	inputLength = (int) strlen (str) +1;
	for (i = 0; i < inputLength; i++)
	{
		if ('(' == str[i])
			parCount++;

		else if (')' == str[i])
			parCount--;

		if (0 == parCount)
		{
			i++;
			while((';' != str[i]) && (! (whiteSpace (str[i]))) && (j < MAX_NAME_LENGTH))
				rootLabel[j++] = str[i++];

			rootLabel[j] = '\0';
			i = inputLength;
		}
		else if (parCount < 0)
			Exit ( (char*)"Generated tree has too many right parentheses.");
	}
	centerNode = decodeNewickSubtree (str, T, &uCount, &nodeCount, &edgeCount);
	snprintf (centerNode->label, MAX_NAME_LENGTH, "%s", rootLabel);
	T->root = centerNode;

	return (T);
}

/*********************************************************/

tree *loadNewickTree (FILE *ifile, int numLeaves)
{
	char c;
	char rootLabel[MAX_NAME_LENGTH];
	char *nextString;
	boolean Comment;
	int i, j, uCount, parCount, inputLength, nodeCount, edgeCount;
	tree *T;
	node *centerNode;

	i = j = uCount = parCount = nodeCount = edgeCount = 0;
	
	T = newTree();
	
	nextString = (char *) mCalloc (numLeaves * INPUT_SIZE, sizeof (char));
	if (NULL == nextString)
		nextString = (char *) mCalloc (MAX_INPUT_SIZE, sizeof(char));

	Comment = FALSE;
	while (1 == fscanf(ifile,"%c",&c))
	{
		if('[' == c)
			Comment = TRUE;

		else if (']' == c)
			Comment = FALSE;

		else if (!(Comment))
		{
			if (whiteSpace(c))
			{
				if (i > 0)
					nextString[i++] = ' ';
			}
			else
				nextString[i++] = c;

			if (';' == c)
				break;
		}
	}

	if ('(' != nextString[0])
		Exit ( (char*)"Invalid input tree file format. Does not start with '('.");

	// Unroot tree in newick string if any
	if (isNwkRootedTree (nextString))
	{
		Message ( (char*)"Input tree is rooted. Unrooting...");
		Unode *tree = UreadNewick (nextString, (int) strlen (nextString) -1);
		UprintTree (nextString, tree);
		UfreeNodes (tree);
	}

	inputLength = (int) strlen (nextString);
	
	for (i = 0; i < inputLength; i++)
	{
		if ('(' == nextString[i])
			parCount++;

		else if (')' == nextString[i])
			parCount--;


		if (0 == parCount)
		{
			i++;
			while((';' != nextString[i]) && (! (whiteSpace (nextString[i]))) && (j < MAX_NAME_LENGTH))
				rootLabel[j++] = nextString[i++];

			rootLabel[j] = '\0';
			i = inputLength;
		}
		else if (parCount < 0)
			Exit ( (char*)"Invalid input tree file format. Too many right parentheses.");
	}

	centerNode = decodeNewickSubtree (nextString, T, &uCount, &nodeCount, &edgeCount);
	snprintf (centerNode->label, MAX_NAME_LENGTH, "%s", rootLabel);
	T->root = centerNode;
	free (nextString);
	if (NULL != T->root->parentEdge)
		Exit ( (char*)"Tree poorly rooted.");

	return(T);
}

/*********************************************************/

void NewickPrintSubtree(tree *T, edge *e, FILE *ofile, const char *format)
{
	if (NULL == e)
		Exit ( (char*)"Newick Printing routine error.");

	if(!(leaf(e->head)))
	{
		fprintf(ofile,"(");
		NewickPrintSubtree(T,e->head->leftEdge,ofile,format);
		fprintf(ofile,",");
		NewickPrintSubtree(T,e->head->rightEdge,ofile,format);
		fprintf(ofile,")");
	}
	fprintf(ofile,"%s",e->head->label);
	//fprintf(ofile,":%f",e->distance);
	fprintf(ofile,format,e->distance);
	
	return;
}

/*********************************************************/

void NewickPrintBinaryTree(tree *T, FILE *ofile, const char *format)
{
	edge *e, *f;
	node *rootchild;
	e = T->root->leftEdge;
	rootchild = e->head;
	fprintf(ofile,"(");
	f = rootchild->leftEdge;
	if (NULL != f)
	{
		NewickPrintSubtree(T,f,ofile,format);
		fprintf(ofile,",");
	}
	f = rootchild->rightEdge;
	if (NULL != f)
	{
		NewickPrintSubtree(T,f,ofile,format);
		fprintf(ofile,",");
	}
	//fprintf(ofile,"%s:%f",T->root->label,e->distance);
	fprintf(ofile,"%s",T->root->label);
	fprintf(ofile,format,e->distance);
	fprintf(ofile,")");
	if (NULL != rootchild->label)
		fprintf(ofile,"%s",rootchild->label);
	fprintf(ofile,";\n");
	
	return;
}

/*********************************************************/

void NewickPrintTrinaryTree(tree *T, FILE *ofile, const char *format)
{
	edge *f;
	f = T->root->leftEdge;
	fprintf(ofile,"(");
	if (NULL != f)
	{
		NewickPrintSubtree(T,f,ofile,format);
		fprintf(ofile,",");
	}
	f = T->root->rightEdge;
	if (NULL != f)
	{
		NewickPrintSubtree(T,f,ofile,format);
		fprintf(ofile,",");
	}
	f = T->root->middleEdge;
	if (NULL != f)
	{
		NewickPrintSubtree(T,f,ofile,format);
		fprintf(ofile,")");
	}
	if (NULL != T->root->label)
		fprintf(ofile,"%s",T->root->label);
	fprintf(ofile,";\n");
	
	return;
}

/*********************************************************/

void NewickPrintTree(tree *T, FILE *ofile, int precision)
{
	char format[8];
	snprintf (format, 8, ":%%.%df", precision);
	
  if (leaf(T->root))
    NewickPrintBinaryTree(T,ofile,format);
  else
    NewickPrintTrinaryTree(T,ofile,format);
    
  return;
}

/*********************************************************/

void NewickPrintSubtreeStr(tree *T, edge *e, char *str, const char *format)
{
	char *tmp;
	if (NULL == e)
		Exit ( (char*)"Newick Printing routine error.");

	if(!(leaf(e->head)))
	{
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, "(", 1);
		NewickPrintSubtreeStr(T,e->head->leftEdge,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ",", 1);
		NewickPrintSubtreeStr(T,e->head->rightEdge,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ")", 1);
	}
	if (strlen (str) < MAX_INPUT_SIZE - strlen (e->head->label) -1)
		strncat (str, e->head->label, strlen (e->head->label));

	//if (strlen (str) < MAX_INPUT_SIZE - 2)
	//	strncat (str, ":", 1);

	tmp = (char *) mCalloc (INPUT_SIZE, sizeof(char));
	if (strlen(tmp))
		strncpy(tmp, "", strlen(tmp));
	//snprintf (tmp, INPUT_SIZE, "%f", e->distance);
	snprintf (tmp, INPUT_SIZE, format, e->distance);
	if (strlen (str) < MAX_INPUT_SIZE - strlen (tmp) -1)
		strncat (str, tmp, strlen (tmp));
	free (tmp);
	
  return;
}

/*********************************************************/

void NewickPrintBinaryTreeStr(tree *T, char *str, const char *format)
{
	edge *e, *f;
	node *rootchild;
	char *tmp;
	e = T->root->leftEdge;
	rootchild = e->head;

	if (strlen (str) < MAX_INPUT_SIZE -2)
		strncat (str, "(", 1);
	f = rootchild->leftEdge;
	if (NULL != f)
	{
		NewickPrintSubtreeStr(T,f,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ",", 1);
	}
	f = rootchild->rightEdge;
	if (NULL != f)
	{
		NewickPrintSubtreeStr(T,f,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ",", 1);
	}
	if (strlen (str) < MAX_INPUT_SIZE - strlen (T->root->label) -1)
		strncat (str, T->root->label, strlen (T->root->label));

	//if (strlen (str) < MAX_INPUT_SIZE - 2)
	//	strncat (str, ":", 1);

	tmp = (char *) mCalloc (INPUT_SIZE, sizeof(char));
	if (strlen(tmp))
		strncpy(tmp, "", strlen(tmp));
	//snprintf (tmp, INPUT_SIZE, "%f", e->distance);
	snprintf (tmp, INPUT_SIZE, format, e->distance);
	if (strlen (str) < MAX_INPUT_SIZE - strlen (tmp) -1)
		strncat (str, tmp, strlen (tmp));
	free (tmp);

	if (strlen (str) < MAX_INPUT_SIZE - 2)
		strncat (str, ")", 1);

	if (NULL != rootchild->label)
		if (strlen (str) < MAX_INPUT_SIZE - strlen (rootchild->label) -1)
		strncat (str, T->root->label, strlen (rootchild->label));

	if (strlen (str) < MAX_INPUT_SIZE - 3)
		strncat (str, ";\n", 2);
	
	return;
}

/*********************************************************/

void NewickPrintTrinaryTreeStr(tree *T, char *str, const char *format)
{
	edge *f;
	f = T->root->leftEdge;
	if (strlen (str) < MAX_INPUT_SIZE -2)
		strncat (str, "(", 1);
	
	if (NULL != f)
    {
		NewickPrintSubtreeStr(T,f,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ",", 1);
	}
	f = T->root->rightEdge;
	if (NULL != f)
	{
		NewickPrintSubtreeStr(T,f,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ",", 1);
	}
	f = T->root->middleEdge;
	if (NULL != f)
	{
		NewickPrintSubtreeStr(T,f,str,format);
		if (strlen (str) < MAX_INPUT_SIZE -2)
			strncat (str, ")", 1);
	}
	if (NULL != T->root->label)
		if (strlen (str) < MAX_INPUT_SIZE - strlen (T->root->label) -1)
			strncat (str, T->root->label, strlen (T->root->label));
	
	if (strlen (str) < MAX_INPUT_SIZE - 3)
		strncat (str, ";\n", 2);
	
	return;
}

/*********************************************************/

void NewickPrintTreeStr(tree *T, char *str, int precision)
{
	char format[8];
	snprintf (format, 8, ":%%.%df", precision);
	
	if (leaf(T->root))
	{
		NewickPrintBinaryTreeStr(T,str,format);
	}
	else
	{
		NewickPrintTrinaryTreeStr(T,str,format);
	}
	return;
}

/*********************************************************/

boolean isNwkRootedTree (char *str)
{
	int i, inputLength, parCount, virCount;
	
	i = parCount = virCount = 0;

	inputLength = (int) strlen (str)+1;

	for (i = 0; i < inputLength; i++)
	{
		if ('(' == str[i])
			parCount++;

		else if (')' == str[i])
			parCount--;

		if ((1 == parCount) && (',' == str[i]))
		{
			virCount++;
		}
	}
	
	if (virCount > 1)
		return (FALSE);
	else
		return (TRUE);
}

/*********************************************************/

Unode *UnewNode (void)
{
	Unode *r = (Unode *) mCalloc (1, sizeof (struct Unode));
	r->father = NULL;
	r->child = NULL;
	r->length = 0;
	r->value = 0.0;
	r->name = NULL;
	
	return r;
}

/*********************************************************/

void UfreeNodes (Unode *tree)
{
	if (NULL==tree)
		return;
	
	for (int i=0 ; i<tree->length ; i++)
		UfreeNodes (tree->child[i]);
	
	free (tree->name);
	free (tree->child);
	free (tree);
	
	return;
}

/*********************************************************/

Unode *UreadNewick (char* line, int length)
{
	if( line==NULL || length==0 )
		return NULL;
	
	Unode *root = UnewNode();
	Unode **tmpNODE = NULL;
	int dblePt = length-1;
	while (dblePt>=0 && line[dblePt]!=':' && line[dblePt]!=')')
		dblePt--;

	int parD = dblePt;
	while (parD>=0 && line[parD]!=')')
		parD--;
	
	if (dblePt > parD)
	{
		root->value = atof (line+dblePt+1);
	}
	else
		dblePt = length;
	
	dblePt = dblePt-parD-1;
	
	if (dblePt>0)
	{
		root->name = (char *) mCalloc (dblePt+1, sizeof(char));
		strncpy (root->name, line+parD+1, (unsigned long) dblePt);
		root->name[dblePt]='\0';
	}
	
	int par=0;
	int start=1;
	
	for (int i=1 ; i<=parD ; i++)
	{
		if (line[i]=='(')
			par++;
		else if (line[i]==')')
			par--;
		
		if ( ( line[i]==',' && par==0 ) || ( i==parD ) )
		{
			tmpNODE = (Unode **) mCalloc (root->length+1, sizeof (struct node*));
			
			for (int j=0 ; j<root->length ; j++)
				tmpNODE[j] = root->child[j];
			
			free (root->child);
			root->child = tmpNODE;
			root->child[root->length] = UreadNewick (line+start, i-start);
			root->child[root->length]->father = root;
			root->length++;
			start = i+1;
		}
	}
	
	return root;
}

/*********************************************************/

void UprintTree (char *f, Unode *tree)
{
	int i, j;
	i = 0;
	j = 1 ;
	
	if (tree->child[0]->length!=2)
	{
		i = 1 ;
		j = 0 ;
	}
	
	sprintf (f, "(" );
	UprintsubTree (f, tree->child[i]->child[0]);
	strncat (f, ",", 1);
	UprintsubTree (f, tree->child[i]->child[1]);
	strncat (f, ",", 1);
	tree->child[j]->father = NULL;
	UprintsubTree (f, tree->child[j]);
	tree->child[j]->father = tree;
	
	char *tmp = (char *) mCalloc (20, sizeof(char));
	snprintf (tmp, 20, ":%f);", tree->child[i]->value + tree->child[j]->value);
	strncat (f, tmp, strlen(tmp));
	free (tmp);

	return;
}

/*********************************************************/

void UprintsubTree (char *f, Unode *tree)
{
	if (tree==NULL )
		return;
	
	if (tree->length!=0)
	{
		strncat (f, "(", 1);
		UprintsubTree (f, tree->child[0]);
		for (int i=1; i<tree->length ; i++)
		{
			strncat( f, ",", 1);
			UprintsubTree (f, tree->child[i]);
		}
		strncat (f, ")", 1);
	}
	
	if (tree->name!=NULL)
		strncat (f, tree->name, strlen (tree->name));
	
	if (tree->father!=NULL)
	{
		char *tmp = (char *) mCalloc (20, sizeof(char));
		snprintf (tmp, 20, ":%f", tree->value);
		strncat (f, tmp, strlen(tmp));
		free (tmp);
	}
	
	return;
}

