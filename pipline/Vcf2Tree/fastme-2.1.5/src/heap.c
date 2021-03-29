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


#include "heap.h"

/*********************************************************/

int *initPerm (int size)
{
	int *p;
	int i;

	p = (int *) mCalloc (size, sizeof (int));

	for (i=0; i<size; i++)
		p[i] = i;

	return (p);
}

/*********************************************************/

void permInverse (int *p, int *q, int length)
{
	int i;

	for (i=0; i<length; i++)
		q[p[i]] = i;

	return;
}

/* swaps two values of a permutation */

/*********************************************************/

void swap (int *p, int *q, int i, int j)
{
	int temp;

	temp = p[i];
	p[i] = p[j];
	p[j] = temp;
	q[p[i]] = i;
	q[p[j]] = j;

	return;
}

/* The usual Heapify function, tailored for our use with a heap of scores
 * will use array p to keep track of indexes
 * after scoreHeapify is called, the subtree rooted at i will be a heap
 * p goes from heap to array, q goes from array to heap */

/*********************************************************/

void heapify (int *p, int *q, double *HeapArray, int i, int n)
{
	int left, right, smallest;
	boolean moreswap = TRUE;

	do
	{
		left = 2 * i;
		right = 2* i + 1;
		if ((left <= n) && (HeapArray[p[left]] < HeapArray[p[i]]))
			smallest = left;

		else
			smallest = i;

		if ((right <= n) && (HeapArray[p[right]] < HeapArray[p[smallest]]))
			smallest = right;

		if (smallest != i)
		{
			swap (p, q, i, smallest);	/* push smallest up the heap */
			i = smallest;	  			/* check next level down */
		}
		else
			moreswap = FALSE;
	} while (moreswap);

	return;
}

/* heap is of indices of elements of v, popHeap takes the index at
 * position i and pushes it out of the heap (by pushing it to the bottom
 * of the heap, where it is not noticed) */

/*********************************************************/

void reHeapElement (int *p, int *q, double *v, int length, int i)
{
	int up, here;

	here = i;
	up = i / 2;

	if ((up > 0) && (v[p[here]] < v[p[up]]))
	{
		while ((up > 0) && (v[p[here]] < v[p[up]]))		/* we push the new
														 * value up the heap */
		{
			swap (p, q, up, here);
			here = up;
			up = here / 2;
		}
	}
	else
		heapify (p, q, v, i, length);

	return;
}

/*********************************************************/

void popHeap (int *p, int *q, double *v, int length, int i)
{
	swap (p, q, i, length);					/* puts new value at the last
											 * position in the heap */
	reHeapElement (p, q, v, length-1, i);	/* put the swapped guy in the
											 * right place */
	return;
}

/*********************************************************/

void pushHeap (int *p, int *q, double *v, int length, int i)
{
	swap (p, q, i, length+1);						/* puts new value at the last
													 * position in the heap */
	reHeapElement (p, q, v, length+1, length+1);	/* put that guy in the
													 * right place */
  return;
}

/*********************************************************/

int makeThreshHeap (int *p, int *q, double *v, int arraySize, double thresh)
{
	int i, heapsize;

	heapsize = 0;

	for (i = 1; i<arraySize; i++)
		if (v[q[i]] < thresh)
			pushHeap (p, q, v, heapsize++, i);

	return (heapsize);
}

