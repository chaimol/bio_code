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


#ifndef HEAP_H_
#define HEAP_H_

#include "utils.h"

int *initPerm (int size);
void permInverse (int *p, int *q, int length);
void swap (int *p, int *q, int i, int j);
void heapify (int *p, int *q, double *HeapArray, int i, int n);
void reHeapElement (int *p, int *q, double *v, int length, int i);
void popHeap (int *p, int *q, double *v, int length, int i);
void pushHeap (int *p, int *q, double *v, int length, int i);
int makeThreshHeap (int *p, int *q, double *v, int arraySize, double thresh);

#endif /*HEAP_H_*/

