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


#ifndef TRAVERSE_H_
#define TRAVERSE_H_

#include "graph.h"

edge *findBottomLeft (edge *e);
edge *findBottomRight (edge *e);
edge *moveRight (edge *e);
edge *moveMiddle (edge *e);
edge *moveRightFirstMiddle (edge *e);
edge *moveLeft (edge *e);
edge *depthFirstTraverse (tree *T, edge *e);
edge *moveUpRight (edge *e);
edge *topFirstTraverse (tree *T, edge *e);
edge *depthRightFirstTraverse (tree *T, edge *e);

#endif /*TRAVERSE_H_*/

