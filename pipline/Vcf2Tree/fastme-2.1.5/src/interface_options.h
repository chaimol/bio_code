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


#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "interface_utilities.h"

#ifndef BOLD
#define BOLD      "\033[00;01m"
#endif

#ifndef FLAT
#define FLAT      "\033[00;00m"
#endif

#ifndef LINE
#define LINE      "\033[00;04m"
#endif


#ifdef _OPENMP
void Usage (int nbthreads);
#else
void Usage (void);
#endif


Options *chooseSettings (int argc, char **argv);
void Set_Defaults_Input (Options *input);
void Get_Input_Interactive (Options *input);
void Get_Input_CommandLine (Options *input, int argc, char **argv);

#endif /*OPTIONS_H_*/
