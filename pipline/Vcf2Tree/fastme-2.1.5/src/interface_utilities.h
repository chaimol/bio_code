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


#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "utils.h"


typedef struct __Options {			/* mostly used in 'interface_options.c' */
	char		*I_data_file;		/* required input file containing user data (sequence alignment or distance matrix) */
	char		*I_tree_file;		/* optional input topology file */
	char		*O_tree_file;		/* output file for the resulting tree */
	char		*O_mat_file;		/* optional output file for the distance matrix */
	char		*O_stat_file;		/* output file for execution informations */
	char		*O_boot_file;		/* output file for bootstrapped trees */
	FILE		*fpI_data_file;
	FILE		*fpI_tree_file;
	FILE		*fpO_tree_file;
	FILE		*fpO_mat_file;
	FILE		*fpO_stat_file;
	FILE		*fpO_boot_file;
	boolean		use_O_mat_file;
	char		*open_mode;
	boolean		is_interleaved;		/* input sequence file format (interleaved or sequential) */
	int			nb_datasets;		/* number of datasets (alignments or matrices, default is 1)*/
	int			nb_bootstraps;		/* number of replicates when bootstrapping */
	int			input_type;			/* input file data type (MATRIX, DNA, PROTEIN, SCOREDIST, default is MATRIX) */
	int			method;				/* method for building initial tree (BAL, OLS, NJ, UNJ, BIONJ or USER, default is BAL) */
	int			model;				/* evolutionary model (F84, TN93, K2P, JC69, Transversion Only, LogDet, or SCOREDIST,
									 * default is NONE which corresponds to distance matrix input) */
	boolean		global_aa_fq;		/* equilibrium frequencies computation
									 * for DNA, nt freq are always counted from the input aln
									 * for protein, the default is to use the model freq. If FALSE, the aa freq are counted from the aln */
	boolean		use_gamma;			/* gamma distributed rates across sites */
	float		gamma;				/* gamma rate variation parameter (alpha) */
	boolean		only_mat;			/* only computes the distance matrix */
	int			precision;			/* number of digits after dot to use on output */
	long		seed;				/* seed for randomization (if doing bootstrapping) */
	boolean		no_gap;				/* remove any site which has a gap in any sequence */
	int			branch;				/* branch lengths to assign to a topology (only when no topology improvement, BAL OLS or NONE) */
	boolean		use_NNI;			/* NNI postprocessing */
	int			NNI;				/* type of NNI tree swapping (BAL (default) or OLS) */
	boolean		use_SPR;			/* SPR postprocessing */
//	boolean		use_TBR;			/* TBR postprocessing */
	int			nb_threads;			/* number of threads used in parallel regions (if any) */
} Options;


void PrintOptions (Options *input);
void Free_Input (Options *input);
int Filexists (char *filename);
void Getstring_Stdin (char *file_name);
FILE *Openfile (char *filename, char *mode);
void askOption (char *question, char *c);
boolean testM (char *c);
int getM (char *c);
boolean testN (char *c);
int getN (char *c);
boolean testW (char *c, boolean none);
int getW (char *c);
boolean testD (char *c);
boolean testP (char *c);
boolean testI (char *c);
int getI (char *c);
boolean testF (char *c);
boolean getF (char *c);
int getModel_DNA (char *c);
int getModel_PROTEIN (char *c);

#endif /*UTILITIES_H_*/

