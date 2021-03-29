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


#include "fastme.h"

boolean isBoostrap;
int verbose;

/*********************************************************/

int main (int argc, char **argv)
{
	Options *options;
	set *species, *species_bk;
	
	species = species_bk = NULL;
	
	time_t t_beg, t_end;

	tree *T = NULL;

	seq **sequences = NULL;

#ifdef _OPENMP
	int i = 0;
#endif

	int numSpecies;
	int setCounter = 0;
	
	int nniCount, sprCount, repCounter, repToPrint, printedRep;
	
	int seqLength = -1;
	
	// Random numbers for bootstraps
	int **rnd = NULL;

	double **D, **A;
	D = A = NULL;

	// DNA sequences from alignment
	char **DNAdata = NULL;
	
	// PROTEIN sequences from alignment
	model *PROTEINdatamodel = NULL;
	allseq *PROTEINdataseq = NULL;

	// Strings containing the best tree and the bootstraped trees
	char *bestTree = NULL;
	char **bootTrees = NULL;
	
	// String to record pseudo-matrices
	char **matStr = NULL;
	
	// Input tree file stream position usefull for bootstraps
	fpos_t fpI_tree_pos;

	isBoostrap = FALSE;

	options = chooseSettings (argc, argv);

//DEBUG
//	PrintOptions (options);

#ifdef _OPENMP
	Message ( (char*)"This analysis will run on %d threads.", options->nb_threads);
#else
	printf ("\n");
	fflush (stdout);
	fflush (stderr);
#endif

	time(&t_beg);

	OpenFiles (options);

	// mem alloc for best tree
	bestTree = (char *) mCalloc (MAX_INPUT_SIZE, sizeof (char));
	if (options->nb_bootstraps > 0)
	{
		bootTrees = (char **) mCalloc (options->nb_bootstraps, sizeof (char *));
	}

	if (PROTEIN == options->input_type && PDIST != options->model)
	{
		PROTEINdatamodel = InitProtModel (options);
	}

	sgenrand ( (unsigned long)options->seed);

	while ( setCounter < options->nb_datasets )
	{
		//nniCount = sprCount = tbrCount = repCounter = 0;
		nniCount = sprCount = repCounter = 0;
		setCounter++;
		species = (set *) mCalloc (1, sizeof(set));
		InitSpeciesAndTrees (options, species, bootTrees, bestTree);
		
		printf ("\n#  Analysing dataset %d\n", setCounter);

		if (setCounter == 1)
			printOptions (options);

		fprintf (options->fpO_stat_file,"Dataset %d\n", setCounter);
		if (options->use_NNI && options->use_SPR)
		{
			fprintf (options->fpO_stat_file,"\tTwo tree searches are performed, with NNIs and SPRs respectively.\n");
			fprintf (options->fpO_stat_file,"\tThe resulting tree is the best (shortest) of both.\n\n");
		}

		if (options->nb_datasets > 1)
		{
			if (options->nb_bootstraps > 0)
			{
				if (setCounter > 1)
					fprintf (options->fpO_boot_file,"\n");

				fprintf (options->fpO_boot_file,"Dataset %d\n", setCounter);
			}
			if (options->use_O_mat_file)
				fprintf (options->fpO_mat_file,"Dataset %d\n", setCounter);
		}

/*********************************************************
		GET MATRIX
**********************************************************/
		if (MATRIX == options->input_type)
		{
			D = loadM (options->fpI_data_file, &numSpecies, species);
			PrintEstimatedMemorySpace (numSpecies, 0, options);
		}
/*********************************************************
		GET DATA FROM SEQUENCES
**********************************************************/
		else
		{
			// read sequences from input file
			sequences = Get_Seq (options->fpI_data_file, options->is_interleaved,
							&numSpecies, &seqLength, options->input_type, species);

			// DNA sequences
			if (DNA == options->input_type || PDIST == options->model || F81LIKE == options->model)
			{
				DNAdata = GetDataFromDNA (numSpecies, sequences);
			}
			// PROTEIN sequences
			else if (PROTEIN == options->input_type)
			{
				PROTEINdataseq = GetDataFromProt (options, numSpecies, sequences, PROTEINdatamodel);
			}

			Free_Seq (sequences, numSpecies);
			PrintEstimatedMemorySpace (numSpecies, seqLength, options);

/*********************************************************
		GET MATRIX FROM ALIGNMENT DATA
**********************************************************/
			Message ( (char*)"Computing pairwise distances...");
			// DNA sequences
			if (DNA == options->input_type || PDIST == options->model || F81LIKE == options->model)
			{
				D = GetMatFromDNA (options, numSpecies, seqLength, DNAdata, NULL);
			}
			// PROTEIN sequences
			else if (PROTEIN == options->input_type)
			{
				D = GetMatFromProt (PROTEINdataseq, PROTEINdatamodel, NULL, options->nb_threads, options->global_aa_fq);
			}
			else D = NULL;
		}

		species_bk = copySet (species);
		
		if (options->use_O_mat_file)
			printMatrix (D, numSpecies, species, options->fpO_mat_file, options->input_type, options->precision);

/*********************************************************
		COMPUTE TREE
**********************************************************/
		if (numSpecies < 4)
			Warning ( (char*)"Cannot compute tree with less than 4 taxa.");
		if (! options->only_mat && numSpecies > 3)
		{
			if (USER == options->method)
				Message ( (char*)"Reading input tree...");
			else
				Message ( (char*)"Computing tree...");

			A = initDoubleMatrix (2*numSpecies-2);
			
			if (USER == options->method)
			{
				if (repCounter == 0)
					fgetpos (options->fpI_tree_file, &fpI_tree_pos);
				else
					fsetpos(options->fpI_tree_file, &fpI_tree_pos);
				T = loadNewickTree (options->fpI_tree_file, numSpecies);
				T = detrifurcate (T);
				compareSets (T, species);
				partitionSizes (T);
			}
			else
			{
				T = ComputeTree (options, D, A, species, numSpecies, options->precision);
			}
		
			T = ImproveTree (options, T, D, A, &nniCount, &sprCount, options->fpO_stat_file);
		
			if (verbose > 0)
			{
				Message ( (char*)"Performed %d NNI(s) and %d SPR(s) on dataset %d.",
						nniCount, sprCount, setCounter);
			}
		}
		
/*********************************************************
		BOOTSTRAPS
**********************************************************/
		if (options->nb_bootstraps > 0 && numSpecies > 3)
		{
			Message ( (char*)"Non parametric bootstrap analysis...");
			printf ("\n  [");
			isBoostrap = TRUE;
			NewickPrintTreeStr (T, bestTree, options->precision);

			repCounter = repToPrint = printedRep = 0;
			
			if (options->use_O_mat_file)
			{
				// mem alloc for pseudo-matrices recorded in strings
				matStr = InitMatStrings (options->nb_bootstraps, numSpecies);
			}
			
			// Generate random numbers for bootstraps before parallel section
			rnd = rndForBootstraps (options, seqLength);

#ifdef _OPENMP
			// Create a proteic model for each thread
			model **protModels = NULL;
			if (PROTEIN == options->input_type)
			{
				protModels = (model **) mCalloc (options->nb_threads, sizeof (model *));
				for (i=0; i<options->nb_threads; i++)
					protModels[i] = Copy_Model (PROTEINdatamodel);
			}

			freeMatrix (A, 2*numSpecies-2);
			
	#pragma omp parallel for private (repCounter, species, D, A, T, nniCount, sprCount) shared (bootTrees, matStr, repToPrint)
#endif
			for (repCounter = 1; repCounter <= options->nb_bootstraps; repCounter++)
			{
				repToPrint++;
				printedRep = PrintBootstrapInfo (options, repToPrint, printedRep);

				//nniCount = sprCount = tbrCount = 0;
				nniCount = sprCount = 0;
				species = copySet (species_bk);


		/***GET MATRIX FROM ALIGNMENT DATA********************/
				if (MATRIX != options->input_type)
				{
					if (DNA == options->input_type || PDIST == options->model || F81LIKE == options->model)
					{
						D = GetMatFromDNA (options, numSpecies, seqLength, DNAdata, rnd[repCounter-1]);
					}
					else if (PROTEIN == options->input_type)
					{
#ifdef _OPENMP
						D = GetMatFromProt (PROTEINdataseq, protModels[omp_get_thread_num()], rnd[repCounter-1], options->nb_threads, options->global_aa_fq);
#else
						D = GetMatFromProt (PROTEINdataseq, PROTEINdatamodel, rnd[repCounter-1], options->nb_threads, options->global_aa_fq);
#endif
					}
					else D = NULL;
				}

#ifdef _OPENMP
				A = initDoubleMatrix (2*numSpecies-2);
#else
				fillZeroMatrix (&A, 2*numSpecies-2);
#endif

		/***COMPUTE TREE*************************************/
				if (USER == options->method)
				{
#ifdef _OPENMP
	#pragma omp critical (readUserTree)
	{
#endif
					fsetpos(options->fpI_tree_file, &fpI_tree_pos);
					T = loadNewickTree (options->fpI_tree_file, numSpecies);
#ifdef _OPENMP
	}	// End readUserTree
#endif
					T = detrifurcate (T);
					compareSets (T, species);
					partitionSizes (T);
				}
				else
				{
					T = ComputeTree (options, D, A, species, numSpecies, options->precision);
				}

				T = ImproveTree (options, T, D, A, &nniCount, &sprCount, options->fpO_stat_file);
				
#ifdef _OPENMP
				freeMatrix (A, 2*numSpecies-2);
#endif

				NewickPrintTreeStr (T, bootTrees[repCounter-1], options->precision);
				
				if (options->use_O_mat_file)
				{
					freeSet (species);
					species = copySet (species_bk);
					printMatrixStr (D, numSpecies, species, matStr[repCounter-1], options->input_type, options->precision);
				}
				
			} // End bootstrap loop

			printedRep = PrintBootstrapInfo (options, options->nb_bootstraps, printedRep);
			printf ("] %d/%d\n", printedRep, options->nb_bootstraps);

#ifdef _OPENMP
			// Delete threaded proteic models
			if (PROTEIN == options->input_type)
			{
				for (i=0; i<options->nb_threads; i++)
					free (protModels[i]);
				free (protModels);
			}
#endif

			isBoostrap = FALSE;
			boot (bestTree, bootTrees, options->nb_bootstraps, options->fpO_tree_file);
			
			printFinalData (options, bootTrees, matStr);
			
			freeIntMat (rnd, options->nb_bootstraps);
	
		}
		else
		{
			if (! options->only_mat && numSpecies > 3)
				NewickPrintTree (T, options->fpO_tree_file, options->precision);
		}

#ifndef _OPENMP
		freeMatrix (A, 2*numSpecies-2);
#endif

		if (numSpecies > 3)
		{
			if (NULL != T)
				freeTree (T);
			if (! options->only_mat)
				fprintf (options->fpO_tree_file, "\n");
		}
		
		if (NULL != matStr)
			free (matStr);
		
		freeMatrix (D, numSpecies);
		freeSet (species);
		freeSet (species_bk);
		if (NULL != DNAdata)
			freeCharMatrix (DNAdata, numSpecies);
		
	} //end datasets loop

	fflush (stdout);
	fflush (stderr);
	
	if (NULL != rnd)		free (rnd);
	if (NULL != bootTrees)	free (bootTrees);
	if (NULL != bestTree)	free (bestTree);

	if (NULL != PROTEINdataseq)
		Free_Cseq (PROTEINdataseq);

	Free_Model (PROTEINdatamodel);
	Free_Input (options);

	time (&t_end);
	PrintTimeInfo (t_beg, t_end);
	
	//exit (EXIT_SUCCESS);

	return 0;
}

/*********************************************************/

void printFinalData (Options *options, char **bootTrees, char **matStr)
{
	int i;
	
	for (i=0; i<options->nb_bootstraps; i++)
	{
		// print pseudo-trees to file
		// free memory allocated for best tree and pseudo-trees strings
		fprintf (options->fpO_boot_file, "%s", bootTrees[i]);
		if (NULL != bootTrees[i])
			free (bootTrees[i]);
		
		// print pseudo-matrices to file
		// free memory allocated for pseudo-matrices strings
		if (options->use_O_mat_file)
		{
			fprintf (options->fpO_mat_file, "%s", matStr[i]);
			if (NULL != matStr[i])
				free (matStr[i]);
		}
	}
	
	return;
}

/*********************************************************/

char **InitMatStrings (int numBoot, int numSpc)
{
	int i;
	
	char **matStr = (char **) mCalloc (numBoot, sizeof (char *));
	
	for (i=0; i<numBoot; i++)
	{
		matStr[i] = (char *) mCalloc ( (numSpc+2) + (13 * numSpc * (numSpc+1)), sizeof (char));
	}
	
	return matStr;
}

/*********************************************************/

int **rndForBootstraps (Options *options, int len)
{
	int **P = NULL;
	
	if (DNA == options->input_type || PDIST == options->model)
	{
		P = bootFilter (options->nb_bootstraps, len);
	}
	else if (PROTEIN == options->input_type)
	{
		P = p_bootPositions (options->nb_bootstraps, len);
	}
	
	return P;
}

/*********************************************************/

int **bootFilter (int nBoot, int len)
{
	int i;
	int **filters;

	filters = (int **) mCalloc (nBoot, sizeof (int *));
	
	for (i=0; i<nBoot; i++)
	{
		filters[i] = initZeroArray (len);
		bootstrapSelect (len, filters[i]);
	}

	return filters;
}

/*********************************************************/

int **p_bootPositions (int nBoot, int len)
{
	int i, j;
	int **pos;

	pos = (int **) mCalloc (nBoot, sizeof (int *));
	
	for (i=0; i<nBoot; i++)
	{
		pos[i] = (int *) mCalloc (len, sizeof (int));
		for (j=0; j<len; j++)
		{
			pos[i][j] = getIntRandom (len);
		}
	}

	return (pos);
}

/*********************************************************/

void freeIntMat (int **mat, int size)
{
	int i;

	for (i=0; i<size; i++)
		if (NULL != mat[i])
			free(mat[i]);
	
	return;
}

/*********************************************************/

allseq *p_bootstraps (allseq *alldata, int ns, int *site_num, int *positions)
{
	int i, position;
	int init_len = 0;

	allseq *boot_data = Copy_Cseq (alldata, alldata->crunch_len, ns);

	for (i=0; i<boot_data->crunch_len; i++)
		boot_data->wght[i] = 0;

	if (isBoostrap && (NULL != positions))
	{
		for (i=0; i<boot_data->init_len; i++)
		{
			position = positions[i];
			boot_data->wght[site_num[position]] += 1;
			init_len++;
		}
	}
	else
	{
		for (i=0; i<boot_data->init_len; i++)
		{
			position = getIntRandom (alldata->init_len);
			boot_data->wght[site_num[position]] += 1;
			init_len++;
		}
	}

	if (init_len != alldata->init_len)
		Exit ( (char*)"Problem when copying sequences for bootstrap.");

	Get_AA_Freqs (boot_data);
	
	return (boot_data);
}

/*********************************************************/

void OpenFiles (Options *options)
{
	options->fpI_data_file = Openfile (options->I_data_file,  (char*)"r");
	options->fpO_stat_file = Openfile (options->O_stat_file, options->open_mode);

	if (! options->only_mat)
		options->fpO_tree_file = Openfile (options->O_tree_file, options->open_mode);
	
	if (options->nb_bootstraps > 0)
		options->fpO_boot_file = Openfile (options->O_boot_file, options->open_mode);

	if (options->use_O_mat_file)
		options->fpO_mat_file = Openfile (options->O_mat_file, options->open_mode);

	if (USER == options->method)
		options->fpI_tree_file = Openfile (options->I_tree_file, (char*)"r");

  return;
}

/*********************************************************/

void printOptions (Options *options)
{
	char *tmp;
	tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));

	FILE *f = options->fpO_stat_file;

	if (options->open_mode[0] == 'w')
	{
		fprintf (f, "\n - FastME %s - \n\n", PACKAGE_VERSION);
		if (!options->only_mat && (options->method == NJ || options->method == UNJ || options->method == BIONJ))
			fprintf (f, "\nPapers to be cited:\n");
		else
			fprintf (f, "\nPaper to be cited:\n");
		
		fprintf (f, "\nFastME 2.0 - A comprehensive, accurate and fast distance-based phylogeny inference program.");
		fprintf (f, "\n\tVincent Lefort, Richard Desper and Olivier Gascuel,");
		fprintf (f, "\n\tMolecular Biology and Evolution 32(10), 2798-800, 2015.");
		
		if (!options->only_mat)
		{
			if (options->method == BIONJ)
			{
				fprintf (f, "\nBIONJ algorithm:");
				fprintf (f, "\n\tGascuel O. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data.");
				fprintf (f, "\n\tMolecular Biology and Evolution, 14(7):685-695");
			}
			if (options->method == NJ)
			{
				fprintf (f, "\nNJ algorithm:");
				fprintf (f, "\n\tSaitou N., Nei M. 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees.");
				fprintf (f, "\n\tMolecular Biology and Evolution, 4(4):406-25");
			}
			if (options->method == UNJ)
			{
				fprintf (f, "\nUNJ algorithm:");
				fprintf (f, "\n\tGascuel O. 1997. Concerning the NJ algorithm and its unweighted version, UNJ.");
				fprintf (f, "\n\tMathematical Hierarchies and Biology,");
				fprintf (f, "\n\tB. Mirkin, F.R. McMorris, F.S. Roberts and A. Rzetsky (eds.),");
				fprintf (f, "\n\tAmerican Mathematical Society, Providence, 149-170");
			}
		}
		
		fprintf (f, "\n\n-------------------------------------------------------------------------------\n");
		fprintf (f, "Settings for this run:\n\n");

//----------------------------------------------------------------------//

		constantToStr (options->input_type, tmp);
		fprintf (f, "  I "
			"                                     Input data type "
			" %-15s \n", tmp);

//----------------------------------------------------------------------//

		if (options->input_type != MATRIX)
		{
			constantToStr (options->model, tmp);
			fprintf (f, "  E "
				"                                  evolutionary model "
				" %-15s \n", tmp);
/*
			if (options->input_type == DNA)
				fprintf (f, "  E "
					"                              DNA evolutionary model "
					" %-15s \n", tmp);

			else if (options->input_type == PROTEIN)
				fprintf (f, "  E "
					"                          PROTEIN evolutionary model "
					" %-15s \n", tmp);
*/

//----------------------------------------------------------------------//

			if (options->use_gamma)
			{
				fprintf (f, "  G "
					"                   Gamma rate variation across sites "
					" %-15f \n", options->gamma);
			}
			else
			{
				fprintf (f, "  G "
					"                   Gamma rate variation across sites "
					" %-15s \n", "no");
			}

//----------------------------------------------------------------------//

			fprintf (f, "  R "
				"                             Remove sites whith gaps "
				" %-15s \n", (options->no_gap ? "yes" : "no"));

//----------------------------------------------------------------------//

			fprintf (f, "  O "
				"                   Output calculated distance matrix "
				" %-15s \n", (options->use_O_mat_file ? "yes" : "no"));

			fprintf (f, "\n");
		}

//----------------------------------------------------------------------//

		if (! options->only_mat)
		{
			fprintf (f, "  D "
				"                                  Number of datasets "
				" %-15d \n", options->nb_datasets);

//----------------------------------------------------------------------//
			constantToStr (options->method, tmp);
			fprintf (f, "  M "
				"                        Initial tree building method "
				" %-15s \n", tmp);

//----------------------------------------------------------------------//

			if (options->use_NNI)
			{
				constantToStr (options->NNI, tmp);
				fprintf (f, "  N "
					"                                  NNI postprocessing "
					" %-15s \n", tmp);
			}
			else
			{
				fprintf (f, "  N "
					"                                  NNI postprocessing "
					" %-15s \n", "no");
			}

//----------------------------------------------------------------------//

			fprintf (f, "  S "
				"                                  SPR postprocessing "
				" %-15s \n", (options->use_SPR ? "yes" : "no"));

//----------------------------------------------------------------------//

			if (!options->use_NNI && !options->use_SPR && options->branch != NONE)
			{
				constantToStr (options->branch, tmp);
				fprintf (f, "  W "
					"             Branch lengths assigned to the topology "
					" %-15s \n", tmp);
			}

//----------------------------------------------------------------------//
/*
			fprintf (f, "  T "
				"                                  TBR postprocessing "
				" %-15s \n", (options->use_TBR ? "yes" : "no"));*/

//----------------------------------------------------------------------//

			if (options->nb_bootstraps > 0)
			{
				fprintf (f, "\n");
				fprintf (f, "  B "
					"                     Bootstrap: number of replicates "
					" %-15d \n", options->nb_bootstraps);
			}

		}

		fprintf (f, "\n-------------------------------------------------------------------------------\n");
	}

	free (tmp);

	return;
}

/*********************************************************/

void printMatrix (double **D, int size, set *nodes, FILE *ofile, int input_type, int precision)
{
	char *tmp;
	char format[8];
	int i;
	node *v, *w;
	set *S, *T;
	double threshold, d;

	tmp = (char *) mCalloc (DECIMAL_DIG, sizeof(char));
	snprintf (format, 8, "%%.%df  ", precision);

	if (PROTEIN == input_type)
		threshold = PROT_DIST_MAX;
	else
		threshold = DNA_DIST_MAX;

	fprintf (ofile, "%d\n", size);
	for (S=nodes; NULL!=S; S=S->secondNode)
	{
		v = S->firstNode;
		fprintf (ofile, "%-10s ", v->label);
		for (T=nodes; NULL!=T; T=T->secondNode)
		{
			w = T->firstNode;
			d = D[v->index2][w->index2];
			if (d > threshold)
			{
				// Build string of length depending on 'precision' to write
				strncat (tmp, "NA", 2);
				// Add [precision - 2] blank spaces
				for (i=0; i<precision-2; i++)
					strncat (tmp, " ", 1);
				// The following code generates warning because 'tmp' is not a string literal
				fprintf (ofile, tmp);
			}
			else
				fprintf (ofile, format, d);
		}
		fprintf (ofile, "\n");
	}
	fprintf (ofile, "\n");
	free (tmp);

	return;
}

/*********************************************************/

void printMatrixStr (double **D, int size, set *nodes, char *str, int input_type, int precision)
{
	int i;
	node *v, *w;
	set *S, *T;
	char *tmp;
	char format[8];
	double threshold, d;

	snprintf (format, 8, "%%.%df  ", precision);

	if (PROTEIN == input_type)
		threshold = PROT_DIST_MAX;
	else
		threshold = DNA_DIST_MAX;
	
	tmp = (char *) mCalloc (INPUT_SIZE, sizeof(char));
	if (strlen (tmp))
		strncpy (tmp, "", strlen(tmp));
	snprintf (tmp, INPUT_SIZE, "%d\n", size);
	strncat (str, tmp, strlen (tmp));

	for (S=nodes; NULL!=S; S=S->secondNode)
	{
		v = S->firstNode;
		snprintf (tmp, INPUT_SIZE, "%-10s ", v->label);
		strncat (str, tmp, strlen (tmp));
		for (T=nodes; NULL!=T; T=T->secondNode)
		{
			w = T->firstNode;
			d = D[v->index2][w->index2];
			if (d > threshold)
			{
				//snprintf (tmp, INPUT_SIZE, "NA          ");
				snprintf (tmp, 2, "NA");
				// Add [precision - 2] blank spaces
				for (i=0; i<precision-2; i++)
					strncat (tmp, " ", 1);
			}
			else
			{
				//snprintf (tmp, INPUT_SIZE, "%7.8f  ", d);
				snprintf (tmp, INPUT_SIZE, format, d);
			}
			strncat (str, tmp, strlen (tmp));
		}
		strcat (str, "\n");
	}
	strcat (str, "\n");
	free (tmp);

	return;
}

/*********************************************************/

void PrintTimeInfo (time_t t_beg, time_t t_end)
{
	ldiv_t hour, min;
	int sec;

	hour = ldiv (t_end - t_beg, 3600);
	min = ldiv (t_end - t_beg, 60);
	min.quot -= hour.quot * 60;
	sec = (int) (t_end-t_beg) % 60;

	if (min.quot > 9)
	{
		if (sec > 9)
			Message ( (char*)"Time used %dh%dm%ds", hour.quot, min.quot, sec);
		else
			Message ( (char*)"Time used %dh%dm0%ds", hour.quot, min.quot, sec);
	}
	else
	{
		if (sec > 9)
			Message ( (char*)"Time used %dh0%dm%ds", hour.quot, min.quot, sec);
		else
			Message ( (char*)"Time used %dh0%dm0%ds", hour.quot, min.quot, sec);
	}

	return;
}

/*********************************************************/

void PrintEstimatedMemorySpace (int nbtax, int nbsites, Options *options)
{
	unsigned long long int schar, sint, ssint, sdouble, sphydbl, nbT, nbS;
	unsigned long long int mem, memMat, memTree, memNNI, memSPR, Matsize, Treesize;
	lldiv_t ko;
	ldiv_t Mo;
	ldiv_t Go;
	
	memMat = memTree = memNNI = memSPR = 0;
	nbT = (unsigned long long int) nbtax;
	nbS = (unsigned long long int) nbsites;
	schar = (unsigned long long int) sizeof (char);
	sint = (unsigned long long int) sizeof (int);
	ssint = (unsigned long long int) sizeof (short int);
	sdouble = (unsigned long long int) sizeof (double);
	sphydbl = (unsigned long long int) sizeof (phydbl);

	Matsize = nbT * nbT;
	Treesize = (nbT * 2) - 2;

	/**
	 ** Estimate memory space needed to compute distance matrix
	 */

	if (options->input_type == DNA)
	{
		memMat += (nbT * sdouble);
		memMat += (nbT * ((schar * nbS) + (10 * sint) + (4 * sdouble)));
		memMat += (Matsize * ( (7 * sdouble) + sint) );
		memMat /= 8;
	}
	else if (options->input_type == PROTEIN)
	{
		// Egein struct
		memMat += (4 * 20 * sdouble);
		memMat += (2 * 20 * sint);
		memMat += (4 * 20 * 20 * sdouble);

		// Pij Qmat
		memMat += (3 * 20 * 20 * sdouble);

		// Sequences
		memMat += ( nbT * MAX_NAME_LENGTH * schar );
		memMat += ( nbT * nbS * (schar + ssint) );
		memMat += ( nbS * 3 * ssint );
		memMat += ( nbS * sint );

		// Matrix
		memMat += ( nbT * MAX_NAME_LENGTH * schar );
		memMat += ( Matsize * 4 * sphydbl );
	}

#ifdef _OPENMP
	memMat = memMat * ( (unsigned long long int) options->nb_threads );
#endif

	// String length for trees
	memMat += ( (1 + (unsigned long long int) options->nb_bootstraps) * MAX_INPUT_SIZE * schar );

	/**
	 ** Estimate memory space needed to compute tree
	 */

	if (options->use_NNI)
	{
		// Size of 1 edge ~ 2 double + 2 int
		memNNI = (Treesize + 1) * 2 * 2 * sdouble;
		memNNI += ( (Treesize + 1) * 2 * sint );
	}
	
	if (options->use_SPR)
		memSPR = 2 * Treesize * Treesize * sdouble;
/*
	if (options->use_TBR)
	{
		//memTree += Treesize * Treesize * Treesize * sdouble;
		memTree = Treesize * Treesize * Treesize * sdouble;
	}
*/
	memTree = (memNNI > memSPR) ? memNNI : memSPR;
	mem = (memMat > memTree) ? memMat : memTree;

	ko = lldiv ( (long long) mem, 1024);
	Mo = ldiv (ko.quot, 1024);
	Go = ldiv (Mo.quot, 1024);

	if (Go.quot > 0)
		Warning ( (char*)"This analysis requires at least %d.%d Go of memory space.", Go.quot, (ldiv (Go.rem, 100)).quot );

	else if (Mo.quot > 100)
		Message ( (char*)"This analysis requires at least %d Mo of memory space.", Mo.quot);

	return;
}

/*********************************************************/

int PrintBootstrapInfo (Options *options, int repCounter, int printedRep)
{
	//int ID = 0;
	
#ifdef _OPENMP
	//ID = omp_get_thread_num();
	#pragma omp master
#endif
	// With parallel version : execute only on master thread
	//if (isBoostrap && ID == 0)
	if (isBoostrap)
	{
		while (printedRep < repCounter)
		{
			if (printedRep > 0 && printedRep < options->nb_bootstraps && 0 == (printedRep) % 20)
				printf ("] %d/%d\n  [", printedRep, options->nb_bootstraps);
			printf (".");
			printedRep++;
		}
	}

	fflush (stdout);
	fflush (stderr);

	return printedRep;
}

/*********************************************************/

model *InitProtModel (Options *options)
{
	model *mod = Make_Model_Basic();
	Set_Defaults_Model (mod);
	mod->whichmodel = options->model;
	Make_Model_Complete (mod, options);
	
	return mod;
}

/*********************************************************/

void InitSpeciesAndTrees (Options *options, set *species, char **bootTrees, char *bestTree)
{
	int i = 0;
	
	species->firstNode = NULL;
	species->secondNode = NULL;

	size_t n = sizeof(bestTree);
	memset (bestTree, '\0', n);
	
	// mem alloc for bootstraped trees strings
	for (i=0; i<options->nb_bootstraps; i++)
	{
		bootTrees[i] = (char *) mCalloc (MAX_INPUT_SIZE, sizeof (char));
	}

	return;
}

/*********************************************************/

char **GetDataFromDNA (int numSpecies, seq **sequences)
{
	int i;

	char **DNAdata = (char **) mCalloc (numSpecies, sizeof (char *));
	
	for (i=0; i<numSpecies; i++)
	{
		DNAdata[i] = (char *) mCalloc (sequences[i]->len, sizeof (char));
		strncpy (DNAdata[i], sequences[i]->state, (unsigned long) sequences[i]->len);
	}

	return DNAdata;
}

/*********************************************************/

double **GetMatFromDNA (Options *options, int numSpecies, int seqLength,
	char **DNAdata, int *filter)
{
	boolean filterIsNULL;
	
	if (NULL == filter)
	{
		filterIsNULL = TRUE;
		if (isBoostrap)
		{
			filter = initZeroArray (seqLength);
			bootstrapSelect (seqLength, filter);
		}
		else
		{
			filter = initOneArray (seqLength);
		}
	}
	else
	{
		filterIsNULL = FALSE;
	}
	double **D = makeDistMatrix (DNAdata, numSpecies, seqLength, options->use_gamma,
						options->gamma, options->model, options->input_type, filter,
						options->no_gap, options->fpO_stat_file, options->use_O_mat_file);
	
	if (filterIsNULL)
	{
		free (filter);
		filter = NULL;
	}
	
	return D;
}

/*********************************************************/

allseq *GetDataFromProt (Options *options, int numSpecies, seq **sequences,
	model *mod)
{
	int i, j;
	allseq *alldata = NULL;
	
	mod->n_otu = numSpecies;
	alldata = Compact_Seq (sequences, mod, options->no_gap);
	Check_Ambiguities (alldata, mod->stepsize);
	Init_Model (alldata, mod, options->global_aa_fq);

	mod->site_num_len = alldata->init_len;
	mod->site_num = (int *) mCalloc (mod->site_num_len, sizeof(int));
	mod->n_site = 0;
	for (i=0; i<alldata->crunch_len; i++)
	{
		for (j=0; j<alldata->wght[i]; j++)
		{
			mod->site_num[mod->n_site] = i;
			mod->n_site++;
		}
	}

	return alldata;
}

/*********************************************************/

double **GetMatFromProt (allseq *alldata, model *mod, int *positions,
	int nbthreads, boolean global_aa_fq)
{
	matrix *mat = NULL;
	double **D = initDoubleMatrix (mod->n_otu);

	if (isBoostrap)
	{
		allseq *boot_data = p_bootstraps (alldata, mod->ns, mod->site_num, positions);
		Init_Model (boot_data, mod, global_aa_fq);
		mat = ML_Dist (boot_data, mod, nbthreads);
		Free_Cseq (boot_data);
	}
	else
	{
		mat = ML_Dist (alldata, mod, nbthreads);
	}
	
	Fill_Missing_Dist (mat);
	D = Copy_PMat_to_DMat (mat);
	Free_Mat (mat);
	
	return (D);
}

/*********************************************************/

tree *ComputeTree (Options *options, double **D, double **A, set *species,
	int numSpecies, int precision)
{
	set *slooper;
	node *addNode;
	tree *T = NULL;
	char format[8];
	snprintf (format, 8, "%%.%df", precision);
	
	switch (options->method)
	{
		case USER:
			break;

		case TaxAddBAL:
			for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
			{
				addNode = copyNode (slooper->firstNode);
				T = BMEaddSpecies (T, addNode, D, A);
			}
			assignBMEWeights (T, A);
			break;

		case TaxAddOLS:
			for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
			{
				addNode = copyNode (slooper->firstNode);
				T = GMEaddSpecies (T, addNode, D, A);
			}
			makeOLSAveragesTable (T, D, A);
			assignOLSWeights (T, A);
			break;
			
		case NJ:
			T = bionj (D, species, numSpecies, TRUE, format);
			compareSets (T, species);
			partitionSizes (T);
			break;

		case BIONJ:
			T = bionj (D, species, numSpecies, FALSE, format);
			compareSets (T, species);
			partitionSizes (T);
			break;

		case UNJ:
			T = unj (D, species, numSpecies, format);
			compareSets (T, species);
			partitionSizes (T);
			break;
	}
	
	return T;
}

/*********************************************************/

tree *ImproveTree (Options *options, tree *T0, double **D, double **A,
	int *nniCount, int *sprCount, FILE *ofile)
{
	//T0 = ME
	//T1 = ME + NNI
	//T2 = ME + SPR
	//T3 = ME + NNI + SPR
	tree *T1, *T2, *T3;
	
	T1 = T2 = T3 = NULL;

	if (!options->use_NNI && !options->use_SPR)
	{
		switch (options->branch)
		{
			case BrBAL:
				if (options->method != TaxAddBAL)
					makeBMEAveragesTable (T0, D, A);

				assignBMEWeights (T0, A);
				break;

			case BrOLS:
				if (options->method != TaxAddOLS)
					assignAllSizeFields (T0);

				makeOLSAveragesTable (T0, D, A);
				assignOLSWeights (T0, A);
				break;

			default:
				break;
		}
	}
	
	if (options->use_NNI)
	{
		if (!isBoostrap)
		{
			if (verbose > 2)
				printf ("\n");

			Message ( (char*)"Performing NNI...");
		}
		
		T1 = copyTree (T0);
		
		switch (options->NNI)
		{
			case BALNNI:
				if (options->method != TaxAddBAL)
					makeBMEAveragesTable (T1, D, A);

				bNNI (T1, A, nniCount, options->fpO_stat_file);
				assignBMEWeights (T1, A);
				break;

			case OLSNNI:
				if (options->method != TaxAddOLS)
					assignAllSizeFields (T1);

				makeOLSAveragesTable (T1, D, A);
				NNI (T1, A, nniCount, options->fpO_stat_file);
				assignOLSWeights (T1, A);
				break;

			default:
				break;
		}
		if (!isBoostrap)
			fprintf (ofile, "\tPerformed %d NNI(s).\n\n", *nniCount);
	}


	if (options->use_SPR)
	{
		if (!isBoostrap)
		{
			if (verbose > 2)
				printf ("\n");

			Message ( (char*)"Performing SPR...");
		}
		
		T2 = copyTree (T0);
		
		makeBMEAveragesTable (T2, D, A);
		SPR (T2, D, A, sprCount, options->fpO_stat_file);
		assignBMEWeights (T2, A);
		
		if (!isBoostrap)
			fprintf (ofile, "\tPerformed %d SPR(s).\n\n", *sprCount);
	}
	
	weighTree (T0);
	
	if (NULL != T1)
		weighTree (T1);

	if (NULL != T2)
		weighTree (T2);

	if ((NULL != T1) && (T0->weight > T1->weight))
	{
		freeTree (T0);
		T0 = T1;	//using T0 as the place to store the minimum evolution tree in all cases
		T1 = NULL;
	}
	else if (NULL != T1)
		freeTree (T1);

	if ((NULL != T2) && (T0->weight > T2->weight))
	{
		freeTree (T0);
		T0 = T2;
		T2 = NULL;
	}
	else if (NULL != T2)
		freeTree (T2);
	
	return T0;
}

