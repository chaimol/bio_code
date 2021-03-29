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


#include "interface_options.h"


/*********************************************************/

Options *chooseSettings (int argc, char **argv)
{
	boolean ask;
	char *tmp;
	char choix;
	int n_trial = 0;

 	Options *input = (Options *) mCalloc (1, sizeof (Options));

	Set_Defaults_Input (input);

	if (argc == 1)
	{
		// Default values (different from CLI)
		input->use_NNI		= TRUE;
		input->use_SPR		= TRUE;
		input->use_gamma	= TRUE;
		input->branch		= BrBAL;
		Get_Input_Interactive (input);
	}
	else
		Get_Input_CommandLine (input, argc, argv);

#ifdef _OPENMP
	if (input->nb_bootstraps > 0 && input->nb_threads > input->nb_bootstraps)
		input->nb_threads = input->nb_bootstraps;
		
	omp_set_num_threads (input->nb_threads);
#endif

	if ((NULL == input->I_data_file) || (! input->I_data_file) || (strlen (input->I_data_file) == 0))
		Exit ( (char*)"You must provide an input file.");

	// ouput files default values
	if ((NULL != input->O_tree_file) && (strlen (input->O_tree_file) == 0))
	{
		strncpy (input->O_tree_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
		strncat (input->O_tree_file, "_fastme_tree.nwk", 16);
	}

	if ((NULL != input->O_stat_file) && (strlen (input->O_stat_file) == 0))
	{
		strncpy (input->O_stat_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
		strncat (input->O_stat_file, "_fastme_stat.txt", 16);
	}

	if ((input->use_O_mat_file) && (NULL != input->O_mat_file) && (strlen (input->O_mat_file) == 0))
	{
		strncpy (input->O_mat_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 16);
		strncat (input->O_mat_file, "_fastme_mat.txt", 15);
	}

	if ((input->nb_bootstraps > 0) && (NULL != input->O_boot_file) && (strlen (input->O_boot_file) == 0)) {
		strncpy (input->O_boot_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
		strncat (input->O_boot_file, "_fastme_boot.txt", 16);
	}

	if ((input->only_mat) && (input->nb_bootstraps > 0)) {
		Warning ( (char*)"Cannot compute bootstraps when only computing distance matrix.");
		input->nb_bootstraps = 0;
	}

	if ((input->input_type == MATRIX) && (input->nb_bootstraps > 0)) {
		Warning ( (char*)"Bootstraps can only be used when inputting aligned sequences.");
		input->nb_bootstraps = 0;
	}

	// No gamma law with p-distance or LogDet
	if (input->model == PDIST || input->model == LOGDET)
		input->use_gamma = FALSE;

	
	// Specific default values for PHYLIP like interface
	if (argc == 1) {
		// BAL or OLS NNI
		if (input->use_NNI)
		{
			if (input->method == TaxAddBAL || input->method == NJ || input->method == BIONJ)
			{
				input->NNI = BALNNI;
			}
			else
			{
				input->NNI = OLSNNI;
			}
		}
		// BAL or OLS Br length
		if (!input->use_NNI && !input->use_SPR)
		{
			if (input->method == TaxAddBAL)
			{
				input->branch = BrBAL;
			}
			else if (input->method == TaxAddOLS)
			{
				input->branch = BrOLS;
			}
			else
				input->branch = NONE;
		}
	}
	// Specific post-treament for CLI
	else
	{
		// Mandatory branch length assignment when no topology improvement with TaxAdd algorithm or user input tree
		if (input->method == TaxAddBAL || input->method == TaxAddOLS || input->method == USER)
		{
			// No tropology improvement choosen
			if (!input->use_NNI && !input->use_SPR)
			{
				// No branch length assignation choosen
				if (input->branch != BrBAL && input->branch != BrOLS)
				{
					// Assign BrOLS branch length with TaxAddOLS
					if (input->method == TaxAddOLS)
					{
						Message ( (char*)"OLS branch lengths will be assigned to the tree.");
						input->branch = BrOLS;
					}
					// Assign BrBAL branch length with TaxAddBAL or user input tree
					else
					{
						Message ( (char*)"BalLS branch lengths will be assigned to the tree.");
						input->branch = BrBAL;
					}
				}
			}
		}
	}

	// Ask question to replace or append files when PHYLIP like interface was used
	if (argc == 1) {
		ask = FALSE;
		tmp = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
		if (Filexists (input->O_tree_file)) {
			strncat (tmp, basename (input->O_tree_file), strlen (basename (input->O_tree_file)));
			ask = TRUE;
		}
		if (Filexists (input->O_stat_file)) {
			if (ask)
				strncat (tmp, "'\n  '", 5);

			strncat (tmp, basename (input->O_stat_file), strlen (basename (input->O_stat_file)));
			ask = TRUE;
		}
		if (Filexists (input->O_boot_file)) {
			if (ask)
				strncat (tmp, "'\n  '", 5);

			strncat (tmp, basename (input->O_boot_file), strlen (basename (input->O_boot_file)));
			ask = TRUE;
		}
		if ((input->use_O_mat_file) && (Filexists (input->O_mat_file))) {
			if (ask)
				strncat (tmp, "'\n  '", 5);

			strncat (tmp, basename (input->O_mat_file), strlen (basename (input->O_mat_file)));
			ask = TRUE;
		}

		if (ask)
		{
			printf ("\n The file(s):\n  '%s'\n already exist(s).\n\n", tmp);
			printf (" Do you want to Replace or Append ? > ");
			n_trial = 0;
			do
			{
				if (n_trial > 0)
					printf ("Please type R or A > ");

				choix = 'X';
				if (scanf ("%c", &choix))
				{
					if (choix == '\n')
						choix = 'X';
				}
				Uppercase (&choix);

		  		if (++n_trial > 10)
		  			Exit ( (char*)"\n");
			}
			while ((choix != 'R') && (choix != 'A'));

			if (choix == 'R')
				strncpy (input->open_mode, "w", 3);

			else
				strncpy (input->open_mode, "a", 3);

		}
		free (tmp);
	}
	
	return input;
}

/*********************************************************/
#ifdef _OPENMP
void Usage (int nbthreads)
#else
void Usage ()
#endif
{
	printf (FLAT"%s\n"BOLD"NAME\n"
		FLAT"\t%s - A comprehensive, accurate and fast distance-based phylogeny inference program.\n\n"
		FLAT"\tVincent Lefort, Richard Desper and Olivier Gascuel,\n"
		FLAT"\tMolecular Biology and Evolution 32(10), 2798-800, 2015.\n"
		FLAT"\tPlease cite this paper if you use this program in your publications.\n", PACKAGE_STRING, PACKAGE_NAME);

	printf (BOLD"\nSYNOPSIS\n"
		BOLD"\t%s"
		FLAT"  ["BOLD"-i "LINE"input data file"FLAT"]"
		FLAT"  ["BOLD"-u "LINE"input user tree file"FLAT"]"
		FLAT"\n\t["BOLD"-o "LINE"output tree file"FLAT"]"
		FLAT"  ["BOLD"-O "LINE"output matrix file"FLAT"]"
		FLAT"  ["BOLD"-I "LINE"output information file"FLAT"]"
		FLAT"\n\t["BOLD"-B "LINE"output bootstrap trees file"FLAT"]"
		FLAT"  ["BOLD"-a"FLAT"]"
		FLAT"\n\t["BOLD"-m "LINE"method"FLAT"]"
		FLAT"  [ "BOLD"-D"LINE"[model]"FLAT" | "BOLD"-P"LINE"[model]"FLAT" ]"
		FLAT"  ["BOLD"-r]"
		FLAT"  ["BOLD"-e]"
		FLAT"  ["BOLD"-g"LINE"[alpha]"FLAT"]"
		FLAT"  ["BOLD"-n"LINE"[NNI]"FLAT"]"
		FLAT"  ["BOLD"-s"FLAT"]"
		FLAT"  ["BOLD"-w "LINE"branch"FLAT"]"
		FLAT"\n\t["BOLD"-d "LINE"datasets"FLAT"]"
		FLAT"  ["BOLD"-b "LINE"replicates"FLAT"]"
		FLAT"  ["BOLD"-z "LINE"seed"FLAT"]"
		FLAT"\n\t["BOLD"-c"FLAT"]"
		FLAT"\n\t["BOLD"-f"FLAT"]"
#ifdef _OPENMP
		FLAT"  ["BOLD"-T "LINE"number of threads"FLAT"]"
#endif
		FLAT"  ["BOLD"-v"FLAT"]"
		FLAT"  ["BOLD"-V"FLAT"]"
		FLAT"  ["BOLD"-h"FLAT"] \n", PACKAGE);

	printf (FLAT"\n\tYou can use %s with no arguments, in this case change the value of\n"
		FLAT"\ta parameter by typing its corresponding character as shown on screen.\n", PACKAGE);

	printf (BOLD"\nOPTIONS\n");
	
	printf (BOLD"\n\t-i "LINE"input data file"BOLD", --input_data="LINE"input data file"
		FLAT"\n\t\tThe "LINE"input data file"FLAT" contains sequence alignment(s)"
		FLAT"\n\t\tor a distance matrix(ces).\n");
		
	printf (BOLD"\n\t-u "LINE"input user tree file"BOLD", --user_tree="LINE"input user tree file"
		FLAT"\n\t\t"BOLD"%s "FLAT"may use an existing topology available in the "LINE"input user tree file"
		FLAT"\n\t\twhich corresponds to the input dataset.\n", PACKAGE_NAME);

	printf (BOLD"\n\t-o "LINE"output tree file"BOLD", --output_tree="LINE"output tree file"
		FLAT"\n\t\t"BOLD"%s "FLAT"will write the infered tree into the "LINE"output tree file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-O "LINE"output matrix file"BOLD", --output_matrix="LINE"output matrix file"
		FLAT"\n\t\tUse this option if you want "BOLD"%s "FLAT"to write the distances"
		FLAT"\n\t\tmatrix computed from the input alignment in the "LINE"output matrix file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-I "LINE"output information file"BOLD", --output_info="LINE"output information file"
		FLAT"\n\t\tUse this option if you want "BOLD"%s "FLAT"to write information"
		FLAT"\n\t\tabout its execution in the "LINE"output information file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-B "LINE"output bootstrap trees file"BOLD", --output_boot="LINE"output bootstrap trees file"
		FLAT"\n\t\tUse this option if you want "BOLD"%s "FLAT"to write bootstrap trees"
		FLAT"\n\t\tin the "LINE"bootstrap trees file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-a, --append"
		FLAT"\n\t\tUse this option to append results to existing output files (if any)."
		FLAT"\n\t\tBy default output files will be overwritten.\n");

	printf (BOLD"\n\t-m "LINE"method"BOLD", --method="LINE"method"
		FLAT"\n\t\t"BOLD"%s "FLAT"computes a tree using a distance algorithm."
		FLAT"\n\t\tYou may choose this "LINE"method"FLAT" from:"
		FLAT"\n\t\t"BOLD"TaxAdd_(B)alME"FLAT", "BOLD"TaxAdd_(O)LSME"FLAT", "BOLD"B(I)ONJ"FLAT" (default),"
		FLAT"\n\t\t"BOLD"(N)J"FLAT" or "BOLD"(U)NJ"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-d"LINE"[model]"BOLD", --dna="LINE"[model]"
		FLAT"\n\t\tUse this option if your input data file contains DNA sequences alignment(s)."
		FLAT"\n\t\tYou may also indicate the evolutionary "LINE"[model]"FLAT" which can be choosen from:"
		FLAT"\n\t\t"BOLD"(p)-distance"FLAT", "BOLD"R(Y) symmetric"FLAT", "BOLD"(R)Y"FLAT", "BOLD"(J)C69"FLAT", "BOLD"(K)2P"
		FLAT", "BOLD"F8(1)"FLAT", "BOLD"F8(4)"FLAT" (default), "BOLD"(T)N93"FLAT", "BOLD"(L)ogDet"FLAT".\n");

	printf (BOLD"\n\t-p"LINE"[model]"BOLD", --protein="LINE"[model]"
		FLAT"\n\t\tUse this option if your input data file contains protein sequences alignment(s)."
		FLAT"\n\t\tYou may also indicate the evolutionary "LINE"[model]"FLAT" which can be choosen from:"
		FLAT"\n\t\t"BOLD"(p)-distance"FLAT", "BOLD"(F)81 like"FLAT", "BOLD"(L)G"FLAT" (default), "BOLD"(W)AG"FLAT", "BOLD"(J)TT"FLAT", "BOLD"Day(h)off"FLAT", "
		FLAT"\n\t\t"BOLD"(D)CMut"FLAT", "BOLD"(C)pRev"FLAT", "BOLD"(M)tREV"FLAT", "BOLD"(R)tREV"FLAT", "BOLD"HIV(b)"FLAT", "BOLD"H(I)Vw"FLAT" or "BOLD"FL(U)"FLAT".\n");

	printf (BOLD"\n\t-r, --remove_gap"
		FLAT"\n\t\tUse this option to completely remove any site which has a gap in"
		FLAT"\n\t\tany sequence. By default, "BOLD"%s "FLAT"is doing pairwise deletion of gaps.\n", PACKAGE_NAME);

	printf (BOLD"\n\t-e, --equilibrium"
		FLAT"\n\t\tThe equilibrium frequencies for DNA are always estimated by counting"
		FLAT"\n\t\tthe occurence of the nucleotides in the input alignment."
		FLAT"\n\t\tFor amino-acid sequences, the equilibrium frequencies are estimated"
		FLAT"\n\t\tusing the frequencies defined by the substitution model."
		FLAT"\n\t\tUse this option if you whish to estimate the amino-acid frequencies"
		FLAT"\n\t\tby counting their occurence in the input alignment.\n");

	printf (BOLD"\n\t-g"LINE"[alpha]"BOLD", --gamma="LINE"[alpha]"
		FLAT"\n\t\tUse this option if you wish to have gamma distributed rates across sites."
		FLAT"\n\t\tBy default, %s runs with no gamma variation."
		FLAT"\n\t\tIf running %s with gamma distributed rates across sites, the "LINE"[alpha]"FLAT" default value is 1.0."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n", PACKAGE_NAME, PACKAGE_NAME);

	printf (BOLD"\n\t-n"LINE"[NNI]"BOLD", --nni="LINE"[NNI]"
		FLAT"\n\t\tUse this option to do "LINE"[NNI]"FLAT" tree topology improvement."
		FLAT"\n\t\tYou may choose the "LINE"[NNI]"FLAT" type from:"
		FLAT"\n\t\t"BOLD"NNI_(B)alME"FLAT" (default) or "BOLD"NNI_(O)LS"FLAT".\n");

	printf (BOLD"\n\t-s, --spr"
		FLAT"\n\t\tUse this option to do "LINE"SPR"FLAT" tree topology improvement.\n");

	printf (BOLD"\n\t-w "LINE"branch"BOLD", --branch_length="LINE"branch"
		FLAT"\n\t\tUse this option to indicate the "LINE"branch"FLAT" length to assign to the tree."
		FLAT"\n\t\tYou may choose the "LINE"branch"FLAT" length from: "
		BOLD"(B)alLS"FLAT" (default), "BOLD"(O)LS"
		FLAT"\n\t\tor "BOLD"(n)one"FLAT". "BOLD"(n)one "FLAT"is only available with BIONJ, NJ or UNJ."
		FLAT"\n\t\tOnly helpful when not improving the tree topology (no NNI nor SPR).\n");

	printf (BOLD"\n\t-D "LINE"datasets"BOLD", --datasets="LINE"datasets"
		FLAT"\n\t\tUse this option to indicate the number of "LINE"datasets"FLAT" in your input"
		FLAT"\n\t\tdata file. Default value is 1.\n");

	printf (BOLD"\n\t-b "LINE"replicates"BOLD", --bootstrap="LINE"replicates"
		FLAT"\n\t\tUse this option to indicate the number of "LINE"replicates"BOLD" %s "FLAT"will"
		FLAT"\n\t\tdo for bootstrapping. Default value is 0."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n", PACKAGE_NAME);

	printf (BOLD"\n\t-z "LINE"seed"BOLD", --seed="LINE"seed"
		FLAT"\n\t\tUse this option to initialize randomization with "LINE"seed"FLAT" value."
		FLAT"\n\t\tOnly helpful when bootstrapping.\n");

	printf (BOLD"\n\t-c"
		FLAT"\n\t\tUse this option if you want %s only to compute distance matrix."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n", PACKAGE_NAME);

	printf (BOLD"\n\t-f "LINE"number of digits"
		FLAT"\n\t\tUse this option to set the number of digits after the dot to use on output."
		FLAT"\n\t\tDefault precision is 12.\n");

#ifdef _OPENMP
	printf (BOLD"\n\t-T "LINE"number of threads"BOLD", --nb_threads="LINE"number of threads"
		FLAT"\n\t\tUse this option to set the number of threads to use."
		FLAT"\n\t\tDefault "LINE"number of threads"FLAT" is %d.\n", nbthreads);
#endif

	printf (BOLD"\n\t-v "LINE"value"BOLD", --verbose="LINE"value"
		FLAT"\n\t\tSets the verbose level to "LINE"value"FLAT" [0-3]."
		FLAT"\n\t\tDefault "LINE"value"FLAT" is 0.\n");

	printf (BOLD"\n\t-V, --version"
		FLAT"\n\t\tPrints the %s version.\n", PACKAGE_NAME);

	printf (BOLD"\n\t-h, --help"
		FLAT"\n\t\tDisplay this usage.\n\n");

/*
	printf (BOLD"\n\t-t, --TBR"
		FLAT"\n\t\tUse this option to do "LINE"TBR"FLAT" postprocessing.\n");
*/
	
	return;
}

/*********************************************************/

void Set_Defaults_Input (Options *input)
{
	input->I_data_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->I_tree_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_tree_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_mat_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_stat_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_boot_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->open_mode	 = (char *) mCalloc (3, sizeof(char));

	input->fpI_data_file	= NULL;
	input->fpI_tree_file	= NULL;
	input->fpO_tree_file	= NULL;
	input->fpO_mat_file		= NULL;
	input->fpO_stat_file	= NULL;
	input->fpO_boot_file	= NULL;

	strncpy (input->open_mode, "w", 3);
	input->use_O_mat_file	= FALSE;
	input->is_interleaved	= TRUE;
	input->nb_datasets		= 1;
	input->nb_bootstraps	= 0;
	input->input_type		= MATRIX;
	input->method			= BIONJ;
	input->model			= NONE;
	input->global_aa_fq		= TRUE;
	input->use_gamma		= FALSE;
	input->gamma			= 1.0;
	input->only_mat			= FALSE;
	input->precision		= 12;
	input->use_NNI			= FALSE;
	input->NNI				= BALNNI;
	input->branch			= NONE;
	input->seed				= time (NULL);
	input->no_gap			= FALSE;
	input->use_SPR			= FALSE;
//	input->use_TBR			= FALSE;
	verbose					= 0;
	
#ifdef _OPENMP
	input->nb_threads		= omp_get_max_threads();
#else
	input->nb_threads		= 1;
#endif

	return;
}

/*********************************************************/

void Get_Input_CommandLine (Options *input, int argc, char **argv)
{
	int c;

	// getopt_long stores the option index here.
	int option_index = 0;
	
	static struct option long_options[] =
	{
		// These options don't set a flag.
		// We distinguish them by their indices.
		{"input_data",     required_argument, 0, 'i'},	// input file (sequence alignment or distance matrix)
		{"user_tree",      required_argument, 0, 'u'},	// input topology file
		{"output_tree",    required_argument, 0, 'o'},	// output file for the resulting tree
		{"output_matrix",  required_argument, 0, 'O'},	// output file for the distance matrix
		{"output_info",    required_argument, 0, 'I'},	// output file for execution informations
		{"output_boot",    required_argument, 0, 'B'},	// output file for bootstrap trees
		{"datasets",       required_argument, 0, 'D'},	// number of datasets (alignments or matrices)
		{"bootstrap",      required_argument, 0, 'b'},	// number of replicates when bootstrapping
		{"method",         required_argument, 0, 'm'},	// method for building initial tree
		{"nni",            optional_argument, 0, 'n'},	// type of NNI
		{"branch_length",  required_argument, 0, 'w'},	// branch lengths to assign to a topology
		{"seed",           required_argument, 0, 'z'},	// seed for randomization
		{"append",         no_argument,       0, 'a'},	// append results to output files
		{"spr",            no_argument,       0, 's'},	// SPR postprocessing
		//{"tbr",            no_argument,       0, 't'},	// TBR postprocessing
		{"verbose",        required_argument, 0, 'v'},
		{"version",        no_argument,       0, 'V'},
		{"help",           no_argument,       0, 'h'},
		{"compute_mat",    no_argument,       0, 'c'},	// only compute the distance matrix
		{"precision",      no_argument,       0, 'f'},	// output number of digits after dot precision
		{"dna",            optional_argument, 0, 'd'},	// evolutionary model for DNA sequence input (sets input data type to DNA)
		{"protein",        optional_argument, 0, 'p'},	// evolutionary model for PROTEIN sequence input (sets input data type to PROTEIN)
		{"gamma",          optional_argument, 0, 'g'},	// gamma variation across sites
		{"remove_gap",     no_argument,       0, 'r'},	// remove any site which has a gap in any sequence
		{"equilibrium",    no_argument,       0, 'e'},	// estimate amino-acids equilibrium frequencies by counting those in the input alignment
#ifdef _OPENMP
		{"nb_threads",     required_argument, 0, 'T'},	// number of threads
#endif

		{0, 0, 0, 0}
	};
		
	while (1)
	{
#ifdef _OPENMP
		c = getopt_long (argc, argv, "i:u:o:O:I:B:D:b:m:n::w:z:asv:Vhcf:d::p::g::reT:", long_options, &option_index);
#else
		c = getopt_long (argc, argv, "i:u:o:O:I:B:D:b:m:n::w:z:asv:Vhcf:d::p::g::re", long_options, &option_index);
#endif

		// Detect the end of the options.
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);

				printf ("\n");

				break;

			case 'a':
				// append results to output files
				strncpy (input->open_mode, "a", 3);

				break;

			case 'B':
				// output bootstrap trees file
				if (NULL != optarg)
					strncpy (input->O_boot_file, optarg, MAX_FILE_NAME_LENGTH);

				break;

			case 'c':
				input->only_mat = TRUE;
				input->use_O_mat_file = TRUE;

				break;

			case 'b':
				// number of replicates when bootstrapping
				input->nb_bootstraps = atoi (optarg);
				if (input->nb_bootstraps < 1)
					Exit ( (char*)"-b option: '%d' invalid value for number of bootstrap replicates.", input->nb_bootstraps);

				break;

			case 'D':
				// number of datasets (alignments or matrices)
				input->nb_datasets = atoi (optarg);
				if (input->nb_datasets < 1)
					Exit ( (char*)"-D option: '%d' invalid value for number of datasets.", input->nb_datasets);

				break;

			case 'd':
				// evolutionary model for DNA sequence input (sets input data type to DNA)
				input->input_type = DNA;
				if (NULL != optarg)
				{
					if (testD (optarg))
						input->model = getModel_DNA (optarg);

					else
						Exit ( (char*)"-d option: '%s' invalid evolutionary model.", optarg);
				}
				else
					input->model = F84;
				
				break;

			case 'e':
				// local equilibrium frequencies computation
				input->global_aa_fq = FALSE;

				break;
			
			case 'f':
				// number of digits after dot for output precision
				input->precision = atoi (optarg);
				if (input->precision < 1 || input->precision > DECIMAL_DIG -3)
					Exit ( (char*)"-f option: '%d' invalid value for number of digits.", input->precision);

				break;

			case 'g':
				// gamma variation across sites
				input->use_gamma = TRUE;
				if (optarg)
				{
					input->gamma = (float) atof(optarg);
					if (input->gamma < 0.)
						Exit ( (char*)"-g option: '%s' invalid value for gamma rate variation parameter (alpha).", optarg);
				}
				else
					input->gamma = 1.0;
				break;

			case 'h':
#ifdef _OPENMP
				Usage (input->nb_threads);
#else
				Usage ();
#endif
				exit (EXIT_SUCCESS);

				break;

			case 'i':
				// input file (sequence alignment or distance matrix)
				if (Filexists (optarg))
					strncpy (input->I_data_file, optarg, MAX_FILE_NAME_LENGTH);

				else
					Exit ( (char*)"-i option: '%s' file does not exist.", optarg);

				break;

			case 'I':
				// output file for execution informations
				if (NULL != optarg)
					strncpy (input->O_stat_file, optarg, MAX_FILE_NAME_LENGTH);

				break;

			case 'm':
				// method for building initial tree
				if (testM (optarg))
					input->method = getM (optarg);

				else
					Exit ( (char*)"-m option: '%s' invalid method.", optarg);

				break;

			case 'n':
				// NNI postprocessing
				input->use_NNI = TRUE;
				if (optarg)
				{
					// type of NNI
					if (testN (optarg))
						input->NNI = getN (optarg);

					else
						Exit ( (char*)"-n option: '%s' invalid NNI type.", optarg);
				}
				else
				{
					if (input->method == TaxAddOLS || input->method == UNJ)
						input->NNI = OLSNNI;
					else
						input->NNI = BALNNI;
				}

				break;

			case 'o':
				// output file for the resulting tree
				if (NULL != optarg)
					strncpy (input->O_tree_file, optarg, MAX_FILE_NAME_LENGTH);

				break;

			case 'O':
				// output file for the distance matrix
				if (NULL != optarg)
					strncpy (input->O_mat_file, optarg, MAX_FILE_NAME_LENGTH);

				input->use_O_mat_file = TRUE;

				break;

			case 'p':
				input->input_type = PROTEIN;
				// evolutionary model for PROTEIN sequence input (sets input data type to PROTEIN)
				if (NULL != optarg)
				{
					if (testP (optarg))
						input->model = getModel_PROTEIN (optarg);

					else
						Exit ( (char*)"-p option: '%s' invalid evolutionary model.", optarg);
				}
				else
					input->model = LG;

				break;

			case 'r':
				// remove any site which has a gap in any sequence
				input->no_gap = TRUE;

				break;

			case 's':
				// SPR postprocessing
				input->use_SPR = TRUE;

				break;

//			case 't':
				// TBR postprocessing
//				input->use_TBR = TRUE;

//				break;
	
#ifdef _OPENMP
			case 'T':
				// number of threads
				input->nb_threads = (atoi (optarg) > omp_get_max_threads()) ? omp_get_max_threads() : atoi (optarg);
				if (input->nb_threads < 1)
					input->nb_threads = 1;
				break;
#endif

			case 'u':
				// input user tree file
				strncpy (input->I_tree_file, optarg, MAX_FILE_NAME_LENGTH);
				input->method = USER;

				break;

			case 'v':
				// verbose
				if (NULL != optarg)
				{
					verbose = atoi (optarg);
					if (verbose > 3)
						verbose = 3;

					if (verbose < 0)
						verbose = 0;
				}
				else
					verbose = 0;

				break;

			case 'V':
				printf ("%s\n", PACKAGE_STRING);
				exit (EXIT_SUCCESS);

				//break;

			case 'w':
				// branch lengths to assign to a topology
				if (testW (optarg, TRUE))
					input->branch = getW (optarg);

				else
					Exit ( (char*)"-w option: '%s' invalid branch length type.", optarg);

				break;

			case 'z':
				// seed for randomization
				input->seed = atoi (optarg);

				break;

			case '?':
				// getopt_long already printed an error message.

				break;

			default:
				abort ();
		}
	}

	return;
}

/*********************************************************/

void Get_Input_Interactive (Options *input)
{
	int n_trial = 0;
	char choix;
	char *tmp, *c;

	printf ("Enter your input data file name > "); fflush (NULL);
	Getstring_Stdin (input->I_data_file);

	while (! Filexists (input->I_data_file))
	{
		if (++n_trial > 10)
			Exit ( (char*)"The file '%s' does not exist.", basename (input->I_data_file));

		printf ("The file '%s' does not exist.\n", basename (input->I_data_file));
		printf ("Enter your input data file name > "); fflush (NULL);
		Getstring_Stdin (input->I_data_file);
	}

	strncpy (input->O_tree_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
	strncat (input->O_tree_file, "_fastme_tree.txt", 16);
	strncpy (input->O_stat_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
	strncat (input->O_stat_file, "_fastme_stat.txt", 16);
	strncpy (input->O_boot_file, input->I_data_file, MAX_FILE_NAME_LENGTH - 17);
	strncat (input->O_boot_file, "_fastme_boot.txt", 16);

	choix = 0;
	do
	{
		printf ("\n - %s %s - \n\n\n", PACKAGE_NAME, PACKAGE_VERSION);
		printf ("Settings for this run:\n\n");

		printf ("  I "
			"         Input data type (distance matrix or sequence alignment) "
			" %-15s \n", (input->input_type == MATRIX ? "Distance matrix" : (input->input_type == DNA ? "DNA sequence alignment" : "Protein sequence alignment")));

		if (input->input_type != MATRIX) {
			tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
			if (input->input_type == DNA) {
				if (input->model == NONE)
					input->model = F84;
				constantToStr (input->model, tmp);
				printf ("  E "
				"                                          DNA evolutionary model\n"
				"                                (P-dist, RY symetric, RY, JC69, K2P,\n"
				"                                             F81, F84, TN93, LogDet) "
				" %-15s \n", tmp);
			}
			else if (input->input_type == PROTEIN) {
				if (input->model == NONE)
					input->model = LG;
				constantToStr (input->model, tmp);
				printf ("  E "
				"                                      Protein evolutionary model\n"
				"                           (P-dist, F81-like, LG, WAG, JTT, Dayhoff,\n"
				"                        DCMut, CpREV, MtREV, RtREV, HIVb, HIVw, FLU) "
				" %-15s \n", tmp);
			}
			free (tmp);

			//if (input->model != PDIST && input->model != LOGDET && input->model != SCOREDIST) {
			if (input->model != PDIST && input->model != LOGDET) {
				printf ("  G "
					"                            Gamma distributed rates across sites "
					" %-15s \n", (input->use_gamma  ? "yes" : "no" ));

				if (input->use_gamma) {
					printf ("  A "
						"                          Gamma rate variation parameter (alpha) "
						" %-15f \n", input->gamma);
				}
			}
			
			printf ("  R "
				"                                         Remove sites whith gaps "
				" %-15s \n", (input->no_gap ? "yes" : "no"));

			printf ("  O "
				"                               Output calculated distance matrix "
				" %-15s \n", (input->use_O_mat_file ? "yes" : "no"));

			printf ("\n");
		}

		printf ("  D "
			"                                              Number of datasets "
			" %-15d \n", input->nb_datasets);

		tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
		constantToStr (input->method, tmp);
		printf ("  M "
			"                                      Initial tree: build method\n"
			"         (TaxAdd_BalME, TaxAdd_OLSME, BIONJ, NJ or UNJ) or user tree "
			" %-15s \n", tmp);
		free (tmp);

		printf ("  N "
			"                                              NNI postprocessing "
			" %-15s \n", (input->use_NNI ? "yes" : "no"));
		
		printf ("  S "
			"                                              SPR postprocessing "
			" %-15s \n", (input->use_SPR ? "yes" : "no"));

//		printf ("  T "
//			"                                              TBR postprocessing "
//			" %-15s \n", (input->use_TBR ? "yes" : "no"));
	
		if (input->input_type != MATRIX) {
			printf ("\n");
			printf ("  B "
				"                                 Bootstrap: number of replicates "
				" %-15d \n", input->nb_bootstraps);
		}

		printf ("\n");
		printf ("\nAre these settings correct? "
			"(type  Y  or letter for one to change)  ");

		choix = 'X';
		if (scanf ("%c", &choix))
		{
			if (choix == '\n')
				choix = 'X';

			else
				getchar();
		}
		Uppercase (&choix);

		if (choix == 'Y')
			break;

		switch (choix)
		{

			case 'A' :
			{
				if (input->use_gamma)
				{
					c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof (char));
					askOption ( (char*)"Gamma rate variation parameter (alpha) > ", c);
					n_trial = 0;
					while (atof (c) < 0.)
					{
						if (++n_trial > 10)
							Exit ( (char*)"Invalid value for gamma rate variation parameter.");

						printf ("\nInvalid value for gamma rate variation parameter\n");
						askOption ( (char*)"Enter a new value > ", c);
					}
					input->gamma = (float) atof (c);
					free (c);
				}
				break;
			}

			case 'B':
			{
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof (char));
				askOption ( (char*)"Number of replicates > ", c);
				n_trial = 0;
				while ((!atoi (c)) || (atoi (c) < 0))
				{
					if (++n_trial > 10)
						Exit ( (char*)"Invalid number of replicates choosen.");

					printf ("\nInvalid number of replicates choosen\n");
					askOption ( (char*)"Enter a new value > ", c);
				}
				input->nb_bootstraps = atoi (c);
				free (c);
				break;
			}

			case 'D':
			{
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof (char));
				askOption ( (char*)"Number of datatsets > ", c);
				n_trial = 0;
				while ((!atoi (c)) || (atoi (c) <= 0))
				{
					if (++n_trial > 10)
						Exit ( (char*)"The number of datasets must be a positive integer.");

					printf ("\nThe number of datasets must be a positive integer\n");
					askOption ( (char*)"Enter a new value > ", c);
				}
				input->nb_datasets = atoi (c);
				free (c);
				break;
			}

			case 'E':
			{
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof (char));
				if (input->input_type == DNA)
				{
					askOption ( (char*)"Choose the DNA evolutionary model:\n  - (P)-distance, R(Y) symetric, (R)Y, (J)C69, (K)2P, F8(1), F8(4), (T)N93, (L)ogDet\n> ", c);
					n_trial = 0;
					while (! testD (c))
					{
						if (++n_trial > 10)
							Exit ( (char*)"Invalid evolutionary model choosen.");

						printf ("\nInvalid evolutionary model choosen\n");
						askOption ( (char*)"Enter a new value > ", c);
					}
					input->model = getModel_DNA (c);
				}
				else if (input->input_type == PROTEIN)
				{
					askOption ( (char*)"Choose the protein evolutionary model:\n  - (P)-distance, (F)81-like, (L)G, (W)AG, (J)TT, Day(h)off\n    (D)CMut, (C)pREV, (M)tREV, (R)tREV, HIV(b), H(I)Vw, FL(U)\n> ", c);
					n_trial = 0;
					while ((! testP (c)) && (*c != 'S') && (*c != 's'))
					{
						if (++n_trial > 10)
							Exit ( (char*)"Invalid evolutionary model choosen.");

						printf ("\nInvalid evolutionary model choosen\n");
						askOption ( (char*)"Enter a new value > ", c);
					}
					input->model = getModel_PROTEIN (c);
				}
				else
					input->model = NONE;

				free (c);
				break;
			}

			case 'G':
			{
				input->use_gamma = (boolean) abs (input->use_gamma - 1);
				break;
			}

			case 'I':
			{
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof (char));
				askOption ( (char*)"Choose your input data type: distance (M)atrix, (D)NA alignment, (P)rotein alignment > ", c);
				n_trial = 0;
				while (! testI (c))
				{
					if (++n_trial > 10)
						Exit ( (char*)"Invalid datatype choosen.");

					printf ("\nInvalid datatype choosen\n");
					askOption ( (char*)"Enter a new value > ", c);
				}
                input->input_type = getI (c);
                if (input->input_type == DNA)
                	input->model = F84;

                else if (input->input_type == PROTEIN)
                	input->model = LG;

                else
                	input->model = NONE;

				free (c);
				break;
			}

			case 'M':
			{
                boolean ask_I_tree_file = FALSE;
				c = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof (char));
				askOption ( (char*)"Choose your method: TaxAdd_(B)alME, TaxAdd_(O)LSME, B(I)ONJ, (N)J, (U)NJ or u(s)er > ", c);
				n_trial = 0;
				while (! testM (c))
				{
					if (++n_trial > 10)
						Exit ( (char*)"Invalid method choosen.");

					printf ("\nInvalid method choosen\n");
					askOption ( (char*)"Enter a new value > ", c);
				}

                input->method = getM (c);
				if (input->method == TaxAddBAL || input->method == TaxAddOLS)
					input->branch = input->method;

				else if (input->method == USER)
				{
					if (strlen (input->I_tree_file) > 0)
					{
						printf ("\nStarting tree topology file: '%s'\n", input->I_tree_file);
						askOption ( (char*)"Do you wish to change ? (Y/n) > ", c);
						if ((*c == 'y') || (*c == 'Y') || (*c == '\n'))
							ask_I_tree_file = TRUE;
					}
					else
						ask_I_tree_file = TRUE;

					if (ask_I_tree_file)
					{
						printf ("\nEnter the starting tree topology file name > "); fflush (NULL);
						Getstring_Stdin (input->I_tree_file);
						n_trial = 0;
						while (! Filexists (input->I_tree_file))
						{
							if (++n_trial > 10)
								Exit ( (char*)"The file '%s' does not exist.", basename (input->I_tree_file));

							printf ("The file '%s' does not exist.\n", basename (input->I_tree_file));
							printf ("Enter the starting tree topology file name > "); fflush (NULL);
							Getstring_Stdin (input->I_tree_file);
						}
					}
				}
				else
					input->branch = NONE;

				free (c);
				break;
			}

			case 'N':
			{
				input->use_NNI = (boolean) abs (input->use_NNI - 1);
				break;
			}

			case 'O':
			{
				input->use_O_mat_file = (boolean) abs (input->use_O_mat_file - 1);
				break;
			}

			case 'R':
			{
				input->no_gap = (boolean) abs (input->no_gap - 1);
				break;
			}

			case 'S':
			{
				input->use_SPR = (boolean) abs (input->use_SPR - 1);
				break;
			}

//			case 'T':
//			{
//				input->use_TBR = abs (input->use_TBR - 1);
//				break;
//			}

			default:
				break;
		}

	} while (1);

	return;
}

