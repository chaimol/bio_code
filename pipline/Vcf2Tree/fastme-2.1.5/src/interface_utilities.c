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


#include "interface_utilities.h"


/*********************************************************/

void PrintOptions (Options *input)
{
	char *tmp = NULL;
	
	Message ( (char*)"Input data file '%s'", input->I_data_file);
	if (strlen (input->I_tree_file) > 0)
	{
		Message ( (char*)"Input tree file '%s'", input->I_tree_file);
	}
	if (strlen (input->O_tree_file) > 0)
	{
		Message ( (char*)"Output tree file '%s'", input->O_tree_file);
	}
	if (input->use_O_mat_file && strlen (input->O_mat_file) > 0)
	{
		Message ( (char*)"Output matrix file '%s'", input->O_mat_file);
	}
	if (strlen (input->O_stat_file) > 0)
	{
		Message ( (char*)"Output stat file '%s'", input->O_stat_file);
	}
	if (strlen (input->O_boot_file) > 0)
	{
		Message ( (char*)"Output bootstrap file '%s'", input->O_boot_file);
	}
	Message ( (char*)"Nb datasets %d", input->nb_datasets);
	Message ( (char*)"Nb bootstraps %d", input->nb_bootstraps);
	switch (input->input_type)
	{
		case MATRIX:
			Message ( (char*)"Input type MATRIX");
			break;
		case DNA:
			Message ( (char*)"Input type DNA");
			break;
		case PROTEIN:
			Message ( (char*)"Input type PROTEIN");
			break;
	}
	if (input->input_type == PROTEIN || input->input_type == DNA)
	{
		if (input->is_interleaved)
		{
			Message ( (char*)"Input data format interleaved");
		}
		else
		{
			Message ( (char*)"Input data format sequential");
		}
		if (input->no_gap)
		{
			Message ( (char*)"Remove sites with gaps");
		}
		else
		{
			Message ( (char*)"Pairwise deletion of gaps");
		}
	}
	if (input->input_type == PROTEIN || input->input_type == DNA)
	{
		switch (input->model)
		{
			case PDIST:
				Message ( (char*)"Model: P distance");
				break;
			case RY:
				Message ( (char*)"Model: RY");
				break;
			case RYSYM:
				Message ( (char*)"Model: RY symetric");
				break;
			case JC69:
				Message ( (char*)"Model: JC69");
				break;
			case F81:
				Message ( (char*)"Model: F81");
				break;
			case F84:
				Message ( (char*)"Model: F84");
				break;
			case K2P:
				Message ( (char*)"Model: K2P");
				break;
			case TN93:
				Message ( (char*)"Model: TN93");
				break;
			case LOGDET:
				Message ( (char*)"Model: LogDet");
				break;
			case F81LIKE:
				Message ( (char*)"Model: F81 like");
				break;
			case WAG:
				Message ( (char*)"Model: WAG");
				break;
			case DAYHOFF:
				Message ( (char*)"Model: Dayhoff");
				break;
			case JTT:
				Message ( (char*)"Model: JTT");
				break;
			case BLOSUM62:
				Message ( (char*)"Model: BLOSUM 62");
				break;
			case MTREV:
				Message ( (char*)"Model: MtRev");
				break;
			case RTREV:
				Message ( (char*)"Model: RtRev");
				break;
			case CPREV:
				Message ( (char*)"Model: CpRev");
				break;
			case DCMUT:
				Message ( (char*)"Model: DCMut");
				break;
			case VT:
				Message ( (char*)"Model: VT");
				break;
			case LG:
				Message ( (char*)"Model: LG");
				break;
			case HIVB:
				Message ( (char*)"Model: HIVb");
				break;
			case HIVW:
				Message ( (char*)"Model: HIVw");
				break;
		}
	}
	
	if (input->input_type == PROTEIN)
	{
		if (input->global_aa_fq)
		{
			Message ( (char*)"Equilibrium frequencies from model");
		}
		else
		{
			Message ( (char*)"Equilibrium frequencies counted from alignment");
		}
	}
	
	if (input->use_gamma)
	{
		Message ( (char*)"Use a gamma law, alpha = %f", input->gamma);
	}
	
	switch (input->method)
	{
		case USER:
			Message ( (char*)"Building tree method: input tree");
			break;
		case TaxAddBAL:
			Message ( (char*)"Building tree method: TaxAdd_BalME");
			break;
		case TaxAddOLS:
			Message ( (char*)"Building tree method: TaxAdd_OLSME");
			break;
		case NJ:
			Message ( (char*)"Building tree method: NJ");
			break;
		case UNJ:
			Message ( (char*)"Building tree method: UNJ");
			break;
		case BIONJ:
			Message ( (char*)"Building tree method: BIONJ");
			break;
	}
	
	if (!input->use_NNI && !input->use_SPR) {
		tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
		constantToStr (input->branch, tmp);
		Message ( (char*)"Assign branch length: %s", tmp);
		free (tmp);
	}
	else if (input->use_NNI && input->use_SPR) 	{
		switch (input->NNI)
		{
			case BALNNI:
				Message ( (char*)"Tree improvement: NNI_BalME & SPR");
				break;
			case OLSNNI:
				Message ( (char*)"Tree improvement: NNI_OLSME & SPR");
				break;
		}
	}
	else if (input->use_NNI) {
		switch (input->NNI)
		{
			case BALNNI:
				Message ( (char*)"Tree improvement: NNI_BalME");
				break;
			case OLSNNI:
				Message ( (char*)"Tree improvement: NNI_OLSME");
				break;
		}
	}
	else if (input->use_SPR) {
		Message ( (char*)"Tree improvement: SPR");
	}
	
//	if (input->use_TBR)
//	{
//		Message("Tree improvement: TBR");
//	}
	
	
	return;
}

/*********************************************************/

void Free_Input (Options *input)
{
	if (NULL != input->I_data_file)
		free (input->I_data_file);
	if (NULL != input->I_tree_file)
		free (input->I_tree_file);/*
	if (NULL != input->I_seqmat_file)
		free (input->I_seqmat_file);*/
	if (NULL != input->O_tree_file)
		free (input->O_tree_file);
	if (NULL != input->O_mat_file)
		free (input->O_mat_file);
	if (NULL != input->O_stat_file)
		free (input->O_stat_file);
	if (NULL != input->O_boot_file)
		free (input->O_boot_file);
	if (NULL != input->open_mode)
		free (input->open_mode);
	if (NULL != input->fpI_data_file)
		fclose (input->fpI_data_file);
	if (NULL != input->fpI_tree_file)
		fclose (input->fpI_tree_file);/*
	if (NULL != input->fpI_seqmat_file)
		fclose (input->fpI_seqmat_file);*/
	if (NULL != input->fpO_tree_file)
		fclose (input->fpO_tree_file);
	if (NULL != input->fpO_mat_file)
		fclose (input->fpO_mat_file);
	if (NULL != input->fpO_stat_file)
		fclose (input->fpO_stat_file);
	if (NULL != input->fpO_boot_file)
		fclose (input->fpO_boot_file);
	if (NULL != input)
		free (input);

	return;
}

/*********************************************************/

int Filexists (char *filename)
{
	FILE *fp;
	fp = fopen (filename, "r");

	if (fp) {
		fclose(fp);
		return 1;
	}
	else
		return 0;
}

/*********************************************************/

void Getstring_Stdin (char *file_name)
{
	if (!(fgets (file_name, MAX_FILE_NAME_LENGTH, stdin)))
		Exit ( (char*)"Cannot read from stdin.");

	if (strchr (file_name, '\n') != NULL)
		*strchr (file_name, '\n') = '\0';

	return;
}

/*********************************************************/

FILE *Openfile (char *filename, char *mode)
{
	FILE *fp;
	fp = NULL;

	fp = fopen(filename, mode);
	if (!fp)
		Exit ( (char*)"Cannot open file '%s'", basename (filename));

	return fp;
}

/*********************************************************/

void askOption (char *question, char *c)
{
	printf ("%s", question);
	Getstring_Stdin (c);

	return;
}

/*********************************************************/

boolean testM (char *c)
{
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 || strncmp (c, "I", 1) == 0 ||
			strncmp (c, "N", 1) == 0 || strncmp (c, "O", 1) == 0 ||
			strncmp (c, "U", 1) == 0 || strncmp (c, "S", 1) == 0)
			ret = TRUE;
		break;
	case 2:
		if (strncmp (c, "NJ", 2) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "BAL", 3) == 0 || strncmp (c, "OLS", 3) == 0 ||
			strncmp (c, "UNJ", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "USER", 4) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "BALME", 5) == 0 || strncmp (c, "OLSME", 5) == 0 ||
			strncmp (c, "BIONJ", 5) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

int getM (char *c)
{
	int ret = TaxAddBAL;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 )
			ret = TaxAddBAL;
		else if (strncmp (c, "I", 1) == 0 )
			ret = BIONJ;
		else if (strncmp (c, "N", 1) == 0 )
			ret = NJ;
		else if (strncmp (c, "O", 1) == 0 )
			ret = TaxAddOLS;
		else if (strncmp (c, "U", 1) == 0 )
			ret = UNJ;
		else if (strncmp (c, "S", 1) == 0)
			ret = USER;
		break;
	case 2:
		if (strncmp (c, "NJ", 2) == 0 )
			ret = NJ;
		break;
	case 3:
		if (strncmp (c, "BAL", 3) == 0 )
			ret = TaxAddBAL;
		else if (strncmp (c, "OLS", 3) == 0 )
			ret = TaxAddOLS;
		else if (strncmp (c, "UNJ", 3) == 0 )
			ret = UNJ;
		break;
	case 4:
		if (strncmp (c, "USER", 4) == 0 )
			ret = USER;
		break;
	case 5:
		if (strncmp (c, "BALME", 5) == 0 )
			ret = TaxAddBAL;
		else if (strncmp (c, "OLSME", 5) == 0 )
			ret = TaxAddOLS;
		else if (strncmp (c, "BIONJ", 5) == 0 )
			ret = BIONJ;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

boolean testN (char *c)
{
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 || strncmp (c, "O", 1) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "BAL", 3) == 0 || strncmp (c, "OLS", 3) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "BALME", 5) == 0 || strncmp (c, "OLSME", 5) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

int getN (char *c)
{
	int ret = BALNNI;

	Uppercase (c);

	if (strncmp (c, "B", 1) == 0 || strncmp (c, "BAL", 3) == 0 ||
		strncmp (c, "BALME", 5) == 0 || strncmp (c, "NNIBALME", 8) == 0 ||
		strncmp (c, "NNI_BALME", 9) == 0)
	{
		ret = BALNNI;
	}
	else if (strncmp (c, "O", 1) == 0 || strncmp (c, "OLS", 3) == 0 ||
		strncmp (c, "OLSME", 5) == 0 || strncmp (c, "NNIOLSME", 8) == 0 ||
		strncmp (c, "NNI_OLSME", 9) == 0)
	{
		ret = OLSNNI;
	}

	return ret;
}

/*********************************************************/

boolean testW (char *c, boolean none)
{
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 || strncmp (c, "O", 1) == 0)
			ret = TRUE;
		else if  (strncmp (c, "N", 1) == 0 && none)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "BAL", 3) == 0 || strncmp (c, "OLS", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "NONE", 4) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "BALLS", 5) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

int getW (char *c)
{
	int ret = BrBAL;

	Uppercase (c);

	if (strncmp (c, "B", 1) == 0 || strncmp (c, "BAL", 3) == 0 ||
		strncmp (c, "BALLS", 5) == 0)
	{
		ret = BrBAL;
	}
	else if (strncmp (c, "O", 1) == 0 || strncmp (c, "OLS", 3) == 0)
	{
		ret = BrOLS;
	}
	else if (strncmp (c, "N", 1) == 0 || strncmp (c, "NONE", 4) == 0)
	{
		ret = NONE;
	}

	return ret;
}

/*********************************************************/

boolean testD (char *c)
{
// (R)Y, F8(1), F8(4), (T)N93, (K)2P, (J)C69, (L)OGDET, (P)DIST, R(Y)SYM
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "R", 1) == 0 || strncmp (c, "1", 1) == 0 ||
			strncmp (c, "4", 1) == 0 || strncmp (c, "T", 1) == 0 ||
			strncmp (c, "K", 1) == 0 || strncmp (c, "J", 1) == 0 ||
			strncmp (c, "L", 1) == 0 || strncmp (c, "P", 1) == 0 ||
			strncmp (c, "Y", 1) == 0)
			ret = TRUE;
		break;
	case 2:
		if (strncmp (c, "RY", 2) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "K2P", 3) == 0 || strncmp (c, "F81", 3) == 0 || strncmp (c, "F84", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "TN93", 4) == 0 || strncmp (c, "JC69", 4) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "PDIST", 5) == 0 || strncmp (c, "RYSYM", 5) == 0)
			ret = TRUE;
		break;
	case 6:
		if (strncmp (c, "LOGDET", 6) == 0 || strncmp (c, "P-DIST", 6) == 0)
			ret = TRUE;
		break;
	case 9:
		if (strncmp (c, "PDISTANCE", 9) == 0)
			ret = TRUE;
		break;
	case 10:
		if (strncmp (c, "P-DISTANCE", 10) == 0)
			ret = TRUE;
		break;
	case 11:
		if (strncmp (c, "RYSYMMETRIC", 11) == 0 )
			ret = TRUE;
		break;
	case 12:
		if (strncmp (c, "RY-SYMMETRIC", 12) == 0 )
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

boolean testP (char *c)
{
//HIV(B), (C)pREV, (D)CMut, (F)81-like, Day(h)off, H(I)Vw, (J)TT, (L)G, (M)tREV, (P)-DIST, (R)tREV, FL(U), (V)T or (W)AG
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "B", 1) == 0 || strncmp (c, "C", 1) == 0 ||
			strncmp (c, "D", 1) == 0 || strncmp (c, "F", 1) == 0 ||
			strncmp (c, "H", 1) == 0 || strncmp (c, "I", 1) == 0 ||
			strncmp (c, "J", 1) == 0 || strncmp (c, "L", 1) == 0 ||
			strncmp (c, "M", 1) == 0 || strncmp (c, "P", 1) == 0 ||
			strncmp (c, "R", 1) == 0 || strncmp (c, "U", 1) == 0 ||
			strncmp (c, "V", 1) == 0 || strncmp (c, "W", 1) == 0)
			ret = TRUE;
		break;
	case 2:
		if (strncmp (c, "LG", 2) == 0 || strncmp (c, "VT", 2) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "JTT", 3) == 0 || strncmp (c, "WAG", 3) == 0 ||
			strncmp (c, "F81", 3) == 0 || strncmp (c, "FLU", 3) == 0)
			ret = TRUE;
		break;
	case 4:
		if (strncmp (c, "HIVB", 4) == 0 || strncmp (c, "HIVW", 4) == 0)
			ret = TRUE;
		break;
	case 5:
		if (strncmp (c, "CPREV", 5) == 0 || strncmp (c, "DCMUT", 5) == 0 ||
			strncmp (c, "MTREV", 5) == 0 || strncmp (c, "RTREV", 5) == 0 ||
			strncmp (c, "PDIST", 5) == 0 )
			ret = TRUE;
		break;
	case 6:
		if (strncmp (c, "P-DIST", 6) == 0)
			ret = TRUE;
		break;
	case 7:
		if (strncmp (c, "DAYHOFF", 7) == 0 || strncmp (c, "F81LIKE", 7) == 0)
			ret = TRUE;
		break;
	case 8:
		if (strncmp (c, "F81-LIKE", 8) == 0)
			ret = TRUE;
		break;
	case 9:
		if (strncmp (c, "PDISTANCE", 9) == 0)
			ret = TRUE;
		break;
	case 10:
		if (strncmp (c, "P-DISTANCE", 10) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}
/*********************************************************/

boolean testI (char *c)
{
//(M)ATRIX, (D)NA, (P)ROTEIN
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "M", 1) == 0 || strncmp (c, "D", 1) == 0 ||
			strncmp (c, "P", 1) == 0)
			ret = TRUE;
		break;
	case 3:
		if (strncmp (c, "DNA", 3) == 0)
			ret = TRUE;
		break;
	case 6:
		if (strncmp (c, "MATRIX", 6) == 0)
			ret = TRUE;
		break;
	case 7:
		if (strncmp (c, "PROTEIN", 7) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

int getI (char *c)
{
	int ret = MATRIX;

	Uppercase (c);

	if (strncmp (c, "M", 1) == 0 || strncmp (c, "MATRIX", 6) == 0)
	{
		ret = MATRIX;
	}
	else if (strncmp (c, "D", 1) == 0 || strncmp (c, "DNA", 3) == 0)
	{
		ret = DNA;
	}
	else if (strncmp (c, "P", 1) == 0 || strncmp (c, "PROTEIN", 7) == 0)
	{
		ret = PROTEIN;
	}

	return ret;
}

/*********************************************************/

boolean testF (char *c)
{
//(I)nterleaved, (S)equential
	boolean ret = FALSE;

	Uppercase (c);

	switch (strlen (c))
	{
	case 1:
		if (strncmp (c, "I", 1) == 0 || strncmp (c, "S", 1) == 0)
			ret = TRUE;
		break;
	case 10:
		if (strncmp (c, "SEQUENTIAL", 10) == 0)
			ret = TRUE;
		break;
	case 11:
		if (strncmp (c, "INTERLEAVED", 11) == 0)
			ret = TRUE;
		break;
	default:
		break;
	}

	return ret;
}

/*********************************************************/

boolean getF (char *c)
{
	boolean ret = TRUE;

	Uppercase (c);

	if (strncmp (c, "I", 1) == 0 || strncmp (c, "INTERLEAVED", 11) == 0)
	{
		ret = TRUE;
	}
	else if (strncmp (c, "S", 1) == 0 || strncmp (c, "SEQUENTIAL", 10) == 0)
	{
		ret = FALSE;
	}

	return ret;
}

/*********************************************************/

int getModel_DNA (char *c)
{
// (R)Y, F8(1), F8(4), (T)N93, (K)2P, (J)C, (L)OGDET, (P)DIST, (S)COREDIST, R(Y)SYM
	int ret = NONE;

	Uppercase (c);

	if (strlen (c) == 1)
	{
		switch (c[0])
		{
			case 'R':
				ret = RY;
				break;
			case '1':
				ret = F81;
				break;
			case '4':
				ret = F84;
				break;
			case 'T':
				ret = TN93;
				break;
			case 'K':
				ret = K2P;
				break;
			case 'J':
				ret = JC69;
				break;
			case 'L':
				ret = LOGDET;
				break;
			case 'P':
				ret = PDIST;
				break;
			case 'Y':
				ret = RYSYM;
				break;
			default:
				break;
		}
	}
	else if (strlen (c) == 2 && strncmp (c, "RY", 2) == 0)
	{
			ret = RY;
	}
	else if (strlen (c) == 3)
	{
		if (strncmp (c, "F81", 3) == 0)
			ret = F81;
	
		else if (strncmp (c, "F84", 3) == 0)
			ret = F84;
	
		else if (strncmp (c, "K2P", 3) == 0)
			ret = K2P;
	}
	else if (strlen (c) == 4)
	{
		if (strncmp (c, "TN93", 4) == 0)
			ret = TN93;
	
		else if (strncmp (c, "JC69", 4) == 0)
			ret = JC69;
	}
	else if (strncmp (c, "RYSYM", 5) == 0 || strncmp (c, "RY-SYM", 6) == 0 || 
			strncmp (c, "RYSYMMETRIC", 10) == 0 || strncmp (c, "RY-SYMMETRIC", 12) == 0)
		ret = RYSYM;
	
	else if (strncmp (c, "LOGDET", 6) == 0)
		ret = LOGDET;
	
	else if (strncmp (c, "PDIST", 5) == 0 || strncmp (c, "P-DIST", 6) == 0 ||
			strncmp (c, "PDISTANCE", 9) == 0 || strncmp (c, "P-DISTANCE", 10) == 0 )
		ret = PDIST;

	
/*	else if (strncmp (c, "S", 1) == 0 || strncmp (c, "SCOREDIST", 9) == 0)
	{
		ret = SCOREDIST;
	}
*/
	return ret;
}

/*********************************************************/

int getModel_PROTEIN (char *c)
{
//HIV(B), (C)pREV, (D)CMut, (F)81-like, Day(h)off, H(I)Vw, (J)TT, (L)G, (M)tREV, (P)-DIST, (R)tREV, FL(U), (V)T or (W)AG

	int ret = NONE;

	Uppercase (c);

	if (strlen (c) == 1)
	{
		switch (c[0])
		{
			case 'B':
				ret = HIVB;
				break;
			case 'C':
				ret = CPREV;
				break;
			case 'D':
				ret = DCMUT;
				break;
			case 'F':
				ret = F81LIKE;
				break;
			case 'H':
				ret = DAYHOFF;
				break;
			case 'I':
				ret = HIVW;
				break;
			case 'J':
				ret = JTT;
				break;
			case 'L':
				ret = LG;
				break;
			case 'M':
				ret = MTREV;
				break;
			case 'P':
				ret = PDIST;
				break;
			case 'R':
				ret = RTREV;
				break;
			case 'U':
				ret = FLU;
				break;
			case 'V':
				ret = VT;
				break;
			case 'W':
				ret = WAG;
				break;
			default:
				break;
		}
	}
	else
	{
		if (strncmp (c, "LG", 2) == 0)
			ret = LG;
	
		else if (strncmp (c, "VT", 2) == 0)
			ret = VT;
	
		else if (strncmp (c, "JTT", 3) == 0)
			ret = JTT;
	
		else if (strncmp (c, "WAG", 3) == 0)
			ret = WAG;
			
		else if (strncmp (c, "FLU", 3) == 0)
			ret = FLU;
			
		else if (strncmp (c, "F81", 3) == 0 || strncmp (c, "F81LIKE", 7) == 0 ||
				strncmp (c, "F81-LIKE", 8) == 0)
			ret = F81LIKE;
	
		else if (strncmp (c, "HIVB", 4) == 0)
			ret = HIVB;
	
		else if (strncmp (c, "HIVW", 4) == 0)
			ret = HIVW;
	
		else if (strncmp (c, "CPREV", 5) == 0)
			ret = CPREV;
	
		else if (strncmp (c, "DCMUT", 5) == 0)
			ret = DCMUT;
	
		else if (strncmp (c, "MTREV", 5) == 0)
			ret = MTREV;
	
		else if (strncmp (c, "RTREV", 5) == 0)
			ret = RTREV;
	
		else if (strncmp (c, "PDIST", 5) == 0 || strncmp (c, "P-DIST", 6) == 0 ||
				strncmp (c, "PDISTANCE", 9) == 0 || strncmp (c, "P-DISTANCE", 10) == 0 )
			ret = PDIST;
	
		else if (strncmp (c, "DAYHOFF", 7) == 0)
			ret = DAYHOFF;
	}
	
	return ret;
}
