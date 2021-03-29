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


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Utilities functions for ditance matrix            ;
;                         computation from amino acids sequences            ;
;                                                                           ;
;                         Vincent Lefort                                    ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         vincent.lefort@lirmm.fr                           ;
;                                                                           ;
;                         code is heavily borrowed from PhyML :             ;
;                                                                           ;
;                         Stephane Guindon                                  ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         guindon@lirmm.fr                                  ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

#include "p_utils.h"

/*********************************************************/

double **Copy_PMat_to_DMat (matrix *mat)
{
	int i, j;
	double **D;

	D = initDoubleMatrix (mat->n_otu);
	for (i=0; i<mat->n_otu; i++)
		for (j=0; j<mat->n_otu; j++)
			D[i][j] = (double) mat->dist[i][j];

	return D;
}

/*********************************************************/

void Free_Seq (seq **d, int n_otu)
{
	int i;
	for (i=0; i<n_otu; i++)
	{
		free (d[i]->name);
		free (d[i]->state);
		if ((NULL != d[i]->is_ambigu) && (d[i]->is_ambigu))
			free (d[i]->is_ambigu);

		free (d[i]);
	}

	free(d);

	return;
}

/*********************************************************/

void Free_Model (model *mod)
{
	int i, j;

	if (NULL != mod)
	{
		if (NULL != mod->Pij_rr && mod->n_catg > 0)
		{
			for (i=0; i<mod->n_catg; i++)
			{
				if (NULL != mod->Pij_rr[i])
				{
					for (j=0; j<mod->ns; j++)
					{
						if (NULL != mod->Pij_rr[i][j])
							free (mod->Pij_rr[i][j]);
					}
					free (mod->Pij_rr[i]);
				}
			}
			free (mod->Pij_rr);
		}

		Free_Eigen (mod->eigen);

		free (mod->pi);
		free (mod->gamma_r_proba);
		free (mod->gamma_rr);
		free (mod->qmat);
		free (mod->qmat_buff);
		free (mod->site_num);
		free (mod);
	}
	
	return;
}

/*********************************************************/

void Free_Cseq (allseq *data)
{
	int i;

	free (data->invar);
	free (data->wght);
	free (data->ambigu);
	free (data->b_frq);
	free (data->sitepatt);

	for (i=0; i<data->n_otu; i++)
	{
		free (data->c_seq[i]->name);
		if (data->c_seq[i]->state)
		{
			free (data->c_seq[i]->state);
			if ((NULL != data->c_seq[i]->is_ambigu) && (data->c_seq[i]->is_ambigu))
				free(data->c_seq[i]->is_ambigu);

		}
		free (data->c_seq[i]);
	}

	free (data->c_seq);
	free (data);

	return;
}

/*********************************************************/

void Free_Prefix_Tree (pnode *n, int size)
{
	int i;

	for (i=0; i<size; i++)
	{
		if (n->next[i])
			Free_Prefix_Tree(n->next[i],size);
	}

	Free_Pnode (n);

	return;
}

/*********************************************************/

void Free_Pnode (pnode *n)
{
	free (n->next);
	free (n);

	return;
}

/*********************************************************/

void Free_Eigen (eigen *eigen_struct)
{
	free (eigen_struct->space_int);
	free (eigen_struct->space);
	free (eigen_struct->e_val);
	free (eigen_struct->e_val_im);
	free (eigen_struct->r_e_vect);
	free (eigen_struct->r_e_vect_im);
	free (eigen_struct->l_e_vect);
	free (eigen_struct->q);
	free (eigen_struct);

	return;
}

/*********************************************************/

seq **Get_Seq (FILE *in, boolean interleaved, int *n_otu, int *len, int itype, set *taxa)
{
	seq **data;
	int i, j;
	char **buff;
	int n_unkn, n_removed,pos;
	int *remove;

	if (interleaved)
		data = Read_Seq_Interleaved (in, n_otu, taxa);

	else
		data = Read_Seq_Sequential (in, n_otu, taxa);

	*len = data[0]->len;
	if (verbose > 0 && !isBoostrap)
		Message ( (char*)"Number of taxa: %d. Sequences length: %d.", *n_otu, *len);

	if (data)
	{
		buff = (char **) mCalloc (*n_otu, sizeof (char *));
		for (i=0; i<*n_otu; i++)
			buff[i] = (char *) mCalloc (*len, sizeof (char));

		remove = (int *) mCalloc (*len, sizeof (int));

		n_removed = 0;

		for (i=0; i<*len; i++)
		{
			for (j=0; j<*n_otu; j++)
			{
				if ((data[j]->state[i] == '?') || (data[j]->state[i] == '-'))
					data[j]->state[i] = 'X';

				if ((itype == DNA) && (data[j]->state[i] == 'N'))
					data[j]->state[i] = 'X';

				if (data[j]->state[i] == 'U')
					data[j]->state[i] = 'T';
			}

			n_unkn = 0;
			for (j=0; j<*n_otu; j++)
				if (data[j]->state[i] == 'X')
					n_unkn++;

			if (n_unkn == *n_otu)
			{
				remove[i] = 1;
				n_removed++;
			}

			for (j=0; j<*n_otu; j++)
				buff[j][i] = data[j]->state[i];
		}

		if (n_removed > 0)
		{
			if (itype == DNA)
				printf("\n . %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n", n_removed);

			else
				printf("\n . %d sites are made from completely undetermined states ('X', '-', '?')...\n", n_removed);
		}

		pos = 0;
		for (i=0; i<*len; i++)
		{
			for (j=0; j<*n_otu; j++)
				data[j]->state[pos] = buff[j][i];

			pos++;
		}

		for (i=0; i<*n_otu; i++)
		{
			data[i]->len = pos;
			if (NULL != buff[i])
				free (buff[i]);
		}

		if (NULL != buff)
			free (buff);

		if (NULL != remove)
			free (remove);
	}

	return data;
}

/*********************************************************/

seq **Read_Seq_Sequential (FILE *in, int *n_otu, set *taxa)
{
	int i;
	node *v;
	char *line, *format;;
	int len, readok;
	seq **data = NULL;
	char c;

	line = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
	format = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
	snprintf (format, MAX_NAME_LENGTH, "%%%ds", MAX_NAME_LENGTH);

	readok = len = 0;
	do
	{
		if (fscanf (in, "%s", line) == EOF)
		{
			free (line);
			return NULL;
		}
		else
		{
			if (strncmp (line, "\n", 1) && strncmp (line, "\r", 1) && strncmp (line, "\t", 1))
			{
				*n_otu = atoi (line);
				if (*n_otu <= 0)
					Exit ( (char*)"Problem with sequence format, invalid number of taxa.");

				data = (seq **) mCalloc (*n_otu, sizeof (seq *));
				if (!(fscanf (in, "%s", line)))
					Exit ( (char*)"Problem with sequence format.");

				len = atoi (line);
				if (len <= 0)
					Exit ( (char*)"Problem with sequence format.");
				else
					readok = 1;
			}
		}
	} while (!readok);

	free (line);

	while (((c = (char) fgetc (in)) != '\n') && (c != ' ') && (c != '\r') && (c != '\t'));
	
	for (i=0; i<*n_otu; i++)
	{
		data[i] = (seq *) mCalloc (1, sizeof (seq));
		data[i]->len = 0;
		data[i]->name = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
		data[i]->state = (char *) mCalloc (len + 1, sizeof (char));
		data[i]->is_ambigu = NULL;
		
		// The following code generates warning because 'format' is not a string literal
		if (!(fscanf (in, format, data[i]->name)))
			Exit ( (char*)"Problem with sequence format, invalid sequence name.");

		v = makeNode (data[i]->name, -1);
		v->index2 = i;
		taxa = addToSet (v, taxa);

		while (data[i]->len < len)
			Read_One_Line_Seq (&data, i, in);

		if (data[i]->len != len)
			Exit ( (char*)"Problem with species %s's sequence (check the format).", data[i]->name);
	}

	return data;
}

/*********************************************************/

void Free_Mat (matrix *mat)
{
	int i;

	for (i=0; i<mat->n_otu; i++)
	{
		free (mat->P[i]);
		free (mat->Q[i]);
		free (mat->dist[i]);
		free (mat->name[i]);
	}

	free (mat->P);
	free (mat->Q);
	free (mat->dist);
	free (mat->name);
	free (mat->on_off);
	free (mat);

	return;
}

/*********************************************************/

seq **Read_Seq_Interleaved (FILE *in, int *n_otu, set *taxa)
{
	char c;
	int i, end, num_block, len, readok;
	char *line, *format;
	node *v;
	seq **data = NULL;

	line = (char *) mCalloc (MAX_INPUT_SIZE, sizeof (char));
	format = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
	snprintf (format, MAX_NAME_LENGTH, "%%%ds", MAX_NAME_LENGTH);

	readok = len = 0;
	do
	{
		if (fscanf (in, "%s", line) == EOF)
		{
			free (format);
			free (line);
			return NULL;
		}
		else
		{
			if (strncmp (line, "\n", 1) && strncmp (line, "\r", 1) && strncmp (line, "\t", 1))
			{
				*n_otu = atoi (line);
				if (*n_otu <= 0)
					Exit ( (char*)"Problem with sequence format, invalid number of taxa.");

				data = (seq **) mCalloc (*n_otu, sizeof(seq *));
				if (!(fscanf (in, "%s", line)))
					Exit ( (char*)"Problem with sequence format.");

				len = atoi (line);
				if (len <= 0)
					Exit ( (char*)"Problem with sequence format.");
				else
					readok = 1;
			}
		}
	} while (!readok);

	while (((c = (char) fgetc (in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));

	end = 0;
	for (i=0; i<*n_otu; i++)
	{
		data[i] = (seq *) mCalloc (1, sizeof (seq));
		data[i]->len = 0;
		data[i]->name = (char *) mCalloc (MAX_NAME_LENGTH, sizeof(char));
		data[i]->state = (char *) mCalloc (len + 1, sizeof(char));
		data[i]->is_ambigu = NULL;
		
		// The following code generates warning because 'format' is not a string literal
		if (!(fscanf (in, format, data[i]->name)))
			Exit ( (char*)"Problem with sequence format, invalid sequence name.");
			
		if (checkLabelExist (taxa, data[i]->name))
			Exit ( (char*)"Duplicated sequence name: '%s'", data[i]->name);
			
		v = makeNode (data[i]->name, -1);
		v->index2 = i;
		taxa = addToSet (v, taxa);

		if (!Read_One_Line_Seq (&data, i, in))
		{
			end = 1;

			if ((i != *n_otu) && (i != *n_otu-1))
				Exit ( (char*)"Problem with species %s's sequence.", data[i]->name);

			break;
		}
	}

	if (data[0]->len == len)
		end = 1;

	if (!end)
	{
		end = 0;
		num_block = 1;

		do
		{
			num_block++;

			// interblock
			if (!fgets (line, MAX_INPUT_SIZE, in))
				break;

			if (line[0] != 13 && line[0] != 10)
				Exit ( (char*)"One or more missing sequences in block %d.", num_block - 1);

			for (i=0; i<*n_otu; i++)
				if (data[i]->len != len)
					break;

			if (i == *n_otu)
				break;

			for (i=0; i<*n_otu; i++)
			{
				if (data[i]->len > len)
					Exit ( (char*)"Problem with species %s's sequence. Observed length=%d expected length=%d.",
						data[i]->name, data[i]->len, len);

				else if (!Read_One_Line_Seq (&data, i, in))
				{
					end = 1;
					if ((i != *n_otu) && (i != *n_otu-1))
						Exit ( (char*)"Problem with species %s's sequence.", data[i]->name);

					break;
				}
			}
		} while (!end);
	}

	for (i=0; i<*n_otu; i++)
	{
		if (data[i]->len != len)
			Exit ( (char*)"Check sequence '%s' length.", data[i]->name);
	}

	free (format);
	free (line);

	return data;
}

/*********************************************************/

boolean Read_One_Line_Seq (seq ***data, int num_otu, FILE *in)
{
	boolean ret;
	char c = ' ';
	const char badSymbol[28] = "ABCDEFGHIKLMNOPQRSTUVWXYZ?-.";

	while (1)
	{
		if ((c == EOF) || (c == 13) || (c == 10))
			break;

		else if ((c==' ') || (c=='\t'))
		{
			c = (char) fgetc (in);
			continue;
		}
		Uppercase (&c);

		if (strchr (badSymbol, c) == NULL)
		{
			Exit ( (char*)"Bad symbol: \"%c\" at position %d of species %s.",
				c, (*data)[num_otu]->len, (*data)[num_otu]->name);
		}

		if (c == '.')
		{
			c = (*data)[0]->state[(*data)[num_otu]->len];
			if (!num_otu)
				Exit ( (char*)"Symbol \".\" should not appear in the first sequence.");

		}

		(*data) [num_otu]->state[(*data) [num_otu]->len] = c;
		(*data) [num_otu]->len++;
		c = (char) fgetc (in);
	}

	if (c == EOF)
		ret = FALSE;
	else
		ret = TRUE;

	return ret;
}

/*********************************************************/

model *Make_Model_Basic ()
{
	model *mod;

	mod	= (model *) mCalloc (1, sizeof (model));

	return mod;
}

/*********************************************************/

void Make_Model_Complete (model *mod, Options *options)
{
	int i, j;

	mod->pi				= (phydbl *) mCalloc (mod->ns, sizeof (phydbl));
	mod->gamma_r_proba	= (phydbl *) mCalloc (mod->n_catg, sizeof (phydbl));
	mod->gamma_rr		= (phydbl *) mCalloc (mod->n_catg, sizeof (phydbl));
	mod->Pij_rr			= (double***) mCalloc (mod->n_catg, sizeof (double **));

	for (i=0; i<mod->n_catg; i++)
	{
		mod->Pij_rr[i] = (double **) mCalloc (mod->ns, sizeof (double *));
		for (j=0; j<mod->ns; j++)
			mod->Pij_rr[i][j] = (double *) mCalloc (mod->ns, sizeof (double));
	}

	mod->qmat		= (double *) mCalloc (mod->ns*mod->ns, sizeof (double));
	mod->qmat_buff	= (double *) mCalloc (mod->ns*mod->ns, sizeof (double));
	mod->eigen		= (eigen *) Make_Eigen_Struct (mod->ns);

	mod->use_gamma = options->use_gamma;
	mod->alpha = options->gamma;

	return;
}

/*********************************************************/

model *Copy_Model (model *in)
{
	int i, j, k;
	model *mod;

	mod	= Make_Model_Basic();
	
	mod->whichmodel = in->whichmodel;
	mod->datatype = in->datatype;
	mod->n_catg = in->n_catg;
	mod->invar = in->invar;
	mod->ns = in->ns;
	mod->stepsize = in->stepsize;
	mod->n_otu = in->n_otu;
	mod->bootstrap = in->bootstrap;
	mod->kappa = in->kappa;
	mod->lambda = in->lambda;
	mod->use_gamma = in->use_gamma;
	mod->alpha = in->alpha;
	mod->pinvar = in->pinvar;
	mod->mr = in->mr;
	mod->n_site = in->n_site;
	
	mod->site_num_len = in->site_num_len;
	mod->site_num = (int *) mCalloc (mod->site_num_len, sizeof(int));
	for (i=0; i<mod->site_num_len; i++)
		mod->site_num[i] = in->site_num[i];
	
	mod->pi				= (phydbl *) mCalloc (in->ns, sizeof (phydbl));
	for (i=0; i<in->ns; i++)
	{
		mod->pi[i] = in->pi[i];
	}
		
	mod->gamma_r_proba	= (phydbl *) mCalloc (in->n_catg, sizeof (phydbl));
	mod->gamma_rr		= (phydbl *) mCalloc (in->n_catg, sizeof (phydbl));
	mod->Pij_rr			= (double***) mCalloc (in->n_catg, sizeof (double **));

	for (i=0; i<in->n_catg; i++)
	{
		mod->gamma_r_proba[i] = in->gamma_r_proba[i];
		mod->gamma_rr[i] = in->gamma_rr[i];
		
		mod->Pij_rr[i] = (double **) mCalloc (in->ns, sizeof (double *));
		for (j=0; j<in->ns; j++)
		{
			mod->Pij_rr[i][j] = (double *) mCalloc (in->ns, sizeof (double));
			for (k=0; k<in->ns; k++)
				mod->Pij_rr[i][j][k] = in->Pij_rr[i][j][k];
		}
	}

	mod->qmat		= (double *) mCalloc (in->ns*in->ns, sizeof (double));
	mod->qmat_buff	= (double *) mCalloc (in->ns*in->ns, sizeof (double));
	for (i=0; i<in->ns; i++)
	{
		for (j=0; j<in->ns; j++)
		{
			mod->qmat[i * in->ns + j] = in->qmat[i * in->ns + j];
			mod->qmat_buff[i * in->ns + j] = in->qmat_buff[i * in->ns + j];
		}
	}
	
	mod->eigen = Copy_Eigen_Struct (in->eigen);
	
	return mod;
}

/*********************************************************/

void Set_Defaults_Model (model *mod)
{
	mod->datatype	= PROTEIN;
	mod->whichmodel	= LG;
	mod->n_catg		= 1;
	mod->kappa		= 4.0;
	mod->use_gamma	= FALSE;
	mod->alpha		= 1.0;
	mod->lambda		= 1.0;
	mod->bootstrap	= 0;
	mod->invar		= 0;
	mod->pinvar		= 0.0;
	mod->stepsize	= 1;
	mod->ns			= 20;

	return;
}

/*********************************************************/

eigen *Make_Eigen_Struct (int ns)
{
	eigen *eig;

	eig					= (eigen *) mCalloc (1, sizeof (eigen));
	eig->size			= ns;
	eig->space			= (double *) mCalloc (2 * ns, sizeof (double));
	eig->space_int		= (int *) mCalloc (2 * ns, sizeof (int));
	eig->e_val			= (double *) mCalloc (ns, sizeof (double));
	eig->e_val_im		= (double *) mCalloc (ns, sizeof (double));
	eig->r_e_vect		= (double *) mCalloc (ns * ns, sizeof (double));
	eig->r_e_vect_im	= (double *) mCalloc (ns * ns, sizeof (double));
	eig->l_e_vect		= (double *) mCalloc (ns * ns, sizeof (double));
	eig->q				= (double *) mCalloc (ns * ns, sizeof (double));

	return eig;
}

/*********************************************************/

eigen *Copy_Eigen_Struct (eigen *in)
{
	int i, j, ns;
	
	eigen *out;
	
	ns = in->size;
	
	out = Make_Eigen_Struct (ns);
	out->size = ns;

	for (i=0; i<ns; i++)
	{
		out->e_val[i]		= in->e_val[i];
		out->e_val_im[i]	= in->e_val_im[i];
		out->space[i]		= in->space[i];
		out->space[ns+i]	= in->space[ns+i];
		out->space_int[i]	= in->space_int[i];
		out->space_int[ns+i]= in->space_int[ns+i];
		for (j=0; j<ns; j++)
		{
			out->r_e_vect[i*ns+j]	= in->r_e_vect[i*ns+j];
			out->r_e_vect_im[i*ns+j]= in->r_e_vect_im[i*ns+j];
			out->l_e_vect[i*ns+j]	= in->l_e_vect[i*ns+j];
			out->q[i*ns+j]			= in->q[i*ns+j];
		}
	}

	return out;
}

/*********************************************************/

allseq *Compact_Seq (seq **data, model *mod, boolean rm_ambigu)
{
	allseq *alldata_tmp,*alldata;
	int i, j, k, site;
	int n_patt, which_patt, n_invar;
	char **sp_names;
	int n_otu, n_sites;
	pnode *proot;
	int compress;
	int n_ambigu, is_ambigu;

	n_otu       = mod->n_otu;
	n_patt      = 0;
	which_patt  = 0;

	sp_names = (char **) mCalloc (n_otu, sizeof(char *));

	for (i=0; i<n_otu; i++)
	{
		sp_names[i] = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
		strncpy (sp_names[i], data[i]->name, MAX_NAME_LENGTH);
	}

	alldata_tmp = Make_Cseq (n_otu, data[0]->len, data[0]->len, sp_names);
	proot = (pnode *) Create_Pnode (MAX_NAME_LENGTH);

	for (i=0; i<n_otu; i++)
		free (sp_names[i]);

	free (sp_names);

	if (data[0]->len%mod->stepsize)
		Exit ( (char*)"Sequence length is not a multiple of %d.", mod->stepsize);

	compress = 1;
	n_ambigu = 0;
	is_ambigu = 0;

	for (site=0; site<data[0]->len; site+=mod->stepsize)
	{
		if (rm_ambigu)
		{
			is_ambigu = 0;
			for (j=0; j<n_otu; j++)
				if (Is_Ambigu (data[j]->state+site))
					break;

			if (j != n_otu)
			{
				is_ambigu = 1;
				n_ambigu++;
			}
		}

		if (!is_ambigu)
		{
			if (compress)
			{
				which_patt = -1;
				Traverse_Prefix_Tree (site, -1, &which_patt, &n_patt, data, mod, proot);
				if (which_patt == n_patt-1)	// New pattern found
				{
					n_patt--;
					k = n_patt;
				}
				else
					k = n_patt-10;

			}
			else
			{
				Warning ( (char*)"Sequences are not compressed.");
				k = n_patt;
			}

			if (k == n_patt)	// Add a new site pattern
			{
				for (j=0; j<n_otu; j++)
					Copy_One_State (data[j]->state+site, alldata_tmp->c_seq[j]->state+n_patt, mod->stepsize);

				for (i=0; i<n_otu; i++)
				{
					for (j=0; j<n_otu; j++)
					{
						if (! (Are_Compatible (alldata_tmp->c_seq[i]->state+n_patt, alldata_tmp->c_seq[j]->state+n_patt)))
							break;
					}
					if (j != n_otu)
						break;

				}

				if ((j == n_otu) && (i == n_otu))	// All characters at that site are compatible -> the site is invariant
				{
					for (j=0; j<n_otu; j++)
					{
						alldata_tmp->invar[n_patt] = (short int) Assign_State (alldata_tmp->c_seq[j]->state+n_patt);
						if (alldata_tmp->invar[n_patt] > -1.)
							break;
					}
				}
				else
					alldata_tmp->invar[n_patt] = -1;

				alldata_tmp->sitepatt[site] = n_patt;
				alldata_tmp->wght[n_patt] += 1;
				n_patt += mod->stepsize;
			}
			else
			{
				alldata_tmp->sitepatt[site] = which_patt;
				alldata_tmp->wght[which_patt] += 1;
			}
		}
	}

	data[0]->len -= n_ambigu;
	alldata_tmp->init_len = data[0]->len;
	alldata_tmp->crunch_len = n_patt;

	for (i=0; i<n_otu; i++)
		alldata_tmp->c_seq[i]->len = n_patt;

	if (verbose > 0 && !isBoostrap)
	{
		Message ( (char*)"%d patterns found (out of a total of %d sites).", n_patt, data[0]->len);
		if ((rm_ambigu) && (n_ambigu))
			Message ( (char*)"Removed %d columns of the alignment as they contain ambiguous characters (e.g., gaps).", n_ambigu);
	}

	n_invar = 0;
	for (i=0; i<alldata_tmp->crunch_len; i++)
		if (alldata_tmp->invar[i] > -1.)
			n_invar += (int) alldata_tmp->wght[i];

	if (verbose > 0 && !isBoostrap)
		Message ( (char*)"%d sites without polymorphism (%.2f%c).", n_invar, 100. * (phydbl) n_invar / data[0]->len, '%');

	alldata_tmp->obs_pinvar = (phydbl) n_invar / data[0]->len;

	n_sites = 0;
	for (i=0; i<alldata_tmp->crunch_len; i++)
		n_sites += alldata_tmp->wght[i];

	if (n_sites != data[0]->len)
		Exit ( (char*)"Number of sites does not correspond to data length.");

	Get_AA_Freqs (alldata_tmp);
	alldata = Copy_Cseq (alldata_tmp, alldata_tmp->crunch_len, mod->ns);

	Free_Cseq (alldata_tmp);
	Free_Prefix_Tree (proot, MAX_NAME_LENGTH);

	return alldata;
}

/*********************************************************/

allseq *Compact_CSeq (allseq *data, model *mod)
{
	allseq *alldata;
	int i, j, k, site;
	int n_patt, which_patt;
	int n_otu;

	n_otu = data->n_otu;

	alldata			= (allseq *) mCalloc (1, sizeof (allseq));
	alldata->n_otu	= n_otu;
	alldata->c_seq	= (seq **) mCalloc (n_otu, sizeof (seq *));
	alldata->wght	= (int *) mCalloc (data->crunch_len, sizeof (int));
	alldata->b_frq	= (phydbl *) mCalloc (mod->ns, sizeof (phydbl));
	alldata->ambigu	= (short int *) mCalloc (data->crunch_len, sizeof (short int));
	alldata->invar	= (short int *) mCalloc (data->crunch_len, sizeof (short int));

	alldata->crunch_len = alldata->init_len = -1;
	for (j=0; j<n_otu; j++)
	{
		alldata->c_seq[j]				= (seq *) mCalloc (1, sizeof (seq));
		alldata->c_seq[j]->name			= (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
		strcpy (alldata->c_seq[j]->name, data->c_seq[j]->name);
		alldata->c_seq[j]->state		= (char *) mCalloc (data->crunch_len, sizeof (char));
		alldata->c_seq[j]->is_ambigu	= (short int *) mCalloc (data->crunch_len, sizeof (short int));
		alldata->c_seq[j]->state[0]		= data->c_seq[j]->state[0];
	}

	n_patt = which_patt =  0;

	for (site=0; site<data->crunch_len; site+=mod->stepsize)
	{
		if (data->wght[site])
		{
			for (k=0; k<n_patt; k+=mod->stepsize)
			{
				for (j=0; j<n_otu; j++)
				{
					if (strncmp (alldata->c_seq[j]->state+k, data->c_seq[j]->state+site, (size_t) mod->stepsize))
						break;
				}
				if (j == n_otu)
				{
					which_patt = k;
					break;
				}
			}

			/*    TO DO
			 *  k = n_patt; */

			if (k == n_patt)
			{
				for (j=0; j<n_otu; j++)
					Copy_One_State (data->c_seq[j]->state+site, alldata->c_seq[j]->state+n_patt, mod->stepsize);

				for (i=0; i<n_otu; i++)
				{
					for (j=0; j<n_otu; j++)
					{
						if (!(Are_Compatible(alldata->c_seq[i]->state+n_patt, alldata->c_seq[j]->state+n_patt)))
							break;
					}
					if (j != n_otu)
						break;
				}

				if ((j == n_otu) && (i == n_otu))
				{
					for (j=0; j<n_otu; j++)
					{
						alldata->invar[n_patt] = (short int) Assign_State (alldata->c_seq[j]->state+n_patt);
						if (alldata->invar[n_patt] > -1.)
							break;
					}
				}
				else
					alldata->invar[n_patt] = -1;

				alldata->wght[n_patt] += data->wght[site];
				n_patt+=mod->stepsize;
			}
			else
				alldata->wght[which_patt] += data->wght[site];
		}
	}

	alldata->init_len	= data->crunch_len;
	alldata->crunch_len	= n_patt;
	for (i=0; i<n_otu; i++)
		alldata->c_seq[i]->len = n_patt;

	Get_AA_Freqs (alldata);

	return alldata;
}

/*********************************************************/

allseq *Make_Cseq (int n_otu, int crunch_len, int init_len, char **sp_names)
{
	allseq *alldata;
	int j;

	alldata				= (allseq *) mCalloc (1, sizeof (allseq));
	alldata->n_otu		= n_otu;
	alldata->c_seq		= (seq **) mCalloc (n_otu, sizeof (seq *));
	alldata->b_frq		= (phydbl *) mCalloc (MAX_NAME_LENGTH, sizeof (phydbl));
	alldata->wght		= (int *) mCalloc (crunch_len, sizeof (int));
	alldata->ambigu		= (short int *) mCalloc (crunch_len, sizeof (short int));
	alldata->invar		= (short int *) mCalloc (crunch_len, sizeof (short int));
	alldata->sitepatt	= (int *) mCalloc (init_len, sizeof (int));
	alldata->crunch_len	= crunch_len;
	alldata->init_len	= init_len;
	alldata->obs_pinvar	= .0;

	for (j=0; j<n_otu; j++)
	{
		alldata->c_seq[j]            = (seq *) mCalloc (1, sizeof (seq));
		alldata->c_seq[j]->name      = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
		strncpy(alldata->c_seq[j]->name, sp_names[j], MAX_NAME_LENGTH);
		alldata->c_seq[j]->state     = (char *) mCalloc (crunch_len, sizeof (char));
		alldata->c_seq[j]->is_ambigu = (short int *) mCalloc (crunch_len, sizeof (short int));
	}

	return alldata;
}

/*********************************************************/

void Traverse_Prefix_Tree (int site, int seqnum, int *patt_num, int *n_patt,
	seq **data, model *mod, pnode *n)
{
	int next_state;

	if (seqnum == mod->n_otu - 1)
	{
		n->weight++;
		if (n->weight == 1)
		{
			n->num = *n_patt;
			(*n_patt) += 1;
		}
		(*patt_num) = n->num;
		return;
	}
	else
	{
		next_state = -1;
		next_state = Assign_State_With_Ambiguity (data[seqnum+1]->state + site);
		if (!n->next[next_state])
			n->next[next_state] = Create_Pnode (MAX_NAME_LENGTH);

		Traverse_Prefix_Tree (site, seqnum + 1, patt_num, n_patt, data, mod, n->next[next_state]);
	}

	return;
}

/*********************************************************/

pnode *Create_Pnode (int size)
{
	pnode *n;
	int i;

	n = (pnode *) mCalloc (1, sizeof (pnode));
	n->next = (pnode **) mCalloc (size, sizeof (pnode *));

	for (i=0; i<size; i++)
		n->next[i] = NULL;

	n->weight = 0;
	n->num = -1;

	return n;
}

/*********************************************************/

int Is_Ambigu (char *state)
{
	if (strchr ("X?-.", state[0]))
		return 1;

	else
		return 0;
}

/*********************************************************/

void Copy_One_State (char *from, char *to, int state_size)
{
	int i;

	for (i=0; i<state_size; i++)
		to[i] = from[i];

	return;
}

/*********************************************************/

int Are_Compatible (char *statea, char *stateb)
{
	char a,b;

	a = statea[0];
	b = stateb[0];

	switch (a)
	{
		case 'A' :
		{
			switch (b)
			{
				case 'A' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'R' :
		{
			switch (b)
			{
				case 'R' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'N' :
		{
			switch (b)
			{
				case 'N' :
				case 'B' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'B' :
		{
			switch (b)
			{
				case 'N' :
				case 'B' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'D' :
		{
			switch (b)
			{
				case 'D' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'C' :
		{
			switch (b)
			{
				case 'C' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'Q' :
		{
			switch (b)
			{
				case 'Q' :
				case 'Z' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'Z' :
		{
			switch (b)
			{
				case 'Q' :
				case 'Z' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'E' :
		{
			switch (b)
			{
				case 'E' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'G' :
		{
			switch (b)
			{
				case 'G' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'H' :
		{
			switch (b)
			{
				case 'H' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'I' :
		{
			switch (b)
			{
				case 'I' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'L' :
		{
			switch (b)
			{
				case 'L' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'K' :
		{
			switch (b)
			{
				case 'K' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'M' :
		{
			switch (b)
			{
				case 'M' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'F' :
		{
			switch (b)
			{
				case 'F' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'P' :
		{
			switch (b)
			{
				case 'P' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'S' :
		{
			switch (b)
			{
				case 'S' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'T' :
		{
			switch (b)
			{
				case 'T' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'W' :
		{
			switch (b)
			{
				case 'W' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'Y' :
		{
			switch (b)
			{
				case 'Y' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'V' :
		{
			switch (b)
			{
				case 'V' :
				case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		case 'X' :
		{
			switch (b)
			{
				case 'A' : case 'R' : case 'N' : case 'B' : case 'D' :
				case 'C' : case 'Q' : case 'Z' : case 'E' : case 'G' :
				case 'H' : case 'I' : case 'L' : case 'K' : case 'M' :
				case 'F' : case 'P' : case 'S' : case 'T' : case 'W' :
				case 'Y' : case 'V' : case 'X' :{break;}
				default :
					return 0;
			}
			break;
		}
		default :
		{
			Exit ( (char*)"Please check that characters `%c` and `%c` correspond to existing amino-acids.", a, b);
			//return 0;
		}
	}

	return 1;
}

/*********************************************************/

int Assign_State (char *c)
{
	int state[3];
	state[0] = state[1] = state[2] = -1;
	switch (c[0])
	{
		case 'A' : {state[0]=0 ; break;}
		case 'R' : {state[0]=1 ; break;}
		case 'N' : {state[0]=2 ; break;}
		case 'D' : {state[0]=3 ; break;}
		case 'C' : {state[0]=4 ; break;}
		case 'Q' : {state[0]=5 ; break;}
		case 'E' : {state[0]=6 ; break;}
		case 'G' : {state[0]=7 ; break;}
		case 'H' : {state[0]=8 ; break;}
		case 'I' : {state[0]=9 ; break;}
		case 'L' : {state[0]=10; break;}
		case 'K' : {state[0]=11; break;}
		case 'M' : {state[0]=12; break;}
		case 'F' : {state[0]=13; break;}
		case 'P' : {state[0]=14; break;}
		case 'S' : {state[0]=15; break;}
		case 'T' : {state[0]=16; break;}
		case 'W' : {state[0]=17; break;}
		case 'Y' : {state[0]=18; break;}
		case 'V' : {state[0]=19; break;}
		case 'B' : {state[0]=2 ; break;}
		case 'Z' : {state[0]=5 ; break;}
		default  : {state[0]=-1; break;}
	}

	return state[0];
}

/*********************************************************/

int Assign_State_With_Ambiguity (char *c)
{
	int state[3];
	state[0] = state[1] = state[2] = -1;
	switch (c[0])
	{
		case 'A' : {state[0]= 0; break;}
		case 'R' : {state[0]= 1; break;}
		case 'N' : {state[0]= 2; break;}
		case 'D' : {state[0]= 3; break;}
		case 'C' : {state[0]= 4; break;}
		case 'Q' : {state[0]= 5; break;}
		case 'E' : {state[0]= 6; break;}
		case 'G' : {state[0]= 7; break;}
		case 'H' : {state[0]= 8; break;}
		case 'I' : {state[0]= 9; break;}
		case 'L' : {state[0]=10; break;}
		case 'K' : {state[0]=11; break;}
		case 'M' : {state[0]=12; break;}
		case 'F' : {state[0]=13; break;}
		case 'P' : {state[0]=14; break;}
		case 'S' : {state[0]=15; break;}
		case 'T' : {state[0]=16; break;}
		case 'W' : {state[0]=17; break;}
		case 'Y' : {state[0]=18; break;}
		case 'V' : {state[0]=19; break;}
		case 'B' : {state[0]= 2; break;}
		case 'Z' : {state[0]= 5; break;}
		case 'X' : case '?' : case '-' :
			state[0]=20;
			break;
		default :
		{
			Exit ( (char*)"Unknown character state : %c. Check the data type.", state[0]);
			//break;
		}
	}

	return state[0];
}

/*********************************************************/


void Get_AA_Freqs (allseq *data)
{
	int i, j, w;
	phydbl A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y;
	phydbl fA, fC, fD, fE, fF, fG, fH, fI, fK, fL, fM, fN, fP, fQ, fR, fS, fT, fV, fW, fY;
	phydbl sum;

	fA = fC = fD = fE = fF = fG = fH = fI = fK = fL =
	fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;

	A = C = D = E = F = G = H = I = K = L =
	M = N = P = Q = R = S = T = V = W = Y = .0;

	for (i=0; i<data->n_otu; i++)
	{
		for (j=0; j<data->crunch_len; j++)
		{
			w = data->wght[j];
			if (w)
			{
				switch (data->c_seq[i]->state[j])
				{
					case 'A' : A+=w; break;
					case 'C' : C+=w; break;
					case 'D' : D+=w; break;
					case 'E' : E+=w; break;
					case 'F' : F+=w; break;
					case 'G' : G+=w; break;
					case 'H' : H+=w; break;
					case 'I' : I+=w; break;
					case 'K' : K+=w; break;
					case 'L' : L+=w; break;
					case 'M' : M+=w; break;
					case 'N' : N+=w; break;
					case 'P' : P+=w; break;
					case 'Q' : Q+=w; break;
					case 'R' : R+=w; break;
					case 'S' : S+=w; break;
					case 'T' : T+=w; break;
					case 'V' : V+=w; break;
					case 'W' : W+=w; break;
					case 'Y' : Y+=w; break;
					case 'Z' : Q+=w; break;
					case 'X' : case '?' : case 'O' : case '-' :
						A+=w*fA;
						C+=w*fC;
						D+=w*fD;
						E+=w*fE;
						F+=w*fF;
						G+=w*fG;
						H+=w*fH;
						I+=w*fI;
						K+=w*fK;
						L+=w*fL;
						M+=w*fM;
						N+=w*fN;
						P+=w*fP;
						Q+=w*fQ;
						R+=w*fR;
						S+=w*fS;
						T+=w*fT;
						V+=w*fV;
						W+=w*fW;
						Y+=w*fY;
						break;
					default :
						break;
				}
			}
		}
	}
		
	if (A < 1./20 || C < 1./20 || D < 1./20 || E < 1./20 || F < 1./20 ||
		G < 1./20 || H < 1./20 || I < 1./20 || K < 1./20 || L < 1./20 ||
		M < 1./20 || N < 1./20 || P < 1./20 || Q < 1./20 || R < 1./20 ||
		S < 1./20 || T < 1./20 || V < 1./20 || W < 1./20 || Y < 1./20)
	// if at least one frequency equals 0 then add a pseudo-count
	// as these are doubles, cannot test equality to 0, then test less than minimum value it can have (1./20)
	{
		A+=1.;	C+=1.;	D+=1.;	E+=1.;	F+=1.;
		G+=1.;	H+=1.;	I+=1.;	K+=1.;	L+=1.;
		M+=1.;	N+=1.;	P+=1.;	Q+=1.;	R+=1.;
		S+=1.;	T+=1.;	V+=1.;	W+=1.;	Y+=1.;
	}
	
	sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
	fA = A/sum;		fC = C/sum;		fD = D/sum;		fE = E/sum;
	fF = F/sum;		fG = G/sum;		fH = H/sum;		fI = I/sum;
	fK = K/sum;		fL = L/sum;		fM = M/sum;		fN = N/sum;
	fP = P/sum;		fQ = Q/sum;		fR = R/sum;		fS = S/sum;
	fT = T/sum;		fV = V/sum;		fW = W/sum;		fY = Y/sum;

	data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
	data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
	data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
	data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
	data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV;

	return;
}

/*********************************************************/

allseq *Copy_Cseq (allseq *ori, int len, int ns)
{
	allseq *new;
	int i,j,n_otu;
	char **sp_names;

	n_otu = ori->n_otu;

	sp_names = (char **) mCalloc (n_otu, sizeof (char *));
	for (i=0; i<n_otu; i++)
	{
		sp_names[i] = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
		strncpy (sp_names[i], ori->c_seq[i]->name, MAX_NAME_LENGTH);
	}

	new = Make_Cseq (n_otu, len, ori->init_len, sp_names);
	new->obs_pinvar = ori->obs_pinvar;

	for (i=0; i<ori->init_len; i++)
		new->sitepatt[i] = ori->sitepatt[i];

	for (j=0; j<ori->crunch_len; j++)
	{
		for (i=0; i<ori->n_otu; i++)
		{
			new->c_seq[i]->state[j]     = ori->c_seq[i]->state[j];
			new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
		}
		new->wght[j]   = ori->wght[j];
		new->ambigu[j] = ori->ambigu[j];
		new->invar[j]  = ori->invar[j];
	}

	for (i=0; i<ori->n_otu; i++)
	{
		new->c_seq[i]->len = ori->c_seq[i]->len;
		strncpy (new->c_seq[i]->name, ori->c_seq[i]->name, MAX_NAME_LENGTH);
	}

	new->init_len	= ori->init_len;
	new->clean_len	= ori->clean_len;
	new->crunch_len	= ori->crunch_len;

	for (i=0; i<ns; i++)
		new->b_frq[i] = ori->b_frq[i];

	new->n_otu = ori->n_otu;

	for (i=0; i<n_otu; i++)
		free (sp_names[i]);

	free (sp_names);

	return new;
}

/*********************************************************/

void Check_Ambiguities (allseq *data, int stepsize)
{
	int i, j;

	for (j=0; j<data->crunch_len; j+=stepsize)
	{
		for (i=0; i<data->n_otu; i++)
		{
			data->ambigu[j]              = 0;
			data->c_seq[i]->is_ambigu[j] = 0;
		}
		for (i=0; i<data->n_otu; i++)
		{
			if (Is_Ambigu(data->c_seq[i]->state+j))
			{
				data->ambigu[j]              = 1;
				data->c_seq[i]->is_ambigu[j] = 1;
			}
		}
	}

	return;
}

/*********************************************************/

void Hide_Ambiguities (allseq *data)
{
	int i;

	for (i=0; i<data->crunch_len; i++)
		if (data->ambigu[i])
			data->wght[i] = 0;

	return;
}

/*********************************************************/

boolean Check_2Sequences_Diff (allseq *data)
{
	int i;
	boolean ret = false;

	for (i=0; i<data->crunch_len; i++)
	{
		if ((! data->ambigu[i]) && (data->c_seq[0]->state[i] !=  data->c_seq[1]->state[i]))
		{
			ret = true;
			break;
		}
	}

	return ret;
}

/*********************************************************/

int Matinv (double *x, int n, int m)
{
/* x[n*m]  ... m>=n */
	int i,j,k;
	int *irow;
	double ee, t, t1, xmax;
	double det;

	ee = 1.0E-10;
	det = 1.0;

	irow = (int *) mCalloc (n, sizeof (int));

	for (i=0; i<n; i++)
	{
		xmax = 0.;
		for (j=i; j<n; j++)
			if (xmax < fabs(x[j*m+i]))
			{
				xmax = fabs (x[j*m+i]);
				irow[i]=j;
			}

		det *= xmax;
		if (xmax < ee)
		{
			free (irow);
			Exit ( (char*)"Cannot invert the matrix of eigen vectors : determinant becomes zero at %3d.", i+1);
		}
		if (irow[i] != i)
		{
			for (j=0; j<m; j++)
			{
				t = x[i*m+j];
				x[i*m+j] = x[irow[i] * m + j];
				x[irow[i] * m + j] = t;
			}
		}
		t = 1./x[i * m + i];

		for (j=0; j<n; j++)
		{
			if (j == i)
				continue;

			t1 = t * x[j * m + i];
			for (k=0; k<m; k++)
				x[j * m + k] -= t1 * x[i * m + k];

			x[j * m + i] = -t1;
		}
		for (j=0; j<m; j++)
			x[i * m + j] *= t;

		x[i * m + i] = t;
	}

	for (i=n-1; i>=0; i--)
	{
		if (irow[i] == i)
			continue;

		for (j=0; j<n; j++)
		{
			t = x[j * m + i];
			x[j * m + i] = x[j * m + irow[i]];
			x[j * m + irow[i]] = t;
		}
	}

	free (irow);

	return (1);
}

/*********************************************************/

matrix *JC69_Dist (allseq *data, model *mod)
{
	int site, i, j, k;
	matrix *mat;
	phydbl **len;

	len = (phydbl **) mCalloc (data->n_otu, sizeof (phydbl *));
	for (i=0; i<data->n_otu; i++)
		len[i] = (phydbl *) mCalloc (data->n_otu, sizeof (phydbl));

	mat = Make_Mat (data->n_otu);
	Init_Mat (mat, data);

	for (site=0; site<data->c_seq[0]->len; site+=mod->stepsize)
	{
		for (j=0; j<data->n_otu-1; j++)
		{
			for (k=j+1; k<data->n_otu; k++)
			{
				if ((!Is_Ambigu (data->c_seq[j]->state+site)) &&
					(!Is_Ambigu (data->c_seq[k]->state+site)))
				{
					len[j][k] += data->wght[site];
					len[k][j] = len[j][k];
					if (strncmp (data->c_seq[j]->state + site, data->c_seq[k]->state + site, (size_t) mod->stepsize))
						mat->P[j][k] += data->wght[site];
				}
			}
		}
	}

	for (i=0; i<data->n_otu-1; i++)
		for (j=i+1; j<data->n_otu; j++)
		{
			if (len[i][j])
				mat->P[i][j] /= len[i][j];

			else
				mat->P[i][j] = 1.;

			mat->P[j][i] = mat->P[i][j];

			if ((1. - (mod->ns) / (mod->ns-1.) * mat->P[i][j]) < .0)
				mat->dist[i][j] = PROT_DIST_MAX;

			else
				mat->dist[i][j] = -(mod->ns-1.) / (mod->ns) *
					(phydbl) log (1.-(mod->ns) / (mod->ns-1.) * mat->P[i][j]);

			if (mat->dist[i][j] > PROT_DIST_MAX)
				mat->dist[i][j] = PROT_DIST_MAX;

			mat->dist[j][i] = mat->dist[i][j];
		}

	for (i=0; i<data->n_otu; i++)
		free (len[i]);

	free (len);

	return mat;
}

/*********************************************************/

matrix *Make_Mat (int n_otu)
{
	matrix *mat;
	int i;

	mat = (matrix *) mCalloc (1, sizeof (matrix));

	mat->n_otu	= n_otu;
	mat->P		= (phydbl **) mCalloc (n_otu, sizeof (phydbl *));
	mat->Q		= (phydbl **) mCalloc (n_otu, sizeof (phydbl *));
	mat->dist	= (phydbl **) mCalloc (n_otu, sizeof (phydbl *));
	mat->on_off	= (int *) mCalloc (n_otu, sizeof (int));
	mat->name	= (char **) mCalloc (n_otu, sizeof (char *));

	for (i=0; i<n_otu; i++)
	{
		mat->P[i]		= (phydbl *) mCalloc (n_otu, sizeof (phydbl));
		mat->Q[i]		= (phydbl *) mCalloc (n_otu, sizeof (phydbl));
		mat->dist[i]	= (phydbl *) mCalloc (n_otu, sizeof (phydbl));
		mat->name[i]	= (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));
	}

	return mat;
}

/*********************************************************/

void Init_Mat (matrix *mat, allseq *data)
{
	int i;

	mat->n_otu = data->n_otu;
	mat->r = mat->n_otu;
	mat->curr_int = mat->n_otu;
	mat->method = 1;

	for (i=0; i<data->n_otu; i++)
	{
		strncpy (mat->name[i], data->c_seq[i]->name, MAX_NAME_LENGTH);
		mat->on_off[i] = 1;
	}

	return;
}

/*********************************************************/

void Fill_Missing_Dist (matrix *mat)
{
	int i, j;

	if (verbose > 1 && !isBoostrap)
		Message ( (char*)"Finalizing distance computation...");

	for (i=0; i<mat->n_otu; i++)
	{
		for (j=i+1; j<mat->n_otu; j++)
		{
			if (i != j)
			{
				if (mat->dist[i][j] < .0)
				{
					Fill_Missing_Dist_XY (i, j, mat);
					mat->dist[j][i] = mat->dist[i][j];
				}
			}
		}
	}

	return;
}

/*********************************************************/

void Fill_Missing_Dist_XY (int x, int y, matrix *mat)
{
	int i, j, cpt, pos_best_estimate;
	phydbl *local_mins, **S1S2;
	phydbl min_crit, curr_crit;

	local_mins	= (phydbl *) mCalloc (mat->n_otu * mat->n_otu, sizeof (phydbl ));
	S1S2		= (phydbl **) mCalloc (mat->n_otu * mat->n_otu, sizeof (phydbl *));

	for (i=0; i<mat->n_otu*mat->n_otu; i++)
		S1S2[i] = (phydbl *) mCalloc (2, sizeof (phydbl));

	cpt = 0;
	for (i=0; i<mat->n_otu; i++)
	{
		if ((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
		{
			for (j=0; j<mat->n_otu; j++)
			{
				if ((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
				{
					if ((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
					{
						S1S2[cpt][0] = MIN (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j],
										mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
						S1S2[cpt][1] = MAX (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j],
										mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
						cpt++;
					}
				}
			}
		}
	}

	Qksort_Matrix (S1S2, 0, 0, cpt-1);

	local_mins[0] = S1S2[0][1];
	for (i=1; i<cpt; i++)
		local_mins[i] = (i * local_mins[i-1] + S1S2[i][1]) / (phydbl)(i+1);

	pos_best_estimate = 0;
	min_crit = curr_crit = DBL_MAX;

	for (i=0; i<cpt-1; i++)
	{
		if ((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
		{
			curr_crit = Least_Square_Missing_Dist_XY (x, y, local_mins[i], mat);
			if (curr_crit < min_crit)
			{
				min_crit = curr_crit;
				pos_best_estimate = i;
			}
		}
	}

	mat->dist[x][y] = local_mins[pos_best_estimate];
	mat->dist[y][x] = mat->dist[x][y];

	for (i=0; i<mat->n_otu*mat->n_otu; i++)
		free (S1S2[i]);

	free (S1S2);
	free (local_mins);

	return;
}

/********************************************************/

void Qksort_Matrix (phydbl **A, int col, int ilo, int ihi)
{
	phydbl pivot;		// pivot value for partitioning array
    int ulo, uhi;		// indices at ends of unpartitioned region
    int ieq;			// least index of array entry with value equal to pivot
    phydbl *tempEntry;	// temporary entry used for swapping

    tempEntry = NULL;

	if (ilo >= ihi)
		return;

	// Select a pivot value.
    pivot = A[(ilo + ihi)/2][col];
    /* Initialize ends of unpartitioned region and least index of entry
	 * with value equal to pivot. */
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi)
	{
		if (A[uhi][col] > pivot)	/* Here, we can reduce the size of
									 * the unpartitioned region and try again. */
			uhi--;

		else	// Here, A[uhi] <= pivot, so swap entries at indices ulo and uhi.
		{
			tempEntry = A[ulo];
			A[ulo] = A[uhi];
			A[uhi] = tempEntry;
			// After the swap, A[ulo] <= pivot.
			if (A[ulo][col] < pivot)
			{
				// Swap entries at indices ieq and ulo.
				tempEntry = A[ieq];
				A[ieq] = A[ulo];
				A[ulo] = tempEntry;
				// After the swap, A[ieq] < pivot, so we need to change ieq.
				ieq++;
				/* We also need to change ulo, but we also need to do
				 * that when A[ulo] = pivot, so we do it after this if statement.*/
			}
			/* Once again, we can reduce the size
			 * of the unpartitioned region and try again.*/
			ulo++;
		}
	}
	/* Now, all entries from index ilo to ieq - 1 are less than the pivot
	 * and all entries from index uhi to ihi + 1 are greater than the
	 * pivot. So we have two regions of the array that can be sorted
	 * recursively to put all of the entries in order.*/
	Qksort_Matrix (A, col, ilo, ieq - 1);
	Qksort_Matrix (A, col, uhi + 1, ihi);

	return;
}

/*********************************************************/

phydbl Least_Square_Missing_Dist_XY (int x, int y, phydbl dxy, matrix *mat)
{
	int i, j;
	phydbl fit;

	fit = .0;
	for (i=0; i<mat->n_otu; i++)
	{
		if ((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
		{
			for (j=0; j<mat->n_otu; j++)
			{
				if ((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
				{
					if ((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
					{
						if (dxy < MIN (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j],
								mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
						{
							fit += pow ((mat->dist[i][x] + mat->dist[j][y]) -
									 (mat->dist[i][y] + mat->dist[j][x]), 2);
						}
						else if ((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
						{
							fit += pow (dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]), 2);
						}
						else
						{
							fit += pow (dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]), 2);
						}
					}
				}
			}
		}
	}

	return fit;
}

/*********************************************************/

void Read_Qmat (double *daa, phydbl *pi, FILE *fp)
{
	int i, j;
	phydbl sum;

	for (i=1; i<20; i++)
	{
		for (j=0; j<19; j++)
		{
			if (!(fscanf (fp, "%lf", &(daa[i * 20 + j]))))
				Exit ( (char*)"The rate matrix format is incorrect.");

			daa[j * 20 + i] = daa[i * 20 + j];
			if(j == i-1)
				break;
		}
	}

	for (i=1; i<20; i++)
	{
		if (!(fscanf (fp, "%lf", pi+i)))
			Exit ( (char*)"The rate matrix format is incorrect.");
	}

	sum = .0;
	for (i=1; i<20; i++)
		sum += pi[i];

	if (fabs (sum - 1.) > 1.E-06)
		Exit ( (char*)"The rate matrix format is incorrect.");

	return;
}
