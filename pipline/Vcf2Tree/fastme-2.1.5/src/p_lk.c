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

#include "p_lk.h"

/*********************************************************/

allseq *Init2Data (allseq *data, int ns)
{
	int i;
	allseq *ret;

	ret			= (allseq *) mCalloc (1, sizeof (allseq));
	ret->c_seq	= (seq **) mCalloc (2, sizeof (seq *));
	ret->b_frq	= (phydbl *) mCalloc (ns, sizeof (phydbl));
	ret->ambigu	= (short int *) mCalloc (data->crunch_len, sizeof (short int));

	ret->n_otu = 2;
	ret->crunch_len = 0;
	ret->init_len = 0;

	for (i=0; i<data->crunch_len; i++)
	{
		if (data->wght[i] > .0)
		{
			ret->crunch_len++;
			ret->init_len += (int) data->wght[i];
		}
	}

	return ret;
}

/*********************************************************/

void Free2Data (allseq *data)
{
	free (data->ambigu);
	free (data->b_frq);
	free (data->c_seq);
	free (data);

	return;
}

/*********************************************************/

matrix *ML_Dist (allseq *data, model *mod, int nbthreads)
{
	int i, j, k, l;
	int state0, state1, len;
	phydbl init, sum;
	boolean warn;
	phydbl d_max;
	matrix *mat;
	allseq **twodata, **tmpdata;
	phydbl **Fs;
	
	mat = JC69_Dist (data, mod);

	for (i=0; i<mod->n_catg; i++)
	{
		mod->gamma_rr[i] = 1.0;
		mod->gamma_r_proba[i] = 1.0;
	}

	warn = FALSE;
	
	// Create F for each thread
	Fs = (phydbl **) mCalloc (nbthreads, sizeof (phydbl *));
	for (i=0; i<nbthreads; i++)
		Fs[i] = (phydbl *) mCalloc (mod->ns*mod->ns, sizeof (phydbl));
	
	// Create data for each thread
	twodata = (allseq **) mCalloc (nbthreads, sizeof (allseq *));
	tmpdata = (allseq **) mCalloc (nbthreads, sizeof (allseq *));
	
	
	if (!isBoostrap)	// here FastME computes distances
						// when the process is NOT computing bootstrap pseudo-trees
						// the code is parallelized with OpenMP for the distances computation
	{
#ifdef _OPENMP
		// Create models for each thread
		model **models = (model **) mCalloc (nbthreads-1, sizeof (model *));
		for (i=0; i<nbthreads-1; i++)
			models[i] = Copy_Model (mod);
	#pragma omp parallel for private (i, j, k, l, init, d_max, len, sum, state0, state1) shared (mat, warn)
#endif
		for (j=0; j<data->n_otu; j++)
		{ // begin for j->n_otu
#ifdef _OPENMP
			tmpdata[omp_get_thread_num()] = Init2Data (data, mod->ns);
			tmpdata[omp_get_thread_num()]->c_seq[0] = data->c_seq[j];
			tmpdata[omp_get_thread_num()]->c_seq[0]->name = data->c_seq[j]->name;
			tmpdata[omp_get_thread_num()]->wght = data->wght;
#else
			tmpdata[0] = Init2Data (data, mod->ns);
			tmpdata[0]->c_seq[0] = data->c_seq[j];
			tmpdata[0]->c_seq[0]->name = data->c_seq[j]->name;
			tmpdata[0]->wght = data->wght;
#endif

			for (k=j+1; k<data->n_otu; k++)
			{ // begin for k->n_otu
#ifdef _OPENMP
				tmpdata[omp_get_thread_num()]->c_seq[1] = data->c_seq[k];
				tmpdata[omp_get_thread_num()]->c_seq[1]->name = data->c_seq[k]->name;
				twodata[omp_get_thread_num()] = Compact_CSeq (tmpdata[omp_get_thread_num()], mod);
#else
				tmpdata[0]->c_seq[1] = data->c_seq[k];
				tmpdata[0]->c_seq[1]->name = data->c_seq[k]->name;
				twodata[0] = Compact_CSeq (tmpdata[0], mod);
#endif

				for (l=0; l<mod->ns; l++)
				{
#ifdef _OPENMP
					twodata[omp_get_thread_num()]->b_frq[l] = data->b_frq[l];
#else
					twodata[0]->b_frq[l] = data->b_frq[l];
#endif
				}

#ifdef _OPENMP
				Check_Ambiguities (twodata[omp_get_thread_num()], 1);
#else
				Check_Ambiguities (twodata[0], 1);
#endif

				// If sequences are different (avoid ambiguities), compute distance
				// Else distance = 0
#ifdef _OPENMP
				if (Check_2Sequences_Diff (twodata[omp_get_thread_num()]))
				{
					Hide_Ambiguities (twodata[omp_get_thread_num()]);
#else
				if (Check_2Sequences_Diff (twodata[0]))
				{
					Hide_Ambiguities (twodata[0]);
#endif

					init = mat->dist[j][k];

					if ((init == PROT_DIST_MAX) || (init < .0))
						init = 0.1;

					d_max = init;

#ifdef _OPENMP
					memset (Fs[omp_get_thread_num()], 0, (size_t)(mod->ns*mod->ns) * sizeof (phydbl));
#else
					memset (Fs[0], 0, (size_t)(mod->ns*mod->ns) * sizeof (phydbl));
#endif

					len = 0;
				
#ifdef _OPENMP
					for (l=0; l<twodata[omp_get_thread_num()]->c_seq[0]->len; l++)
					{
						state0 = Assign_State (twodata[omp_get_thread_num()]->c_seq[0]->state + l);
						state1 = Assign_State (twodata[omp_get_thread_num()]->c_seq[1]->state + l);
						if ((state0 > -1) && (state1 > -1))
						{
							Fs[omp_get_thread_num()][mod->ns*state0+state1] += twodata[omp_get_thread_num()]->wght[l];
							len += (int) twodata[omp_get_thread_num()]->wght[l];
						}
					}
#else
					for (l=0; l<twodata[0]->c_seq[0]->len; l++)
					{
						state0 = Assign_State (twodata[0]->c_seq[0]->state + l);
						state1 = Assign_State (twodata[0]->c_seq[1]->state + l);
						if ((state0 > -1) && (state1 > -1))
						{
							Fs[0][mod->ns*state0+state1] += twodata[0]->wght[l];
							len += (int) twodata[0]->wght[l];
						}
					}
#endif

					if (len > .0)
					{
						for (i=0; i<mod->ns*mod->ns; i++)
						{
#ifdef _OPENMP
							Fs[omp_get_thread_num()][i] /= (phydbl) len;
#else
							Fs[0][i] /= (phydbl) len;
#endif
						}
					}

					sum = 0.;
					for (i=0; i<mod->ns*mod->ns; i++)
					{
#ifdef _OPENMP
						sum += Fs[omp_get_thread_num()][i];
#else
						sum += Fs[0][i];
#endif
					}

					if (sum < .001)
						d_max = -1.;

					else if ((sum > 1. - .001) && (sum < 1. + .001))
					{
#ifdef _OPENMP
						if (omp_get_thread_num() == 0)
							d_max = Opt_Dist_F (d_max, Fs[omp_get_thread_num()], mod, nbthreads);
						else
							d_max = Opt_Dist_F (d_max, Fs[omp_get_thread_num()], models[omp_get_thread_num()-1], nbthreads);
#else
						d_max = Opt_Dist_F (d_max, Fs[0], mod, nbthreads);
#endif
					}

					else
						Exit ( (char*)"Invalid value when computing distance. sum = %f.", sum);

					if (d_max >= PROT_DIST_MAX)
					{
						warn = TRUE;
						d_max = PROT_DIST_MAX;
					}
					// Do not correct for dist < BL_MIN,
					// otherwise Fill_Missing_Dist will not be called
				}
				else
					d_max = 0.;

				mat->dist[j][k] = d_max;
				mat->dist[k][j] = mat->dist[j][k];

#ifdef _OPENMP
				Free_Cseq (twodata[omp_get_thread_num()]);
#else
				Free_Cseq (twodata[0]);
#endif

			} // end for k->n_otu

#ifdef _OPENMP
		Free2Data (tmpdata[omp_get_thread_num()]);
#else
		Free2Data (tmpdata[0]);
#endif

		}  // end for j->n_otu

#ifdef _OPENMP		
	for (i=0; i<nbthreads-1; i++)
		free (models[i]);
	free (models);
#endif

	}
	else 	// here FastME computes distances
			// when the process is computing bootstrap pseudo-trees
			// the code is the same for the mono or multi -threaded process
			// the parallelization was already done with OpenMP for the bootstraps computation
	{
		for (j=0; j<data->n_otu; j++)
		{ // begin for j->n_otu
			tmpdata[0] = Init2Data (data, mod->ns);
			tmpdata[0]->c_seq[0] = data->c_seq[j];
			tmpdata[0]->c_seq[0]->name = data->c_seq[j]->name;
			tmpdata[0]->wght = data->wght;
			
			for (k=j+1; k<data->n_otu; k++)
			{ // begin for k->n_otu
				tmpdata[0]->c_seq[1] = data->c_seq[k];
				tmpdata[0]->c_seq[1]->name = data->c_seq[k]->name;
				twodata[0] = Compact_CSeq (tmpdata[0], mod);
				
				for (l=0; l<mod->ns; l++)
				{
					twodata[0]->b_frq[l] = data->b_frq[l];
				}
				
				Check_Ambiguities (twodata[0], 1);
				
				// If sequences are different (avoid ambiguities), compute distance
				// Else distance = 0
				if (Check_2Sequences_Diff (twodata[0]))
				{
					Hide_Ambiguities (twodata[0]);
				
					init = mat->dist[j][k];
					
					if ((init == PROT_DIST_MAX) || (init < .0))
						init = 0.1;
					
					d_max = init;
					
					memset (Fs[0], 0, (size_t)(mod->ns*mod->ns) * sizeof (phydbl));
					
					len = 0;
					
					for (l=0; l<twodata[0]->c_seq[0]->len; l++)
					{
						state0 = Assign_State (twodata[0]->c_seq[0]->state + l);
						state1 = Assign_State (twodata[0]->c_seq[1]->state + l);
						if ((state0 > -1) && (state1 > -1))
						{
							Fs[0][mod->ns*state0+state1] += twodata[0]->wght[l];
							len += (int) twodata[0]->wght[l];
						}
					}
					
					if (len > .0)
					{
						for (i=0; i<mod->ns*mod->ns; i++)
						{
							Fs[0][i] /= (phydbl) len;
						}
					}
					
					sum = 0.;
					for (i=0; i<mod->ns*mod->ns; i++)
					{
						sum += Fs[0][i];
					}
					
					if (sum < .001)
						d_max = -1.;
					
					else if ((sum > 1. - .001) && (sum < 1. + .001))
					{
						d_max = Opt_Dist_F (d_max, Fs[0], mod, nbthreads);
					}
					
					else
						Exit ( (char*)"Invalid value when computing distance. sum = %f.", sum);
					
					if (d_max >= PROT_DIST_MAX)
					{
						warn = TRUE;
						d_max = PROT_DIST_MAX;
					}
					// Do not correct for dist < BL_MIN,
					// otherwise Fill_Missing_Dist will not be called
				}
				else
					d_max = 0.;
				
				mat->dist[j][k] = d_max;
				mat->dist[k][j] = mat->dist[j][k];
				
				Free_Cseq (twodata[0]);
				
			} // end for k->n_otu
			
			Free2Data (tmpdata[0]);	
			
		} // end for j->n_otu
	
	} // end if/else isBootstrap
	
	for (i=0; i<nbthreads; i++)
		free (Fs[i]);
	free (Fs);
	free (tmpdata);
	free (twodata);

	if (warn && !isBoostrap)
		Warning ( (char*)"Give up this dataset because at least one distance exceeds %.2f.", PROT_DIST_MAX);

	return mat;
}

/*********************************************************/

phydbl Lk_Dist (phydbl *F, phydbl dist, model *mod)
{
	int i, j;
	phydbl len, lnL;

	len = -1.;
	
	for (i=0; i<mod->n_catg; i++)
	{
		len = dist * mod->gamma_rr[i];
			
		if (len < BL_MIN)
			len = BL_MIN;

		else if (len > BL_MAX)
			len = BL_MAX;

		PMat (len, mod);
	}
	
	lnL = .0;
	
	for (i=0; i<mod->ns; i++)
	{
		for (j=0; j<mod->ns; j++)
		{
			lnL += F[mod->ns*i+j] * (phydbl)log (partialLK (mod, i, j));
		}
	}

	return lnL;
}

/*********************************************************/

phydbl partialLK (model *mod, int i, int j)
{
	int k;
	phydbl lk = .0;
	
	for (k=0; k<mod->n_catg; k++)
	{
		lk += mod->gamma_r_proba[k] * mod->pi[i] * 
				(phydbl)(mod->Pij_rr[k][i][j]);
	}
	
	return lk;
}


