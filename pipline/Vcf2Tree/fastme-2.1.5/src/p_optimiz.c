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

#include "p_optimiz.h"


/*********************************************************/

phydbl Opt_Dist_F (phydbl dist, phydbl *F, model *mod, int nbthreads)
{
	// required to avoid warning with mono-thread compilation
	nbthreads = nbthreads +0;
	
	phydbl ax, bx, cx, ret;
	phydbl *optdist;
	
	if (dist < BL_MIN)
		dist = BL_MIN;
	
	ax = BL_MIN;
	bx = dist;
	cx = BL_MAX;

	optdist = (phydbl *) mCalloc (1, sizeof (phydbl));

	if (!isBoostrap)
	{
#ifdef _OPENMP
		int i;
		phydbl **optdists = (phydbl **) mCalloc (nbthreads, sizeof (phydbl *));
		for (i=0; i<nbthreads; i++)
			optdists[i] = (phydbl *) mCalloc (1, sizeof (phydbl));
	
		*optdists[omp_get_thread_num()] = dist;
		Dist_F_Brent (ax, bx, cx, 1.E-10, 1000, optdists[omp_get_thread_num()], F, mod);
		ret = *optdists[omp_get_thread_num()];
		
		for (i=0; i<nbthreads; i++)
			free (optdists[i]);
		free (optdists);
#else
		*optdist = dist;
		Dist_F_Brent (ax, bx, cx, 1.E-10, 1000, optdist, F, mod);
		ret = *optdist;
#endif
	}
	else
	{
		*optdist = dist;
		Dist_F_Brent (ax, bx, cx, 1.E-10, 1000, optdist, F, mod);
		ret = *optdist;
	}

	free (optdist);

	return ret;
}

/*********************************************************/

phydbl Dist_F_Brent (phydbl ax, phydbl bx, phydbl cx, phydbl tol,
	int n_iter_max, phydbl *param, phydbl *F, model *mod)
{
	int iter;
	phydbl a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	phydbl init_lnL, curr_lnL;
	phydbl e=0.0;

	//optimize distance, not likelihood
	phydbl old_param, cur_param;

	d = 0.0;
	a = ((ax < cx) ? ax : cx);
	b = ((ax > cx) ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = -Lk_Dist (F, fabs (bx), mod);
	curr_lnL = init_lnL = -fw;

	old_param = cur_param = fabs (bx);

	for (iter=1; iter<=BRENT_ITMAX; iter++)
	{
		xm = 0.5 * (a + b);

		tol1 = tol * fabs (x) + BRENT_ZEPS;
		tol2 = 2.0 * tol1;

		if ((iter > 1) && fabs (old_param-cur_param) < 1.E-06)
		{
			*param = x;
			curr_lnL = Lk_Dist (F, *param, mod);
			return -curr_lnL;
		}

		if (fabs (e) > tol1)
		{
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0)
				p = -p;

			q = fabs (q);
			etemp = e;
			e = d;
			
			if (fabs (p) >= fabs (0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
			{
				d = BRENT_CGOLD * (e = (x >= xm ? a-x : b-x));
			}
			else
			{
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN (tol1, xm - x);
			}
		}
		else
		{
			d = BRENT_CGOLD * (e = (x >= xm ? a-x : b-x));
		}

		u = (fabs (d) >= tol1 ? x + d : x + SIGN (tol1, d));
		if (u < BL_MIN)
			u = BL_MIN;

		(*param) = fabs (u);
		fu = -Lk_Dist (F, fabs (u), mod);
		curr_lnL = -fu;

		if (fu <= fx)
		{
			if (iter > n_iter_max)
				return -fu;

			if (u >= x)
				a = x;

			else
				b = x;

			SHFT (v, w, x, u)
			SHFT (fv, fw, fx, fu)
		}
		else
		{
			if (u < x)
				a = u;

			else
				b = u;

			if (fu <= fw || ( fabs (w - x) < DBL_EPSILON ) )
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || ( fabs (v - x) < DBL_EPSILON ) || ( fabs (v - w) < DBL_EPSILON) )
			{
				v = u;
				fv = fu;
			}
		}
		old_param = cur_param;
		cur_param = *param;
	}

	Exit ( (char*)"Too many iterations in BRENT.");

	//return (-1);
}



