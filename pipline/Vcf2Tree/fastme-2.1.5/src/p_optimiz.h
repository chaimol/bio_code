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

#ifndef P_OPTIMIZ_H
#define P_OPTIMIZ_H

#include "p_models.h"
#include "p_lk.h"

#ifndef UNLIKELY
#define UNLIKELY -1.e10
//#define UNLIKELY -DBL_EPSILON
#endif

#ifndef BRENT_ZEPS
#define BRENT_ZEPS 1.e-10
//#define BRENT_ZEPS DBL_EPSILON
#endif

#ifndef BRENT_ITMAX
#define BRENT_ITMAX 10000
#endif

#ifndef BRENT_CGOLD
#define BRENT_CGOLD 0.3819660
#endif


phydbl Opt_Dist_F (phydbl dist, phydbl *F, model *mod, int nbthreads);

phydbl Dist_F_Brent (phydbl ax, phydbl bx, phydbl cx, phydbl tol,
	int n_iter_max, phydbl *param, phydbl *F, model *mod);


#endif /*P_OPTIMIZ_H_*/
