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

#ifndef P_MODELS_H
#define P_MODELS_H

#include "p_eigen.h"

void Init_Model (allseq *data, model *mod, boolean global_aa_fq);
int Init_Qmat_Dayhoff (double *daa, phydbl *pi);
int Init_Qmat_JTT (double *daa, phydbl *pi);
int Init_Qmat_MtREV (double *daa, phydbl *pi);
int Init_Qmat_LG (double *daa, phydbl *pi);
int Init_Qmat_WAG (double *daa, phydbl *pi);
int Init_Qmat_DCMut (double *daa, phydbl *pi);
int Init_Qmat_RtREV (double *daa, phydbl *pi);
int Init_Qmat_CpREV (double *daa, phydbl *pi);
int Init_Qmat_VT (double *daa, phydbl *pi);
int Init_Qmat_HIVb (double *daa, phydbl *pi);
int Init_Qmat_HIVw (double *daa, phydbl *pi);
int Init_Qmat_FLU (double *daa, phydbl *pi);
void PMat (phydbl l, model *mod);
void PMat_Zero_Br_Len (model  *mod, double ***Pij);
void PMat_Empirical (phydbl l, model *mod, double ***Pij);

#endif /*P_MODELS_H_*/

