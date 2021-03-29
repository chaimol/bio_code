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


#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <libgen.h>
#include <assert.h>


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


typedef char boolean;

extern boolean isBoostrap;
extern int verbose;


#ifndef true
#define true 1
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef false
#define false 0
#endif
#ifndef FALSE
#define FALSE 0
#endif

// constant for Newick strings
#ifndef INPUT_SIZE
#define INPUT_SIZE 100
#endif
#ifndef MAX_INPUT_SIZE
#define MAX_INPUT_SIZE 1000000
#endif
#ifndef MAX_NAME_LENGTH
#define	MAX_NAME_LENGTH	64
#endif

// maximum distance matrix size
#ifndef MAXSIZE
#define MAXSIZE 65536
#endif

#ifndef MAX_FILE_NAME_LENGTH
#define MAX_FILE_NAME_LENGTH 256
#endif

#ifndef DNA_DIST_MAX
#define DNA_DIST_MAX 5.00
#endif
#ifndef PROT_DIST_MAX
#define PROT_DIST_MAX 20.00
#endif
#ifndef NONE
#define NONE 0
#endif
#ifndef UP
#define UP 1
#endif
#ifndef DOWN
#define DOWN 2
#endif
#ifndef LEFT
#define LEFT 3
#endif
#ifndef RIGHT
#define RIGHT 4
#endif
#ifndef SKEW
#define SKEW 5
#endif

#ifndef ReadOpenParenthesis
#define ReadOpenParenthesis 0
#endif
#ifndef ReadSubTree
#define ReadSubTree 1
#endif
#ifndef ReadLabel
#define ReadLabel 2
#endif
#ifndef ReadWeight
#define ReadWeight 3
#endif
#ifndef AddEdge
#define AddEdge 4
#endif
#ifndef ReadSize
#define ReadSize 5
#endif
#ifndef ReadEntries
#define ReadEntries 6
#endif
#ifndef Done
#define Done 7
#endif

#ifndef TaxAddBAL
#define TaxAddBAL 1
#endif
#ifndef TaxAddOLS
#define TaxAddOLS 2
#endif
#ifndef BALNNI
#define BALNNI 3
#endif
#ifndef OLSNNI
#define OLSNNI 4
#endif
#ifndef NJ
#define NJ 5
#endif
#ifndef UNJ
#define UNJ 6
#endif
#ifndef BIONJ
#define BIONJ 7
#endif
#ifndef BrBAL
#define BrBAL 8
#endif
#ifndef BrOLS
#define BrOLS 9
#endif
#ifndef USER
#define USER 10
#endif
#ifndef NONE
#define NONE 11
#endif

#ifndef DNA_ALPHABET_SIZE
#define DNA_ALPHABET_SIZE 4
#endif
#ifndef DNA_ALPHABET
#define DNA_ALPHABET "ACGT"
#endif

#ifndef PROTEIN_ALPHABET_SIZE
#define PROTEIN_ALPHABET_SIZE 26
#endif
#ifndef PROTEIN_ALPHABET
#define PROTEIN_ALPHABET "ABCDEFGHIKLMNPQRSTVWYZX*?-"
#endif

// constants for DNA alphabet
#ifndef ADENINE
#define ADENINE 0
#endif
#ifndef CYTOSINE
#define CYTOSINE 1
#endif
#ifndef GUANINE
#define GUANINE 2
#endif
#ifndef THYMINE
#define THYMINE 3
#endif

// nucleotides model constants
#ifndef PDIST
#define PDIST 11
#endif
#ifndef RYSYM
#define RYSYM 12
#endif
#ifndef RY
#define RY 13
#endif
#ifndef JC69
#define JC69 14
#endif
#ifndef F81
#define F81 15
#endif
#ifndef F84
#define F84 16
#endif
#ifndef TN93
#define TN93 17
#endif
#ifndef K2P
#define K2P 18
#endif
#ifndef LOGDET
#define LOGDET 20
#endif

// input data type constants
#ifndef MATRIX
#define MATRIX 21
#endif
#ifndef DNA
#define DNA 22
#endif
#ifndef PROTEIN
#define PROTEIN 23
#endif
//#ifndef SCOREDIST
//#define SCOREDIST 24
//#endif

// protein models constants
#ifndef F81LIKE
#define F81LIKE 30
#endif
#ifndef WAG
#define WAG 31
#endif
#ifndef DAYHOFF
#define DAYHOFF 32
#endif
#ifndef JTT
#define JTT 33
#endif
#ifndef BLOSUM62
#define BLOSUM62 34
#endif
#ifndef MTREV
#define MTREV 35
#endif
#ifndef RTREV
#define RTREV 36
#endif
#ifndef CPREV
#define CPREV 37
#endif
#ifndef DCMUT
#define DCMUT 38
#endif
#ifndef VT
#define VT 39
#endif
#ifndef LG
#define LG 40
#endif
#ifndef HIVB
#define HIVB 41
#endif
#ifndef HIVW
#define HIVW 42
#endif    
#ifndef FLU
#define FLU 43
#endif


void constantToStr (int c, char *str);
int *initZeroArray (int l);
int *initOneArray (int l);
double **initDoubleMatrix (int d);
void fillZeroMatrix (double ***A, int d);
boolean whiteSpace (char c);
void Uppercase (char *str);
void Exit (char *message, ...) __attribute__ ((noreturn));
void Warning (char *message, ...);
void Message (char *message, ...);
void Debug (char *message, ...);
void *mCalloc (int nb, size_t size);
boolean isNumeric (const char *p);
char *getLine (FILE *file, char* line, const int len);
char *str_replace (const char *s, char ch, const char *repl);
int countFields (char *str, const char c);
char **str_split (char* str, const char delim);

#endif /*UTILS_H_*/
