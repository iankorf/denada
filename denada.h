/******************************************************************************\
denada.h

Copyright (C) 2003-2005 Ian Korf

\******************************************************************************/

#define dnDURATION_LIMIT 100
#define dnSTATE_LIMIT 100

#include "toolbox.h"
#include "sequence.h"


void dnLogarithm_set (int);


struct dnState {
	char   * name;
	double   init;
	double   term;
	int      transitions;
	ik_tvec  state;
	ik_fvec  score;
	int      order;
	int      count;
	double * emit;
	int      durations;
	double * duration;
};
typedef struct dnState * dnState;
void    dnState_free (dnState);
dnState dnState_read (FILE *);
void    dnState_write (FILE *, const dnState);


struct dnHMM {
	char    * name;
	int       states;
	dnState   state[dnSTATE_LIMIT];
	double    transition[dnSTATE_LIMIT][dnSTATE_LIMIT];
};
typedef struct dnHMM * dnHMM;
void  dnHMM_free (dnHMM);
dnHMM dnHMM_read (FILE *);
void  dnHMM_write (FILE *, const dnHMM);


struct dnMatrix {
	ik_dna     dna;
	dnHMM      hmm;
	double  ** vscore;
	int     ** vtrace;
	int     ** vjumps;
	double  ** fscore;
	double  ** bscore;
    ik_ivec    vpath;
    double     vmax;
    ik_ivec    ppath;
    double     forward;
    double     backward;
};
typedef struct dnMatrix * dnMatrix;
void     dnMatrix_free (dnMatrix);
dnMatrix dnMatrix_new (const dnHMM, const ik_dna);
void     dnMatrix_decode (dnMatrix);
void     dnMatrix_output (const dnMatrix);
void     dnMatrix_browse (const dnMatrix);
void     dnMatrix_posterior_states (const dnMatrix);
void     dnMatrix_posterior_coords (const dnMatrix);
void     dnMatrix_testing (const dnMatrix);

