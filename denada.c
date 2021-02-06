/******************************************************************************\
denada.c

Copyright (C) 2003-2006 Ian Korf

\******************************************************************************/

#include "denada.h"

static int LOGARITHMS = 1; /* use by default */
void dnLogarithm_set (int i) {
	LOGARITHMS = i;
}

static double LOG2 = 0.6931471806;

static double dn_log2 (double x) {
	if (x <= 0) ik_exit(1, "attempt to take log of %g", x);
	return (log(x) / LOG2);
}

static double p2s (double p) {
	if (!LOGARITHMS) return p;
	if (p == 0) return MIN_SCORE;
	if (p < 0 || p > 1) ik_exit(1, "prob not in range 0..1 (%f)", p);
	return dn_log2(p);
}

static double s2p (double s) {
	if (!LOGARITHMS) return s;
	if (s == MIN_SCORE) return 0;
	return pow(2, s);
}

/*
static double sum_log2p (double v1, double v2) {
	double larger, smaller, p1, p2, sum;
	
	if (v1 > v2) {
		larger = v1;
		smaller = v2;
	} else {
		larger = v2;
		smaller = v1;
	}
	
	if (larger == MIN_SCORE) return MIN_SCORE;
	
	if (larger - smaller > 20) return larger;
	
	p1 = exp2(smaller - v1);
	p2 = exp2(smaller - v2);
	sum = p1 + p2;
	
	return (log(sum)/log(2)) + larger;
}
*/

void dnState_free (dnState state) {
	ik_free(state->name);
}

dnState dnState_read (FILE * stream) {
	int		i;
	char	tag[64], name[64];
	double	value;
	dnState s = ik_malloc(sizeof(struct dnState));
		
	if (fscanf(stream, "%s %s %lf %lf %d %d %d", tag, name, &s->init, &s->term,
		&s->transitions, &s->order, &s->durations) != 7)
		ik_exit (1, "denada-State parse error 1");
	s->init = p2s(s->init);
	s->term = p2s(s->term);
			
	/* name */
	s->name = ik_malloc(strlen(name) +1);
	strcpy(s->name, name);
	
	/* transitions */
	s->state = ik_tvec_new();
	s->score = ik_fvec_new();
	for (i = 0; i < s->transitions; i++) {
		if (fscanf(stream, "%s", name) != 1)
			ik_exit(2, "denada-State parse error 2");
		ik_tvec_push(s->state, name);
		if (fscanf(stream, "%lf", &value) != 1)
			ik_exit(3, "denada-State parse error 3");
		ik_fvec_push(s->score, p2s(value));
	}
	
	/* emissions */
	s->count = pow(4, (double)s->order +1);
	s->emit = ik_malloc(sizeof(double) * s->count);
	for (i = 0; i < s->count; i++) {
		if (fscanf(stream, "%lf", &value) != 1)
			ik_exit(4, "denada-State parse error 4");
		s->emit[i] = p2s(value);
	}
	
	/* durations */
	if (s->durations) s->duration = ik_malloc(sizeof(double) * s->durations);
	for (i = 0; i < s->durations; i++) {
		if (fscanf(stream, "%lf", &value) != 1) 
			ik_exit(5, "denada-State parse error 5");
		s->duration[i] = p2s(value);
	}
	
	/* explicit length states may not start or end */
	if (s->durations) {
		if (LOGARITHMS) {
			if (s->init != MIN_SCORE) ik_exit(1, "explicit states can't init");
			if (s->term != MIN_SCORE) ik_exit(1, "explicit states can't term");
		} else {
			if (s->init != 0) ik_exit(1, "explict states can't init");
			if (s->term != 0) ik_exit(1, "explict states can't term");
		}
	}
	
	return s;
}

void dnState_write (FILE * stream, const dnState s) {
	int i;
	
	/* header */
	fprintf(stream, "\tdenada-State %s %.3f %.3f %d %d\n\n", 
			s->name, s2p(s->init), s2p(s->term), s->transitions, s->durations);
	
	/* transitions */
	for (i = 0; i < s->transitions; i++) {
		fprintf(stream, "\t\t%s %.3f\n", s->state->elem[i], s2p(s->score->elem[i]));
	}
	
	/* emissions */
	for (i = 0; i < s->count; i++) {
		if (i % 4 == 0) fprintf(stream, "\n\t\t");
		fprintf(stream, "%.3f ", s2p(s->emit[i]));
	}
	if (s->durations) fprintf(stream, "\n");
	
	/* durations */
	for (i = 0; i < s->durations; i++) {
		if (i % 5 == 0) fprintf(stream, "\n\t\t");
		fprintf(stream, "%.3f ", s2p(s->duration[i]));
	}
	fprintf(stream, "\n\n");
}

void dnHMM_free (dnHMM hmm) {
	int i;
	ik_free(hmm->name);
	for (i = 0; i < hmm->states; i++) dnState_free(hmm->state[i]);
}

dnHMM dnHMM_read (FILE * stream) {
	int		 i, j, from, to, count;
	double	 score;
	char	 line[1024], tag[256], name[256];
	dnHMM	 hmm = ik_malloc(sizeof(struct dnHMM));
	ik_map	 name_to_number;
	
	/* clear out */
	hmm->name = NULL;
	hmm->states = 0;
	for (i = 0; i < dnSTATE_LIMIT; i++) hmm->state[i] = NULL;
	for (i = 0; i < dnSTATE_LIMIT; i++)
		for (j = 0; j < dnSTATE_LIMIT; j++) hmm->transition[i][j] = MIN_SCORE;
	
	/* find denada-HMM header, skip comments (lines starting with #) */
	while (fgets(line, sizeof(line), stream) != NULL) {
		if (line[0] == '#') continue;
		
		if (sscanf(line, "%s %s %d", tag, name, &hmm->states) == 3) {
			if (strcmp(tag, "denada-HMM") != 0) {
				ik_exit(1, "unrecognized HMM type: %s", line);
			} else {
				break;
			}
		}
	}
	
	/* name */
	hmm->name = ik_malloc(strlen(name) +1);
	strcpy(hmm->name, name);
	
	/* states */
	for (i = 0; i < hmm->states; i++) hmm->state[i] = dnState_read(stream);
	
	/* state names and numbers */
	name_to_number = ik_map_new();
	for (i = 0; i < hmm->states; i++) {
		ik_map_set(name_to_number, hmm->state[i]->name, (void*)i);
	}

	/* transition */
	for (i = 0; i < hmm->states; i++) {
		from = i;
		for (j = 0; j < hmm->state[i]->transitions; j++) {
			to = (size_t)ik_map_get(name_to_number, hmm->state[i]->state->elem[j]);
			score = hmm->state[i]->score->elem[j];
			hmm->transition[from][to] = score;
		}
	}

	/*
	   explicit length state restrictions:
		  - not allowed to be adjacent
		  - not allowed more than one connection on either side
	*/
	for (i = 0; i < hmm->states; i++) {
		if (hmm->state[i]->durations == 0) continue;
		count = 0;
		for (j = 0; j < hmm->states; j++) {
			if (hmm->transition[i][j] != MIN_SCORE) count++;
		}
		if (count > 1) {
			ik_exit(1, "explicit duration states allowed only 1 link");
		}
	}
				
	for (i = 0; i < hmm->states; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (i == j) continue;
			if (hmm->transition[i][j] == MIN_SCORE) continue;
			if (hmm->state[i]->durations && hmm->state[j]->durations)
				ik_exit(1, "explicit duration states must not be adjacent");
		}
	}
	
	
	/* clean up */
	ik_map_free(name_to_number);
	
	return hmm;
}

void dnHMM_write (FILE *stream, const dnHMM hmm) {
	int i;
	
	fprintf(stream, "denada-HMM %s %d\n\n", hmm->name, hmm->states);
	for (i = 0; i < hmm->states; i++) {
		dnState_write(stream, hmm->state[i]);
	}
}

static double emission (const dnState state, const char * seq, int pos) {
	int i, p, c, index;
	
	if (LOGARITHMS) {
		if (pos < state->order) return -2;
		if (seq[pos] == 4)		return -2;
	} else {
		if (pos < state->order) return 0.25;
		if (seq[pos] == 4)		return 0.25;
	}
	
	index = 0;
	for (i = 0; i <= state->order; i++) {
		c = seq[pos -i];
		if (c == 4) {
			if (LOGARITHMS) return -2;
			else			return 0.25;
		}
		p = ik_POWER[4][i];
		index += p * c;
	}
	
	return state->emit[index];
}

struct dnMax {
	double score;
	int	   state;
	int	   jumps;
};

static struct dnMax viterbi_calc (dnMatrix m, int this_coor, int this_state) {
	int		prev_state, prev_coor = this_coor -1;
	struct	dnMax max;
	double	emit, trans, prev, total_score, max_score, max_state;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	max.score = MIN_SCORE;
	max.state = this_state;
	max.jumps = 0; /* normal states do not jump */
	
	emit = emission(hmm->state[this_state], dna->num, this_coor);
	if (emit == MIN_SCORE) return max;
	
	max_score = MIN_SCORE;
	max_state = -1;
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		trans = hmm->transition[prev_state][this_state];
		prev = m->vscore[prev_state][prev_coor];
		
		if (trans == MIN_SCORE) continue;
		if (prev == MIN_SCORE) continue;
		
		total_score = emit + trans + prev;
		if (total_score > max_score) {
			max_score = total_score;
			max_state = prev_state;
		}
	}
	
	if (max_score == MIN_SCORE) return max;
	max.score = max_score;
	max.state = max_state;
	
	return max;
}


static struct dnMax viterbiX_calc (
								  dnMatrix m,
								  const int this_coor,
								  const int this_state)
{
	int		i, prev_state, limit;
	double	emit_score[dnDURATION_LIMIT], state_score[dnDURATION_LIMIT];
	double	prev_score, trans_score, this_score, total_score;
	struct	dnMax max;
	dnHMM	hmm = m->hmm;
	ik_dna	dna = m->dna;
	dnState state = m->hmm->state[this_state];
	
	max.score = MIN_SCORE;
	max.state = -1;
	max.jumps = 0;
	
	/* pre-calculate emit_score for all distances as a cumulative score */
	emit_score[0] = emission(state, dna->num, this_coor);	
	limit = state->durations;
	if (limit > this_coor) limit = this_coor;
	for (i = 1; i < limit; i++) {
		this_score = emission(state, dna->num, this_coor - i);
		prev_score = emit_score[i-1];
		if (this_score == MIN_SCORE) emit_score[i] = MIN_SCORE;
		else if (prev_score == MIN_SCORE) emit_score[i] = MIN_SCORE;
		else emit_score[i] = this_score + prev_score;
	}
	
	/* pre-calculate state_score (emission + duration) */
	for (i = 0; i < limit; i++) {
		if (state->duration[i] == MIN_SCORE) state_score[i] = MIN_SCORE;
		else if (emit_score[i] == MIN_SCORE) state_score[i] = MIN_SCORE;
		else state_score[i] = emit_score[i] + state->duration[i];
	}
	
	/* calculate max for each state and duration */
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		if (prev_state == this_state) continue;
		trans_score = hmm->transition[prev_state][this_state];
		if (trans_score == MIN_SCORE) continue;
		for (i = 0; i < limit; i++) {
			if (state_score[i] == MIN_SCORE) continue;
			prev_score = m->vscore[prev_state][this_coor -i];
			if (prev_score == MIN_SCORE) continue;
			total_score = prev_score + state_score[i] + trans_score;
			if (total_score > max.score) {
				max.score = total_score;
				max.state = prev_state;
				max.jumps = i;
			}
		}
	}
	
	return max;
}


static struct dnMax viterbiP_calc (
								   dnMatrix m, 
								   int this_coor, 
								   int this_state) 
{
	int		prev_state, prev_coor = this_coor -1;
	struct	dnMax max;
	double	emit, total_score, max_score, max_state;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
		
	max.score = 0;
	max.state = this_state;
	max.jumps = 0; /* normal states do not jump */
	
	emit = emission(hmm->state[this_state], dna->num, this_coor);
	
	max_score = 0;
	max_state = -1;
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		total_score = emit * hmm->transition[prev_state][this_state]
			* m->vscore[prev_state][prev_coor];
		if (total_score > max_score) {
			max_score = total_score;
			max_state = prev_state;
		}
	}
	
	if (max_score == 0) return max;
	max.score = max_score;
	max.state = max_state;
	
	return max;
}

static struct dnMax viterbiXP_calc (
								   dnMatrix m, 
								   const int this_coor, 
								   const int this_state) 
{
	int		i, prev_state, limit;
	double	emit_score[dnDURATION_LIMIT], state_score[dnDURATION_LIMIT];
	double	prev_score, this_score, total_score;
	struct	dnMax max;
	dnHMM	hmm = m->hmm;
	ik_dna	dna = m->dna;
	dnState state = m->hmm->state[this_state];
		
	max.score = 0;
	max.state = -1;
	max.jumps = 0;
	
	/* pre-calculate emit_score for all distances as a cumulative score */
	emit_score[0] = emission(state, dna->num, this_coor);	
	limit = state->durations;
	if (limit > this_coor) limit = this_coor;
	for (i = 1; i < limit; i++) {
		this_score = emission(state, dna->num, this_coor - i);
		prev_score = emit_score[i-1];
		emit_score[i] = this_score * prev_score;
	}
	
	/* pre-calculate state_score (emission * duration) */
	for (i = 0; i < limit; i++) {
		state_score[i] = emit_score[i] * state->duration[i];
	}
	
	/* calculate max for each state and duration */
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		for (i = 0; i < limit; i++) {
			prev_score = m->vscore[prev_state][this_coor -i];
			total_score = 
				prev_score 
				* state_score[i] 
				* hmm->transition[prev_state][this_state];
			if (total_score > max.score) {
				max.score = total_score;
				max.state = prev_state;
				max.jumps = i;
			}
		}
	}
	
	return max;
}

static void dnMatrix_viterbi (dnMatrix m) {
	int		i, j, k, max_state;
	double	max_score;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	struct	dnMax max;
		
	/* initialization */
	i = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->init == MIN_SCORE) {
			m->vscore[j][i] = MIN_SCORE;
		} else if (emission(hmm->state[j], dna->num, i) == MIN_SCORE) {
			m->vscore[j][i] = MIN_SCORE;
		} else {
			m->vscore[j][i] = hmm->state[j]->init 
				+ emission(hmm->state[j], dna->num, 0);
		}
		m->vtrace[j][i] = j;
		m->vjumps[j][i] = 0;
	}
	
	/* induction */
	for (i = 1; i < dna->length; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (hmm->state[j]->durations) max = viterbiX_calc(m, i, j);
			else						  max = viterbi_calc(m, i, j);
			m->vscore[j][i] = max.score;
			m->vtrace[j][i] = max.state;
			m->vjumps[j][i] = max.jumps;
		}
	}
	
	/* termination */
	i = dna->length -1;
	max_state = -1;
	max_score = MIN_SCORE;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->term == MIN_SCORE) {
			m->vscore[j][i] = MIN_SCORE;
		} else {
			m->vscore[j][i] += hmm->state[j]->term;
		}
		if (m->vscore[j][i] > max_score) {
			max_score = m->vscore[j][i];
			max_state = j;
		}
	}
	m->vmax = max_score;
		
	/* trace */
	i = dna->length -1;
	j = max_state;			   /* from termination step above */
	ik_ivec_push(m->vpath, j); /* start with max state */
	while (i > 0) {
		if (m->vtrace[j][i] == -1) {
			ik_exit(1, "fatal traceback error");
		} else if (m->vtrace[j][i] == j) {
			ik_ivec_push(m->vpath, j);
			i--;
		} else {
			if (m->vjumps[j][i] == 0) {
				ik_ivec_push(m->vpath, m->vtrace[j][i]);
				j = m->vtrace[j][i];
				i--;
			} else {
				max_state = m->vtrace[j][i];
				for (k = 0; k < m->vjumps[j][i]; k++) {
					ik_ivec_push(m->vpath, j);
				}
				i -= m->vjumps[j][i];
				j = max_state;
			}
		}
	}
}

static void dnMatrix_viterbiP (dnMatrix m) {
	int		i, j, k, max_state;
	double	max_score;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	struct	dnMax max;
	
	/* initialization */
	i = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->init == 0) {
			m->vscore[j][i] = 0;
		} else if (emission(hmm->state[j], dna->num, i) == 0) {
			m->vscore[j][i] = 0;
		} else {
			m->vscore[j][i] = 
			hmm->state[j]->init * emission(hmm->state[j], dna->num, 0);
		}
		m->vtrace[j][i] = j;
		m->vjumps[j][i] = 0;
	}
	
	/* induction */
	for (i = 1; i < dna->length; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (hmm->state[j]->durations) max = viterbiXP_calc(m, i, j);
			else						  max = viterbiP_calc(m, i, j);
			m->vscore[j][i] = max.score;
			m->vtrace[j][i] = max.state;
			m->vjumps[j][i] = max.jumps;
		}
	}
	
	/* termination */
	i = dna->length -1;
	max_state = -1;
	max_score = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->term == 0) {
			m->vscore[j][i] = 0;
		} else {
			m->vscore[j][i] *= hmm->state[j]->term;
		}
		if (m->vscore[j][i] > max_score) {
			max_score = m->vscore[j][i];
			max_state = j;
		}
	}
	m->vmax = max_score;
	
	if (max_state == -1) {
		ik_exit(1, "Viterbi underflowed, -no-logs was a bad idea");
	}
	
	/* trace */
	i = dna->length -1;
	j = max_state;			   /* from termination step above */
	ik_ivec_push(m->vpath, j); /* start with max state */
   
	while (i > 0) {
		if (m->vtrace[j][i] == -1) {
			ik_exit(1, "fatal traceback error");
		} else if (m->vtrace[j][i] == j) {
			ik_ivec_push(m->vpath, j);
			i--;
		} else {
			if (m->vjumps[j][i] == 0) {
				ik_ivec_push(m->vpath, m->vtrace[j][i]);
				j = m->vtrace[j][i];
				i--;
			} else {
				max_state = m->vtrace[j][i];
				for (k = 0; k < m->vjumps[j][i]; k++) {
					ik_ivec_push(m->vpath, j);
				}
				i -= m->vjumps[j][i];
				j = max_state;
			}
		}
	}
}


static double forward_calc (const dnMatrix m, const int coor, const int s1, const int s2) {
	dnHMM  hmm = m->hmm;
	ik_dna dna = m->dna;
	double prev, emit, dura, tran, P;
	int	   i, limit, xduration;
	
	xduration = hmm->state[s1]->durations;
	tran = hmm->transition[s2][s1];
	
	if (!xduration) {
		prev = m->fscore[s2][coor-1];
		emit = emission(hmm->state[s1], dna->num, coor);
		P = prev * emit * tran;
		if (P > 1) ik_exit(1, "P > 1 at f2calc geo");
	} else {
		limit = (coor < xduration) ? coor : xduration;
		emit = 1;
		P = 0;
		for (i = 1; i <= limit; i++) {
			emit *= emission(hmm->state[s1], dna->num, coor -i +1);
			dura = hmm->state[s1]->duration[i-1];
			prev = m->fscore[s2][coor -i];
			P += prev * emit * dura * tran;
		}
		if (P > 1) ik_exit(1, "P > 1 at f2calc exp");
	}
	
	return P;
}

double backward_calc (const dnMatrix m, const int coor, const int s1, const int s2) {
	dnHMM  hmm = m->hmm;
	ik_dna dna = m->dna;
	double prev, emit, dura, tran, P;
	int	   i, limit, exp_length, from_end, jump_state;
	
	exp_length = hmm->state[s2]->durations;
	tran = hmm->transition[s1][s2];
	
	if (!exp_length) {
		prev = m->bscore[s2][coor+1];
		emit = emission(hmm->state[s2], dna->num, coor+1);
		P = prev * emit * tran;
		if (P > 1) ik_exit(1, "P > 1 at b2calc geo");
	} else {
		jump_state = -1;
		for (i = 0; i < hmm->states; i++) {
			if (hmm->transition[s2][i] != MIN_SCORE) {
				jump_state = i;
				break;
			}
		}
		if (jump_state == -1) ik_exit(1, "exp state no not possible");
		
		from_end = dna->length - coor -1;
		limit = (exp_length > from_end) ? from_end : exp_length;
		emit = 1;
		P = 0;
		for (i = 1; i <= limit; i++) {
			emit *= emission(hmm->state[s2], dna->num, coor +i);
			dura = hmm->state[s2]->duration[i-1];
			prev = m->bscore[jump_state][coor+i]; // coor +i +1 ?
			P += prev * emit * dura * tran;
		}
		if (P > 1) ik_exit(1, "P > 1 at b2calc exp");
	}
	
	return P;
}

static void dnMatrix_forwardP (dnMatrix m) {
	int		i, j, k;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	/* initialization */
	i = 0;
	for (j = 0; j < hmm->states; j++) m->fscore[j][i] = 
		hmm->state[j]->init * emission(hmm->state[j], dna->num, i);
	
	/* induction */
	for (i = 1; i < dna->length; i++) {
		for (j = 0; j < hmm->states; j++) {
			m->fscore[j][i] = 0;
			for (k = 0; k < hmm->states; k++) {
				if (hmm->transition[k][j] == MIN_SCORE) continue;
				m->fscore[j][i] += forward_calc(m, i, j, k);
			}
		}
	}
	
	/* termination */
	i = dna->length -1;
	m->forward = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->durations) continue;
		m->forward += m->fscore[j][i] * hmm->state[j]->term;
	}
	
}

static void dnMatrix_backwardP (const dnMatrix m) {
	int		i, j, k;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	/* initialization */
	i = dna->length -1;
	for (j = 0; j < hmm->states; j++) m->bscore[j][i] = hmm->state[j]->term;
	
	/* induction */
	for (i = dna->length -2; i >= 0; i--) {
		for (j = 0; j < hmm->states; j++) {
			m->bscore[j][i] = 0;
			for (k = 0; k < hmm->states; k++) {
				if (hmm->transition[j][k] == MIN_SCORE) continue;
				m->bscore[j][i] += backward_calc(m, i, j, k);
			}
		}
	}
	
	/* termination */
	i = 0;
	m->backward = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->durations) continue;
		m->backward += m->bscore[j][i] * 
			hmm->state[j]->init * emission(hmm->state[j], dna->num, i);
	}
}


static void dnMatrix_posteriorP (const dnMatrix m) {
	int		i, j;
	double	posterior, maxP, maxS;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
		
	for (i = 0; i < dna->length; i++) {
		maxP = 0;
		maxS = -1;
		for (j = 0; j < hmm->states; j++) {
			posterior = m->fscore[j][i] * m->bscore[j][i];
			if (posterior > maxP) {
				maxP = posterior;
				maxS = j;
			}
		}
		if (maxS == -1) ik_exit(1, "dnMatrix_posteriorP");
		ik_ivec_push(m->ppath, maxS);
	}
}

void dnMatrix_free (dnMatrix m) {
	int i;
	
	for (i = 0; i < m->hmm->states; i++) {
		ik_free(m->vscore[i]);
		ik_free(m->vtrace[i]);
		ik_free(m->vjumps[i]);
		ik_free(m->bscore[i]);
		ik_free(m->fscore[i]);
	}
	
	ik_free(m->vscore);
	ik_free(m->vtrace);
	ik_free(m->vjumps);
	ik_free(m->bscore);
	ik_free(m->fscore);
	ik_ivec_free(m->vpath);
	ik_ivec_free(m->ppath);
}

dnMatrix dnMatrix_new (const dnHMM hmm, const ik_dna dna) {
	int		 i, j;
	dnMatrix m = ik_malloc(sizeof(struct dnMatrix));
	
	m->dna = dna;
	m->hmm = hmm;
	
	m->vscore = ik_malloc(sizeof(double) * hmm->states);
	m->vtrace = ik_malloc(sizeof(int)	* hmm->states);
	m->vjumps = ik_malloc(sizeof(int)	* hmm->states);
	m->bscore = ik_malloc(sizeof(double) * hmm->states);
	m->fscore = ik_malloc(sizeof(double) * hmm->states);
	for (i = 0; i < hmm->states; i++) {
		m->vscore[i] = ik_malloc(sizeof(double) * dna->length);
		m->vtrace[i] = ik_malloc(sizeof(int)   * dna->length);
		m->vjumps[i] = ik_malloc(sizeof(int)   * dna->length);
		m->bscore[i] = ik_malloc(sizeof(double) * dna->length);
		m->fscore[i] = ik_malloc(sizeof(double) * dna->length);
	}
	
	for (i = 0; i < hmm->states; i++) {
		for (j = 0; j < dna->length; j++) {
			m->vscore[i][j] = 0;
			m->vtrace[i][j] = -1;
			m->bscore[i][j] = 0;
			m->fscore[i][j] = 0;
		}
	}
	
	m->vpath = ik_ivec_new();
	m->vmax = 0;
	m->ppath = ik_ivec_new();
	
	return m;
}

void dnMatrix_decode (dnMatrix m) {
	if (LOGARITHMS) {
		dnMatrix_viterbi(m);
	} else {
		dnMatrix_viterbiP(m);		 
		dnMatrix_forwardP(m);
		dnMatrix_backwardP(m);
		dnMatrix_posteriorP(m);
		
	}
}

void dnMatrix_output (const dnMatrix m) {
	int i, pos, prev;
	
	printf(">%s\n", m->dna->def);
	printf("Viterbi decoding\n");
	prev = 1;
	for (i = m->vpath->size -1; i > 0; i--) {
		pos = m->vpath->size -i;
		if (m->vpath->elem[i -1] != m->vpath->elem[i]) {
			printf("%d\t%d\t%s\n",
				   prev, pos, m->hmm->state[m->vpath->elem[i]]->name);
			prev = pos +1;
		}
	}
	printf("%d\t%d\t%s\n", prev, m->dna->length,
		   m->hmm->state[m->vpath->elem[0]]->name);
	
	if (LOGARITHMS) return;
	
	printf("Posterior decoding\n");
	prev = 0;
	for (i = 0; i < m->ppath->size -1; i++) {
		if (m->ppath->elem[i+1] != m->ppath->elem[i]) {
			printf("%d\t%d\t%s\n", prev +1, i+1,
				   m->hmm->state[m->ppath->elem[i]]->name);
			prev = i + 1;
		}
	}
	printf("%d\t%d\t%s\n", prev+1, m->dna->length,
		   m->hmm->state[m->ppath->elem[m->ppath->size-1]]->name);
}

void dnMatrix_browse (const dnMatrix m) {
	int		i, j;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	printf("\nMatrix\n");
	printf("DNA ");
	for (j = 0; j < hmm->states; j++) {
		printf(" %1d:%-22s ", j, hmm->state[j]->name);
	}
	printf("\n");
	for (i = 0; i < dna->length; i++) {
		printf("%c%2d  ", dna->seq[i], i);
		for (j = 0; j < hmm->states; j++) {
			printf("%1d", m->vtrace[j][i]);
			printf("/");
			printf("%1d", m->vjumps[j][i]);
			printf("/");
			
			if (m->vscore[j][i] == 0) printf("EMPTY");
			else					  printf("%4.0e", (m->vscore[j][i]));
			printf("/");
			if (m->fscore[j][i] == 0) printf("EMPTY");
			else					  printf("%4.0e", (m->fscore[j][i]));
			 printf("/");
			if (m->bscore[j][i] == 0) printf("EMPTY");
			else					  printf("%4.0e", (m->bscore[j][i]));
			printf("  ");
		}
		printf("\n");
	}
}


void dnMatrix_posterior_states (const dnMatrix m) {
	int		i, j;
	double	posterior, fraction;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	double	sum[dnDURATION_LIMIT];
	
	for (j = 0; j < hmm->states; j++) sum[j] = 0;
	for (i = 0; i < dna->length; i++) {
		for (j = 0; j < hmm->states; j++) {
			posterior = m->fscore[j][i] * m->bscore[j][i] / m->forward;
			sum[j] += posterior;
		}
	}
	
	printf("\nPosterior States\n");
	printf("Pos\t");
	for (j = 0; j < hmm->states; j++) {
		printf("%s\t", hmm->state[j]->name);
	}
	printf("\n");
	
	for (i = 0; i < dna->length; i++) {
		printf("%d\t", i+1);
		
		for (j = 0; j < hmm->states; j++) {
			posterior = m->fscore[j][i] * m->bscore[j][i] / m->forward;
			fraction = posterior / sum[j];
			printf("%.3f\t", fraction);
		}
		printf("\n");
	}
}

void dnMatrix_posterior_coords2 (const dnMatrix m) {
	int		i, j;
	double	 fraction, total;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	printf("\nPosterior Coords\n");
	printf("Pos\t");
	for (j = 0; j < hmm->states; j++) {
		printf("%s\t", hmm->state[j]->name);
	}
	printf("\n");
	
	for (i = 0; i < dna->length; i++) {
		printf("%d\t", i+1);
		
		total = 0;
		for (j = 0; j < hmm->states; j++) {
			total += m->fscore[j][i] * m->bscore[j][i];
		}
		
		for (j = 0; j < hmm->states; j++) {
			fraction = m->fscore[j][i] * m->bscore[j][i] / total;
			printf("%.3f\t", fraction);
		}
		printf("\n");
	}
}

void dnMatrix_posterior_coords (const dnMatrix m) {
	int		i, j;
	double	posterior, fraction, total;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	printf("\nPosterior Coords\n");
	printf("Pos\t");
	for (j = 0; j < hmm->states; j++) {
		printf("%s\t", hmm->state[j]->name);
	}
	printf("\n");
	
	for (i = 0; i < dna->length; i++) {
		printf("%d\t", i+1);
		
		total = 0;
		for (j = 0; j < hmm->states; j++) {
			posterior = m->fscore[j][i] * m->bscore[j][i] / m->forward;
			total += posterior;
		}
		
		for (j = 0; j < hmm->states; j++) {
			posterior = m->fscore[j][i] * m->bscore[j][i] / m->forward;
			fraction = posterior / total;
			printf("%.3f\t", fraction);
		}
		printf("\n");
	}
}

void dnMatrix_testing (dnMatrix m) {
	int i, j;
	double p, sum;
	ik_dna dna = m->dna;
	dnHMM  hmm = m->hmm;
	
	printf("\nchecking, foward=%g, backward=%g\n", m->forward, m->backward);
	for (i = 0; i < dna->length; i++) {
		sum = 0;
		for (j = 0; j < hmm->states; j++) {
			if (hmm->state[j]->durations) continue;
			p = m->fscore[j][i] * m->bscore[j][i];
			sum += p;
		}
		printf("%3d %6e %g\n", i, sum, sum / m->forward);
	}
}



