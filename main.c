/*****************************************************************************\
 denada.c

Duration Explicit Nucleic Acid Decoding Algorithm

Copyright (C) 2003-2006 Ian Korf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "denada.h"
#include "sequence.h"
#include "toolbox.h"

int main (int argc, char *argv[]) {
	ik_pipe   hmm_file;
	ik_pipe   fasta_file;
	ik_fasta  fasta;
	ik_dna    dna;
	dnHMM     hmm;
	dnMatrix  matrix;
	
	/* set the program name */
	ik_set_program_name(argv[0]);
    
	/* unadvertised options */
	ik_register_option("-show-dna", 0);
	ik_register_option("-show-hmm", 0);
	ik_register_option("-show-matrix", 0);
	ik_register_option("-posterior-states", 0);
	ik_register_option("-posterior-coords", 0);
	ik_register_option("-no-logs", 0);
	ik_register_option("-test-posterior", 0);
	ik_parse_options(&argc, argv);
	
	if (argc != 3) {
		printf("denada - Duration Explicit Nucleic Acid Decoding Algorithm\n\n");
		printf("usage: %s <HMM file> <FASTA file>\n", argv[0]);
		exit(1);
	}
	
	hmm_file   = ik_pipe_open(argv[1], "r");
	fasta_file = ik_pipe_open(argv[2], "r");
    
	if (ik_option("-no-logs")) {
		/*fprintf(stderr, "logarithms turned off, beware underflow\n");*/
		dnLogarithm_set(0);
	}
	hmm = dnHMM_read(hmm_file->stream);
	ik_pipe_close(hmm_file);
	
	if (ik_option("-show-hmm")) dnHMM_write(stdout, hmm);
	
	/* process files */
	while ((fasta = ik_fasta_read(fasta_file->stream)) != NULL) {
		dna = ik_dna_new(fasta->def, fasta->seq);
		if (ik_option("-show-dna")) printf("%s\n", dna->seq);
		matrix = dnMatrix_new(hmm, dna);
		dnMatrix_decode(matrix);
        dnMatrix_output(matrix);
        
        /*if (ik_option("-no-logs")) {*/
            if (ik_option("-show-matrix")) dnMatrix_browse(matrix);
            if (ik_option("-posterior-states")) dnMatrix_posterior_states(matrix);
            if (ik_option("-posterior-coords")) dnMatrix_posterior_coords(matrix);
            if (ik_option("-test-posterior")) dnMatrix_testing(matrix);
		/*}*/
        		
		/* clean up */
		dnMatrix_free(matrix);
		ik_dna_free(dna);
		ik_fasta_free(fasta);
	}
	
	/* clean up */
	dnHMM_free(hmm);
	ik_pipe_close(fasta_file);
	return 0;
}
