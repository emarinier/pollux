/*

Pollux
Copyright (C) 2014  Eric Marinier

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

---

Contributions: Oliver Grant

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Encoding.h"

unsigned int encode_Nucleotide(char n1) {
    if (n1 == 'A' || n1 == 'a') { 
        return 0;
    } else if (n1 == 'G' || n1 == 'g') { 
        return 1;
    } else if (n1 == 'C' || n1 == 'c') { 
        return 2;
    } else if (n1 == 'T' || n1 == 't') {
        return 3;
    } else { // used for filling 
        return 0;
    }    
}

void fill_int(unsigned long long int *nucl_seq, char n1) {
    *nucl_seq <<= 2;
    *nucl_seq += encode_Nucleotide(n1);
}

void encode_sequence(struct read* rd, char* seq) {
    if (strchr(seq, '\n')) {
        seq[strlen(seq) - 1] = 0;
    }
    int seq_words = ceil((double)strlen(seq)/(double)32);
    rd->sequence       = malloc(seq_words * sizeof(unsigned long long int));
    // initialize the array to be 0's
    for (int i = 0; i < seq_words; i++) {
        rd->sequence[i] = 0;
    }

    rd->length = strlen(seq);

    
    // for each character in the sequence, encode it
    for (int i = 0; i < rd->length; i++) {
        fill_int(&rd->sequence[i / 32], seq[i]);
    }

    // pad the rest with 0's
    for (int i = rd->length; i < (seq_words*32); i++) {
        fill_int(&rd->sequence[i / 32], 'X');
    }
}
