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

*/

#ifndef ERROROUTPUT_H
#define	ERROROUTPUT_H

#include <stdio.h>
#include "KMerHashTable.h"
#include "Utility.h"
#include "Reads.h"

#ifdef	__cplusplus
extern "C" {
#endif
    
void writeKMerCounts(FILE* output,
        unsigned long long int* sequence, unsigned int length,
        KMerHashTable* kmers,  unsigned int kmerSize);

void outputKMerCounts(FILE* output,
        Reads** reads, unsigned int numInputFiles,
        KMerHashTable* kmers,  unsigned int kmerSize);

void outputReads(FILE* output, Reads* reads);

void outputRead(FILE* output, struct read* read);

void outputReadsFASTK(FILE* output, Reads* reads, KMerHashTable* kmers, unsigned int kmerSize);

#ifdef	__cplusplus
}
#endif

#endif	/* ERROROUTPUT_H */

