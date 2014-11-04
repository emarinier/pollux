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

#ifndef ERRORPROCESSING_H
#define	ERRORPROCESSING_H

#include "KMerHashTable.h"
#include "Globals.h"
#include "Reads.h"

#ifdef	__cplusplus
extern "C" {
#endif  

extern unsigned int KMER_SIZE;
    
/**
 * This function will initiate error correcting.
 */    
int processCorrection(int numInputFiles, char* inputFileNames, char* output,
        bool paired, bool substitutions, bool insertions, bool deletions, bool homopolymers,
        bool filtering, bool qualityUpdating);

int convertFASTQToFASTK(int numInputFiles, char* inputFileNames, char* outputDirectory);

void hashSequence(unsigned long long int* sequence, unsigned int sequenceLength,
        KMerHashTable* kmers, unsigned int kmerSize);

void processReads(Reads** reads, unsigned int numInputFiles, 
        KMerHashTable* kmers,  unsigned int kmerSize);

#ifdef	__cplusplus
}
#endif

#endif	/* ERRORPROCESSING_H */

