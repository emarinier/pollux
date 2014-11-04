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

#include "Reads.h"
#include "KMerHashTable.h"
#include "Utility.h"

#ifndef CORRECTION_H
#define	CORRECTION_H

#ifdef	__cplusplus
extern "C" {
#endif
    
typedef struct Correction Correction;
typedef bool (*CorrectionFunction)(struct read* read, Correction* correction);

struct Correction
{
    Reads** reads;
    unsigned int numReadSets;
    
    KMerHashTable* kmers;
    unsigned int kmerSize;   
    unsigned int lowKMerThreshold;
    
    // Enabled Corrections:
    bool substitutions;
    bool insertions;
    bool deletions;
    bool homopolymers;
    
    bool filtering;
    bool qualityUpdating;
    
    char* outputDirectory;
    
    CorrectionFunction correctionFunction;    // Correction function pointer.
    
};

Correction* createCorrection(Reads** reads, unsigned int numReadSets,
        KMerHashTable* kmers, unsigned int kmerSize, unsigned int lowKMerThreshold,
        char* outputDirectory, CorrectionFunction correctionFunction);
        
Reads** correctionGetReads(Correction* correction);
unsigned int correctionGetNumReadSets(Correction* correction);
KMerHashTable* correctionGetKMers(Correction* correction);
unsigned int correctionGetKMerSize(Correction* correction);
unsigned int correctionGetLowThreshold(Correction* correction);
CorrectionFunction correctionGetFunction(Correction* correction);
char* correctionGetOutputDirectory(Correction* correction);

#ifdef	__cplusplus
}
#endif

#endif	/* CORRECTION_H */

