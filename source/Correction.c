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

#include <stdlib.h>
#include "Correction.h"

Correction* createCorrection(Reads** reads, unsigned int numReadSets,
        KMerHashTable* kmers, unsigned int kmerSize, unsigned int lowKMerThreshold,
        char* outputDirectory, CorrectionFunction correctionFunction)
{
    Correction* correction = (Correction*)malloc(sizeof(Correction));
    
    correction->reads = reads;
    correction->numReadSets = numReadSets;
    
    correction->kmers = kmers;
    correction->kmerSize = kmerSize;
    correction->lowKMerThreshold = lowKMerThreshold;
    
    correction->substitutions = true;
    correction->insertions = true;
    correction->deletions = true;
    correction->homopolymers = true;
    
    correction->outputDirectory = outputDirectory;
    
    correction->correctionFunction = correctionFunction;
    
    return correction;
}

Reads** correctionGetReads(Correction* correction)
{
    return correction->reads;
}

unsigned int correctionGetNumReadSets(Correction* correction)
{
    return correction->numReadSets;
}

KMerHashTable* correctionGetKMers(Correction* correction)
{
    return correction->kmers;
}

unsigned int correctionGetKMerSize(Correction* correction)
{
    return correction->kmerSize;
}

unsigned int correctionGetLowThreshold(Correction* correction)
{
    return correction->lowKMerThreshold;
}

CorrectionFunction correctionGetFunction(Correction* correction)
{
    return correction->correctionFunction;
}

char* correctionGetOutputDirectory(Correction* correction)
{
    return correction->outputDirectory;
}
