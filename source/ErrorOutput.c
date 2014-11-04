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

#include "ErrorOutput.h"

#include <time.h>
#include <stdlib.h>
#include "Utility.h"
#include "Counting.h"
#include "Reads.h"
#include "Correction.h"

void writeKMerCounts(FILE* output,
        unsigned long long int* sequence, unsigned int length,
        KMerHashTable* kmers,  unsigned int kmerSize)
{
    unsigned long long int kmer;
    unsigned long long int count;
    
    // Iterate over all k-mers within the read:
    for(int j = 0; j <= (int)length - (int)kmerSize; j++)
    {
        // Get the next k-mer:
        kmer = getKMer(sequence, j, j + kmerSize);
        count = KMerTableLookup(kmers, kmer);
        
        fprintf(output, "%2llu ", count % 100);
    }     
}

void outputKMerCounts(FILE* output,
        Reads** reads, unsigned int numInputFiles,
        KMerHashTable* kmers,  unsigned int kmerSize)
{
    // Variables:
    unsigned int readLength;
    unsigned long long int* sequence;
    unsigned long long int kmer;
    unsigned long long int count;
    
    // Iterate over all files:
    for(int file = 0; file < numInputFiles; file++)
    {
        readsReset(reads[file]);
        
        // Iterate over all reads:
        while(readsHasNext(reads[file]))
        {
            struct read* current = readsGetNext(reads[file]);

            fprintf(output, "%d (%c) : \n", current->number, '!');

            readLength = current->length;
            sequence = current->sequence;

            // Iterate over all k-mers within the read:
            for(int j = 0; j <= (int)readLength - (int)kmerSize; j++)
            {
                // Get the next k-mer:
                kmer = getKMer(sequence, j, j + kmerSize);

                count = KMerTableLookup(kmers, kmer);

                // Safety check for existence of k-mer:
                // Exists:
                if(count != 0)
                {
                    fprintf(output, "%2llu ", count % 100);
                }
                // Doesn't exist -- ERROR:
                else
                {
                    fprintf(output, "\nERROR - KMER NOT RECOGNIZED!\n");
                    fprintf(output, "%llu\n", kmer);
                    writeAsNucleotides(output, &kmer, 0, kmerSize);
                    fprintf(output, "\n");

                    return;
                }
            }        

            fprintf(output, "\n");

            // Nucleotides:
            writeAsNucleotidesSpaced(output, sequence, 0, readLength);
            fprintf(output, "\n");

            // Indices:
            for(int j = 0; j < (int)readLength; j++)
            {
                fprintf(output, "%2d ", j % 100);
            }

            fprintf(output, "\n");
        }        
    }        
}

void outputRead(FILE* output, struct read* read)
{
    fprintf(output, "%s", read->seqName1);
    writeAsNucleotides(output, read->sequence, 0, read->length); fprintf(output, "\n"); 
    fprintf(output, "%s", read->seqName2);
    fprintf(output, "%s\n", read->quality);
}

void outputReads(FILE* output, Reads* reads)
{   
    readsReset(reads);
    
    // Iterate over all reads:
    while(readsHasNext(reads))
    {   
        outputRead(output, readsGetNext(reads));
    }   
}

void outputReadFASTK(FILE* output, struct read* read, KMerHashTable* kmers, unsigned int kmerSize)
{
    fprintf(output, "%s", read->seqName1);
    writeAsNucleotides(output, read->sequence, 0, read->length); fprintf(output, "\n"); 
    fprintf(output, "%s", read->seqName2);
    fprintf(output, "%s\n", read->quality);
    
    // Is the read long enough to have k-mers of this size?
    // NO:
    if(read->length < kmerSize)
    {
        fprintf(output, "0");
    }
    // YES:
    else
    {
       // K-mers:
        int total = read->length - kmerSize + 1;
        unsigned int readKMers[total];

        // Get kmer counts.
        getKMerCounts(read->sequence, read->length, kmers, kmerSize, readKMers);

        for(int i = 0; i < total; i++)
        {
            fprintf(output, "%d ", readKMers[i]);
        } 
    }    
    
    fprintf(output, "\n");
}

void outputReadsFASTK(FILE* output, Reads* reads, KMerHashTable* kmers, unsigned int kmerSize)
{   
    readsReset(reads);
    
    // Iterate over all reads:
    while(readsHasNext(reads))
    {   
        outputReadFASTK(output, readsGetNext(reads), kmers, kmerSize);
    }   
}

