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

#include "Counting.h"
#include "Utility.h"

void getKMerCounts(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, unsigned int* counts)
{
    // Variables:
    unsigned long long int kmer;
    int total = length - kmerSize + 1;
    
    // Iterate over all k-mers within the read:
    for(int i = 0; i < total; i++)
    {
        // Get the next k-mer:
        kmer = getKMer(sequence, i, i + kmerSize);

        // Get the count:     
        counts[i] = KMerTableLookup(kmers, kmer);
    }
}

unsigned int areCountsBelowThreshold(unsigned int* counts, unsigned int start, 
        unsigned int end, const unsigned int THRESHOLD)
{
    for(int i = start; i < end; i++)
    {        
        if(counts[i] > THRESHOLD)
        {
            return 0;
        }
    }
    
    return 1;
}
