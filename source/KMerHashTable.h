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

#include "hash-table.h"

#ifndef KMERHASHTABLE_H
#define	KMERHASHTABLE_H

#ifdef	__cplusplus
extern "C" {
#endif
    
typedef struct Correction Correction;
    
typedef struct
{
    HashTable* table;
} KMerHashTable;

/**
 * Creates a new KMerHashTable.
 * 
 * @return 
 */
KMerHashTable* newKMerHashTable();

/**
 * Insters count at the location associated with kmer. 
 * 
 * @param kmerTable The kmer table to work with.
 * @param kmer The kmer value.
 * @param count The count value to store.
 * @return Whether or not the insert was successful.
 */
int KMerTableInsert(KMerHashTable* kmerTable, unsigned long long int kmer, unsigned long long int count);

/**
 * Does a lookup and returns the count value stored at the location associated 
 * with the kmer.
 * 
 * @param kmerTable The kmer table to work with.
 * @param kmer The kmer value.
 * @return The count associated with the kmer.
 */
unsigned long long int KMerTableLookup(KMerHashTable* kmerTable, unsigned long long int kmer);

void addKMersToTable(KMerHashTable* table, unsigned long long int* sequence, 
        unsigned int sequenceLength, unsigned int kmerSize);

/**
 * This function returns the maximum k-mer count for the k-mer hash table.
 * 
 * ex:
 * 
 * kmerA -> 1 instances
 * kmerB -> 5 instances
 * kmerC -> 3 instances
 * 
 * result: 5
 * 
 * @param kmerTable The k-mer table to work with.
 * @return The maximum k-mer counts.
 */
unsigned int getMaxKMerCount(KMerHashTable* kmerTable);

/**
 * Creates a distribution of k-mer counts, such that the indices correspond to 
 * the number of k-mers that have that count.
 * 
 * ex:
 * 
 * kmerA -> 1 instances
 * kmerB -> 5 instances
 * kmerC -> 3 instances
 * kmerD -> 3 instances
 * 
 * result: [0, 1, 0, 2, 0, 1]
 * 
 * @param kmerTable The k-mer table to work with.
 * @param max The maximum k-mer counts.
 * @return The k-mer count distribution.
 */
unsigned int* createDistribution(KMerHashTable* kmerTable, unsigned int max);

/**
 * This function determines the number of repeats k-mers in the k-mer hash table.
 * 
 * @param kmerTable The k-mer table to work with.
 * @return The number of repeats.
 */
unsigned int getNumRepeats(KMerHashTable* kmerTable);

void preprocessKMers(KMerHashTable* kmerTable, Correction* correction);


#ifdef	__cplusplus
}
#endif

#endif	/* KMERHASHTABLE_H */

