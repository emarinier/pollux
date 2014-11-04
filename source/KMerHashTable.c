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
#include "KMerHashTable.h"
#include "Utility.h"
#include "Correction.h"

static inline void printProgress(int x, int n, int r)
{
    // Only update r times.
    if (n/r == 0 || x % (n/r) != 0) return; // Mod 0 causes problems.
    
    double result;
    
    if(n / r != 1)
        result = x / (n/r) * (100/r);
    else
        result = (x + 1) / ((double)n/(double)r) * ((double)100/(double)r);
 
    printf("%d%% ", (int)result);
    fflush(stdout);
}

void addKMersToTable(KMerHashTable* table, unsigned long long int* sequence, 
        unsigned int sequenceLength, unsigned int kmerSize)
{
    //Variables:
    unsigned long long int kmer;
    unsigned long long int count;
    
    for(int i = 0; i <= (int)sequenceLength - (int)kmerSize; i++)
    {
        // Get the next k-mer:
        kmer = getKMer(sequence, i, i + kmerSize);
        
        // Does the k-mer exist?
        count = (unsigned long long int)hash_table_lookup(table->table, (HashTableKey)(kmer));

        // The k-mer does exist:            
        if(count != 0)
        {                    
            // Update:
            KMerTableInsert(table, kmer, (count + 1));                
        }
        // The k-mer doesn't exist:
        else
        {
            // Initialize:
            KMerTableInsert(table, kmer, 1);
        }
    }    
}

void preprocessKMers(KMerHashTable* kmerTable, Correction* correction)
{
    // Data structures:
    HashTable* hashTable = kmerTable->table;
    HashTableIterator iterator;
    
    unsigned long long int MAX_KMER_COUNT = (1024 + 1);
    unsigned long long int counts[MAX_KMER_COUNT];
    
    unsigned long long int kmer;    
    unsigned long long int count;
    
    unsigned long long int current = 0;
    unsigned long long int total = hash_table_num_entries(hashTable);
    
    unsigned long long int unique = 0;
    
    hash_table_iterate(hashTable, &iterator);
    
    // Initialize Counts:
    for (int i = 0; i < MAX_KMER_COUNT; i++)
    {
        counts[i] = 0;
    }
    
    // Iterate Over K-Mers:
    while(hash_table_iter_has_more(&iterator))
    {
        printProgress(current, total, 20);
        current++;
        
        kmer = (unsigned long long int)hash_table_iter_next_key(&iterator);
        count = (unsigned long long int)hash_table_lookup(hashTable, (HashTableKey)(kmer));

        // Tally Counts:
        if(count <= MAX_KMER_COUNT)
        {
            counts[count] += 1;
        }
        
        // Remove Unique:
        if(count == 1)
        {
            hash_table_remove(hashTable, (HashTableKey)(kmer), false);
            unique++;
        }
    }
    
    // Unique K-Mer Information:
    printf("\n");    
    printf("Removed %llu unique k-mers from the set of %llu total k-mers.\n", unique, total);
    
    printf("Resizing...\n");
    hash_table_resize(hashTable);
    printf("Finished resizing...\n");
    
    unsigned int currentKMerCount = 1;
        
    // Loop until we find a low-to-high number transitions:
    while (currentKMerCount <= MAX_KMER_COUNT && counts[currentKMerCount] > counts[currentKMerCount + 1])
    {
        currentKMerCount++;
    }
    
    if(currentKMerCount < MAX_KMER_COUNT)
    {
        correction->lowKMerThreshold = currentKMerCount;
    }
    
   printf("Low k-mer count value was observed to be %d.\n", currentKMerCount);
}

int KMerHashFunctionEquals(HashTableValue value1, HashTableValue value2)
{
    return (value1 == value2);
}

unsigned long KMerHashFunctionHash(HashTableKey kmer)
{
    return (unsigned long)kmer;
}

KMerHashTable* newKMerHashTable()
{
    HashTable* hashTable = hash_table_new(KMerHashFunctionHash, KMerHashFunctionEquals);    
    KMerHashTable* kmerTable;
    
    if((kmerTable = malloc(sizeof *kmerTable)) != NULL)
    {
        kmerTable->table = hashTable;
    }

    return kmerTable;
}

int KMerTableInsert(KMerHashTable* kmerTable, unsigned long long int kmer, unsigned long long int count)
{
    return hash_table_insert(kmerTable->table, (HashTableKey)(kmer), (HashTableValue)(count));
}

/**
 * 
 * BE EXTREMELY CAREFUL USING THIS LOCALLY!!
 * 
 * SAFER: count = (unsigned long long int)hash_table_lookup(kmerTable->table, (HashTableKey)(kmer));
 *      - Might miss singletons.
 * 
 */
unsigned long long int KMerTableLookup(KMerHashTable* kmerTable, unsigned long long int kmer)
{
    unsigned long long int count;    
    HashTableValue value = hash_table_lookup(kmerTable->table, (HashTableKey)(kmer));
    
    count = (unsigned long long int)value;
    
    if(count == 0)
    {
        count = 1;
    }
    
    return count;
}

unsigned int getMaxKMerCount(KMerHashTable* kmerTable)
{
    HashTable* hashTable = kmerTable->table;
    HashTableIterator iterator;
    unsigned long long int kmer;
    
    unsigned long long int current = 0;
    unsigned long long int max = 0;    
    
    hash_table_iterate(hashTable, &iterator);
    
    while(hash_table_iter_has_more(&iterator))
    {
        kmer = (unsigned long long int)hash_table_iter_next_key(&iterator);
        current = KMerTableLookup(kmerTable, kmer);
        
        max = getMax(current, max);
    }
    
    return max;
}

unsigned int* createDistribution(KMerHashTable* kmerTable, unsigned int max)
{    
    // Data structures:
    HashTable* hashTable = kmerTable->table;
    HashTableIterator iterator;
    unsigned int* distribution = (unsigned int*)malloc((max + 1) * sizeof(unsigned int*));
    
    unsigned long long int kmer;    
    unsigned long long int current;
    
    // Initialize:
    for(int i = 0; i <= max; i++)
    {
        distribution[i] = 0;
    }
    
    hash_table_iterate(hashTable, &iterator);
    
    while(hash_table_iter_has_more(&iterator))
    {
        kmer = (unsigned long long int)hash_table_iter_next_key(&iterator);
        current = KMerTableLookup(kmerTable, kmer);

        distribution[current]++;
    }
    
    return distribution;
}

unsigned int getNumRepeats(KMerHashTable* kmerTable)
{   
    // Variables:
    unsigned int repeats = 0;
    unsigned int max = getMaxKMerCount(kmerTable);
    
    // Data Structures:
    unsigned int* distribution = createDistribution(kmerTable, max);    

    for(int i = 2; i <= max; i++)
    {
        repeats += distribution[i] * i;
    }
    
    free(distribution);
    
    return repeats;
}


