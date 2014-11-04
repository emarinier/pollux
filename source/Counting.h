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

#include "KMerHashTable.h"

#ifndef COUNTING_H
#define	COUNTING_H

#ifdef	__cplusplus
extern "C" {
#endif

/**
 * This function fills the passed array with the number of the occurances of a 
 * given k-mer over the entire length of the sequence.
 * 
 * @param sequence The sequence to get the k-mer counts for.
 * @param length The length of the sequence.
 * @param kmers The k-mer hash table data structure.
 * @param kmerSize The length of the kmers.
 * @param counts The counts array to fill. There will be (length - kmerSize + 1) 
 *      entries expected to be filled.
 */    
void getKMerCounts(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, unsigned int* counts);

/**
 * This function determines whether or not all the entries between the specified 
 * range in the passed array are below the given threshold.
 * 
 * @param counts The array to examine.
 * @param start The starting index.
 * @param end The ending index.
 * @param THRESHOLD The threshold. Less than or equal to the threshold!
 * 
 * @return Whether or not the values are less than or equal to the threshold.
 */
unsigned int areCountsBelowThreshold(unsigned int* counts, unsigned int start, 
        unsigned int end, const unsigned int THRESHOLD);


#ifdef	__cplusplus
}
#endif

#endif	/* COUNTING_H */

