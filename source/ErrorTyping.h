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

#ifndef ERRORTYPING_H
#define	ERRORTYPING_H

#include "KMerHashTable.h"
#include "Utility.h"
#include "Reads.h"

#ifdef	__cplusplus
extern "C" {
#endif

#define UNKNOWN '?'
#define HIGH_QUALITY 'H'
#define LOW_COVERAGE 'L'
#define HOMOPOLYMER 'P'
#define ERROR 'E'
#define INTERNAL_ERROR 'I'
#define EXTERNAL_ERROR 'X'
#define CORRECTED 'C'
#define BAD 'B'

/**
 * Determines whether or not the distance between values is large enough to be 
 * considered a jump or a discontinuity.
 * 
 * @param value1
 * @param value2
 * @return Zero if not a jump, non-zero if a jump.
 */
unsigned int isJump(unsigned int value1 , unsigned int value2);

/**
 * Determines whether or not the values in the array in the specified range are 
 * below the threshold or not.
 * 
 * @param counts The array of counts.
 * @param start The starting index.
 * @param end The ending index.
 * @param THRESHOLD The treshold value.
 * @return Zero if at least a single item is strictly greater than the threshold,
 *      non-zero if all the items are less than or equal to the threshold.
 */
unsigned int areCountsBelowThreshold(unsigned int* counts, unsigned int start, 
        unsigned int end, const unsigned int THRESHOLD);

/**
 * Determines whether or not the given sequence will be typed as an external 
 * single nucleotide error.
 * 
 * Individual typings are not necessarily mutually exclusive from other typings.
 * 
 * @param sequence
 * @param length
 * @param kmers
 * @param kmerSize
 * @param THRESHOLD
 * @return 
 */
unsigned int typeExternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

/**
 * Determines whether or not the given sequence will be typed as an internal 
 * single nucleotide error.
 * 
 * Individual typings are not necessarily mutually exclusive from other typings.
 * 
 * @param sequence
 * @param length
 * @param kmers
 * @param kmerSize
 * @param THRESHOLD
 * @return 
 */
unsigned int typeInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

/**
 * Determines whether or not the given sequence will be typed as high quality 
 * read.
 * 
 * Individual typings are not necessarily mutually exclusive from other typings.
 * 
 * @param sequence
 * @param length
 * @param kmers
 * @param kmerSize
 * @param THRESHOLD
 * @return 
 */
unsigned int typeHighQuality(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

/**
 * Determines whether or not the given sequence will be typed as low coverage 
 * read.
 * 
 * Individual typings are not necessarily mutually exclusive from other typings.
 * 
 * @param sequence
 * @param length
 * @param kmers
 * @param kmerSize
 * @param THRESHOLD
 * @return 
 */
unsigned int typeLowCoverage(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

/**
 * Types a given sequence.
 * 
 * @param sequence The sequence to type as a pointer to the first item.
 * @param sequenceLength The length of the sequence in nucleotide bases.
 * @param kmers The kmer data structure.
 * @param kmerSize The length of the kmers.
 * @return A single character representation of the type.
 */
char typeSequence(unsigned long long int* sequence, unsigned int sequenceLength,
        KMerHashTable* kmers, unsigned int kmerSize);

/**
 * Types all the given reads.
 * 
 * @param reads The reads to type.
 * @param numReads The number of reads in the reads data structure.
 * @param kmers The kmer data structure.
 * @param kmerSize The length of the kmers.
 */
void typeReads(struct read** reads, unsigned int numReads,
        KMerHashTable* kmers,  unsigned int kmerSize);

int getInternalErrorPosition(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

int getExternalErrorPosition(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

void printDistributionOfTypes(Reads** reads, unsigned int numInputFiles);

void printDistributionOfInternalErrorLengths(Reads** reads, unsigned int numInputFiles, 
        KMerHashTable* kmers,  unsigned int kmerSize);

int getStartOfInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

int getLengthOfInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD);

char typeSequence(unsigned long long int* sequence, unsigned int sequenceLength,
        KMerHashTable* kmers, unsigned int kmerSize);

#ifdef	__cplusplus
}
#endif

#endif	/* ERRORTYPING_H */

