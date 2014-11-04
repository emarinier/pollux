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

#include <stdio.h>

#ifndef UTILITY_H
#define	UTILITY_H

#ifdef	__cplusplus
extern "C" {
#endif
    
typedef int bool;
#define true 1
#define false 0

#define MAX_CORRECTIONS 30                      // 'Minimum' max corrections performed.
                                                // - More important for short reads.

struct Sequence
{
    unsigned long long int* sequence;           // Pointer to the actual sequence.
    char* quality;                              // Pointer to the quality scores.
    
    unsigned int length;                        // Number of nucleotides in sequence.
    unsigned int blocks;                        // The number of blocks of memory used.
    unsigned int numCorrections;    
    char* corrections;
    int* homopolymerSize;                       // Associated with corrections.
    char type;                                  // Sequence (error) type.
};

/**
 * This function returns a k-mer from the passed sequence array as a single 
 * 64-bit integer. Consequently, the means that the k-mer size cannot be larger 
 * than 32 nucleotide bases.
 * 
 * @param sequence The sequence array from which to pull the k-mer.
 * @param startNucleotidePosition The starting position of the k-mer in 
 *      nucleotide bases, relative to the entire sequence (with sequence[0] 
 *      base 0 being the start). The starting base is included in the k-mer.
 * @param endNucleotidePosition  The ending position of the k-mer in 
 *      nucleotide bases, relative to the entire sequence (with sequence[0] 
 *      base 0 being the start). The ending base is not included in the k-mer.
 * @return A 64-bit integer representation of the k-mer with the first k high 
 *      bits representing the k-mer and the 64 - (2 * k) lower bits set to 0.
 */    
unsigned long long int getKMer(unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition);

/**
 * This function will print to standard output a sequence of nucleotides 
 * defined by the function call. It will not print an end of line character.
 * 
 * @param sequence The sequence array from which to print the nucleotides.
 * @param startNucleotidePosition The starting position of the printing in 
 *      nucleotide bases, relative to the entire sequence (with sequence[0] 
 *      base 0 being the start). The starting base is included in the print.
 * @param endNucleotidePosition The ending position of the printing in 
 *      nucleotide bases, relative to the entire sequence (with sequence[0] 
 *      base 0 being the start). The ending base is not included in the print.
 */
void printAsNucleotides(unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition);

void writeAsNucleotides(FILE* file, unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition);
    
/**
 * This function will print to standard output a sequence of nucleotides 
 * defined by the value passed to the function call. It will interpret the 
 * entire value as being valid nucleotide data. This will print 32 nucleotide 
 * bases to standard output.
 * 
 * @param value The 64-bit representation of 32 nucleotide bases.
 */
void printValueAsNucleotides(unsigned long long int value);

/**
 * This function will print to the specified file a sequence of nucleotides 
 * defined by the value passed to the function call. It will interpret the 
 * entire value as being valid nucleotide data. This will print 32 nucleotide 
 * bases to standard output.
 * 
 * @param file The file to write to. Use "stdout" for standard output.
 * @param value The 64-bit representation of the 32 nucleotide bases.
 */
void writeValueAsNucleotides(FILE* file, unsigned long long int value);

/**
 * This is a helper function to get max of two numbers.
 * 
 * @param value1
 * @param value2
 * @return The larger of the two numbers.
 */
int getMax(int value1, int value2);

/**
 * This is a helper function to get min of two numbers.
 * 
 * @param value1
 * @param value2
 * @return The smaller of the two numbers.
 */
int getMin(int value1, int value2);

/**
 * This function returns the reverse for a single 64bit sequence block.
 * 
 * @param sequence The single 64bit sequence VALUE (not pointer) to reverse.
 * @return The VALUE of the reverse.
 */
unsigned long long int getReverse(unsigned long long int sequence);

/**
 * This function returns the reverse compliment for the passed sequence. This 
 * function will allocate memory. The original sequence will be unchanged. The 
 * remaining bits in the last block will be masked out with 0's.
 * 
 * @param sequence The sequence to reverse compliment.
 * @param length The length of the sequence (>= 1).
 * @return 
 */
unsigned long long int* createReverseCompliment(unsigned long long int* sequence,
        unsigned int length);

char getBase(unsigned long long int* nucleotideSequence, 
        unsigned int nucleotidePosition);

void setBase(unsigned long long int* nucleotideSequence, 
        unsigned int nucleotidePosition, char nucleotide);

void changeBase(struct Sequence* sequence, unsigned int nucleotidePosition, 
        char nucleotide, char quality);

void deleteBase(struct Sequence* sequence, unsigned int nucleotidePosition);

void insertBase(struct Sequence* sequence, unsigned int nucleotidePosition, 
        char nucleotide, char quality);

unsigned int getNumMemoryBlocks(unsigned int length);

void writeAsNucleotidesSpaced(FILE* file, unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition);

int getHomopolymerLength(unsigned long long int* nucleotideSequence, 
        int sequenceLength, unsigned int nucleotidePosition);

int getHomopolymerLeftmostNucleotide(unsigned long long int* nucleotideSequence, 
        unsigned int nucleotidePosition);

void setHomopolymerLength(struct Sequence* sequence, 
        unsigned int nucleotidePosition, unsigned int size, char quality);

#ifdef	__cplusplus
}
#endif

#endif	/* UTILITY_H */

