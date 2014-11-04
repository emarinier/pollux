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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Utility.h"

unsigned int getNumMemoryBlocks(unsigned int length)
{
    const int NUCLEOTIDES_PER_BLOCK = 32;
    unsigned int blocks;
    
    
    // Determine number of blocks:
    blocks = length / NUCLEOTIDES_PER_BLOCK; // Number of full blocks.

    if(length % NUCLEOTIDES_PER_BLOCK > 0)
    {
        blocks++;   // Partially filled block.
    }
    
    return blocks;
}

unsigned long long int getKMer(unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition)
{
    const unsigned long long int MASK = 0xFFFFFFFFFFFFFFFF;
    
    unsigned int startSequenceIndex = startNucleotidePosition / 32;
    unsigned int endNucleotideIndex = (endNucleotidePosition - 1) / 32;
    
    unsigned long long int startSequence = sequence[startSequenceIndex];
    unsigned long long int endSequence = sequence[endNucleotideIndex];
    
    unsigned long long int result;
    
    result = startSequence << ((startNucleotidePosition * 2) % 64);
    result = result & (MASK << ((startNucleotidePosition * 2) % 64));
    
    // Do we need to grab the bits at the start of the second sequence and
    // append it to the end of the result?
    if(startSequenceIndex != endNucleotideIndex)
    {
        result = result + (endSequence >> (64 - (startNucleotidePosition * 2) % 64));
    }
    
    //Mask away the extra bits:
    result = result & (MASK << (64 - (endNucleotidePosition - startNucleotidePosition) * 2));
    
    return result;
}

unsigned long long int getReverse(unsigned long long int sequence)
{
    unsigned long long int MASK = 0x3;
    
    unsigned long long int nucleotide;
    unsigned long long int reverse = 0x0;

    
    for(int i = 0; i < 32; i++)
    {
        // Get the last nucleotide:
        nucleotide = (sequence >> (i * 2)) & MASK;        
                
        // Add the nucleotide to the reverse:
        reverse = reverse | (nucleotide << ((32 - i - 1) * 2));
    }  
    
    return reverse;
}

unsigned long long int* createReverseCompliment(unsigned long long int* sequence,
        unsigned int length)
{
    unsigned int endNucleotidePosition = getMax(0, (length - 1));
    unsigned int endNucleotideIndex = endNucleotidePosition / 32;
    unsigned int extraBases = (32 - length) % 32;
     
    unsigned long long int reverse;
    
    unsigned long long int* reverseCompliment = 
        (unsigned long long int*)malloc((endNucleotideIndex + 1) * sizeof(unsigned long long int*));
    
    for(int i = 0; i <= endNucleotideIndex; i++)
    {
        reverse = getReverse(sequence[endNucleotideIndex - i]);
        reverseCompliment[i] = ~(reverse);
    }
    
    //Collapse:
    for(int i = 0; i <= endNucleotideIndex; i++)
    {
        //Shuffle:
        reverseCompliment[i] = reverseCompliment[i] << (extraBases * 2);
        
        //Fill trailing bits if necessary:
        if(i < endNucleotideIndex && extraBases > 0)
        {            
            reverseCompliment[i] = reverseCompliment[i] 
                | (reverseCompliment[i + 1] >> ((32 - extraBases) * 2));
        }
    }
    
    return reverseCompliment;
}

void printAsNucleotides(unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition)
{
    writeAsNucleotides(stdout, sequence, startNucleotidePosition, 
            endNucleotidePosition);
}

void writeAsNucleotides(FILE* file, unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition)
{
    const unsigned long long int MASK = 0x3;

    // Iterate over all nucleotides:
    for(int i = startNucleotidePosition; i < endNucleotidePosition; i++)
    {
        // Determine the current sequence:
        unsigned long long int currentSequence = sequence[i / 32];
        unsigned int position = (i * 2) % 64;     // bit position w/i sequence
        
        unsigned long long int nucleotide = (currentSequence >> (62 - position)) & MASK;
        
        // Print appropriate character:
        if(nucleotide == 0x0)
        {
            fprintf(file, "A");
        }
        else if(nucleotide == 0x1)
        {
            fprintf(file, "G");
        }
        else if(nucleotide == 0x2)
        {
            fprintf(file, "C");
        }
        else if(nucleotide == 0x3)
        {
            fprintf(file, "T");
        }
        // Error:
        else
        {
            fprintf(file, "E");
        }
    }
}

void writeAsNucleotidesSpaced(FILE* file, unsigned long long int* sequence, 
        unsigned int startNucleotidePosition, 
        unsigned int endNucleotidePosition)
{
    const unsigned long long int MASK = 0x3;

    // Iterate over all nucleotides:
    for(int i = startNucleotidePosition; i < endNucleotidePosition; i++)
    {
        // Determine the current sequence:
        unsigned long long int currentSequence = sequence[i / 32];
        unsigned int position = (i * 2) % 64;     // bit position w/i sequence
        
        unsigned long long int nucleotide = (currentSequence >> (62 - position)) & MASK;
        
        // Print appropriate character:
        if(nucleotide == 0x0)
        {
            fprintf(file, " A ");
        }
        else if(nucleotide == 0x1)
        {
            fprintf(file, " G ");
        }
        else if(nucleotide == 0x2)
        {
            fprintf(file, " C ");
        }
        else if(nucleotide == 0x3)
        {
            fprintf(file, " T ");
        }
        // Error:
        else
        {
            fprintf(file, " E ");
        }
    }
}

void printValueAsNucleotides(unsigned long long int value)
{
    writeValueAsNucleotides(stdout, value);
}

void writeValueAsNucleotides(FILE* file, unsigned long long int value)
{
    unsigned long long int mask = 0x3;
    
    unsigned long long int nucleotide;
    
    for(int i = 0; i < 32; i++)
    {
        nucleotide = (value >> ((32 - i - 1) * 2)) & mask;
                
        if(nucleotide == 0x0)
        {
            fprintf(file, "A");
        }
        else if(nucleotide == 0x1)
        {
            fprintf(file, "G");
        }
        else if(nucleotide == 0x2)
        {
            fprintf(file, "C");
        }
        else if(nucleotide == 0x3)
        {
            fprintf(file, "T");
        }
        else
        {
            fprintf(file, "[E]");
        }
    }
}

int getMax(int value1, int value2)
{
    if(value1 > value2)
    {
        return value1;
    }
    else
    {
        return value2;
    }
}

int getMin(int value1, int value2)
{
    if(value1 < value2)
    {
        return value1;
    }
    else
    {
        return value2;
    }
}

char getBase(unsigned long long int* nucleotideSequence, 
        unsigned int nucleotidePosition)
{
    unsigned long long int mask = 0x3;
    unsigned long long int result;
    char base = 'E';

    unsigned int sequenceIndex = nucleotidePosition / 32;    
    unsigned long long int sequence = nucleotideSequence[sequenceIndex];
    
    unsigned int shift = (31 - (nucleotidePosition % 32)) * 2;
    
    result = (sequence >> shift) & mask;
    
    if(result == 0x0)
    {
        base = 'A';
    }
    else if(result == 0x1)
    {
       base = 'G';
    }
    else if(result == 0x2)
    {
        base = 'C';
    }
    else if(result == 0x3)
    {
        base = 'T';
    }

    return base;
}

void setQuality(struct Sequence* sequence, int location, char quality)
{
    // Within range?
    if(0 <= location && location < sequence->length)
    {
        sequence->quality[location] = quality;
    }
}

void insertQuality(struct Sequence* sequence, int location, char quality)
{
    // Within range?
    if(0 <= location && location < sequence->length)
    {
        // Make room:
        sequence->quality = (char*) realloc (sequence->quality, strlen(sequence->quality) + 2);
        
        // Make a temporary string:
        char * temp = (char*) malloc (strlen(sequence->quality) + 1);
        strcpy(temp, &(sequence->quality[location]));   // Copy end bit in.
        
        // Push end bit over by one after the insert location:
        strcpy(&(sequence->quality[location + 1]), temp);
        
        // Make insertion:
        sequence->quality[location] = quality;
        
        // Remove temporary string:
        free(temp);
    }
}

void deleteQuality(struct Sequence* sequence, int location)
{
    // Within range?
    if(0 <= location && location < sequence->length)
    {
        // Make a temporary string:
        char * temp = (char*) malloc (strlen(sequence->quality) + 1);
        strcpy(temp, &(sequence->quality[location + 1]));   // Copy end bit in.
        
        // Remove:
        strcpy(&(sequence->quality[location]), temp);
        
        // Shrink:
        // THE STRING HAS ALREADY BEEN SHRUNK!
        // The terminating character will be in the correct place. 
        // Adjust the memory to match the new string length.
        // strlen(x) + 1 is the length plus 1 for the terminating character.
        sequence->quality = realloc (sequence->quality, strlen(sequence->quality) + 1);      
        
        // Remove temporary string:
        free(temp);
    }
}

void setBase(unsigned long long int* nucleotideSequence, 
        unsigned int nucleotidePosition, char nucleotide)
{
    unsigned long long int mask = 0x3;
    unsigned long long int base;

    unsigned int sequenceIndex = nucleotidePosition / 32;    
    unsigned int shift = (31 - (nucleotidePosition % 32)) * 2;
    
    if (nucleotide == 'A')
    {
        base = 0x0;
    }
    else if(nucleotide == 'G')
    {
        base = 0x1;
    }
    else if(nucleotide == 'C')
    {
        base = 0x2;
    }
    else if(nucleotide == 'T')
    {
        base = 0x3;
    }
    
    //white out the bits to be changed
    nucleotideSequence[sequenceIndex] &= ~(mask << shift);
    
    //write the new bits:
    nucleotideSequence[sequenceIndex] |= (base << shift);
}

void changeBase(struct Sequence* sequence, unsigned int nucleotidePosition, 
        char nucleotide, char quality)
{    
    setBase(sequence->sequence, nucleotidePosition, nucleotide);
    setQuality(sequence, nucleotidePosition, quality);
}

void deleteBase(struct Sequence* sequence, unsigned int nucleotidePosition)
{
    const unsigned long long int MASK = 0xFFFFFFFFFFFFFFFF;
    
    // Indices involved in deletion:
    unsigned int startIndex = nucleotidePosition / 32;
    unsigned int endNucleotidePosition = getMax(0, (sequence->length - 1));
    unsigned int endIndex = endNucleotidePosition / 32;
    
    unsigned long long int temp;
    char base;
    
    // Handle the block of memory where the deletion occurs:
    if(nucleotidePosition % 32 != 31)   
    {
        // Not needed if deleting last nucleotide in block.
        temp = sequence->sequence[startIndex] & (MASK >> ((nucleotidePosition * 2 % 64) + 2));
        sequence->sequence[startIndex] = sequence->sequence[startIndex] & ~(MASK >> (nucleotidePosition * 2 % 64));
        sequence->sequence[startIndex] = sequence->sequence[startIndex] | (temp << 2);
    }
    
    base = getBase(sequence->sequence, ((startIndex + 1) * 32));
    setBase(sequence->sequence, ((startIndex + 1) * 32 - 1), base);
    
    // Shift the rest of the memory to the left by two:
    for(int i = startIndex + 1; i < endIndex; i++)
    {
        sequence->sequence[i] = sequence->sequence[i] << 2;
        
        base = getBase(sequence->sequence, ((i + 1) * 32));
        setBase(sequence->sequence, ((i + 1) * 32 - 1), base);
    }
    
    if(startIndex != endIndex)  // Don't want to shift twice!
    {
        // Shift the last block of memory:
        sequence->sequence[endIndex] = sequence->sequence[endIndex] << 2;
    }
    
    // White out the last base:
    setBase(sequence->sequence, endNucleotidePosition, 'A');
    
    // Update the qualities:
    // NOTE: Needs to be done BEFORE length changes.
    deleteQuality(sequence, nucleotidePosition);
    
    // Update length:
    sequence->length = sequence->length - 1;
}

void insertBase(struct Sequence* sequence, unsigned int nucleotidePosition, 
        char nucleotide, char quality)
{
    const unsigned long long int MASK = 0xFFFFFFFFFFFFFFFF;
    char base;

    unsigned int sequenceIndex = nucleotidePosition / 32;
    
    int blocksNeeded = getNumMemoryBlocks(sequence->length + 1);  
        
    // Do we need to allocate another block?
    if(blocksNeeded > sequence->blocks)
    {
        sequence->sequence = realloc(sequence->sequence, (blocksNeeded * sizeof(unsigned long long int)));
        sequence->blocks = blocksNeeded;
    }
    
    // Push everything after insert block over once.    
    for(int i = sequence->blocks - 1; i > sequenceIndex; i--)
    {
        sequence->sequence[i] = sequence->sequence[i] >> 2;
        
        base = getBase(sequence->sequence, (i * 32 - 1));
        setBase(sequence->sequence, i * 32, base);
    }
    
    // Handle the block of memory where the insertion occurs:
    unsigned long long int temp = sequence->sequence[sequenceIndex] & (MASK >> (nucleotidePosition * 2 % 64));
    sequence->sequence[sequenceIndex] = sequence->sequence[sequenceIndex] & ~(MASK >> (nucleotidePosition * 2 % 64));
    sequence->sequence[sequenceIndex] = sequence->sequence[sequenceIndex] | (temp >> 2);
    
    // Set the base in the newly created gap.
    setBase(sequence->sequence, nucleotidePosition, nucleotide);
    
    // Update the sequence.
    sequence->length = sequence->length + 1;
    
    // Update the qualities:
    // NOTE: Needs to be done AFTER length changes.
    insertQuality(sequence, nucleotidePosition, quality);
}

int getHomopolymerLength(unsigned long long int* nucleotideSequence, 
        int sequenceLength, unsigned int nucleotidePosition)
{
    char base = getBase(nucleotideSequence, nucleotidePosition);
    
    int count = 1;
    
    // Left
    for(int i = nucleotidePosition - 1; i >= 0 && getBase(nucleotideSequence, i) == base; i--)
    {
        count++;
    }
    
    // Right
    for(int i = nucleotidePosition + 1; i < sequenceLength && getBase(nucleotideSequence, i) == base; i++)
    {
        count++;
    }
    
    return count;
}

int getHomopolymerLeftmostNucleotide(unsigned long long int* nucleotideSequence, 
        unsigned int nucleotidePosition)
{
    char base = getBase(nucleotideSequence, nucleotidePosition);
    
    int leftmostNucleotide = nucleotidePosition;
    
    // Left
    for(int i = nucleotidePosition - 1; i >= 0 && getBase(nucleotideSequence, i) == base; i--)
    {
        leftmostNucleotide--;
    }
    
    return leftmostNucleotide;
}

void setHomopolymerLength(struct Sequence* sequence, 
        unsigned int nucleotidePosition, unsigned int size,
        char quality)
{
    int originalSize = getHomopolymerLength(sequence->sequence, sequence->length, nucleotidePosition);
    int difference = originalSize - size;
    
    char base = getBase(sequence->sequence, nucleotidePosition);
    char current;
    int position;
    
    // Currently too long.
    if(difference > 0)
    {
        current = base;
        
        // Locate leftmost nucleotide in homopolymer.
        for(position = nucleotidePosition; 
            position >= 0 && getBase(sequence->sequence, position) == base;
            position--);
        
        // Shift back one if necessary.
        
        if(position < 0)
        {
            position = 0;
        }
        
        if(getBase(sequence->sequence, position) != base)
        {
            position++;
        }      
        
        // Delete bases.
        for(int i = difference; i > 0; i--)
        {
            deleteBase(sequence, position);
        }
    }
    // Currently too short.
    else if(difference < 0)
    {
        // Insert bases.
        for(int i = difference; i < 0; i++)
        {
            insertBase(sequence, nucleotidePosition, base, quality);
        }
    }
}
