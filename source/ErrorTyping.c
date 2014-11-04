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

#include "ErrorTyping.h"
#include "Counting.h"

unsigned int isJump(unsigned int value1 , unsigned int value2)
{
    const unsigned int JUMP_VALUE_THRESHOLD = 3;
    const float JUMP_PERCENT_THRESHOLD = 0.2;
    
    // The bigger value should be value1
    if(value2 > value1)
    {
        return isJump(value2, value1);
    }

    unsigned int difference = value1 - value2;   // positive

    // Is the difference less than a small percentage of the larger value?
    if(difference <= (value1 * JUMP_PERCENT_THRESHOLD))
    {
        return 0;       // Not a jump.
    }
    // Is the difference less than a small numeric threshold?
    else if(difference <= JUMP_VALUE_THRESHOLD)
    {
        return 0;       // Not a jump.
    }
    else
    {
        return 1;       // A jump.
    }

}

int getExternalErrorPosition(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    unsigned int total = length - kmerSize + 1;
    
    // Data structures:
    unsigned int counts[total];
    getKMerCounts(sequence, length, kmers, kmerSize, counts);
    
    // Beginning?
    if(counts[0] <= THRESHOLD)
    {
        // Find a jump from low to high:
        for(int i = 1; i < kmerSize; i++)
        {
            //Is there a low to high jump?
            //(CHECK JUMP FIRST!)
            if(isJump(counts[i], counts[i - 1]) 
                && counts[i - 1] <= THRESHOLD
                && counts[i] > THRESHOLD)
            {
                // Return position:
                return (i - 1);
            }
            //Is it still below the threshold?
            //(CHECK AFTER JUMP CHECK!)
            else if(counts[i] > THRESHOLD)
            {
                // Went above threshold.
                break;
            }
        }
    }
    
    // End?
    if(counts[total - 1] <= THRESHOLD)
    {
        // Find a jump from low to high:
        // (working backwards)
        for(int i = total - 2; i >= total - kmerSize; i--)
        {
            //Is there a low to high jump?
            //(CHECK JUMP FIRST!)
            if(isJump(counts[i], counts[i + 1]) 
                && counts[i + 1] <= THRESHOLD
                && counts[i] > THRESHOLD)
            {
                // Return position:
                return (i + kmerSize);
            }
            //Is it still below the threshold?
            //(CHECK AFTER JUMP CHECK!)
            else if(counts[i] > THRESHOLD)
            {
                // Went above threshold.
                break;
            }
        }
    }
    
    return -1;
}

unsigned int typeExternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Analyze:    
    // Non-negative:
    if(getExternalErrorPosition(sequence, length, kmers, kmerSize, THRESHOLD) >= 0)
    {
        return 1;
    }
    
    return 0;  
}

int getInternalErrorPosition(unsigned long long int* sequence, unsigned int sequenceLength, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    int result = -1;
    int length = -1;
    int start = - 1;
    
    length = getLengthOfInternalError(sequence, sequenceLength, kmers, kmerSize, THRESHOLD);
    
    // Is there a reasonably sized internal error?
    if(length >= (kmerSize / 2))
    {
        // Where does it start?
        start = getStartOfInternalError(sequence, sequenceLength, kmers, kmerSize, THRESHOLD);
    
        if(start > 0)
        {
            result = start + kmerSize - 1;      // Pinpoint error location.
        }
    }
    
    return result;
}

int getStartOfInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    unsigned int total = length - kmerSize + 1;
    int start = -1;
    
    // Data structures:
    unsigned int counts[total];
    getKMerCounts(sequence, length, kmers, kmerSize, counts);
    
    // Analyze:
    // Find start:
    for(int i = 1; i < (int)total - (int)kmerSize; i++)
    {
        // Find the first "drop" from high to low:
        // The value before the discontinuity should be above threshold.
        // The value after the discontinuity should be below threshold.
        if(isJump(counts[i], counts[i - 1]) 
                && counts[i - 1] > THRESHOLD
                && counts[i] <= THRESHOLD)
        {
            start = i;
            break;
        }
    }
    
    return start;
}

int getLengthOfInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    unsigned int total = length - kmerSize + 1;
    int start = -1;
    int end = -1;
    int lengthOfError = -1;
    
    // Data structures:
    unsigned int counts[total];
    getKMerCounts(sequence, length, kmers, kmerSize, counts);
    
    // Analyze:
    // Find start:
    start = getStartOfInternalError(sequence, length, kmers, kmerSize, THRESHOLD);
    
    if(start > 0)
    {
        // Find end:
        for(int i = start; i < total - 1 && counts[i] <= THRESHOLD; i++)
        {        
            // Does it jump back up to above the threshold?
            if(isJump(counts[i], counts[i + 1]) 
                && counts[i] <= THRESHOLD
                && counts[i + 1] > THRESHOLD)
            {       
                end = i;
                break;
            }
        }
    }
    
    if(start > 0 && end > 0)
    {
        lengthOfError = end - start + 1;
    }
    
    return lengthOfError;
}

int getEndOfInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    int start = getStartOfInternalError(sequence, length, kmers, kmerSize, THRESHOLD);
    int lengthOfError = getLengthOfInternalError(sequence, length, kmers, kmerSize, THRESHOLD);
    int end = -1;
    
    if(start > 0 && lengthOfError > 0)
    {
        end = start + length;
    }
    
    return end;
}

unsigned int typeInternalError(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    int lengthOfError = getLengthOfInternalError(sequence, length, 
        kmers, kmerSize, THRESHOLD);
        
    // Analyze:    
    if(lengthOfError > 0)
    {
        return 1;
    }
    
    return 0;    
}

unsigned int typeHighQuality(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    unsigned int total = length - kmerSize + 1;
    
    // Data structures:
    unsigned int counts[total];
    getKMerCounts(sequence, length, kmers, kmerSize, counts);
    
    // Analyze:
    for(int i = 0; i < total; i++)
    {
        if(counts[0] <= THRESHOLD || (i > 0 && isJump(counts[i], counts[i - 1])))
        {
            return 0;
        }   
    }    
    
    return 1;    
}

unsigned int typeLowCoverage(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize, const unsigned int THRESHOLD)
{
    // Variables:
    unsigned int total = length - kmerSize + 1;
    
    int min = 0x7FFFFFFF;
    
    // Data structures:
    unsigned int counts[total];
    getKMerCounts(sequence, length, kmers, kmerSize, counts);
    
    // Analyze:
    for(int i = 0; i < total; i++)
    {
        if(counts[i] < min)
        {
            min = counts[i];
        }
        
        // Is the value less than the threshold?
        // Or is there a jump?
        if(i > 0 && isJump(counts[i], counts[i - 1]))
        {
            return 0;
        }
    }
    
    //No jumps..
    // Is there some point in the read with low coverage?
    if(min <= THRESHOLD)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

unsigned int typeHomopolymer(unsigned long long int* sequence, unsigned int length, 
        KMerHashTable* kmers,  unsigned int kmerSize)
{
    // Variables:
    unsigned int total = length - kmerSize + 1;
    char current = 0;
    char previous = getBase(sequence, 0);
    int homopolymerLength = 1;
    
    // Data structures:
    unsigned int counts[total];
    getKMerCounts(sequence, length, kmers, kmerSize, counts);
    
    // Analyze:
    for(int i = 1; i < length; i++)
    {
        current = getBase(sequence, i);
        
        // Continuing homopolymer
        if(current == previous)
        {
            homopolymerLength++;
        }
        // Ending homopolymer
        else
        {   
            // Have we found an error profile that matches homopolymers?
            if(homopolymerLength >= 4)
            {
                // Not located in the start?
                if(!(i - homopolymerLength < kmerSize))
                {
                    //printf("NOT IN START\n");
                    //printf("s = %d, e = %d, l = %d, b=%d, a=%d, c[b]=%d, c[a]=%d\n", i - homopolymerLength, i, homopolymerLength, (i - homopolymerLength) - kmerSize, i + 1 - kmerSize, counts[(i - homopolymerLength) - kmerSize], counts[i + 1 - kmerSize]);
                    //printf("homopolymer=%c\n", previous);
                    
                    if(counts[(i - homopolymerLength) - kmerSize] * 0.35 >= counts[i + 1 - kmerSize])
                    {
                        //printf("FLAGGED homopolymer=%c\n", previous);
                        return 1;
                    }
                }
                // Not located in the end?
                else if(!(i > total))
                {
                    //printf("NOT IN END\n");
                    //printf("s = %d, e = %d, l = %d, b=%d, a=%d, c[b]=%d, c[a]=%d\n", i - homopolymerLength, i, homopolymerLength, i - homopolymerLength - 1, i, counts[i - homopolymerLength - 1], counts[i]);
                    //printf("homopolymer=%c\n", previous);
                    
                    if(counts[i] * 0.35 >= counts[i - homopolymerLength - 1])
                    {
                        //printf("FLAGGED homopolymer=%c\n", previous);
                        return 1;
                    }
                }
            }
                
            previous = current;
            homopolymerLength = 1;
        }
    }
    
    return 0;  
}

/*
char typeSequence(unsigned long long int* sequence, unsigned int sequenceLength,
        KMerHashTable* kmers, unsigned int kmerSize)
{
    char result = UNKNOWN;
    
    if(sequenceLength > kmerSize)
    {
        if(typeHighQuality(sequence, sequenceLength, kmers, kmerSize, LOW_COVERAGE_THRESHOLD))
        {
            result = HIGH_QUALITY;
        }
        else if(typeLowCoverage(sequence, sequenceLength, kmers, kmerSize, LOW_COVERAGE_THRESHOLD))
        {
            result = LOW_COVERAGE;
        }
        else if(typeInternalError(sequence, sequenceLength, kmers, kmerSize, 1))
        {
            result = INTERNAL_ERROR;
        }
        else if(typeExternalError(sequence, sequenceLength, kmers, kmerSize, 1))
        {
            result = EXTERNAL_ERROR;
        }
    }
    
    return result;
}
*/

void printDistributionOfTypes(Reads** reads, unsigned int numInputFiles)
{
    // Counts:
    unsigned int highQuality = 0;
    unsigned int lowCoverage = 0;
    unsigned int internal = 0;
    unsigned int external = 0;
    unsigned int corrected = 0;
    unsigned int unknown = 0;
    
    // Iterate over all files:
    for(int file = 0; file < numInputFiles; file++)
    {    
        // Iterate over all reads:
        readsReset(reads[file]);
        
        while(readsHasNext(reads[file]))
        {   
            struct read* current = readsGetNext(reads[file]); 
            
            switch(current->type)
            {
                case HIGH_QUALITY: highQuality++; break;
                case LOW_COVERAGE: lowCoverage++; break;
                case UNKNOWN: unknown++; break;
                case INTERNAL_ERROR: internal++; break;
                case EXTERNAL_ERROR: external++; break;
                case CORRECTED: corrected++; break;
            }
        }
    }
    
    printf("High Quality: %d\n", highQuality);
    printf("Low Coverage: %d\n", lowCoverage);
    printf("Unknown: %d\n", unknown);
    printf("Internal: %d\n", internal);
    printf("External: %d\n", external);
    printf("Corrected: %d\n", corrected);
}
