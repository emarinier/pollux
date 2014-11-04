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

#include "ErrorCorrection.h"
#include "Utility.h"
#include "ErrorTyping.h"
#include "Counting.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "Reads.h"
#include "Correction.h"
#include "ErrorOutput.h"

typedef enum 
{
    SUB_A = 0,
    SUB_T,
    SUB_C,
    SUB_G,
    INS,
    DEL_L_A,
    DEL_L_T,
    DEL_L_C,
    DEL_L_G,
    DEL_R_A,
    DEL_R_T,
    DEL_R_C,
    DEL_R_G,
    CorrectionType_Count       // Needs to be last! Count of items in enum!
} CorrectionType;

unsigned int substitutionErrors = 0;
unsigned int insertionErrors = 0;
unsigned int deletionErrors = 0;
unsigned int multipleErrors = 0;
unsigned int homopolymerErrors = 0;

#define HOMOPOLYMER_SIZE_RANGE 10        // [< -RANGE, -RANGE, .. , RANGE, > RANGE]
#define HOMOPOLYMER_SIZE_OUTSIDE_NEGATIVE 0     // < -RANGE
#define HOMOPOLYMER_SIZE_NEGATIVE_START 1       // -RANGE
#define HOMOPOLYMER_SIZE_NEGATIVE_END (HOMOPOLYMER_SIZE_NEGATIVE_START + HOMOPOLYMER_SIZE_RANGE - 1) // (-1)
#define HOMOPOLYMER_SIZE_POSITIVE_START (HOMOPOLYMER_SIZE_NEGATIVE_END + 1)     // (1))
#define HOMOPOLYMER_SIZE_POSITIVE_END (HOMOPOLYMER_SIZE_POSITIVE_START + HOMOPOLYMER_SIZE_RANGE - 1)    // + RANGE
#define HOMOPOLYMER_SIZE_OUTSIDE_POSITIVE (HOMOPOLYMER_SIZE_POSITIVE_END + 1)   // > + RANGE

#define MINIMUM_HOMOPOLYMER_SIZE 1

unsigned int homopolymerSize[HOMOPOLYMER_SIZE_OUTSIDE_POSITIVE + 1] = {0};

void printCorrectionResults()
{
    unsigned int homopolymerInsertions = 0;
    unsigned int homopolymerDeletions = 0;
    
    int currentSize;
    
    printf("\n");
    
    printf("Corrected...\n");
    printf("Reads with Multiple Corrections: %d\n", multipleErrors);
    printf("\n");
    printf("Substitution Corrections: %d\n", substitutionErrors);
    printf("Single Insertion Corrections: %d\n", insertionErrors);
    printf("Single Deletion Corrections: %d\n", deletionErrors);
    printf("\n");
    printf("Total Homopolymer Corrections: %d\n", homopolymerErrors);
    printf("\n");
    printf("Breakdown:\n");
    
    printf("< -%d: %d\n", HOMOPOLYMER_SIZE_RANGE, homopolymerSize[HOMOPOLYMER_SIZE_OUTSIDE_NEGATIVE]);
            
    for(int i = HOMOPOLYMER_SIZE_NEGATIVE_START; i <= HOMOPOLYMER_SIZE_NEGATIVE_END; i++)
    {
        currentSize = (i - HOMOPOLYMER_SIZE_RANGE - 1);
        homopolymerDeletions += homopolymerSize[i] * -currentSize;
        printf("%d : %d\n", currentSize, homopolymerSize[i]);
    }
    
    for(int i = HOMOPOLYMER_SIZE_POSITIVE_START; i <= HOMOPOLYMER_SIZE_POSITIVE_END; i++)
    {
        currentSize = (i - HOMOPOLYMER_SIZE_RANGE);
        homopolymerInsertions += homopolymerSize[i] * currentSize;
        printf("%d : %d\n", currentSize, homopolymerSize[i]);
    }
    
    printf("> %d: %d\n", HOMOPOLYMER_SIZE_RANGE, homopolymerSize[HOMOPOLYMER_SIZE_OUTSIDE_POSITIVE]);
    
    printf("\n");
    printf("Total Insertions (Single + Homopolymer): %d\n", insertionErrors + homopolymerInsertions);
    printf("Total Deletions (Single + Homopolymer): %d\n", deletionErrors + homopolymerDeletions);
}

bool isHighToLow(unsigned int* kmers, unsigned int kmerLocation)
{
    return (kmers[kmerLocation] > kmers[kmerLocation + 1]);
}

int getSequenceLocation(struct Sequence* sequence, int kmerLocation, 
        KMerHashTable* kmers, unsigned int kmerSize)
{
    // K-mers.
    int total = sequence->length - kmerSize + 1;
    unsigned int kmerCounts[total];
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmerCounts);
    
    int sequenceLocation;
    
    // Order of discrepancy?
    if(isHighToLow(kmerCounts, kmerLocation))
    {
        // High -> Low
        sequenceLocation = kmerLocation + kmerSize;
    }
    else
    {
        // Low -> High
        sequenceLocation = kmerLocation;
    }
    
    return sequenceLocation;
}

char getAverageQuality(struct Sequence* sequence, int first, int second)
{
    const char ERROR_VALUE = 33;
    
    int total = 0;
    int count = 0;
    
    // Within range:
    if(0 <= first  && first < sequence->length)
    {
        total += (int)sequence->quality[first];
        count++;
    }
    
    // Within range:
    if(0 <= second && second < sequence->length)
    {
        total += (int)sequence->quality[second];
        count++;
    }
    
    if(count != 0)
    {
        return (char)(total / count);
    }
    else
    {
        return ERROR_VALUE;
    }
}

int evaluateCorrection(bool highToLow, unsigned int* kmerCounts, int total, int kmerLocation)
{
    int count = 0;
    
    // HIGH -> LOW
    if(highToLow)
    {
        for(int i = kmerLocation; i < total - 1; i++)
        {
            if(isJump(kmerCounts[i], kmerCounts[i + 1]))
            {
                break;
            }
            else
            {
                count++;
            }
        }
    }
    // LOW -> HIGH
    else
    {
        for(int i = kmerLocation; i >= 0; i--)
        {
            if(isJump(kmerCounts[i], kmerCounts[i + 1]))
            {
                break;
            }
            else
            {
                count++;
            }
        }
    }
    
    return count;
}

bool doSubstitutionCorrection(
        struct Sequence* sequence, unsigned int kmerLocation, char substitution,
        KMerHashTable* kmers, unsigned int kmerSize)
{
    int sequenceLocation = getSequenceLocation(sequence, kmerLocation, kmers, kmerSize);    
    char originalBase = getBase(sequence->sequence, sequenceLocation);
    char originalQuality = sequence->quality[sequenceLocation];
    char newQuality = getAverageQuality(sequence, sequenceLocation - 1, sequenceLocation + 1);
    
    // K-mers:
    int total = sequence->length - kmerSize + 1;
    unsigned int kmersInitial[total];
    unsigned int kmersCorrection[total];
    
    // Get initial counts.
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmersInitial);
    
    // Try correction.   
    changeBase(sequence, sequenceLocation, substitution, newQuality);
    
    // Get new counts.
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmersCorrection);
    
    // Note changes:
    int kmersImproved = evaluateCorrection(isHighToLow(kmersInitial, kmerLocation), 
            kmersCorrection, sequence->length - kmerSize + 1, kmerLocation);
    
    // Revert:
    changeBase(sequence, sequenceLocation, originalBase, originalQuality);
    
    return kmersImproved;
}

// Deletion Error: Need to insert a base.
bool doDeletionCorrection(struct Sequence* sequence, unsigned int kmerLocation, 
        char nucleotide, bool insertLeft,
        KMerHashTable* kmers, unsigned int kmerSize)
{    
    int sequenceLocation = getSequenceLocation(sequence, kmerLocation, kmers, kmerSize);
    char newQuality = getAverageQuality(sequence, sequenceLocation - 1, sequenceLocation + 1);
    
    // K-mers.
    int total = sequence->length - kmerSize + 1;
    unsigned int kmersInitial[total];
    unsigned int kmersCorrection[total];
    
    // Get initial counts.
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmersInitial);
    
    // Try correction.
    if(insertLeft)
    {
        insertBase(sequence, sequenceLocation, nucleotide, newQuality);
    }
    else
    {
        insertBase(sequence, sequenceLocation + 1, nucleotide, newQuality);
    }
    
    // Get new counts.
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmersCorrection);
    
    bool highToLow = isHighToLow(kmersInitial, kmerLocation);
    int kmersImproved;
    
    // Note changes:
    if(highToLow)
    {
        kmersImproved = evaluateCorrection(highToLow, kmersCorrection, sequence->length - kmerSize + 1, kmerLocation);
    }
    else
    {
        kmersImproved = evaluateCorrection(highToLow, kmersCorrection + 1, sequence->length - kmerSize + 1, kmerLocation);
            // (kmersCorrection + 1) because the 'good' bases has moved right
            // we don't want to look at garbage and think it's gold
    }
    
    // Revert:
    if(insertLeft)
    {
        deleteBase(sequence, sequenceLocation);   
    }
    else
    {
        deleteBase(sequence, sequenceLocation + 1);   
    }
    
    return kmersImproved - 1;
}

bool correctErrorAtLocationInsertion(struct Sequence* sequence,
        KMerHashTable* kmers, unsigned int kmerSize, int kmerLocation)
{   
    // Safety:
    if(sequence->length <= kmerSize + 1)
    {
        // We will be making the read even shorter.
        return false;
    }
 
    // Sequence:
    int sequenceLocation = getSequenceLocation(sequence, kmerLocation, kmers, kmerSize);
    char originalBase = getBase(sequence->sequence, sequenceLocation);
    char originalQuality = sequence->quality[sequenceLocation];
    
    // K-mers.
    int total = sequence->length - kmerSize + 1;
    unsigned int kmersInitial[total];   // Original
    unsigned int kmersCorrection[total - 1];    // Original with deletion. (-1)
    
    // Get initial counts.
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmersInitial);
    
    // Delete:
    deleteBase(sequence, sequenceLocation);    
    
    // Get new counts.
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmersCorrection);
    
    bool highToLow = isHighToLow(kmersInitial, kmerLocation);
    int kmersImproved;
    
    // Note changes:
    if(highToLow)
    {
        kmersImproved = evaluateCorrection(highToLow, kmersCorrection, sequence->length - kmerSize + 1, kmerLocation);
    }
    else
    {
        kmersImproved = evaluateCorrection(highToLow, kmersCorrection - 1, sequence->length - kmerSize + 1, kmerLocation);
            // (kmersCorrection - 1) because 'good' bases have fallen into the position
    }
    
    // Revert:
    insertBase(sequence, sequenceLocation, originalBase, originalQuality);
    
    return kmersImproved;
}

int getNextKMerDiscrepancy(double* kmerDiscrepancies, int total)
{
    double score = -1;
    int kmerLocation = -1;
    
    for(int i = 0; i < total; i++)
    {
        if(kmerDiscrepancies[i] > 0 && kmerDiscrepancies[i] > score)
        {
            score = kmerDiscrepancies[i];
            kmerLocation = i;
        }
    }
    
    return kmerLocation;
}

double scoreKMerDiscrepancy(int value1, int value2)
{
    int high = getMax(value1, value2);
    int low = getMin(value1, value2);
    
    high = getMax(high, 1);
    low = getMax(low, 1);
    
    //double score = ((double)low * (double)low) / (double)high;
    double score = (double)high - (double)low;
    
    return score;
}

void findKMerDiscrepancies(struct Sequence* sequence, 
        KMerHashTable* kmers, unsigned int kmerSize, double* kmerDiscrepancies)
{    
    // K-mers.
    int total = sequence->length - kmerSize + 1;
    unsigned int kmerCounts[total];
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmerCounts);
    
    // Analyze possible discrepancies
    for(int i = 0; i < total - 1; i++)
    {
        // Is there something we would consider a discrepancy?
        if(isJump(kmerCounts[i], kmerCounts[i + 1]))
        {
            kmerDiscrepancies[i] = scoreKMerDiscrepancy(kmerCounts[i], kmerCounts[i + 1]);
        }
        // Not a discrepancy.
        else
        {
            kmerDiscrepancies[i] = -1;
        }
    }
}

double getAverageKMerCount(struct Sequence* sequence,
        KMerHashTable* kmers, unsigned int kmerSize, int start, int end)
{    
    double average = 0;
    
    // K-mers:
    int total = sequence->length - kmerSize + 1;
    unsigned int kmerCounts[total];
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmerCounts);
    
    // Averaged k-mers:
    for(int i = start; i < end; i++)
    {
        average += kmerCounts[i];
    }
    
    average = (double)average / (double)(end - start);
    
    return average;
}

bool correctErrorAtLocationHomopolymer(struct Sequence* sequence,
        KMerHashTable* kmers, unsigned int kmerSize, int kmerLocation)
{
    bool highToLow;
    int initialHomopolymerLength;
    int currentHomopolymerLength;
    int homopolymerLeftmostNucleotide;
    double currentAverage;
    double maxAverage = 0;
    int bestLength;
    
    int start, end; // Used in averaging.
    
    // K-mers:
    int total = sequence->length - kmerSize + 1;
    unsigned int kmerCounts[total];
    getKMerCounts(sequence->sequence, sequence->length, kmers, kmerSize, kmerCounts);
    
    int sequenceLocation = getSequenceLocation(sequence, kmerLocation, kmers, kmerSize);
    char newQuality = getAverageQuality(sequence, sequenceLocation - 1, sequenceLocation + 1);  //TODO: THIS SHOULD BE TEMP
    
    // High -> Low
    if(isHighToLow(kmerCounts, kmerLocation))
    {   
        highToLow = true;
        homopolymerLeftmostNucleotide = getHomopolymerLeftmostNucleotide(sequence->sequence, sequenceLocation - 1);
        
        // Initialize Average (Single):
        start = (homopolymerLeftmostNucleotide - kmerSize + 1) + getHomopolymerLength(sequence->sequence, sequence->length, homopolymerLeftmostNucleotide);
        end = start + 1;
            
        maxAverage = getAverageKMerCount(sequence, kmers, kmerSize, start, end);
    }
    // Low -> High
    else
    {   
        highToLow = false;
        homopolymerLeftmostNucleotide = getHomopolymerLeftmostNucleotide(sequence->sequence, sequenceLocation + 1);
        
        // Initialize Average (Single):
        start = homopolymerLeftmostNucleotide - 1;
        end = start + 1;
            
        maxAverage = getAverageKMerCount(sequence, kmers, kmerSize, start, end);
    }   
    
    // Initialize:
    initialHomopolymerLength = getHomopolymerLength(sequence->sequence, sequence->length, homopolymerLeftmostNucleotide);
    
    // Homopolymer must be at least minimum size to try operating on it.
    if (initialHomopolymerLength < MINIMUM_HOMOPOLYMER_SIZE)
    {
        return false;
    }
    
    currentHomopolymerLength = initialHomopolymerLength;
    bestLength = initialHomopolymerLength;
    
    for (int i = getMax(initialHomopolymerLength / 2, 1); 
            i <= initialHomopolymerLength * 2 && currentHomopolymerLength < kmerSize; i++)
    {        
        setHomopolymerLength(sequence, homopolymerLeftmostNucleotide, i, newQuality);

        if(sequence->length <= kmerSize)
        {
            // The read is too short to gain any information.
            continue;
        }

        currentHomopolymerLength = getHomopolymerLength(sequence->sequence, sequence->length, homopolymerLeftmostNucleotide);
        
        // FIND K-MER RANGE:
        // Idea is to use as few k-mers as needed. 
        // Have the k-mers be mostly the better part of the read.
        // High -> Low
        if(highToLow)
        {
            start = (homopolymerLeftmostNucleotide - kmerSize + 1) + (currentHomopolymerLength);
            end = start + 2;
        }
        // Low -> High
        else
        {
            start = homopolymerLeftmostNucleotide - 2;
            end = start + 2;
        }
        
        if(start < 0 || end > sequence->length - kmerSize + 1)
        {
            break;      // We're outside the k-mer range and will always be.
        }    
        
        currentAverage = getAverageKMerCount(sequence, kmers, kmerSize, start, end);            
        
        if(currentAverage > maxAverage)
        {
            maxAverage = currentAverage;
            bestLength = i;
        }

    }
    
    setHomopolymerLength(sequence, homopolymerLeftmostNucleotide, bestLength, newQuality);
    
    // Did we find something better?
    if(bestLength != initialHomopolymerLength)
    {        
        // Update information.
        sequence->corrections[sequence->numCorrections] = 'H';
        sequence->homopolymerSize[sequence->numCorrections] = bestLength - initialHomopolymerLength;
                
        sequence->numCorrections = sequence->numCorrections + 1;
        
        return true;
    }
    else
    {
        return false;
    }
}

void assignCorrection(CorrectionType correctionType, struct Sequence* sequence, 
        unsigned int kmerLocation, KMerHashTable* kmers, unsigned int kmerSize)
{
    int sequenceLocation = getSequenceLocation(sequence, kmerLocation, kmers, kmerSize);    
    char newQuality = getAverageQuality(sequence, sequenceLocation - 1, sequenceLocation + 1);
    
    switch(correctionType)
    {
        case SUB_A:
            changeBase(sequence, sequenceLocation, 'A', newQuality);
            break;
            
        case SUB_T:
            changeBase(sequence, sequenceLocation, 'T', newQuality);
            break;
            
        case SUB_C:
            changeBase(sequence, sequenceLocation, 'C', newQuality);
            break;
            
        case SUB_G:
            changeBase(sequence, sequenceLocation, 'G', newQuality);
            break;
            
        case INS:
            deleteBase(sequence, sequenceLocation);  
            break;
            
        case DEL_L_A:
            insertBase(sequence, sequenceLocation, 'A', newQuality);
            break;
            
        case DEL_L_T:
            insertBase(sequence, sequenceLocation, 'T', newQuality);
            break;
            
        case DEL_L_C:
            insertBase(sequence, sequenceLocation, 'C', newQuality);
            break;
            
        case DEL_L_G:
            insertBase(sequence, sequenceLocation, 'G', newQuality);
            break;
            
        case DEL_R_A:
            insertBase(sequence, sequenceLocation + 1, 'A', newQuality);
            break;
            
        case DEL_R_T:
            insertBase(sequence, sequenceLocation + 1, 'T', newQuality);
            break;
            
        case DEL_R_C:
            insertBase(sequence, sequenceLocation + 1, 'C', newQuality);
            break;
            
        case DEL_R_G:
            insertBase(sequence, sequenceLocation + 1, 'G', newQuality);
            break;
            
        default:
            printf("UNRECOGNIZED ERROR CORRECTION TYPE!\n");
            exit(0);
    }
    
    if(correctionType == SUB_A || correctionType == SUB_T 
            || correctionType == SUB_C || correctionType == SUB_G)
    {
        sequence->corrections[sequence->numCorrections] = 'S';
    }
    else if(correctionType == INS)
    {
        sequence->corrections[sequence->numCorrections] = 'I';
    }
    else if(correctionType == DEL_L_A || correctionType == DEL_L_T 
            || correctionType == DEL_L_C || correctionType == DEL_L_G
            || correctionType == DEL_R_A || correctionType == DEL_R_T
            || correctionType == DEL_R_C || correctionType == DEL_R_G)
    {
        sequence->corrections[sequence->numCorrections] = 'D';
    }
    else
    {
        printf("UNRECOGNIZED ERROR CORRECTION TYPE!\n");
        exit(0);
    }
    
    sequence->numCorrections = sequence->numCorrections + 1;
}

bool correctErrorAtLocation(struct Sequence* sequence, int kmerLocation, 
        Correction* correction)
{
    // KMers:
    KMerHashTable* kmers = correctionGetKMers(correction);
    unsigned int kmerSize = correctionGetKMerSize(correction);            
    int total = sequence->length - kmerSize + 1;
    
    // Safety:
    if (kmerLocation <= 0 || kmerLocation >= (total - 2))
    {
        // Out of bounds.
        return false;
    }   
    
    int correctionAttempt[CorrectionType_Count];
    
    // Initialize:
    for(int i = 0; i < CorrectionType_Count; i++)
    {
        correctionAttempt[i] = 0;
    }
    
    // Substitution Error
    if(correction->substitutions)
    {
        correctionAttempt[SUB_A] = doSubstitutionCorrection(sequence, kmerLocation, 'A', kmers, kmerSize);
        correctionAttempt[SUB_T] = doSubstitutionCorrection(sequence, kmerLocation, 'T', kmers, kmerSize);
        correctionAttempt[SUB_C] = doSubstitutionCorrection(sequence, kmerLocation, 'C', kmers, kmerSize);
        correctionAttempt[SUB_G] = doSubstitutionCorrection(sequence, kmerLocation, 'G', kmers, kmerSize);
    }
    
    // Insertion Error
    if(correction->insertions)
    {
        correctionAttempt[INS] = correctErrorAtLocationInsertion(sequence, kmers, kmerSize, kmerLocation);
    }
    
    // Deletion Error
    if(correction->deletions)
    {
        correctionAttempt[DEL_L_A] = doDeletionCorrection(sequence, kmerLocation, 'A', true, kmers, kmerSize);
        correctionAttempt[DEL_L_T] = doDeletionCorrection(sequence, kmerLocation, 'T', true, kmers, kmerSize);
        correctionAttempt[DEL_L_C] = doDeletionCorrection(sequence, kmerLocation, 'C', true, kmers, kmerSize);
        correctionAttempt[DEL_L_G] = doDeletionCorrection(sequence, kmerLocation, 'G', true, kmers, kmerSize);
        correctionAttempt[DEL_R_A] = doDeletionCorrection(sequence, kmerLocation, 'A', false, kmers, kmerSize);
        correctionAttempt[DEL_R_T] = doDeletionCorrection(sequence, kmerLocation, 'T', false, kmers, kmerSize);
        correctionAttempt[DEL_R_C] = doDeletionCorrection(sequence, kmerLocation, 'C', false, kmers, kmerSize);
        correctionAttempt[DEL_R_G] = doDeletionCorrection(sequence, kmerLocation, 'G', false, kmers, kmerSize);
    }
    
    // Find the best correction:
    CorrectionType bestCorrectionIndex = 0;
    int bestCorrectionValue = 0;
    
    for(int i = 0; i < CorrectionType_Count; i++)
    {
        if(correctionAttempt[i] > bestCorrectionValue)
        {
            bestCorrectionIndex = i;
            bestCorrectionValue = correctionAttempt[i];
        }
    }
    
    if(bestCorrectionValue >= 2)
    {
        assignCorrection(bestCorrectionIndex, sequence, kmerLocation, kmers, kmerSize);
        return true;
    }
    
    // Homopolymers
    if(correction->homopolymers
            && correctErrorAtLocationHomopolymer(sequence, kmers, kmerSize, kmerLocation))
    {
        return true;
    }
            
    return false;
}

void recordHomopolymerSize(int size)
{
    if(size > HOMOPOLYMER_SIZE_RANGE)
    {
        homopolymerSize[HOMOPOLYMER_SIZE_RANGE * 2 + 1]++;
    }
    else if(size < -HOMOPOLYMER_SIZE_RANGE)
    {
        homopolymerSize[0]++;
    }
    else if(size > 0)
    {
        homopolymerSize[HOMOPOLYMER_SIZE_RANGE + size]++;
    }
    else if(size < 0)
    {
        homopolymerSize[(HOMOPOLYMER_SIZE_RANGE + 1) + size]++;
    }
    else
    {
        printf("Homopolymer size recording error!");
    }
}

void applyCorrection(struct Sequence* sequence, struct read* read)
{
    free(read->sequence);
    read->sequence = sequence->sequence;      // NOT QUITE CORRECT!

    free(read->quality);
    read->quality = sequence->quality;

    read->length = sequence->length;

    // Update information.
    for(int i = 0; i < sequence->numCorrections; i++)
    {
        switch(sequence->corrections[i])
        {
            case 'S': substitutionErrors++; break;
            case 'I': insertionErrors++; break;
            case 'D': deletionErrors++; break;
            
            // Homopolymers:
            case 'H':  homopolymerErrors++;
                       recordHomopolymerSize(sequence->homopolymerSize[i]);
                       break;
        }
    }

    if(sequence->numCorrections > 1)
    {
        multipleErrors++;
    }
    
    free(sequence->corrections);
    free(sequence->homopolymerSize);
}

void revertCorrection(struct Sequence* sequence, struct read* read)
{
    free(sequence->sequence);
    free(sequence->quality);
    free(sequence->corrections);
    free(sequence->homopolymerSize);
}

void initializeSequence(struct Sequence* sequence, struct read* read, 
        Correction* correction, int maxCorrections)
{
    // KMers:
    KMerHashTable* kmers = correctionGetKMers(correction);
    unsigned int kmerSize = correctionGetKMerSize(correction);
    
    // Creating temporary date structure.   
    sequence->length = read->length;
    sequence->blocks = getNumMemoryBlocks(sequence->length);
    
    // Copy sequence:
    sequence->sequence = (unsigned long long int*) malloc(sequence->blocks * sizeof(unsigned long long int));
    memcpy(sequence->sequence, read->sequence, sequence->blocks * sizeof(unsigned long long int));
    
    // Copy quality scores:
    sequence->quality = (char*) malloc(strlen(read->quality) + 1);
    strcpy(sequence->quality, read->quality);
    
    // Correction information:
    sequence->numCorrections = 0;
    
    // Recorded corrections:
    sequence->corrections = (char*) malloc(maxCorrections * sizeof(char));
    sequence->homopolymerSize = (int*) malloc(maxCorrections * sizeof(int));
        // The size associated with the homopolymer correction record.
}

bool isHighQuality(struct read* read, Correction* correction)
{
    // K-mers:
    KMerHashTable* kmers = correctionGetKMers(correction);
    unsigned int kmerSize = correctionGetKMerSize(correction);
    
    int total = read->length - kmerSize + 1;
    double numUniqueKMers = 0;
    
    unsigned int kmerCounts[total];
    getKMerCounts(read->sequence, read->length, kmers, kmerSize, kmerCounts);
    
    // Scan k-mers:
    for(int i = 0; i < total; i++)
    {
        if(kmerCounts[i] == 1)
        {
            numUniqueKMers++;
        }
    }
    
    bool result;
    double fraction = (double)numUniqueKMers / (double)total;
    
    if (fraction < 0.5)
    {
        result = true;
    }
    else
    {
        result = false;
    }
    
    return result;
}

bool correctRead(struct read* read, Correction* correction)
{
    // KMers:
    KMerHashTable* kmers = correctionGetKMers(correction);
    unsigned int kmerSize = correctionGetKMerSize(correction);
    int maxCorrections = getMax(MAX_CORRECTIONS, (int)(0.2 * read->length));
    
    struct Sequence sequence;
    
    if(read->length <= kmerSize)
    {
        return false;
    }
    
    initializeSequence(&sequence, read, correction, maxCorrections);
    
    //-----------------------------//
    
    // Get initial discrepancies:
    int total = sequence.length - kmerSize + 1;    
    double* discrepancies = (double*)malloc((total - 1) * sizeof(double*));    
    findKMerDiscrepancies(&sequence, kmers, kmerSize, discrepancies);
    
    int kmerLocation = getNextKMerDiscrepancy(discrepancies, total - 1);
    
    // No discrepancies?
    if(kmerLocation < 0)
    {
        // Read is high quality or low coverage.
        if(typeHighQuality(sequence.sequence, sequence.length, kmers, kmerSize, correction->lowKMerThreshold))
        {
            sequence.type = HIGH_QUALITY;
        }
        else
        {
            sequence.type = LOW_COVERAGE;
        }
        
        // Differentiate between low coverage and bad.
        if(!isHighQuality(read, correction))
        {
            read->type = BAD;
        }
        
        free(sequence.sequence);
        free(sequence.quality);
        free(discrepancies);
        
        return true;
    }
    
    // While there is a location of something to correct:
    while(kmerLocation >= 0 && kmerLocation < sequence.length)
    {        
        // Safety -- In case we're stuck in a correction loop or read is terrible.
        if(sequence.numCorrections >= maxCorrections)
        {
            break;
        }
        
        if(correctErrorAtLocation(&sequence, kmerLocation, correction))
        {
            // Correction successful!
            // Get new discrepancies:
            total = sequence.length - kmerSize + 1;
            
            free(discrepancies);
            discrepancies = (double*)malloc((total - 1) * sizeof(double*));            
            findKMerDiscrepancies(&sequence, kmers, kmerSize, discrepancies); 
        }
        else
        {
            // Correction unsuccessful.
            discrepancies[kmerLocation] = -1;     // Move to next discrepancy.
        }
        
        kmerLocation = getNextKMerDiscrepancy(discrepancies, total - 1);
    }
 
    // Were we successful?
    if(sequence.numCorrections < MAX_CORRECTIONS)
    {
        // Correction worked.
        applyCorrection(&sequence, read);
    }
    else
    {
        // Couldn't correct.
        revertCorrection(&sequence, read);
    }
    
    free(discrepancies);
    
    // Check to see if the read is considered bad.
    if(isHighQuality(read, correction))
    {
        read->type = CORRECTED;
    }
    else
    {
        read->type = BAD;
    }
    
    return true;
}


