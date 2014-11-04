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
#include <stdlib.h>
#include <string.h>
#include "Globals.h"
#include "Utility.h"
#include "ErrorProcessing.h"
#include "ErrorOutput.h"
#include "ErrorTyping.h"
#include "ErrorCorrection.h"
#include "Reads.h"
#include "Correction.h"
#include <libgen.h>

unsigned int KMER_SIZE = 31;

const int LEFT = 0;
const int RIGHT = 1;

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

void hashSequence(unsigned long long int* sequence, unsigned int sequenceLength,
        KMerHashTable* kmers, unsigned int kmerSize)
{
    unsigned long long int* reverse;
    
    addKMersToTable(kmers, sequence, sequenceLength, kmerSize);

    reverse = createReverseCompliment(sequence, sequenceLength);       
    addKMersToTable(kmers, reverse, sequenceLength, kmerSize);
    
    free(reverse);
}

void hashReads(Correction* correction)
{
    // Reads:
    Reads** reads = correctionGetReads(correction);
    unsigned int numReadSets = correctionGetNumReadSets(correction);
    struct read* current;
    
    // Read:
    unsigned long long int* sequence;
    unsigned int length;
    
    // KMers:
    KMerHashTable* kmers = correctionGetKMers(correction);
    unsigned int kmerSize = correctionGetKMerSize(correction);
    
    // Iterate over all files:    
    for(int file = 0; file < numReadSets; file++)
    {
        printf("Processing file %d/%d...\n", file + 1, numReadSets);        
        
        readsReset(reads[file]);
        
        // Iterate over all reads:
        for(int i = 0; readsHasNext(reads[file]); i++)
        {            
            printProgress(i, readsGetCount(reads[file]), 20);
            
            current = readsGetNext(reads[file]);  
            sequence = current->sequence;
            length = current->length;

            hashSequence(sequence, length, kmers, kmerSize);
        }
        
        printf("\n");
        
        printf("Preprocessing k-mers...\n");
        preprocessKMers(kmers, correction);
        printf("Finished preprocessing k-mers!\n\n");
    }
}

void executePairedCorrection(Correction* correction,
        FILE* leftCorrectedFile, FILE* rightCorrectedFile, 
        FILE* leftGarbageFile, FILE* rightGarbageFile, FILE* extraFile)
{    
    // Reads:
    Reads** reads = correctionGetReads(correction);
    
    // Correction Function:
    CorrectionFunction correctionFunction = correctionGetFunction(correction);
    
    int maxReads = getMax(readsGetCount(reads[LEFT]), readsGetCount(reads[RIGHT]));
    int current = 0;
    
    printf("Correcting paired files.\n");

    readsReset(reads[LEFT]);
    readsReset(reads[RIGHT]);
    
    // Safety check:
    if(!(readsHasNext(reads[LEFT]) && readsHasNext(reads[RIGHT])))
    {
        printf("No reads to correct!\n");        
        return;
    }
    
    struct read* leftRead = readsGetNext(reads[LEFT]);
    struct read* rightRead = readsGetNext(reads[RIGHT]);
    
    // Always correct the read immediately after loading.
    correctionFunction(leftRead, correction);
    correctionFunction(rightRead, correction);
    
    while(readsHasNext(reads[LEFT]) && readsHasNext(reads[RIGHT]))
    {   
        printProgress(current, maxReads, 20);
        current++;

        // The current pair matches.
        // This is to accommodate read filtering which might have been done at
        //      the input level, which existed in the past.
        if(leftRead->number == rightRead->number)
        {
            // Garbage the pair of reads only if BOTH are bad.
            if (correction->filtering && leftRead->type == BAD && rightRead->type == BAD)
            {
               // write both to garbage pair files
                outputRead(leftGarbageFile, leftRead);
                outputRead(rightGarbageFile, rightRead);
            }
            else
            {
                // write both to corrected pair files
                outputRead(leftCorrectedFile, leftRead);
                outputRead(rightCorrectedFile, rightRead);
            }
            
            leftRead = readsGetNext(reads[LEFT]);
            rightRead = readsGetNext(reads[RIGHT]);
            
            correctionFunction(leftRead, correction);
            correctionFunction(rightRead, correction);
        }
        // The left read is ahead of the right.
        else if(leftRead->number > rightRead->number)
        {
            // write right to extra file
            outputRead(extraFile, rightRead);

            rightRead = readsGetNext(reads[RIGHT]);
            correctionFunction(rightRead, correction);
        }
        // The right read is ahead of the left.
        else if(leftRead->number < rightRead->number)
        {
            // write left to extra file
            outputRead(extraFile, leftRead);

            leftRead = readsGetNext(reads[LEFT]);
            correctionFunction(leftRead, correction);
        }
    }

    // Write leftover reads to extra:
    while(readsHasNext(reads[LEFT]))
    {
        outputRead(extraFile, leftRead);
        
        leftRead = readsGetNext(reads[LEFT]);
        correctionFunction(leftRead, correction);
    }

    while(readsHasNext(reads[RIGHT]))
    {
        outputRead(extraFile, rightRead);
        
        rightRead = readsGetNext(reads[RIGHT]);
        correctionFunction(rightRead, correction);
    }

    printf("\n");
    printCorrectionResults();
    printf("\n");   
}

int processSimpleCorrection(Correction* correction)
{
    // Reads:
    Reads** reads = correctionGetReads(correction);
    unsigned int numReadSets = correctionGetNumReadSets(correction);
    
    // Correction Function:
    CorrectionFunction correctionFunction = correctionGetFunction(correction);
    
    // Files:
    char correctedFileName[1024];
    char garbageFileName[1024];
    FILE* correctedFile;
    FILE* garbageFile;
    
    char* outputDirectory = correctionGetOutputDirectory(correction);
    char baseName[1024];
    
    // Iterate over all files:
    for(int file = 0; file < numReadSets; file++)
    {
        printf("Correcting file %d/%d...\n", file + 1, numReadSets);
        
        // Prepare output:
        strcpy(baseName, basename(readsGetFileName(reads[file])));
            // We don't want the original to get changed.
        
        strcpy(correctedFileName, outputDirectory);
        strcpy(garbageFileName, outputDirectory);
        
        if(correctedFileName[strlen(correctedFileName) - 1] != '/')
                strcat(correctedFileName, "/");
        
        if(garbageFileName[strlen(correctedFileName) - 1] != '/')
                strcat(garbageFileName, "/");
        
        strcat(correctedFileName, baseName);
        strcat(correctedFileName, ".corrected");
        correctedFile = fopen(correctedFileName, "w");
        
        if(correction->filtering)
        {
            strcat(garbageFileName, baseName);
            strcat(garbageFileName, ".low");
            garbageFile = fopen(garbageFileName, "w");
        }
        
        // Safety check:
        if(correctedFile == 0)
        {
            printf("Could not create file: %s\n", correctedFileName);
            return 1;
        }
        
        if(correction->filtering)
        {
            if(garbageFile == 0)
            {
                printf("Could not create file: %s\n", garbageFileName);
                return 1;
            }
        }
        
        readsReset(reads[file]);
        
        // Iterate over all reads:
        for(int i = 0; readsHasNext(reads[file]); i++)
        {            
            printProgress(i, readsGetCount(reads[file]), 20);
            
            struct read* current = readsGetNext(reads[file]);
            
            // Correction:
            correctionFunction(current, correction);
            
            // Output:
            if (current->type != BAD || correction->filtering == false)
            {
                outputRead(correctedFile, current);
            }
            else
            {
                outputRead(garbageFile, current);
            }             
        }
        
        // Close file:
        fclose(correctedFile);
        
        if(correction->filtering)
        {
            fclose(garbageFile);
        }
        
        printf("\n");
        printCorrectionResults();
        printf("\n");
    }     
}

void processPairedCorrection(Correction* correction)
{   
    char leftCorrectedFileName[1024];
    char rightCorrectedFileName[1024];
    char leftGarbageFileName[1024];
    char rightGarbageFileName[1024];    
    char extraFileName[1024];
    
    char* outputDirectory = correctionGetOutputDirectory(correction);
    char baseName[1024];

    FILE* leftCorrectedFile;
    FILE* rightCorrectedFile;
    FILE* leftGarbageFile;
    FILE* rightGarbageFile;
    FILE* extraFile;
    
    Reads** reads = correctionGetReads(correction);
    
    // Left Corrected:
    strcpy(baseName, basename(readsGetFileName(reads[LEFT])));
        // We don't want the original to get changed.
    
    strcpy(leftCorrectedFileName, outputDirectory);
    if(leftCorrectedFileName[strlen(leftCorrectedFileName) - 1] != '/')
                strcat(leftCorrectedFileName, "/");
    strcat(leftCorrectedFileName, baseName);
    strcat(leftCorrectedFileName, ".corrected");
    leftCorrectedFile = fopen(leftCorrectedFileName, "w");

    // Right Corrected:
    strcpy(baseName, basename(readsGetFileName(reads[RIGHT])));
        // We don't want the original to get changed.
    
    strcpy(rightCorrectedFileName, outputDirectory);
    if(rightCorrectedFileName[strlen(rightCorrectedFileName) - 1] != '/')
                strcat(rightCorrectedFileName, "/");
    strcat(rightCorrectedFileName, baseName);
    strcat(rightCorrectedFileName, ".corrected");
    rightCorrectedFile = fopen(rightCorrectedFileName, "w");
    
    // Left Garbage:
    if(correction->filtering)
    {
        strcpy(baseName, basename(readsGetFileName(reads[LEFT])));
                // We don't want the original to get changed.

        strcpy(leftGarbageFileName, outputDirectory);
        if(leftGarbageFileName[strlen(leftGarbageFileName) - 1] != '/')
                    strcat(leftGarbageFileName, "/");
        strcat(leftGarbageFileName, baseName);
        strcat(leftGarbageFileName, ".low");
        leftGarbageFile = fopen(leftGarbageFileName, "w");

        // Right Garbage:
        strcpy(baseName, basename(readsGetFileName(reads[RIGHT])));
            // We don't want the original to get changed.

        strcpy(rightGarbageFileName, outputDirectory);
        if(rightGarbageFileName[strlen(rightGarbageFileName) - 1] != '/')
                    strcat(rightGarbageFileName, "/");
        strcat(rightGarbageFileName, baseName);
        strcat(rightGarbageFileName, ".low");
        rightGarbageFile = fopen(rightGarbageFileName, "w");
    }

    // Extra:
    strcpy(extraFileName, outputDirectory);
    strcat(extraFileName, "/extra.corrected");
    extraFile = fopen(extraFileName, "w");

    executePairedCorrection(correction, leftCorrectedFile, rightCorrectedFile, 
            leftGarbageFile, rightGarbageFile, extraFile);

    // Close files:
    fclose(leftCorrectedFile);
    fclose(rightCorrectedFile);
    
    if(correction->filtering)
    {
        fclose(leftGarbageFile);
        fclose(rightGarbageFile);
    }
    
    fclose(extraFile);    
}

Correction* preprocessing(int numInputFiles, char* inputFileNames, char* outputDirectory) 
{
    unsigned int LOW_COVERAGE_THRESHOLD_DEFAULT = 3;

    // Data Structures:
    KMerHashTable* kmers = newKMerHashTable();
    Reads** reads = (Reads**)malloc(sizeof(Reads*) * numInputFiles);
    Correction* correction = (Correction*)malloc(sizeof(Correction));
    
    // CREATE READS:      
    printf("Creating read objects...\n");    
    for(int i = 0; i < numInputFiles; i++)
    {
        char* inputFileName = &(inputFileNames[i * 200]);
        printf("Reading file: %s\n", inputFileName);
        
        reads[i] = createReads(inputFileName);
    }    
    printf("Finished creating read objects!\n\n");    

    // META OBJECT:
    correction = createCorrection(reads, numInputFiles, 
            kmers, KMER_SIZE, LOW_COVERAGE_THRESHOLD_DEFAULT,
            outputDirectory, NULL);
    
    // CONSTRUCT KMERS:
    printf("Constructing k-mers...\n");       
    hashReads(correction);
    printf("Finished constructing k-mers!\n\n");
    
    return correction;
}

int convertFASTQToFASTK(int numInputFiles, char* inputFileNames, char* outputDirectory) 
{
    printf("FASTK CONVERSION\n\n");
    
    // META OBJECT:
    Correction* correction = preprocessing(numInputFiles, inputFileNames, outputDirectory);
    
    // # Reads:
    Reads** reads = correctionGetReads(correction);
    unsigned int numReadSets = correctionGetNumReadSets(correction);
    
    // KMers:
    KMerHashTable* kmers = correctionGetKMers(correction);
    unsigned int kmerSize = correctionGetKMerSize(correction);
    
    // Files:
    char convertedFileName[1024];
    FILE* convertedFile;  
    char baseName[1024];
    
    // Iterate over all files:
    for(int file = 0; file < numReadSets; file++)
    {        
        printf("Converting file %d/%d...\n", file + 1, numReadSets);
        
        // Prepare output:
        strcpy(baseName, basename(readsGetFileName(reads[file])));
            // We don't want the original to get changed.
        
        strcpy(convertedFileName, outputDirectory);
        if(convertedFileName[strlen(convertedFileName) - 1] != '/')
                strcat(convertedFileName, "/");
        strcat(convertedFileName, baseName);
        strcat(convertedFileName, ".fastk");
        convertedFile = fopen(convertedFileName, "w");
        
        // Safety check:
        if(convertedFile == 0)
        {
            printf("Could not create file: %s\n", convertedFileName);
            return 1;
        }
        
        outputReadsFASTK(convertedFile, reads[file], kmers, kmerSize);
    }
}

int processCorrection(int numInputFiles, char* inputFileNames, char* outputDirectory, 
        bool paired, bool substitutions, bool insertions, bool deletions, bool homopolymers,
        bool filtering, bool qualityUpdating) 
{
    printf("ERROR CORRECTION\n\n");
    
    // META OBJECT:
    Correction* correction = preprocessing(numInputFiles, inputFileNames, outputDirectory);
    
    // CORRECTION FUNCTION:    
    correction->correctionFunction = &correctRead;
    
    // ENABLED CORRECTIONS:
    correction->substitutions = substitutions;
    correction->insertions = insertions;
    correction->deletions = deletions;
    correction->homopolymers = homopolymers;
    
    correction->filtering = filtering;
    correction->qualityUpdating = qualityUpdating;
    
    // CORRECT READS:
    printf("Correcting reads...\n");  
    if(paired)
    {
        processPairedCorrection(correction);
    }    
    else
    {
        processSimpleCorrection(correction);
    }
    printf("Finished correcting reads!\n\n");
    
    return 0;
}
