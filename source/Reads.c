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
#include <string.h>

#include <ctype.h>
#include "Reads.h"
#include "Encoding.h" 
#include "Utility.h"

int BATCH_SIZE = 200000;        // Number of reads loaded in memory.
int NUCLEOTIDE = 0;             // [0, 1, 2, 4] : replaces N's deterministically 

Reads* createReads(char* fileName)
{
    Reads* reads = (Reads*)malloc(sizeof(Reads));
    
    FILE* file = fopen(fileName, "r");
    
    if (file == 0) 
    {
        printf("Could not open file location: %s for reading.\n", fileName);
        return 0;
    }
    
    char ch;
    int lines = 0;
    
    // Count the number of lines in the file:
    while (EOF != (ch = getc(file))) 
    {
        if ('\n' == ch) {
            ++lines;
        }
    }
    
    fclose(file);
    file = fopen(fileName, "r");
    
    if (file == 0) 
    {
        printf("Could not open file location: %s for reading.\n", fileName);
        return 0;
    }
    
    reads->fileName = fileName;
    reads->file = file;
    reads->total = lines / 4;
    
    reads->current = 0;
    reads->readData = 0;
    reads->ID = 0;    
    
    return reads;
}

void freeReads(Reads* reads)
{
    // Do we already have reads open?
    if(reads->readData != 0)
    {
        int limit;
        
        // We don't want to free memory we haven't allocated.
        // How many reads do we need to free in this block? The last block 
        //      might not need all freed.
        
        // The entire block needs to be freed.
        // This is the case when we're not in the last, or the last block was
        //      entirely full.
        if(reads->current < reads->total || reads->total % BATCH_SIZE == 0)
        {
            limit = BATCH_SIZE;
        }
        // Only part of the block needs to be freed.
        else
        {
            limit = reads->total % BATCH_SIZE;
        }
        
        // Free all the open reads.
        for(int i = 0; i < limit; i++)
        {
            struct read* current = &(reads->readData[i]);
            
            free(current->seqName1);
            free(current->sequence);
            free(current->seqName2);
            free(current->quality);
            
            free(current->basecontig);
        }
        
        free(reads->readData);
    }
}

char getNextReplacementNucleotide()
{
    char result = 'A';
    
    switch(NUCLEOTIDE)
    {
        case 0:
            result = 'A';
            break;
        case 1:
            result = 'C';
            break;
        case 2:
            result = 'G';
            break;
        case 3:
            result = 'T';
            break;
    }
    
    NUCLEOTIDE = (NUCLEOTIDE + 1) % 4;
    
    return result;
}

// Replaces internal N's with other nucleotides.
void replaceN(char* string)
{
    for(int i = 0; i < strlen(string); i++)
    {
        if(string[i] == 'N' || string[i] == 'n')
        {
            string[i] = getNextReplacementNucleotide();
        }
    }
}

void deleteCharacter(char* string, int pos)
{
    int length = strlen(string);
    
    if (pos < 0 || pos >= length)
        return;
    
    for(int i = pos; i < length; i++)
    {
        string[i] = string[i + 1];
    }
    
    string[length] = '\0';
}

void trimNs(char* sequence, char* quality)
{
    // Delete leading N's.
    while(strlen(sequence) >= 1 && 
            (sequence[0] == 'N' || sequence[0] == 'n'))
    {
        deleteCharacter(sequence, 0);
        deleteCharacter(quality, 0);
    }
    
    // Delete trailing N's.
    while(strlen(sequence) >= 1 && 
            (sequence[strlen(sequence) - 1] == 'N' || sequence[strlen(sequence) - 1] == 'n'))
    {
        deleteCharacter(sequence, strlen(sequence) - 1);
        deleteCharacter(quality, strlen(quality) - 1);
    }    
}

void trimSpaces(char* string)
{
    // Delete leading spaces.
    while(strlen(string) >= 1 && isspace(string[0]))
    {
        deleteCharacter(string, 0);
    }
    
    // Delete trailing N's.
    while(strlen(string) >= 1 && isspace(string[strlen(string) - 1]))
    {
        deleteCharacter(string, strlen(string) - 1);
    }
}

void loadReads(Reads* reads)
{
    freeReads(reads);   // FREE FIRST!    
    reads->readData = (struct read*)malloc(BATCH_SIZE * sizeof(struct read));      // LOAD SECOND!
      
    for(int i = 0; i < BATCH_SIZE && (reads->current + i < reads->total); i++)
    {
        struct read* current = &(reads->readData[i]);
        
        reads->ID = reads->ID + 1;

        // TODO: THESE LENGTHS SHOULD BE VARIABLE!!!
        char seqName1[2048];
        char sequence[2048];
        char seqName2[2048];
        char quality[2048];
        
        // Get the 4 lines of a FASTQ file:
        if(     fgets(seqName1, 2048, reads->file) == NULL ||
                fgets(sequence, 2048, reads->file) == NULL ||
                fgets(seqName2, 2048, reads->file) == NULL ||
                fgets(quality, 2048, reads->file) == NULL )
        {
            printf("CRITICAL: FAILED TO READ INPUT!\n");
            exit(1);
        }
       
        trimSpaces(sequence);
        trimSpaces(quality);
    
        trimNs(sequence, quality);
        replaceN(sequence);
        
        // Sequence name 1:
        current->seqName1 = malloc(strlen(seqName1) + 1);
        strcpy(current->seqName1, seqName1);
        
        // Encode sequences:
        encode_sequence(current, sequence);
        
         // Sequence name 2:
        current->seqName2 = malloc(strlen(seqName2) + 1);
        strcpy(current->seqName2, seqName2);       
        
        //Quality:
        current->quality = malloc(strlen(quality) + 1);
        strcpy(current->quality, quality);
        
        current->basecontig = NULL;
        current->correct_pos = 0;
        current->number = reads->ID;        
    }
}

struct read* readsGetNext(Reads* reads)
{
    struct read* result;
    
    // Do we need to load more reads?
    if(reads->current % BATCH_SIZE == 0)
    {
        loadReads(reads);
    }
    
    result = &(reads->readData[reads->current % BATCH_SIZE]);
    reads->current = reads->current + 1;    
    
    return result;
}

bool readsHasNext(Reads* reads)
{
    return (reads->current < reads->total);
}

int readsReset(Reads* reads)
{   
    freeReads(reads);
    fclose(reads->file);    
    
    reads->file = fopen(reads->fileName, "r");
    
    if (reads->file == 0) 
    {
        printf("Could not open file location: %s for reading.\n", reads->fileName);
        return 1;
    }
    
    reads->current = 0;
    reads->readData = 0;
    reads->ID = 0;
    
    return 0;
}

void readsDestroy(Reads* reads)
{
    freeReads(reads);
    fclose(reads->file);
    
    free(reads);
}

int readsGetCount(Reads* reads)
{
    return reads->total;
}

char* readsGetFileName(Reads* reads)
{
    return reads->fileName;
}
