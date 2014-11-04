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
#include "Utility.h"

#ifndef READS_H
#define	READS_H

#ifdef	__cplusplus
extern "C" {
#endif
    
extern int BATCH_SIZE;  // Batch size in reads.

// Single read:
struct read {
    char* seqName1;
    unsigned long long int* sequence;
    char* seqName2;
    char* quality;   
    
    struct contig* basecontig;
    unsigned long int startpos;
    int length;
    int number;
    short forward;
    short correct_pos;
    char type;
};

// Collection of reads:
typedef struct
{
    char* fileName;
    FILE* file;
    int total;
    
    int current;
    int ID;
    
    struct read* readData;
    
} Reads;

Reads* createReads(char* fileName);
struct read* readsGetNext(Reads* reads);
bool readsHasNext(Reads* reads);
int readsReset(Reads* reads);
void readsDestroy(Reads* reads);
int readsGetCount(Reads* reads);
char* readsGetFileName(Reads* reads);

#ifdef	__cplusplus
}
#endif

#endif	/* KMERHASHTABLE_H */

