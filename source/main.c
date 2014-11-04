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

#include "Utility.h"
#include "ErrorProcessing.h"
#include <string.h>
#include <stdlib.h>
#include "Globals.h"
#include "Reads.h"
#include <unistd.h>

bool checkInput(int numInputFiles, char* inputFileNames, char* outputFileName, 
        bool paired, enum SEQUENCING_TECHNOLOGY type, bool fastk, unsigned int KMER_SIZE)
{
    if (numInputFiles < 1)
    {
        printf("ERROR: Need at least one input file.\n");
        return false;
    }
            
    if (paired && numInputFiles != 2)
    {
        printf("ERROR: Provide only two input files when specifying paired input.\n");
        return false;
    }
    
    if (KMER_SIZE < 4)
    {
        printf("ERROR: k-mer size is too small.\n");
        return false;
    }
    
    if (KMER_SIZE > 31)
    {
        printf("ERROR: k-mer size is too large.\n");
        return false;
    }
    
    return true;
}

void help()
{
    printf("\n");
    printf("USAGE: \n");
    printf("\n");
    
    printf("ERROR CORRECTION\n");
    printf("Required: \n");
    printf("\n");
    printf("\t-i \t[file] \tSpecify one or many FASTQ input files.\n");
    printf("\n");
    
    printf("Optional: \n");
    printf("\n");
    printf("\t-o \t \tOutput directory.\n");
    printf("\t-p \t \tSpecify input should be treated as paired.\n");
    printf("\n");
    printf("\t-s \t[bool] \tSubstitution corrections. \"true\" or \"false\".\n");
    printf("\t-n \t[bool] \tInsertion corrections. \"true\" or \"false\".\n");
    printf("\t-d \t[bool] \tDeletion corrections. \"true\" or \"false\".\n");
    printf("\t-h \t[bool] \tHomopolymer corrections. \"true\" or \"false\".\n");
    printf("\t-f \t[bool] \tLow k-mer read filtering. \"true\" or \"false\".\n");
    //printf("\t-q \t[bool] \tQuality score updating. \"true\" or \"false\".\n");
    printf("\n");
    printf("\t-k \t[int] \tSpecify the k-mer size.\n");
    printf("\t-b \t[int] \tSpecify the input batch size.\n");
    printf("\n");
    
    printf("FASTK CONVERSION\n");
    printf("Required: \n");
    printf("\t-fastk \t\tFASTQ to FASTK file conversion.\n");
    printf("\t-i \t[file] \tSpecify one or many input files.\n");
    printf("\n");
    
    printf("EXAMPLES: \n");
    printf("\n");
    printf("./error -i file1.fastq\n");
    printf("\tCorrect a single file.\n");
    printf("\n");
    printf("./error -i file1.fastq file2.fastq file3.fastq\n");
    printf("\tCorrect multiple file together.\n");
    printf("\n");
    printf("./error -i frag_1.fastq frag_2.fastq -p\n");
    printf("\tCorrect two paired files.\n");
    printf("\n");
    
    printf("./error -fastk -i file1.fastq\n");
    printf("\tConvert file1.fastq to FASTK format.\n");
    
    printf("\n");    
}

int main(int argc, char**argv) 
{
    char outputDirectory[1024];
    enum SEQUENCING_TECHNOLOGY type = UNKNOWN;
    
    int numInputFiles = 0;
    char* inputFileNames;
    
    bool paired = false;
    bool fastk = false;
    
    // Corrections:
    bool substitutions = true;
    bool insertions = true;
    bool deletions = true;
    bool homopolymers = true;
    bool qualityUpdating = true;
    bool filtering = true;
    
    printf("\n");
    printf("Pollux 1.00\n");
    printf("Source compiled on %s at %s.\n", __DATE__, __TIME__);
    printf("\n");
    
    if(argc == 1)
    {
        help();
        return 1;
    }
    
    // Set default output directory:
    getcwd(outputDirectory, sizeof(outputDirectory));
    
    for(int i = 1; i < argc; i++)
    {
        printf("Considering argument: %s ", argv[i]);
        
        // HELP
        if(strcmp("-help", argv[i]) == 0 || strcmp("--help", argv[i]) == 0)
        {
            help();
            return 1;
        }        
        // INPUT FILES
        else if(strcmp("-i", argv[i]) == 0 && i < (argc - 1))
        {
            // Determine the number of input files.
            for(int arg = i + 1; arg < argc && argv[arg][0] != '-'; arg++)
            {
                numInputFiles++;
            }
            
            // Allocate space for names.
            inputFileNames = (char*) malloc (sizeof(char*) * numInputFiles * 200);
            
            // Grab file names.
            for(int file = 0; file < numInputFiles; file++)
            {
                strcpy(&(inputFileNames[file * 200]), argv[i + 1 + file]);
            }
            
            i += numInputFiles;
            
            printf(": number of input files is %d\n", numInputFiles);
        }
        // OUTPUT FILE(S)
        else if(strcmp("-o", argv[i]) == 0 && i < (argc - 1))
        {
            strcpy(outputDirectory, argv[i + 1]);
            i++;
            
            printf(": output directory is %s\n", outputDirectory);
        }
        // TYPE OF DATA
        else if(strcmp("-t", argv[i]) == 0)
        {
            // Illumina
            if(strcmp("illumina", argv[i + 1]) == 0)
            {
                type = ILLUMINA;
                printf(": input type is Illumina\n");
            }
            else if(strcmp("ion", argv[i + 1]) == 0)
            {
                type = ION;
                printf(": input type is Ion Torrent\n");
            }
            else if(strcmp("454", argv[i + 1]) == 0)
            {
                type = ROCHE454;
                printf(": input type is 454\n");
            }
            else
            {
                type = UNKNOWN;
                printf(": input type is unknown!\n");
            }
            
            i++;            
        }
        // PAIRED INPUT
        else if(strcmp("-p", argv[i]) == 0)
        {
            paired = true;
            
            printf(": paired input\n");
        }
        // BATCH SIZE
        else if(strcmp("-b", argv[i]) == 0)
        {
            BATCH_SIZE = atoi(argv[i + 1]);
            
            printf(": batch size is %d\n", BATCH_SIZE);
            
            i++; 
        }
        // FASTK
        else if(strcmp("-fastk", argv[i]) == 0)
        {
            fastk = true;
            
            printf(": FASTK conversion\n");
        }
        // SUBSTITUTIONS
        else if(strcmp("-s", argv[i]) == 0)
        {
            // Enabled
            if(strcmp("true", argv[i + 1]) == 0)
            {
                substitutions = true;
                printf(": substitutions enabled\n");
            }
            else if(strcmp("false", argv[i + 1]) == 0)
            {
                substitutions = false;
                printf(": substitutions disabled\n");
            }
            else
            {
                printf(": input not understood!\n");
            }
            
            i++; 
        }
        // INSERTIONS
        else if(strcmp("-n", argv[i]) == 0)
        {
            // Enabled
            if(strcmp("true", argv[i + 1]) == 0)
            {
                insertions = true;
                printf(": insertions enabled\n");
            }
            else if(strcmp("false", argv[i + 1]) == 0)
            {
                insertions = false;
                printf(": insertions disabled\n");
            }
            else
            {
                printf(": input not understood!\n");
            }
            
            i++; 
        }
        // DELETIONS
        else if(strcmp("-d", argv[i]) == 0)
        {
            // Enabled
            if(strcmp("true", argv[i + 1]) == 0)
            {
                deletions = true;
                printf(": deletions enabled\n");
            }
            else if(strcmp("false", argv[i + 1]) == 0)
            {
                deletions = false;
                printf(": deletions disabled\n");
            }
            else
            {
                printf(": input not understood!\n");
            }
            
            i++; 
        }
        // HOMOPOLYMERS
        else if(strcmp("-h", argv[i]) == 0)
        {
            // Enabled
            if(strcmp("true", argv[i + 1]) == 0)
            {
                homopolymers = true;
                printf(": homopolymers enabled\n");
            }
            else if(strcmp("false", argv[i + 1]) == 0)
            {
                homopolymers = false;
                printf(": homopolymers disabled\n");
            }
            else
            {
                printf(": input not understood!\n");
            }
            
            i++; 
        }
        // FILTERING
        else if(strcmp("-f", argv[i]) == 0)
        {
            // Enabled
            if(strcmp("true", argv[i + 1]) == 0)
            {
                filtering = true;
                printf(": read filtering enabled\n");
            }
            else if(strcmp("false", argv[i + 1]) == 0)
            {
                filtering = false;
                printf(": read filtering disabled\n");
            }
            else
            {
                printf(": input not understood!\n");
            }
            
            i++; 
        }
        // QUALITY UPDATING
        /*
        else if(strcmp("-q", argv[i]) == 0)
        {
            // Enabled
            if(strcmp("true", argv[i + 1]) == 0)
            {
                qualityUpdating = true;
                printf(": quality updating enabled\n");
            }
            else if(strcmp("false", argv[i + 1]) == 0)
            {
                qualityUpdating = false;
                printf(": quality updating disabled\n");
            }
            else
            {
                printf(": input not understood!\n");
            }
            
            i++; 
        }
        */
        // KMER SIZE
        else if(strcmp("-k", argv[i]) == 0)
        {
            KMER_SIZE = atoi(argv[i + 1]);
            
            printf(": k-mer size is %d\n", KMER_SIZE);
            
            i++; 
        }
        
        // PROBLEM
        else
        {
            printf("\nProblem with argument: %s\n", argv[i]);
            return 1;
        }
    }
    
    printf("\n");
    
    for(int file = 0; file < numInputFiles; file++)
    {
        printf("File name: %s\n", &(inputFileNames[file * 200]));
    }
    
    printf("\n");
    
    if(checkInput(numInputFiles, inputFileNames, outputDirectory, paired, type, fastk, KMER_SIZE))
    {
        // FASTK CONVERSION
        if(fastk)
        {
            convertFASTQToFASTK(numInputFiles, inputFileNames, outputDirectory);
        }
        // ERROR CORRECTION
        else
        {
            processCorrection(numInputFiles, inputFileNames, outputDirectory, paired,
                    substitutions, insertions, deletions, homopolymers, filtering, qualityUpdating);
        }        
    }
    
    return 0;    
}
