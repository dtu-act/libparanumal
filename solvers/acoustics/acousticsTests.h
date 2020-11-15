#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string>
#include <unistd.h>


/**
 * Function to compare two files.
 * Returns 0 if both files are equivalent, otherwise returns
 * -1 and sets line and col where both file differ.
 */
int compareFiles(FILE * fPtr1, FILE * fPtr2, int * line, int * col);

/**
 * Function to compare pressures in two files outputted from libParanumal.
 * Returns 0 if both files are equivalent, otherwise returns
 * -1 and sets line where both file differ.
 */
int comparePressureFiles(FILE * fPtr1, FILE * fPtr2, int * line, double eps);

/**
 * Function to open file given dir path and filename. Sets filehandle iFP and resulting filepath
 * Returns 0 is success, -1 otherwise
 */
int openFile(std::string filepath, FILE **iFP);