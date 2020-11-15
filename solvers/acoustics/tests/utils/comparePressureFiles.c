#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * Function to compare pressures in two files outputted from libParanumal.
 * Returns 0 if both files are equivalent, otherwise returns
 * -1 and sets line where both file differ.
 */
int comparePressureFiles(FILE * fPtr1, FILE * fPtr2, int * line, double eps)
{
    *line  = 0;

    float ts1;
    float ts2;
    float pres1;
    float pres2;
    
    int res1 = 0;
    int res2 = 0;
    
    res1 = fscanf(fPtr1, "%f %f", &ts1, &pres1);
    res2 = fscanf(fPtr2, "%f %f", &ts2, &pres2);

    if (res1 == EOF && res2 == EOF) {
        // file is empty
        return 0;
    }        
    else if (res1 != res2) {
        // one of the files is empty
        return -1;
    }

    do
    {
        // if time vector or pressure vector is not same then return -1

        if (fabs(ts1 - ts2) > eps)
            return -1;

        if (fabs(pres1 - pres2) > eps)
            return -1;

        res1 = fscanf(fPtr1, "%f %f", &ts1, &pres1);
        res2 = fscanf(fPtr2, "%f %f", &ts2, &pres2);

        *line  += 1;

    } while (res1 != EOF && res2 != EOF);


    /* If both files have reached end */
    if (res1 == EOF && res2 == EOF)
        return 0;
    else
        return -1;
}