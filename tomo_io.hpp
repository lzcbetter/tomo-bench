#ifndef _TOMO_IO_
#define _TOMO_IO_

#include <stdio.h>


/*
***************************************************************************************************
                    Global Variables
                    all capital variables denotes parameter
***************************************************************************************************
*/
unsigned int P;      // number of projections
unsigned int S;      // number of sinograms
unsigned int C;      // number of columns
float *pin_buf;      // pointer to input buffer
float *pout_buf;     // pointer to output buffer

#endif