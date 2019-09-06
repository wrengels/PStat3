//
//  PStatC.c
//  PStat3
//
//  Created by Bill Engels on 6/21/17.
//  Copyright © 2017 Bill Engels. All rights reserved.
//


#ifdef DEBUG    // This stuff can be removed in the R-only version
#include "PStatC.h"

void rmultinom(int n, double* prob, int K, int* rN) {
    int i, j;
    double x, cp;
    for (i = 0; i < K; i++) rN[i] = 0;
    for (j = 0; j < n ; j++) {
        x = (double)rand()/RAND_MAX;
        i = 0;
        cp=prob[0];
        while (x > cp) {
            cp += prob[++i];
        }
        (rN[i])++;
    }
}

#else
#include <R.h>
#include <Rmath.h>
#endif

#include <time.h>
#include <string.h>

//This call (MWCranx) generates a random unsigned long named ranx
#define MWCranx  t=a*Q[++ii]+c; c=(t>>32);  ranx=(unsigned)(t+c); if(ranx<c){ranx++;c++;} Q[ii]= ranx;
// Turn it into a random integer in the range 0 to k-1 with (unsigned)((unsigned short)ranx * (k)) >> 16;
#define MWCrank(k) (unsigned)(((unsigned long long)ranx * (unsigned long long)(k)) >> 32)

void sayit() {
    double x = 5;
    x = 3 * x;
    Rprintf("x is %4.7f",sqrt(x));
};
