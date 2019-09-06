// ;
// ;
//  esum.c
//  PStat3
//
//  Created by Bill Engels on 7/6/17.
//  Copyright © 2017 Bill Engels. All rights reserved.
//
/*
Copyright 2009-2014  William R. Engels

This file used in computing an exact statistical tests.
It is part of the 'PStat3' package for R
and is made available under the terms of the GNU General Public
License, version 2, or at your option, any later version,
incorporated herein by reference.

This program is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License
along with this program; if not, a copy is available at
http://www.r-project.org/Licenses/

*/


//
/************************************************************
 *                  For R Package `HWxtest`
 *
 * DESCRIPTION
 *
 *Function to perform an exact test for Hardy-
 *Weinberg proportions. The method uses recursion as described
 *by
 *
 *           Engels, 2009, Genetics 183, pp1431-1441
 *
 *************************************************************/

#define COUNTTYPE unsigned short

#ifdef DEBUG    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "esum.h"

    #define R_pow_di pow  // The following defines convert R calls to normal C
    #define Calloc(x,y) calloc((x),sizeof(y))
    #define Free free
    #define Rprintf printf
    #define lgammafn lgamma
#else
    #include <R.h>
    #include <Rmath.h>
#endif


#include <time.h>

// Globals

time_t start;
int timeLimit;
time_t start;
uint64_t nextCheck, checkGap;
int HN;
unsigned hdex;
double histoLeft, histoRight, statLeft, statRight, statSpan;
double * x; // The actual values -- n of them
uint64_t * hdata;
int64_t  ntrials = 0;
uint64_t nlower = 0, nupper= 0;


void subsetof(int k, int m, double stat)       // (unsigned) k from: (unsigned) m with: (FPTYPE) stat
{
    int bigusd, x1,x2, x3,x4, hdex;
    double stat1, stat2, stat3, stat4;
    if(ntrials > nextCheck) {
        nextCheck += checkGap;
        // See if time limit has been exceeded
//        Rprintf("\n%lld subsets so far", ntrials);
        if(time(NULL) - start > timeLimit) ntrials = -llabs(ntrials);
    }
    if (ntrials < 0) return;

    if (k<5) {
        switch (k) {
            case 4:
                for (x4=0; x4<m-2; x4++) {
                    stat4 = stat + x[x4];
                    for (x3=x4+1; x3<m-1; x3++) {
                        stat3 = stat4 + x[x3];
                        for (x2=x3+1; x2<m; x2++) {
                            stat2 = stat3 + x[x2];
                            for (x1=x2+1; x1<= m; x1++) {
                                stat1 = stat2 + x[x1];
                                ntrials ++;
//                                if (groupSwap) stat1 = swConst - stat1;
                                if (stat1 < statLeft) {
                                    nlower++;
                                } else if (stat1 > statRight) {
                                    nupper++;
                                }
                                    hdex = (unsigned) (statSpan * (stat1 - histoLeft));
                                    if (hdex >= 0 && hdex < HN) hdata[hdex]++;
                            }
                        }
                    }
                }
                break;
            case 3:
                for (x3=0; x3<m-1; x3++) {
                    stat3 = stat + x[x3];
                    for (x2=x3+1; x2<m; x2++) {
                        stat2 = stat3 + x[x2];
                        for (x1=x2+1; x1<= m; x1++) {
                            stat1 = stat2 + x[x1];
    //                        if (groupSwap) stat1 = swConst - stat1;
                            ntrials ++;
                            if (stat1 < statLeft) {
                               nlower++;
                            } else if (stat1 > statRight) {
                                nupper++;
                            }
                                hdex = (unsigned) (statSpan * (stat1 - histoLeft));
                                if (hdex >= 0 && hdex < HN) hdata[hdex]++;
                        }
                    }
                }
                break;
            case 2:
                for (x1=0; x1<m; x1++) {
                    stat2 = stat + x[x1];
                    for (x2=x1+1; x2<=m; x2++) {
                    	stat1 = 7;
                        stat1 = stat2 + x[x2];
//                        if (groupSwap) stat1 = swConst - stat1;
                        ntrials ++;
                        if (stat1 < statLeft) {
                            nlower++;
                        } else if (stat1 > statRight) {
                            nupper++;
                        }
                            hdex = (unsigned) (statSpan * (stat1 - histoLeft));
                            if (hdex >= 0 && hdex < HN) hdata[hdex]++;
                    }
                }
                break;
            case 1:
                for (x1=0; x1<= m; x1++) {
                    stat1 = stat + x[x1];
//                    if (groupSwap) stat1 = swConst - stat1;
                    ntrials ++;
                    if (stat1 < statLeft) {
                        nlower++;
                    } else if (stat1 > statRight) {
                        nupper++;
                    }
                        hdex = (unsigned) (statSpan * (stat1 - histoLeft));
                        if (hdex >= 0 && hdex < HN) hdata[hdex]++;
                }
            default:
                break;
        }
        if (0) {
            return;
        }
        
    } else {
        for (bigusd = m; bigusd >= k-1; bigusd--) {
            subsetof(k-1, bigusd-1, stat + x[bigusd]);
//            [self subsetof:k-1 from:bigusd-1 with:stat + x[bigusd]];
        }
    }
}


void esumTest(
              double * rx,  // n values to sum up k-subsets of
              int * rn, // The total number of values
              int * rk, // The number in each subset
              double * rnTodo, // number to examine - rsum only
              double * rstatLeft, // left boundary of sum
              double * rstatRight, // right boundary of sum for P value
              int * rhistobins, // number of bins for histogram from r
              double * rhistoLeft, // left end of range of histogram from r
              double * rhistoRight, // right end of range of histogram from r
              int * rsafeSecs, // Stop after this many seconds
              // Returned results
              double * rnlower, // Returned: number of subsets in left tail
              double * rnupper, // Returned: number of subsets in right tail
              double * rhistoData, // Returned: histogram data. length = histobins
              double * subsetsExamined,
              double * timespent // how much time it took
              )
{
    // Set up globals
    HN = *rhistobins;
    timeLimit = *rsafeSecs;
    time_t startclock = clock();
    statLeft = *rstatLeft;
    statRight = *rstatRight;
    statSpan = (HN-1)/(*rhistoRight - *rhistoLeft);
    histoLeft = *rhistoLeft;
    x = rx;
    ntrials = 0;
    nlower = 0;
    nupper = 0;
    checkGap = 1e6;
    nextCheck = checkGap;

    hdata = Calloc(HN + 1, uint64_t);
    for (int i = 0; i < HN; i++) hdata[i] = 0;
    start = time(NULL);
    
    subsetof(*rk, *rn - 1, 0);
    
    *rnlower = nlower;
    *rnupper = nupper;
    *subsetsExamined = ntrials;
   // Rprintf("\nSubsets examined = %1.0f\n", *subsetsExamined);
    for (int i = 0; i<HN; i++) {
        rhistoData[i] = hdata[i];
    }
    *timespent = ((double)clock() - startclock)/CLOCKS_PER_SEC;
   Free(hdata);
}
