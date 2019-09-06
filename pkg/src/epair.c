// ;
//  esum.c
//  PStat3
//
//  Created by Bill Engels on 7/6/17.
//  Copyright © 2017 Bill Engels. All rights reserved.
//
/*
Copyright 2009-2014  William R. Engels

This file used in computing an exact statistical test.
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
    #include "epair.h"

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

// Globals (We need global variables for use by recursive routines.)

extern time_t start;
extern int timeLimit;
extern time_t start;
extern uint64_t nextCheck, checkGap;
extern int HN;
extern unsigned hdex;
extern double histoLeft, histoRight, statLeft, statRight, statSpan;
extern double * x; // The actual values -- n of them
extern uint64_t * hdata;
extern int64_t  ntrials;
extern uint64_t nlower, nupper;
int npairs, ng8, counter;
double * x256; // Sums of 8-sets. total of ng8 * 256 of them.
double * x256last;
int isCancelled;
unsigned remaindr, remainderGap;

/*
 This recursive function computes the 256 possible sums of differences for the
 8 differences starting with the initial k
 */
//-(void) addDif: (int) k stat: (double) sum
void addDif(int k, double sum)
{
    double xk = x[k];
    if(k >= npairs) xk = 0; // same as padding x with zeroes.
    if((k+1)%8) {
        addDif(k+1, sum);
        addDif(k+1, sum + xk);
    } else {
        x256[counter ++] = sum * 2;
        x256[counter ++] = (sum + xk) * 2;
        
    }
}

void doStats(double * x256point, double sum) {
    double stat;
    int hdex;
    if(isCancelled) return;
    x256last = & (x256[256 * (ng8-1)]);
    if (x256point < x256last) {
        for (int i = 0; i < 256; i++) {
            doStats(x256point + 256, sum + x256point[i]);
        }
    } else {
        for (int i = 0; i < 256; i += remainderGap) {
            (ntrials)++;
            stat = sum + x256point[i];
            if (stat > statRight) nupper++;
            if (stat < statLeft) nlower++;
            hdex = (unsigned) (statSpan * (stat - histoLeft));
            if (hdex >= 0 && hdex < HN) hdata[hdex]++;
            if (ntrials > nextCheck) {
                nextCheck += checkGap;
                if (time(NULL) - start > timeLimit) {
                    isCancelled = 1;
                }
            }
        }
    }
    if (isCancelled) return;
}

void epairTest(
              double * rx,  // absolute values of pair differences
              int * rn, // The total number of pairs
//              int * rk, // The number in each subset
              double * rnTodo, // number to examine - rpair only
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
              double * timespent
              )
{
    // Set up globals
    HN = *rhistobins;
    timeLimit = *rsafeSecs;
    
    statLeft = *rstatLeft;
    statRight = *rstatRight;
    double maxsum = 0;
    for(int i=0;i<(*rn); i++) maxsum += rx[i];
    statSpan = (HN-1)/(*rhistoRight - *rhistoLeft);
    histoLeft = *rhistoLeft;
    x = rx;
    ntrials = 0;
    nupper = 0;
    nlower = 0;
    checkGap = 1e6;
    nextCheck = checkGap;
    npairs = *rn;
    ng8 = npairs/8;
    if(npairs % 8) ng8++;
    remaindr =  8- 8 * ng8 + *rn;
    remainderGap = R_pow_di(2, 8-remaindr);
    counter = 0;

    hdata = Calloc(HN + 1, uint64_t);
    for (int i = 0; i < HN; i++) hdata[i] = 0;
    x256 = Calloc(ng8 * 256, double);
    for(int i=0; i< ng8*256; i++) x256[i] = 0;
    start = time(NULL);
    time_t startclock = clock();
    // Set up x256 array using a recursive routine.
    //This array pre-calculates all 256 possible sums for each group of 8 x's
    
    for(int i=0; i<ng8; i++)  {
        addDif(8 * i, 0);
    }
    
    //Rprintf("\nHere are the difs, recursively computed:\n");
   // for(int i=0; i< ng8*256; i++) Rprintf("%f, ", x256[i]);
        
    isCancelled = 0;
    x256last = & (x256[256 * (ng8-1)]);
    
    // Start Main Loop
    
    doStats(x256, -maxsum);  // This is where most of the work happens
        
    *rnlower = nlower;
    *rnupper = nupper;
    *subsetsExamined = ntrials;
    for (int i = 0; i<HN; i++) {
        rhistoData[i] = hdata[i];
    }
    double timespent_double = ((double)clock() - startclock)/CLOCKS_PER_SEC;
    *timespent = timespent_double;
    
    Free(hdata);
    Free(x256);
}
