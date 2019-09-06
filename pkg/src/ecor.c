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

time_t startt;
extern int timeLimit;
extern uint64_t nextCheck, checkGap;
extern int HN;
extern unsigned hdex;
extern double histoLeft, histoRight, statLeft, statRight, statSpan;
extern double * x; // The actual values -- n of them
double * y, * xy;
extern uint64_t * hdata;
extern int64_t  ntrials;
extern uint64_t nlower, nupper;
extern int npairs, counter;
extern int isCancelled;
int * d;
int safeSecs;
double *xn1, *xn2, *xn3;

void permute(int start, double stat) {
    
    /*
     finds all permutations of an array of integers d[0] to d[n-1] starting with 'start'
     
     After each permutation is completed, the floating point value 'stat' contains the
     corresponding correlation coefficient associated with that permutation. Computing
     this value requires that an array 'xy' be set up in advance. The used of the passed
     parameter 'stat'  to distribute the computation efficiently is original as far as I know.
     
     Here is a web site discussing various permutation algorithms. They refer to the
     one used here as "exeter" but it is a common method.
     http://www.bearcave.com/random_hacks/permute.html
     */
    
    // done with a permutation -- but for the final six.
    //Doing the final six together increases the overall speed.
    //xn1, xn2 and xn3 point to the last three rows of xy
    int tmp;
    double r, r0;
    
    if(isCancelled) return;
    if (ntrials > nextCheck){// It's time to check
        nextCheck += checkGap;
        if(time(NULL) - startt >= safeSecs) {//Time to stop
            isCancelled = 1;
            Rprintf("\nProcess timed out\ntime(NULL)= %ld, start= %d", time(NULL), start);
        }
    }
    int n = npairs;
    if(start == npairs - 3) {
        ntrials += 6;
        
        // do ....... 3-2-1
        r0 = stat + xn3[d[n-3]];
        r  = r0 + xn2[d[n-2]] + xn1[d[n-1]];
        if (r < statLeft) nlower++;
        else if (r > statRight) nupper++;
        hdex = (unsigned) (statSpan * (r - histoLeft));
        if (hdex < HN) hdata[hdex]++;
        
        // do ....... 3-1-2
        r  = r0 + xn2[d[n-1]] + xn1[d[n-2]];
        if (r < statLeft) nlower++;
        else if (r > statRight) nupper++;
        hdex = (unsigned) (statSpan * (r - histoLeft));
        if (hdex < HN) hdata[hdex]++;

        // do ....... 1-2-3
        r0 = stat + xn3[d[n-1]];
        r  = r0 + xn2[d[n-2]] + xn1[d[n-3]];
        if (r < statLeft) nlower++;
        else if (r > statRight) nupper++;
        hdex = (unsigned) (statSpan * (r - histoLeft));
        if (hdex < HN) hdata[hdex]++;

        // do ....... 1-3-2
        r  = r0 + xn2[d[n-3]] + xn1[d[n-2]];
        if (r < statLeft) nlower++;
        else if (r > statRight) nupper++;
        hdex = (unsigned) (statSpan * (r - histoLeft));
        if (hdex < HN) hdata[hdex]++;

        // do ....... 2-1-3
        r0 = stat + xn3[d[n-2]];
        r  = r0 + xn2[d[n-1]] + xn1[d[n-3]];
        if (r < statLeft) nlower++;
        else if (r > statRight) nupper++;
        hdex = (unsigned) (statSpan * (r - histoLeft));
        if (hdex < HN) hdata[hdex]++;
        
        // do ....... 2-3-1
        r  = r0 + xn2[d[n-3]] + xn1[d[n-1]];
        if (r < statLeft) nlower++;
        else if (r > statRight) nupper++;
        hdex = (unsigned) (statSpan * (r - histoLeft));
        if (hdex < HN) hdata[hdex]++;

        // just completed a group of six trials
        if(isCancelled) return;
        
    } else {
        if(isCancelled) return;
        for(int i=start; i<npairs; i++) {
            tmp = d[i];
            d[i] = d[start];
            d[start] = tmp;
            permute(start + 1, stat + xy[tmp + npairs * start]);
            d[start] = d[i];
            d[i] = tmp;
        }
    }
}

void ecorTest(
              double * rx,  // x values
              double * ry,  // y values
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
    safeSecs = * rsafeSecs;
    statLeft = *rstatLeft;
    statRight = *rstatRight;
    double maxsum = 0;
    for(int i=0;i<(*rn); i++) maxsum += rx[i];
    statSpan = (HN-1)/(*rhistoRight - *rhistoLeft);
    histoLeft = *rhistoLeft;
    x = rx; y = ry;
    ntrials = 0;
    nupper = 0;
    nlower = 0;
    checkGap = 1e6;
    nextCheck = checkGap;
    npairs = *rn;
    double aconst = rhistoData[2];
    double cconst = rhistoData[3];
//    double maxCorrelation = rhistoData[0];
//    double minCorrelation = rhistoData[1];
    
    hdata = Calloc(HN + 1, uint64_t);
    xy = Calloc((*rn) * (*rn), double);
    d = Calloc(npairs, int);
    
    for(int i = 0; i < (*rn); i++) {
        for(int j = 0; j < (*rn); j++) {
            xy[i * npairs + j] = x[i] * y[j] * aconst + cconst;
        }
    }
    double obcor = 0;
    for(int i = 0; i < npairs; i++) obcor += xy[i*npairs + i];
    for (int i = 0; i < HN; i++) hdata[i] = 0;
    startt = time(NULL);
    time_t startclock = clock();
 
    for(int i = 0; i< npairs; i++) d[i] = i;
    int nn1, nn2, nn3;
    // just so nn1, nn2 and nn3 don't have to be recalculated a zillion times in permute
    nn1 = npairs * npairs - npairs;
    nn2 = npairs * (npairs-2);
    nn3 = npairs * (npairs-3);
    xn1 = &(xy[nn1]);
    xn2 = &(xy[nn2]);
    xn3 = &(xy[nn3]);
    
    isCancelled = 0;
    
    // Start Main Loop
    permute(0, 0);
        
    *rnlower = nlower;
    *rnupper = nupper;
    *subsetsExamined = ntrials;
    for (int i = 0; i<HN; i++) {
        rhistoData[i] = hdata[i];
    }
    double timespent_double = ((double)clock() - startclock)/CLOCKS_PER_SEC;
    *timespent = timespent_double;
    Free(hdata);
    Free(xy);
    Free(d);
}
