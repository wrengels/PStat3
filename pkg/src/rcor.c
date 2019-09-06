// ;
//  rsum.c
//  PStat2.5
//
//  Created by Bill Engels on 7/16/17.
//  Copyright © 2017 Bill Engels. All rights reserved.
//


#define COUNTTYPE unsigned short

#ifdef DEBUG    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <string.h>  // needed for memset
    #include "rpair.h"

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


void rcorTest(
              double * rx,  // x values
              double * ry,  // y values
              int * rn,     // The total number of pairs
//              int * rk,     // The number of values in x group
              double * rnTodo, // requested number of subsets to examine - for rpair only
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
              double * subsetsExamined,  // the number actually generated
              double * timespent // how much time it took
              )
{
               int i,j,r;
               time_t start, startclock;
               double *x, *y, *xy;
               int counter;
               int npairs;
               int timeLimit = *rsafeSecs;
               uint64_t checkGap = 1e4;
               uint64_t nextCheck = checkGap;
               int HN;
               double histoLeft, histoRight, statLeft, statRight, statSpan;
               int64_t ntrials;
               uint64_t * hdata;
               uint64_t nlower, nupper;
               unsigned hdex;
               double stat;
               int isCancelled;
               double maxsum = 0;
//               double maxCorrelation = rhistoData[0];
 //              double minCorrelation = rhistoData[1];
               double aconst = rhistoData[2];
               double cconst = rhistoData[3];
               unsigned rseed = (unsigned)(rhistoData[4] + .5);
                  
                  
    // for MWC256
    uint32_t Q[256];
    uint64_t t=0, c = 362436;
    uint8_t ii = 255;
    Q[0] = 0;
    uint64_t a=1540315826LL;
    uint32_t ranx;
    srandom((unsigned)time(NULL)/2 + rseed/2);
    for (i = 0; i < 256; i++) Q[i] = (uint32_t)random();
    
    //This call (MWCranx) generates a random unsigned long named ranx
#define MWCranx  t=a*Q[++ii]+c; c=(t>>32);  ranx=(uint32_t)(t+c); if(ranx<c){ranx++;c++;} Q[ii]= ranx;
    // Turn it into a random integer in the range 0 to k-1 with (unsigned)((unsigned short)ranx * (k)) >> 16;
#define MWCrank(k) (unsigned)(((uint64_t)ranx * (uint64_t)(k)) >> 32)
    
    //Run off some random MWCs just to mix things up
    for (i=0; i<1000; i++) MWCranx;

 
    
                  // Set up globals
 
                    int nreps = 10;
                    if(*rnupper > .1) nreps = (int)(fabs(*rnupper) + .1); // to set nreps programmatically
                  HN = *rhistobins;
                  timeLimit = *rsafeSecs;
                  
                  statLeft = *rstatLeft;
                  statRight = *rstatRight;
                  for(int i=0;i<(*rn); i++) maxsum += rx[i];
                  histoLeft = *rhistoLeft;
                  histoRight = *rhistoRight;
                  statSpan = (HN-1)/(histoRight - histoLeft);
                  x = rx;
                  y = ry;
                  ntrials = 0;
                  uint64_t nTodo = *rnTodo;
                  nupper = 0;
                  nlower = 0;
                  checkGap = 1e6;
                  nextCheck = checkGap;
                  npairs = *rn;
//                  ng8 = (npairs + 7)/8;
//                  ng32 = (ng8 + 3)/4;
//                  remaindr =  8- 8 * ng8 + *rn;
//                  remainderGap = R_pow_di(2, 8-remaindr);
                  counter = 0;
                  
                  hdata = Calloc(HN + 1, uint64_t);
                  xy = Calloc((npairs) * (npairs), double);
    
                  for(int i = 0; i < (npairs); i++) {
                      for(int j = 0; j < (npairs); j++) {
                          xy[i * npairs + j] = x[i] * y[j] * aconst + cconst;
                      }
                  }
                  
                  for (int i = 0; i < HN; i++) hdata[i] = 0;
                  start = time(NULL);
                  startclock = clock();

    // Make and scramble a collection of nreps sets of integers
    unsigned *dc, *di, * d = Calloc(nreps * npairs, unsigned);
    dc = d;
    // fill with the integers
    for (i=0; i<nreps; i++) {
        for (j=0; j<npairs ; j++) {
            dc[j] = j;
        }
        dc+= npairs;
    }
    unsigned temp, *wi,  *wlast, randomRange;
    // scramble each of the NREPS sets
    dc = d;
    for (i=0; i<nreps; i++) {
        for (j=0; j<npairs; j++) {
            MWCranx r = MWCrank(npairs - j);
            temp = dc[r+j];
            dc[r+j] = dc[j];
            dc[j] = temp;
        }
        dc += npairs;
    }
    
    // make another set of integers to be scrambled for EACH iteration
    unsigned *wc, * w = Calloc(npairs, unsigned);
    for(i=0; i<npairs; i++) w[i] = i;
    isCancelled = 0;
    
    // Now, here comes the main loop;!
    while(ntrials < nTodo && !isCancelled) {
        randomRange = npairs;
        wlast = &w[npairs-1];
        // randomly permute the w integers
        for(wi=w; wi <= wlast; wi++) {
            MWCranx r = MWCrank(randomRange--);
            wc = wi + r;
            temp = *wc;
            *wc = *wi;
            *wi = temp;
        } // for wi
        
        //Find correlation for each of the nreps sets
        dc = d;
        for (i = 0; i < nreps; i++) {
            stat = 0; di = dc;
            for (wi = w ; wi <= wlast; wi++) {
                stat += xy[npairs * (*wi) + *di];
                di++;
            }
            //         NSLog(@"\nstat = %12.9f\n\n**************************\n", stat);
            dc += npairs;
            if (stat > statRight) {
                nupper++;
            } else if (stat < statLeft) {
                nlower++;
            }
                hdex = (unsigned) (statSpan * (stat - histoLeft));
                if (hdex < HN) hdata[hdex]++;
        }
        ntrials += nreps;
        if (ntrials > nextCheck){// It's time to check
            nextCheck += checkGap;
            if(time(NULL) - start >= *rsafeSecs) {//Time to stop
                isCancelled = 1;
                Rprintf("\nProcess timed out\n");
            }
        }
    } // while ntrials
    
                  
                  *rnlower = nlower;
                  *rnupper = nupper;
                  *subsetsExamined = ntrials;
 //                 Rprintf("\nSubsets examined = %1.0f\n", *subsetsExamined);
                  for (int i = 0; i<HN; i++) {
                      rhistoData[i] = hdata[i];
                  }
                  *timespent = ((double)clock() - startclock)/CLOCKS_PER_SEC;
                  
 //                 Rprintf ("\nB is set to %f \nAnd there were %lld trials taking %f seconds\n", *rnTodo, ntrials, *timespent);
                  
                    Free(hdata);
                    Free(xy);
                    Free(d);
                    Free(w);
                  
 //                 Rprintf("\nJust finished running rcor\n");

}
