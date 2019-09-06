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
    #include "rsum.h"

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

extern time_t start;
//int timeLimit;
//time_t start;
//uint64_t nextCheck, checkGap;
//int HN;
//unsigned hdex;
//double histoLeft, histoRight, statLeft, statRight, statSpan;
//double * x; // The actual values -- n of them
//uint64_t * hdata;
//int64_t  ntrials = 0;
//uint64_t nlower = 0, nupper= 0;


void rsumTest(
              double * rx,  // n values to sum up k-subsets of
              int * rn,     // The total number of values, including both x and y
              int * rk,     // The number of values in x group
              double * rnTodo, // requested number of subsets to examine - for rsum only
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
              ) {
    int i,j,r;
                  unsigned nreps = 10;
                  if(*rnupper > .1) nreps = (int)(fabs(*rnupper) + .1); // to set nreps programmatically
                  double * x;
                  x = rx;
                 start = time(NULL);
                  time_t startclock = clock();
                  uint64_t checkGap = 1e4;
                  uint64_t nextCheck = checkGap;
    
                  
                  
    // for MWC256
      uint32_t Q[256];
      uint64_t t=0, c = 362436;
      uint8_t ii = 255;
      Q[0] = 0;
      uint64_t a=1540315826LL;
      uint32_t ranx;
      unsigned rseed = (unsigned)(rhistoData[4] + .4);
      srandom((unsigned)time(NULL)/2 + rseed/2);
    for (i = 0; i < 256; i++) Q[i] = (uint32_t)random();
    
    //This call (MWCranx) generates a random unsigned long named ranx
#define MWCranx  t=a*Q[++ii]+c; c=(t>>32);  ranx=(uint32_t)(t+c); if(ranx<c){ranx++;c++;} Q[ii]= ranx;
    // Turn it into a random integer in the range 0 to k-1 with (unsigned)((unsigned short)ranx * (k)) >> 16;
#define MWCrank(k) (unsigned)(((uint64_t)ranx * (uint64_t)(k)) >> 32)
    
    //Run off some random MWCs just to mix things up
    for (i=0; i<1000; i++) {
        MWCranx
    }
    
                  unsigned *dc, *d;
                  unsigned ncount = *rn;
                  d = Calloc(nreps * ncount , unsigned);
                  
                  dc = d;
                  // fill with the integers
                  for (i=0; i<nreps; i++) {
                      for (j=0; j<ncount; j++) {
                          dc[j] = j;
                      }
                      dc+= ncount;
                  }
                  
                  uint64_t nTrials = 0, nTodo = *rnTodo;
                  uint64_t nupper = 0, nlower = 0;
                  uint64_t * h = Calloc(*rhistobins, uint64_t);
                  unsigned temp, *wi, *wlast, randomRange;
                  
                  // scramble each of the NREPS sets
                  dc = d;
                  for (i=0; i<nreps; i++) {
                      for (j=0; j<ncount; j++) {
                          MWCranx r = MWCrank(ncount - j);
                          temp = dc[r+j];
                          dc[r+j] = dc[j];
                          dc[j] = temp;
                      }
                      dc += ncount;  // move pointer to the next set
                  }

                  unsigned *wc, *w; // = malloc(ncount * sizeof(ITYPE));
                  unsigned HN = *rhistobins;
                  
                  double statSpan = (HN-1)/(*rhistoRight - *rhistoLeft);
                  
                  w = Calloc(ncount, unsigned);
                  double sum1;
                  int hdex;

                  // Make a set to scramble. w contains index to x
                  for (i=0; i<ncount; i++) {
                      w[i] = i;
                  }
                  
                  wlast = &(w[*rk-1]); // so wlast points to the last of the x values, and before the y values
                  
                  // This is the main loop
                  while(nTrials < nTodo) {
                      // Scramble the first ncount1 values in w
                      randomRange = ncount;
                      for (wi = w; wi <= wlast; wi++) {
                          MWCranx r = MWCrank(randomRange--);
                          wc = wi + r;
                          temp = *wc;
                          *wc = *wi;
                          *wi = temp;
                      }
                      // get stats
                      for (i=0; i<nreps; i++) {
                          dc = d + i * ncount;  // dc now points to one of the scrambled sets of integers
                          sum1 = 0;
                          for (wi=w; wi<=wlast; wi++) {
                              sum1 += x[(dc[*wi])];
                          }
                          if (sum1 < *rstatLeft) nlower++;
                          if (sum1 > *rstatRight) nupper++;
                          hdex = (unsigned) (statSpan * (sum1 - *rhistoLeft));
                          if (hdex >= 0 && hdex < HN) h[hdex]++;
                          
                      } // for i
                      nTrials += nreps; // because each loop performs nreps trials
                      if (nTrials > nextCheck) { // it's time to check
                          nextCheck += checkGap;  // kick the can down the road
                          if(time(NULL) - start > *rsafeSecs) nTodo = nTrials;  // It's time to stop
                      }
                  } // while nTrials
                  
                  // Return the computed values to R
                  *rnupper = nupper;
                  *rnlower = nlower;
                  for (i=0; i<HN; i++) {
                      rhistoData[i] = h[i];
                  }
                  *subsetsExamined = nTrials;
                  *timespent = ((double)clock() - startclock)/CLOCKS_PER_SEC;
                  Free(d); Free(w); Free(h);
}
