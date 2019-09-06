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

// globals for recursive diff
extern double * x;
extern int counter;
extern double * x256;
extern int npairs;

/*
 This recursive function computes the 256 possible sums of differences for the
 8 differences starting with the initial k
 */
//-(void) addDif: (int) k stat: (double) sum
void addDifr(int k, double sum)
{
    double xk = x[k];
    if(k >= npairs) xk = 0; // same as padding x with zeroes.
    if((k+1)%8) {
        addDifr(k+1, sum);
        addDifr(k+1, sum + xk);
    } else {
        x256[counter ++] = sum * 2;
        x256[counter ++] = (sum + xk) * 2;
    }
}


void rpairTest(
              double * rx,  // n differences
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
              ) {
                  int i;
                  x = rx;
                  start = time(NULL);
                  int timeLimit = *rsafeSecs;
                  uint64_t checkGap = 1e4;
                  uint64_t nextCheck = checkGap;
                  int HN;
                  double histoLeft, histoRight, statLeft, statRight, statSpan;
                  int64_t ntrials;
                  uint64_t * hdata;
                  uint64_t nlower, nupper;
                  unsigned hdex;
                  double sum1;
                  int ng8, ng32;
                  unsigned remaindr, remainderGap;
                  int isCancelled;
                  double * x256last;
    
                  
                  
    // for MWC256
    uint32_t Q[256];
    uint64_t t, c = 362436;
    uint8_t ii = 255;
    Q[0] = 0;
    uint64_t a=1540315826LL;
    uint32_t ranx;
    unsigned rseed = (unsigned)(rhistoData[4] + .5); // optional use of seed passed in rhistoData array

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
 
    
                  // Set up globals
                  HN = *rhistobins;
                  timeLimit = *rsafeSecs;
                  
                  statLeft = *rstatLeft;
                  statRight = *rstatRight;
                  double maxsum = 0;
                  for(int i=0;i<(*rn); i++) maxsum += rx[i];
                  histoLeft = *rhistoLeft;
                  histoRight = *rhistoRight;
                  statSpan = (HN-1)/(histoRight - histoLeft);
                  x = rx;
                  ntrials = 0;
                  uint64_t nTodo = *rnTodo;
                  nupper = 0;
                  nlower = 0;
                  checkGap = 1e6;
                  nextCheck = checkGap;
                  npairs = *rn;
                  ng8 = (npairs + 7)/8;
                  ng32 = (ng8 + 3)/4;
                  remaindr =  8- 8 * ng8 + *rn;
                  remainderGap = R_pow_di(2, 8-remaindr);
                  counter = 0;
                  double minsum = 0;
                  for(int i=0; i<npairs; i++) minsum -= fabs(x[i]);
                  hdata = Calloc(HN + 1, uint64_t);
                  for (int i = 0; i < HN; i++) hdata[i] = 0;
                  x256 = Calloc(ng8 * 256, double);
                  uint32_t * ranxx;
                  ranxx = Calloc(ng32, uint32_t);
                  uint8_t * ranByte;
                  ranByte = ranxx; // so we can treat the ng32 vals as ng8 bytes
                  for(int i=0; i< ng8*256; i++) x256[i] = 0;
                  start = time(NULL);
                  time_t startclock = clock();

                  // Set up x256 array using a recursive routine.
                  //This array pre-calculates all 256 possible sums for each group of 8 x's
                  counter = 0;
                  for(int i=0; i<ng8; i++)  {
                      addDifr(8 * i, 0);
                  }
                  
                  //Rprintf("\nHere are the difs, recursively computed:\n");
                  // for(int i=0; i< ng8*256; i++) Rprintf("%f, ", x256[i]);
                  
                  isCancelled = 0;
                  x256last = & (x256[256 * (ng8-1)]);
                  
                  // Start Main Loop
                  while(ntrials < nTodo) {
                      for(i=0; i<ng32; i++) {
                          MWCranx
                          ranxx[i] = ranx;
                      }
                      sum1 = minsum; //*rhistoLeft;
                      ranByte = ranxx;
                      for(int i = 0; i<ng8; i++) sum1 += x256[i*256 + ranByte[i]];
                      if (sum1 < *rstatLeft) nlower++;
                      if (sum1 > *rstatRight) nupper++;
                      hdex = (unsigned) (statSpan * (sum1 - *rhistoLeft));
                      if (hdex >= 0 && hdex < HN) hdata[hdex]++;
 //                     if(ntrials < 10) Rprintf("\ntrial %lld: sum1= %f : hdex = %d", ntrials, sum1, hdex);
//                      Rprintf("\n\nTrial %lld   : ", ntrials);
//                      for(i=0; i<ng8; i++) Rprintf("%d, ", ranByte[i]);
                      
                      
                      ntrials++;
                      if (ntrials > nextCheck){// It's time to check
                          nextCheck += checkGap;
                          if(time(NULL) - start >= *rsafeSecs) {//Time to stop
                              nTodo = ntrials;
                              Rprintf("\nProcess timed out\n");
                      }
                    }
                  }
                  
                  
                  *rnlower = nlower;
                  *rnupper = nupper;
                  *subsetsExamined = ntrials;
                  for (int i = 0; i<HN; i++) {
                      rhistoData[i] = hdata[i];
                  }
                  *timespent = ((double)clock() - startclock)/CLOCKS_PER_SEC;
                  
                  Free(x256);
                  Free(ranxx);
                  Free(hdata);

}
