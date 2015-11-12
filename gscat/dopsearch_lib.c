/*
$Log: dopsearch_lib.c,v $

Revision 1.0  2009/06/05 aj
implementation for various functions

*/

#include <errno.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
#include "rtypes.h"
#include "dmap.h"
#include "option.h"
#include "rtime.h"
#include "radar.h"
#include "limit.h"
#include "rprm.h"
#include "invmag.h"
#include "fitdata.h"
#include "cfitdata.h"
#include "scandata.h"
#include "fitread.h"
#include "fitscan.h"
#include "fitindex.h"
#include "fitseek.h"
#include "oldfitread.h"
#include "oldfitscan.h"
#include "cfitread.h"
#include "cfitscan.h"
#include "fitscan.h"
#include "dopsearch_lib.h"
#include "ribmath.h"

/*
double ribfabs(double a)
{
  if(a < 0)
    a = (a)*(-1.);

  return a;
}
*/
void initEvent(struct eventNode *myNode)
{
  myNode->ranges = malloc(150*sizeof(int));
  myNode->rangemax = -1;
  myNode->rangemin = 999;
}

void initStack(struct coordStack *myStack)
{
  if(myStack == NULL)
  {
    myStack = malloc(sizeof(struct coordStack));
  }
  myStack->time = malloc(250000*sizeof(int));
  myStack->rng = malloc(250000*sizeof(int));
  if (myStack->time == NULL)
  {
    fprintf(stderr, "Insufficient memory to initialize stack.\n");
    exit(1);
  }
  myStack->top = -1;
}

void freeStack(struct coordStack *myStack)
{
  /* Get rid of array. */
  free(myStack->time);
  free(myStack->rng);

  myStack->time = NULL;
  myStack->rng = NULL;
  myStack->top = -1; 
}


void push(struct coordStack *myStack, int time, int rng)
{
  /* Put information in array*/
  myStack->top++;
  myStack->time[myStack->top] = time;
  myStack->rng[myStack->top] = rng;
}


void pop(struct coordStack *myStack)
{
  myStack->top--;
}

/*this function was written by RJ Barnes and Kile Baker
int dbl_cmp(const void *x,const void *y) {
  double *a,*b;
  a=(double *) x;
  b=(double *) y;
  if (*a > *b) return 1;
  else if (*a == *b) return 0;
  else return -1;
}
*/   
/*double stdDev(double *arr, int n)
{
  double mean=0., sum=0., dev=0.;
  int i=0;
  for(i=0; i<n; i++) 
  {
    sum += arr[i];
  }
  if(sum == 0.)
    return 0;
  else
    mean = sum/((double)n);

  sum=0.;
  for(i=0; i<n; i++)
  {
    sum += (arr[i]-mean)*(arr[i]-mean);
  }
  dev = (double)(sum/((double)n));

  return dev;
}*/
void do_search(struct eventNode *myEvent, int k, int j, struct coordStack *myStack,
                int *toprng, int *botrng, double *eventStart, double *eventEnd,
                int *n_col, int *n_gray, int n_soundings, int *startk,
                int *endk, double *event_w, double *v_tot, double *w_l_tot,
                double **scatterp)
{
  int temp;
  struct Beam myBeam = myEvent[k].beam;

  if(j < myEvent[k].rangemin)
    myEvent[k].rangemin = j;
  if(j > myEvent[k].rangemax)
    myEvent[k].rangemax = j;

  if(myEvent[k].ranges[j].checked == 0)
  {
    scatterp[*n_col+*n_gray][0] = myBeam.rng[j].v;
    scatterp[*n_col+*n_gray][1] = myBeam.rng[j].w_l;

    if(ribfabs(myBeam.rng[j].v) >= 15.)
    {
      temp = *n_col + 1;
      *n_col = temp;
    }
    else
    {
      temp = *n_gray + 1;
      *n_gray = temp;
    }
    if(myBeam.cpid == 3310)
    {
      *event_w += .125;
      *v_tot += myBeam.rng[j].v * .125;
      *w_l_tot += myBeam.rng[j].w_l * .125;
    }
    else
    {
      *event_w += 1.;
      *v_tot += myBeam.rng[j].v;
      *w_l_tot += myBeam.rng[j].w_l;
    }
  }
  myEvent[k].ranges[j].checked = 1;
  /*the target beam sounding*/

  if(myBeam.time < *eventStart)
  {
    *eventStart = myBeam.time;
    *startk     = k;
  }
  if(myBeam.time > *eventEnd) 
  {
    *eventEnd = myBeam.time;
    *endk     = k;
  }
  if(j < *botrng)
    *botrng = j;
  if(j > *toprng)
    *toprng = j;

  int ind = 1;
  int time_lim = 361.;

  while((k>ind)&&(k<n_soundings-ind-1)&&(ribfabs(myEvent[k+ind].beam.time - myEvent[k].beam.time) < time_lim)
      &&(ribfabs(myEvent[k-ind].beam.time - myEvent[k].beam.time) < time_lim))
  {
    if((k>ind+1)&&(myEvent[k-ind].beam.sct[j] != 0)
        &&(myEvent[k-ind].ranges[j].checked==0))
    {
      push(myStack,k-ind,j);
      return;
    }
    if((myEvent[k].beam.sct[j-1] != 0)&&(myEvent[k].ranges[j-1].checked==0)
            &&(j>6))
    {
      push(myStack,k,j-1);
      return;
    }
    if((k<n_soundings-ind-1)&&(myEvent[k+ind].beam.sct[j] != 0)
        &&(myEvent[k+ind].ranges[j].checked==0))
    {
      push(myStack,k+ind,j);
      return;
    }
    if((myEvent[k].beam.sct[j+1] != 0)&&(myEvent[k].ranges[j+1].checked==0)
            &&(j<myBeam.nrang))
    {
      push(myStack,k,j+1);
      return;
    }
    if((k>ind+1)&&(j>6)&&(myEvent[k-ind].beam.sct[j-1] != 0)
        &&(myEvent[k-ind].ranges[j-1].checked==0))
    {
      push(myStack,k-ind,j-1);
      return;
    }
    if((k>ind+1)&&(j<myBeam.nrang)&&(myEvent[k-ind].beam.sct[j+1] != 0)
        &&(myEvent[k-ind].ranges[j+1].checked==0))
    {
      push(myStack,k-ind,j+1);
      return; 
    }
    if((k<n_soundings-ind-1)&&(j<myBeam.nrang)&&(myEvent[k+ind].beam.sct[j+1] != 0)
        &&(myEvent[k+ind].ranges[j+1].checked==0))
    {
      push(myStack,k+ind,j+1);
      return;
    }
    if((k<n_soundings-ind-1)&&(j>6)&&(myEvent[k+ind].beam.sct[j-1] != 0)
        &&(myEvent[k+ind].ranges[j-1].checked==0))
    {
      push(myStack,k+ind,j-1);
      return;
    }
    ind++;

  }

  pop(myStack);
  return;
}

void search_for_event(struct eventNode *myEvent, int k, int j, struct coordStack *myStack,
                int *toprng, int *botrng, double *eventStart, double *eventEnd,
                int *n_col, int *n_gray, int n_soundings, int *startk,
                int *endk, double *event_w, double *v_tot, double *w_l_tot,
                double **scatterp)
{
  push(myStack,k,j);
  int temp1,temp2, temp3, temp4;
  while(myStack->top > -1)
  {
    temp1 = myStack->time[myStack->top];
    temp3 = myStack->rng[myStack->top];
    do_search(myEvent,myStack->time[myStack->top],myStack->rng[myStack->top],
                myStack, toprng, botrng, eventStart, eventEnd,
                n_col, n_gray, n_soundings, startk, endk, event_w, v_tot, w_l_tot,
                scatterp);
    temp2 = myStack->time[myStack->top];
    temp4 = myStack->rng[myStack->top];
    /*if(((ribfabs(temp1-temp2) > 3)||(ribfabs(temp3-temp4) > 2))&&(myStack->top > -1))
    {
      fprintf(stderr, "we have a problem\n");
      exit(-1);
    }*/

  }
  return;
}

void initBeam(struct FitData *fit, struct RadarParm *prm, struct Scan *myScan, int num)
{
  int i,j;
  myScan->stid = prm[0].stid;

  myScan->st_time = TimeYMDHMSToEpoch(prm[0].time.yr,prm[0].time.mo,prm[0].time.dy,
                                      prm[0].time.hr,prm[0].time.mt,(double)prm[0].time.sc);
  myScan->ed_time = TimeYMDHMSToEpoch(prm[num-1].time.yr,prm[num-1].time.mo,prm[num-1].time.dy,
                                      prm[num-1].time.hr,prm[num-1].time.mt,(double)prm[num-1].time.sc);
  myScan->num = num;
  for(i=0;i<num;i++)
  {
    myScan->bm[i].scan = prm[i].scan;
    myScan->bm[i].bm = prm[i].bmnum;
    myScan->bm[i].bmazm = prm[i].bmazm;
    myScan->bm[i].time = TimeYMDHMSToEpoch(prm[i].time.yr,prm[i].time.mo,prm[i].time.dy,
                                      prm[i].time.hr,prm[i].time.mt,(double)prm[i].time.sc);
    myScan->bm[i].cpid = prm[i].cp;
    myScan->bm[i].intt.sc = prm[i].intt.sc;
    myScan->bm[i].intt.us = prm[i].intt.us;
    myScan->bm[i].nave = prm[i].nave;
    myScan->bm[i].frang = prm[i].frang;
    myScan->bm[i].rsep = prm[i].rsep;
    myScan->bm[i].rxrise = prm[i].rxrise;
    myScan->bm[i].freq = prm[i].tfreq;
    myScan->bm[i].noise = prm[i].noise.search;
    myScan->bm[i].atten = prm[i].atten;
    myScan->bm[i].channel = prm[i].channel;
    myScan->bm[i].nrang = prm[i].nrang;
    for(j=0;j<prm[i].nrang;j++)
    {
      myScan->bm[i].sct[j] = fit[i].rng[j].qflg;
      myScan->bm[i].rng[j].gsct = fit[i].rng[j].gsct;
      myScan->bm[i].rng[j].p_0 = fit[i].rng[j].p_0;
      myScan->bm[i].rng[j].v = fit[i].rng[j].v;
      myScan->bm[i].rng[j].v_e = fit[i].rng[j].v_err;
      myScan->bm[i].rng[j].w_l = fit[i].rng[j].w_l;
      myScan->bm[i].rng[j].w_l_e = fit[i].rng[j].w_l_err;
      myScan->bm[i].rng[j].p_l = fit[i].rng[j].p_l;
      myScan->bm[i].rng[j].p_l_e = fit[i].rng[j].p_l_err;
      myScan->bm[i].rng[j].elv = fit[i].elv[j].normal;
      myScan->bm[i].rng[j].phi0 = fit[i].rng[j].phi0;
      myScan->bm[i].rng[j].out = 0;
    }
  }
}

void copyBeam(struct eventNode * myNode, struct FitData myFit, struct RadarParm myPrm,
          int scanind, int beamind)
{
  int i;

  myNode->beam.scan = myPrm.scan;
  myNode->beam.bm = myPrm.bmnum;
  myNode->beam.bmazm = myPrm.bmazm;
  myNode->beam.time = TimeYMDHMSToEpoch(myPrm.time.yr,myPrm.time.mo,myPrm.time.dy,
                                      myPrm.time.hr,myPrm.time.mt,(double)myPrm.time.sc);
  myNode->beam.cpid = myPrm.cp;
  myNode->beam.intt.sc = myPrm.intt.sc;
  myNode->beam.intt.us = myPrm.intt.us;
  myNode->beam.nave = myPrm.nave;
  myNode->beam.frang = myPrm.frang;
  myNode->beam.rsep = myPrm.rsep;
  myNode->beam.rxrise = myPrm.rxrise;
  myNode->beam.freq = myPrm.tfreq;
  myNode->beam.noise = myPrm.noise.search;
  myNode->beam.atten = myPrm.atten;
  myNode->beam.channel = myPrm.channel;
  myNode->beam.nrang = myPrm.nrang;
  myNode->beam.scanind = scanind;
  myNode->beam.beamind = beamind;
  for(i=0;i<myPrm.nrang;i++)
  {
    myNode->beam.sct[i] = myFit.rng[i].qflg;
    myNode->beam.rng[i].v = myFit.rng[i].v;
    myNode->beam.rng[i].p_0 = myFit.rng[i].p_0;
    myNode->beam.rng[i].w_l = myFit.rng[i].w_l;
    myNode->beam.rng[i].p_l = myFit.rng[i].p_l;
    myNode->beam.rng[i].phi0 = myFit.rng[i].phi0;
  }

}

void initFit(struct FitData *fit, struct RadarParm *prm, struct Beam myBeam)
{
  int j,yr,mo,dy,hr,mt;
  double sc;

  prm->stid = 33;

  prm->scan = myBeam.scan;
  prm->bmnum = myBeam.bm;
  prm->bmazm = myBeam.bmazm;
  TimeEpochToYMDHMS(myBeam.time,&yr,&mo,&dy,&hr,&mt,&sc);
  prm->time.yr = yr;
  prm->time.mo = mo;
  prm->time.dy = dy;
  prm->time.hr = hr;
  prm->time.mt = mt; 
  prm->time.sc = sc;
  prm->cp = myBeam.cpid;
  prm->intt.sc = myBeam.intt.sc;
  prm->intt.us = myBeam.intt.us;
  prm->nave = myBeam.nave;
  prm->frang = myBeam.frang;
  prm->rsep = myBeam.rsep; 
  prm->rxrise = myBeam.rxrise;
  prm->tfreq = myBeam.freq;
  prm->noise.search = myBeam.noise;
  prm->atten = myBeam.atten;
  prm->channel = myBeam.channel;
  prm->nrang = myBeam.nrang;
  for(j=0;j<prm->nrang;j++)
  {
    fit->rng[j].qflg = myBeam.sct[j];
    fit->rng[j].gsct = myBeam.rng[j].gsct;
    fit->rng[j].p_0 = myBeam.rng[j].p_0;
    fit->rng[j].v = myBeam.rng[j].v;
    fit->rng[j].v_err = myBeam.rng[j].v_e;
    fit->rng[j].w_l = myBeam.rng[j].w_l;
    fit->rng[j].w_l_err = myBeam.rng[j].w_l_e;
    fit->rng[j].p_l = myBeam.rng[j].p_l;
    fit->rng[j].p_l_err = myBeam.rng[j].p_l_e;
    fit->elv[j].normal = myBeam.rng[j].elv;
    fit->rng[j].phi0 = myBeam.rng[j].phi0;
  }

}
