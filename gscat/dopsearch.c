
/*
$Log: dopsearch.c,v $

Revision 1.0  2009/06/05 aj
this code reads a fit file and searches for plasmapause events

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
#include "rpos.h"
#include "fitdata.h"
#include "cfitdata.h"
#include "scandata.h"
#include "fitread.h"
#include "fitwrite.h"
#include "fitscan.h"
#include "fitindex.h"
#include "fitseek.h"
#include "oldfitread.h"
#include "oldfitscan.h"
#include "cfitread.h"
#include "cfitscan.h"
#include "fitscan.h"
#include "dopsearch_lib.h"
#include "aacgm.h"
#include "mlt.h"
#include "ribmath.h"
#include "scan2.h"


struct CFitdata cfit;
struct Scan src;
struct RadarNetwork *network;
struct Radar *radar;
struct RadarSite *site;
struct OptionData opt;

int nbox;
int ebmno=0;
int ebm[32*3];
int minrng=-1;
int maxrng=-1;


/*double myAbs(f)
{
  if(f < 0.)
    f *= -1.;
  return f;
}*/

int main(int argc,char *argv[])
{
  /*declarations*/
  extern int errno;
  int toprng, botrng, n_col, n_gray, n_soundings, col, gray, offset=0;
  double eventStart, eventEnd;
  struct coordStack myStack;
  int concat = 0, mode = 0, farg = 0, tgtcol, tgtgray,pflg,aflg;
  int startk, endk,last_k = 0;
  int i, j, k, t;
  double nsc, sc, jTime, myTime;
  int s = 0, new = 0, old = 0;
  int yr=0,mo=0,dy=0,hr,mt; 
  struct Beam myBeam;
  struct RadarParm tempprm,tempprmr;
  struct FitData tempfit,tempfitr;
  double *med_vel;
  int tgtbeam;
  FILE *fp;
  FILE *fitfp = NULL;
  int **w, gflg;  
  char *envstr;
  char *dname = NULL,*iname = NULL;
  unsigned char help = 0;
  int first = 1;
  unsigned char option = 0;
  struct scan * scans = malloc(1500*sizeof(struct scan));
  fprintf(stderr,"%ld\n",sizeof(struct scan));
  struct scan tempscan;
  if(scans == NULL)
    fprintf(stderr,"NO MEMORY for scans!!\n");
  double grho,glat,glon,gazm,mlat,mlon,mazm,myMlt,elv;
  int yr_sec;


  double **scatterp = malloc(2000000*sizeof(double *));
  if(scatterp == NULL)
    fprintf(stderr,"NO MEMORY for scatterp!!\n");
  for(i=0; i<2000000; i++)
  {
    scatterp[i]=malloc(2*sizeof(double));
    if(scatterp[i] == NULL)
      fprintf(stderr,"NO MEMORY!!\n");
  }


  double event_w, v_tot, w_l_tot;
  int p,q;



  struct eventNode *tgtscans=malloc(35000*sizeof(struct eventNode));
  if(tgtscans == NULL)
    fprintf(stderr,"NO MEMORY!!\n");
  int n_scns=0, tgt_scns=0;



  envstr=getenv("SD_RADAR");
  if(envstr==NULL)
  {
    fprintf(stderr,"Environment variable 'SD_RADAR' must be defined.\n");
    return -1;
  }
  fp=fopen(envstr,"r");

  if(fp==NULL) 
  {
    fprintf(stderr,"Could not locate radar information file.\n");
    exit(-1);
  }

  network=RadarLoad(fp);
  fclose(fp);
  if(network==NULL)
  {
    fprintf(stderr,"Failed to read radar information.\n");
    exit(-1);
  }

  envstr=getenv("SD_HDWPATH");
  if(envstr==NULL)
  {
    fprintf(stderr,"Environment variable 'SD_HDWPATH' must be defined.\n");
    exit(-1);
  }

  RadarLoadHardware(envstr,network);

  /*add options to control performance*/
  OptionAdd(&opt,"-help",'x',&help);
  OptionAdd(&opt,"-option",'x',&option);
  OptionAdd(&opt,"new",'x',&new);
  OptionAdd(&opt,"concat",'x',&concat);                 /*MAYBE*/
  OptionAdd(&opt,"year",'i',&yr);
  OptionAdd(&opt,"mon",'i',&mo);
  OptionAdd(&opt,"day",'i',&dy);
  OptionAdd(&opt,"beam",'i',&tgtbeam);

  /*process the options*/
  farg=OptionProcess(1,argc,argv,&opt,NULL);
  if(option == 1)
  {
    OptionDump(stdout,&opt);
    exit(0);
  }
  if(mode > 0)
        mode--;

  old = !new;

  double stEpoch = TimeYMDHMSToEpoch(yr,mo,dy,0,0,0);
  double edEpoch = TimeYMDHMSToEpoch(yr,mo,dy,23,59,59);



  if(argc-farg > 1)
  {
    dname=argv[argc-2];
    iname=argv[argc-1];
  }
  else
    dname=argv[argc-1];


  /*go through all of the day's scans
  and find beam 7 soundings*/
  for(tgtbeam=0;tgtbeam<30;tgtbeam++)
  {
    n_soundings=0;
    tgt_scns=0;
    fitfp=fopen(dname,"r");
    if(fitfp==NULL)
    {
      fprintf(stderr,"File not found.\n");
      exit(-1);
    }

    n_scns = 0;
    /*read the first entry in the fit file*/
    s=FitFread(fitfp,&tempprmr,&tempfitr);
    if(s==-1)
    {
      fprintf(stderr,"Error reading file\n");
      exit(-1);
    }
    s = fitReadtoScan(fitfp,&tempscan,0,&tempfitr,&tempprmr,7);
    radar=RadarGetRadar(network,tempscan.prm[0].stid);
    if (radar==NULL)
    {
      fprintf(stderr,"Failed to get radar information.\n");
      exit(-1);
    }
    site=RadarYMDHMSGetSite(radar,tempscan.prm[0].time.yr,tempscan.prm[0].time.mo,
                            tempscan.prm[0].time.dy,tempscan.prm[0].time.hr,
                            tempscan.prm[0].time.mt,tempscan.prm[0].time.sc);
    myTime = TimeYMDHMSToEpoch(tempscan.prm[0].time.yr,tempscan.prm[0].time.mo,
                            tempscan.prm[0].time.dy,tempscan.prm[0].time.hr,
                            tempscan.prm[0].time.mt,tempscan.prm[0].time.sc);


    do
    {
      /*go through each beam sounding in the scan*/
      for(i=0; i<tempscan.nbeams; i++)
      {
        /*the target beam sounding*/
        tempfit = tempscan.fit[i];
        tempprm = tempscan.prm[i];
        /*we have got our target beam*/
        if(tempprm.bmnum == tgtbeam && (tempprm.channel == 0 || tempprm.channel == 1))
        {
          n_soundings++;
          initEvent(&tgtscans[tgt_scns]);
          if(myTime >= stEpoch && myTime <= edEpoch) k = n_scns;
          else k = -1;
          copyBeam(&tgtscans[tgt_scns],tempfit,tempprm,k,i);
          tgtscans[tgt_scns].rangemax = -1;
          tgtscans[tgt_scns].rangemin = 300;
          /*go through range gates 0-75, assigning the beamsoundings'
          data to our new structure*/
          for(j=0; j<tempprm.nrang; j++)
          {
            tgtscans[tgt_scns].ranges[j].checked = 0;
            /*normalize velocity and scatter*/
            if(tgtscans[tgt_scns].beam.sct[j] == 0)
              tgtscans[tgt_scns].beam.rng[j].v = 0.0;
            else
              tgtscans[tgt_scns].beam.sct[j] = 1;
          }
          tgt_scns++;
        }
      }

      if(myTime >= stEpoch && myTime <= edEpoch && (tempprm.channel == 0 || tempprm.channel == 1))
      {
        if(first)
        {
          scans[n_scns] = tempscan;
          for(i=0; i<tempscan.nbeams; i++)
            for(j=0;j<tempscan.prm[i].nrang;j++)
            {
              scans[n_scns].fit[i].rng[j].gsct=0.;/*
              scans[n_scns].fit[i].rng[j].p_l=1.;*/
            }
        }
        n_scns++;
      }

      s = fitReadtoScan(fitfp,&tempscan,0,&tempfitr,&tempprmr,7);
      if(s != -1)
        myTime = TimeYMDHMSToEpoch(tempscan.prm[0].time.yr,tempscan.prm[0].time.mo,
                                tempscan.prm[0].time.dy,tempscan.prm[0].time.hr,
                                tempscan.prm[0].time.mt,tempscan.prm[0].time.sc);

    } while (s!=-1);

    fprintf(stderr,"n_scnas: %d\n",n_scns);

    fclose(fitfp);
    first = 0;

    w = malloc(tgt_scns*sizeof(int *));
    for(i=0; i<tgt_scns; i++)
      w[i] = malloc(150*sizeof(int));
    /*output some data to ascii files that will be used for plotting*/

    for(k=0; k<tgt_scns; k++)
    {
      myBeam = tgtscans[k].beam;
      for(j=5; j<myBeam.nrang; j++)
      {
        /*do the search*/
        if((myBeam.sct[j] == 1)&&(tgtscans[k].ranges[j].checked==0))
        {
          /*initialize some vars*/
          eventStart = 999999999999.;
          eventEnd = -99999999999.;
          startk = 999999;
          endk = -999999;
          toprng = -1;
          botrng = 250;
          n_col = 0;
          n_gray = 0;
          event_w = 0.;
          v_tot = 0.;
          w_l_tot = 0.;
          initStack(&myStack);
          for(i=0;i<35000;i++)
          {
            tgtscans[i].rangemin = 999;
            tgtscans[i].rangemax = -1;
          }
          for(i=0;i<2000000;i++)
            for(t=0;t<2;t++)
              scatterp[i][t] = 0.;
          /*actually do the search*/
          search_for_event(tgtscans, k, j, &myStack, &toprng, &botrng,
                            &eventStart, &eventEnd, &n_col, &n_gray,
                            n_soundings, &startk, &endk, &event_w, &v_tot, &w_l_tot,
                            scatterp);
          freeStack(&myStack);

          TimeEpochToYMDHMS(eventStart,&yr,&mo,&dy,&hr,&mt,&sc);
          int nhr, nmin;
          TimeEpochToYMDHMS(eventEnd,&yr,&mo,&dy,&nhr,&nmin,&nsc);

          med_vel = malloc((n_col+n_gray)*sizeof(double));
          if(med_vel == NULL)
            fprintf(stderr,"NO MEM FOR MED_vel\n");
          for(i=0;i<n_col+n_gray;i++)
            med_vel[i] = scatterp[i][0];
          qsort(med_vel, n_col+n_gray, sizeof(double), dbl_cmp);

          pflg = 0;
          aflg = 0;

          /*if we have found a plasmapause event*/
          if(1/*((eventEnd-eventStart > 9800)
              &&(toprng-botrng > 4)
              &&(((double)n_col/(double)n_gray) > .2))
              ||((eventEnd-eventStart > 7200) 
              &&(toprng-botrng > 4)
              &&(((double)n_col/(double)n_gray) > .33))
              ||((eventEnd-eventStart > 3499)
              &&(toprng-botrng > 4)
              &&(((double)n_col/(double)n_gray) > .475))*/)
          {
            /*try to eliminate leading ground scatter*/
            offset = 0;
            for(i=startk; i<endk; i++)
            {
              col = 0;
              gray = 0;
              tgtcol = 0;
              tgtgray = 0;
              for(t=tgtscans[i].rangemin; t<tgtscans[i].rangemax+1; t++)
              {
                if(tgtscans[i].beam.sct[t])
                {
                  if(ribfabs(tgtscans[i].beam.rng[t].v) < 15.)
                  {
                    gray++;
                    tgtgray++;
                  }
                  else
                  {
                    col++;
                    tgtcol++;
                  }
                }
              }
              if(i<tgt_scns-1)
              {
                for(t=tgtscans[i+1].rangemin; t<tgtscans[i+1].rangemax+1; t++)
                {
                  if((i<tgt_scns-1)&&(tgtscans[i+1].beam.sct[t]))
                  {
                    if(ribfabs(tgtscans[i+1].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
                }
              }
              if(i < tgt_scns-2)
              {
                for(t=tgtscans[i+2].rangemin; t<tgtscans[i+2].rangemax+1; t++)
                {
                  if(tgtscans[i+2].beam.sct[t])
                  {
                    if(ribfabs(tgtscans[i+2].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
                }
              }
              if(i<tgt_scns-3)
              {
                for(t=tgtscans[i+3].rangemin; t<tgtscans[i+3].rangemax+1; t++)
                {
                  if((i<tgt_scns-1)&&(tgtscans[i+3].beam.sct[t]))
                  {
                    if(ribfabs(tgtscans[i+3].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
                }
              }
              if(i < tgt_scns-4)
              {
                for(t=tgtscans[i+4].rangemin; t<tgtscans[i+4].rangemax+1; t++)
                {
                  if(tgtscans[i+4].beam.sct[t])
                  {
                    if(ribfabs(tgtscans[i+4].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
                }
              }
              if(((double)col/(double)gray > .75)||(col+gray < 10))
              {
                startk++;
                eventStart = tgtscans[i+1].beam.time;
                n_col -= tgtcol;
                n_gray -= tgtgray;
                if(n_gray < 1)
                  n_gray = 1;
              }
              else
                break;
            }

            /*try to eliminate following ground scatter*/
            offset = 0;
            for(i=endk; i>startk; i--)
            {
              col = 0;
              gray = 0;
              tgtcol = 0;
              tgtgray = 0;
              for(t=tgtscans[i].rangemin; t<tgtscans[i].rangemax+1; t++)
              {
                if(tgtscans[i].beam.sct[t])
                {
                  if(ribfabs(tgtscans[i].beam.rng[t].v) < 15.)
                  {
                    gray++;
                    tgtgray++;
                  }
                  else
                  {
                    col++;
                    tgtcol++;
                  }
                }
              }
              if(i > 0)
              {
                for(t=tgtscans[i-1].rangemin; t<tgtscans[i-1].rangemax+1; t++)
                  if(tgtscans[i-1].beam.sct[t])
                  {
                    if(ribfabs(tgtscans[i-1].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
              }
              if(i>1)
              {
                for(t=tgtscans[i-2].rangemin; t<tgtscans[i-2].rangemax+1; t++)
                  if(tgtscans[i-2].beam.sct[t])
                  {
                    if(ribfabs(tgtscans[i-2].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
              }
              if(i > 2)
              {
                for(t=tgtscans[i-3].rangemin; t<tgtscans[i-3].rangemax+1; t++)
                  if(tgtscans[i-3].beam.sct[t])
                  {
                    if(ribfabs(tgtscans[i-3].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
              }
              if(i>3)
              {
                for(t=tgtscans[i-4].rangemin; t<tgtscans[i-4].rangemax+1; t++)
                  if(tgtscans[i-4].beam.sct[t])
                  {
                    if(ribfabs(tgtscans[i-4].beam.rng[t].v) < 15.)
                      gray++;
                    else
                      col++;
                  }
              }

              if(((double)col/(double)gray > .75) || (col+gray < 10))
              {
                endk--;
                eventEnd = tgtscans[i-1].beam.time;
                n_col -= tgtcol;
                n_gray -= tgtgray;
                if(n_gray < 1)
                  n_gray = 1;
              }
              else
                break;
            }

            int cnt = 0;
            for(q=startk; q<endk; q++)
              for(p=tgtscans[q].rangemin; p<tgtscans[q].rangemax+1; p++)
              {
                myBeam = tgtscans[q].beam;
                if(myBeam.sct[p])
                {
                  scatterp[cnt][0] = myBeam.rng[p].v;
                  cnt++;
                }
              }

            free(med_vel);
            med_vel = malloc((cnt)*sizeof(double));
            if(med_vel == NULL)
              fprintf(stderr,"NO MEM FOR MED_vel\n");
            for(i=0;i<cnt;i++)
              med_vel[i] = scatterp[i][0];
            qsort(med_vel, cnt, sizeof(double), dbl_cmp);

						gflg = 0;

					if(((eventEnd-eventStart > 32400)
            &&(toprng-botrng > 4)
            &&(((double)n_gray/(double)n_col) > 4.)))
						gflg = 1;
          else
              continue;

            last_k = endk;
            if(gflg)
            {
              TimeEpochToYMDHMS(eventStart,&yr,&mo,&dy,&hr,&mt,&sc);
              int nhr, nmin;
              TimeEpochToYMDHMS(eventEnd,&yr,&mo,&dy,&nhr,&nmin,&nsc);
              fprintf(stderr,"a gs event found!  ");
              fprintf(stderr,"%d:%d - %d:%d  ranges %d-%d  %d\n",
                            hr, mt, nhr, nmin, botrng, toprng,tgtbeam);

              for(q=startk; q<endk; q++)
              {
                for(p=tgtscans[q].rangemin; p<tgtscans[q].rangemax+1; p++)
                {
                  myBeam = tgtscans[q].beam;
                  if(myBeam.sct[p] && myBeam.scanind != -1)
                  {
                    scans[myBeam.scanind].fit[myBeam.beamind].rng[p].gsct=1.;
                    TimeEpochToYMDHMS(myBeam.time,&yr,&mo,&dy,&hr,&mt,&sc);
                    RPosGeo(1,myBeam.bm,p,site,myBeam.frang,myBeam.rsep,myBeam.rxrise,200.,&grho,&glat,&glon);
                    RPosRngBmAzmElv(myBeam.bm,p,yr,site,myBeam.frang,myBeam.rsep,myBeam.rxrise,200.,&gazm,&elv);
                    RPosInvMag(myBeam.bm,p,yr,site,myBeam.frang,myBeam.rsep,myBeam.rxrise,200.,&mlat,&mlon,&mazm);
                    yr_sec = TimeYMDHMSToYrsec(yr,mo,dy,hr,mt,sc);
                    myMlt = AACGMConvertMLT(yr,yr_sec,mlon);
                    /*fprintf(stdout,"%d  %d  %d  %d  %d  %d  %d  %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
                                    myBeam.bm,p,yr,mo,dy,hr,mt,(int)sc,myBeam.time,
                                    myBeam.rng[p].v,myBeam.rng[p].w_l,glat,glon,gazm,
                                    mlat,mlon,mazm,myMlt);
                     scans[myBeam.scanind].fit[myBeam.beamind].rng[p].p_l=17.;*/
                  }
                }
              }
            }
          }
          free(med_vel);
        }
      }
    }
    free(w);
  }

  for(k=0; k<=n_scns; k++)
    for(j=0; j<scans[k].nbeams; j++)
    {
      myTime = TimeYMDHMSToEpoch(scans[k].prm[j].time.yr,scans[k].prm[j].time.mo,
                          scans[k].prm[j].time.dy,scans[k].prm[j].time.hr,
                          scans[k].prm[j].time.mt,scans[k].prm[j].time.sc);
      /*
        if(myTime >= stEpoch && myTime <= edEpoch)
      */
        s=FitFwrite(stdout,&scans[k].prm[j],&scans[k].fit[j]);
    }

  free(scatterp);
  free(scans);
  free(tgtscans);
  return 0;

}
