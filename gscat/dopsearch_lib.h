/*
$Log: dopsearch_lib.h,v $

Revision 1.0  2009/06/05  aj
definitions for various functions in dopsearch_lib

*/

struct coordStack{
  int *time;
  int *rng;
  int top;
}; 

struct RangeCell {
  int gsct;
  double p_0;
  double p_0_e;
  double v;
  double v_e;
  double w_l;
  double w_l_e;
  double p_l;
  double p_l_e;
  double elv;
  double phi0;
  short out;
};

struct Beam {
  int scan;
  int bm;
  float bmazm;
  double time;
  int cpid;
  struct {
    int sc;
    int us;
  } intt;
  int nave;
  int frang;
  int rsep;
  int rxrise;
  int freq;
  int noise;
  int atten;
  int channel;
  int nrang;
  unsigned char sct[150];
  struct RangeCell rng[150];
  int beamind;
  int scanind;
};

struct Scan {
  int stid;
  struct {
    int major;
    int minor;
  } version;

  double st_time;
  double ed_time;
  int num;
  struct Beam bm[50];
};

struct cellNode{
  int checked;
};

struct eventNode{
  struct Beam beam;
  struct cellNode *ranges;
  int rangemax;
  int rangemin;
};

void initEvent(struct eventNode *myNode);
void initStack(struct coordStack *myStack);
void freeStack(struct coordStack *myStack);
void push(struct coordStack *myStack, int time, int rng);
void pop(struct coordStack *myStack);
int dbl_cmp(const void *x,const void *y);
double stdDev(double *arr, int n);
void do_search(struct eventNode *myEvent, int k, int j, struct coordStack *stack,
                int *toprng, int *botrng, double *eventStart, double *eventEnd,
                int *n_col, int *n_gray, int n_soundings, int *startk,
                int *endk, double *event_w, double *v_tot, double *w_l_tot,
                double **scatterp);
void search_for_event(struct eventNode *myEvent, int k, int j, struct coordStack *stack,
                int *toprng, int *botrng, double *eventStart, double *eventEnd,
                int *n_col, int *n_gray, int n_soundings, int *startk,
                int *endk, double *event_w, double *v_tot, double *w_l_tot,
                double **scatterp);
void initBeam(struct FitData *fit, struct RadarParm *prm, struct Scan *myScan, int num);
void initFit(struct FitData *fit, struct RadarParm *prm, struct Beam myBeam);
void copyBeam(struct eventNode * myNode, struct FitData myFit, struct RadarParm myPrm,
              int scanind, int beamind);

