#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "rna1_SSA.h"

#define NSPE 3 /* Number of species */
#define NPAR 4 /* Number of parameters */
#define NRXN 4 /* Number of reactions */
#define NTPS 2 /* Number of time points */
#define NREP 50 /* Number of replicas */

int main () {
  double x0[NSPE]; /* Initial State */
  double xt0[NSPE]; /* Temporal initial state */
  double p[NPAR]; /* Parameters */
  double tps[NTPS]; /* Time points */

  double maxE = 10000000000; /* Maximal number of events */
  double maxt = 100000; /* Maximal time */
  double rd = 0;
  int i = 0, qd = 0;
  FILE *ofp;  char outputF[255];

  /* Time points */
  tps[0] =   0.00000e+00;
  tps[1] =   1.00000e+05;
  x0[0] =   1.00000e+01; /* G0 */
  x0[1] =   0.00000e+00; /* G1 */
  x0[2] =   0.00000e+00; /* M */
  p[0] =   1.00000e-02; /* k01 */
  p[1] =   1.00000e-02; /* k10 */
  p[2] =   1.00000e-02; /* km */
  p[3] =   1.00000e-04; /* gm */
  /* Set random seed */
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned char *po = (unsigned char *)&tv;
  unsigned seed = 0;
  size_t z;
  for (z = 0; z < sizeof tv; z++){
    seed = seed * (UCHAR_MAX + 2U) + po[z];
  }
  srand (seed);

  sprintf(outputF,"rna1_doc5_scan_ss.out");
  ofp = fopen(outputF, "w");
  for(i=0;i<NREP;i++){
    SSAs(maxE,maxt,x0,p,NTPS,tps,ofp);
   }
  fclose(ofp);
}
