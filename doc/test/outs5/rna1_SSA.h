#define NSPE 3 /* Number of species */
#define NPAR 4 /* Number of parameters */
#define NRXN 4 /* Number of reactions */

void v_mltadd(int n, double x1[], double x2[], double s, double x3[]){
 /* x3 <- x1 + s*x2 */
  int q = 0;
  for (q=0;q<n;q++){
    x3[q] = x1[q] + s*x2[q];
  }
}
double v_sum (int n, double x[]){
 /* Sums the elements of x up to n */
  int q = 0;
  double sum = 0;
  for (q=0;q<n;q++){
    sum += x[q];
  }
  return sum;
}
void SSAs (double maxE, double maxt, double x0[], double p[], int NTPS, double tps[], FILE *s) {
  /* Set random seed */
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned char *po = (unsigned char *)&tv;
  unsigned seed = 0;
  size_t i;
  for (i = 0; i < sizeof tv; i++){
    seed = seed * (UCHAR_MAX + 2U) + po[i];
  }
  srand (seed);

  
  /* Variables declarations */
  double x[NSPE]; /* Current state */
  double R[NRXN]; /* Current Propensities */
  double event = 0; /* Number of events */
  double t = 0; /* Current time */
  double rsum = 0; /* For sum of propensities */
  double lsum = 0; /* For cumsum of propensities */
  double r1 = 0; /* Random number */
  double r2 = 0; /* Random number */
  double dt = 0; /* Time updates */
  int j, q;
  int ct = 0;
  /* Initialise values */
  for(q=0;q<NSPE;q++){
    x[q] = x0[q];
  }
            R[0] = p[0]*x[0]; /* k01*G0 */
            R[1] = p[1]*x[1]; /* k10*G1 */
            R[2] = p[2]*x[1]; /* km*G1 */
            R[3] = p[3]*x[2]; /* gm*M */

  do {
    rsum = v_sum(NRXN,R); /* sum of propensities*/
    r2 = rand() / (double) RAND_MAX; /* rand, 0 and 1*/
    dt = log(1/r2)/rsum; /* Increase in time */
    t += dt; /* time update */
    while ((NTPS > ct) && (tps[ct] <= t)){
      fprintf(s,"%e",tps[ct]);
      for (q=0;q<NSPE;q++){
        fprintf(s," %d",(int)x[q]);
      }
      fprintf(s,"\n");
      ct++;
    }
    if (t<maxt){
      r1 = (rand() / (double) RAND_MAX)*rsum;
      lsum = 0;
      for (j=1;j<=NRXN;j++){
        lsum += R[j-1];
        if (lsum >= r1){
          switch(j){
          case 1: /* G0 => G1 */
            x[0] -= 1; /* G0 */
            x[1] += 1; /* G1 */
            R[0] = p[0]*x[0]; /* k01*G0 */
            R[1] = p[1]*x[1]; /* k10*G1 */
            R[2] = p[2]*x[1]; /* km*G1 */
            break;
          case 2: /* G1 => G0 */
            x[0] += 1; /* G0 */
            x[1] -= 1; /* G1 */
            R[0] = p[0]*x[0]; /* k01*G0 */
            R[1] = p[1]*x[1]; /* k10*G1 */
            R[2] = p[2]*x[1]; /* km*G1 */
            break;
          case 3: /* G1 => G1 + M */
            x[2] += 1; /* M */
            R[3] = p[3]*x[2]; /* gm*M */
            break;
          case 4: /* M =>  */
            x[2] -= 1; /* M */
            R[3] = p[3]*x[2]; /* gm*M */
            break;
          default:
            break;
          }
          break;
	}
      }
      event++;
    }
  }while((maxE>event) && (maxt>t));
  return;
}

void SSAt (double maxE, double maxt, double x0[], double p[], int NTPS, double tps[], FILE *s) {
  /* Set random seed */
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned char *po = (unsigned char *)&tv;
  unsigned seed = 0;
  size_t i;
  for (i = 0; i < sizeof tv; i++){
    seed = seed * (UCHAR_MAX + 2U) + po[i];
  }
  srand (seed);

  
  /* Variables declarations */
  double x[NSPE]; /* Current state */
  double R[NRXN]; /* Current Propensities */
  double event = 0; /* Number of events */
  double t = 0; /* Current time */
  double rsum = 0; /* For sum of propensities */
  double lsum = 0; /* For cumsum of propensities */
  double r1 = 0; /* Random number */
  double r2 = 0; /* Random number */
  double dt = 0; /* Time updates */
  int j, q;
  int ct = 0;
  /* Initialise values */
  for(q=0;q<NSPE;q++){
    x[q] = x0[q];
  }
            R[0] = p[0]*x[0]; /* k01*G0 */
            R[1] = p[1]*x[1]; /* k10*G1 */
            R[2] = p[2]*x[1]; /* km*G1 */
            R[3] = p[3]*x[2]; /* gm*M */
  double x1[NSPE]; /* State first moments */
  double x2[NSPE]; /* State second moments */
  double xx[NSPE];
  double it = 0; /* InterIntervals time */
  do {
    rsum = v_sum(NRXN,R); /* sum of propensities*/
    r2 = rand() / (double) RAND_MAX; /* rand, 0 and 1*/
    dt = log(1/r2)/rsum; /* Increase in time */
    t += dt; /* time update */
    /* Average trajectories over time intervals */
    if ((NTPS > ct) && (tps[ct] <= t)){
      /* Add reminder */
      it += dt-(t-tps[ct]);
      v_mltadd(NSPE,x1,x,dt-(t-tps[ct]),x1); /* First moment */
      for (q=0;q<NSPE;q++){
         xx[q] = x[q]*x[q];
      }
      v_mltadd(NSPE,x2,xx,dt,x2); /* Second moment */
      if (it>0){ /* Normalise by interval length*/
       for (q=0;q<NSPE;q++){
         x1[q] = x1[q]/it;
         x2[q] = x2[q]/it;
       }
      }else{
       for (q=0;q<NSPE;q++){
         x1[q] = x[q];
         x2[q] = xx[q];
       }
      }
      fprintf(s,"%e %e",tps[ct]-it,tps[ct]); /* Print initial and final times*/
      for (q=0;q<NSPE;q++){
        fprintf(s," %e",x1[q]); /* Print first moments */
      }
      for (q=0;q<NSPE;q++){
        fprintf(s," %e",x2[q]); /* Print second moments */
      }
      fprintf(s,"\n");
      ct++;
      /* Possibly we cross more than one interval */
      while((NTPS > ct) && (tps[ct] <= t)){
       for (q=0;q<NSPE;q++){
         x1[q] = x[q];
         x2[q] = xx[q];
       }
       fprintf(s,"%e %e",tps[ct-1],tps[ct]); /* Print initial and final times*/
       for (q=0;q<NSPE;q++){
         fprintf(s," %e",x1[q]); /* Print first moments */
       }
       for (q=0;q<NSPE;q++){
         fprintf(s," %e",x2[q]); /* Print second moments */
       }
       fprintf(s,"\n");
	ct++;
      }
      /* Reset and add reminder */
       for (q=0;q<NSPE;q++){
         x1[q] = 0;
         x2[q] = 0;
       }
      it = t-tps[ct-1];
      v_mltadd(NSPE,x1,x,it,x1); /* First moment */
      for (q=0;q<NSPE;q++){
         xx[q] = x[q]*x[q];
      }
      v_mltadd(NSPE,x2,xx,it,x2); /* Second moment */
    }else{
      it += dt;
      v_mltadd(NSPE,x1,x,dt,x1); /* First moment */
      for (q=0;q<NSPE;q++){
         xx[q] = x[q]*x[q];
      }
      v_mltadd(NSPE,x2,xx,dt,x2); /* Second moment */
    }
    if (t<maxt){
      r1 = (rand() / (double) RAND_MAX)*rsum;
      lsum = 0;
      for (j=1;j<=NRXN;j++){
        lsum += R[j-1];
        if (lsum >= r1){
          switch(j){
          case 1: /* G0 => G1 */
            x[0] -= 1; /* G0 */
            x[1] += 1; /* G1 */
            R[0] = p[0]*x[0]; /* k01*G0 */
            R[1] = p[1]*x[1]; /* k10*G1 */
            R[2] = p[2]*x[1]; /* km*G1 */
            break;
          case 2: /* G1 => G0 */
            x[0] += 1; /* G0 */
            x[1] -= 1; /* G1 */
            R[0] = p[0]*x[0]; /* k01*G0 */
            R[1] = p[1]*x[1]; /* k10*G1 */
            R[2] = p[2]*x[1]; /* km*G1 */
            break;
          case 3: /* G1 => G1 + M */
            x[2] += 1; /* M */
            R[3] = p[3]*x[2]; /* gm*M */
            break;
          case 4: /* M =>  */
            x[2] -= 1; /* M */
            R[3] = p[3]*x[2]; /* gm*M */
            break;
          default:
            break;
          }
          break;
	}
      }
      event++;
    }
  }while((maxE>event) && (maxt>t));
  return;
}

