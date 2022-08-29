/* mcmcMSC2s.c
MCMC program for estimating tau and theta from multiple loci of 2 sequences from 2 species 
under the multispecies coalescent mode.  Mutation model is JC.
ziheng yang on 2 February 2019

   cl -Ox mcmcMSC2s.c
   gcc -o mcmcMSC2s -O3 mcmcMSC2s.c -lm
   ./mcmcMSC2s <datafile>
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <stdarg.h>

#define NLOCI 15000
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define swap2(a,b,y) { y=a; a=b; b=y; }
#define Pi  3.1415926535897932384626433832795

int nloci = 14663;
//int nloci = 1000;
int ni[NLOCI], xi[NLOCI];

int zero(double x[], int n)
{
   int i;
   for (i = 0; i < n; i++) x[i] = 0;
   return (0);
}

double sum(double x[], int n)
{
   int i; double t = 0; 
   for (i = 0; i < n; i++) t += x[i];
   return(t);
}

unsigned int z_rndu = 1237;
double rndu(void)
{
/* 32-bit integer assumed, from Ripley (1987) p. 46 or table 2.4 line 2. */
   z_rndu = z_rndu * 69069 + 1;
   if (z_rndu == 0 || z_rndu == 4294967295)  z_rndu = 13;
   return z_rndu / 4294967295.0;
}

void SetSeed(int seed)
{
   /* Note seed is of type int with -1 meaning "please find a seed".
   */
   if (seed <= 0)
      seed = abs(2 * (int)time(NULL) + 1);
   z_rndu = (unsigned int)seed;
}

void zerror(const char* format, ...)
{
   va_list argptr;
   va_start(argptr, format);
   fprintf(stderr, "\nerror: ");
   vfprintf(stderr, format, argptr);
   va_end(argptr);
   exit(1);
}

static time_t time_start;
void starttimer(void)
{
   time_start = time(NULL);
}

char* printtime(char timestr[])
{
   /* print time elapsed since last call to starttimer()
   */
   time_t t;
   int h, m, s;

   t = time(NULL) - time_start;
   h = (int)t / 3600;
   m = (int)(t % 3600) / 60;
   s = (int)(t - (t / 60) * 60);
   if (h) sprintf(timestr, "%d:%02d:%02d", h, m, s);
   else   sprintf(timestr, "%2d:%02d", m, s);
   return(timestr);
}

/* log density of tau & theta priors and of ti: log f(tau, theta)f(ti|tau, theta) */
#if(1)
double logprior(double tau, double theta, double ti[])
{
   double a_tau = 2, b_tau = 400, a_theta = 2, b_theta = 1000;
   int i;
   double lnp = 0, p, sumt = 0;

   /* gamma prior on tau and gamma prior on theta */
   lnp = (a_tau - 1) * log(tau) - b_tau * tau + (a_theta - 1) * log(theta) - b_theta * theta;
   /* exp densities for t */
   for (i = 0; i < nloci; i++)
      sumt += ti[i];
   lnp += nloci * log(2 / theta) - 2 / theta * sumt;
   return (lnp);
}
#else
double logprior(double tau, double theta, double ti[])
{
   int i;
   double mu_tau = 0.005, mu_theta = 0.001;
   double lnp = 0, p, sumt = 0;

   lnp = -tau / mu_tau - theta / mu_theta;
   for (i = 0; i < nloci; i++) sumt += ti[i];
   lnp += nloci * log(2 / theta) - 2 / theta * sumt;
   return(lnp);
}
#endif


/* log likelihood */
double loglikelihood(double tau, double ti[])
{
   int use_jc = 0, i;   /* use_jc = 0: infinite-sites; 1: jukes-cantor */
   double lnp = 0, p, lambda;

   if (use_jc)
      for (i = 0; i < nloci; i++) {
         p = 0.75 - 0.75 * exp(-8.0 * (tau + ti[i]) / 3);
         lnp += xi[i] * log(p) + (ni[i] - xi[i]) * log(1 - p);
      }
   else
      for (i = 0; i < nloci; i++) {
         lambda = 2 * ni[i] * (tau + ti[i]);
         lnp += xi[i] * log(lambda) - lambda;
      }

   return(lnp);
}

double lnacceptanceratio_locus(int locus, double ti, double tinew, double tau, double theta)
{
   double lnR, p, pnew;

   lnR = -2 / theta*(tinew - ti);
   p = 0.75 - 0.75 * exp(-8.0 * (tau + ti) / 3);
   pnew = 0.75 - 0.75 * exp(-8.0 * (tau + tinew) / 3);
   lnR += xi[locus] * log(pnew / p) + (ni[locus] - xi[locus]) * log((1 - pnew) / (1 - p));
   return(lnR);
}



/* PjumpOptimum = 0.30; optimum step length for Bactrian moves. */

int ResetStepLengths(double Pjump[], double finetune[], int nsteps, double PjumpOptimum)
{
   /* this abjusts the MCMC proposal step lengths, using equation 9 in
      Yang, Z. & Rodríguez, C. E. 2013 Searching for efficient Markov chain Monte Carlo proposal kernels. 
      Proc. Natl .Acad. Sci. U.S.A. 110, 19307–19312.  
      PjumpOptimum = 0.3 is also from that paper.
   */
   int j, verybadstep = 0;
   double maxstep = 99;  /* max step length */
   int noisy = 3;

   if (noisy >= 3) {
      printf("\n(nsteps = %d)\nCurrent Pjump:    ", nsteps);
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", Pjump[j]);
      printf("\nCurrent finetune: ");
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", finetune[j]);
   }

   for (j = 0; j < nsteps; j++) {
      if (Pjump[j] < 0.001) {
         finetune[j] /= 100;
         verybadstep = 1;
      }
      else if (Pjump[j] > 0.999) {
         finetune[j] = min2(maxstep, finetune[j] * 100);
         verybadstep = 1;
      }
      else {
         finetune[j] *= tan(Pi / 2 * Pjump[j]) / tan(Pi / 2 * PjumpOptimum);
         finetune[j] = min2(maxstep, finetune[j]);
      }
   }

   if (noisy >= 3) {
      printf("\nNew     finetune: ");
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", finetune[j]);
      printf("\n\n");
   }
   return(verybadstep);
}

double Eff_IntegratedCorrelationTime(double x[], int n, double* mx, double* vx, double* rho1)
{
   /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
      sequence method.
      Note that this destroys x[].
   */
   double Tint = 1, rho0 = 0, rho, m = 0, s = 0;
   int  i, ir, minNr = 10, maxNr = 2000;

   /* if(n<1000) puts("chain too short for calculating Eff? "); */
   for (i = 0; i < n; i++) m += x[i];
   m /= n;
   for (i = 0; i < n; i++) x[i] -= m;
   for (i = 0; i < n; i++) s += x[i] * x[i];
   s = sqrt(s / n);
   for (i = 0; i < n; i++) x[i] /= s;

   if (mx) { *mx = m; *vx = s * s; }
   if (s / (fabs(m) + 1) < 1E-9)
      Tint = n;
   else {
      for (ir = 1; ir < min2(maxNr, n - minNr); ir++) {
         for (i = 0, rho = 0; i < n - ir; i++)
            rho += x[i] * x[i + ir];
         rho /= (n - ir);
         if (ir == 1) *rho1 = rho;
         if (ir > minNr && rho + rho0 < 0) break;
         Tint += rho * 2;
         rho0 = rho;
      }
   }
   return (1 / Tint);
}

int main(int argc, char *argv[])
{
   FILE* fin;
   char datafile[1024] = "HC.SitesDiffs.txt", timestr[64], line[1024];
   int i, j, ir, nround = 0;
   double t[NLOCI], ti, tinew;
   double lnprior, lnpriornew, lnlike, lnlikenew, lnacceptance;
   double mx[2] = { 0 }, vx[2 * 2] = {0}, d[2];
   double naccept[3] = { 0 }, Pjump[3] = { 0 }, ntotal = 0, xtotal = 0;
   double tau = 0.01, theta = 0.01, taunew, thetanew;
   double steplength[3] = { 0.0005, 0.0005, 0.02 };  /* wtau wtheta wt */
   int nr = 100000, burnin = 8000;
   int resetsteps = 1;
   double* sample[2];   /* sample[0] has tau and sample[1] has theta */
   double  Eff[2], rho1[2];
   
   printf("Usage:\n   mcmcMSC2s <datafile>\n");
   if (argc > 1) strcpy(datafile, argv[1]);
   fin = (FILE*)fopen(datafile, "r");
   if (fin == NULL)
      zerror("data file (%s) open error\n", datafile);
   /*
   printf("\ninput L (number of loci, between 1 and 14663) (e.g., type  1000 Enter)? ");
   scanf("%d", &nloci);
   printf("input N  tau_initial  theta_inital  w_tau  w_theta  w_ti: ");
   scanf("%d%lf%lf%lf%lf%lf", &nr, &tau, &theta, steplength+0, steplength+1, steplength+2);
   */

   sample[0] = (double*)malloc(nr * sizeof(double) * 2);
   sample[1] = sample[0] + nr;
   if (sample[0] == NULL)  zerror("oom");

   /* get a seed for the random number generator */
   SetSeed(-1);
   fgets(line, 1024, fin);  /* pop up the header line */
   for (i = 0; i < nloci; i++) {
      fscanf(fin, "%d%d", &ni[i], &xi[i]);
      ntotal += ni[i];
      xtotal += xi[i];
   }
   fclose(fin);
   printf("read in %4d loci\n", nloci);
   printf("average sequence divergence = %6.0f/%6.0f = %9.6f\n", xtotal, ntotal, xtotal / ntotal);

   starttimer();
   /* initialize */
   for (i = 0; i < nloci; i++)
      t[i] = theta / 2 * (0.5 + rndu());   /* initial values, 0.5-1.5 times theta0/2 */

   if (nr < 1 || tau <= 0 || theta <= 0 || steplength[0] <= 0 || steplength[1] <= 0 || steplength[2] <= 0) {
      printf("bad initial values.  Try again.");
      exit(-1);
   }
   lnprior = logprior(tau, theta, t);
   lnlike = loglikelihood(tau, t);
   printf("Window sizes (tau theta ti): %9.6f%9.6f%9.6f\n", steplength[0], steplength[1], steplength[2]);
   printf("\nRunning MCMC to collect %d samples in mcmc.txt...", nr);
   printf("\nSampling starts after burnin (%d)...\n", burnin);
   printf("initials (tau theta): %9.6f%9.6f  lnp0 = %9.3f %9.3f\n", tau, theta, lnprior, lnlike);

   printf("\n  Pjump: tau theta     t    mean_tau mean_theta   logprior    loglike\n");
   nround = 0;
   for (ir = -burnin; ir < nr; ir++) {   /* loop over nr iterations */
      nround++;
      if (ir == 0 || (resetsteps && nround >= 100 && burnin >= 200 && ir < 0 && ir % (burnin / 4) == 0)) {
         /* reset step lengths.  Do this four times. */
         if (resetsteps && burnin >= 200) {
            for (i = 0; i < 3; i++) Pjump[i] = naccept[i] / nround;
            ResetStepLengths(Pjump, steplength, 3, 0.4);
         }
         zero(Pjump, 3);
         nround = 1;
         mx[0] = mx[1] = 0;
         vx[0] = vx[1] = vx[2] = vx[3] = 0;
         naccept[0] = naccept[1] = naccept[2] = 0;
      }

      /* change tau */
      taunew = tau + (rndu() - 0.5) * steplength[0];
      if (taunew < 0) taunew = -taunew;
      lnpriornew = logprior(taunew, theta, t);
      lnlikenew = loglikelihood(taunew, t);
      lnacceptance = lnpriornew - lnprior + lnlikenew - lnlike;
      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {  /* accept */
         tau = taunew;
         lnprior = lnpriornew;
         lnlike = lnlikenew;
         naccept[0]++;
      }
      /* change theta, no need to calculate likelihood */
      thetanew = theta + (rndu() - 0.5) * steplength[1];
      if (thetanew < 0) thetanew = -thetanew;
      lnpriornew = logprior(tau, thetanew, t);
      lnacceptance = lnpriornew - lnprior;
      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) { /* accept */
         theta = thetanew;
         lnprior = lnpriornew;
         naccept[1]++;
      }

      /* Change ti for all loci.  This uses lnacceptanceratio_locus to save computation. */
      for (j = 0; j < nloci; j++) {
         ti = t[j];
         tinew = ti + (rndu() - 0.5) * steplength[2];
         if (tinew < 0) tinew = -tinew;
         lnacceptance = lnacceptanceratio_locus(j, ti, tinew, tau, theta);  /* this is diff in logposterior  */
         if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) { /* accept */
            t[j] = tinew;
            naccept[2] += 1.0 / nloci;
         }
      }
      /* re-calculate log prior and log likelihood so that they stay up-to-date */
      lnprior = logprior(tau, theta, t);
      lnlike = loglikelihood(tau, t);

      if (ir >= 0) {
         sample[0][ir] = tau; sample[1][ir] = theta;
      }

      d[0] = tau - mx[0];
      d[1] = theta - mx[1];
      vx[0] += d[0] * d[0] * (nround - 1.0) / nround;
      vx[1] += d[0] * d[1] * (nround - 1.0) / nround;
      vx[3] += d[1] * d[1] * (nround - 1.0) / nround;
      mx[0] += d[0] / nround;
      mx[1] += d[1] / nround;

      if ((ir + 1) % 1000 == 0) {
         for (i = 0; i < 3; i++) Pjump[i] = naccept[i] / nround;
         printf("\r%4.0f%% %6.3f%6.3f%6.3f ", (ir + 1.0)/nr*100, Pjump[0], Pjump[1], Pjump[2]);
         printf(" %10.6f%10.6f %11.3f %10.3f", mx[0], mx[1], lnprior, lnlike);
         if ((ir + 1) % (nr / 10) == 0) 
			printf(" %s\n", printtime(timestr));
      }
   }
   vx[0] /= nround - 1;
   vx[2] = vx[1] /= nround - 1;
   vx[3] /= nround - 1;
   if (vx[0] > 0.0 && vx[3] > 0.0) {
      vx[0] = sqrt(vx[0]);
      vx[3] = sqrt(vx[3]);
      vx[1] = vx[2] /= vx[0] * vx[3];
   }
   printf("\nSD and correlation:\n%12.6f %12.6f\n%12.6f %12.6f\n", vx[0], vx[1], vx[2], vx[3]);
   Eff[0] = Eff_IntegratedCorrelationTime(sample[0], nr, &mx[0], &vx[0], &rho1[0]);
   Eff[1] = Eff_IntegratedCorrelationTime(sample[1], nr, &mx[1], &vx[3], &rho1[1]);

   printf("Eff:  %12.6f%12.6f \n", Eff[0], Eff[1]);
   printf("rho1: %12.6f%12.6f \n", rho1[0], rho1[1]);
}
