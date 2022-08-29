/* mcmcMSC2s.c
MCMC program for estimating tau and theta from multiple loci of 2 sequences from 2 species
under the multispecies coalescent mode.  Mutation model is JC.
ziheng yang on 2 February 2019

   cl -Ox mcmcMSC2s-GammaPriorTheta.c
   gcc -o mcmcMSC2s-GammaPriorTheta -O3 mcmcMSC2s-GammaPriorTheta.c -lm
   ./mcmcMSC2s-GammaPriorTheta <datafile>
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>

#define NLOCI 15000
//int nloci = 14663, ni[NLOCI], xi[NLOCI];
int nloci = 1000, ni[NLOCI], xi[NLOCI];
int tau_move_t = 1;          /* 0: basic window; 1: mirror gibbs without t; 2: : mirror gibbs with t change */
int tau_move_tau_kernel = 1; /* 0: window no mirror; 1: mirror exp-normal; 2: : mirror exp-cauchy */
int theta_move = 2;          /* 0: sliding window; 1: gibbs with gamma; 2: gibbs with exp-cauchy */
//int tau_move_A1 = 1;    /* A1: mirror-gibbs for tau & theta; A2: mirror on tau, t, gibbs on theta */

double a_tau = 2, b_tau = 400, a_theta = 2, b_theta = 1000;  /* gamma priors on tau & theta */
int use_jc = 0; /* use_jc = 0: infinite-sites; 1: jukes-cantor */
int print_mcmc = 1;  /* 1: print mcmc sample into file; 0: don't print  */
char mcmcf[1024] = "mcmc.txt";

#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))
#define Pi 3.1415926535897932384626433832795
#define rndexp(mean) (-(mean)*log(rndu()))
#define square(a) ((a) * (a))


int zero(double x[], int n)
{
   int i;
   for (i = 0; i < n; i++)
      x[i] = 0;
   return (0);
}

double sum(double x[], int n)
{
   int i;
   double t = 0;
   for (i = 0; i < n; i++)
      t += x[i];
   return (t);
}

unsigned int z_rndu = 1237;
double rndu(void)
{
   /* 32-bit integer assumed, from Ripley (1987) p. 46 or table 2.4 line 2. */
   z_rndu = z_rndu * 69069 + 1;
   if (z_rndu == 0 || z_rndu == 4294967295)
      z_rndu = 13;
   return z_rndu / 4294967295.0;
}

double rndNormal(void)
{
   /* Standard normal variate, using the Box-Muller method (1958), improved by
      Marsaglia and Bray (1964).  The method generates a pair of N(0,1) variates,
      but only one is used.
      Johnson et al. (1994), Continuous univariate distributions, vol 1. p.153.
   */
   double u, v, s;

   for (;;)
   {
      u = 2 * rndu() - 1;
      v = 2 * rndu() - 1;
      s = u * u + v * v;
      if (s > 0 && s < 1)
         break;
   }
   s = sqrt(-2 * log(s) / s);
   return (u * s); /* (v*s) is the other N(0,1) variate, wasted. */
}

double rndgamma(double a)
{
   /* This returns a random variable from gamma(a, 1).
      Marsaglia and Tsang (2000) A Simple Method for generating gamma variables",
      ACM Transactions on Mathematical Software, 26 (3): 363-372.
      This is not entirely safe and is noted to produce zero when a is small (0.001).
    */
   double a0 = a, c, d, u, v, x;

   if (a < 1)
      a++;

   d = a - 1.0 / 3.0;
   c = (1.0 / 3.0) / sqrt(d);

   for (;;)
   {
      do
      {
         x = rndNormal();
         v = 1.0 + c * x;
      } while (v <= 0);

      v *= v * v;
      u = rndu();

      if (u < 1 - 0.0331 * x * x * x * x)
         break;
      if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
         break;
   }
   v *= d;

   if (a0 < 1)
      v *= pow(rndu(), 1 / a0);
   if (v == 0)
      printf("\a\nrndgamma returning 0.\n");
   return v;
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
   if (h)
      sprintf(timestr, "%d:%02d:%02d", h, m, s);
   else
      sprintf(timestr, "%2d:%02d", m, s);
   return (timestr);
}

/* log prior */
double logprior(double tau, double theta, double ti[])
{
   int i;
   double lnp = 0, sumt = 0;

   /* gamma prior on tau and gamma prior on theta */ 
   lnp = (a_tau - 1) * log(tau) - b_tau * tau
       + (a_theta - 1) * log(theta) - b_theta * theta;
   /* exp densities for t */
   for (i = 0; i < nloci; i++)
      sumt += ti[i];
   lnp += nloci * log(2 / theta) - 2 / theta * sumt;
   return (lnp);
}

/* log likelihood */
double loglikelihood(double tau, double ti[])
{
   int i;
   double lnp = 0, p, lambda;

   if (use_jc)
      for (i = 0; i < nloci; i++)
      {
         p = 0.75 - 0.75 * exp(-8.0 * (tau + ti[i]) / 3);
         lnp += xi[i] * log(p) + (ni[i] - xi[i]) * log(1 - p);
      }
   else
      for (i = 0; i < nloci; i++)
      {
         lambda = 2 * ni[i] * (tau + ti[i]);
         lnp += xi[i] * log(lambda) - lambda;
      }

   return (lnp);
}

double lnacceptanceratio_locus(int locus, double ti, double tinew, double tau, double theta)
{
   double lnR, p, pnew;

   lnR = -2 / theta * (tinew - ti);
   p = 0.75 - 0.75 * exp(-8.0 * (tau + ti) / 3);
   pnew = 0.75 - 0.75 * exp(-8.0 * (tau + tinew) / 3);
   lnR += xi[locus] * log(pnew / p) + (ni[locus] - xi[locus]) * log((1 - pnew) / (1 - p));
   return (lnR);
}

double logdnorm(double x, double mu, double sigma)
{
   return -0.5 / (sigma * sigma) * (x - mu) * (x - mu) - log(sqrt(2 * Pi) * sigma);
}

double logdLN(double x, double mu, double sigma)
{ /* the log density for the log-normal distribuion */
   return -0.5 * square((log(x) - mu) / sigma) - log(sqrt(2 * Pi) * sigma * x);
}


double rndCauchy(void)
{
/* Standard Cauchy variate, generated using inverse CDF */
   return tan(Pi * (rndu() - 0.5));
}

double PDFCauchy(double x, double m, double sigma)
{
   double z = (x - m) / sigma;
   return 1 / (Pi * sigma * (1 + z * z));
}


/* PjumpOptimum = 0.30; optimum step length for Bactrian moves. */
int ResetStepLengths(double Pjump[], double finetune[], int nsteps, double PjumpOptimum)
{
   /* this abjusts the MCMC proposal step lengths, using equation 9 in
      Yang, Z. & Rodr�guez, C. E. 2013 Searching for efficient Markov chain Monte Carlo proposal kernels.
      Proc. Natl .Acad. Sci. U.S.A. 110, 19307�19312.
      PjumpOptimum = 0.3 is also from that paper.
   */
   int j, verybadstep = 0;
   double maxstep = 99; /* max step length */
   int noisy = 3;

   if (noisy >= 3)
   {
      printf("\n(nsteps = %d)\nCurrent Pjump:    ", nsteps);
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", Pjump[j]);
      printf("\nCurrent finetune: ");
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", finetune[j]);
   }

   for (j = 0; j < nsteps; j++)
   {
      if (Pjump[j] < 0.001)
      {
         finetune[j] /= 100;
         verybadstep = 1;
      }
      else if (Pjump[j] > 0.999)
      {
         finetune[j] = min2(maxstep, finetune[j] * 100);
         verybadstep = 1;
      }
      else
      {
         finetune[j] *= tan(Pi / 2 * Pjump[j]) / tan(Pi / 2 * PjumpOptimum);
         finetune[j] = min2(maxstep, finetune[j]);
      }
   }

   if (noisy >= 3)
   {
      printf("\nNew     finetune: ");
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", finetune[j]);
      printf("\n\n");
   }
   return (verybadstep);
}


int HPDinterval(double x[], int n, double HPD[2], double alpha)
{
   /* This calculates the HPD interval at the alpha level.
   */
   int jL0 = (int)(n * alpha / 2), jU0 = (int)(n * (1 - alpha / 2)), jL, jU, jLb = jL0;
   double w0 = x[jU0] - x[jL0], w = w0;

   HPD[0] = x[jL0];
   HPD[1] = x[jU0];
   if (n < 3) return(-1);
   for (jL = 0, jU = jL + (jU0 - jL0); jU < n; jL++, jU++) {
      if (x[jU] - x[jL] < w) {
         jLb = jL;
         w = x[jU] - x[jL];
      }
   }
   HPD[0] = x[jLb];
   HPD[1] = x[jLb + jU0 - jL0];
   return(0);
}

double Eff_IntegratedCorrelationTime(double x[], int n, double* mx, double* vx, double* rho1)
{
   /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
      sequence method.
      Note that this destroys x[].
   */
   double Tint = 1, rho0 = 0, rho, m = 0, s = 0;
   int i, ir, minNr = 10, maxNr = 2000;

   /* if(n<1000) puts("chain too short for calculating Eff? "); */
   for (i = 0; i < n; i++)
      m += x[i];
   m /= n;
   for (i = 0; i < n; i++)
      x[i] -= m;
   for (i = 0; i < n; i++)
      s += x[i] * x[i];
   s = sqrt(s / n);
   for (i = 0; i < n; i++)
      x[i] /= s;

   if (mx)
   {
      *mx = m;
      *vx = s * s;
   }
   if (s / (fabs(m) + 1) < 1E-9)
      Tint = n;
   else
   {
      for (ir = 1; ir < min2(maxNr, n - minNr); ir++)
      {
         for (i = 0, rho = 0; i < n - ir; i++)
            rho += x[i] * x[i + ir];
         rho /= (n - ir);
         if (ir == 1)
            *rho1 = rho;
         if (ir > minNr && rho + rho0 < 0)
            break;
         Tint += rho * 2;
         rho0 = rho;
      }
   }
   return (1 / Tint);
}

int comparedouble(const void* a, const void* b)
{
   double aa = *(double*)a, bb = *(double*)b;
   return (aa > bb ? 1 : (aa < bb ? -1 : 0));
}

void summarize(double* data, int p, int n)
{
   /* this summarizes the mcmc sample, data[p*n].
   sorting destroys the original sample so a copy is made into y.  */
   int i, j, k;
   double *y, rho1, median, EqualTail[2], HPD[2], Tint;
   double *mean, *variance;

   y = (double*)malloc((n + p + p * p) * sizeof(double));
   if ((y == NULL)) zerror("oom");
   mean = y + n;
   variance = mean + p;

   /* mean, median, Eff, rho1 */
   for (j = 0; j < p; j++) {
      memmove(y, data + j * n, n * sizeof(double));
      Tint = 1 / Eff_IntegratedCorrelationTime(y, n, &mean[j], &variance[j * p + j], &rho1);
      memmove(y, data + j * n, n * sizeof(double));
      qsort(y, (size_t)n, sizeof(double), comparedouble);

      median = (n % 2 == 0 ? (y[n / 2 - 1] + y[n / 2]) / 2 : y[n / 2]);
      EqualTail[0] = y[(int)(n * .025)];
      EqualTail[1] = y[(int)(n * .975)];

      HPDinterval(y, n, HPD, 0.05);
      printf("  var %2d: mean&sd = %9.6f%9.6f ", j + 1, mean[j], sqrt(variance[j * p + j]));
      printf("equal-tail: (%8.6f, %8.6f)  HPD (%8.6f, %8.6f) ", EqualTail[0], EqualTail[1], HPD[0], HPD[1]);
      printf("eff: %6.4f rho1: %6.4f\n", 1 / Tint, rho1);
   }

   /* variance-covariance matrix */
   memset(variance, 0, p * p * sizeof(double));
   for (j = 0; j < p; j++)
      for (k = 0; k <= j; k++)
         for (i = 0; i < n; i++)
            variance[j * p + k] += (data[j * n + i] - mean[j]) * (data[k * n + i] - mean[k]);
   for (j = 0; j < p; j++)
      for (k = 0; k < j; k++)
         variance[k * p + j] = variance[j * p + k] /= sqrt(variance[j * p + j] * variance[k * p + k]);
   for (j = 0; j < p; j++)
      variance[j * p + j] = sqrt(variance[j * p + j] / (n - 1.0));
   printf("\nCorrelation matrix (SD on diagonal)\n");
   for (j = 0; j < p; j++) {
      for (k = 0; k < p; k++)
         printf(" %12.6f", variance[j * p + k]);
      printf("\n");
   }

   printf("\nDebugging info: mean & sd of log-tau and log-theta by delta technique\n");
   for (j = 0; j < p; j++) {   /* variance[j * p + j] has the SD */
      double sln = variance[j*p+j] / mean[j], mln = log(mean[j]) + sln * sln / 2;
      printf("%10s:%12.6f%12.6f\n", (j==0?"tau":"theta"), mln, sln);
   }
   free(y);
}

int main(int argc, char* argv[])
{
   FILE* fin;
   char datafile[1024] = "HC.SitesDiffs.txt", timestr[64], line[1024];
   int i, j, ir, nround = 0;
   double t[NLOCI], ti, tinew, tnew[NLOCI];
   double lnprior, lnpriornew, lnlike, lnlikenew, lnacceptance;
   double mtau=0, mtheta=0, d[2];
   double naccept[3] = { 0 }, Pjump[3] = { 0 }, ntotal = 0, xtotal = 0;
   double tau = 0.005, theta = 0.005, taunew, thetanew;
   double steplength[3] = { 0.0005, 0.0005, 0.02 }; /* wtau wtheta wt */
   int nr = 100000, burnin = 16000;
   int resetsteps = 1;
   double *sample[2]; /* sample[0] has tau and sample[1] has theta */
   FILE* fmcmc;
   double x, xnew, y, ynew, mu_lntau, mu_lntheta, sigma_lntau, sigma_lntheta, sigma_lntautheta, rho;
   double A, B, a1, b1, a = a_theta, b = b_theta;

   printf("Usage:\n   mcmcMSC2s <datafile>\n");
   if (argc > 1)
      strcpy(datafile, argv[1]);
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
   if (sample[0] == NULL)
      zerror("oom");
   if (print_mcmc) {
      fmcmc = fopen(mcmcf, "w");
      if (fmcmc == NULL) zerror("%s file open error", mcmcf);
      fprintf(fmcmc, "gen\ttau\ttheta\tlnprior\tlnlike\n");
   }

   /* get a seed for the random number generator */
   SetSeed(-123);
   fgets(line, 1024, fin); /* pop up the header line */
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
      t[i] = rndexp(theta / 2);  /* initial t[i] from the exponential with mean theta/2 */

   if (nr < 1 || tau <= 0 || theta <= 0 || steplength[0] <= 0 || steplength[1] <= 0 || steplength[2] <= 0) {
      printf("bad initial values.  Try again.");
      exit(-1);
   }
   lnprior = logprior(tau, theta, t);
   lnlike = loglikelihood(tau, t);
   printf("Window sizes (tau theta ti): %9.6f%9.6f%9.6f\n", steplength[0], steplength[1], steplength[2]);
   printf("\nRunning MCMC to collect %d samples in mcmc.txt...", nr);
   printf("\nSampling starts after burnin (%d)...\n", burnin);
   printf("initials (tau theta): %9.6f%9.6f  lnprior & lnlike: %9.3f %9.3f\n", tau, theta, lnprior, lnlike);

   printf("\n Pjump:   tau  theta      t    mean_tau mean_theta   mean_t    logprior    loglike\n");

   mu_lntau = 0, mu_lntheta = 0, sigma_lntau = 0, sigma_lntheta = 0, sigma_lntautheta = 0, rho = 0;
   for (ir = -burnin; ir < nr; ir++) { /* loop over nr iterations */
      nround++;   /*** 2022.8.21 ziheng edits ***/
      if (ir == 0 || (resetsteps && nround >= 100 && burnin >= 200 && ir < 0 && ir % (burnin / 4) == 0)) {
         /* reset step lengths.  Do this four times. */
         if (resetsteps && burnin >= 200) {
            for (i = 0; i < 3; i++)
               Pjump[i] = naccept[i] / nround;
            ResetStepLengths(Pjump, steplength, 3, 0.4);
         }

         sigma_lntau = sqrt(sigma_lntau / (nround - 1));
         sigma_lntheta = sqrt(sigma_lntheta / (nround - 1));
         sigma_lntautheta /= nround - 1;
         rho = sigma_lntautheta / (sigma_lntau * sigma_lntheta);
         printf("mean_x, mean_y, SD_x, SD_y, Corr:  %12.6f%12.6f%12.6f%12.6f%12.6f \n", mu_lntau, mu_lntheta, sigma_lntau, sigma_lntheta, rho);
         /*** 2022.8.21 ziheng edits ***/
         if (ir < 0) { mu_lntau = mu_lntheta = sigma_lntau = sigma_lntheta = sigma_lntautheta = 0; }
         zero(Pjump, 3);
         nround = 1;
         mtau = mtheta = 0;
         naccept[0] = naccept[1] = naccept[2] = 0;
      }

      /*** (A) change tau: sliding window during burnin, mirror-gibbs after ***/
      lnacceptance = 0;
      thetanew = theta;
      memmove(tnew, t, nloci * sizeof(double));
      if (ir < 0 || tau_move_t == 0) { /* change tau: use sliding window during the burnin or if tau_move_t = 0*/
         taunew = tau + (rndu() - 0.5) * steplength[0];
         if (taunew < 0) taunew = -taunew;
      }
      else { /* change tau by mirror gibbs: x = log tau, y = log theta */
         /* propose tau' */
         x = log(tau);
         if (tau_move_tau_kernel==0)         /* window no mirror, calculate step length? */
            xnew = x + (rndu() - 0.5) * steplength[0]*200;
         else if (tau_move_tau_kernel == 1)  /* exp-normal mirror */
            xnew = rndNormal() * sigma_lntau / sqrt(2) + 2 * mu_lntau - x;
         else if (tau_move_tau_kernel == 2)  /* exp-cauchy mirror  */
            xnew = rndCauchy() * sigma_lntau / sqrt(2) + 2 * mu_lntau - x;  
         taunew = exp(xnew);
         lnacceptance = xnew - x;

         /* propose theta' */
         y = log(theta);
         double s = sqrt(1 - rho * rho) * sigma_lntheta;
         double mnew = mu_lntheta + rho * sigma_lntheta * (xnew - mu_lntau) / sigma_lntau;  /* x=log tau*/
         double znew = rndNormal();
         double m    = mu_lntheta + rho * sigma_lntheta * (x - mu_lntau) / sigma_lntau;
         double z = (y - m) / s;
         ynew = znew * s + mnew;
         thetanew = exp(ynew);
         lnacceptance += (ynew - y)    /* log(thetanew/theta) from the exp transform */
                      + (znew*znew - z*z)/2;

         /* propose t_i' */
         if (tau_move_t==2) {  /* if (tau_move_t==2) change t[] */
            lnacceptance += (ynew - y) * nloci;
            for (j = 0; j < nloci; j++) {
               tnew[j] = rndexp(thetanew / 2);
               lnacceptance += -t[j] * 2 / theta + tnew[j] * 2 / thetanew;
            }
         }
      }  /* ifelse(tau_move_t) */
      lnpriornew = logprior(taunew, thetanew, tnew);
      lnlikenew = loglikelihood(taunew, tnew);
      lnacceptance += lnpriornew - lnprior + lnlikenew - lnlike;
      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
         tau = taunew;
         theta = thetanew;
         memmove(t, tnew, nloci * sizeof(double));
         lnprior = lnpriornew;
         lnlike = lnlikenew;
         naccept[0]++;
      }

      /*** (B) change theta, no need to calculate likelihood ***/
      if (theta_move == 0) { /* theta_move == 0: sliding window */
         thetanew = theta + (rndu() - 0.5) * steplength[1];
         if (thetanew < 0)
            thetanew = -thetanew;
         lnpriornew = logprior(tau, thetanew, t);
         lnacceptance = 0;
      }
      else if (theta_move) {  /* theta_move == 1: gibbs-gamma; 2: gibbs-exp-cauchy */
         A = nloci;
         B = 2 * sum(t, nloci);
         a1 = a + A;
         b1 = b * (a + A) * (a + A + 1) / (a * (a + 1) + B * b);
         if (theta_move == 1) { /* gibbs-gamma */
            thetanew = rndgamma(a1) / b1;
            lnacceptance = (a1 - 1) * log(theta / thetanew) - b1 * (theta - thetanew);
         }
         else {                 /* gibbs-exp-cauchy */
            /* m and s for log theta approximated using the delta technique */
            double m = a1/b1, s = sqrt(a1)/b1, my = log(m) + 0.5*s*s/(m*m), sy = s/m;
            double y = log(theta), z = (y - my) / sy, ynew, znew;
            znew = rndCauchy();
            ynew = znew * sy + my;
            thetanew = exp(ynew);
            lnacceptance = log((1 + znew * znew) / (1 + z * z)) + ynew - y;
         }
      }
      lnpriornew = logprior(tau, thetanew, t);
      lnacceptance += lnpriornew - lnprior;
      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) { /* accept */
         theta = thetanew;
         lnprior = lnpriornew;
         naccept[1]++;
      }

      /*** (C) Change ti for all loci by sliding window.  This uses lnacceptanceratio_locus. ***/
      for (j = 0; j < nloci; j++) {
         ti = t[j];
         tinew = ti + (rndu() - 0.5) * steplength[2];
         if (tinew < 0)
            tinew = -tinew;
         lnacceptance = lnacceptanceratio_locus(j, ti, tinew, tau, theta); /* this is diff in logposterior  */
         if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) { /* accept */
            t[j] = tinew;
            naccept[2] += 1.0 / nloci;
         }
      }
      /* re-calculate log prior and log likelihood so that they stay up-to-date */
      lnprior = logprior(tau, theta, t);
      lnlike = loglikelihood(tau, t);

      if (ir >= 0) {
         sample[0][ir] = tau;
         sample[1][ir] = theta;
         if (print_mcmc) {
            fprintf(fmcmc, "%d\t%.6f\t%.6f\t%.2f\t%.2f\n", ir+1, tau, theta, lnprior, lnlike);
            if ((ir + 1) % 100 == 0) fflush(fmcmc);
         }
      }
      if (ir < 0) {
         d[0] = log(tau) - mu_lntau;   d[1] = log(theta) - mu_lntheta;
         sigma_lntau += d[0] * d[0] * (nround - 1.0) / nround;
         sigma_lntautheta += d[0] * d[1] * (nround - 1.0) / nround;
         sigma_lntheta += d[1] * d[1] * (nround - 1.0) / nround;
         mu_lntau += d[0] / nround;    mu_lntheta += d[1] / nround;
      }
      mtau += (tau - mtau) / nround;
      mtheta += (theta - mtheta) / nround;

      if ((ir + 1) % 1000 == 0) {
         for (i = 0; i < 3; i++)
            Pjump[i] = naccept[i] / nround;
         printf("\r%4.0f%% %7.4f%7.4f%7.4f ", (ir + 1.0) / nr * 100, Pjump[0], Pjump[1], Pjump[2]);
         printf(" %10.6f%10.6f%10.6f %11.3f %10.3f", mtau, mtheta, sum(t,nloci)/nloci, lnprior, lnlike);
         if ((ir + 1) % (nr / 10) == 0) printf(" %s\n", printtime(timestr));
      }
   }  /** for(ir) **/

   printf("\n\nSummarizing the mcmc sample\n");
   summarize(&sample[0][0], 2, nround);
   free(sample[0]);
   fclose(fmcmc);
   printf("\nlikelihood model: %s\n", (use_jc ? "jukes-cantor" : "infinite-sites"));
   printf("tau move:         %s\n", (tau_move_t == 0 ? "basic window" : (tau_move_t == 1 ? "mirror gibbs no t" : "mirror gibbs with t")));
   printf("tau move tau kernel: %s\n", (tau_move_tau_kernel == 0 ? "window" : (tau_move_tau_kernel == 1 ? "exp-normal" : "exp-cauchy")));
   printf("theta move:       %s\n", (theta_move == 0 ? "sliding-window" : (theta_move == 1 ? "gibbs gamma" : "gibbs exp-cauchy")));
}
