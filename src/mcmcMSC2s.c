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

#define NLOCI 15000
//int nloci = 14663;
int nloci = 1000;
int ni[NLOCI], xi[NLOCI];

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

/* log prior */
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

/* log likelihood */
double loglikelihood(double tau, double ti[])
{
   int i;
   double lnp = 0, p, sumt = 0;

   for (i = 0; i < nloci; i++) {
      p = 0.75 - 0.75 * exp(-8.0 * (tau + ti[i]) / 3);
      lnp += xi[i] * log(p) + (ni[i] - xi[i]) * log(1 - p);
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

int main(int argc, char *argv[])
{
   FILE *fin, *fout = (FILE*)fopen("mcmc.txt", "w");
   char datafile[1024] = "HC.SitesDiffs.txt", timestr[64];
   char line[1024];
   int i, j, ir, nround = 0;
   double t[NLOCI], ti, tinew;
   double lnprior, lnpriornew, lnlike, lnlikenew, lnacceptance, meantau = 0, meantheta = 0;
   double naccept[3] = { 0 }, ntotal = 0, xtotal = 0;
   double tau = 0.01, theta = 0.01, taunew, thetanew;
   double wtau = 0.0005, wtheta = 0.0005, wt = 0.02;
   int nr = 100000, burnin = 8000;

   printf("Usage:\n   mcmcMSC2s <datafile>\n");
   if (argc > 1) strcpy(datafile, argv[1]);
   fin = (FILE*)fopen(datafile, "r");
   if (fin == NULL) {
      printf("data file (%s) open error\n", datafile);
      exit(-1);
   }
   if (fout == NULL) {
      printf("mcmc sample file (mcmc.txt) open error\n");
      exit(-1);
   }

   printf("\ninput L (number of loci, between 1 and 14663) (e.g., type  1000 Enter)? ");
   scanf("%d", &nloci);
   printf("input N  tau_initial  theta_inital  w_tau  w_theta  w_ti: ");
   scanf("%d%lf%lf%lf%lf%lf", &nr, &tau, &theta, &wtau, &wtheta, &wt);

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
      t[i] = -theta / 2 * (0.5 - rndu());   /* initial values, 0.5-1.5 times theta0/2 */

   if (nr < 1 || tau <= 0 || theta <= 0 || wtau <= 0 || wtheta <= 0 || wt <= 0) {
      printf("bad initial values.  Try again.");
      exit(-1);
   }
   lnprior = logprior(tau, theta, t);
   lnlike = loglikelihood(tau, t);
   printf("Window sizes (tau theta ti): %9.6f%9.6f%9.6f\n", wtau, wtheta, wt);
   printf("\nRunning MCMC to collect %d samples in mcmc.txt...", nr);
   printf("\nSampling starts after burnin (%d)...\n", burnin);
   printf("initials (tau theta): %9.6f%9.6f  lnp0 = %9.3f %9.3f\n", tau, theta, lnprior, lnlike);

   printf("\n  Pjump: tau theta     t    mean_tau mean_theta   logprior    loglike\n");
   nround = 0;
   for (ir = -burnin; ir < nr; ir++) {   /* loop over nr iterations */
      if (ir == 0) {
         nround = 0;
         meantau = meantheta = 0;
         naccept[0] = naccept[1] = naccept[2] = 0;
      }
      nround++;
      /* change tau */
      taunew = tau + (rndu() - 0.5)*wtau;
      if (taunew < 0) taunew = -taunew;
      lnpriornew = logprior(taunew, theta, t);
      lnlikenew = loglikelihood(taunew, t);
      lnacceptance = lnpriornew - lnprior + lnlikenew - lnlike;
      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {  /* accept */
         tau = taunew;
         lnprior = lnpriornew;  lnlike = lnlikenew;
         naccept[0]++;
      }
      /* change theta, no need to calculate likelihood */
      thetanew = theta + (rndu() - 0.5)*wtheta;
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
         tinew = ti + (rndu() - 0.5)*wt;
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

      if (ir >= 0) fprintf(fout, "%.6f\t%.6f\n", tau, theta);
      meantau += tau;
      meantheta += theta;
      if ((ir + 1) % 1000 == 0) {
         printf("\r%4.0f%% %6.3f%6.3f%6.3f ", (ir + 1.0)/nr*100, naccept[0]/nround, naccept[1]/nround, naccept[2]/nround);
         printf(" %10.6f%10.6f %11.3f %10.3f", meantau/nround, meantheta/nround, lnprior, lnlike);
         if ((ir + 1) % (nr / 10) == 0) printf(" %s\n", printtime(timestr));
      }
   }
   fclose(fout);
}
