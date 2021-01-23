/* mcmcMSC2s-HPC.c
   MCMC program for estimating tau and theta from multiple loci of 2 sequences from 2 species
   under the multispecies coalescent mode.  Mutation model is JC.
   ziheng yang on 2 February 2019

   gcc -o mcmcMSC2s-HPC -Wno-unused-result -O3 mcmcMSC2s-HPC.c -lm -lpthread
   icc -o mcmcMSC2s-HPC -O3 mcmcMSC2s-HPC.c -lm -lpthread
   cl -Ox mcmcMSC2s-HPC.c -link libpthreadVC3.lib
   ./mcmcMSC2s <datafile>
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<pthread.h>

#if(defined(__linux__) && 0)
#define _GNU_SOURCE
#include<sched.h>
#include<unistd.h>
//#include <sys/resource.h>
//#include <sys/sysinfo.h>

void pin_to_core(int t)
{
   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(t, &cpuset);
   if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
      puts("Error while pinning thread to core. ");
      exit(-1);
   }
}
#endif

#define NTHREADS 4
#define NLOCI 15000  /* This should be larger than 14663, the number of loci in total */

struct data_s {
   int nloci, ni[NLOCI], xi[NLOCI];
   double wtau, wtheta, wt, tau, theta, t[NLOCI];
   int nr, burnin;
}  data;

struct thread_data_s {
   pthread_t thread;
   pthread_mutex_t mutex;
   pthread_cond_t cond;

   volatile int work;   /* -1: end; 0: idle waiting for work; 1 2 3 etc.: work to do */
   int id, locus_start, nloci, naccept_age;
}  thread_data[NTHREADS];

int nthreads = NTHREADS;
unsigned int z_rndu[NTHREADS*4] = {1,3,5,7};   /* if(nthreads > NTHREADS) */
int debug = 0;

double rndu(int thread_id)
{
   /* From Ripley (1987) p. 46 or table 2.4 line 2.  32-bit integer assumed.  */
   z_rndu[thread_id] = z_rndu[thread_id] * 69069 + 1;
   if (z_rndu[thread_id] == 0)  z_rndu[thread_id] = 12345671;
   return ldexp((double)(z_rndu[thread_id]), -32);
}

void SetSeed(int seed, int thread_id)
{
/* Note seed is of type int with -1 meaning "please find a seed". */
   if (seed <= 0)
      seed = abs(2 * (int)time(NULL) + 1);
   z_rndu[thread_id] = (unsigned int)seed;
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

void read_data(char datafile[], double* theta, double* tau, double t[])
{
   FILE* fin;
   char line[1024];
   int i;
   double ntotal = 0, xtotal = 0;

   // data.nloci = 1000;
   data.nloci = 14663;
   data.wtau = 0.0005;  data.wtheta = 0.0005;  data.wt = 0.0123;
   data.nr = 100000; data.burnin = 8000;
   *tau = 0.01;  *theta = 0.01;
   /* data.nloci = 14663; */
   fin = (FILE*)fopen(datafile, "r");
   if (fin == NULL) puts("data file open error");
   fgets(line, 1024, fin);  /* pop up the header line */
   for (i = 0; i < data.nloci; i++) {
      fscanf(fin, "%d%d", &data.ni[i], &data.xi[i]);
      ntotal += data.ni[i];
      xtotal += data.xi[i];
   }
   fclose(fin);
   printf("read in %5d loci\n", data.nloci);
   printf("average sequence divergence = %6.0f/%6.0f = %9.6f\n", xtotal, ntotal, xtotal / ntotal);

   /* read_data */
   for (i = 0; i < data.nloci; i++)
      t[i] = - *theta / 2 * (0.5 - rndu(0));   /* initial values, 0.5-1.5 times theta0/2 */

   printf("\ninput L (number of loci, between 1 and 14663) (e.g., type  1000 Enter)?\n");
   scanf("%d", &data.nloci);
   printf("input N  tau_0 theta_0 w_tau  w_theta  w_ti: ");
   scanf("%d%lf%lf%lf%lf%lf", &data.nr, tau, theta, &data.wtau, &data.wtheta, &data.wt);
   if (data.nr < 1 || *tau <= 0 || *theta <= 0 || data.wtau <= 0 || data.wtheta <= 0 || data.wt <= 0) {
      printf("bad initial values.  Try again.");
      exit(-1);
   }
   printf("Window sizes (tau theta ti): %9.6f%9.6f%9.6f\n", data.wtau, data.wtheta, data.wt);
   printf("\nRunning MCMC to collect %d samples in mcmc.txt...", data.nr);
   printf("\nSampling starts after burnin (%d)...\n", data.burnin);
}

/* log prior */
double logprior(double tau, double theta, double ti[])
{
   int i;
   double mu_tau = 0.005, mu_theta = 0.001;
   double lnp = 0, p, sumt = 0;

   lnp = -tau / mu_tau - theta / mu_theta;
   for (i = 0; i < data.nloci; i++) sumt += ti[i];
   lnp += data.nloci * log(2 / theta) - 2 / theta * sumt;
   return(lnp);
}

/* log likelihood */
double loglikelihood(double tau, double ti[])
{
   int i;
   double lnp = 0, p, sumt = 0;

   for (i = 0; i < data.nloci; i++) {
      p = 0.75 - 0.75 * exp(-8.0 * (tau + ti[i]) / 3);
      lnp += data.xi[i] * log(p) + (data.ni[i] - data.xi[i]) * log(1 - p);
   }
   return(lnp);
}

double logacceptanceratio_locus(int locus, double ti, double tinew, double tau, double theta)
{
   double lnR, p, pnew;

   lnR = -2 / theta * (tinew - ti);
   p = 0.75 - 0.75 * exp(-8.0 * (tau + ti) / 3);
   pnew = 0.75 - 0.75 * exp(-8.0 * (tau + tinew) / 3);
   lnR += data.xi[locus] * log(pnew / p)
       + (data.ni[locus] - data.xi[locus]) * log((1 - pnew) / (1 - p));
   return(lnR);
}

/* update coalescent times (t_i) at nloci, starting from locus_start. */
int update_times(int locus_start, int nloci, int thread_id)
{
   int i, naccept = 0;
   double tinew, lnR;

   for (i = locus_start; i < locus_start + nloci; i++) {
      tinew = data.t[i] + (rndu(thread_id) - 0.5) * data.wt;
      if (tinew < 0) tinew = -tinew;
      lnR = logacceptanceratio_locus(i, data.t[i], tinew, data.tau, data.theta);
      if (lnR >= 0 || rndu(thread_id) < exp(lnR)) {
         data.t[i] = tinew;
         naccept++;
      }
   }
   return (naccept);
}

void* thread_worker(void* arg)
{
   int id = (int)arg;
   int l0 = thread_data[id].locus_start, l1 = l0 + thread_data[id].nloci;
   int threads_start = 0; 

#if (defined(__linux__))
   // pin_to_core(threads_start + id);
#endif
   pthread_mutex_lock(&thread_data[id].mutex);
   /* loop until signalled to quit */
   while (thread_data[id].work >= 0) {
      /* wait for work available */
      if (thread_data[id].work == 0)
         pthread_cond_wait(&thread_data[id].cond, &thread_data[id].mutex);
      if (thread_data[id].work > 0) {
         if (debug) printf("Thread %d working, loci %5d -%5d\n", id, l0 + 1, l1 + 1);
         thread_data[id].naccept_age = update_times(thread_data[id].locus_start, thread_data[id].nloci, id);
         if (debug) printf("\t\t\t\tThread %d finished, loci %5d -%5d\n", id, l0 + 1, l1 + 1);
         thread_data[id].work = 0;
         pthread_cond_signal(&thread_data[id].cond);
      }
   }
   pthread_mutex_unlock(&thread_data[id].mutex);
   pthread_exit(NULL);
   return(NULL);
}

void threads_init(void)
{
   int id;
   for (id = 0; id < nthreads; ++id) {
      thread_data[id].work = 0;
      pthread_mutex_init(&thread_data[id].mutex, NULL);
      pthread_cond_init(&thread_data[id].cond, NULL);
      if (pthread_create(&thread_data[id].thread, NULL, thread_worker, (void*)id))
         puts("Cannot create thread");
   }
}

void threads_wakeup(int work, void* data)
{
   int id;
   /* put threads to work */
   for (id = 0; id < nthreads; ++id) {
      pthread_mutex_lock(&thread_data[id].mutex);
      thread_data[id].work = work;
      pthread_cond_signal(&thread_data[id].cond);
      pthread_mutex_unlock(&thread_data[id].mutex);
   }
   /* wait for threads to finish work */
   for (id = 0; id < nthreads; ++id) {
      pthread_mutex_lock(&thread_data[id].mutex);
      while (thread_data[id].work > 0)
         pthread_cond_wait(&thread_data[id].cond, &thread_data[id].mutex);
      pthread_mutex_unlock(&thread_data[id].mutex);
   }
   /* TODO: threads have now finished, master thread collects results */
}

void threads_exit()
{
   for (int id = 0; id < nthreads; ++id) {
      pthread_mutex_lock(&thread_data[id].mutex);
      thread_data[id].work = -1;
      pthread_cond_signal(&thread_data[id].cond);
      pthread_mutex_unlock(&thread_data[id].mutex);
      /* wait for worker to quit */
      if (pthread_join(thread_data[id].thread, 0))
         puts("Cannot join thread");
      pthread_cond_destroy(&thread_data[id].cond);
      pthread_mutex_destroy(&thread_data[id].mutex);
   }
}

int main(int argc, char* argv[])
{
   FILE* fout = (FILE*)fopen("mcmc.txt", "w");
   char datafile[1024] = "HC.SitesDiffs.txt", timestr[64];
   int i, ir, nround = 0;
   double lnprior, lnpriornew, lnlike, lnlikenew, lnacceptance, taunew, thetanew, meantau = 0, meantheta = 0;
   double naccept[3] = { 0 };

   printf("Usage:\n   mcmcMSC2s <datafile>\n");
   if (argc > 1) strcpy(datafile, argv[1]);
   if (fout == NULL) puts("outfile open error");   
   read_data(datafile, &data.theta, &data.tau, data.t);
   starttimer();
   for (i = 0; i < nthreads; i++)  SetSeed(-1, i);
   lnprior = logprior(data.tau, data.theta, data.t);
   lnlike = loglikelihood(data.tau, data.t);
   printf("N = %5d w_tau w_theta w_t = %9.6f%9.6f%9.6f\n", data.nr, data.wtau, data.wtheta, data.wt);
   printf("initials (tau theta): %9.6f%9.6f  lnp0 = %9.3f %9.3f\n", data.tau, data.theta, lnprior, lnlike);

   if (nthreads > 1) {
      int start = 0, remaining = data.nloci % nthreads, per_thread = data.nloci / nthreads;
      for (i = 0; i < nthreads; ++i) {
         thread_data[i].locus_start = start;
         thread_data[i].nloci = per_thread + (remaining > 0 ? 1 : 0);
         if (remaining) remaining--;
         start += thread_data[i].nloci;
         printf("Thread %d : loci %5d -- %5d\n", i, thread_data[i].locus_start + 1, thread_data[i].locus_start + thread_data[i].nloci);
      }
      threads_init();
   }

   printf("\n  Pjump: tau theta     t  mean_tau mean_theta   logprior     loglike\n");
   nround = 0;
   for (ir = -data.burnin; ir < data.nr; ir++) {   /* loop over nr iterations */
      if (debug) printf("\nmcmc round %2d\n", ir);
      if (ir == 0) {
         nround = 0;
         meantau = meantheta = 0;
         naccept[0] = naccept[1] = naccept[2] = 0;
      }
      nround++;
      /* change tau, serial move */
      taunew = data.tau + (rndu(0) - 0.5) * data.wtau;
      if (taunew < 0) taunew = -taunew;
      lnpriornew = logprior(taunew, data.theta, data.t);
      lnlikenew = loglikelihood(taunew, data.t);
      lnacceptance = lnpriornew - lnprior + lnlikenew - lnlike;
      if (lnacceptance >= 0 || rndu(0) < exp(lnacceptance)) {  /* accept */
         data.tau = taunew;
         lnprior = lnpriornew;  lnlike = lnlikenew;
         naccept[0]++;
      }
      /* change theta, serial move, no need to calculate likelihood */
      thetanew = data.theta + (rndu(0) - 0.5) * data.wtheta;
      if (thetanew < 0) thetanew = -thetanew;
      lnpriornew = logprior(data.tau, thetanew, data.t);
      lnacceptance = lnpriornew - lnprior;
      if (lnacceptance >= 0 || rndu(0) < exp(lnacceptance)) { /* accept */
         data.theta = thetanew;
         lnprior = lnpriornew;
         naccept[1]++;
      }
      /* Change ti for all loci.  This can be done in parrallel.  */
      if (nthreads <= 1)
         naccept[2] += update_times(0, data.nloci, 0);
      else {
         threads_wakeup((int)1, NULL);
         for (i = 0; i < nthreads; i++)
            naccept[2] += (double)thread_data[i].naccept_age / data.nloci;
      }
      lnprior = logprior(data.tau, data.theta, data.t);
      lnlike = loglikelihood(data.tau, data.t);

      if (ir >= 0) fprintf(fout, "%.6f\t%.6f\n", data.tau, data.theta);
      meantau += data.tau;
      meantheta += data.theta;
      if ((ir + 1) % 1000 == 0) {
         printf("\r%4.0f%% %6.3f%6.3f%6.3f ", (ir+1.0)/data.nr*100, naccept[0]/nround, naccept[1]/nround, naccept[2]/nround);
         printf(" %9.6f%9.6f %11.3f %9.3f", meantau / nround, meantheta / nround, lnprior, lnlike);
         if ((ir + 1) % (data.nr / 10) == 0) printf(" %s\n", printtime(timestr));
      }
   }
   if (nthreads > 1) threads_exit();
   fclose(fout);
   printf("\nTime used: %s\n", printtime(timestr));
}
