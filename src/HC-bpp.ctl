          seed = -1

       seqfile = HC-neutral-A.txt
      Imapfile = HC-Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

  species&tree = 2  H  C
                    1  1
                 (H, C);

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 1000 * 14663

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.002 e   # invG(a, b) for theta
      tauprior = 3 0.02      # invG(a, b) for root tau & Dirichlet(a) for other tau's

      finetune = 1: 0.2 0.002 0.00008 0.00004 0.006 0.25 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0 0 
        burnin = 1000
      sampfreq = 2
       nsample = 100000
       threads = 8 1 1
