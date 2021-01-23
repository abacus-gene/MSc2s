# MSc2s

This is a simple Bayesian inference project using MCMC, to estimate tau and theta under the multispecies coalescent model for two species. The data file HC.SitesDiffs.txt has two columns, listing the number of sites (n_i) and the number of differences (x_i) between the human and the chimpanzee for the 14,663 neutral autosomal loci, from table 1 of Burgess and Yang (2008). The same model was analyzed by Takahata (1986) using a method of moments and by Takahata et al. (1995) using ML under the infinite-sites model.

Under the MSC model, the sequence divergence at locus i is 2(\tau + t_i), where t_i has the exponential distribution with mean \theta/2. Given t_i or sequence divergence, the data x_i has a binomial distribution from the Jukes-Cantor mutation model. The objective of the inference is to derive the posterior of \tau and \theta. The MCMC takes a sample from the joint conditional of \tau, \theta, {t_i}.

For the first L = 1000 loci, the posterior mean and 95% CI are tau: 0.00401 (0.00372, 0.00430), theta: 0.00396 (0.00335, 0.00457), with a correlation of -0.76. Window sizes are around .001 .001 .02, for tau, theta and ti, giving acceptance rates of 0.229 0.388 0.340.

For all L = 14663 loci, the results are tau: 0.00417 (0.00409, 0.00424), theta: 0.00436 (0.00419, 0.00452), with the correlation -0.75. Window sizes are in the order of .0002 .0005 .02, for tau, theta and ti, giving acceptance rates of around 0.320 0.232 0.360.

## References
Burgess R, Yang Z. 2008. Estimation of hominoid ancestral population sizes under Bayesian coalescent models incorporating mutation rate variation and sequencing errors. Mol Biol Evol 25:1979-1994.

Jukes TH, Cantor CR. 1969. Evolution of protein molecules. Pp. 21-123 in Munro HN, ed. Mammalian Protein Metabolism. Academic Press, New York.

Takahata N. 1986. An attempt to estimate the effective size of the ancestral species common to two extant species from which homologous genes are sequenced. Genet Res 48:187-190.

Takahata N, Satta Y, Klein J. 1995. Divergence time and population size in the lineage leading to modern humans. Theor Popul Biol 48:198-221.
