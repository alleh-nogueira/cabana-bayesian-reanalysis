# A Bayesian (Re-)Analysis of CABANA Trial

We provide the R code used in our independent Bayesian (re-)analysis of the [CABANA trial](10.1001/jama.2019.0693). Randomization and outcome data were extracted directly from the original publication.

Please, reproduce the results presented in the paper by running the following files:
- [*Figure 1.*](figure-01_priors-used.pdf) [Priors used for the Bayesian reanalysis of the CABANA trial.](code-01_formulate-priors.R)
- [*Figure 2.*](figure-02_posteriors-non-informative.pdf) [Bayesian reanalysis of CABANA under a noninformative prior.](code-02_non-informative-analysis.R)
- [*Figure 3.*](figure-03_posteriors-empirical.pdf) [Bayesian reanalysis of CABANA under evidence-based priors.](code-03_empirical-analysis.R)
- [*Figure 4.*](figure-04_posteriors-standard-primary.pdf) [Bayesian reanalysis of the primary outcome under standardized theoretical priors.](code-04_standard-analysis-primary.R)

The results presented in the supplement may also be reproduced by running the following files:
- [*Figure S1.*](figure-05_posteriors-standard-death.pdf) [Bayesian reanalysis of mortality under standardized theoretical priors.](code-05_standard-analysis-death.R)
- [*Table S1.*](table-01_benefit-probabilities-primary.xlsx) [Probability of treatment effects on the primary outcome according to varying prior beliefs.](code-06_benefit-probabilities-primary.R)
- [*Table S2.*](table-02_benefit-probabilities-death.xlsx) [Probability of treatment effects on the primary outcome according to varying prior beliefs.](code-07_benefit-probabilities-death.R)