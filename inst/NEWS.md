# MANOVA.RM 0.5.1
* improve RM function such that
  - no interaction term necessary
  - within-subjects factor can be explicitely specified (instead of no.subf)
* speed up computation time by removing parallel computing

# MANOVA.RM 0.4.3
* fix ordering of data for multRM
* update documentation
* fix bug in simCI without interactions
* improve error messages for missing data

# MANOVA.RM 0.4.1
* include new function multRM for multivariate repeated measures designs
* update documentation and vignette
* update output of RM (now also states the names of the within-subject factors)
* improve post-hoc comparisons

# MANOVA.RM 0.3.4
* include warning for incorrectly labelled subjects (in RM designs with more than one whole-plot factor)
* include graphics for example documentation
* update package description
* implement checks for missing values in RM

# MANOVA.RM 0.3.1
* implementation of post-hoc comparisons
* formula without interaction term possible (in MANOVA)
* warning message for singular covariance matrices (RM)
* improve runtime for permutation
* improve output
* allow for univariate calculations in MANOVA.wide

# MANOVA.RM 0.2.3
* add documentation for S3 methods
* improve documentation in vignette
* fix handling of variables coded as factors in nested designs

# MANOVA.RM 0.2.2
* warning message for singular covariance matrices
* implement output of confidence intervals for interactions (RM design, option: CI.info in plot(...))

# MANOVA.RM 0.2.1
* improved plotting routine for RM models, now allowing for more user-specified parameters
* parametric bootstrap now also implemented for the ATS in RM-Designs
* excluded calculation for ATS in multivariate settings (no sensible test statistic in this context)
* fixed some bugs in ordering of data for multivariate higher-way layouts
and setting of random seed in permutation procedure
* CIs can now also be calculated using the resampling quantiles of the WTS
* include factor AgeGroup in EEGwide data


# MANOVA 0.1.2
* included EEG data example in wide format
* fixed a bug in the calculation of mean values and covariances in the MANOVA() function (one-way layout)

# MANOVA 0.1.1
* included new function MANOVA.wide() for data provided in wide format
* fixed a bug in the calculation of mean values and covariance matrices in two- and higher-way layouts
* re-structured the analysis of nested designs

# MANOVA 0.0.5
* included calculation and plotting options for confidence ellipsoids for contrasts in multivariate factorial designs
* fixed some bugs in Wild bootstrap routine and one-way designs
* included modified ATS (MATS) and parametric and wild bootstrap thereof in MANOVA()

# MANOVA 0.0.4
* fixed small bug in the parametric bootstrap routine in the RM() function


# MANOVA 0.0.3
* fixed some bugs in the MANOVA() function, especially in the parametric bootstrap
