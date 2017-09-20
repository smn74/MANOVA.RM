# MANOVA.RM 0.1.3
* improved plotting routine for RM models, now allowing for more user-specified parameters
* excluded calculation for ATS in multivariate settings (no sensible test statistic in this context)
* include factor AgeGroup in EEGwide data
* fixed some bugs in ordering of data for multivariate higher-way layouts
and setting of random seed in permutation procedure
* CIs can now also be calculated using the resampling quantiles

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
