<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{An Introduction to MANOVA.RM}
-->


## Introduction

This vignette documents the use of the `MANOVA.RM` package for the analysis of semi-parametric repeated measures designs and multivariate data. 
The package consists of two parts - one for repeated measurements and one for multivariate data - which will be explained in detail below. Both functions calculate the Wald-type statistic 
(WTS) and the ANOVA-type statistic (ATS) based on means. These test statistics can be used for arbitrary semi-parametric designs, even with unequal covariance matrices among groups
and small sample sizes.
Furthermore, different resampling approaches are provided in order to improve the
small sample behavior of the test statistics. 

## The `RM` function

The `RM` function calculates the above mentioned test statistics in a repeated measures design with an arbitrary number of crossed whole-plot and 
sub-plot factors.
The resampling methods provided are a permutation procedure, a parametric bootstrap approach and a
wild bootstrap using Rademacher weights. The wild bootstrap is also implemented for the ATS.


### Data Example 1 (One whole-plot and two sub-plot factors)

For illustration purposes, we consider the data set `o2cons`, which is included in `MANOVA.RM`. 


```r
library(MANOVA.RM)
data(o2cons)
```

The data set contains measurements on the oxygen consumption of leukocytes in the presence and absence of inactivated staphylococci 
at three consecutive time points. More details on the study can be found in Friedrich et al. (2017).
Due to the study design, both time and staphylococci are sub-plot factors while the treatment (Verum vs. Placebo) is a whole-plot
factor.


```r
head(o2cons)
```

```
##     O2 Staphylococci Time Group Subject
## 1 1.48             1    6     P       1
## 2 2.81             1   12     P       1
## 3 3.56             1   18     P       1
## 4 1.04             0    6     P       2
## 5 2.07             0   12     P       2
## 6 2.81             0   18     P       2
```

We will now analyze this data using the `RM` function.
The `RM` function takes as arguments:
* `formula`: A formula consisting of the outcome variable on the left hand side of a \~ operator and the factor
variables of interest on the right hand side. An interaction term must be specified. The time variable must be the last factor in the formula.
* `data`: A data.frame, list or environment containing the variables in `formula`.
* `subject`: The column name of the subjects variable in the data frame.
* `no.subf`: The number of sub-plot factors, default is 1.
* `iter`: The number of iterations for the resampling approach. Default value is 10,000.
* `alpha`: The significance level, default is 0.05.
* `resampling`: The resampling method, one of 'Perm', 'paramBS' or 'WildBS'. Default is
set to 'Perm'.
* `CPU`: The number of cores used for parallel computing. If omitted, cores are detected via \code{\link[parallel]{detectCores}}.
* `seed`: A random seed for the resampling procedure. If omitted, no reproducible seed is set.


```r
model1 <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons, 
             subject = "Subject", no.subf = 2, iter = 10000, 
             resampling = "Perm", CPU = 1, seed = 1234)
summary(model1)
```

```
## Call: 
## O2 ~ Group * Staphylococci * Time
## 
## Descriptive:
##    Group Staphylococci Time  n    Means Lower 95 % CI Upper 95 % CI
## 1      P             0    6 12 1.321667      1.150408      1.492926
## 5      P             0   12 12 2.430000      2.195672      2.664328
## 9      P             0   18 12 3.425000      3.123488      3.726512
## 3      P             1    6 12 1.618333      1.478675      1.757991
## 7      P             1   12 12 2.434167      2.164383      2.703950
## 11     P             1   18 12 3.526667      3.272821      3.780512
## 2      V             0    6 12 1.394167      1.200641      1.587692
## 6      V             0   12 12 2.570000      2.355010      2.784990
## 10     V             0   18 12 3.676667      3.374016      3.979317
## 4      V             1    6 12 1.655833      1.471327      1.840340
## 8      V             1   12 12 2.799167      2.500171      3.098162
## 12     V             1   18 12 4.029167      3.801804      4.256529
## 
## Wald-Type Statistic (WTS):
##                          Test statistic df      p-value
## Group                         11.167304  1 8.325153e-04
## Staphylococci                 20.400635  1 6.280894e-06
## Group:Staphylococci            2.554304  1 1.099942e-01
## Time                        4113.057018  2 0.000000e+00
## Group:Time                    24.105270  2 5.829176e-06
## Staphylococci:Time             4.334106  2 1.145146e-01
## Group:Staphylococci:Time       4.302876  2 1.163168e-01
## 
## ANOVA-Type Statistic (ATS):
##                          Test statistic      df1      df2      p-value
## Group                         11.167304 1.000000 316.2776 9.323112e-04
## Staphylococci                 20.400635 1.000000      Inf 6.280894e-06
## Group:Staphylococci            2.554304 1.000000      Inf 1.099942e-01
## Time                         960.208241 1.524477      Inf 0.000000e+00
## Group:Time                     5.393468 1.524477      Inf 9.237194e-03
## Staphylococci:Time             2.365958 1.982999      Inf 9.434742e-02
## Group:Staphylococci:Time       2.147250 1.982999      Inf 1.172659e-01
## 
## p-values resampling:
##                          Perm (WTS) Perm (ATS)
## Group                        0.0032         NA
## Staphylococci                0.0003         NA
## Group:Staphylococci          0.1236         NA
## Time                         0.0000         NA
## Group:Time                   0.0004         NA
## Staphylococci:Time           0.1519         NA
## Group:Staphylococci:Time     0.1507         NA
```

The output consists of four parts: `model1$Descriptive` gives an overview of the descriptive statistics: The number of observations, 
mean and confidence intervals (based on quantiles of the t-distribution) are displayed for each factor level combination.
`model1$WTS` contains the results for the Wald-type test: The test statistic, degree of freedom and p-values based on the asymptotic \(\chi^2\) distribution
are displayed. Note that the $\chi^2$ approximation is very liberal for small sample sizes.
`model1$ATS` contains the corresponding results based on the ATS. 
This test statistic tends to rather
conservative decisions in case of small sample sizes and is even asymptotically only an approximation, thus not providing an asymptotic level $\alpha$ test.
Finally, `model1$resampling` contains the p-values based on 
the chosen resampling approach. For the ATS, only the wild bootstrap procedure can be applied.
Due to the above mentioned issues for small sample sizes, the resampling procedure is recommended for such situations.

### Data Example 2 (Two sub-plot and two whole-plot factors)

We consider the data set `EEG` from the `MANOVA.RM` package: At the Department of Neurology, University Clinic of Salzburg, 160 patients were diagnosed
with either AD, MCI, or SCC, based on neuropsychological diagnostics (Bathke et al.(2015)). This data set contains z-scores for 
brain rate and Hjorth complexity,
each measured at frontal, temporal and central electrode positions and averaged across hemispheres. In addition to standardization, complexity
values were multiplied by -1 in order to make them more easily comparable to brain rate
values: For brain rate we know that the values decrease with age and pathology, while
Hjorth complexity values are known to increase with age and pathology.
The three between-subjects factors considered were sex (men vs. women), diagnosis (AD
vs. MCI vs. SCC), and age ($< 70$ vs. $>= 70$ years). Additionally, the within-subjects factors region (frontal, temporal, central) and 
feature (brain rate, complexity) structure the response vector.


```r
data(EEG)
set.seed(987)
EEG_model <- RM(resp ~ sex * diagnosis * feature * region, 
                     data = EEG, subject = "id", no.subf = 2, resampling = "WildBS",
                     iter = 1000,  alpha = 0.01, CPU = 1)
summary(EEG_model)
```

```
## Call: 
## resp ~ sex * diagnosis * feature * region
## 
## Descriptive:
##    sex diagnosis    feature   region  n       Means Lower 99 % CI
## 1    M        AD  brainrate  central 12 -1.00974303   -4.88050009
## 13   M        AD  brainrate  frontal 12 -1.00676081   -4.99056034
## 25   M        AD  brainrate temporal 12 -0.98728648   -4.49344448
## 7    M        AD complexity  central 12 -1.48789095  -10.05319599
## 19   M        AD complexity  frontal 12 -1.08562580   -6.90620385
## 31   M        AD complexity temporal 12 -1.32044015   -7.20316997
## 3    M       MCI  brainrate  central 27 -0.44742815   -1.59051263
## 15   M       MCI  brainrate  frontal 27 -0.46371997   -1.64609499
## 27   M       MCI  brainrate temporal 27 -0.50628015   -1.58436017
## 9    M       MCI complexity  central 27 -0.25680596   -1.13870013
## 21   M       MCI complexity  frontal 27 -0.45922121   -1.99708310
## 33   M       MCI complexity temporal 27 -0.49000571   -1.79618083
## 5    M       SCC  brainrate  central 20  0.45927248   -0.41381197
## 17   M       SCC  brainrate  frontal 20  0.24296492   -0.67023546
## 29   M       SCC  brainrate temporal 20  0.40875170   -1.21026961
## 11   M       SCC complexity  central 20  0.34866657   -0.06998078
## 23   M       SCC complexity  frontal 20  0.09532784   -1.03651243
## 35   M       SCC complexity temporal 20  0.31400312   -0.59757724
## 2    W        AD  brainrate  central 24 -0.29374893   -1.97838450
## 14   W        AD  brainrate  frontal 24 -0.15918975   -1.81317962
## 26   W        AD  brainrate temporal 24 -0.28539451   -1.77644839
## 8    W        AD complexity  central 24 -0.12774393   -1.37182217
## 20   W        AD complexity  frontal 24  0.02573385   -1.21242033
## 32   W        AD complexity temporal 24 -0.19371702   -1.67017970
## 4    W       MCI  brainrate  central 30 -0.10647305   -1.07581790
## 16   W       MCI  brainrate  frontal 30 -0.07356190   -1.03211025
## 28   W       MCI  brainrate temporal 30 -0.06924853   -1.06377356
## 10   W       MCI complexity  central 30  0.09398357   -0.46428563
## 22   W       MCI complexity  frontal 30  0.13130291   -0.76790741
## 34   W       MCI complexity temporal 30  0.12144023   -0.65193792
## 6    W       SCC  brainrate  central 47  0.53736580   -0.04918933
## 18   W       SCC  brainrate  frontal 47  0.54829110   -0.06244934
## 30   W       SCC  brainrate temporal 47  0.55891259   -0.01510047
## 12   W       SCC complexity  central 47  0.38428656    0.10967160
## 24   W       SCC complexity  frontal 47  0.40347289   -0.03793270
## 36   W       SCC complexity temporal 47  0.50641224    0.13237975
##    Upper 99 % CI
## 1      2.8610140
## 13     2.9770387
## 25     2.5188715
## 7      7.0774141
## 19     4.7349522
## 31     4.5622897
## 3      0.6956563
## 15     0.7186551
## 27     0.5717999
## 9      0.6250882
## 21     1.0786407
## 33     0.8161694
## 5      1.3323569
## 17     1.1561653
## 29     2.0277730
## 11     0.7673139
## 23     1.2271681
## 35     1.2255835
## 2      1.3908866
## 14     1.4948001
## 26     1.2056594
## 8      1.1163343
## 20     1.2638880
## 32     1.2827457
## 4      0.8628718
## 16     0.8849864
## 28     0.9252765
## 10     0.6522528
## 22     1.0305132
## 34     0.8948184
## 6      1.1239209
## 18     1.1590315
## 30     1.1329256
## 12     0.6589015
## 24     0.8448785
## 36     0.8804447
## 
## Wald-Type Statistic (WTS):
##                              Test statistic df      p-value
## sex                              9.97296133  1 1.588558e-03
## diagnosis                       42.38284398  2 6.261558e-10
## sex:diagnosis                    3.77699775  2 1.512988e-01
## feature                          0.08646315  1 7.687226e-01
## sex:feature                      2.16726982  1 1.409763e-01
## diagnosis:feature                5.31686865  2 7.005782e-02
## sex:diagnosis:feature            1.73538341  2 4.199197e-01
## region                           0.06961597  2 9.657908e-01
## sex:region                       0.87591221  2 6.453541e-01
## diagnosis:region                 6.12148594  4 1.902575e-01
## sex:diagnosis:region             1.53151907  4 8.210440e-01
## feature:region                   0.65247268  2 7.216346e-01
## sex:feature:region               0.42264779  2 8.095118e-01
## diagnosis:feature:region         7.14232065  4 1.285557e-01
## sex:diagnosis:feature:region     2.27379685  4 6.855437e-01
## 
## ANOVA-Type Statistic (ATS):
##                              Test statistic      df1      df2      p-value
## sex                              9.97296133 1.000000 657.4156 1.661290e-03
## diagnosis                       13.12350327 1.343095 657.4156 5.927116e-05
## sex:diagnosis                    1.90378291 1.343095 657.4156 1.635926e-01
## feature                          0.08646315 1.000000      Inf 7.687226e-01
## sex:feature                      2.16726982 1.000000      Inf 1.409763e-01
## diagnosis:feature                1.43695004 1.561987      Inf 2.380325e-01
## sex:diagnosis:feature            1.03088198 1.561987      Inf 3.416028e-01
## region                           0.01784709 1.610689      Inf 9.650363e-01
## sex:region                       0.37086704 1.610689      Inf 6.440475e-01
## diagnosis:region                 1.09114816 2.045826      Inf 3.368500e-01
## sex:diagnosis:region             0.37621401 2.045826      Inf 6.912290e-01
## feature:region                   0.12552423 1.420760      Inf 8.099026e-01
## sex:feature:region               0.07709516 1.420760      Inf 8.636365e-01
## diagnosis:feature:region         0.82935062 1.624110      Inf 4.146809e-01
## sex:diagnosis:feature:region     0.61120071 1.624110      Inf 5.098169e-01
## 
## p-values resampling:
##                              WildBS (WTS) WildBS (ATS)
## sex                                 0.000        0.000
## diagnosis                           0.000        0.000
## sex:diagnosis                       0.132        0.146
## feature                             0.797        0.797
## sex:feature                         0.141        0.141
## diagnosis:feature                   0.062        0.254
## sex:diagnosis:feature               0.462        0.385
## region                              0.972        0.982
## sex:region                          0.665        0.697
## diagnosis:region                    0.172        0.329
## sex:diagnosis:region                0.870        0.810
## feature:region                      0.787        0.921
## sex:feature:region                  0.900        0.960
## diagnosis:feature:region            0.091        0.498
## sex:diagnosis:feature:region        0.732        0.662
```

We find significant effects at level $\alpha = 0.01$ of the whole-plot factors sex and diagnosis, while none of the sub-plot factors or interactions become significant.

### Plotting

The `RM()` function is equipped with a plotting option, displaying the calculated means along with $(1-\alpha)$ confidence intervals.
The `plot` function takes an `RM` object as an argument. In addition, the factor of interest may be specified. If this argument is 
omitted in a two- or higher-way layout, the user is asked to specify the factor for plotting. Furthermore, additional graphical parameters
can be used to customize the plots. The optional argument `legendpos` specifies the position of the legend in higher-way layouts.


```r
plot(EEG_model, factor = "sex", main = "Effect of sex on EEG values")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
plot(EEG_model, factor = "sex:diagnosis", legendpos = "topleft", col = c(4, 2))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
plot(EEG_model, factor = "sex:diagnosis:feature", legendpos = "center")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)



## The `MANOVA` function

The `MANOVA` function calculates the above mentioned test statistics for multivariate data in a design with crossed or nested factors.
Additionally, a modified ANOVA-type statistic (MATS) is calculated which has the additional advantage of being applicable to designs
involving singular covariance matrices and is invariant under scale transformations of the data, see Friedrich and Pauly (2017) for details.
The resampling methods provided are a parametric bootstrap approach and a
wild bootstrap using Rademacher weights. The wild bootstrap is also implemented for the ATS, while the parametric approach works only for the WTS and MATS,
see Konietschke et al. (2015) for details. 
Note that only balanced nested designs (i.e., the same number of factor levels $b$ for each level of the factor $A$) with up to three factors are 
implemented. Designs involving both crossed and nested factors are not implemented.


### Data Example MANOVA (two crossed factors)

We again consider the data set `EEG` from the `MANOVA.RM` package, but now we ignore the sub-plot factor structure. Therefore, we are
now in a multivariate setting with 6 measurements per patient and three crossed factors sex, age and diagnosis. Due to the small number of subjects in
some groups (e.g., only 2 male patients aged $<$ 70 were diagnosed with AD) we restrict our analyses to two factors at a time.
The analysis of this example is shown below.

The `MANOVA` function takes as arguments:
* `formula`: A formula consisting of the outcome variable on the left hand side of a \~ operator and the factor
variables of interest on the right hand side. An interaction term must be specified.
* `data`: A data.frame, list or environment containing the variables in `formula`.
* `subject`: The column name of the subjects variable in the data frame.
* `iter`: The number of iterations for the resampling procedure. Default value is 10000.
* `alpha`: The significance level, default is 0.05.
* `resampling`: The resampling method, one of 'paramBS' and 'WildBS'. Default is
set to 'paramBS'.
* `CPU`: The number of cores used for parallel computing. If omitted, cores are detected automatically.
* `seed`: A random seed for the resampling procedure. If omitted, no reproducible seed is set.



```r
data(EEG)
EEG_MANOVA <- MANOVA(resp ~ sex * diagnosis, 
                     data = EEG, subject = "id", resampling = "paramBS", 
                     iter = 1000,  alpha = 0.01, CPU = 1, seed = 987)
summary(EEG_MANOVA)
```

```
## Call: 
## resp ~ sex * diagnosis
## 
## Descriptive:
##   sex diagnosis  n      Mean 1     Mean 2     Mean 3     Mean 4
## 1   M        AD 12 -0.98728648 -1.0067608 -1.0097430 -1.3204402
## 3   M       MCI 27 -0.50628015 -0.4637200 -0.4474281 -0.4900057
## 5   M       SCC 20  0.40875170  0.2429649  0.4592725  0.3140031
## 2   W        AD 24 -0.28539451 -0.1591898 -0.2937489 -0.1937170
## 4   W       MCI 30 -0.06924853 -0.0735619 -0.1064731  0.1214402
## 6   W       SCC 47  0.55891259  0.5482911  0.5373658  0.5064122
##        Mean 5      Mean 6
## 1 -1.08562580 -1.48789095
## 3 -0.45922121 -0.25680596
## 5  0.09532784  0.34866657
## 2  0.02573385 -0.12774393
## 4  0.13130291  0.09398357
## 6  0.40347289  0.38428656
## 
## Wald-Type Statistic (WTS):
##               Test statistic df      p-value
## sex                12.604176  6 4.977046e-02
## diagnosis          55.158000 12 1.695621e-07
## sex:diagnosis       9.790162 12 6.343637e-01
## 
## ANOVA-Type Statistic (ATS):
##               Test statistic      df1 df2      p-value
## sex                 7.333559 1.796816 Inf 1.047708e-03
## diagnosis           9.624691 2.323563 Inf 2.244646e-05
## sex:diagnosis       1.536449 2.323563 Inf 2.116279e-01
## 
## modified ANOVA-Type Statistic (MATS):
##               Test statistic
## sex                 45.26305
## diagnosis          194.16517
## sex:diagnosis       18.40144
## 
## p-values resampling:
##               paramBS (WTS) paramBS (ATS) paramBS (MATS)
## sex                   0.124            NA          0.003
## diagnosis             0.000            NA          0.000
## sex:diagnosis         0.748            NA          0.223
```

The output consists of several parts: First, some descriptive statistics of the data set are displayed, namely the sample size and mean for each factor level combination and each dimension (dimensions occur in the same order as in the original data set). In this example, Mean 1 to Mean 3 correspond to the brainrate (temporal, frontal, central) while Mean 4--6 correspond to complexity. Second, the results based on the WTS are displayed. For each factor, the test statistic, degree of freedom and p-value is given. The output for the ATS is similar. For the MATS, only the value of the test statistic is given, since inference is here based on the resampling procedure only. The resampling-based p-values are finally displayed for all test statistics (where applicable).

### Confidence regions

The `MANOVA` function is equipped with a function for calculating and plotting of confidence regions.
Details on the methods can be found in Friedrich and Pauly (2017). 

#### Confidence regions

Confidence regions can be calculated using the `conf.reg` function. Note that confidence regions can only be plotted in designs with 2 dimensions.
The `conf.reg` function takes as arguments:
* `object`: A `MANOVA` object calculated via `MANOVA()`.
* `nullhypo`: In designs involving more than one factor, it is necessary to specify the null hypothesis, i.e., the contrast of interest.

As an example, we consider the data set `water` from the package `HSAUR`. The data set contains measurements of mortality and drinking water hardness for 61 cities in England and Wales. Suppose we want to analyse whether these measurements differ between northern and southern towns. Since the data set is in wide format, we first need to transform it into long format before applying the `MANOVA` function.


```r
library(HSAUR)
```

```
## Loading required package: tools
```

```r
data(water)
water$subject <- 1:dim(water)[1]
library(tidyr)
data_long <- gather(water, measurement, response, mortality:hardness, factor_key=TRUE)
head(data_long)
```

```
##   location       town subject measurement response
## 1    South       Bath       1   mortality     1247
## 2    North Birkenhead       2   mortality     1668
## 3    South Birmingham       3   mortality     1466
## 4    North  Blackburn       4   mortality     1800
## 5    North  Blackpool       5   mortality     1609
## 6    North     Bolton       6   mortality     1558
```

```r
test <- MANOVA(response ~ location, data = data_long, subject = "subject", iter = 1000, resampling = "paramBS", CPU = 1, seed = 123)
summary(test)
```

```
## Call: 
## response ~ location
## 
## Descriptive:
##   location  n    Means    Means
## 1    North 35 1633.600 30.40000
## 2    South 26 1376.808 69.76923
## 
## Wald-Type Statistic (WTS):
## Test statistic             df        p-value 
##   5.158438e+01   2.000000e+00   6.289191e-12 
## 
## ANOVA-Type Statistic (ATS):
## Test statistic            df1            df2        p-value 
##   4.909672e+01   1.089611e+00            Inf   3.267386e-13 
## 
## modified ANOVA-Type Statistic (MATS):
##          [,1]
## [1,] 69.88179
## 
## p-values resampling:
##  paramBS (WTS)  paramBS (ATS) paramBS (MATS) 
##              0             NA              0
```

```r
cr <- conf.reg(test)
cr
```

```
## Center: 
##           [,1]
## [1,] 256.79231
## [2,] -39.36923
## 
## Scale:
## [1] 10.852716  2.736354
## 
## Eigenvectors:
##      [,1] [,2]
## [1,]   -1    0
## [2,]    0   -1
```

The output consists of the necessary parameters specifying the ellipsoid: the center, the eigenvectors which determine the axes of the ellipsoid as well as the scaling factors for the eigenvectors, which are calculated based on the eigenvalues, the bootstrap quantile and the total sample size. For more information on the construction of the confidence ellipses see Friedrich and Pauly (2017). 
For observations with two dimensions, the confidence ellipse can be plotted using the generic `plot` function. The usual plotting parameters can be used to customize the plots.


```r
plot(cr, col = 2, lty = 2, xlab = "Difference in mortality", ylab ="Difference in water hardness")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)


## optional GUI

The `MANOVA.RM` package is equipped with an optional graphical user interface, which is based on `RGtk2`. The GUI may be started in `R` (if `RGtk2` is installed) using the
command `GUI.RM()` and `GUI.MANOVA()` for repeated measures designs and multivariate data, respectively. 


```r
GUI.MANOVA()
```

The user can specify the data location
(either directly or via the "load data" button), the formula, the number of iterations for the resampling approach and
the significance level. Furthermore, one needs to specify the number of sub-plot factors (for the repeated measures design only),
the 'subject' variable in the 
data frame and the resampling method.
Additionally, one can specify whether or not headers are
included in the data file, and which separator (e.g., ',' for *.csv files) and character symbols are used for decimals
in the data file. 
The GUI for `RM` also provides a plotting option, which generates a new window
for specifying the factors to be plotted (in higher-way layouts) along with a few plotting
parameters.

## References

* Bathke, A. et al. (2016).
  Using EEG, SPECT, and Multivariate Resampling Methods
  to Differentiate Between Alzheimer's and other Cognitive Impairments. arXiv preprint arXiv:1606.09004.
      
* Friedrich, S., Brunner, E. and Pauly, M. (2017). Permuting longitudinal data
  in spite of the dependencies. Journal of Multivariate Analysis, 153, 255-265.
 
* Friedrich, S., and Pauly, M. (2017). MATS: Inference for potentially singular and
  heteroscedastic MANOVA. In preparation.
 
* Konietschke, F., Bathke, A. C., Harrar, S. W. and Pauly, M. (2015). 
  Parametric and nonparametric bootstrap methods for general MANOVA. Journal 
  of Multivariate Analysis, 140, 291-301.  
