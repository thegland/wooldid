# wooldid: Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way Fixed Effects Regressions

This program offers a suite of tools in STATA for implementing difference-in-differences style analyses with staggered treatment onset using the two-way fixed effects approach proposed in Wooldridge (2021) and the high dimensional fixed effects estimators developed by Correia (2017). Features include:

- Estimation of multiple types of overall average treatment effects
- Estimation of event study style estimates of treatment effects in particular years relative to treatment using the user's choice of a pooled or fixed reference period
- Estimation of various pre-trend tests: various types of average pre-treatment effects, testing the hypothesis that all pre-treatment relative time period specific effects are zero, testing the hypothesis that all relative time period specific effects (pre-and post-treatment) are on a single line
- Estimation of subgroup/treatment-arm specific effects, including comparisons of effects across subgroups, tests for homogeneity of treatment effects across subgroups, and subgroup specific pre-trend tests
- Estimation of the average marginal treatment effect of a continuous treatment variable (experimental) 
- Production of event study plots, histograms of treatment effects by cohorts, and plots of treatment effects across the distribution of a continuous treatment variable
- Inference using clustered standard errors, using unconditional standard errors, and using Ibragimov and Muller (2010) style cluster-robust inference
- Estimation of treatment effect elasticities and of treatment effects using poisson estimation 

&nbsp;

### Author Information and Contact

*Author:* Thomas A. Hegland, Agency for Healthcare Research and Quality

*Contact:* thomas.hegland@ahrq.hhs.gov, thomashegland.com, @thomas_hegland  

*Citation:* Hegland, Thomas A. wooldid: Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way
Fixed Effects Regressions. Statistical Software Components, 2023.

*Disclaimer:* This document reflects only the views of the author and does not necessarily represent the views of the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government. The program with which this file is associated is the work of the author and does not come with any endorsement by or warranty from the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government.  

&nbsp;

### How to Install
 ```
 ssc install wooldid
```
or
  ```
 net install wooldid, from(https://raw.githubusercontent.com/thegland/wooldid/master/wooldidinstall/) replace
  ```

&nbsp;
    
### Short Guide to Syntax

_For more detail on syntax, please refer to wooldid's included stata help file. This is only a brief overview and does not contain as much detail as the help file._

```
wooldid y i t ttre [if] [aw/pw=weights], options
```

- y is the outcome variable
- i identifies cohorts, which may contain one or more observations per cohort-time period
- t is the time variable
- ttre is a variable containing cohort-specific dates of treatment onset; missing indicates never-treated cohorts
&nbsp;

#### Selected Options by Purpose:

##### Choice of Estimand(s): 
- att - average treatment effect within the treated group
- att_i - simple average of all cohort-specific treatment effects 
- att_it - simple average of all cohort-period-specific treatment effects
- att_itime - simple average of all cohort-specific treatment effects, where cohort-specific treatment effects are simple averages of each cohort's cohort-period-specific treatment effects  
- semielasticity - compute estimates as semielasticities (% change in y in response to treatment)

##### Standard Errors and Inference: 
- cluster(_varlist_) - specify the clustering variable(s)
- robust - use heteroskedasticity robust SEs; relevant only if cluster() is not specified
- unconditionalse - compute _unconditional_ standard errors, which account for sampling variation in relevant covariates
- customvce(_custom_text_) - specify a custom VCE option
- impvalues - implement Ibragimov and Muller-style inference

##### Event Study: 
- esfixedbaseperiod - specifies use of a fixed reference period, rather than a pooled one
- esprelength(_+int_) - number of pre-treatment coefficients
- espostlength(_+int_) - number of post-treatment coefficients
- jointtests - test hypothesis that all pre-treatment coefficients are 0
- adversarial - test hypothesis that all event study coefficients fall on a single line
- makeplots - produce event study plots
- esrelativeto(_-int_) - choice of reference period if a fixed period is being used; default: period immediately prior to treatment
- oostreatedcontrols(_+int_) - determines how cohorts treated after the end of the sample are handled

##### Subgroup Effects: 
- subgroups(_variable_) - variable whose levels specify subgroups for which estimates are to be produced
- suppresssgprimaryeffects - do not produce main results for each subgroup, only compare subgroups to eachother

##### Controls and Fixed Effects: controls() fe() timetrends coarsecohortcontrols()
- controls(_fv varlist_) - specify control variables
- fe(_fv varlist_) - specify additional fixed effects
- timetrends - include cohort-specific time trends, and make adjustments to estimation process and identification checks to accommodate these trends
- coarsecohortcontrols(_fv-required varlist_) - include interactions between specified controls and indicator variables for cohort-periods where treatment effects are of interest

##### Supporting Features: 
- poisson - switch to poisson estimation
- poisexpresults - present poisson results in terms of differences in exponentiated linear predictions
- histogramcohorteffects - present histograms of cohort-specific effects
- summarystats - produce summary statistics on the outcome variable and any included continuous treatment variables
- saveplots(_prefix_) - save all graphs produced to files with a specified prefix
- replace - overwrite graphs on disk when saving graph files

##### Continuous Treatments: 
- contreat(_variable_) - specifies the continuous treatment variable
- contreatelasticitytype(_dydx_, _eydx_, _dyex_, or _eyex_) - specifies what type of continuous treatment effects are to be estimated: standard or some time of elasticity
- lattice(), contreatpoly(), and latticeignoreweights - tools for enabling the continuous treatment variable to affect the outcome variable within cohort-periods using a flexible, stepwise function (see help file for details)
- contreatcontrols(_untreatedZero OR one or more of: basic, byTreated, byTtre, byCohort, byTime_) - specifies what controls for the effect of the continuous treatment outside treated-cohort periods are included; akin to specifying the fixed effects to be included in the binary treatment case
- contreatwithin - enables estimation of marginal effect of treatment using only within-treated cohort variation in the continuous treatment variable
-  makecfxplots - produce plots showing the effect of treatment across the distribution of the continuous treatment variable
-  cfxlattice() - specifies a variable containing the bins in which treatment effects are to be produced for makecfxplots; if not specified, the variable in lattice() will be used
-  cfxplottypes(_post, pre, and/or integers specifying relative time periods_) - specifies whether the plots showing the effect of treatment across the distribution of the continuous treatment variable will show post-treatment effects, pre-treatment effects, or effects in specific years relative to treatment

##### Technical Features: 
- update - update your local installation of wooldid using the most recent version of wooldid available on github
- verbose - print updates as estimation progresses and save the results from the underlying regression estimated to a model called ___wooldidfullmodel
- safety() - if set to "off", turns off various safety features that block estimation when wooldid detects a situation where it may yield unreliable results
- Other features documented in the help file.
 

&nbsp;

### Description of Estimation Approach

wooldid estimates various types of difference-in-differences treatment effects in the vein of Wooldridge (2021), with extensions to handle continuous treatments and estimation of effects by subgroup or treatment arm. The program does not currently accommodate cases where treatment is reversible -- once treated, cohorts must remain treated.

Broadly speaking, wooldid estimates difference-in-differences treatment effects using a two-step procedure. In step 1, wooldid estimates the underlying two-way fixed effects regression, which consists of (at baseline) a regression of the outcome variable on cohort fixed effects, time fixed effects, and a set of indicator variables used to flexibly capture the effect of treatment. Each indicator variable is 1 in a specific cohort-period (i-t) cell, and 0 otherwise. The slate of indicator variables spans all cohort-periods for which treatment effects are to be estimated. In step 2, wooldid completes a set of calculations using the {help margins} command that convert the large slate of coefficients on the treatment indicator variables from the underlying regression in step 1 into treatment effects of interest. These margins calculations proceed by computing, for all treated observations, the difference between the predicted values from the underlying regression and the counterfactual predicted values that would be observed given each observation was not treated. Various averages of these differences are then presented to obtain various types of user-specified average treatment effects. To flesh this out, the underlying regression used is as follows: 

Step 1 (Underlying Regression): Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it) 

where i indexes cohorts, t indexes time periods, j indexes individual units (possibly equivalent to cohorts, possibly more numerous than cohorts), and where Z is the set of all (i,t) pairs where i is a treated cohort and t is a treated period for cohort i. In the above regression, Tz_it is a dummy variable that is 1 if, for a given z=(iz,tz) in Z, i==iz and t==tz, TreatedFlag_it is a variable that is 1 if an observation is part of an (i,t) cell in Z and 0 otherwise, fe_i is a cohort fixed effect, fe_t is a time fixed effect,
and Y_jt is the outcome variable. 

The output of the underlying Step 1 regression is then used for margins calculations in step 2. The exact margins calculation performed depends on the estimands specified by the user. The overall average treatment effect would be obtained by computing fitted values from the underlying regression among all treated observations (those where TreatedFlag_it is equal to 1), subtracting the counterfactual fitted values for these observations that set TreatedFlag_it to 0, and then taking the average of these differences (applying any sample weights). To target estimands other than the overall average treatment effect, this procedure can be modified by either introducing a set of weights (or modifying the pre-existing sample weights) when averaging the computed differences in the margins calculation or by restricting the sample so that the average only contains observations from a subgroup or time period of interest. The margins calculation can also be modified here to compute effects in terms of semi-elasticities. 

With respect to inference, there are several approaches to standard error estimation available. The default estimates standard errors in the underlying regression and then obtains standard errors on the treatment effects produced by the margins calculation using the delta method. These errors are conditional standard errors in the sense that they are conditional on the data: the difference in fitted values potentially depends on the value of various covariates, and this approach takes these values as fixed when performing the margins calculation. Given these standard errors, p-values are then obtained using a t-distribution of same degrees of freedom as the residual number of degrees of freedom in the underlying regression model. The alternative is to have margins produce unconditional standard errors, which account for sampling variation in covariates.  

This overall approach can be extended to cover many additional use cases with little or no modification. The user can switch from standard regression to poisson regression without requiring any essential modification to the regression equation or margins procedure. Additional controls and fixed effects can be added to the underlying regression without requiring adjustments to the margins calculation. Effects specific to particular subgroups may be estimated simply by performing the margins calculation within a specified subgroup. Where multiple subgroups are present, wooldid can use margins to calculate pairwise comparisons of the subgroup-specific effects and compute joint tests of the hypothesis that all subgroups have the same effect. This can be useful for tests of differences in treatment efficacy across treatment arms and can be used to replicate functionality that interaction terms in a traditional difference-in-differences would have provided. 

In a similar vein, event study style estimates that present treatment effects specific to various periods relative to treatment can be obtained by performing the margins calculation restricted to the subsample of observations corresponding with particular periods relative to treatment. If effects in periods prior to treatment are of interest, this requires modifying the underlying regression with additional treatment indicator variables for cohort-periods that correspond with the relative time period of interest in treated cohorts. The regression equation itself remains the same: only the set of (i,t) pairs appearing in Z changes for this approach. The margins calculations involved also remain similar, though may be performed on varying subsets of the data by relative period. When computing event study estimates, wooldid can supplement the standard treatment effect estimates with joint tests of the hypothesis that all pre-treatment relative time period effects are zero, as well as of the hypothesis that all relative time period effects fall on a single line (namely, the line of best fit through the estimated relative time period effects). 

In terms of implementing event study-style estimates, wooldid offers two approaches. The default pooled reference period approach is conceptually similar to estimating results like normal, but with the treatment onset date for all treated cohorts pulled earlier in time by just enough periods to estimate all requested relative time period effects. With this approach, we add no more (i,t) pairs to Z than is necessary and we allow the reference period to potentially pool multiple periods together for each cohort. Event study estimates derived from this approach always have a normalized effect of 0 in the earliest period reported (the pooled reference period). The alternative approach is to use a fixed reference period, in which case Z is modified to include all (i,t) pairs for treated cohorts, except the pair corresponding with the period immediately prior to treatment (or corresponding with some other user specified period). In this case, the regression is saturated with i-t indicator variables and all treatment effects are estimated relative to the single specified reference period. 

With greater modification, the approach above can be extended to handle continuous treatment variables and treatments that come with an associated measure of treatment intensity or dosage. When cohorts are coarse -- i.e., there are multiple units within each cohort-period cell -- and when a continuous treatment variable varies within i-t cells, estimates of the effect of the continuous treatment variable can be computed within each treated i-t cell. To ensure treatment effect heterogeneity is adequately captured, a flexible function of the continuous treatment variable can be used to measure its effect within each i-t cell if needed. Note that this approach is a bit speculative, and is not proposed in Wooldridge (2021). 

The particular modification to the underlying regression used for estimation in the continuous treatment context is as follows: 
Step 1 (Continuous): Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it * f(Contreat_jt)) + g(f(Contreat_jt)) 

where Contreat_jt is a continuous treatment variable that potentially varies arbitrarily across observations, where f() is a flexible function of Contreat_jt, and where
g() consists of interactions between various fixed effects-like indicator variables and f(Contreat_jt). Average treatment effects can be computed using the same margins calculation as before. The average marginal effect of the continuous treatment variable can be computed similar to the main margins calculation, but comparing differences in fitted values after modifying the variable Contreat_jt rather than modifying the variable TreatedFlag_it. In addition to overall and event study style estimates, wooldid can also compute average treatment effect estimates at varying points along the distribution of the continuous treatment variable. It can also compute various types of elasticity estimates in lieu of standard estimates of the effect of the continuous treatment variable. For more detail on how f(Contreat_jt) and g(Contreat_jt) are constructed, please refer to the help file's discussion of this in its section continuous treatments.

A final feature available in wooldid is its implementation of Ibragimov and Muller (2010) style cluster robust inference. This approach estimates the same underlying regression as normal, but then modifies the margins calculation significantly. In particular, for any treatment effect requested (be it an average treatment effect or an effect of a continuous treatment variable), this approach has the margins calculation compute that effect separately within each cluster. The final estimate presented is then the simple average of these cluster-specific estimates, with inference being conducted by t-testing the cluster-level estimates. 

_Note: greater detail on the program's implementation can be found in wooldid's help file, including comparisons with other estimators, greater detail on how constinuous treatments are handled, tips for working with wooldid output, etc._

&nbsp;

### Citations

Borusyak, Kirill, Xavier Jaravel, and Jann Spiess. Revisiting Event Study Designs: Robust and Efficient Estimation. Working Paper, 2023. 

Correia, Sergio. reghdfe: Stata module for linear and instrumental-variable/gmm regression absorbing multiple levels of fixed effects. Statistical Software Components s457874, 2017.
Boston College Department of Economics. https://ideas.repec.org/c/boc/bocode/s457874.html 

Ibragimov, Rustam and Muller, Ulrich K. t-Statistic Based Correlation and Heterogeneity Robust Inference. Journal of Business & Economic Statistics, 2010. 

Wooldridge, Jeffrey M. Two-Way Fixed Effects, the Two-Way Mundlak Regression, and Difference-in-differences Estimators. Working Paper, 2021. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3906345 

