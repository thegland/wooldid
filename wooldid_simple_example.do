
****** Implement a Simple Example - The Baseline Diff-in-Diff Implementation

clear
set seed 1234
set sortseed 1234
set obs 100

* Create 10 cohorts, observed for 10 periods each
gen i = mod(_n,10) + 1
bysort i: gen t = _n

* Create the treatment onset variable:
* cohorts 1, 2, and 3 are treated in period 3
* cohorts 4, 5, and 6 are treated in period 7
* cohorts 7, 8, 9, and 10 are never treated
gen ttre = 3 if i <= 3
replace ttre = 7 if i >= 4 & i <= 6

* Create a simple treatment variable with modest heterogenous effects
gen y = rnormal(0,.11)
replace y = y + ttre if t >= ttre

* Estimate ATT including cohort (i) fixed effects and time period (t) fixed effects
* and allowing for treatment effects to vary by cohort and period.
* Cluster standard errors at the cohort level
* and compute subgroup effects by treatment onset groups.
wooldid y i t ttre, att cluster(i) subgroup(ttre)

* Re-do the above, but only allowing treatment effects to vary by treatment onset group and year,
* instead of by cohort i and year. Still include i FE in the specification, however.
* Note: requires creating an ancillary variable that consists of ttre, but with the missing values
* replaced with an arbitrary value.
gen ttre_groups = ttre
replace ttre_groups = 9999 if ttre == .
wooldid y ttre_groups t ttre, att cluster(i) subgroup(ttre) fe(i)

* Return to using i as the cohort variable, but now produce event study estimates,
* an event study plot, and using the period immediately prior to treatment as the reference period.
wooldid y i t ttre, att cluster(i) subgroup(ttre) espre(3) espost(3) esfixedbaseperiod makeplots




****** Now do an example with a continuous treatment

* Create a simple dataset where observations within cohorts have a latent intensity level that affects their outcome in all periods
* and, where the teratment causes this intensity to impact the outcome variable y even more strongly.
clear
set obs 10

gen i = _n
expand 10 
bysort i: gen t = _n 

gen ttre = 5 if i <= 5
gen treated = t >= ttre

expand 10 
gen intensity = rnormal(10,.5)

gen y = rnormal(0,.01) + intensity*2 + 3*i + t*-2 + treated*intensity

* Estimate wooldid with a continuous treatment, diffing out the effect of intensity from the treatment period using an estimate of intensity
* drawing upon both the control group and the treatment group's untreated period collectively.
* The Average Marginal effect should be about equal to 1 (the coefficient on treated*intensity) and 
* the ATT should be about 10 (the mean of the intensity variable * the AME).
wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) esf espre(5) espost(5)



