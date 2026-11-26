{smcl}
{* *! version 1.3 2024-07-31}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ftools" "help ftools"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{vieweralsosee "gtools" "help gtools"}{...}
{vieweralsosee "margins" "help margins"}{...}
{vieweralsosee "hdidregress" "help hdidregress"}{...}
{vieweralsosee "did_imputation" "help did_imputation"}{...}
{viewerjumpto "Syntax" "wooldid##syntax"}{...}
{viewerjumpto "Description" "wooldid##description"}{...}
{viewerjumpto "All Options" "wooldid##options"}{...}
{viewerjumpto "Specifying the Estimand(s)" "wooldid##specifyestimand"}{...}
{viewerjumpto "Standard Errors" "wooldid##standarderrors"}{...}
{viewerjumpto "Bootstrap" "wooldid##bootstrap"}{...}
{viewerjumpto "Randomization Inference" "wooldid##randinf"}{...}
{viewerjumpto "Supporting Features" "wooldid##supporting"}{...}
{viewerjumpto "Event Studies" "wooldid##eventstudy"}{...}
{viewerjumpto "Subgroup Effects" "wooldid##interactions"}{...}
{viewerjumpto "Continuous Treatments" "wooldid##continuous"}{...}
{viewerjumpto "Controls and Fixed Effects" "wooldid##controls"}{...}
{viewerjumpto "Technical Features" "wooldid##technical"}{...}
{viewerjumpto "Working with Output" "wooldid##output"}{...}
{viewerjumpto "Examples" "wooldid##examples"}{...}
{viewerjumpto "Bugs and Assistance" "wooldid##bugs"}{...}
{viewerjumpto "Citations" "wooldid##citations"}{...}
{viewerjumpto "Author Info and How to Cite" "wooldid##author"}{...}
{viewerjumpto "Disclaimer" "wooldid##disclaimer"}{...}
{title:Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment }{...}
{title:Onset Using Heterogeneity-Robust Two-Way Fixed Effects Regressions}

{pstd} {bf:wooldid} -- Suite for implementing difference-in-differences style analyses with
staggered treatment onset using the two-way fixed effects approach proposed in Wooldridge (2021)
and the high dimensional fixed effects estimators developed by Correia (2017). Features include
estimation of various types of average treatment effects, comparison of treatment effects
across subgroups, estimation of event study style estimates, parallel trends assumption tests
following Roth (2022), and estimation of average and marginal treatment effects in the presence of a
continuous treatment (or treatment intensity) variable. This suite also can produce event study
plots, histograms of treatment effects by cohort, and plots of treatment
effects across the distribution of a specified continuous treatment variable. Options for inference
include analytic SEs, bootstrap (standard and Bayesian), randomization inference,
simultaneous bands, and Ibragimov and Muller (2010) style cluster-robust inference.  {p_end}

{pstd}{ul:Shortcuts:} {help wooldid##description:Description}, {help wooldid##options:Options},
{help wooldid##output:Working with Output}, {help wooldid##examples:Examples},
{help wooldid##bugs:Bugs and Assistance}, {help wooldid##citations:Citations},
{help wooldid##author:Author Info and How to Cite}, {help wooldid##disclaimer:Disclaimer} {p_end}

{pstd}{ul:Dependencies:} {cmd:wooldid} requires users have STATA 16 or higher installed,
{help reghdfe} version 5.7.3, 
{help ftools}, and {help gtools}. Poisson estimation also requires the package
{help ppmlhdfe}. {p_end}

{pstd}{ul:Github Link:} {browse "https://github.com/thegland/wooldid/"} - the command
{it:wooldid, update} will cause this program to update itself using the materials
at this link. {p_end}



{marker syntax}{...}
{title:Syntax}

{phang}{cmd: wooldid} {it:y i t ttre} [if] [aw/pw=weights] [{cmd:,}
{help wooldid##specifyestimand:estimands} {help wooldid##standarderrors:standard_errors}
{help wooldid##bootstrap:bootstrap} {help wooldid##randinf:randomization_inference}
{help wooldid##supporting:supporting_features} {help wooldid##eventstudy:event_study}
{help wooldid##interactions:subgroup_effects} {help wooldid##continuous:continuous_treatment}
{help wooldid##controls:controls_and_FE} {help wooldid##technical:technical_features}]{p_end}

{synoptset 22 tabbed}{...}
{synopt : {it:y}}Outcome variable. {p_end}
{synopt : {it:i}}Cohort identifier; represents collections of units allowed to have different
treatment effects by year. Cohorts must either group units 
by year of treatment onset, or must be subsets of such treatment onset groups. Coarser cohort
definitions allow for faster estimation and more stable 
estimation of SEs with clusters (with analytic clustered SEs, cohorts should be coarser than
clusters). Setting cohorts equal to the unit-level identifier, 
in general, yields poor results (unit FE can be included separately). The cohort variable i must
consist only of integers greater than 0. {p_end}
{synopt : {it:t}}Time period identifier; must consist of integers greater than 0.{p_end}
{synopt : {it:ttre}}Variable containing cohort-specific dates of treatment onset. Missing indicates
never-treated cohorts; cohorts treated before the sample begins are dropped. {p_end}

{phang} {ul:Options by Purpose:} {p_end}
{synoptset 1 tabbed }{...}
{synopt : {help wooldid##specifyestimand:Specifying the Estimand(s):}} {opt att} {opt att_it}
{opt att_i} {opt att_itime} {p_end}

{synopt : {help wooldid##standarderrors:Standard Errors:}}
{opt clus:ter(variable(s)[numerical])}  {opt imp:values} 
{opt simult:aneousbands} {opt level(#)}  {p_end}

{synopt : {help wooldid##bootstrap:Bootstrap Inference:}} {opt boot:strap(standard|bayes)}
{opt bootci:type(percentile|bc|bca|wald)} {opt bootdf(normal|clusters|#)} 
{opt reps(#)} {opt bootsave(name[, replace])}    {p_end}

{synopt : {help wooldid##randinf:Randomization Inference:}} {opt randinf(shuffle|probabilities)}
{opt ristr:ata(stratavariable)} {opt riprobs(spec)} 
{opt riside(left|right)} {opt riabs} {opt rit:stat} {opt reps(#)} {opt riex:act}
{opt risave(name[, replace])}  {p_end}

{synopt : {help wooldid##supporting:Supporting Features:}} {opt hist:ogramcohorteffects}
{opt save:plots(filename_stem[, replace])}  {opt sum:marystats}
{opt pois:son}    {p_end}

{synopt : {help wooldid##eventstudy:Event Study Estimation:}}  {opt esf:ixedbaseperiod}
{opt esr:elativeto(-int)} {opt espre:length(int>=0)} {opt espost:length(int>=0)}
 {opt pretrend:test(powerlevels)}  {opt makep:lots} {opt oos:treatedcontrols(int>=0)}   {p_end}

{synopt : {help wooldid##interactions:Subgroup Effects:}}
{opt sub:groups(variable[>=0 int])}  {p_end}

{synopt : {help wooldid##continuous:Continuous Treatments:}} {opt contreat(variable[numerical])}
{opt contreatelas:ticity}  {opt contreatw:ithin}
{opt contreatp:oly(+int)} {opt lattice(variable[+int] OR latticetype latticenum)}
{opt latticeig:noreweights} {opt makec:fxplots} {opt cfxlattice(variable[+int])}
{opt cfxplott:ypes(post, pre, and/or integers specifying relative time periods)}
{opt makecfxplotsbysg}  {opt cfxplotnose}  {opt contreatc:ontrols(controltype)}  {p_end}
{phang3} {it:-latticetype}: {it:all}, {it:byTime}, {it:byTreated}, {it:byTtre}, {it:byCohort},
{it:byTreatedByTime}, {it:byTtreByTime}, {it:byCohortByTime}, {it:baseTreated}, {it:baseTreatedRef},
{it:baseTreatedPost}, or {it:baseTreatedPre}  {p_end}
{phang3} {it:-latticenum}: integer > 1  {p_end}
{phang3} {it:-contreatcontroltype}: {it:untreatedZero}, or one or more of {it:basic},
{it:byTreated}, {it:byTtre}, {it:byCohort}, {it:byTime}  {p_end}

{synopt : {help wooldid##controls:Controls and Fixed Effects:}} {opt control:s(fv varlist)}
{opt fe(fv varlist)} {opt timetrend:s} {opt coarse:cohortcontrols(fv-required varlist)}
{opt inter:activecontrols(varlist)} {p_end}

{synopt : {help wooldid##technical:Technical Features:}}
{opt update} {opt fast:nose} {opt regtol(+int)} {opt ccno:absorb}
{opt ccc_absorb(varlist[numerical])}
{opt nofastcollapse} {opt rseed(#)}  {opt ver:bose} {opt iteffects} {opt safety(off)}    {p_end}



{marker description}{...}
{title:Description of Estimation Approach}

{pstd}
{bf:wooldid} estimates various types of difference-in-differences treatment effects in the vein of
Wooldridge (2021), with extensions to handle continuous treatments and estimation of effects by
subgroup or treatment arm. The program does not currently accommodate cases where treatment is
reversible -- once treated, cohorts must remain treated.{p_end}

{pstd}
Broadly speaking, {it:wooldid} estimates difference-in-differences treatment effects using a
two-step procedure. In step 1, {it:wooldid} estimates the underlying two-way fixed effects
regression,
which consists of (at baseline) a regression of the outcome variable on cohort fixed effects, time
fixed effects, and a set of indicator variables used to flexibly capture the effect of treatment.
Each indicator variable is 1 in a specific cohort-period (i-t) cell, and 0 otherwise. The slate of
indicator variables spans all cohort-periods for which treatment effects are to be estimated. In
step 2, {it:wooldid} completes a set of calculations using the {help margins} command that convert
the large slate of coefficients on the treatment indicator variables from the underlying
regression in step 1 into treatment effects of interest. These margins calculations proceed by
computing, for all treated observations, the difference between the predicted values from
the underlying regression and the counterfactual predicted values that would be observed given each
observation was not treated. Various averages of these differences are then presented to obtain
various types of user-specified average treatment effects. To flesh this out, the underlying
regression used is as follows: {p_end}

{p2col 5 8 8 0 : Step 1 (Underlying Regression):} Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it *
TreatedFlag_it * Tz_it) {p_end}

{pstd} where i indexes cohorts, t indexes time periods, j indexes individual units (possibly
equivalent to cohorts, possibly more numerous than cohorts), and where Z is the set of all (i,t)
pairs where i is a treated cohort and t is a treated period for cohort i. In the above regression,
Tz_it is a dummy variable that is 1 if, for a given z=(iz,tz) in Z, i==iz and t==tz,
TreatedFlag_it is a variable that is 1 if an observation is part of an (i,t) cell in Z and 0
otherwise, fe_i is a cohort fixed effect, fe_t is a time fixed effect,
and Y_jt is the outcome variable. {p_end}

{pstd} The results from the above step 1 regression are then used in step 2 to produce the final
average treatment effects. The exact calculations used depend on the estimand 
specified by the user, but for an average treatment effect on the treated group (ATT), the margins
calculation proceeds by computing fitted values for each treated observation 
(those where TreatedFlag_it is equal to 1), subtracting the counterfactual fitted values for these
observations when they are untreated (obtained by setting TreatedFlag_it to 0), 
and then taking the average of these differences (applying any sample weights). To target other
estimands, such as simple averages of the i-t level treatment effects, averages of 
cohort-specific treatment effects, or relative-time-period-specific treatment effects, this
procedure is modified by either computing the average only over the relevant subsample 
(e.g., for period-specific effects) or by reweighting the sample (e.g., for averages of cohort
effects). For more on the different estimand-types available, 
see {help wooldid##specifyestimand:estimands}. {p_end}

{pstd} With respect to inference, there are several approaches available. The default estimates
heteroskedasticity robust standard errors in the underlying regression and then obtains
 standard errors on the treatment effects produced by the margins calculation using the delta
method. Given these standard errors, p-values are then obtained using a t-distribution of 
the same degrees of freedom as the residual number of degrees of freedom in the underlying
regression model. Inference is also available using analytic clustered standard errors, the 
bootstrap, the Bayesian bootstrap, and randomization inference. When treatment effects can be
estimated independently within clusters, inference via the method of Ibragimov and Muller 
(2010) is also available, in which case inference is done by t-testing the set of cluster-specific
estimates against 0 (possibly changing the requested estimand). Finally, note that to 
speed estimation of analytic SEs, the delta method calculations are done at the cohort-period (i-t)
level, rather than at the observation-level - except in cases where the two methods 
might be reasonably expected to diverge. This speed-up trick can be disabled using the option
{opt nofastcollapse} if one's specification is likely to produce differences between i-t 
level averages and observation-level averages.
For more on standard error computation, see {help wooldid##standarderrors:standard_errors}. {p_end}


{pstd} This overall approach can be extended to cover many additional use cases with little or no
modification. The user can switch from standard regression to poisson regression,
and can add controls and fixed effects to the underlying regression, without any essential
modification to the margins calculation. Effects specific to particular subgroups may be estimated
simply by
performing the margins calculation within a specified subgroup, and {it:wooldid} can compare the
subgroup-specific effects pairwise and jointly. This can be useful for tests of differences
in treatment efficacy across treatment arms and can be used to replicate functionality that
interaction terms in a traditional difference-in-differences would have provided.
See {help wooldid##interactions:subgroups} for more. {p_end}

{pstd} In a similar vein, event study style estimates are obtained by performing the margins
calculation restricted to the subsample of observations corresponding with particular periods
relative to treatment. Effects in periods prior to treatment require adding to Z the (i,t) pairs of
treated cohorts that correspond with the relative time period of interest; the regression
equation itself remains the same. {it:wooldid} offers two approaches. The default pooled reference
period approach is conceptually similar to estimating results like normal, but with the
treatment onset date for all treated cohorts pulled earlier in time by just enough periods to
estimate all requested relative time period effects: no more (i,t) pairs are added to Z than
necessary, the reference period may pool multiple periods for each cohort, and the earliest period
reported (the pooled reference period) always has a normalized effect of 0. With a fixed
reference period, Z includes all (i,t) pairs for treated cohorts except the pair corresponding with
the period immediately prior to treatment (or with some other user specified period), so the
regression is saturated with i-t indicator variables and all treatment effects are estimated
relative to that single reference period. See {help wooldid##eventstudy:event_study}. {p_end}

{pstd} With greater modification, the approach extends to continuous treatment variables and to
treatments that come with an associated measure of treatment intensity or dosage. When
cohorts are coarse -- i.e., there are multiple units within each cohort-period cell -- and the
continuous treatment variable varies within i-t cells, its effect can be estimated within each
treated i-t cell, using a flexible function of it if needed: {p_end}

{p2col 5 8 8 0 : Step 1 (Continuous):} Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it *
TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it * f(Contreat_jt)) + g(f(Contreat_jt)) {p_end}

{pstd} where Contreat_jt is a continuous treatment variable that potentially varies arbitrarily
across observations, f() is a flexible function of Contreat_jt, and g() consists of
interactions between various fixed effects-like indicator variables and f(Contreat_jt). Average
treatment effects are computed using the same margins calculation as before. The average
marginal effect of the continuous treatment variable is computed similarly, but comparing
differences in fitted values after modifying Contreat_jt rather than TreatedFlag_it. {it:wooldid}
can also compute average treatment effect estimates at varying points along the distribution of the
continuous treatment variable. See {help wooldid##continuous:continuous_treatment}. {p_end}


{pstd}
{ul:Comment on Sample Formation, Cohorts Treated Prior to Sample Start:} The estimation sample
{it:wooldid} uses may be more restricted than your overall sample: {it:wooldid}
will drop observations with missing data, observations with 0-valued weights, observations from
treated cohorts not observed prior to treatment (or not observed in the user-specified
reference year if {it:esfixedbaseperiod} is specified), and observations from time periods where no
untreated
(including not yet treated) units are observed. Under continuous treatment designs, treated i-t
cells with no variation in the continuous treatment variable
are also dropped. With {it:timetrends}, so are cohorts not observed enough times to estimate a time
trend using the pre-treatment data. {p_end}

{pstd}
{ul:Comment on Cohorts Treated After Sample End:} By default, {it:wooldid} handles cohorts treated
after the end of the sample identically to cohorts that have never been treated;
see {it:oostreatedcontrols()}. {p_end}



{marker options}{...}
{title:Options}

{marker specifyestimand}{...}
{dlgtab:Specifying the Estimand(s)}

{phang}{opt att}: Estimate average treatment on the treated (ATT) effects; this option is the
default and may be combined with options att_it, att_i, and att_itime to estimate multiple different
effect types. ATT effects are computed by performing the margins calculation (see step 2 in
{help wooldid##description:Description} for more on the margins calculation)
using whatever weights (if any) are used in the main regression command.{p_end}

{phang}{opt att_it}: Estimate simple averages of the i-t cell-specific treatment effects
(ATT_it). The treatment effect in each cohort-period is calculated via the margins calculation
using any weights that may apply to the observations falling in each i-t cell. The final ATT_it
estimate is the simple, unweighted mean of all of the cohort-period-specific estimates.
This is done by reweighting the data used in the margins calculation (after estimation of the
underlying regression) by w/W_it where w is an observation's sample weight and W_it is the
total sample weight associated with all observations appearing in a given i-t cell. If no weights
are specified, w is taken to have a value of 1 for all observations.   {p_end}

{phang}{opt att_i}: Estimate the simple average of cohort-specific treatment effects (ATT_i). The
treatment effect for each cohort is calculated using the margins calculation and any sample
weights. The final ATT_i estimate is the simple, unweighted mean of the cohort-specific treatment
effects. ATT_i effects are computed by performing the margins calculation, but reweighting
observations by either w/W_ipre or w/W_ipost, where w is an observation's sample weight, W_ipre is
the total sample weight associated with cohort-i in the pre-treatment period
(but outside the reference period), and W_ipost is the total sample weight associated with cohort-i
in the post-treatment period. Each observation is then reweighted based on whether it
falls in the pre- or post-treatment period. Equivalent to {it:att_it} when producing event study
horizon-specific estimates.  {p_end}

{phang}{opt att_itime}: Estimate the simple average of cohort-specific simple average treatment
effects (ATT_itime). Treatment effect calculation here begins like for {it:att_it}, obtaining
treatment effects for individual cohort-periods of interest using weights. Then, cohort-specific
effects are obtained by taking the simple mean of the cohort-period-specific effects.
The overall presented ATT_itime estimate is then equal to the simple, unweighted mean of these
cohort-specific effects. ATT_itime effects are calculated by performing the margins calculation,
but reweighting observations by w/(W_it * T_iprepost), where w is an observation's sample weight,
W_it is the amount of sample weight falling in a given i-t cell, and T_iprepost contains
(depending on the time period) either the number of post-treatment time periods in which a cohort-i
is observed or the number of pre-treatment, non-reference period time periods in which a
cohort i is observed. Equivalent to {it:att_it} when producing event study horizon-specific
estimates. {p_end}


{marker standarderrors}{...}
{dlgtab:Standard Errors}

{phang}{opt clus:ter(clusvars)}: Variable(s) specifying clustering level. Clustered standard errors
are computed analytically.  {p_end}

{phang}{opt imp:values}: Implements Ibragimov and Muller (2010) style cluster-robust inference for
all estimates. The idea behind the Ibragimov and Muller approach is to compute
cluster-specific treatment effect estimates, present their average as one's final treatment effect
estimate, and conduct inference by t-testing the cluster-specific estimates against 0.
IM-style inference implicitly changes {it:wooldid}'s
target estimand. For example, if clustering at the cohort level, IM-style inference will cause
{it:att} to yield the same effect estimates as {it:att_i}. Note that {it:wooldid} does not adjust
the
weights used for computing different estimands when IM-style inference is requested. IM-style
inference requires a single cluster variable, with each treated cohort inside one cluster. IM-style
inference is not currently allowed with the use of subgroups or the option
{it:makecfxplots}. {p_end}

{phang}{opt simult:aneousbands}: Compute a sup-t simultaneous confidence band for the event study
relative-year effect estimates. Obtained via bootstrap, if used. Otherwise, it 
is simulated from the analytic event study covariance. {p_end}

{phang}{opt level(#)}: Confidence level for reported confidence intervals.   {p_end}



{pstd} {ul:Discussion of Ibragimov and Muller-Style Inference:}
The virtue of the Ibragimov and Muller approach is that it can significantly improve on analytic
clustered standard errors when the number of clusters is small and when the clusters are
very heterogeneous in size. A key cost to adopting this approach is the required assumption that
cluster-specific treatment effect estimates are independent. This assumption seems most likely to
be satisfied when estimating the effect of continuous treatment variables using only within-cluster
variation, or in a bespoke diff-in-diff set up where clusters contain both treated and control
units and where time fixed effects are not shared across clusters. In the standard diff-in-diff
setting, {it:wooldid} will still present IM-style estimates by computing cluster-specific treatment
effect estimates, despite the fact that the existence of a shared control group (and perhaps a
shared reference period, if {it:esfixedbaseperiod} is not being used) likely causes mischief with
the
independence assumption.  {p_end}



{marker bootstrap}{...}
{dlgtab:Bootstrap Inference}

{phang}{opt boot:strap(standard|bayes)}: Bootstrap standard errors and inference using either the
standard (pairs cluster, when clustering) bootstrap or the Bayesian bootstrap. 
The default behavior is to show standard errors calculated as the standard deviation of the
bootstrap replicate result distribution, coupled with percentile-based intervals and p-values.
The Bayesian bootstrap proceeds like the standard bootstrap, except that instead of resampling
replicates, each replicate involves re-estimation after multiplying each observation's 
(or cluster of observations') sample weight with a Dirichlet weight. Note: the Bayesian bootstrap is
not available with {it:att_it}, {it:att_i} or {it:att_itime}. {p_end}

{phang}{opt bootci:type(percentile|bc|bca|wald)}: Select the type of confidence interval and p-value
shown. By default, the CI and p-value are based on percentiles of the bootstrap replicate 
distribution ({it:percentile}). Specify wald to obtain CIs and p-values based on the bootstrapped
standard error (the standard deviation of the bootstrap distribution) coupled with either a normal
or 
t-distribution approximation, depending on {it:bootdf()}. Bias-corrected percentile intervals and
p-values are also available ({it:bc}), along with bias-corrected and accelerated ({it:bca})
intervals (only 
when clusters are specified). {p_end}

{phang}{opt bootdf(normal|clusters|#)}: Specifies the reference distribution for bootstrap inference
using just the bootstrap standard error ({it:bootcitype(wald)}). The default is {it:normal} without
clustering, and 
{it:clusters} with clustering (which specifies a t-distribution with number of clusters - 1 degrees
of freedom). An integer can also be specified to use a t-distribution with the specified number of
degrees of freedom. 
The dof set here also is used as an input into other downstream tests when bootstrapping (e.g.,
{it:lincom}, {it:pretrendtest()}). {p_end}

{phang}{opt reps(#)}: Number of replicates used; default 100.  {p_end}

{phang}{opt bootsave(name[, replace])}: Save the full set of bootstrap replicate results to the
specified file.  {p_end}




{marker randinf}{...}
{dlgtab:Randomization Inference}

{phang}{opt randinf(shuffle|probabilities)}: Produce confidence intervals and p-values based on
randomization inference (SE and t remain analytic). 
The default behavior is to show a two-tailed randomization inference p-value (testing the sharp null
hypothesis of no effects for any cohorts in 
any period), coupled with confidence intervals based on quantiles of the randomization inference
estimate distribution subtracted from the main estimate (i.e., shift intervals).
If {it:shuffle} is specified, treatment onset dates are rerandomized across clusters, holding
constant the number of clusters per treatment onset date 
and the number of untreated clusters. If {it:probabilities} is specified, clusters are rerandomized
based on user-specified probabilities of falling in 
any given treatment onset group, allowing the number of clusters per onset date to
change. Note: {it:randinf()} requires users 
specify a clustering scheme where each cluster has one treatment onset date (or is untreated) and
either (a) all cohorts are subsets of clusters (including the case where the cluster variable is
the 
cohort variable {it:i}), or (b) all clusters are subsets of cohorts and cohorts group observations
by treatment onset date.  {p_end}

{phang}{opt ristr:ata(variable)}: Specify strata containing clusters such that, with
{it:randinf(shuffle)}, treatment onset dates and 
treatment status are rerandomized across clusters separately within each stratum. This ensures
treatment onset dates maintain their composition 
with respect to strata, which can be useful for ensuring rerandomizations have the same allocation
of clusters with respect to size or 
expected treatment-effect variance (a condition generally required for randomization inference to be
valid). {p_end}

{phang}{opt riprobs(spec)}: Specify the probability of assigning each cluster to each treatment
onset date (or untreated status), when using 
{it:randinf(probabilities)}. The syntax for {it:spec}
is: {it:(ttre level level ...) (cluster1 p1 p2 ...) ... (clusterN p1 p2 ...)}, 
where {it:ttre} is the treatment onset variable, each level is the numerical level specifying a
treatment onset period (including '.' for 
untreated), and where the ensuing cluster-specific groups give the assignment probabilities for
placing a given cluster into the 
respective listed level of the treatment onset variable. {p_end}

{phang}{opt riside(left|right)}: Swap from producing two-sided randomization inference p-values to
the given type of one-sided p-value. {p_end}

{phang}{opt riabs}: Swap from producing two-sided p-values based on percentiles of the randomization
inference distribution to calculating 
p-values based on the share of the randomization inference distribution larger than the observed
estimate in absolute value. The confidence 
intervals are also adjusted based on this, and made symmetric. {p_end}

{phang}{opt rit:stat}: Compute randomization inference test statistics based on the randomization
inference for t-statistics rather than 
that for coefficients. Comes at great cost with respect to estimation speed. {p_end}

{phang}{opt reps(#)}: Number of rerandomizations used; default 100.  {p_end}

{phang}{opt riex:act}: Compute every possible rerandomization, in place of {it:reps()}. {p_end}

{phang}{opt risave(name[, replace])}: Save the full set of rerandomization results to the specified
file.  {p_end}




{marker supporting}{...}
{dlgtab:Supporting Features}

{phang}{opt hist:ogramcohorteffects}: Produce a histogram of cohort i-specific treatment effects;
return those effects in a matrix for the user to manipulate. {p_end}

{phang}{opt save:plots(string[, replace])}: Save all graphs generated to .gph files beginning with
the name specified in the option. Also saves .dta files beginning with the same
name containing the data used to generate any requested event study and continuous treatment effects
plots. Graphs and datasets generated have a suffix appended to the name specifying the exact
graph/dataset
generated. {p_end}

{phang}{opt sum:marystats}: Compute and present summary statistics of the outcome variable y and, if
present, of the continuous treatment variable contreat.
Summary statistics are computed overall, within the treatment group by time period (post-treatment
vs pre-treatment), and within the control group. If a
fixed reference period is not being used and pre-treatment relative time period effects are
requested, statistics on the size of the reference group relative to the rest of the
pre-treatment period will also be shown. If event study effects are requested, the composition of
each event study horizon will also be shown. {p_end}

{phang}{opt pois:son}: Implement the specified {it:wooldid} command using poisson
regression. Treatment effects
are presented in terms of differences in linear predictions (the poisson {it:predict(xb)}
option). Poisson regression is especially useful in situations where one would like to use
log(y) as an outcome, but where y contains natural 0s, as well as when treatment and control groups
are only on parallel trends in proportions (e.g., one group is always 2x the other). {p_end}



{marker eventstudy}{...}
{dlgtab:Event Study Estimation}

{phang}{opt esf:ixedbaseperiod}: Use a fixed reference period for estimation, rather than a pooled
reference period.   {p_end}

{phang}{opt esr:elativeto(-integer)}: The relative time period used as the fixed reference period,
if a fixed reference period is being used. The default is -1, the time period immediately
prior to treatment.  {p_end}

{phang}{opt espre:length(integer)}: The number of pre-treatment relative time period-specific
effects to be estimated. If any pre-treatment effects are requested and
{it:esfixedbaseperiod} is specified, all relative time period effects between the
{it:esrelativeto()} period and treatment onset will be presented.{p_end}

{phang}{opt espost:length(integer)}: The number of post-treatment relative time period-specific
effects to be estimated. {p_end}

{phang}{opt pretrend:test(powerlevels)}: Reports Roth (2022) pre-trend power diagnostics for each
requested power level. Requires {it:esfixedbaseperiod} and {it:esprelength()}. Reported results
include:
the slope of the line fit through the pre-treatment period-specific event horizon coefficients, the
joint test of the hypothesis that all pre-treatment coefficients are 
zero, the probability (given the observed variance-covariance structure of the pre-treatment
coefficients) that at least one pre-treatment coefficient 
would be statistically significant (at the 5% level) in the absence of a pre-trend, and a table
showing - given the observed variance-covariance structure of the 
pre-treatment coefficients - (a) the minimum detectable trend (mdt) that causes the joint test of
the pre-treatment coefficient null to reject (at the 5% level) at 
the user-specified power level, (b) the mdt that causes at least one pre-treatment coefficient to be
rejected (at the 5% level) at the user-specified power level, and 
(c)/(d) the biases in the main treatment effect implied by those smallest detectable
pre-trends.  {p_end}

{phang}{opt makep:lots}: Produce event study plots. {p_end}

{phang}{opt oos:treatedcontrols(integer)}: By default, {it:wooldid} uses cohorts treated after the
end of the sample the same way as never-treated control cohorts (i.e., such cohorts are
"out-of-sample treated controls"). When this option is specified with a value of T, cohorts treated
within T periods of the end of sample are instead handled the same as other treated cohorts.
Such cohorts are potentially allowed to contribute to pre-treatment relative time period-specific
effects when this option is specified,
but not otherwise. Note that if {it:esfixedbaseperiod} is specified, the out-of-sample-treated
cohorts will be dropped if the specified reference period is not observed
within the sample for these cohorts. {p_end}


{pstd}{ul on}Two Approaches to Estimation -- Pooled Reference Period vs Fixed Reference
Period:{ul off} Difference-in-differences treatment effect estimates involve comparing average
outcomes during
treated cohorts' treatment periods to treated cohorts' average outcomes during a pre-treatment
reference period. By default, {it:wooldid} does comparisons using a pooled reference period,
meaning that treated periods for a given cohort are compared to a reference period containing all
pre-treatment observations for the cohort. If the user specifies the option
{it:esfixedbaseperiod}, {it:wooldid} will switch to employing a fixed reference period, meaning that
all treatment effect estimates will be formed by comparing treated cohorts' average outcomes
in post-treatment periods to their outcomes in a particular pre-treatment reference period
(specified in relative time period terms). By default, this pre-treatment reference period is the
time period immediately prior to treatment (relative period -1), though the user can specify
alternative behavior.  {p_end}

{pstd} With respect to which type of reference period should be used, fixed reference periods offer
the virtue of transparency: all treated cohorts are compared to one reference period,
which is the same for all cohorts in relative time terms. By contrast, pooled reference periods can
be more opaque, with different cohorts having reference periods that contain
varying numbers of time periods. In exchange for reduced transparency, pooled reference periods
offer the possibility of greater precision in treatment effect estimates, owing
to the fact that they pool together more data and thus have greater opportunity to average out
idiosyncratic variation in outcomes within the reference period. {p_end}

{pstd} While there is no general solution to this transparency-precision tradeoff, use of a pooled
reference period is likely disadvantageous whenever the outcome variable
is evolving over time within the pre-treatment period. It also tends to be disadvantageous when one
is seeking to estimate more pre-treatment event study style treatment effects
than can be estimated for all cohorts, as can occur if one has staggered treatment onset and wants
to estimate the largest number of pre-treatment coefficients possible. In this case,
cohorts with fewer observed pre-treatment periods will have their earliest observed period set as
their reference period, even if this period is a period for which treatment effects are
reported for other cohorts. If this issue combines with cohorts not being on the same outcome
trends, the resulting event study plots may show treatment effect time paths that evolve in
counterintuitive and difficult to interpret ways, possibly obscuring even severe parallel trends
assumption violations from detection. {p_end}

 

{marker interactions}{...}
{dlgtab:Subgroup Effects}

{phang}{opt sub:groups()}: Specifies the subgroup (SG) variable, whose levels specify subgroups for
which {it:wooldid} will compute subgroup specific estimates.
The SG variable must consist of non-negative integer values (or missing values) and its levels must
specify collections of cohorts i. For estimates to be valid, a given cohort i cannot
be split across multiple subgroups. Cohorts with all missing values are not placed in any subgroup
and are ignored when performing subgroup computations.   {p_end}


{pstd} {ul:Description of Approach to Subgroup Effect Estimation:} Given a variable dividing cohorts
i into different subgroups, {it:wooldid} will estimate treatment effects separately
for each subgroup (i.e., each level of the SG variable). This is achieved without modification to
the underlying
regression specification in step 1 (see {help wooldid##description:here}). Rather, it is done by
carrying out the margins calculations in step (2) within the subsamples of the data
specified by the subgroup variable.
Additionally, {it:wooldid} computes all possible pairwise comparisons of the overall
effects (and pre-treatment effects) estimates by SG (e.g., SG 1 vs 2, 1 vs 3, 2 vs 3, etc.),
presenting for each comparison the difference in effect estimates and a p-value for each difference.
A joint test of the hypothesis that all SGs share the same effect estimate (normalized to be the
effect estimate of the numerically lowest numbered SG) is also computed.  {p_end}



{marker continuous}{...}
{dlgtab:Continuous Treatments}

{phang}{opt contreat(varname)}: The continuous treatment variable. Must be defined for all
observations (set to 0 if not naturally defined within untreated cohort-periods).{p_end}

{phang}{opt contreatelas:ticity}: Compute AME-type estimates as dyex semielasticities, i.e. the
level change in y/% change in contreat. Achieved by modifying the calculation done at the 
margins calculation stage. Successful estimation of these elasticities tends to require a larger
sample, as well as y and contreat variables that are not sparsely distributed. 
Note: with {it:poisson}, yields regular elasticities rather than semielasticities. {p_end}

{phang}{opt contreatw:ithin}: Specify this when you wish to use only within-treated cohort variation
in the continuous treatment variable. This option disables a number of estimability checks 
that would otherwise cause {it:wooldid} to stop estimation in the absence of any untreated cohorts,
and causes {it:wooldid} to skip sample cleaning routines that drop time
periods for which a control cohort is not present.   {p_end}

{phang}{opt contreatp:oly(integer)}: Includes polynomial terms in {it:contreat} of the specified
degree. Be cautious using higher degree polynomials. {p_end}

{phang}{opt lattice(latticevar OR latticetype latticenum)}: Specifies a lattice variable, a variable
specifying groups of observations across which the effect of {it:contreat} is
allowed to vary, even within cohort-periods and within the reference period. If just the name of a
variable is passed to {it:lattice()}, this variable will be used as the lattice variable.
If instead a lattice type and a lattice number are passed, this option will generate a lattice
variable using {it:latticenum} quantiles, with the quantiles being formed within groups
according to the {it:latticetype} (see below). Any specified sample weights will be used when
constructing the lattice. If {it:lattice} is not specified, {it:wooldid} will proceed without
using a lattice variable (as though you specified a lattice variable that is always 1). {p_end}

{phang} {it:latticetype:} Specifies groups in which to estimate {it:latticenum} quantiles of
{it:contreat}. Choice of lattice type can be important for both interpretive reasons
(e.g., when interpreting the output of {it:makecfxplots}) and for estimation feasibility
reasons: clustered standard error computation will often fail if multiple clusters are not
represented at each lattice point. Lattice type options include: {p_end}
{synoptset 20 tabbed}{...}
{synopt : {it:all}} Quantiles are formed looking at the full, within-sample distribution of the
continuous treatment variable, including values for untreated observations. {p_end}
{synopt : {it:byTime}} Quantiles are formed using the distribution of the continuous treatment
variable within each time period t separately. Observations will thus be placed in the
lowest quantile bin if their value of {it:contreat} is in the lowest 1/{it:latticenum} of the
distribution for their time period, even if their value would be high in the sample
overall.  {p_end}
{synopt : {it:byTreated}} Quantiles are formed separately in the treated and untreated
groups. {p_end}
{synopt : {it:byTtre}} Quantiles are formed separately for the control group and for each treatment
onset cohort. {p_end}
{synopt : {it:byCohort}} Quantiles are formed separately within each cohort i. {p_end}
{synopt : {it:byTreatedByTime}} Quantiles are formed by time period t, separately for the treated
and untreated groups.  {p_end}
{synopt : {it:byTtreByTime}} Quantiles are formed by time period t, separately for the control group
and each treatment onset cohort. {p_end}
{synopt : {it:byCohortByTime}} Quantiles are formed separately by cohort and time period t. {p_end}
{synopt : {it:baseTreated}} Quantiles are formed within the treated group. These quantile range
definitions are then used to define lattice points for the rest of the sample. {p_end}
{synopt : {it:baseTreatedRef}} Quantiles are formed within the reference period for the treated
group. These quantile range definitions are then used to define lattice points
for the rest of the sample. {p_end}
{synopt : {it:baseTreatedPost}} Quantiles are formed within the treated period for the treatment
group. These quantile range definitions are then used to define lattice points
for the rest of the sample. {p_end}
{synopt : {it:baseTreatedPre}} Quantiles are formed within the full pre-treatment period for the
treatment group. These quantile range definitions are then used to define lattice
points for the rest of the sample. {p_end}

{phang}{opt latticeig:noreweights}: Ignore sample weights when computing the quantiles for a
generated lattice. {p_end}

{phang}{opt makec:fxplots}: Produce plots presenting the average marginal effect of {it:contreat}
and the average treatment effect within each bin of the lattice variable or, if
specified, of the variable in {it:cfxlattice()}. These plots are only sure to be valid for the ATT
estimand type, but may be valid for others if lattice groupings
do not crosscut cohorts (or, less stringently, cohort-periods, for ATT_it). {p_end}

{phang}{opt cfxlattice(varname)}: The bins in which effects are calculated for the continuous
treatment effects plots; by default, the lattice variable. {p_end}

{phang}{opt cfxplott:ypes(pre post integers)}: The cfxplots to be produced: {it:post} (the default)
for the post-treatment period, {it:pre} for the pre-treatment period, and any integers
for just the specified periods relative to treatment. {p_end}

{phang}{opt makecfxplotsbysg}: Produce continuous treatment effects plots for each group specified
in the subgroup variable.   {p_end}

{phang}{opt cfxplotnose}: Do not estimate standard errors for the cfxplots; speeds
estimation. {p_end}

{phang}{opt contreatc:ontrols(controltype)}: Specifies what controls for {it:contreat} will be
added: the specified controls are always equal to f(Contreat_jt) (the function used to capture
  the effect of contreat within treated i-t cells) interacted with a set of fixed effects-like
indicator variables (time period, cohort, etc.). The specified controls will partial out
the effect of {it:contreat} in a particular way, affecting the nature of the reference effect of
{it:contreat} to which estimated effects of {it:contreat} in treated cohort-periods are
compared. This option must be specified to proceed with a continuous treatment, and care must be
taken when specifying it: different options correspond with very different
substantive assumptions. See more on this topic in the description section below. {p_end}

{phang} {it:contreatcontroltype}: untreatedZero; or one or more of the following: basic, byTreated,
byTtre, byCohort, and byTime  {p_end}
{synopt : {it:basic}} Controls for the effect of contreat (i.e., for f(Contreat_jt)) are added
without any interactions.  {p_end}
{synopt : {it:byTreated}} Controls for the effect of contreat are interacted with indicators for
membership in the treatment group and in the control group. {p_end}
{synopt : {it:byTtre}} Controls for the effect of contreat are interacted with indicators for
membership in each treatment onset group and in the control group.  {p_end}
{synopt : {it:byCohort}} Controls for the effect of contreat are interacted with cohort (i)
membership indicators.  {p_end}
{synopt : {it:byTime}} Controls for the effect of contreat are interacted with time period (t)
indicator variables. {p_end}
{synopt : {it:untreatedZero}} Special option to formulate the regression to work when a continuous
treatment variable is always equal to 0 outside of treated cohort-periods. Setting
the continuous treatment variable to always equal 0 outside of treated cohort-periods can be an
appropriate way to handle cases where the variable is unknowable or undefined outside
of treated cohort-periods. If this option is specified when the always 0 condition is not met,
{it:wooldid} will return poorly identified estimates. {p_end}

{pstd} {ul:On Specifying contreatcontrols() and lattice:} The regression equation {it:wooldid}
adopts when estimating
the effects of a continuous treatment variable is:
Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it *
f(Contreat_jt)) + g(f(Contreat_jt)) {p_end}

{pstd} The idea behind this approach is that, just as each i-t cell is allowed its own treatment
effect when a binary treatment is specified, here, each i-t cell is allowed its own treatment
effect intercept and its own marginal effect of the continuous treatment variable. If cohort-period
cells are large enough that concerns remain that the effect of Contreat_jt may be non-linear
within those cells, the function f() can be used to allow for use of a more complicated function of
the effect of Contreat_jt within i-t cells. Note that continuous treatment designs require 
that cohort-period cells contain multiple clusters for use of analytic clustered standard error
estimation. {p_end}

{pstd} The options {it:contreatpoly()} and {it:lattice()} determine the contents of the function
f(), which in general is a stepwise function of contreat. If no lattice is specified,
f(Contreat_jt) = beta_it1*Contreat_jt^1 + beta_it2*Contreat_jt^2 + ... + beta_itp*Contreat_jt^p,
where p is the polynomial order specified in {it:contreatpoly(p)}. If a lattice variable is
specified, this is modified to f(Contreat_jt) = sum_over_l_in_L[L_l + L_l * Poly(Contreat_jt)],
where L_l is an indicator variable that is one if a given observation falls in
 lattice bin l and is 0 otherwise and where Poly(Contreat_jt) is a term capturing the user specified
number of polynomial terms of Contreat_jt. {p_end}

{pstd} It is important to note that specifying a lattice has two effects on the function f(). First,
it allows the marginal effect of Contreat_jt to vary within i-t cells, whereas before
it could not. Second, it adds additional treatment effect intercept terms that exist within i-t
cells. A similar effect to including a lattice thus could be achieved via an alternative,
finer definition of the cohort variable i. This has some implications for interpreting average
marginal effect estimates in the presence of a lattice. AME-type estimates do not capture
the impact of discontinuities in the effect of Contreat_jt -- more precisely, they do not capture
the effect of the intercept terms a lattice variable might add to f() (nor that of the
intercept terms associated with the cohort-period indicator variables themselves). When a specified
lattice has many groups and almost all of the effect of Contreat_jt is captured by the
intercept lattice terms alone, the identification of marginal effects of Contreat_jt can be
threatened. AME and ATT estimates can also be discordant under certain circumstances: for example,
if within each i-t cell, Y = -Contreat + 100*floor(Contreat/10), the entire positive effect of
Contreat could be loaded onto the intercept terms and yield a positive ATT coupled with a
negative AME. In such cases, it may be wise to either use a coarser lattice or to focus on average
treatment effect estimates, which do capture the impact of those intercept terms --
this is especially the case when considering continuous treatment effects plots. {p_end}
 {pstd} The function g(f(Contreat_jt)) is a control function of Contreat_jt intended to serve as the
continuous treatment-equivalent of the binary treatment's fixed effects. As such, g()
 determines to a significant degree what source of variation is used to identify the effect of
Contreat_jt. This leads to a complex set of choices: since average marginal effects of
 Contreat_jt can be identified using only within cohort-period variation in treatment intensity,
many approaches to identification exist. {p_end}
 {pstd} The first set of options for specifying g() include the control types {it:basic},
{it:byTreated}, {it:byTtre}, and {it:byCohort}. These options are analogous to specifying different
 types of unit fixed effects, and are implemented by including interactions between f(Contreat_jt)
and then either a constant ({it:basic}), a treatment vs control group indicator ({it:byTreated}),
 a set of treatment onset group indicators ({it:byTtre}), or a set of cohort indicators
({it:byCohort}). For example, when {it:basic} is specified, g(f(Contreat_jt)) = f(Contreat_jt),
while
 when {it:byCohort} is specified, g(f(Contreat_jt)) = sum_over_i[fe_i * f(Contreat_jt)]. Since the
treatment effects of Contreat_jt are estimated using interactions between treated i-t cell
 indicator variables and f(Contreat_jt), these different g() options set the reference effect of
Contreat_jt relative to which Contreat_jt effects are estimated in treated cells. {p_end}

{pstd} In addition to the above control types, {it:byTime} can be specified to add interactions
between time period indicator variables and f(Contreat_jt) to g(). These time period
interactions can serve a role analogous to that of the second diff in the binary
difference-in-differences case, partialling out time period specific effects of Contreat_jt. The
control type
{it:byTime} can be specified on its own, or in conjunction with another control type. For example,
if both of the control types {it:byCohort} and {it:byTime} are specified,
one would have g() = sum_over_i[fe_i*f(Contreat_jt)] + sum_over_t[fe_t*f(Contreat_jt)].  {p_end}

{pstd} If the continuous treatment variable is not meaningfully defined outside of treated
cohort-periods, the control type {it:untreatedZero} can be used to specify g() = 0. By
not including a g() function at all, identification of marginal effects of Contreat_jt will occur
using only within treated cohort-period variation. Note: specifying {it:untreatedZero}
makes other technical changes to the estimation process in order to facilitate this case and should
not be specified outside of this narrow situation. {p_end}

{pstd} If attempting to identify the average marginal effects of a continuous treatment variable
without a strictly untreated control group (i.e., by exploiting within treated cohort
variation in treatment intensity), the option {it:contreatwithin} should be specified. This will
prevent {it:wooldid} from objecting to the presence of time periods without a control group
and to the presence of i-t fixed effects (if specified). When specified, {it:wooldid} will still
present ATT-type average treatment effects estimates, which in this case should essentially
be the integral of the estimated average marginal effects across the distribution of the continuous
treatment variable. Note: one should not specify the control type {it:byTime} in this case,
since in the absence of an untreated group, the time period interactions will force {it:wooldid} to
partial out genuine treatment effects for some group, threatening identification. {p_end}


{marker controls}{...}
{dlgtab:Controls and Fixed Effects}

{phang}{opt control:s()}: The contents of this option are passed directly into the regression
specification as controls, without any modification. {p_end}

{phang}{opt fe()}: Specify fixed effects. The contents of this option are passed directly to
{it:reghdfe} as fixed effects (into its absorb statement) without modification. 
Variables given will be interpreted as categorical variables, though i.X##c.X2 type notation is also
permitted for passing heterogeneous slope controls and the like.  {p_end}

{phang}{opt timetrend:s}: Control for cohort specific time trends in the regression (i.e., for
i.i##c.t). Importantly, this option is required for heterogeneous time trend controls to work,
even when they are specified manually, as it also updates {it:wooldid}'s various estimability
checks, drops cohorts without sufficient pre-treatment data to estimate a time trend,
and, in the fixed reference period case, ensures that the two reference periods that must be
normalized to 0 to estimate the time trend are chronologically adjacent. This option does not
automatically include any time trends that are a function of any included continuous treatment
variable.  {p_end}

{phang}{opt coarse:cohortcontrols()}: Includes the specified controls as well as their interactions
with all Tz_it dummies included in the model: one interaction with each
cohort-time period for which a treatment effect is estimated (if {it:esfixedbaseperiod} is
specified, every estimable treated cohort by relative time period effect). These interactions use
versions of the controls that have been de-meaned within i-t cells (and set to 0 outside of treated,
non-reference period i-t cells). The effect of these controls is partialled out from the
reported treatment effects. Unlike simply specifying cohort by time period interactions with a given
control, this pools the effect of the controls in the reference period and the control
group. All coarse cohort controls must be specified using factor notation ("c." for continuous, "i."
for categorical), so that {it:wooldid} knows how to interact them. {p_end}

{phang}{opt inter:activecontrols()}: Includes the specified controls and their interactions with the
Tz_it dummies, as in {it:coarsecohortcontrols()}, except rather than partialling out the effect of
the controls, the treatment effects are permitted to vary linearly with the de-meaned controls. The
stage 1 regression differs: the controls enter as i.t#c.X rather than pooled, their interactions
with the Tz_it dummies use the controls de-meaned within cohort i rather than within i-t cells, and
the controls are given without factor notation. In the stage 2 margins calculation, when treatment
is turned on,
{it:coarsecohortcontrols()} ignores the interaction between the i-t cell treatment dummies and these
controls, whereas {it:interactivecontrols()} includes the interaction term as part of
the ATT (etc.) calculation. Relative to {it:coarsecohortcontrols()}, this specification most closely
mimics the approach to including covariates that appears in Wooldridge (2021). When using only one
interactive control, ATTs
should be similar to those obtained when specifying that variable as a continuous treatment variable
and then separately including controls that interact it with time dummy variables. {p_end}


{pstd} {ul:Caution Regarding Continuous Treatment Variables and Controls:} If manually specifying
additional controls as a function of a continuous treatment variable, do not include the
continuous treatment variable itself among your controls, as this can adversely affect estimates
produced in {it:wooldid}'s margins calculation stage; a control with identical values but a
different variable name is okay. Such controls are best specified with {it:contreatcontrols(basic)},
or any option other than {it:untreatedZero}. {p_end}


{marker technical}{...}
{dlgtab:Technical Features}

{phang}{opt update}: If specified, the command {it:wooldid,} {it:update} will update your local
installation of {it:wooldid} based on the most recent edition of the program available on
github. {p_end}

{phang}{opt fast:nose}: Forgo standard error estimation for an increase in estimation speed. {p_end}

{phang}{opt regtol(#)}: Specifies the convergence tolerance to be used by {it:reghdfe} and
{it:ppmlhdfe}. By default, {it:reghdfe} sets a value of 8; {it:wooldid} has a default of 9.
You may wish to set this higher or lower; higher may be advisable in settings where you ask
{it:reghdfe} to absorb interaction terms with many continuous variables -- as tends
to occur when using continuous treatment variables. {p_end}

{phang}{opt ccno:absorb}: Forces continuous treatment controls to be included in the main regression
statement, rather than in the {it:absorb()} option, which is the default
whenever using the commands {it:reghdfe} or {it:ppmlhdfe} to estimate the underlying
regression. Under certain conditions, this can improve estimation stability and avoid problems
with variables being inappropriately dropped from the regression. {p_end}

{phang}{opt ccc_absorb()}: An alternate version of coarsecohortcontrols() that admits only
continuous variables and passes the controls to the {it:absorb()} statement when using
a high dimensional fixed effects estimator. This is often less efficient than the alternative
approach in terms of computation time, but can be important in terms of memory
efficiency in the main regression or under certain other circumstances. {p_end}

{phang}{opt nofastcollapse}: Forces all stage 2 margins calculations of treatment effects and
standard errors to occur at the observation-level, rather than at the cohort-period (i-t) 
level. This significantly slows estimation time, but may be necessary if the ATT and related
estimands cannot be estimated just based on i-t level summaries. This occurs automatically when  
{it:interactivecontrols()}, {it:contreatpoly()}, or controls with the continuous treatment variable
inside them are specified. However, manually specifying this option may be necessary 
in some cases.  {p_end}

{phang}{opt rseed(#)}: Specify the random number seed, for replicability of bootstrapping,
randomization inference, and simultaneous bands results.  {p_end}

{phang}{opt ver:bose}: Prints the regression command {it:wooldid} uses and prints progress reports
throughout the estimation process.
This option also stores the underlying regression model used by {it:wooldid} to
{it:___wooldidfullmodel}; this model can be accessed using {it:estimates restore}. {p_end}

{phang}{opt iteffects}: Produce the full table of cohort-period (i-t) level effects, labelled by
cohort, period, etc. Not available with continuous treatments. {p_end}

{phang}{opt safety(string)}: When {it:safety(off)} is not explicitly specified, {it:wooldid} will
caution the user or halt estimation when {it:wooldid} is likely to yield unreliable estimates.
Most important among these checks are tests for whether or not a variable intended to capture the
treatment effect in a given i-t cell was dropped from the underlying regression for
multicollinearity, as may occur if one's fixed effects or continuous treatment control variables
partial out one of the main effects. Identification failures of this sort will generally be
detected,
though not necessarily when the problem concerns continuous treatments being analyzed with a lattice
or with polynomial terms.   {p_end}




{marker output}{...}
{title:Working With Output}

{pstd}
{it:wooldid} generates output primarily in the form of matrices and scalars stored in e(), though
all main results can also be accessed using standard utilities like {it:esttab}.
The commands {it:est store}, {it:est restore}, and {it:est save} are compatible with {it:wooldid}
for saving its output. The full set of {it:wooldid} output in e() can be viewed using the command
{it:ereturn list}. {p_end}

{pstd} The primary results are stored in the following matrices: {p_end}
{synoptset 30 tabbed}{...}
{synopt:{it:e(wdidmainresults)}} All primary and event study results.{p_end}

{synopt:{it:e(wdidsgresults_{it:sgcase})}} Primary and event study results for the subgroup
{it:sgcase}, if subgroups are used.{p_end}

{synopt:{it:e(pwsgs)}} Pairwise differences between the subgroup-specific treatment effects, if
subgroups are used.{p_end}

{synopt:{it:e(b)}, {it:e(V)}} The estimates (e.g., main_att, ES_att_post0) and their
covariance. Without {it:bootstrap()}, e(V) is 0 off the diagonal except within one estimand type's
event study coefficients for one (sub)sample, so {it:test} and {it:lincom} are not valid across
other estimates.{p_end}

{pstd} In addition to these core results, the following materials are also returned in e(): {p_end}

{synopt: SG Effect Homogeneity Tests} For each estimate type (ATT, ATT_it, AME, etc.), the p-value
of the joint test that all subgroups have the same such post-treatment effect
(namely, the same effect as the lowest numbered SG) is in the scalar {it:e(joint_att_allsgssame)}
and analogously named scalars; for pre-treatment effects, in scalars like
{it:e(joint_pre_att_allsgssame)}.{p_end}

{synopt:{it:e(N)}, {it:e(N_dropped)}, {it:e(Nclus)}} Estimation sample observation count, number of
observations dropped from initially specified sample, and cluster count. {p_end}

{synopt:R^2 measures} For standard estimation, e(r2), e(r2adj), e(r2_within), and e(r2_withinadj);
for poisson estimation, the pseudo-R^2 e(r2pseud). {p_end}

{synopt:{it:e(histogramestimates)}} The cohort i-specific treatment effects used to generate
histograms, when requested.{p_end}

{synopt:{it:e(iteffects)}} The cohort-period (i-t) level effects, when requested.{p_end}

{synopt:{it:e(summarystats)}} Summary statistics, when requested. Some of the most commonly
referenced are also returned directly as scalars into e(). {p_end}

{synopt:{it:e(roth_power_att)}} The {it:pretrendtest()} table for ATT (likewise for att_it, ame, and
ame_it); the slope and joint test p-value are in {it:e(roth_slope_att)} and
{it:e(roth_prep_att)}.{p_end}

{synopt: CFX Plot Info} The information used to produce each continuous treatment effects plot, in a
matrix e(cfx[item]_[type]): [item] is one of att, att_i, att_it, or att_itime
(average treatment effects) or ame, ame_i, ame_it, or ame_itime (marginal effects of the continuous
treatment); [type] is cfxpost (post-treatment period), cfxpre (pre-treatment
period), or relyrZ (relative time period Z, where the letter n represents a negative sign). {p_end}

{synopt:{it:e(sample)}} The estimation sample used, useful for identifying which observations
{it:wooldid} drops during its sample formation and sanitation procedures.{p_end}

{synopt:{it:e(cmdline)}} The input you fed to {it:wooldid} on the command line.{p_end}

{synopt:randinf() results} {it:e(ri_p)}, {it:e(ri_ci)} (the shift intervals), and {it:e(ridist_b)}
(one row per rerandomization), aligned with {it:e(b)}.{p_end}

{synopt:{it:e(estimate_type)}} Whether your treatment effect estimates use a pooled reference period
or a fixed reference period. {p_end}

{synopt:{it:e(b_es_att)}, {it:e(V_es_att)}} The ATT event study coefficients and their covariance
(e.g., for {it:honestdid}); {it:e(numpre_att)} is the number of pre-treatment coefficients. Likewise
for att_it, ame, and ame_it, and for subgroups (e.g., e(b_es_att_s1) for subgroup 1).{p_end}

{synopt:Bootstrap results} {it:e(bootdist)} (one row per replicate) and the bootstrap intervals and
p-values in {it:e(boot_ci_*)} and {it:e(boot_p_*)}.{p_end}

{synopt:{it:Other Items}} {it:e(ll)}, {it:e(ll_0)}, {it:e(rmse)}, {it:e(rss)}, {it:e(F)}, and
{it:e(chi2)} (poisson only) are passed through directly
from the underlying regression estimator, wherever available. {p_end}

{pstd} Outside of e(), {it:saveplots()} writes its files to STATA's working directory, and
{it:verbose} stores the stage 1 underlying regression as the model {it:___wooldidfullmodel}
(see those options). {p_end}

{pstd} {ul:Tips for accessing output and remaking graphs:} {it:matlist} prints one of {it:wooldid}'s
output matrices. {it:svmat} adds one to the data in memory as variables that you can
directly manipulate (e.g., to remake the event study graphs): for example,
{it:svmat e(histogramestimates), names(col)} adds the cohort i-specific treatment effects, named
after the
matrix's columns. The program {help event_plot} can also generate event study plots from {it:e()} if
pointed to the results of interest (e.g., for ATT-type event study coefficients,
{it:stub_lag(ES_att_post#)} and {it:stub_lead(ES_att_pre#)}); it does not automatically insert 0s to
represent the missing reference period. {p_end}




{marker examples}{...}
{title:Examples}

{phang}0) Conventional diff-in-diff, targeting the ATT and clustering standard errors by
unit:{p_end}

{p 10 10} {cmd:wooldid y i t ttre, cluster(unit) att}

{phang}1) Conventional diff-in-diff, targeting both ATT and ATT_it, with production of a histogram
showing cohort-i specific effects: {p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att att_it histogramcohorteffects}

{phang}2) Conventional diff-in-diff, with added controls and fixed effects:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, controls(controlvar1) fe(i.fe_bonus)}

{phang}3) Conventional diff-in-diff, using coarse cohort controls and including cohort-specific time
trends:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre,}{...}
{cmd: coarsecohortcontrols(c.controlvar1 i.categoricalcontrolvar) timetrends }

{phang}4) Event study with 3 pre-treatment and 3 post-treatment estimates using a fixed reference
period and generating event study plots, as well as 
implementing the full pre-trend testing suite for power levels 80 and 90: {p_end}

{p 10 10}   {cmd:wooldid y i t ttre, espre(3) espost(3) esfixedbaseperiod makeplots}{...}
{cmd: pretrendtest(80 90)  }

{phang}5) Event study with 3 pre-treatment and 3 post-treatment estimates using a fixed reference
period set to 2 periods before treatment, rather than the default
of the period immediately prior to treatment:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, espre(3) espost(3) esfixedbaseperiod esrelativeto(-2)}

{phang}6) Diff-in-diff with separate treatment effects calculated by treatment onset cohort:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, subgroups(ttre)}

{phang}7) Difference-in-differences with a treatment intensity variable (one that is 0 outside the
treatment period), estimated using a stepwise-linear function and a lattice with 5 quantiles
(formed based on the treated group's treatment period), and generating continuous treatment effects
plots for the post treatment period:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(intensity) contreatcontrols(untreatedZero)}{...}
{cmd: makecfxplot lattice(baseTreatedPost 5)}

{phang}8) Difference-in-differences with a continuous treatment variable that varies within all
cohorts, with a lattice containing 5 quantiles (formed separately by treatment group and
by period), controlling for the continuous treatment variable by treatment group and by time period,
with continuous treatment effects plots for the pre-treatment period:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(contreatvar) contreatcontrols(byTreated byTime)}{...}
{cmd: makecfxplot lattice(byTreatedByTime 5) esprelength(3) cfxplottypes(pre)}



{marker bugs}{...}
{title:Bugs and Assistance}

{pstd}{ul:Known Issues:} {p_end}

{p 10 10}   Estimation with cohort-specific time trends sometimes encounters numerical stability and
related issues. 
It is advised that you do not manually specify time-trends variables for inclusion
in the fixed effects (i.e., in the {it:absorb()} statement when using {it:reghdfe}). Estimate
stability seems greater when the option {it:esfixedbaseperiod} is specified.

{p 10 10}   The automatic tests for model identification that {it:wooldid} implements are not fully
comprehensive: {it:wooldid} will almost always succeed in identifying
when a treatment effect has been partialled out by a fixed effect in the binary treatment case and
will often succeed in doing so when a continuous treatment variable is specified
without a lattice variable or any polynomial terms. However, when a continuous treatment is
specified with a lattice and/or a polynomial term, and especially in the presence of an
unbalanced panel, {it:wooldid} will be unlikely to detect problems related to important treatment
effects being partialled out by controls or fixed effects. As such, one should
be very careful when specifying one's fixed effects and {it:contreatcontrols} when using a
continuous treatment.

{pstd}If you encounter a bug - be it the program crashing or just the program delivering a weird
result - please reach out to me about it. This program is still
under development, so assistance finding and correcting bugs is much appreciated. I am also happy to
help and offer guidance on how to use this program.{p_end}

{pstd} Note: if you encounter unusual, hard to track down errors, try uninstalling and reinstalling
{it:reghdfe} and {it:ftools} first. This can fix some problems. {p_end}



{marker citations}{...}
{title:Citations}

{phang} Correia, Sergio. reghdfe: Stata module for linear and instrumental-variable/gmm regression
absorbing multiple levels of fixed effects. Statistical Software Components s457874, 2017.
Boston College Department of Economics. https://ideas.repec.org/c/boc/bocode/s457874.html {p_end}

{phang} Ibragimov, Rustam and Muller, Ulrich K. t-Statistic Based Correlation and Heterogeneity
Robust Inference. Journal of Business & Economic Statistics, 2010. {p_end}

{phang} Roth, Jonathan. Pretest with Caution: Event-Study Estimates after Testing for Parallel
Trends. American Economic Review: Insights, 2022. {p_end}

{phang} Wooldridge, Jeffrey M. Two-Way Fixed Effects, the Two-Way Mundlak Regression, and
Difference-in-differences Estimators. Working Paper,
2021. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3906345 {p_end}



{marker author}{...}
{title:Author & How to Cite}

{pstd} Author: Thomas A. Hegland (Agency for Healthcare Research and Quality)  {p_end}

{pstd} Contact: thomas.hegland@ahrq.hhs.gov, thomashegland.com, @thomas_hegland  {p_end}

{pstd} Citation: Hegland, Thomas A. wooldid: Estimation of Difference-in-Differences Treatment
Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way
Fixed Effects Regressions. Statistical Software Components s459238, 2023.{p_end}



{marker disclaimer}{...}
{title:Disclaimer}

{pstd}This document reflects only the views of the author and does not necessarily represent the
views of the Agency for Healthcare Research and Quality, the Department
of Health and Human Services, or the United States government. The program with which this help file
is associated is the work of the author and does not come
with any endorsement by or warranty from the Agency for Healthcare Research and Quality, the
Department of Health and Human Services, or the United States government.  {p_end}
