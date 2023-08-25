{smcl}
{* *! version 1.0 2023-08-24}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "regress" "help regress"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{vieweralsosee "poisson" "help poisson"}{...}
{vieweralsosee "gtools" "help gtools"}{...}
{vieweralsosee "margins" "help margins"}{...}
{vieweralsosee "hdidregress" "help hdidregress"}{...}
{vieweralsosee "did_imputation" "help did_imputation"}{...}
{viewerjumpto "Syntax" "wooldid##syntax"}{...}
{viewerjumpto "Description" "wooldid##description"}{...}
{viewerjumpto "All Options" "wooldid##options"}{...}
{viewerjumpto "Specifying the Estimand(s)" "wooldid##specifyestimand"}{...}
{viewerjumpto "Standard Errors" "wooldid##standarderrors"}{...}
{viewerjumpto "Supporting Features" "wooldid##supporting"}{...}
{viewerjumpto "Event Studies" "wooldid##eventstudy"}{...}
{viewerjumpto "Subgroup Effects" "wooldid##interactions"}{...}
{viewerjumpto "Continuous Treatments" "wooldid##continuous"}{...}
{viewerjumpto "Controls and Fixed Effects" "wooldid##controls"}{...}
{viewerjumpto "Technical Features" "wooldid##technical"}{...}
{viewerjumpto "Working with Output" "wooldid##output"}{...}
{viewerjumpto "Comparison with Related Estimators" "wooldid##comparison"}{...}
{viewerjumpto "Examples" "wooldid##examples"}{...}
{viewerjumpto "Bugs and Assistance" "wooldid##bugs"}{...}
{viewerjumpto "Citations" "wooldid##citations"}{...}
{viewerjumpto "Author Info How to Cite" "wooldid##author"}{...}
{viewerjumpto "Disclaimer" "wooldid##disclaimer"}{...}
{title:Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way Fixed Effects Regressions}

{pstd} {bf:wooldid} -- Suite for implementing difference-in-differences style analyses with staggered treatment onset using the two-way fixed effects approach proposed in Wooldridge (2021)
and the high dimensional fixed effects estimators developed by Correia (2017). Features include estimation of various types of average treatment effects, comparison of treatment effects
across subgroups, estimation of event study style estimates and related parallel trends assumption tests, and estimation of average and marginal treatment effects in the presence of a
continuous treatment (or treatment intensity) variable [experimental]. This suite also can produce event study plots, histograms of treatment effects by cohort, and plots of treatment
effects across the distribution of a specified continuous treatment variable. Options for inference include inference using standard analytic clustered standard errors, inference using
standard errors on average treatment effects that account for sampling variation in covariates (unconditional standard errors),
and Ibragimov and Muller (2010) style cluster-robust inference.  {p_end}

{pstd}{ul:Shortcuts:} {help wooldid##description:Description}, {help wooldid##options:Options}, {help wooldid##output:Working with Output},
{help wooldid##comparison:Comparison with Related Estimators}, {help wooldid##examples:Examples},
{help wooldid##bugs:Bugs and Assistance}, {help wooldid##citations:Citations}, {help wooldid##author:Author Info and How to Cite}, {help wooldid##disclaimer:Disclaimer} {p_end}

{pstd}{ul:Dependencies:} {cmd:wooldid} requires users have STATA 16 or higher installed, {help reghdfe} version 5.7.3, {help ftools} (for {it:reghdfe}),
and (optionally, for speed) {help gtools}. Poisson estimation also requires the package {help ppmlhdfe}. {p_end}

{pstd}{ul:Github Link:} {browse "https://github.com/thegland/wooldid/"} - the command {it:wooldid, update} will cause this program to update itself using the materials 
at this link. {p_end}



{marker syntax}{...}
{title:Syntax}

{phang}{cmd: wooldid} {it:y i t ttre} [if] [aw/pw=weights] [{cmd:,} {help wooldid##specifyestimand:estimands} {help wooldid##standarderrors:standard_errors}
{help wooldid##supporting:supporting_features} {help wooldid##eventstudy:event_study} {help wooldid##interactions:subgroup_effects} {help wooldid##continuous:continuous_treatment}
{help wooldid##controls:controls_and_FE} {help wooldid##technical:technical_features} ]{p_end}

{synoptset 22 tabbed}{...}
{synopt : {it:y}}Outcome variable. {p_end}
{synopt : {it:i}}Cohort identifier; must consist of integers greater than 0. Cohorts can be fine (i = individual units) or coarse (i = collections of individual units), but should be no coarser
than to consist of treatment onset cohorts. Coarser cohort definitions tend to enable faster estimation, and you can specify coarse cohorts
 alongside unit-level fixed effects. Simulations suggest the estimator is most statistically efficient when {it:i} is set at the level at which you cluster. {p_end}
{synopt : {it:t}}Time period identifier; must consist of integers greater than 0.{p_end}
{synopt : {it:ttre}}Variable containing cohort-specific dates of treatment onset; missing indicates never-treated cohorts, cohorts treated before the sample begins are dropped. {p_end}

{phang} {ul:Options by Purpose:} {p_end}
{synoptset 1 tabbed }{...}
{synopt : {help wooldid##specifyestimand:Specifying the Estimand(s):}} {opt att} {opt att_it} {opt att_i} {opt att_itime} {opt semi:elasticitiy} {p_end}

{synopt : {help wooldid##standarderrors:Standard Errors:}} {opt clus:ter(variable(s)[numerical])}  {opt rob:ust} {opt uncond:itionalse} {opt customvce(customstring)} {opt imp:values}  {p_end}

{synopt : {help wooldid##supporting:Supporting Features:}} {opt hist:ogramcohorteffects} {opt save:plots(filename_stem)} {opt replace} {opt sum:marystats}
{opt pois:son} {opt poisexp:results}  {p_end}

{synopt : {help wooldid##eventstudy:Event Study Estimation:}}  {opt esf:ixedbaseperiod} {opt espre:length(int>=0)} {opt espost:length(int>=0)} {opt esr:elativeto(-int)}
 {opt joint:tests}   {opt adv:ersarial} {opt makep:lots} {opt oos:treatedcontrols(int>=0)}  {p_end}

{synopt : {help wooldid##interactions:Subgroup Effects:}} {opt subgroups(variable[>=0 int])} {opt suppress:sgprimaryeffects}  {p_end}

{synopt : {help wooldid##continuous:Continuous Treatments:}} {opt contreat(variable[numerical])} {opt contreatelas:ticitytype(dydx, eydx, dyex, or eyex)}  {opt contreatw:ithin}
{opt contreatp:oly(+int)} {opt lattice(variable[+int] OR latticetype latticenum)} {opt latticeig:noreweights} {opt makec:fxplots} {opt cfxlattice(variable[+int])}
{opt cfxplott:ypes(post, pre, and/or integers specifying relative time periods)} {opt makecfxplotsbysg}  {opt cfxplotnose}  {opt contreatc:ontrols(conttreatcontroltype)}  {p_end}
{phang3} {it:-latticetype}: {it:all}, {it:byTime}, {it:byTreated}, {it:byTtre}, {it:byCohort}, {it:byTreatedByTime}, {it:byTtreByTime}, {it:byCohortByTime}, {it:baseTreated}, {it:baseTreatedRef},
{it:baseTreatedPost}, or {it:baseTreatedPre}  {p_end}
{phang3} {it:-latticenum}: integer > 1  {p_end}
{phang3} {it:-conttreatcontroltype}: {it:untreatedZero}, or one or more of {it:basic}, {it:byTreated}, {it:byTtre}, {it:byCohort}, {it:byTime}  {p_end}

{synopt : {help wooldid##controls:Controls and Fixed Effects:}} {opt control:s(fv varlist)} {opt fe(fv varlist)} {opt timetrend:s} {opt coarse:cohortcontrols(fv-required varlist)}   {p_end}

{synopt : {help wooldid##technical:Technical Features:}}
{opt update} {opt customw:eightmultiplier(variable[>=0 numerical])}  {opt clean:matrices} {opt regtol(+int)} {opt ccno:absorb} {opt ccc_absorb(varlist[numerical])}
{opt ver:bose} {opt safety(off)}  {p_end}



{marker description}{...}
{title:Description of Estimation Approach}

{pstd}
{bf:wooldid} estimates various types of difference-in-differences treatment effects in the vein of Wooldridge (2021), with extensions to handle continuous treatments and estimation of effects by
subgroup or treatment arm. The program does not currently accommodate cases where treatment is reversible -- once treated, cohorts must remain treated.{p_end}

{pstd}
Broadly speaking, {it:wooldid} estimates difference-in-differences treatment effects using a two-step procedure. In step 1, {it:wooldid} estimates the underlying two-way fixed effects regression,
which consists of (at baseline) a regression of the outcome variable on cohort fixed effects, time fixed effects, and a set of indicator variables used to flexibly capture the effect of treatment.
Each indicator variable is 1 in a specific cohort-period (i-t) cell, and 0 otherwise. The slate of indicator variables spans all cohort-periods for which treatment effects are to be estimated. In
step 2, {it:wooldid} completes a set of calculations using the {help margins} command that convert the large slate of coefficients on the treatment indicator variables from the underlying
regression in step 1 into treatment effects of interest. These margins calculations proceed by computing, for all treated observations, the difference between the predicted values from
the underlying regression and the counterfactual predicted values that would be observed given each observation was not treated. Various averages of these differences are then presented to obtain
various types of user-specified average treatment effects. To flesh this out, the underlying regression used is as follows: {p_end}

{p2col 5 8 8 0 : Step 1 (Underlying Regression):} Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it) {p_end}

{pstd} where i indexes cohorts, t indexes time periods, j indexes individual units (possibly equivalent to cohorts, possibly more numerous than cohorts), and where Z is the set of all (i,t)
pairs where i is a treated cohort and t is a treated period for cohort i. In the above regression, Tz_it is a dummy variable that is 1 if, for a given z=(iz,tz) in Z, i==iz and t==tz,
TreatedFlag_it is a variable that is 1 if an observation is part of an (i,t) cell in Z and 0 otherwise, fe_i is a cohort fixed effect, fe_t is a time fixed effect,
and Y_jt is the outcome variable. {p_end}

{pstd} The output of the underlying Step 1 regression is then used for margins calculations in step 2. The exact margins calculation performed depends on the estimands
specified by the user. The overall average treatment effect would be obtained by computing fitted values from the underlying regression among all treated observations (those where
  TreatedFlag_it is equal to 1), subtracting the counterfactual fitted values for these observations that set TreatedFlag_it to 0, and then taking the average of these differences
  (applying any sample weights). To target estimands other than the overall average treatment effect, this procedure can be modified by either
  introducing a set of weights (or modifying the pre-existing sample weights) when averaging the computed differences in the margins calculation or by restricting the sample so that the
  average only contains observations from a subgroup or time period of interest. The margins calculation can also be modified here to compute effects in terms of semi-elasticities.
For more on the different estimand-types available, see {help wooldid##estimands:estimands}. {p_end}
 {pstd} With respect to inference, there are several approaches to standard error estimation available. The default estimates standard errors in the underlying regression and then obtains
 standard errors on the treatment effects produced by the margins calculation using the delta method. These errors are conditional standard errors in the sense that they are conditional
 on the data: the difference in fitted values potentially depends on the value of various covariates, and this approach takes these values as fixed when performing the margins calculation.
 Given these standard errors, p-values are then obtained using a t-distribution of same degrees of freedom as the residual number of degrees of freedom in the underlying regression model.
 The alternative is to have {it:margins} produce unconditional standard errors (see the option {it:unconditional} in {help margins}), which account for sampling variation in covariates.
 For more on standard error computation, see {help wooldid##standarderrors:standard_errors}. {p_end}

{pstd} This overall approach can be extended to cover many additional use cases with little or no modification. The user can switch from standard regression to poisson regression without
requiring any essential modification to the regression equation or margins procedure. Additional controls and fixed effects can be added to the underlying regression without requiring
adjustments to the margins calculation. Effects specific to particular subgroups may be estimated simply by performing the margins calculation within a specified subgroup. Where multiple
subgroups are present, {it:wooldid} can use {help margins} to calculate pairwise comparisons of the subgroup-specific effects and compute joint tests of the hypothesis that all subgroups
have the same effect. This can be useful for tests of differences in treatment efficacy across treatment arms and can be used to replicate functionality that interaction terms in a
traditional difference-in-differences would have provided. See {help wooldid##interactions:subgroups} for more. {p_end}

{pstd} In a similar vein, event study style estimates that present treatment effects specific to various periods relative to treatment can be obtained by performing the margins calculation
restricted to the subsample of observations corresponding with particular periods relative to treatment. If effects in periods prior to treatment are of interest, this requires modifying
the underlying regression with additional treatment indicator variables for cohort-periods that correspond with the relative time period of interest in treated cohorts. The regression equation
itself remains the same: only the set of (i,t) pairs appearing in Z changes for this approach. The margins calculations involved also remain similar, though may be performed on varying
subsets of the data by relative period. When computing event study estimates, {it:wooldid} can supplement the standard treatment effect estimates with joint tests of the hypothesis that all
pre-treatment relative time period effects are zero, as well as of the hypothesis that all relative time period effects fall on a single line (namely, the line of best fit through the
estimated relative time period effects). For more on this topic, see {help wooldid##eventstudy:event_study}.{p_end}

{pstd} In terms of implementing event study-style estimates, {it:wooldid} offers two approaches. The default pooled reference period approach is conceptually similar to estimating results like
normal, but with the treatment onset date for all treated cohorts pulled earlier in time by just enough periods to estimate all requested relative time period effects. With this approach, we
add no more (i,t) pairs to Z than is necessary and we allow the reference period to potentially pool multiple periods together for each cohort. Event study estimates derived from this approach
always have a normalized effect of 0 in the earliest period reported (the pooled reference period). The alternative approach is to use a fixed reference period, in which case Z is modified to
include all (i,t) pairs for treated cohorts, except the pair corresponding with the period immediately prior to treatment (or corresponding with some other user specified period). In this case,
the regression is saturated with i-t indicator variables and all treatment effects are estimated relative to the single specified reference period. {p_end}

{pstd} With greater modification, the approach above can be extended to handle continuous treatment variables and treatments that come with an associated measure of treatment intensity
or dosage. When cohorts are coarse -- i.e., there are multiple units within each cohort-period cell -- and when a continuous treatment variable varies
within i-t cells, estimates of the effect of the continuous treatment variable can be computed within each treated i-t cell. To ensure treatment effect
heterogeneity is adequately captured, a flexible function of the continuous treatment variable can be used to measure its effect within each i-t cell if needed. 
Note that this approach is a bit speculative and is not covered in Wooldridge (2021), though Wooldridge has 
{browse "https://twitter.com/jmwooldridge/status/1695115782739435922?s=61&t=pXC6R4euc7LaN6MMB5VApQ":informally suggested} 
handling continuous treatments using interactions between treated cell indicators and the continuous treatment variable. {p_end}

{pstd} The particular modification to the underlying regression used for estimation in the continuous treatment context is as follows: {p_end}
{p2col 5 8 8 0 : Step 1 (Continuous):} Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it * f(Contreat_jt)) + g(f(Contreat_jt)) {p_end}

{pstd} where Contreat_jt is a continuous treatment variable that potentially varies arbitrarily across observations, where f() is a flexible function of Contreat_jt, and where
g() consists of interactions between various fixed effects-like indicator variables and f(Contreat_jt). Average treatment effects can be computed using the same margins calculation
as before. The average
marginal effect of the continuous treatment variable can be computed similar to the main margins calculation, but comparing differences in fitted values after modifying the variable
Contreat_jt rather than modifying the variable TreatedFlag_it. In addition to overall and event study style estimates, {it:wooldid} can also compute average treatment effect estimates
at varying points along the distribution of the continuous treatment variable. It can also compute various types of elasticity estimates in lieu of standard estimates of the effect
of the continuous treatment variable. For more detail on options related to continuous treatments, see {help wooldid##continuous:continuous_treatment}. {p_end}

{pstd} A final feature available in {it:wooldid} is its implementation of Ibragimov and Muller (2010) style cluster robust inference. This approach estimates the same underlying
regression as normal, but then modifies the margins calculation significantly. In particular, for any treatment effect requested (be it an average treatment effect or an effect of
a continuous treatment variable), this approach has the margins calculation compute that effect separately within each cluster. The final estimate presented is then the simple
average of these cluster-specific estimates, with inference being conducted by t-testing the cluster-level estimates. For more on this approach and when it is appropriate to apply,
see {help wooldid##standarderrors:standard_errors}. {p_end}

{pstd}
{ul: Comment on Sample Formation, Cohorts Treated Prior to Sample Start:} The estimation sample {it:wooldid} uses may be more restricted than your overall sample: {it:wooldid}
will drop observations with missing data, observations with 0-valued weights, observations from treated cohorts not observed prior to treatment (or not observed in the user-specified
reference year if {it:esfixedbaseperiod} is specified), and observations from time periods no untreated
(including not yet treated) units are observed. If the user requests that cohort-specific time trends be included, this restriction will extend to dropping cohorts that are
not observed enough times to estimate a time trend using the pre-treatment data. {p_end}

{pstd}
{ul: Comment on Cohorts Treated After Sample End:} By default, {it:wooldid} handles cohorts treated after the end of the sample identically to cohorts that have never been treated. The
option {it:oostreatedcontrols()} allows the user to specify that cohorts of this type should instead be handled more like cohorts treated within-sample
and thus that these cohorts should be allowed to contribute to pre-treatment relative time period-specific estimates in the event study context.
{p_end}



{marker options}{...}
{title:Options}

{marker specifyestimand}{...}
{dlgtab:Specifying the Estimand(s)}

{phang}{opt att}: Estimate average treatment on the treated (ATT) effects; this option is the default and may be combined with options att_it, att_i, and att_itime to estimate multiple different
effect types. ATT effects are computed by performing the margins calculation (see step 2 in {help wooldid##description:Description} for more on the margins calculation)
using whatever weights (if any) are used in the main regression command.{p_end}

{phang}{opt att_it}: Estimate simple averages of the i-t cell-specific treatment effects (ATT_it). The treatment effect in each cohort-period is calculated via the margins calculation
using any weights that may apply to the observations falling in each i-t cell. The final ATT_it estimate is the simple, unweighted mean of all of the cohort-period-specific estimates.
This is done by reweighting the data used in the margins calculation (after estimation of the underlying regression) by w/W_it where w is an observation's sample weight and W_it is the
total sample weight associated with all observations appearing in a given i-t cell. If no weights are specified, w is taken to have a value of 1 for all observations. If cohorts are
fine and thus each i-t cell contains only one observation, this is equivalent to ignoring any specified sample weights when performing margins calculations.   {p_end}

{phang}{opt att_i}: Estimate the simple average of cohort-specific treatment effects (ATT_i). The treatment effect for each cohort is calculated using the margins calculation and any sample
weights. The final ATT_i estimate is the simple, unweighted mean of the cohort-specific treatment effects. ATT_i effects are computed by performing the margins calculation, but reweighting
observations by either w/W_ipre or w/W_ipost, where w is an observation's sample weight, W_ipre is the total sample weight associated with cohort-i in the pre-treatment period
(but outside the reference period), and W_ipost is the total sample weight associated with cohort-i in the post-treatment period. Each observation is then reweighted based on whether it
falls in the pre- or post-treatment period. If weights are not used, w is taken to have a value of 1 for all observations.{p_end}

{phang}{opt att_itime}: Estimate the simple average of cohort-specific simple average treatment effects (ATT_itime). Treatment effect calculation here begins like for {it:att_it}, obtaining
treatment effects for individual cohort-periods of interest using weights. Then, cohort-specific effects are obtained by taking the simple mean of the cohort-period-specific effects.
The overall presented ATT_itime estimate is then equal to the simple, unweighted mean of these cohort-specific effects. ATT_itime effects are calculated by performing the margins calculation,
but reweighting observations by w/(W_it * T_iprepost), where w is an observation's sample weight, W_it is the amount of sample weight falling in a given i-t cell, and T_iprepost contains
(depending on the time period) either the number of post-treatment time periods in which a cohort-i is observed or the number of pre-treatment, non-reference period time periods in which a
cohort i is observed. If weights are not used, w is taken to have a value of 1 for all observations.{p_end}

{phang}{opt semi:elasticity}: Compute treatment effect estimates in terms of semi-elasticities (% change in y / unit change in x). This is achieved by adjusting the margins calculation
to present the average percentage change, rather than the average level difference, between the two sets of predicted outcomes. This is implemented using the
{it:eydx} option of the {help margins} command. Note that for estimation to succeed, this option tends to require large samples and a y-variable that is not sparsely distributed. {p_end}


{marker standarderrors}{...}
{dlgtab:Standard Errors}

{phang}{opt clus:ter(clusvars)}: Variable(s) specifying clustering level. Clustered standard errors are computed analytically. By default,
regression estimation using this option will underlying rely on either {help reghdfe} or, for poisson, {help ppmlhdfe}. {p_end}

{phang}{opt rob:ust}: Heteroskedasticity robust standard errors; relevant only if cluster() is not specified. By default, regression estimation using this
option will underlyingly rely on either {it:reghdfe} or, for poisson, {it:ppmlhdfe}. {p_end}

{phang}{opt uncond:itionalse}: Compute standard errors that account for sampling variation in covariates. This option also forces {it:wooldid} to abandon
high dimensional fixed effect estimation packages for the underlying two-way fixed effect regression in favor of the commands {it:regress} and {it:poisson}.
As such, this option should be slower to use on average. This option requires specifying either robust or clustered standard errors. {p_end}

{phang}{opt customvce()}: Forces {it:wooldid} to switch to using either the {help regress} or {help poisson} command for underlying regression estimation,
and then passes the user specified text without modification into the {it:vce()} option for these commands. This allows specifying a broader set of
standard error calculation methods than allowed by the default high dimensional fixed effects estimation packages used by {it:wooldid}. {p_end}

{phang}{opt imp:values}: Implements Ibragimov and Muller (2010) style cluster-robust inference for all estimates. The idea behind the Ibragimov and Muller approach is to compute
cluster-specific treatment effect estimates, present their average as one's final treatment effect estimate, and conduct inference by t-testing the cluster-specific estimates against 0.
Using IM-style inference requires specifying a single cluster variable and is not compatible with the {it:unconditionalse} option. IM-style inference implicitly changes {it:wooldid}'s
target estimand. For example, if clustering at the cohort level, IM-style inference will cause {it:att} to yield the same effect estimates as {it:att_i}. Note that {it:wooldid} does not adjust the
weights used for computing different estimands when IM-style inference is requested: in general, IM-style inference will yield incorrect results for non-att estimands if clusters cross-cut cohorts
or (worse) cohort-time period cells. IM-style inference is not currently allowed with the use of subgroups, the option {it:jointtests}, and the option {it:makecfxplots}. {p_end}

{pstd} {ul: Discussion of Unconditional Standard Errors:}
The overall treatment effect estimates presented by {it:wooldid} depend on not just the coefficients from the underlying regression, but also on the observed values of any included covariates
or continuous treatment variables. The standard errors {it:wooldid} presents for these final estimates are produced during the margins calculation stage
using the delta method, taking the observed values of the relevant variables as given (i.e., conditional on those values). This contrasts with approaches that aim to account for sampling variation
in those observed values (i.e., unconditional standard errors; see the discussion of {it:vce(unconditional)} in {help margins}). Whether one should wish to take the observed values of one's
covariates as fixed likely varies from context to context. {p_end}

{pstd} Does using unconditional standard errors make a large difference? Wooldridge (2021) expresses the view that conditional and unconditional standard errors often should be similar.
Intuitively, this ought to be true whenever samples are large enough that sampling variation in covariates is ignorable (more precisely, ignorable within cohort-time period cells).
Unconditional standard errors may be more impactful, however, when samples are small, or when sample weights make particular observations very influential. Additionally, when
using unconditional standard errors, subgroup-specific estimates (as obtained using {it:subgroups()} or when computing event study estimates) can be sensitive to data from observations
outside the focal subgroup. On a technical basis, this is because {it:wooldid} uses the {it:subpopulation} features in {it:margins} when producing these subgroup estimates. {p_end}

{pstd} {ul: Caution Regarding Unconditional Standard Errors When Used with Certain Weights and Estimands:} When using aweights and estimands other than the ATT, tests using unconditional
standard errors sometimes appear to be incorrectly sized, rejecting the null hypothesis at extraordinarily high rates. The issue with non-ATT estimands likely is related to the margins
calculation doing computations on predicted values using different weights (weights modified to hit non-ATT estimands) than used in the initial regression. The issue with using aweights does
not always occur, though when it proves to be an issue it tends to result in poorly sized test statistics for both {it:wooldid} and {help hdidregress}. Given these issues, if {it:safety(off)}
is not specified, {it:wooldid} will warn against using unconditional standard errors with aweights and stop estimation if using them with a non-ATT estimand. Broadly speaking, unconditional
standard errors tend to be highly sensitive to one's choice of weights. {p_end}

{pstd} {ul: Discussion of Ibragimov and Muller-Style Inference:}
The virtue of the Ibragimov and Muller approach is that it can significantly improve on analytic clustered standard errors when the number of clusters is small and when the clusters are
very heterogeneous in size. A key cost to adopting this approach is the required assumption that cluster-specific treatment effect estimates are independent. This assumption seems most likely to
be satisfied when estimating the effect of continuous treatment variables using only within-cluster variation, or in a bespoke diff-in-diff set up where clusters contain both treated and control
units and where time fixed effects are not shared across clusters. In the standard diff-in-diff setting, {it:wooldid} will still present IM-style estimates by computing cluster-specific treatment
effect estimates, despite the fact that the existence of a shared control group (and perhaps a shared reference period, if {it:esfixedbaseperiod} is not being used) likely causes mischief with the
independence assumption.  {p_end}


{marker supporting}{...}
{dlgtab:Supporting Features}

{phang}{opt hist:ogramcohorteffects}: Produce a histogram of cohort i-specific treatment effects; return those effects in a matrix for the user to manipulate. {p_end}

{phang}{opt save:plots(string)}: Save all graphs generated to .gph files beginning with the name specified in the option. Also saves .dta files beginning with the same
name containing the data used to generate any requested event study plots. Graphs and datasets generated have a suffix appended to the name specifying the exact graph/dataset
generated. By default, everything is saved to STATA's working directory, though you should be able to pass a directory through this string along with your filename.{p_end}

{phang}{opt replace}: When saving materials to the file names specified in {it:saveplots()}, overwrite any existing files with the specified names. {p_end}

{phang}{opt sum:marystats}: Compute and present summary statistics of the outcome variable y and, if present, of the continuous treatment variable contreat.
Summary statistics are computed overall, within the treatment group by time period (post-treatment vs pre-treatment), and within the control group. If a
fixed reference period is not being used and pre-treatment relative time period effects are requested, statistics on the size of the reference group relative to the rest of the
pre-treatment period will also be shown. {p_end}

{phang}{opt pois:son}: Implement the specified {it:wooldid} command using poisson regression, achieved by switching from the regression command {it:reghdfe} to {it:ppmlhdfe}. Treatment effects
are presented in terms of differences in linear predictions (the poisson {it:predict(xb)} option). Poisson regression is especially useful in situations where one would like to use
log(y) as an outcome, but where y contains natural 0s, as well as when treatment and control groups are only on parallel trends in proportions (e.g., one group is always 2x the other). {p_end}

{phang}{opt poisexp:results}: When implementing poisson, present treatment effects in terms of differences in exponentiated linear predictions (the poisson "number of events" or {it:predict(n)}
option). Results of this sort will thus not have a log-type interpretation. {p_end}


{marker eventstudy}{...}
{dlgtab:Event Study Estimation}

{phang}{opt esf:ixedbaseperiod}: Use a fixed reference period for estimation, rather than a pooled reference period.   {p_end}

{phang}{opt esr:elativeto(-integer)}: The relative time period used as the fixed reference period, if a fixed reference period is being used. The default is -1, the time period immediately
prior to treatment.  {p_end}

{phang}{opt espre:length(integer)}: The number of pre-treatment relative time period-specific effects to be estimated. If any pre-treatment effects are requested and
{it:esfixedbaseperiod} is specified, all relative time period effects between the {it:esrelativeto()} period and treatment onset will be presented.{p_end}

{phang}{opt espost:length(integer)}: The number of post-treatment relative time period-specific effects to be estimated. {p_end}

{phang}{opt joint:tests}: Calculate joint statistical significance tests of all requested pre-treatment relative time period-specific effects. Caution: simulations suggest these
tests tend to be particularly sensitive to whether you have a sufficient number of clusters in your sample for standard error estimation. {p_end}

{phang}{opt adv:ersarial}: Tests the joint hypothesis that all estimated relative time period-specific effects fall on a single line, where the line is chosen rather adversarially.
Specifically, the line against which estimates are tested is that chosen by the OLS regression of {it:coef_rt} on {it:time_rt} where {it:coef_rt} is the estimated coefficient in
relative time period
{it:rt} and where {it:time_rt} is equal to the relative time period. Note that the variable {it:coef_rt} includes values of 0 for the reference period (or periods, if {it:timetrends} and
{it:esfixedbaseperiod} are specified). The final p-values presented are from the F-test of the the joint null hypothesis that all non-reference period relative-time period specific effects
are equal to the fitted values from the OLS best fit line. These tests will not be conducted if an insufficient number of relative time period-specific treatment effects are requested
(either overall, or within specific subgroups). {p_end}

{phang}{opt makep:lots}: Produce event study plots. {p_end}

{phang}{opt oos:treatedcontrols(integer)}: By default, {it:wooldid} uses cohorts treated after the end of the sample the same way as never-treated control cohorts (i.e., such cohorts are
"out-of-sample treated controls"). When this option is specified with a value of T, cohorts treated within T periods of the end of sample are instead handled the same as other treated cohorts.
Such cohorts are potentially allowed to contribute to pre-treatment relative time period-specific effects when this option is specified,
but not otherwise. Note that if {it:esfixedbaseperiod} is specified, the out-of-sample-treated cohorts will be dropped if the specified reference period is not observed
within the sample for these cohorts. {p_end}

{pstd}{ul: Two Approaches to Estimation -- Pooled Reference Period vs Fixed Reference Period:} Difference-in-differences treatment effect estimates involve comparing average outcomes during
treated cohorts' treatment periods to treated cohorts' average outcomes during a pre-treatment reference period. By default, {it:wooldid} does comparisons using a pooled reference period,
meaning that treated periods for a given cohort are compared to a reference period containing all pre-treatment observations for the cohort. If the user specifies the option
{it:esfixedbaseperiod}, {it:wooldid} will switch to employing a fixed reference period, meaning that all treatment effect estimates will be formed by comparing treated cohorts' average outcomes
in post-treatment periods to their outcomes in a particular pre-treatment reference period (specified in relative time period terms). By default, this pre-treatment reference period is the
time period immediately prior to treatment (relative period -1), though the user can specify alternative behavior. For more on how this is implemented, see {help wooldid##discussion:discussion}.{p_end}

{pstd} With respect to which type of reference period should be used, fixed reference periods offer the virtue of transparency: all treated cohorts are compared to one reference period,
which is the same for all cohorts in relative time terms. By contrast, pooled reference periods can be more opaque, with different cohorts having reference periods that contain
varying numbers of time periods. In exchange for reduced transparency, pooled reference periods offer the possibility of greater precision in treatment effect estimates, owing
to the fact that they pool together more data and thus have greater opportunity to average out idiosyncratic variation in outcomes within the reference period. {p_end}

{pstd} While there is no general solution to this transparency-precision tradeoff, use of a pooled reference period is likely disadvantageous whenever the outcome variable
is evolving over time within the pre-treatment period. It also tends to be disadvantageous when one is seeking to estimate more pre-treatment event study style treatment effects
than can be estimated for all cohorts, as can occur if one has staggered treatment onset and wants to estimate the largest number of pre-treatment coefficients possible. In this case,
cohorts with fewer observed pre-treatment periods will have their earliest observed period set as their reference period, even if this period is a period for which treatment effects are
reported for other cohorts. If this issue combines with cohorts not being on the same outcome trends, the resulting event study plots may show treatment effect time paths that evolve in
counterintuitive and difficult to interpret ways, possibly obscuring even severe parallel trends assumption violations from detection. {p_end}

{pstd}{ul: Invalid Estimands:} For relative-time period specific estimates, distinctions between the estimand types {it:att_i}, {it:att_itime}, and {it:att_it} collapse. The latter estimand
will yield estimates that weight each cohort's contribution to the relative-time period specific estimate equally; the former two are differentiated from this only in terms of how they treat i-t
cells not relevant to a given relative-time period specific estimate. As such, only {it:att} and {it:att_it} type results will be presented when relative-period specific estimates are requested.
For the {it:att_i} and {it:att_itime} concepts, {it:att_it} type results will be shown. {p_end}


{marker interactions}{...}
{dlgtab:Subgroup Effects}

{phang}{opt sub:groups()}: Specifies the subgroup (SG) variable, whose levels specify subgroups for which {it:wooldid} will compute subgroup specific estimates.
The SG variable must consist of non-negative integer values (or missing values) and its levels must specify collections of cohorts i. For estimates to be valid, a given cohort i cannot
be split across multiple subgroups. Observations given missing values are not placed in any subgroup and are ignored when performing subgroup computations.   {p_end}

{phang}{opt suppress:sgprimaryeffects}: Skips estimation of subgroup-specific treatment effects and computes only the differences in effects between SGs. This option should
speed estimation if these treatment effects are not otherwise of interest. {p_end}

{pstd} {ul: Description of Approach to Subgroup Effect Estimation:} Given a variable dividing cohorts i into different subgroups, {it:wooldid} will estimate treatment effects separately
for each subgroup (i.e., each level of the SG variable). This is achieved without modification to the underlying
regression specification in step 1 (see {help wooldid##description:here}). Rather, it is done by carrying out the margins calculations in step (2) within the subsamples of the data
specified by the subgroup variable.
Additionally, {it:wooldid} computes all possible pairwise comparisons of the overall
effects (and pre-treatment effects) estimates by SG (e.g., SG 1 vs 2, 1 vs 3, 2 vs 3, etc.), presenting for each comparison the difference in effect estimates and a p-value for each difference.
A joint test of the hypothesis that all SGs share the same effect estimate (normalized to be the effect estimate of the numerically lowest numbered SG) is also computed. Note: standard errors
for subgroup-specific estimates can fail when using clustered standard errors in conjunction with subgroups that contain few (or only one) cluster. {p_end}

{pstd} A special case to consider is the difference-in-differences-in-differences (triple diff) case. In the triple diff case, all cohorts are categorized into two groups, with both groups
existing among treated cohorts and control cohorts. The target estimand may be thought of as the difference in the diff-in-diff estimates obtained separately for each of the two groups. To
implement a triple diff design in {it:wooldid}, one can specify the variable containing the two triple diff groups as a SG variable and then manually specify inclusion of SG by period
fixed effects (i.e., specifying {it:subgroups(sgvarname) fe(i.sgvarname##i.t)}). The difference in SG level estimates returned will be the triple diff estimate. {p_end}


{marker continuous}{...}
{dlgtab:Continuous Treatments}

{phang}{opt contreat(varname)}: The continuous treatment variable. Must be defined for all observations (set to 0 if not naturally defined within untreated cohort-periods).{p_end}

{phang}{opt contreatelas:ticitytype(dydx/eyex/eydx/dyex)}: Specifies which type of treatment effect to compute. {it:eyex} corresponds with computing the elasticity
(% change in y/% change in contreat), {it:eydx} corresponds with semi-elasticities comparable to what can be produced for the main effect of treatment
(% change in y/unit change in continuous treatment), {it:dyex} gives the other semi-elasticity (change in y/% change in continuous treatment), and {it:dydx} is the default
(not an elasticity). Note that successful estimation of these elasticities tends to require a larger sample, as well as y and contreat variables that are not sparsely distributed. {p_end}

{phang}{opt contreatw:ithin}: Specify this when you wish to use only within-treated cohort variation in the continuous treatment variable. This option sets the option {it:safety} to off, turns off
estimability checks that would otherwise cause {it:wooldid} to stop estimation in the absence of any untreated cohorts, and causes {it:wooldid} to skip sample cleaning routines that drop time
periods for which a control cohort is not present.   {p_end}

{phang}{opt contreatp:oly(integer)}: Includes polynomial terms in {it:contreat} of the specified degree. Be cautious using higher degree polynomials. {p_end}

{phang}{opt lattice(latticevar OR latticetype latticenum)}: Specifies a lattice variable, a variable specifying groups of observations across which the effect of {it:contreat} is
allowed to vary, even within cohort-periods and within the reference period. If just the name of a variable is passed to {it:lattice()}, this variable will be used as the lattice variable.
If instead a lattice type and a lattice number are passed, this option will generate a lattice variable using {it:latticenum} quantiles, with the quantiles being formed within groups
according to the {it:latticetype} (see below). Any specified sample weights will be used when constructing the lattice. If {it:lattice} is not specified, {it:wooldid} will proceed without
using a lattice variable (as though you specified a lattice variable that is always 1). {p_end}

{phang} {it:latticetype:} Specifies groups in which to estimate {it:latticenum} quantiles of {it:contreat}. Choice of lattice type can be important for both interpretive reasons
(e.g., when interpreting the output of {it:makecfxplots}) and for estimation feasibility reasons: clustered standard error computation will often fail if multiple clusters are not
represented at each lattice point. Lattice type options include: {p_end}
{synoptset 20 tabbed}{...}
{synopt : {it:all}} Quantiles are formed looking at the full, within-sample distribution of the continuous treatment variable, including values for untreated observations. {p_end}
{synopt : {it:byTime}} Quantiles are formed using the distribution of the continuous treatment variable within each time period t separately. Observations will thus be placed in the
lowest quantile bin if their value of {it:contreat} is in the lowest {it:latticenum}-th percent for their time period, even if their value would be high in the sample overall.  {p_end}
{synopt : {it:byTreated}} Quantiles are formed separately in the treated and untreated groups. {p_end}
{synopt : {it:byTtre}} Quantiles are formed separately for the control group and for each treatment onset cohort. {p_end}
{synopt : {it:byCohort}} Quantiles are formed separately within each cohort i. {p_end}
{synopt : {it:byTreatedByTime}} Quantiles are formed by time period t, separately for the treated and untreated groups.  {p_end}
{synopt : {it:byTtreByTime}} Quantiles are formed by time period t, separately for the control group and each treatment onset cohort. {p_end}
{synopt : {it:byCohortByTime}} Quantiles are formed separately by cohort and time period t. {p_end}
{synopt : {it:baseTreated}} Quantiles are formed within the treated group. These quantile range definitions are then used to define lattice points for the rest of the sample. {p_end}
{synopt : {it:baseTreatedRef}} Quantiles are formed within the reference period for the treated group. These quantile range definitions are then used to define lattice points
for the rest of the sample. {p_end}
{synopt : {it:baseTreatedPost}} Quantiles are formed within the treated period for the treatment group. These quantile range definitions are then used to define lattice points
for the rest of the sample. {p_end}
{synopt : {it:baseTreatedPre}} Quantiles are formed within the full pre-treatment period for the treatment group. These quantile range definitions are then used to define lattice
points for the rest of the sample. {p_end}

{phang}{opt latticeig:noreweights}: When constructing a lattice, this option forces {it:wooldid} to ignore sample weights when computing quantiles for lattice production. {p_end}

{phang}{opt makec:fxplots}: Produce plots presenting the average marginal effect of {it:contreat} and the average treatment effect within each bin specified by either the lattice
variable or, if specified, the variable in {it:cfxlattice()}. These plots are only sure to be valid for the ATT estimand type, but may be valid for others if lattice groupings
do not crosscut cohorts (or, less stringently, cohort-periods, for ATT_it). {p_end}

{phang}{opt cfxlattice(varname)}: A variable containing the bins in which effects are to be calculated for the continuous treatment effects plots. If not specified,
the lattice variable will be used. {p_end}

{phang}{opt cfxplott:ypes(pre post integers)}: A list of types of cfxplots to be produced. If {it:post} is among the types specified, plots will be produced giving
treatment effects across the continuous treatment variable within the post-treatment period. If {it:pre} is among the types specified, plots will be produced within the pre-treatment period.
If any integers are specified, plots will be produced just within the specified periods relative to treatment. If nothing is specified, the option post is assumed. {p_end}

{phang}{opt makecfxplotsbysg}: Produce continuous treatment effects plots for each group specified in the subgroup variable.   {p_end}

{phang}{opt cfxplotnose}: Do not estimate standard errors for the cfxplots; speeds estimation. {p_end}

{phang}{opt contreatc:ontrols(controltype)}: Specifies what controls for {it:contreat} will be added: the specified controls are always equal to f(Contreat_jt) (the function used to capture
  the effect of contreat within treated i-t cells) interacted with a set of fixed effects-like indicator variables (time period, cohort, etc.). The specified controls will partial out
the effect of {it:contreat} in a particular way, affecting the nature of the reference effect of {it:contreat} to which estimated effects of {it:contreat} in treated cohort-periods are
compared. This option must be specified to proceed with a continuous treatment, and care must be chosen when specifying it: different options correspond with very different
substantive assumptions. See more on this topic in the description section below. {p_end}

{phang} {it:conttreatcontroltype}: untreatedZero; or one or more of the following: basic, byTreated, byTtre, byCohort, and byTime  {p_end}
{synopt : {it:basic}} Controls for the effect of contreat (i.e., for f(Contreat_jt)) are added without any interactions.  {p_end}
{synopt : {it:byTreated}} Controls for the effect of contreat are interacted with indicators for membership in the treatment group and in the control group. {p_end}
{synopt : {it:byTtre}} Controls for the effect of contreat are interacted with indicators for membership in each treatment onset group and in the control group.  {p_end}
{synopt : {it:byCohort}} Controls for the effect of contreat are interacted with cohort (i) membership indicators.  {p_end}
{synopt : {it:byTime}} Controls for the effect of contreat are interacted with time period (t) indicator variables. {p_end}
{synopt : {it:untreatedZero}} Special option to formulate the regression to work when a continuous treatment variable is always equal to 0 outside of treated cohort-periods. Setting
the continuous treatment variable to always equal 0 outside of treated cohort-periods can be an appropriate way to handle cases where the variable is unknowable or undefined outside
of treated cohort-periods. If this option is specified when the always 0 condition is not met, {it:wooldid} will return poorly identified estimates. {p_end}

{pstd} {ul:On Specifying contreatcontrols() and lattice:} As discussed {help wooldid##discussion:here}, the regression equation {it:wooldid} adopts when estimating
the effects of a continuous treatment variable is:
Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it * f(Contreat_jt)) + g(f(Contreat_jt)) {p_end}

{pstd} The idea behind this approach is that, just as each i-t cell is allowed its own treatment effect when a binary treatment is specified, here, each i-t cell is allowed its own treatment
effect intercept and its own marginal effect of the continuous treatment variable. If cohort-period cells are large enough that concerns remain that the effect of Contreat_jt may be non-linear
within those cells, the function f() can be used to allow for use of a more complicated function of the effect of Contreat_jt within i-t cells.{p_end}

{pstd} The options {it:contreatpoly()} and {it:lattice()} determine the contents of the function f(), which in general is a stepwise function of contreat. If no lattice is specified,
f(Contreat_jt) = beta_it1*Contreat_jt^1 + beta_it2*Contreat_jt^2 + ... + beta_itp*Contreat_jt^p, where p is the polynomial order specified in {it:contreatpoly(p)}. If a lattice variable is
specified, this is modified to f(Contreat_jt) = sum_over_l_in_L[L_l + L_l * Poly(Contreat_jt)], where L_l is an indicator variable that is one if a given observation falls in
 lattice bin l and is 0 otherwise and where Poly(Contreat_jt) is a term capturing the user specified number of polynomial terms of Contreat_jt. {p_end}

{pstd} It is important to note that specifying a lattice has two effects on the function f(). First, it allows the marginal effect of Contreat_jt to vary within i-t cells, whereas before
it could not. Second, it adds additional treatment effect intercept terms that exist within i-t cells. A similar effect to including a lattice thus could be achieved via an alternative,
finer definition of the cohort variable i. This has some implications for interpreting average marginal effect estimates in the presence of a lattice. AME-type estimates do not capture
the impact of discontinuities in the effect of Contreat_jt -- more precisely, they do not capture the effect of the intercept terms a lattice variable might add to f() (nor that of the
intercept terms associated with the cohort-period indicator variables themselves). When a specified lattice has many groups and almost all of the effect of Contreat_jt is captured by the
intercept lattice terms alone, the identification of marginal effects of Contreat_jt can be threatened. AME and ATT estimates can also be discordant under certain circumstances: for example,
if within each i-t cell, Y = -Contreat + 100*floor(Contreat/10), the entire positive effect of Contreat could be loaded onto the intercept terms and yield a positive ATT coupled with a
negative AME. In such cases, it may be wise to either use a coarser lattice or to focus on average treatment effect estimates, which do capture the impact of those intercept terms --
this is especially the case when considering continuous treatment effects plots. {p_end}
 {pstd} The function g(f(Contreat_jt)) is a control function of Contreat_jt intended to serve as the continuous treatment-equivalent of the binary treatment's fixed effects. As such, g()
 determines to a significant degree what source of variation is used to identify the effect of Contreat_jt. This leads to a complex set of choices: since average marginal effects of
 Contreat_jt can be identified using only within cohort-period variation in treatment intensity, many approaches to identification exist. {p_end}
 {pstd} The first set of options for specifying g() include the control types {it:basic}, {it:byTreated}, {it:byTtre}, and {it:byCohort}. These options are analogous to specifying different
 types of unit fixed effects, and are implemented by including interactions between f(Contreat_jt) and then either a constant ({it:basic}), a treatment vs control group indicator ({it:byTreated}),
 a set of treatment onset group indicators ({it:byTtre}), or a set of cohort indicators ({it:byCohort}). For example, when {it:basic} is specified, g(f(Contreat_jt)) = f(Contreat_jt), while
 when {it:byCohort} is specified, g(f(Contreat_jt)) = sum_over_i[fe_i * f(Contreat_jt)]. Since the treatment effects of Contreat_jt are estimated using interactions between treated i-t cell
 indicator variables and f(Contreat_jt), these different g() options set the reference effect of Contreat_jt relative to which Contreat_jt effects are estimated in treated cells. {p_end}

{pstd} In addition to the above control types, {it:byTime} can be specified to add interactions between time period indicator variables and f(Contreat_jt) to g(). These time period
interactions can serve a role analogous to that of the second diff in the binary difference-in-differences case, partialling out time period specific effects of Contreat_jt. The control type
{it:byTime} can be specified on its own, or in conjunction with another control type. For example, if the control types specified are {it:byCohort} and {it:byTime} are both specified,
one would have g() = sum_over_i[fe_i*f(Contreat_jt)] + sum_over_t[fe_t*f(Contreat_jt)].  {p_end}

{pstd} If the continuous treatment variable is not meaningfully defined outside of treated cohort-periods, the control type {it:untreatedZero} can be used to specify g() = 0. By
not including a g() function at all, identification of marginal effects of Contreat_jt will occur using only within treated cohort-period variation. Note: specifying {it:untreatedZero}
makes other technical changes to the estimation process in order to facilitate this case and should not be specified outside of this narrow situation. {p_end}

{pstd} If attempting to identify the average marginal effects of a continuous treatment variable without a strictly untreated control group (i.e., by exploiting within treated cohort
variation in treatment intensity), the option {it:contreatwithin} should be specified. This will prevent {it:wooldid} from objecting to the presence of time periods without a control group
and to the presence of i-t fixed effects (if specified). When specified, {it:wooldid} will still present ATT-type average treatment effects estimates, which in this case should essentially
be the integral of the estimated average marginal effects across the distribution of the continuous treatment variable. Note: one should not specify the control type {it:byTime} in this case,
since in the absence of an untreated group, the time period interactions will force {it:wooldid} to partial out genuine treatment effects for some group, threatening identification. {p_end}


{marker controls}{...}
{dlgtab:Controls and Fixed Effects}

{phang}{opt control:s()}: The contents of this option are passed directly to {it:reghdfe} as controls, without any modification. {p_end}

{phang}{opt fe()}: The contents of this option are passed directly to {it:reghdfe} as fixed effects (into its absorb statement) without modification. Variables given will be interpreted as
categorical variables, though i.X##c.X2 type notation is also permitted for passing heterogeneous slope controls and the like.  {p_end}

{phang}{opt timetrend:s}: Control for cohort specific time trends in the regression (i.e., for i.i##c.t). Importantly, this option is required for heterogeneous time trend controls to work,
even when they are specified manually, as it also updates {it:wooldid}'s various estimability checks, drops cohorts without sufficient pre-treatment data to estimate a time trend,
and, in the fixed reference period case, ensures that the two reference periods that must be normalized to 0 to estimate the time trend are chronologically adjacent. This option does not
automatically include any time trends that are a function of any included continuous treatment variable.  {p_end}

{phang}{opt coarse:cohortcontrols()}: Includes the specified controls as well as interactions between them with all Tz_it dummies included in the model, meaning one interaction with each
cohort-time period for which a treatment effect is estimated (if {it:esfixedbaseperiod} is specified, that includes every estimable treated cohort by relative time period effect).
This differs from simply specifying cohort
by time period interactions with a given control in that it pools the effect of the controls in the reference period and the control group. Note that all specified coarse cohort controls must be
specified using factor notation, so that {it:wooldid} knows if they are to be interacted as a continuous variable (prefixed with "c.") or as a categorical variable (prefixed with "i.").

{pstd} {ul:Caution Regarding Continuous Treatment Variables and Controls:} If seeking to manually specify additional controls as a function of a continuous treatment variable, be sure not
to include the continuous treatment variable itself among your controls, as this can adversely affect estimates produced in {it:wooldid}'s margins calculation stage. It is okay to include
a control with values identical to that of the continuous treatment variable: it just needs to have a different variable name. It is recommended that any custom such controls be specified
in conjunction with the option {it:contreatcontrols(basic)}, or any option other than {it:untreatedZero}. {p_end}


{marker technical}{...}
{dlgtab:Technical Features}

{phang}{opt update}: If specified, the command {it:wooldid,} {it:update} will update your local installation of {it:wooldid} based on the most recent edition of the program available on github. {p_end}

{phang}{opt customw:eightmultiplier(varname)}: When computing various treatment effects using the margins calculation (see step 2 {help wooldid##description:here}), this variable is used to
weight the observations (or if a weight would already be used, that weight is multiplied by this variable). This allows production of certain bespoke types of treatment effects. {p_end}

{phang}{opt clean:matrices}: If specified, {it:wooldid} will clear out all matrices storied in STATA's main memory and in mata. Occasionally, {it:wooldid} enters a state where it repeatedly
crashes for difficult to determine reasons (relating to {it:reghdfe} crashing and leaving certain content in mata). Doing this can sometimes recover {it:wooldid} from this problem. {p_end}

{phang}{opt regtol}: Specifies the convergence tolerance to be used by {it:reghdfe} and {it:ppmlhdfe}. By default, {it:reghdfe} sets a value of 8; {it:wooldid} has a default of 9.
You may wish to set this higher or lower; higher may be advisable in settings where you ask {it:reghdfe} to absorb interaction terms with many continuous variables -- as tends
to occur when using continuous treatment variables. {p_end}

{phang}{opt ccno:absorb}: Forces continuous treatment controls to be included in the main regression statement, rather than in the {it:absorb()} option, which is the default
whenever using the commands {it:reghdfe} or {it:ppmlhdfe} to estimate the underlying regression. Under certain conditions, this can improve estimation stability and avoid problems
with variables being inappropriately dropped from the regression. {p_end}

{phang}{opt ccc_absorb()}: An alternate version of coarsecohortcontrols() that admits only continuous variables and passes the controls to the {it:absorb()} statement when using
a high dimensional fixed effects estimator. This is often less efficient than the alternative approach in terms of computation time, but can be important in terms of memory
efficiency in the main regression or under certain other circumstances. {p_end}

{phang}{opt ver:bose}: Prints the regression command {it:wooldid} uses, displays raw output from the regression estimated, and prints progress reports throughout the estimation process.
This option also stores the underlying regression model used by {it:wooldid} to {it:___wooldidfullmodel}; this model can be accessed using {it:estimates restore}. {p_end}

{phang}{opt safety(string)}: When {it:safety(off)} is not explicitly specified, {it:wooldid} will caution the user or halt estimation when {it:wooldid} is likely to yield unreliable estimates.
Most important among these checks are tests for whether or not a variable intended to capture the treatment effect in a given i-t cell was dropped from the underlying regression for
multicolinearity, as may occur if one's fixed effects or continuous treatment control variables partial out one of the main effects. Identification failures of this sort will generally be detected,
though not necessarily when the problem concerns continuous treatments being analyzed with a lattice or with polynomial terms.   {p_end}



{marker output}{...}
{title:Working With Output}

{pstd}
{it:wooldid} generates output primarily in the form of matrices and scalars stored stored in e(), though all main results can also be accessed using standard utilities like {it:esttab}.
The commands {it:est restore} and {it:est save} are compatible with {it:wooldid} for saving its output. The full set of {it:wooldid} output in e() can be viewed using the command
{it:ereturn list}. {p_end}

{pstd} The primary results are stored in the following matrices: {p_end}
{synoptset 30 tabbed}{...}
{synopt:{it:e(wdidmainresults)}} All primary and event study results.{p_end}

{synopt:{it:e(wdidsgresults_{it:sgcase})}} Primary and event study results for the subgroup {it:sgcase}, if subgroups are used and presentation of their
primary results has not been suppressed.{p_end}

{synopt:{it:e(pwsgs)}} Pairwise differences between the subgroup-specific treatment effects, if subgroups are used.{p_end}


{pstd} In addition to these core results, the following materials are also returned in e(): {p_end}
{synopt: Joint Tests by Horizon} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented in the scalars e(joint_horiz)
(ATT), e(joint_horiz_i) (ATT_i), e(joint_Choriz) (AME; continuous treatments only), and e(joint_Choriz_i) (AME_i; continuous treatments only). {p_end}

{synopt: Joint Tests by Subgroup & Time Period} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented by
subgroup, when used, in the matrices: e(joint_horizon_bysg), e(joint_horizon_i_byssg), e(joint_cont_horizon_bysg), and e(joint_cont_horizon_i_bysg). {p_end}

{synopt: Adversarial Line Tests:} P-values from joint tests of the hypothesis that all (non-reference) relative-time period treatment effect estimates fall on the same line
 are found in the local macros {it:e(adversarial_pval_att)}, {it:e(adversarial_pval_ame)}, {it:e(adversarial_pval_att_s[SG])}, etc. The coefficient and intercept of the line used in
 these tests may be found in matrices named {it:e(adversarial_line_att)}, etc. {p_end}

{synopt: SG Effect Homogeneity Tests} For each estimate type (ATT, ATT_it, AME, etc.), the joint test that all subgroups have the same such post-treatment effect (namely, the
same effect as the lowest numbered SG) may be found in the scalar {it:e(joint_att_allsgssame)} and analogously named scalars for other treatment effect types. For pre-treatment effects,
the test can be found in objects with names like {it:e(joint_pre_att_allsgssame)}.{p_end}

{synopt:{it:e(N)}, {it:e(N_dropped)}, {it:e(Nclus)}} Estimation sample observation count, number of observations dropped from initially specified sample, and cluster count. {p_end}

{synopt:R^2 measures} For standard estimation, these include e(r2), e(r2adj), e(r2_within), e(r2_withinadj). For poisson estimation, these include e(r2pseud), a pseudo-R^2 measure. {p_end}

{synopt:{it:e(histogramestimates)}} The cohort i-specific treatment effects used to generate histograms, when requested.{p_end}

{synopt:{it:e(summarystats)}} Summary statistics on the outcome variable and (if present) continuous treatment variable, when requested.{p_end}

{synopt: CFX Plot Info} Information used to generate the continuous treatment effects plots. The information used to produce each plot will be stored in a matrix with the following naming
pattern: e(cfx[item]_[type]), where [item] is either one of att, att_i, att_it, or att_itime (for overall average treatment effects) or one of ame, ame_i, ame_it, or ame_itime (for marginal
effects of the continuous treatment). Meanwhile, [type] is either cfxpost (estimated in the post-treatment period), cfxpre (estimated in the pre-treatment period), or relyrZ (estimated in
relative time period Z, where the letter n represents a negative sign). {p_end}

{synopt:{it:e(sample)}} A function characterizing the estimation sample used; useful for identifying which observations are dropped by {it:wooldid} during its sample formation and sanitation
procedures.{p_end}

{synopt:{it:e(cmdline)}} The input you fed to {it:wooldid} on the command line.{p_end}

{synopt:{it:e(estimate_type)}} A macro recording whether your treatment effect estimates use a pooled reference period or a fixed reference period. {p_end}

{synopt:{it:Other Items:}} {it:e(ll)}, {it:e(ll_0)}, {it:e(rmse)}, {it:e(rss)}, {it:e(F)}, {it:e(chi2)} (poisson only), and {it:e(converged)} (poisson only) are passed through directly
from the underlying regression estimator, wherever available. {p_end}

{pstd} In addition to these objects returned in e(), if saveplots() is specified, a dataset containing the information used to produce event study plots along with .gph files containing
all graphs generated will be saved to a file with your specified name in STATA's working directory. {p_end}

{pstd} If the option {it:verbose} is specified, the stage 1 underlying regression that {it:wooldid} estimates will also be stored to a model called {it:___wdidfullmodel}, which can be
accessed using standard {it:estimates restore} type syntax. {p_end}

{pstd} {ul:Tips for accessing output:} Use the {it:matlist} command to print the contents of one of {it:wooldid}'s output matrices, for example the command {it:matlist e(wdidmainresults)} will
display the contents of the main results matrix. The command {it:dis} can be used to print the content of a particular scalar macro. To convert one of {it:wooldid}'s output matrices
into a dataset that you can directly manipulate (e.g., to remake the automatically generated graphs), use the {it:svmat} command. For example, the command
{it:svmat e(histogramestimates), names(col)} will create in memory a dataset containing the contents of {it:e(histogramestimates)}
(i.e., the cohort i-specific treatment effects), using as variable names the names of the matrix's columns. {p_end}

{pstd} {ul:Compatible Post-Estimation Commands:} {it:estimates save}, {it:estimates store}, {it:estimates restore}, {it:esttab}, and similar commands. {p_end}

{pstd} {ul:Incompatible Post-Estimation Commands:} {it:lincom}, {it:test}, and all other post-estimation commands that rely on receiving a valid variance-covariance matrix from
an estimator. While {it:wooldid} does return a coefficient vector ({it:e(b)}) and a variance-covariance matrix ({it:e(V)}), these are for display purposes only. The variance-covariance
matrix has all off-diagonal elements set to 0, meaning that while {it:esttab} and the like will present correct estimates using these inputs, efforts to combine or otherwise test different
estimates will not, in general, yield valid results {p_end}



{marker comparison}{...}
{title:Comparison to Related Estimators}

{pstd} {ul:Comparing Output to hdidregress twfe:}
The STATA 18 command {help hdidregress} has an option for implementing Wooldridge (2021) style estimates. Its implementation takes the Mundlak approach to estimation
and uses the command {help regress}, in conjunction with {help margins}, to compute estimates. Estimates from {it:hdidregress twfe} with the option {it:notyet}
(i.e., allowing use of not-yet treated cohorts as controls) should be comparable to estimates from {it:wooldid} using the {it:att} estimand-type, a pooled
reference period (i.e., not using {it:esfixedbaseperiod}), unconditional standard errors, and sample weights which (if present) do not vary within units.
Estimates produced by the two will still differ, however, when certain types of panel imbalance are present. For example, when a treatment onset group contains
multiple cohorts and some of those cohorts are not observed in periods when others are observed, {it:wooldid} and {it:hddiregress twfe} will yield different results.
In such cases, {it:wooldid} approaches the situation by presenting the average treatment effect among only the observed cohort-periods. For example, suppose cohorts 1 and 2
are in the same treatment onset group, cohort 1 is observed for the first 3 periods after treatment onset, and cohort 2 is only observed for the first 2. In this case, {it:wooldid}
and {it:hdidregress} will report the same treatment effects for the treatment onset group for the first 2 periods after treatment. The two programs will disagree on the third
effect: {it:wooldid} will report the same 3rd period effect that it (and {it:hdidregress}) would report if cohort 2 were not incldued in the dataset at all; {it:hdidregress} does not.
Finally, in order to get estimates using controls to
match between the two methods, a {it:wooldid} user will have to manually adjust the controls as needed (e.g., given a control that is constant within cohorts, a
{it:wooldid} user will need to specify i.t#c.control to replicate the behavior of including that control in {it:hdidregress}).  {p_end}

{pstd} {ul:Comparing Output to did_imputation:}
The Borusyak, Jaravel, and Spiess (2023) difference-in-differences imputation estimator is implemented in the package {help did_imputation}. Wooldridge (2021)
proves some equivalence results between their estimator and the Wooldridge estimator. For binary treatments, {it:wooldid} will yield the same estimates as
{it:did_imputation} given that one uses the {it:att} estimand-type, uses a pooled reference period (i.e., not using {it:esfixedbaseperiod}), does not
use unconditional standard errors, and does not use controls. With controls, {it:wooldid} and {it:did_imputation} will still yield the same results provided
that each cohort-time period cell contains only 1 observation. In this 'fine cohorts' case, the inclusion of i-t cell indicators ensures that the included controls
have their coefficients determined solely from the same portion of the sample as is used to estimate the regression in the BJS estimator. However, outside of special
cases, results with controls will not match if cohorts are coarse (contain multiple observations per i-t cell). To assist in this case, the {it:wooldid} option
{it:coarsecohortcontrols()}
will add to the {it:wooldid} regression your controls plus interactions between your controls and the i-t cell indicators associated with periods for which treatment
effects are requested. While the {it:coarsecohortcontrols()} option will not yield the same estimates as {it:did_imputation}, it arguably reflects an approach closer to
the spirit of BJS in that it increases the degree to which the effect of the controls is independent between the control/reference portion of the sample and
the portion of the sample where treatment effects are of interest. {p_end}



{marker examples}{...}
{title:Examples}

{phang}0) Conventional diff-in-diff, targeting the ATT and clustering standard errors by cohort:{p_end}

{p 10 10} {cmd:wooldid y i t ttre, cluster(i) att}

{phang}1) Conventional diff-in-diff, targeting both ATT and ATT_it, with production of histogram showing cohort-i specific effects: {p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att att_it histogramcohorteffects}

{phang}2) Conventional diff-in-diff, with added controls and fixed effects:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, controls(controlvar1) fe(i.fe_bonus)}

{phang}3) Conventional diff-in-diff, using coarse cohort controls and including cohort-specific time trends:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, coarsecohortcontrols(c.controlvar1 i.categoricalcontrolvar) timetrends }

{phang}4) Event study with 3 pre-treatment and 3 post-treatment estimates using a pooled reference period and generating event study plots, as well as testing the hypotheses
that all pre-treatment effects are 0 and that all time-period specific effects fall on the same line: {p_end}

{p 10 10}   {cmd:wooldid y i t ttre, espre(3) espost(3)  makeplots jointtests adversarial}

{phang}5) Event study with 3 pre-treatment and 3 post-treatment event study estimates, calculating all estimates relative a fixed reference period (fixed to the
period 2 periods before treatment):{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, espre(3) espost(3) esfixedbaseperiod esrelativeto(-2)}

{phang}6) Diff-in-diff with separate treatment effects calculated by treatment onset cohort:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, subgroups(ttre)}

{phang}8) Difference-in-differences with a treatment intensity variable (one that is 0 outside the treatment period), estimated using a stepwise-linear function and a lattice with 5 quantiles
(formed based on the treated group's treatment period), and generating continuous treatment effects plots for the post treatment period:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(intensity) contreatcontrols(untreatedZero) makecfxplot lattice(baseTreatedPost 5))}

{phang}9) Difference-in-differences with a continuous treatment variable that varies within all cohorts, with a lattice containing 5 quartiles (formed separately by treatment group and
by period), controlling for the continuous treatment variable by treatment group and by time period, with continuous treatment effects plots for the pre-treatment period:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(contreatvar) contreatcontrols(byTreated byTime) makecfxplot lattice(byTreatedByTime 5) cfxplottypes(pre))}



{marker bugs}{...}

{title:Bugs and Assistance}

{pstd}{ul:Known Issues:} {p_end}

{p 10 10}   Estimation with cohort-specific time trends sometimes encounters numerical stability and related issues. Identical specifications can yield slightly different underlying
coefficients depending on the underlying program used for estimation ({it:regress} vs {it:reghdfe}). It is advised that you do not manually specify time-trends variables for inclusion
in the fixed effects (i.e., in the {it:absorb()} statement when using {it:reghdfe}). Estimate stability seems greater when the option {it:esfixedbaseperiod} is specified.

{p 10 10}   Estimation with continuous treatment controls sometimes faces numerical stability issues as well. For more consistent results, specify either a higher {it:regtol()} or
            specify the option {it:ccnoabsorb} to prevent the continuous treatment controls from entering into {it:reghdfe}'s {it:absorb()} statement.

{p 10 10}   The automatic tests for model identification that {it:wooldid} implements are not fully comprehensive: {it:wooldid} will almost always succeed in identifying
when a treatment effect has been partialled out by a fixed effect in the binary treatment case and will often succeed in doing so when a continuous treatment variable is specified
without a lattice variable or any polynomial terms. However, a continuous treatment is specified with a lattice and/or a polynomial term, and especially in the presence of an
unbalanced panel, {it:wooldid} will be unlikely to detect problems related to important treatment effects being partialled out by controls or fixed effects. As such, one should
be very careful when specifying one's fixed effects and {it:contreatcontrols} when using a continuous treatment.

{p 10 10}   Under various circumstances, {it:wooldid} cannot estimate standard errors for all requested objects. In this case, {it:wooldid} will return point estimates without 
standard errors. This can cause the functionality of certain post-estimation tools to break since, in this situation, {it:wooldid} will not return an {it:e(V)} matrix. 

{pstd}If you encounter a bug - be it the program crashing or just the program delivering a weird result - please reach out to me about it. This program is still
under development, so assistance finding and correcting bugs is much appreciated. I am also happy to help and offer guidance on how to use this program.{p_end}

{pstd} Note: if you encounter unusual, hard to track down errors, try uninstalling and reinstalling {it:reghdfe} and {it:ftools} first. This can fix some problems. {p_end}



{marker citations}{...}

{title:Citations}

{phang} Borusyak, Kirill, Xavier Jaravel, and Jann Spiess. Revisiting Event Study Designs: Robust and Efficient Estimation. Working Paper, 2023. {p_end}

{phang} Correia, Sergio. reghdfe: Stata module for linear and instrumental-variable/gmm regression absorbing multiple levels of fixed effects. Statistical Software Components s457874, 2017.
Boston College Department of Economics. https://ideas.repec.org/c/boc/bocode/s457874.html {p_end}

{phang} Ibragimov, Rustam and Muller, Ulrich K. t-Statistic Based Correlation and Heterogeneity Robust Inference. Journal of Business & Economic Statistics, 2010. {p_end}

{phang} Wooldridge, Jeffrey M. Two-Way Fixed Effects, the Two-Way Mundlak Regression, and Difference-in-differences Estimators. Working Paper, 2021. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3906345 {p_end}



{marker author}{...}
{title:Author & How to Cite}

{pstd} Author: Thomas A. Hegland (Agency for Healthcare Research and Quality)  {p_end}

{pstd} Contact: thomas.hegland@ahrq.hhs.gov, thomashegland.com, @thomas_hegland  {p_end}

{pstd} Citation: Hegland, Thomas A. wooldid: Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way
Fixed Effects Regressions. Statistical Software Components, 2023.{p_end}



{marker disclaimer}{...}
{title:Disclaimer}
{pstd}This document reflects only the views of the author and does not necessarily represent the views of the Agency for Healthcare Research and Quality, the Department
of Health and Human Services, or the United States government. The program with which this help file is associated is the work of the author and does not come
with any endorsement by or warranty from the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government.  {p_end}
