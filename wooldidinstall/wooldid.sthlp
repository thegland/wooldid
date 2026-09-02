{smcl}
{* *! version 1.4 2026-08-19}{...}
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
{viewerjumpto "Runtime and What Drives It" "wooldid##runtime"}{...}
{viewerjumpto "Working with Output" "wooldid##output"}{...}
{viewerjumpto "Comparison with Related Estimators" "wooldid##comparison"}{...}
{viewerjumpto "Examples" "wooldid##examples"}{...}
{viewerjumpto "Bugs and Assistance" "wooldid##bugs"}{...}
{viewerjumpto "Citations" "wooldid##citations"}{...}
{viewerjumpto "Author Info How to Cite" "wooldid##author"}{...}
{viewerjumpto "Disclaimer" "wooldid##disclaimer"}{...}
{title:Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way Fixed Effects Regressions}

{pstd} {bf:wooldid} -- Suite for implementing difference-in-differences style analyses with staggered treatment onset using the two-way fixed effects approach proposed in Wooldridge (2025)
and the high dimensional fixed effects estimators developed by Correia (2017). Features include estimation of various types of average treatment effects, comparison of treatment effects
across subgroups, estimation of event study style estimates and related parallel trends assumption tests, and estimation of average and marginal treatment effects in the presence of a
continuous treatment (or treatment intensity) variable [experimental]. This suite also can produce event study plots, histograms of treatment effects by cohort, and plots of treatment
effects across the distribution of a specified continuous treatment variable. Options for inference include inference using standard analytic clustered standard errors, inference using
standard errors on average treatment effects that account for sampling variation in covariates (unconditional standard errors),
Ibragimov and Muller (2010) style cluster-robust inference, and randomization (permutation) inference for settings - most importantly few treated clusters - where cluster-asymptotic
inference is unreliable ({it:permute()}). A fast exploratory mode ({it:fastnose}) skips standard error estimation entirely for specification search.  {p_end}

{pstd}{ul:Shortcuts:} {help wooldid##description:Description}, {help wooldid##options:Options}, {help wooldid##output:Working with Output},
{help wooldid##runtime:Runtime and What Drives It}, {help wooldid##comparison:Comparison with Related Estimators}, {help wooldid##examples:Examples},
{help wooldid##bugs:Bugs and Assistance}, {help wooldid##citations:Citations}, {help wooldid##author:Author Info and How to Cite}, {help wooldid##disclaimer:Disclaimer} {p_end}

{pstd}{ul:Dependencies:} {cmd:wooldid} requires users have STATA 16 or higher installed, {help reghdfe} version 5.7.3, {help ftools} (for {it:reghdfe}),
and {help gtools}. Poisson estimation also requires the package {help ppmlhdfe}. {p_end}

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
{synopt : {it:t}}Time period identifier; must be numeric and on the same grid as {it:ttre} ({it:t} minus {it:ttre} must be an integer wherever {it:ttre} is nonmissing - see the {it:ttre} entry). Integer values are recommended, but the value itself only labels periods on the default estimation path: zero, negative, and non-consecutive values (gaps) are all fine. With {opt uncond:itionalse}, {opt customvce()}, or {opt inter:activecontrols()}, {it:t} enters the model as a native Stata factor variable and must then consist of non-negative integers no larger than 2,147,483,647 (Stata's factor-variable cap): recode large calendar codes such as %tc datetimes to small integers first.{p_end}
{synopt : {it:ttre}}Variable containing cohort-specific dates of treatment onset; missing indicates never-treated cohorts, cohorts treated before the sample begins are dropped. Binary designs require a comparison group: some units must be never-treated or first treated after the estimation window (an informative error is raised otherwise; see {it:oostreatedcontrols}). Store {it:ttre} (and {it:t}) as {cmd:double} when onset codes are large: 8-digit calendar dates exceed float precision, so a float-stored onset variable is rounded before {it:wooldid} sees it, silently merging adjacent cohorts. {it:ttre} must take a single value (or be missing on every row) within each group defined by {it:i} over the estimation sample: {it:wooldid} builds relative time row by row from {it:ttre}, so an {it:i} group carrying different onsets on different rows - including an onset recorded only from the onset period onward, with missing values pre-treatment - would be silently split across cohorts. {it:wooldid} verifies this and exits with an informative error if it is violated (bypass with {it:safety(off)}). Onset times must also sit on the same grid as {it:t} - {it:t} minus {it:ttre} must be an integer wherever {it:ttre} is nonmissing (integer spacing is not required: biennial or five-year panels are fine). Fractional offsets previously survived to the estimation stage, where they either failed with opaque errors or silently mislabeled the onset period as pre-treatment; {it:wooldid} now rejects them at parse time (bypass with {it:safety(off)}). If every offending onset falls outside the usable estimation window (pre-sample, or beyond the last observed period), the rejection softens to a Caution, since such onsets are dropped or treated as never-treated downstream anyway. {p_end}

{phang} {ul:Options by Purpose:} {p_end}
{synoptset 1 tabbed }{...}
{synopt : {help wooldid##specifyestimand:Specifying the Estimand(s):}} {opt att} {opt att_it} {opt att_i} {opt att_itime} {opt semi:elasticity} {p_end}

{synopt : {help wooldid##standarderrors:Standard Errors:}} {opt clus:ter(variable(s)[numerical])}  {opt rob:ust} {opt uncond:itionalse} {opt customvce(customstring)} {opt imp:values} {opt fast:nose} {opt fastcollapse} {opt fullprecisionslow} {opt perm:ute(method)} {opt permreps(int)} {opt permex:act} {opt permstr:ata(variable)} {opt permpr:obs(spec)} {opt permts:tat} {opt permtc:i} {opt boot:strap(mode)} {opt reps(int)} {opt ties} {opt level(#)} {opt dots(int)} {opt rseed(#)} {opt bootalp:ha(spec)} {opt bootd:f(spec)} {p_end}

{synopt : {help wooldid##supporting:Supporting Features:}} {opt hist:ogramcohorteffects} {opt save:plots(filename_stem)} {opt replace} {opt sum:marystats}
{opt pois:son} {opt poisexp:results}  {p_end}

{synopt : {help wooldid##eventstudy:Event Study Estimation:}}  {opt esf:ixedbaseperiod} {opt espre:length(int>=0)} {opt espost:length(int>=0)} {opt esr:elativeto(-int)}
 {opt joint:tests}   {opt adv:ersarial} {opt altmethod_adversarial} {opt makep:lots} {opt oos:treatedcontrols(int>=0)} {opt antic:ipation(int>=0)} {opt pretrendpow:er(numlist)}  {p_end}

{synopt : {help wooldid##interactions:Subgroup Effects:}} {opt sub:groups(variable[>=0 int])} {opt suppress:sgprimaryeffects}  {p_end}

{synopt : {help wooldid##continuous:Continuous Treatments:}} {opt contreat(variable[numerical])} {opt contreatelas:ticitytype(dydx, eydx, dyex, or eyex)}  {opt contreatw:ithin}
{opt contreatp:oly(+int)} {opt lattice(variable[+int] OR latticetype latticenum)} {opt latticeig:noreweights} {opt makec:fxplots} {opt cfxlattice(variable[+int])}
{opt cfxplott:ypes(post, pre, and/or integers specifying relative time periods)} {opt makecfxplotsbysg}  {opt cfxplotnose}  {opt contreatc:ontrols(conttreatcontroltype)}  {p_end}
{phang3} {it:-latticetype}: {it:all}, {it:byTime}, {it:byTreated}, {it:byTtre}, {it:byCohort}, {it:byTreatedByTime}, {it:byTtreByTime}, {it:byCohortByTime}, {it:baseTreated}, {it:baseTreatedRef},
{it:baseTreatedPost}, or {it:baseTreatedPre}  {p_end}
{phang3} {it:-latticenum}: integer > 1  {p_end}
{phang3} {it:-conttreatcontroltype}: {it:untreatedZero}, or one or more of {it:basic}, {it:byTreated}, {it:byTtre}, {it:byCohort}, {it:byTime}  {p_end}

{synopt : {help wooldid##controls:Controls and Fixed Effects:}} {opt control:s(fv varlist)} {opt fe(fv varlist)} {opt timetrend:s} {opt coarse:cohortcontrols(fv-required varlist)}
{opt inter:activecontrols(varlist)} {p_end}

{synopt : {help wooldid##technical:Technical Features:}}
{opt update} {opt customw:eightmultiplier(variable[>=0 numerical])}  {opt clean:matrices} {opt regtol(+int)} {opt ccno:absorb} {opt ccc_absorb(varlist[numerical])}
{opt ver:bose} {opt gte:ffects} {opt emptycellsoverride} {opt safety(on|off)}   {opt oldsyntax}  {p_end}



{marker description}{...}
{title:Description of Estimation Approach}

{pstd}
{bf:wooldid} estimates various types of difference-in-differences treatment effects in the vein of Wooldridge (2025), with extensions to handle continuous treatments and estimation of effects by
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
For more on the different estimand-types available, see {help wooldid##description:estimands}. {p_end}
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
Note that this approach is a bit speculative and is not covered in Wooldridge (2025), though Wooldridge has
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
and thus that these cohorts should be allowed to contribute to pre-treatment relative time period-specific estimates in the event study context. Relatedly, when units may respond
in anticipation of a known upcoming treatment, the option {it:anticipation(#)} treats the effective onset as {it:ttre} − #, so anticipation windows are not misclassified as clean
pre-treatment periods.
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
As such, this option should be slower to use on average. This option requires specifying either robust or clustered standard errors (a robust or clustered vce supplied via {it:customvce()} also satisfies this requirement, as do pweights, which imply robust standard errors). {p_end}

{phang}{opt customvce()}: Forces {it:wooldid} to switch to using either the {help regress} or {help poisson} command for underlying regression estimation,
and then passes the user specified text without modification into the {it:vce()} option for these commands. This allows specifying a broader set of
standard error calculation methods than allowed by the default high dimensional fixed effects estimation packages used by {it:wooldid}. Not compatible with {it:robust} or
{it:cluster()}: specify the clustering or robustness desired inside the {it:vce()} text itself. {p_end}

{phang}{opt imp:values}: Implements Ibragimov and Muller (2010) style cluster-robust inference for all estimates. The idea behind the Ibragimov and Muller approach is to compute
cluster-specific treatment effect estimates, present their average as one's final treatment effect estimate, and conduct inference by t-testing the cluster-specific estimates against 0.
Using IM-style inference requires specifying a single cluster variable and is not compatible with the {it:unconditionalse} option. IM-style inference implicitly changes {it:wooldid}'s
target estimand. For example, if clustering at the cohort level, IM-style inference will cause {it:att} to yield the same effect estimates as {it:att_i}. Note that {it:wooldid} does not adjust the
weights used for computing different estimands when IM-style inference is requested: in general, IM-style inference will yield incorrect results for non-att estimands if clusters cross-cut cohorts
or (worse) cohort-time period cells. IM-style inference is not currently allowed with the use of subgroups, the option {it:jointtests}, and the option {it:makecfxplots}. {p_end}

{phang}{opt fast:nose}: Fast exploratory mode. Every {it:margins} call runs with {it:nose}, so no standard errors, test statistics, p-values, or confidence intervals are computed
for any estimand: only point estimates are produced, the corresponding cells of the results tables display as missing (.), and the output carries the line "SEs and related materials
are not estimated in fast mode." Intended for specification search on large panels, where the delta-method variance calculation inside {it:margins} dominates run time (measured
speedups of roughly 2-4x on typical specifications, and far larger on poisson event-study runs, where fast mode completes in seconds on configurations that otherwise run for tens of
minutes); re-run the chosen specification without {it:fastnose} before reporting results. Under {it:fastnose}: {it:e(V)} is not posted ({it:e(properties)} is {it:b}), so {it:test},
{it:lincom}, {it:esttab} standard errors, and anything else needing a variance-covariance matrix will refuse; {it:e(b_es_)}, {it:e(V_es_)}, and {it:e(numpre_)} are not returned, so
{it:honestdid} and Roth-style diagnostics have no inputs; no joint tests of any kind are computed (including the subgroup homogeneity tests); the continuous-treatment-effects plot
matrices carry only the point-estimate column, exactly as under {it:cfxplotnose}; and event-study plots are drawn without confidence bands. Not compatible with
{it:unconditionalse}, {it:adversarial}, {it:jointtests}, {it:pretrendpower()}, {it:impvalues}, {it:permtstat}, or {it:permtci} (each is refused with an error). Point estimates can differ from a full run in
roughly the last digit of precision (absolute differences on the order of 1e-15; relative differences are typically similar but can reach roughly 1e-13 on subgroup pairwise
contrasts, whose differencing amplifies rounding, up to roughly 1e-14 under {it:timetrends}, and larger still in relative terms where an estimate is itself near zero), because
{it:margins} follows a different computational path when the variance is skipped. {p_end}

{phang}{opt fastcollapse}: The cell collapse. {it:wooldid} speeds up its {it:margins} calls by evaluating them on one row per cohort-by-period-by-treatment cell instead of on every
observation. {it:wooldid}'s estimands are weighted averages of cell-level treatment effects, and within the conditions listed below the treatment derivative is constant inside such
a cell, so the collapsed row - carrying the cell's summed weight - represents the whole cell exactly. The regression itself is untouched and still runs on the full data, so the
coefficients, the variance matrix and the clustering are unaffected; only the averaging step that {it:margins} performs is handed a smaller dataset. The full data is restored
before any plotting or display work. {bf:Continuous treatments are covered too}: with {it:contreat()} the collapsed row additionally carries the cell's weighted mean of the
treatment variable, which reproduces every reported effect exactly because each of them is linear in the dose inside a cell - the treatment-dummy derivative is
{it:delta + beta*dose}, the dose derivative ({it:ame}) is constant, and {it:contreatelasticitytype(dyex)} is proportional to the dose. The forms that are {ul:not} linear in the
dose are listed below and turn the collapse off. {bf:The collapse is on by default} whenever the request permits it, so this option is never needed to obtain it: specifying it asks {it:wooldid}
to {bf:insist} on the collapse, and turns the silent skip described under {it:What turns the collapse off} into an error. Use {opt fullprecisionslow} to turn the collapse off. The
collapse composes with {it:fastnose}, is applied inside {cmd:permute()} and {cmd:bootstrap()} replications, and is recorded in the hidden macro {it:e(wdid_fcol)} (1 when it was
applied, 0 when it was not); {it:verbose} reports the decision as the run proceeds. {p_end}

{phang2}{ul:When it pays.} The saving grows with the number of observations per cell and with the size of the regression, and is unaffected by the choice of standard errors. On a
30,000-observation administrative panel with roughly 500 rows per cell it removed 35-60% of total run time; on a 475,000-row version, 56%. It pays for itself inside a single
{it:margins} call at that scale - a bare {it:att} is already enough, with no event-study horizons requested: a vanilla {it:att} run with {it:contreat()} was 1.24 times faster on a
4,000-observation sample, 1.25 times on a 29,000-observation one and 1.33 times on a 145,000-observation one, so building the cell means does not eat the saving. Adding
event-study horizons makes it a better deal rather than a worse one, and so does asking for more estimands: the same 29,000-observation run went to 2.0 times faster with all four. It gains little when there
are only a handful of rows per cell (roughly four or fewer), where the collapse can cost slightly more than it saves, and it helps least on designs with very many cohorts and
periods, where the absorbing regression rather than {it:margins} dominates the run. See {help wooldid##runtime:Runtime and what drives it}. {p_end}

{phang2}{ul:Results agree to machine precision but are not bit-identical.} Collapsing changes the order in which {it:margins} accumulates its averages, and the controls' contribution
cancels algebraically between the treated and untreated predictions without cancelling bit for bit. If you diff output across runs, expect low-order digits to differ between a
collapsed and an uncollapsed run; if you need byte-identical reproductions of a result produced by an earlier {it:wooldid} version, which had no collapse, run with
{opt fullprecisionslow}. {p_end}

{phang2}{ul:How large are the differences?} The difference is bounded, in {bf:absolute} terms, by about one unit in the last place of the linear prediction - the size of the
prediction divided by 2 to the 52nd power - and every difference measured was a small integer multiple of that one quantum. On a panel whose predictions are of order 10, that is about 2e-15. Note this scales with the prediction, not with the
estimate: an outcome with a very large mean raises {it:|xb|} and the absolute differences with it. {p_end}

{phang2}Because the bound is absolute, what it means in {bf:relative} terms depends on how big the estimate is. An ordinary coefficient - say an ATT of order 0.5 - agrees to roughly
15 significant digits (relative differences around 3e-15). An estimate sitting near zero, such as a small pre-treatment coefficient of order 0.01, can show a relative difference
around 1e-13, purely because the denominator is small rather than because the error is larger. Standard errors and variances do not amplify at all: they were measured at about
1e-16 relative regardless of the scale of the predictions. P-values in a far tail move more in relative terms - a p-value of 2e-07 moved by about 3e-11 relative - which is a
property of tail probabilities, not of this option. {p_end}

{phang2}{ul:Continuous treatments carry a second rounding step, and a looser bound.} The figures above are for the binary case, where the collapsed row only re-orders the
summation. With {it:contreat()} the row also carries the cell's weighted mean of the dose, which is itself a rounded quantity, so the bound is set by that mean rather than by one
unit in the last place of the prediction. Measured over 305 paired continuous-treatment runs comparing every stored result cell by cell: no cell moved by more than 1e-9 in
relative terms at any magnitude down to 1e-12; the largest relative move on a quantity of ordinary size was 1.9e-10, on a p-value of 1.2e-8; the largest absolute move anywhere
was 4.8e-10, in the tenth significant digit of a t statistic of 39.2; coefficients topped out at 1.1e-10 relative, on a pre-treatment coefficient of 3.2e-05, and variances at
2.6e-11. Two decimal orders of that headroom are spent on near-zero cells: an ordinary coefficient still agrees to about 15 significant digits. As in the binary case, no digit
{it:wooldid} displays is affected. {p_end}

{phang2}{ul:When should you care?} For almost all work, not at all: no digit {it:wooldid} displays is affected, and the differences are several orders of magnitude below any
sampling uncertainty you are reporting. Two situations warrant attention. First, if you compare output byte-for-byte across runs, or archive results for exact reproduction - then
be deliberate about whether the collapse is on or off ({opt fullprecisionslow}) in every run you compare, and bear in mind that a request outside the eligible set is run
uncollapsed without notice, and that a {it:subgroups()} run can switch to the full data part-way through (see below; {it:e(wdid_fcguard)} records it). Second, if your analysis differences two nearly equal estimates, or otherwise divides by a quantity that is near zero, since that amplifies the relative
size of any rounding difference. {p_end}

{phang2}A useful calibration for the second case: merely permuting the order of the rows in your dataset - leaving the data otherwise untouched, with the collapse off -
moves the same near-zero coefficient about 19 times {bf:more} than the collapse does. So if a conclusion is sensitive at this scale, it is already sensitive to the order your
observations happen to be stored in, which is worth knowing whether or not the collapse is in use. {p_end}

{phang2}{ul:What turns the collapse off: the complete list.} There are three ways the collapse does not run, and this paragraph is the single place that enumerates them.
{bf:(1) You turn it off}, with {opt fullprecisionslow}. {bf:(2) The request is ineligible}, in which case the collapse is {bf:skipped} and the run proceeds on the full data
exactly as it would have before the collapse existed - silently, unless {it:verbose} is on, which names the option responsible. If {it:fastcollapse} was specified explicitly the
same combination is instead refused with an error, and {it:safety(off)} does {ul:not} override that refusal, because each entry below is a demonstrated wrong answer rather than
an accepted risk. {bf:(3) A {it:subgroups()} run sets it aside part-way through}, which is a screen rather than an eligibility rule and is described in the next paragraph. {p_end}

{phang2}The ineligible list, grouped by the condition each one breaks. {ul:The treatment effect is no longer constant inside a cell:} {it:interactivecontrols()}, which interacts
the treatment dummy with a variable that varies within the cell; and, under {it:contreat()}, the treatment variable {ul:also} appearing inside {it:controls()},
{it:coarsecohortcontrols()}, {it:fe()} or {it:ccc_absorb()}, which makes the dose derivative vary within the cell on a variable the collapse does not average. (A control term
built from a {ul:copy} of the treatment variable under a different name is not detected and stays eligible; it was also measured exact, {it:margins} not differentiating through
the copy.) {ul:The reported quantity is not linear in the coefficients, or not linear in the dose:} {it:semielasticity} and the continuous-treatment elasticity forms
{it:contreatelasticitytype(eyex)} and {it:(eydx)}, which divide by the prediction and so re-introduce the within-cell variation that a plain marginal effect cancels;
{it:contreatpoly()} above 1, which makes the effect a polynomial in the dose that the cell mean does not reproduce; {it:poisexpresults}, whose exponentiated reporting scale is
not linear in the coefficients; and {it:poisson} together with {it:contreatelasticitytype(dyex)}, which {it:margins} reports through an {it:expression()} it differentiates
numerically, so that replacing the dose by its cell mean changes the round-off of the derivative rather than cancelling out of it.
{ul:The read splits a cell on something the cell key does not contain:}
{it:lattice()}, {it:makecfxplots} and {it:cfxlattice()}, which bin the dose; and {it:impvalues}, which computes effects over the clusters.
{ul:{it:margins} objects to the collapsed row count itself:} {it:unconditionalse}, whose variance is built from observation-level influence functions and whose row count
{it:margins} checks against {it:e(N)}. {ul:Unverified rather than known-broken:} {it:customvce()}, refused pending measurement. {p_end}

{phang2}Everything else is eligible, including several combinations that are easy to assume are not. {it:contreat()} is eligible on the default {it:dydx} scale and under
{it:contreatelasticitytype(dyex)}, with every flavour of {it:contreatcontrols()} and with {it:ccnoabsorb}. Poisson estimation on the default (log) reporting scale is eligible,
including with {it:contreat()}. So are {it:fastnose}, {it:jointtests}, {it:adversarial} and {it:altmethod_adversarial}, {it:esfixedbaseperiod} and {it:esrelativeto()},
{it:pretrendpower()}, {it:histogramcohorteffects}, {it:summarystats}, {it:gteffects} (which {it:wooldid} refuses alongside {it:contreat()} for a separate reason), {it:makeplots}
and {it:saveplots()}, and the replications of {cmd:permute()} and {cmd:bootstrap()}, which are evaluated collapsed like any other run. {p_end}

{phang2}{ul:{bf:subgroups()} is eligible, but screened, for a subtler reason worth understanding.} When {it:margins} judges a subgroup-by-horizon effect not estimable it does not
error - it computes the number anyway and flags it only by marking the name it returns, and that flag rests on a rank decision made at the last digit. wooldid's subgroup code is
the one place that reads the flag, and it does so three separate ways: to drop the affected horizon with a caution, to skip the adversarial pre-trend test, and inside the subgroup
joint tests. So it is the one place where the collapse's last-digit change in summation order can change {bf:what} is reported rather than its low-order digits: a horizon that
should be dropped can come back as a meaningless near-zero, or a reported joint test can turn to missing. This can only happen where the variance matrix has a direction that is
singular in truth - a subgroup effect, or a contrast between subgroups, whose true variance is zero, as when a fixed effect absorbs every cell ({it:fe()} on a cohort-by-period
variable), the fit is exact, or the clustered variance is rank-deficient - so that the number {it:margins} computes there is cancellation residue rather than an estimate. wooldid
therefore screens each subgroup {it:margins} call after the fact, inspecting the returned variance matrix for a not-estimable marker, a subgroup effect or a between-subgroup
contrast whose variance is below 1e-9 of the largest, or a contrast direction singular at that tolerance (the smallest genuine subgroup contrast measured across the test inventory
is 3.5e-5 and the residue is around 1e-16, so the band has ample room on both sides). The screen is a few matrix operations per call. If it finds nothing, the collapsed result
stands. If it finds a singular direction, then under the default safety setting wooldid restores the full data at that point and evaluates that call and every remaining subgroup
call on it, so everything downstream of the flag is exactly what {opt fullprecisionslow} would have produced (the main results, and any subgroup horizons already read before the
flag, were computed collapsed and agree with their full-data values to the last digit as described above); {it:verbose} reports when this happens. Under {it:safety(off)} wooldid
instead keeps the collapsed rows, as requested, and prints a caution naming what it found, since in that situation which subgroup horizons and joint tests are reported can depend on
the summation order. The outcome is recorded in the hidden macro {it:e(wdid_fcguard)}: empty when nothing was flagged, {it:fallback}, or {it:caution}. An explicit
{it:fastcollapse} is a request the screen may set aside in the same way; and a fallback carries into {it:bootstrap()} replications, which are otherwise evaluated collapsed. {p_end}

{phang2}{ul:When does the screen fire?} In practice when a contrast between two subgroups has {bf:exactly} zero variance in the full data too, so the "significance" of that
contrast rests on rounding residue whichever way it is computed. The common route is clustering on a variable that nests the cohort variable - {it:cluster(state)} where the cohorts
are state-by-group cells, or clustering on the cohort itself - so that subgroups of identical cohort-and-period composition are perfectly correlated; a handful of coarse clusters,
a subgroup that is an exact shifted copy of another, and a subgroup nested inside a unit reach the same place. wooldid's own identical-standard-error, all-zero-variance or
not-estimable caution usually accompanies such a run but does not always (a robust variance on a saturated design, or clustering on the period, can be singular there without it),
so the screen adds coverage. Every configuration on which the collapse had been measured to change reported subgroup output is caught by it, and across an adversarial probe of
over 400 further subgroup configurations - noise shrunk toward the tolerance, weights driven toward zero, near-balanced compositions, rescaled and heteroskedastic outcomes,
estimability edges, replications, every downstream reporting path - no result read while the collapse was engaged differed from its full-data value by more than amplified
rounding: at most 1e-7 relative, on between-subgroup contrasts whose variance sat within a factor of ten of the tolerance, never a point estimate, standard error or reported
set. Healthy
subgroup runs pass untouched: realistic heteroskedasticity leaves the screened variance ratios between 0.1 and 1, seven to eight orders of magnitude from the tolerance. The
one way a healthy variance matrix reaches it is outcome scales differing between subgroups by more than about a factor of 1e4 (a variance ratio of 1e-9): the cost is then the
speedup under the default safety, and under {it:safety(off)} a caution that overstates the situation. Non-subgroup runs are not exposed: no code outside the subgroup path reads a
{it:margins} result name, and 123 non-subgroup comparisons on aggregating data showed no structural difference at all. {p_end}

{phang2}{ul:Memory.} The full data is held aside while the collapsed window runs, so peak memory includes an extra copy of your dataset. This is immaterial for typical panels but
worth knowing if you are already near your memory limit, particularly if you have issued your own {help preserve} before calling {it:wooldid}; {opt fullprecisionslow} removes
the extra copy. {p_end}

{phang}{opt fullprecisionslow}: Turn the default cell collapse off, so that every {it:margins} call is evaluated on the full estimation sample as it was before the collapse
existed. Results are then bit-identical to those of earlier {it:wooldid} versions given the same data and options, at the cost of the run time the collapse would have saved
(see {opt fastcollapse} for the magnitudes on both sides). Use it when you need byte-for-byte reproduction of an archived result, when you are comparing output across runs at
the last digit, or when memory is tight and the extra copy of the data held during the collapsed window is unwelcome. Forwarded to {cmd:permute()} and {cmd:bootstrap()}
replications. Cannot be combined with {it:fastcollapse} (the pair is refused as contradictory). {p_end}

{phang}{opt permute(method)}: Randomization (permutation) inference in the Conley-Taber spirit, for settings - most importantly few treated clusters - where cluster-asymptotic standard errors are unreliable. The observed
estimation runs exactly as usual (full inference, display, plots); afterwards, wooldid re-randomizes which clusters carry which treatment-onset dates, re-estimates the requested estimands under each counterfactual
assignment (quiet fast-mode runs, since only point estimates are needed), and compares the observed estimate to the resulting null distribution. Requires {cmd:cluster()} with a single variable - the randomization
unit (the randomization design consists of the clusters having at least one row where the {it:if} expression evaluates to exactly 1 with a nonmissing cluster value; rows {it:wooldid}
later removes by its deterministic markout, e.g. missing outcomes or zero weights, still count toward the design and its possible-assignment count) - and every cluster must sit inside exactly one treatment-onset group (onset constant within cluster; never-treated is its own group; never-treated onsets must be coded as plain missing {cmd:.} - tagged missing codes {cmd:.a}-{cmd:.z} are refused, recode them to {cmd:.} first). Two methods: {cmd:permute(shuffle)} permutes the observed onset labels
across clusters, holding the number of clusters at each onset fixed ({cmd:permstrata(varname)} confines the shuffling to within strata, so each onset group keeps its composition with respect to the strata;
clusters must be strict subsets of strata); {cmd:permute(probabilities)} draws each cluster's onset from user-supplied probabilities via {cmd:permprobs()}. Specify {cmd:permreps(#)} (draws with replacement) or,
for shuffle only, {cmd:permexact} (enumerate every distinct assignment; a caution prints above 1,000 assignments and enumeration is refused above 10,000). The number of possible assignments is always computed and
returned, with a caution when {cmd:permreps()} exceeds 10 percent of it; a design admitting only one possible assignment (e.g. every stratum onset-homogeneous, or every {cmd:permprobs()} row deterministic) is refused - there is nothing to randomize. Outputs per estimand: a two-sided permutation p-value (sampled draws use the add-one Monte-Carlo form; exact enumerations report the plain
share, with the observed assignment a member - note that in symmetric designs the mirror-image assignment ties the observed estimate in absolute value, so exact p-values have a floor of 2/M); the 5th-95th
percentile band of the permuted-effect distribution (a null band under no effect, not a confidence interval for the parameter); and the standard deviation of that distribution (labelled SD; exact enumerations report the population SD - denominator M - of the complete enumerated distribution, while sampled draws report the n-1 sample SD). The raw draws are in
{it:e(permdist)}. {cmd:permreps()} must be at least 20. Replications that fail to estimate are excluded and counted ({it:e(perm_failed)}), with a caution above 10 percent and a refusal when fewer
than max(20, half) succeed (exact enumerations are exempt from that floor, but any exact-mode failure draws a caution: the enumeration is then incomplete and the reported exact
p-value is a share over completed assignments only, so the plain-share / observed-is-a-member / 2/M-floor reading above no longer applies). A degenerate null distribution (every completed draw equal) draws a caution. A missing or unreadable observed comparison statistic (e.g. a degenerate observed standard error under
{it:permtstat}) is refused loudly with r(459) rather than silently zeroing the tail counts. Draws consume the session
random-number stream: {cmd:set seed} beforehand for reproducibility. Under {cmd:subgroups()} the permutation statistics cover the pooled estimands only (rep runs skip the subgroup machinery). With
{cmd:contreat()} the dose stays with its unit and only onset timing is re-randomized: the null being tested is that timing is as-if random given dose. Permutation statistics cover the att-family estimands
(the continuous-treatment ame quantities are not permuted in this version). Note also that reassigned onsets can change which periods are estimable in a given draw: each draw is the estimate wooldid
would compute under that counterfactual assignment, estimability screens included. Under {it:permtstat}, {it:e(perm_lb/ub/sd_*)} and the displayed band/SD are on the t-ratio scale. Not compatible with {it:impvalues},
{it:unconditionalse}, {it:contreatwithin}, {it:contreatcontrols(untreatedZero)} (replications reassign onset while the dose path stays fixed, contradicting untreatedZero's dose-zero-outside-treatment premise), or {it:oldsyntax}. {p_end}

{phang}{opt permprobs(spec)}: The assignment-probability design for {cmd:permute(probabilities)}, structured as {cmd:(onsetvar level level ...)} followed by one {cmd:(clustervalue p p ...)} group per cluster.
The first group names the actual timing variable and lists every in-sample onset level (never-treated is the level {cmd:.}) in the column order the probabilities use; each later group gives one cluster's
probabilities over those levels (nonnegative, summing to 1; rows within the acceptance tolerance of 1e-6 are rescaled to sum exactly 1 before drawing). Every in-sample cluster value must appear exactly once. {cmd:permreps()} is required. The possible-assignment count reported is the product of the
per-cluster support sizes; if the observed assignment has zero probability under the design, a caution notes that the p-value is a tail probability under that design distribution rather than a randomization
test. {p_end}

{phang}{opt permtstat}: Randomize on t-ratios instead of raw estimates. Rep runs then go through the full inference path (slower: fast mode cannot be exploited, and requesting {it:fastnose} alongside is refused);
draws with undefined t-ratios count as failed. Because standard errors vary across assignments, the exact-test floor is 1/M rather than 2/M. {p_end}

{phang}{opt permtci}: Additionally report an inverted-t 95% confidence interval: the observed estimate plus or minus the 95th percentile of the permuted absolute t-ratios times the observed standard error.
This is a fixed-zero-null pivotality approximation - it relies on the permuted t-ratio distribution approximating the t-ratio's distribution under shifted nulls - and is not a
sharp-null randomization-test inversion. Implies {cmd:permtstat}. Unavailable (with a caution) when the observed run has no usable standard error, and refused loudly when the
observed comparison statistic itself is missing. {p_end}

{phang}{opt boot:strap(mode)}: Bootstrap standard errors and inference. When specified, the reported standard errors, z statistics, p-values, confidence intervals, joint tests and subgroup contrasts are bootstrap quantities, replacing the analytic ones; point estimates are unchanged, and the analytic covariance is preserved in the hidden {it:e(V_analytic)} (and {it:e(V_es_analytic_*)}). Each replication re-estimates the full model, so runtime scales with {cmd:reps()}. Two modes. {cmd:bootstrap(standard)} is the cluster bootstrap: clusters - the {cmd:cluster()} variable, or observations when no {cmd:cluster()} is given - are resampled with replacement. For the observation-weighted estimands ({it:att}, {it:ame}) replications are computed by cluster frequency reweighting, which is the same statistic as physical resampling and matches Stata's own {helpb bootstrap} prefix replication for replication under a common seed. When a cell-weighted estimand ({it:att_it}, {it:att_i}, {it:att_itime}) is requested, replications instead physically resample and re-identify duplicated clusters, so a twice-drawn cluster contributes its treatment cells twice; that leg requires {cmd:cluster()} to nest {it:i} and is refused otherwise. Re-identification roughly doubles the cost of each replication, because the resampled data carries distinct cohort identities for duplicated clusters and so presents a larger cell count to the regression. That cost is inherent to making those estimands' bootstrap correct rather than an inefficiency, but it is worth budgeting for when choosing {cmd:reps()}. {cmd:bootstrap(bayes)} is the Bayesian bootstrap: rather than resampling, each replication draws a total Dirichlet weight per cluster and splits it within the cluster in proportion to any existing weights, so no cluster ever leaves the sample and replication estimability is stable. Because reweighting cannot resample cells, {cmd:bootstrap(bayes)} combined with a cell-weighted estimand omits the cell-sampling variance component and can understate those standard errors when treatment effects vary across clusters: it is refused under default safety and proceeds with a caution under {cmd:safety(off)}, while {cmd:bootstrap(standard)} carries no such limitation. In simulations on staggered designs with 24 clusters both modes were close to calibrated for {it:att}, with the resampling mode nearer nominal than the reweighting mode. Requires {cmd:reps()}. Incompatible with {cmd:permute()}, {cmd:impvalues}, {cmd:fastnose}, {cmd:unconditionalse}, {cmd:customvce()}, {cmd:pretrendpower()}, {cmd:suppresssgprimaryeffects}, {cmd:oldsyntax}, and multi-way {cmd:cluster()}. With {cmd:jointtests} the joint tests become chi-squared Wald tests on the bootstrap covariance, and that combination is exempt from the analytic joint-test screen described under {cmd:jointtests} - bootstrap inference is what that screen recommends. Plots, saved plot data, {it:gteffects} cell standard errors, histogram output and adversarial line-test p-values are computed from the analytic covariance before the replications run and remain analytic; a caution says so. Stata's own {helpb bootstrap} prefix can also be used on an unweighted wooldid call ({cmd:bootstrap _b, reps(#) seed(#) cluster(c) idcluster(nc): wooldid y i t ttre, ... fastnose}), but it rejects weighted commands, cannot deliver the Bayesian flavor, and discards replications whose coefficient vector lost a horizon - the built-in {cmd:bootstrap()} option avoids all three. {p_end}

{phang}{opt reps(#)}: Number of bootstrap replications; required with {cmd:bootstrap()}. No minimum is imposed, so small values can be used to time a specification - but bootstrap inference needs many replications: hundreds for standard errors, and a thousand or more before the percentile and bias-corrected intervals are worth reading. Replications that fail to estimate - for example a resample retaining no comparison group - are excluded and counted in {it:e(boot_failed)}, the convention Stata's own {helpb bootstrap} follows. {p_end}

{phang}{opt level(#)}: Confidence level for the bootstrap confidence intervals, between 10 and 99.99; the default is 95. Applies to bootstrap output only. {p_end}

{phang}{opt rseed(#)}: Seed for the bootstrap draws, passed to {cmd:set seed}; as with the built-in prefixes this advances the session's random-number state. Without it the draws consume the current state and the run is not reproducible. {p_end}

{phang}{opt dots(#)}: Show a progress dot every {it:#} replications, with {cmd:dots(0)} suppressing the display; the default is 1. {p_end}

{phang}{opt ties}: Adjust the bias-corrected confidence intervals' median-bias constant for ties between the replications and the observed statistic, exactly as the built-in {helpb bootstrap}'s {cmd:ties} option does. It has no effect on continuous statistics, which do not tie. {p_end}

{phang}{opt bootalp:ha(spec)}: The Dirichlet convention for {cmd:bootstrap(bayes)}. {cmd:bootalpha(equal)}, the default, draws Dirichlet(1,...,1) across clusters - the Rubin convention. {cmd:bootalpha(size)} instead gives cluster {it:c} the Gamma shape {it:n_c}, its observation count, reproducing what {helpb bayesboot} does when passed {cmd:id(}{it:clustervar}{cmd:)}; it is available for exact comparison against that command, but it is not a sound default for panel data, because tying the Dirichlet's concentration to the observation count rather than the cluster count shrinks the cluster-weight dispersion by roughly the square root of the mean cluster size and the resulting standard errors are badly understated. {p_end}

{phang}{opt bootd:f(spec)}: The reference distribution for the reported bootstrap statistics. {cmd:bootdf(normal)} uses standard normal critical values and reports z statistics - the convention Stata's own {helpb bootstrap} prefix follows. {cmd:bootdf(clusters)} instead uses the t distribution with {it:G}-1 degrees of freedom, where {it:G} is the number of resampling units ({cmd:cluster()} groups, or observations when no {cmd:cluster()} is given), and a positive integer sets the degrees of freedom directly. {bf:The default depends on whether you clustered}: {cmd:bootdf(clusters)} when {cmd:cluster()} was given, and {cmd:bootdf(normal)} when it was not. The bootstrap standard error is estimated rather than known, so where a cluster count exists the t({it:G}-1) correction is the appropriate one; with no {cmd:cluster()} the resampling unit is the observation, and a {it:G}-1 rule would be {it:N}-1 and indistinguishable from the normal. An explicit {cmd:bootdf()} always overrides the default. This changes only the reference distribution, never the bootstrap standard errors themselves: the p-values, the normal-theory confidence intervals and the joint tests widen, with the joint tests reported as F({it:k},{it:df}) rather than chi-squared({it:k}); the percentile and bias-corrected intervals are quantile-based and are unaffected. The t option is worth considering when the number of clusters is modest: simulations on 24-cluster designs found coverage running below nominal even where the bootstrap standard error itself was calibrated, because the standard error is estimated rather than known, and a t({it:G}-1) critical value is the standard correction. Under every setting {it:e(df_r)} is aligned with the reference distribution actually reported - the same degrees of freedom under {cmd:bootdf(clusters)} or an explicit integer, and missing under {cmd:bootdf(normal)} so that the fallback is chi-squared and z - so a {cmd:test} or {cmd:lincom} run by hand afterwards reproduces the reported p-values rather than diverging from them. The analytic value is preserved as hidden {it:e(df_r_analytic)}. {p_end}

{pstd} {ul: Discussion of Unconditional Standard Errors:}
The overall treatment effect estimates presented by {it:wooldid} depend on not just the coefficients from the underlying regression, but also on the observed values of any included covariates
or continuous treatment variables. The standard errors {it:wooldid} presents for these final estimates are produced during the margins calculation stage
using the delta method, taking the observed values of the relevant variables as given (i.e., conditional on those values). This contrasts with approaches that aim to account for sampling variation
in those observed values (i.e., unconditional standard errors; see the discussion of {it:vce(unconditional)} in {help margins}). Whether one should wish to take the observed values of one's
covariates as fixed likely varies from context to context. {p_end}

{pstd} Does using unconditional standard errors make a large difference? Wooldridge (2025) expresses the view that conditional and unconditional standard errors often should be similar.
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
generated. Everything is saved to Stata's working directory; the argument must be a plain name stem (paths are rejected by the option's syntax).{p_end}

{phang}{opt replace}: When saving materials to the file names specified in {it:saveplots()}, overwrite any existing files with the specified names. {p_end}

{phang}{opt sum:marystats}: Compute and present summary statistics of the outcome variable y and, if present, of the continuous treatment variable contreat. Also displays the event-study composition matrix's coverage scalars (the reference-normalization pool's N, treated-unit and onset-cohort counts, and att-style weight mass, plus the same for the entire treated estimation sample; see {it:e(escomp_*)} under Stored Results).
Summary statistics are computed overall, within the treatment group by time period (post-treatment vs pre-treatment), and within the control group. If a
fixed reference period is not being used and pre-treatment relative time period effects are requested, statistics on the size of the reference group relative to the rest of the
pre-treatment period will also be shown. {p_end}

{phang}{opt pois:son}: Implement the specified {it:wooldid} command using poisson regression, achieved by switching from the regression command {it:reghdfe} to {it:ppmlhdfe}. Treatment effects
are presented in terms of differences in linear predictions (the poisson {it:predict(xb)} option). Poisson regression is especially useful in situations where one would like to use
log(y) as an outcome, but where y contains natural 0s, as well as when treatment and control groups are only on parallel trends in proportions (e.g., one group is always 2x the other).
Poisson estimation accepts pweights but not aweights (refused with an error). If {it:ppmlhdfe} fails on the initial specification, {it:wooldid} retries once with any time-trend,
coarse-cohort-control, and continuous-treatment terms moved out of the absorbed-effects list into the regression itself (a caution prints); results are equivalent up to numerical
implementation. A pre-existing variable named {it:_ppmlhdfe_d} is overwritten by {it:ppmlhdfe}'s own scaffolding on poisson runs and then removed by {it:wooldid}'s cleanup - treat
that name as reserved in your data. {p_end}

{phang}{opt poisexp:results}: When implementing poisson, present treatment effects in terms of differences in exponentiated linear predictions (the poisson "number of events" or {it:predict(n)}
option). Results of this sort will thus not have a log-type interpretation. Specifying {it:poisexpresults} without {it:poisson} is an error, since the option has no effect outside of poisson estimation. {p_end}


{marker eventstudy}{...}
{dlgtab:Event Study Estimation}

{phang}{opt esf:ixedbaseperiod}: Use a fixed reference period for estimation, rather than a pooled reference period.   {p_end}

{phang}{opt esr:elativeto(-integer)}: The relative time period used as the fixed reference period, if a fixed reference period is being used. The default is -1, the time period immediately
prior to treatment. Values above -1 (the onset period or later) are rejected with an error: the fixed reference period must lie strictly before onset. Without {it:esfixedbaseperiod}, this option has no effect and a Caution is printed. {p_end}

{phang}{opt espre:length(integer)}: The number of pre-treatment relative time period-specific effects to be estimated. If any pre-treatment effects are requested and
{it:esfixedbaseperiod} is specified, all relative time period effects between the {it:esrelativeto()} period and treatment onset will be presented.{p_end}

{phang}{opt espost:length(integer)}: The number of post-treatment relative time period-specific effects to be estimated. With {it:subgroups()} and {it:esprelength()} > 0, a value of
0 is raised to 1 (a caution prints): subgroup estimation requires at least one post-treatment horizon. {p_end}

{phang}{opt joint:tests}: Calculate joint statistical significance tests of all requested pre-treatment relative time period-specific effects. Caution: simulations suggest these
tests tend to be particularly sensitive to whether you have a sufficient number of clusters in your sample for standard error estimation. Two further reliability conditions apply.
First, {it:jointtests} is not compatible with wooldid's analytic clustered SEs when clustering is at the level of the cohort variable {it:i} or above: the joint test p-values are
severely oversized in that configuration (the cluster-robust covariance cannot capture the within-cluster cross-horizon variation the tests depend on, a problem that does not
diminish with more cohorts or larger cells), so {it:wooldid} refuses the combination with an error under default safety ({it:safety(off)} downgrades the error to a caution and
proceeds, reporting the unreliable p-values as-is); if joint-type inference is needed in such settings, use bootstrap-based clustered standard errors or randomization inference
(see {it:permute()}). Second, in the continuous treatment case the joint tests are only correctly sized when the i-by-t cells each contain many observations: the tests aggregate
within-cell dose-slope variances that carry a downward small-sample bias which compounds across horizons, so designs with only a handful of observations per cohort-period cell
yield oversized joint tests even under robust or finely clustered SEs (note that {it:lattice()} refines cells smaller and {it:contreatpoly} adds within-cell parameters, both of
which sharpen this concern). {p_end}

{phang}{opt adv:ersarial}: Tests the composite null hypothesis that all estimated relative time period-specific effects fall on a single line - some line, with slope and intercept
unknown - a null that is adversarial toward finding pre-trend violations in that the line is chosen to fit the estimates as well as possible. The p-value comes from a calibrated
minimum-distance test on the event-study coefficients and their full variance-covariance matrix (the same objects returned in {it:e(b_es_type)} and {it:e(V_es_type)}): the nuisance
line is estimated by generalized least squares under a generalized inverse of that covariance, and the minimized quadratic distance is referred to an F distribution whose numerator
degrees of freedom equal the residual rank (the rank of the event-study covariance minus the number of line parameters identified), with the regression's residual degrees of freedom
in the denominator (a chi-squared reference is used when no residual degrees of freedom are available). The p-value requires {it:esfixedbaseperiod}: under the pooled reference regime
no single omitted period exists, so a caution is printed and the fitted line is reported without a p-value. If the event-study covariance is degenerate (as is typical under
{it:timetrends}), a caution is printed and no p-value is computed. For display, wooldid also reports a line fit by the unweighted OLS regression of {it:coef_rt} on {it:time_rt},
where {it:coef_rt} is the estimated coefficient in relative time period {it:rt} and {it:time_rt} is the relative time period; the variable {it:coef_rt} includes values of 0 for the
reference period (or periods, if {it:timetrends} and {it:esfixedbaseperiod} are specified). These tests will not be conducted if an insufficient number of relative time
period-specific treatment effects are requested (either overall, or within specific subgroups).
{p_end}

{phang}{it:Size of the adversarial test.} Simulation evidence indicates that these line tests are particularly sensitive to small cluster count issues, and that they tend to be
severely oversized when each i-by-t cell does not contain {it:many} clusters. Rejection rates under a true null were close to nominal when every treated cohort-by-period cell was
spanned by six or more clusters, but rose sharply as that number fell - roughly threefold when cells were spanned by only two clusters, and by an order of magnitude when each cell
lay entirely within a single cluster (as happens automatically when clustering on the cohort variable {it:i} itself). Enlarging cells does not substitute for this: cell size was
found to be immaterial to the test's size once the cluster requirement was met, and did not repair matters when it was not. Small overall cluster counts degrade the test
independently. Two further conditions are worth keeping in view. Under a continuous treatment the line test additionally requires cells containing many observations - on the order
of fifty - because each cell's dose slope is otherwise too imprecisely estimated for the test's covariance to be trustworthy. And where errors contain a common shock at the
cohort-by-period level, no analytic variance estimator can capture it: such a shock is absorbed into the cell's own coefficient and leaves no trace in the residuals, so the test
becomes badly oversized regardless of the clustering chosen, increasingly so as cells grow. Where these conditions are not met, prefer randomization inference (see {it:permute()})
or bootstrap-based standard errors for pre-trend inference. {it:wooldid} does not test for these conditions at run time; they are a matter for the user's judgement about the
design. {p_end}

{phang}{opt altmethod_adversarial}: Alternate method of implementing the test above, which projects a trend forward using the pre-treatment coefficients only (plus the planted
zero-valued reference row or rows), and then tests whether the post-treatment coefficients collectively fall in line with the projected values. The test propagates the full
forecast covariance: writing the post-treatment coefficients as {it:b_Q}, the pre-treatment fit set as {it:b_P}, and the linear operator taking {it:b_P} to the forecasts as {it:A},
the forecast error {it:b_Q - A b_P} is tested against its complete covariance {it:V_QQ + A V_PP A' - A V_PQ - V_QP A'}, so uncertainty in the fitted pre-period line and the
pre/post covariance both count against the null. Like the primary implementation, the p-value requires {it:esfixedbaseperiod}. {p_end}

{phang}{opt pretrendpow:er(numlist)}: Reports Roth (2022)-style pre-trend power diagnostics for each requested power level (specified in percent as integers strictly between 5 and 100,
e.g. {it:pretrendpower(50 80 90 95)}): a GLS estimate of the linear pre-trend slope with its standard error and confidence interval (fit through the reference period on the pre-treatment
event-study coefficients using their full variance-covariance matrix from {it:e(V_es_{it:type})} — a more efficient and validly-inferenced line than the display-oriented
{it:e(adversarial_line_{it:type})}, which is an unweighted fit over all horizons, or over the pre-period horizons only under {it:altmethod_adversarial}); the minimum linear pre-trend each pre-test would detect with the requested power ({it:mdt_wald} powers
the joint pre-trend test, an F test with the regression's residual degrees of freedom; {it:mdt_nis} powers the test that any pre-period coefficient is individually significant — the
convention of Roth's {it:pretrends} package, whose size is also reported since it exceeds its nominal level when there are several pre-periods); and the bias each such trend would imply
for that block's pooled estimand {ul:if it continued into the post-treatment period} (the minimum detectable trend itself depends only on the design's precision, never on the estimated
coefficients). When only {it:att_i} or {it:att_itime} is requested, the reported block is the shared {it:att_it} event-study curve ({it:ame_it} under {cmd:contreat()}) standing in for them: its
hbar and bias columns describe {it:att_it}-family post-period weighting only and do {ul:not} quantify the bias such a trend would imply for the requested {it:att_i}/{it:att_itime}
(or {it:ame_i}/{it:ame_itime}) estimands, whose post-period weights differ. With {cmd:subgroups()}, diagnostics are computed for the overall event-study types only ({it:att}, {it:att_it},
{it:ame}, {it:ame_it}); subgroup event-study objects ({it:e(b_es_*_s#)}/{it:e(V_es_*_s#)}) are not analyzed and no {it:e(roth_*_s#)} results are returned. Requires {cmd:esfixedbaseperiod} (under the pooled reference the single omitted reference period the linear-trend normalization needs is not defined) and at least one
pre-treatment event-study coefficient; not compatible with {cmd:impvalues}. Results are returned in {it:e(roth_power_{it:type})} and companion scalars (see Stored Results) and skipped
with a caution when the pre-period variance-covariance block is unavailable (e.g. under {cmd:timetrends}) or carries missing or non-positive variance entries; rank-deficient blocks are
handled with a generalized inverse and reduced degrees of freedom, with a caution (note that on an exactly rank-deficient block every reported statistic reflects only the linearly
independent subset of pre-period coefficients the generalized inverse retains - not a summary of all pre-period evidence - and depends on which columns are retained, so treat those
numbers as one representative of an underdetermined family and interpret with care). A requested power level at or below the NIS test's actual size is achievable by a zero trend, so {it:mdt_nis} and {it:bias_nis} are exactly 0 there (a missing cell always means the
computation was unavailable, never a boundary answer). The exact-quadrature {it:mdt_nis} column is computed only for up to 4 pre-periods (its cost grows
steeply with dimension); for larger designs, feed {it:e(b_es_{it:type})}/{it:e(V_es_{it:type})}/{it:e(numpre_{it:type})} to Roth's {it:pretrends} package (subject to the same positional-geometry conditions as the honestdid handoff — see the {it:e(b_es)} Stored Results entry), which implements the same
convention with a normal critical value. When the underlying estimator posts no residual degrees of freedom (e.g. the
ppmlhdfe path), the joint test and minimum detectable trends use chi-squared rather than F distributions and the slope CI uses a normal critical value. Under {cmd:poisson} without
{cmd:poisexpresults}, the diagnostics are on the linear-index (xb) scale like the event-study coefficients themselves. {p_end}

{phang}{opt makep:lots}: Produce event study plots. Note on graph management (applies to all of {it:wooldid}'s plotting options): {it:wooldid} names the graphs it draws
(e.g. att_fullmodel, hist_att_i, cfxatt_cfxpost) and creates them with the graph-name {it:replace} suboption, so re-running overwrites {it:wooldid}'s own prior graphs of the
same name; it does not clear other graphs from memory. Two consequences: a user graph whose name exactly matches a {it:wooldid} graph name will be overwritten when that plot
is recreated, and {it:wooldid} graphs from an earlier, differently-configured run linger in memory until that configuration is run again (or the user drops them). {p_end}

{phang}{opt oos:treatedcontrols(integer)}: By default, {it:wooldid} uses cohorts treated after the end of the sample the same way as never-treated control cohorts (i.e., such cohorts are
"out-of-sample treated controls"). When this option is specified with a value of T, cohorts treated within T periods of the end of sample are instead handled the same as other treated cohorts.
Such cohorts are potentially allowed to contribute to pre-treatment relative time period-specific effects when this option is specified,
but not otherwise. Note that if {it:esfixedbaseperiod} is specified, the out-of-sample-treated cohorts will be dropped if the specified reference period is not usable for these cohorts
(not observed within the sample, or removed because no comparison unit exists in that calendar period). When {it:anticipation()} is also specified, the out-of-sample window applies to the {ul:effective} (shifted) onset, so a cohort whose recorded onset
lies just beyond the window can become usable once its effective onset falls within it. {p_end}

{phang}{opt antic:ipation(integer)}: Treats the effective treatment onset as {it:ttre} − #, so that periods in which units may act in anticipation of a known upcoming treatment are counted
as treated periods rather than clean pre-treatment periods (which would otherwise contaminate the reference period and pre-trend estimates). The shift is applied once, before relative time
is constructed, and is exactly equivalent to re-running {it:wooldid} after replacing {it:ttre} with {it:ttre} − # (this is also the convention of {it:did_imputation}'s {it:shift()} option):
all event-study horizons are counted from the {ul:effective} onset, so relative period 0 is # periods before the recorded onset, and you may wish to raise {it:espostlength()} by # to keep
the same number of genuine post-treatment horizons. Note an inherent asymmetry: specifying too SMALL an anticipation length leaves a loud residual pre-trend signal, but too LARGE a value attenuates the estimated effects with clean pre-trends and no diagnostic - prefer the smallest defensible value. Notes: # is expressed in the units of {it:t} (not in panel-period counts — with, e.g., biennial data an odd # lands between observed
periods, and {it:wooldid} prints a caution); units whose effective onset falls at or before their first usable period are dropped (with a caution reporting the count), and under
{cmd:esfixedbaseperiod}, units whose {ul:shifted} reference period falls on a period unusable for them (not observed, or dropped for lack of a comparison unit) are dropped, with a caution when the drop is attributable to the shift
itself; with {cmd:subgroups({it:ttre})}, subgroup labels remain the recorded (unshifted) onset values. When active, the shift is recorded in {it:e(anticipation)}. Negative values are
not allowed. {p_end}

{pstd}{ul: Two Approaches to Estimation -- Pooled Reference Period vs Fixed Reference Period:} Difference-in-differences treatment effect estimates involve comparing average outcomes during
treated cohorts' treatment periods to treated cohorts' average outcomes during a pre-treatment reference period. By default, {it:wooldid} does comparisons using a pooled reference period,
meaning that treated periods for a given cohort are compared to a reference period containing all pre-treatment observations for the cohort. If the user specifies the option
{it:esfixedbaseperiod}, {it:wooldid} will switch to employing a fixed reference period, meaning that all treatment effect estimates will be formed by comparing treated cohorts' average outcomes
in post-treatment periods to their outcomes in a particular pre-treatment reference period (specified in relative time period terms). By default, this pre-treatment reference period is the
time period immediately prior to treatment (relative period -1), though the user can specify alternative behavior. For more on how this is implemented, see {help wooldid##description:discussion}.{p_end}

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
The SG variable must contain only non-negative integers (or missing values) on in-sample observations of treated units - both pre- and post-onset rows; values on never-treated observations and on rows wooldid's deterministic screens discard (e.g. missing outcomes) are ignored by the check, since no subgroup computation uses them and its levels must specify collections of cohorts i. A given cohort i cannot be split across multiple
subgroups: {it:wooldid} exits with an informative error if it is. Under {it:safety(off)} such a run proceeds using share weighting (see the subgroup discussion below). Observations given missing values are not placed in any subgroup and are ignored when performing subgroup computations. Because subgroup values are
embedded in the names of stored results and plots, values 8 or more digits long can hit Stata's 31/32-character name limits and crash some option combinations with r(198); recode
ID-style subgroup codes to small integers first (e.g., {cmd:egen newsg = group(oldsg)}).   {p_end}

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

{pstd} Each treated unit of {it:i} should lie wholly within one subgroup. Because subgroups never enter the regression itself (they only select and reweight the unit-by-period effect cells),
a subgroup variable that varies {ul:within} a unit-by-period cell makes subgroup contrasts structurally unable to reflect within-cell subgroup differences, and {it:wooldid} exits with an
informative error (bypass with {it:safety(off)}; the usual remedy is defining {it:i} at the finer level over which the subgroup varies, e.g. {it:egen newi = group(i sg)}). A subgroup variable
that changes over time within a unit is likewise refused under default safety. With {it:safety(off)} such a run proceeds with share weighting - {it:att_i} and {it:att_itime} enter a split
unit into each subgroup weighted by its share of sample weight there rather than with equal unit weighting, while by-subgroup {it:att} and {it:att_it} are unaffected - and a caution prints. Observations with missing subgroup values remain allowed and belong to no subgroup. {p_end}

{pstd} A special case to consider is the difference-in-differences-in-differences (triple diff) case. In the triple diff case, all cohorts are categorized into two groups, with both groups
existing among treated cohorts and control cohorts. The target estimand may be thought of as the difference in the diff-in-diff estimates obtained separately for each of the two groups. To
implement a triple diff design in {it:wooldid}, one can specify the variable containing the two triple diff groups as a SG variable and then manually specify inclusion of SG by period
fixed effects (i.e., specifying {it:subgroups(sgvarname) fe(i.sgvarname##i.t)}). The difference in SG level estimates returned will be the triple diff estimate. {p_end}


{marker continuous}{...}
{dlgtab:Continuous Treatments}

{phang}{opt contreat(varname)}: The continuous treatment variable. Must be defined for all observations (set to 0 if not naturally defined within untreated cohort-periods).
With continuous treatments, reliable analytic inference additionally requires that cohort-by-period cells contain a large number of observations relative to the number of
within-cell parameters estimated. If the continuous treatment carries no identifying variation for any treated cohort-period cell's dose slope (for example, a dose constant
within each cohort - every treated-cell slope term is then dropped for multicollinearity), the within-cell average marginal effects are not identified: under default safety
{it:wooldid} names this situation, reports the ame results as missing rather than as the regression's structural zeros, and continues with the unaffected att results
(expected zero doses in untreated periods, as under {it:contreatcontrols(untreatedZero)}, do not trigger this).{p_end}

{phang}{opt contreatelas:ticitytype(dydx/eyex/eydx/dyex)}: Specifies which type of treatment effect to compute. {it:eyex} corresponds with computing the elasticity
(% change in y/% change in contreat), {it:eydx} corresponds with semi-elasticities comparable to what can be produced for the main effect of treatment
(% change in y/unit change in continuous treatment), {it:dyex} gives the other semi-elasticity (change in y/% change in continuous treatment), and {it:dydx} is the default
(not an elasticity). A supplied value outside this list is rejected with an error. Note that successful estimation of these elasticities tends to require a larger sample, as well as y and contreat variables that are not sparsely distributed. {p_end}

{phang}{opt contreatw:ithin}: Specify this when you wish to use only within-treated cohort variation in the continuous treatment variable. This option disables some of {it:safety}'s features (the estimability-relaxation checks; the data-integrity screens, such as onset constancy within units, still run unless {it:safety(off)} is given explicitly), turns off
estimability checks that would otherwise cause {it:wooldid} to stop estimation in the absence of any untreated cohorts, and causes {it:wooldid} to skip sample cleaning routines that drop time
periods for which a control cohort is not present. Requires {opt contreat(varname)}: because this option disables those checks, requesting it without a continuous treatment is rejected with an error rather than silently ignored.   {p_end}

{phang}{opt contreatp:oly(integer)}: Includes polynomial terms in {it:contreat} of the specified degree, which must be a positive integer (default 1; 0 and negative values are rejected). Be cautious using higher degree polynomials. {p_end}

{phang}{opt lattice(latticevar OR latticetype latticenum)}: Specifies a lattice variable, a variable specifying groups of observations across which the effect of {it:contreat} is
allowed to vary, even within cohort-periods and within the reference period. If just the name of a variable is passed to {it:lattice()}, this variable will be used as the lattice variable.
If instead a lattice type and a lattice number are passed, this option will generate a lattice variable using {it:latticenum} quantiles, with the quantiles being formed within groups
according to the {it:latticetype} (see below). Any specified sample weights will be used when constructing the lattice, unless the option {it:latticeignoreweights} is specified (which
constructs the lattice quantiles without using weights). 
If {it:lattice} is not specified, {it:wooldid} will proceed without
using a lattice variable (as though you specified a lattice variable that is always 1). A user-supplied lattice variable must be numeric and contain only non-negative
integers (it is used as a factor variable; the value check can be bypassed with {it:safety(off)}); a lattice type passed without a quantile count is rejected with a
specific error at parse time. Be aware that a lattice is the most expensive scaling axis in {it:wooldid}: runtime grows at roughly the 1.9th power of the
cohort-by-period cell count once a lattice is present (versus about the 1.3rd to 1.7th power without one), and the absorbing regression rather than
{it:margins} becomes the dominant cost. Each additional lattice point multiplies the terms the regression must carry, so prefer the smallest
{it:latticenum} that answers the question. See {help wooldid##runtime:Runtime and what drives it}. {p_end}

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
the lattice variable will be used. The variable must be numeric and, on observations with a nonmissing treatment onset within the sample restriction, contain only
non-negative integers ({it:margins} accepts only non-negative integers in {it:over()}; the value check can be bypassed with {it:safety(off)}). Values on never-treated
observations are ignored: no continuous treatment effects plot uses those rows. If the variable has missing values on observations with a nonmissing onset, {it:wooldid}
prints a Caution: {it:margins} silently omits those observations from the plotted bin averages. These screens apply only when {it:makecfxplots} is actually requested; {it:cfxlattice()} without {it:makecfxplots} has no effect and draws a Caution.
Note that {it:makecfxplots} requires a continuous treatment ({it:contreat()}); previously, requesting it with {it:cfxlattice()} alone was silently ignored. {p_end}

{phang}{opt cfxplott:ypes(pre post integers)}: A list of types of cfxplots to be produced. If {it:post} is among the types specified, plots will be produced giving
treatment effects across the continuous treatment variable within the post-treatment period. If {it:pre} is among the types specified, plots will be produced within the pre-treatment period.
If any integers are specified, plots will be produced just within the specified periods relative to treatment. If nothing is specified, the option post is assumed. {it:pre} and
relative-period integers outside the reported window require an event-study window: without {it:esfixedbaseperiod} or {it:esprelength()} > 0 they are dropped from the request. {p_end}

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

{pstd} {ul:On Specifying contreatcontrols() and lattice:} As discussed {help wooldid##description:here}, the regression equation {it:wooldid} adopts when estimating
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
cohort-time period for which a treatment effect is estimated (if {it:esfixedbaseperiod} is specified, that includes every estimable treated cohort by relative time period effect). Note that
reported treatment effects will not vary as a direct function of these controls; the effect of these controls are partialled out from the reported treatment effects.
This option differs from simply specifying cohort
by time period interactions with a given control in that it pools the effect of the controls in the reference period and the control group. Note that all specified coarse cohort controls must be
specified using factor notation, so that {it:wooldid} knows if they are to be interacted as a continuous variable (prefixed with "c.") or as a categorical variable (prefixed with "i."). {p_end}

{phang}{opt inter:activecontrols()}: Includes the specified variables in the regression in two capacities. First, the variables are included as controls interacted with time period-specific
indicator variables. Second, the variables are de-meaned by cohort (i) and then treatment effects are permitted to vary linearly with the de-meaned variables. This is achieved by including
interactions between the de-meaned controls and all Tz_it dummies included in the model. The specification used is similar to as with {it:coarsecohortcontrols}, but differs in that these terms
are marked so that variation associated with them is reported as part of the treatment effects (as opposed to being partialled out from the treatment effects). This specification most closely
mimics the approach to including covariates that appears in Wooldridge (2025). When using only one interactive control, ATTs obtained using this option should be similar to those obtained
when specifying that variable as a continuous treatment variable and then separately including controls that interact the continuous treatment variable with time dummy variables. {p_end}


{pstd} {ul:Caution Regarding Continuous Treatment Variables and Controls:} If seeking to manually specify additional controls as a function of a continuous treatment variable, be sure not
to include the continuous treatment variable itself among your controls, as this can adversely affect estimates produced in {it:wooldid}'s margins calculation stage. It is okay to include
a control with values identical to that of the continuous treatment variable: it just needs to have a different variable name. It is recommended that any custom such controls be specified
in conjunction with the option {it:contreatcontrols(basic)}, or any option other than {it:untreatedZero}. {p_end}

{pstd} {ul:Caution Regarding Continuous Treatments and Cluster-Robust Standard Errors:} Valid clustered standard errors cannot be obtained analytically for the ame-family
quantities when clustering at the cohort ({it:i}) level or above, necessitating use of bootstrapping or other methods (detailed below).
With {cmd:contreat()}, the regression gives every treated cohort-by-period cell its own
treatment-intensity slope. When any such cell lies entirely within a single cluster of the specified {cmd:cluster()} (which is automatic when clustering on the cohort identifier
{it:i} itself), the cluster-robust variance estimator is structurally blind to within-cell variation in the continuous treatment: the "standard errors" it would deliver for the ame-family
quantities do not measure their sampling variability (in the {it:contreatcontrols(basic)}/{it:byTime} cases they collapse to the variance of the absorbed baseline dose control, one number
repeated across every ame quantity — including configurations where the covariance appears full rank). {it:wooldid} detects this design, prints an explanatory caution, sets the ame-family
standard errors, tests, and joint diagnostics (adversarial ame p-values, ame joint pre-trend tests, subgroup ame comparisons, ame pre-trend power) to missing, and withholds
{it:e(V_es_ame*)} while still returning the unaffected point estimates and {it:e(b_es_ame*)}. Valid ame inference under such designs: bootstrap standard errors, IM p-values
({it:impvalues} — cross-cluster dispersion inference, which remains valid), randomization inference, {it:vce(robust)} when errors are independent across observations (not appropriate
under within-unit serial dependence), or clustering on a variable strictly finer than the cohort-by-period cells. Note that binary-treatment (att) standard errors under cohort-level clustering are also somewhat understated on such saturated designs (a milder, documented
phenomenon); wooldid does not suppress them. {p_end}


{marker technical}{...}
{dlgtab:Technical Features}

{phang}{opt update}: If specified, the command {it:wooldid,} {it:update} will update your local installation of {it:wooldid} based on the most recent edition of the program available on github. {p_end}

{phang}{opt customw:eightmultiplier(varname)}: When computing various treatment effects using the margins calculation (see step 2 {help wooldid##description:here}), the weights normally 
used to produce the requested type of treatment effect (ATT, ATT_i, etc.) will be multiplied by this variable. When requesting an ATT without specifying sample weights, this is 
equivalent to estimating all regressions without weights, but then using these weights when averaging treatment effects to obtain the ATT. This allows production of certain 
bespoke types of treatment effects. Values must be non-negative on in-sample observations (negative values would silently corrupt the weight masses and estimates; they are rejected
with an error, bypassable with {it:safety(off)}); missing values drop the observation. Under {cmd:impvalues} the reported event effects are equal-weight averages of per-cluster
margins, so the {it:e(escomp)} weight-share columns describe sample weight composition, not the Ibragimov-Muller cross-cluster aggregation. {p_end}

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
efficiency in the main regression or under certain other circumstances. Note: be extra careful to ensure that specifications featuring ccc_absorb() are identified; not all cases
of identification failure in this case are bound to be automatically detected. {p_end}

{phang}{opt ver:bose}: Prints the regression command {it:wooldid} uses, displays raw output from the regression estimated, and prints progress reports throughout the estimation process.
This option also stores the underlying regression model used by {it:wooldid} to {it:___wooldidfullmodel}; this model can be accessed using {it:estimates restore}. {p_end}

{phang}{opt gte:ffects}: Returns (and, for tables of up to 200 rows, prints) the individual cohort-by-relative-period effect estimates underlying {it:wooldid}'s aggregates
(the exact building blocks for the plain binary designs the identity notes below describe; under transformations, time-varying interactive controls, or {cmd:impvalues} the aggregates
are not simple reweightings of these rows), in {it:e(gteffects)} — one row per estimated {it:i}-by-relative-period cell, read directly off the underlying regression's coefficient vector at essentially zero additional
computation. Note that a coefficient dropped by the regression (collinearity omission) simply has no row: the table does not distinguish a partially identified horizon from a complete
one on its own - under default safety the dropped-coefficient screen errors first, but under {it:safety(off)} compare the row count against the expected cohort-by-period cells. Columns: {it:i} (the cohort/unit identifier value), {it:onset} (that unit's effective treatment onset, after any {cmd:anticipation()} shift), {it:caltime} (the calendar
period the cell covers), {it:relyr} (period relative to onset), {it:estimate}, {it:se}, {it:p}, {it:lb_95ci}, {it:ub_95ci}, {it:inwindow} (1 if the cell falls inside the reported
event-study window), {it:nobs} (estimation-sample observations in the cell), and the cell's total estimation weight under each requested estimand ({it:w_att}, {it:w_att_it},
{it:w_att_i}, {it:w_att_itime}; missing for estimands not requested). The weighted mean of {it:estimate} under any of these weight columns over the relevant cells reproduces the
corresponding aggregate exactly, so the table supports arbitrary custom aggregations and is directly comparable to the group-time tables of {it:csdid} or {it:xthdidregress}
(the cells are ATT(g,t)-style group-time effects when {it:i} indexes cohorts; with a fine unit-level {it:i} they are unit-by-period effects and the table can be large). Notes:
standard errors are the regression's coefficient-level (conditional) VCE — {it:unconditionalse} and {it:impvalues} do not alter them; under {it:poisson} the estimates are on the
linear-index (xb) scale ({it:poisexpresults} and {it:semielasticity} are not applied); under {it:timetrends} the coefficient-level variances are typically not identified and the
{it:se} column is missing (degenerate zero variances are always reported as missing, never zero); with {cmd:subgroups()} a single full-sample table is returned; with time-varying
{cmd:interactivecontrols()} the interaction coefficients absorb control variation and the cell-level reading is unreliable. Not available with {cmd:contreat()}, where the
cohort-by-period interaction coefficients are dose-zero intercepts rather than treatment effects. Fully available under {it:fastnose} (the regression's own VCE is intact even
though the {it:margins} stage skips variances), making {it:gteffects} the one inference surface that survives fast mode. If no cohort-by-period coefficients can be read off the
fit (heavily collinear designs), no table is returned and a caution prints. Two reading notes: the weighted-mean identity itself
remains exact even under time-varying {cmd:interactivecontrols()} (the caveat there concerns interpreting individual cells, not the aggregation); and when the reported
event-study window has no post-treatment horizons (e.g. {cmd:espostlength(0)}), {it:inwindow} is 0 on every post-treatment row by construction. {p_end}

{phang}{opt emptycellsoverride}: By default, {it:wooldid} tells Stata to drop interaction terms that correspond with empty cells, which is achieved by the command
{it:set emptycells drop}, which improves estimation efficiency in this setting. This option overrides this behavior. This should only be of interest in conjunction
with the option {it:verbose}. Note that it is an expensive option: dropping empty cells was measured to be worth a factor of roughly 2.24 in runtime on a
representative specification, and overriding it forgoes that saving. {p_end}

{phang}{opt safety(on|off)}: When {it:safety(off)} is not explicitly specified, {it:wooldid} will caution the user or halt estimation when {it:wooldid} is likely to yield unreliable estimates.
Only {it:on} and {it:off} are accepted values; anything else is rejected with an error (formerly, unrecognized values were silently treated as {it:on}).
Most important among these checks are tests for whether or not a variable intended to capture the treatment effect in a given i-t cell was dropped from the underlying regression for
multicolinearity, as may occur if one's fixed effects or continuous treatment control variables partial out one of the main effects. Identification failures of this sort will generally be detected,
though not necessarily when the problem concerns continuous treatments being analyzed with a lattice or with polynomial terms. Default safety also screens for very granular {it:i} variables:
when {it:i} takes more than 100 distinct values, {it:wooldid} cautions that runtime grows faster than linearly in the cohort-by-period cell count built from {it:i} and that the specification
can also exceed Stata's variable limits (consider unit-level fixed effects via {it:fe()} plus a coarser {it:i}, such as treatment onset groups; see
{help wooldid##runtime:Runtime and what drives it}), and if the count of variables the specification requires exceeds the current
{help maxvar} setting, {it:wooldid} halts with an informative error before estimation rather than failing mid-run with r(907). Default safety additionally refuses {it:jointtests} combined
with analytic clustered SEs clustered at the level of the cohort variable {it:i} or above, where the joint test p-values are severely oversized (see {it:jointtests}); {it:safety(off)}
downgrades that error to a caution and proceeds.   {p_end}

{phang}{opt oldsyntax}: In older versions, {it:wooldid} used the suffix {it:contemp} to flag the event study coefficient containing the period contemporaneous with treatment. 
This turned out to make {it:wooldid} relatively incompatible with the popular event study plot production program {help event_plot}, so the {it:contemp} suffix has been 
replaced with the suffix {it:post0}. This option reverts {it:wooldid} to using the older coefficient syntax, for compatibility with any older code that relies on it. {p_end}




{marker runtime}{...}
{title:Runtime and what drives it}

{pstd}Users are often surprised by which specifications are slow, because the intuitive lever - the number of observations - is not the binding one. The guidance
below comes from direct measurement across a range of panel sizes and shapes. {p_end}

{pstd}{bf:Cost is driven by the cohort-by-period cell count, not by sample size.} Write K for the number of interaction cells {it:wooldid} builds, which is
roughly the number of distinct values of {it:i} times the number of relative periods in the event-study window. Runtime grows at about the 1.3rd to 1.7th power
of K, while sample size, subgroup count and event-study width are each merely linear. The effect is large enough to invert expectations: in one measurement a
29,067-observation panel with K=29 ran in 1.4 seconds, while a 4,152-observation panel with K=2,257 took 86.8 seconds - seven times fewer rows and sixty-three
times more time. {p_end}

{pstd}{bf:What this implies in practice.} Coarsening {it:i} is normally the single largest speedup available: move to cohorts or treatment-onset groups and
recover unit-level variation with {it:fe()}. Halving the number of observations, by contrast, buys little, because it leaves K untouched - this is also why
subsampling is a poor way to shorten a long {cmd:bootstrap()} or {cmd:permute()} run, where shrinking the window or coarsening {it:i} works and dropping rows
largely does not. Narrowing the event-study window helps in proportion to the horizons removed. {p_end}

{pstd}{bf:The expensive options.} {it:lattice()} is the steepest axis, at roughly the 1.9th power of K, and it shifts the dominant cost from {it:margins} into
the absorbing regression; prefer the smallest {it:latticenum} that answers the question. {it:emptycellsoverride} forgoes a measured factor of about 2.24.
The {cmd:bootstrap()} re-identification leg required by {it:att_it}, {it:att_i} and {it:att_itime} roughly doubles per-replication cost. On the other side,
{it:fastnose} skips the variance calculation that accounts for most of a {it:margins} call and is correspondingly worth several times over on runs where
standard errors are not needed. The default cell collapse (see {it:fastcollapse}) attacks the same cost from
the other side - it shrinks the data each {it:margins} call sees rather than skipping the variance - and the two compose. The collapse is already on unless
{it:fullprecisionslow} or an ineligible option turns it off, so a run that must be evaluated on the full data ({it:unconditionalse}, {it:contreatpoly()} above 1, ...,
or a {it:subgroups()} run whose screen has set the collapse aside) is the one paying the full per-row {it:margins} cost described here. {p_end}

{pstd}{bf:Where the time goes.} Nearly all of a typical run is spent inside Stata's {it:margins}, and most of each {it:margins} call is the delta-method
variance calculation rather than the point estimates. Very little time is spent in {it:wooldid}'s own code, which is why the levers that matter are the ones
that change the size of the problem handed to {it:margins} - K above all - rather than anything about how the data are stored or passed in. {p_end}

{marker output}{...}
{title:Working With Output}

{pstd}
{it:wooldid} generates output primarily in the form of matrices and scalars stored stored in e(), though all main results can also be accessed using standard utilities like {it:esttab}.
The commands {it:est restore} and {it:est save} are compatible with {it:wooldid} for saving its output. The full set of {it:wooldid} output in e() can be viewed using the command
{it:ereturn list}. {p_end}

{pstd} The primary results are stored in the following matrices, scalars, and macros: {p_end}
{synoptset 30 tabbed}{...}
{synopt:{it:e(wdidmainresults)}} All primary and event study results.{p_end}

{synopt:{it:e(wdidsgresults_{it:sgcase})}} Primary and event study results for the subgroup {it:sgcase}, if subgroups are used and presentation of their
primary results has not been suppressed.{p_end}

{synopt:{it:e(pwsgs)}} Pairwise differences between the subgroup-specific treatment effects, if subgroups are used.{p_end}

{synopt:{it:e(b_es_{it:type})}, {it:e(V_es_{it:type})}} The event study coefficient vector and its {ul:full} variance-covariance matrix (carrying genuine off-diagonal covariances, unlike {it:e(V)}), relative-period-labeled as {it:pre#}/{it:post#} with the reference period omitted. Returned per estimand type: {it:att} (binary) and, with {cmd:contreat()}, {it:ame}, plus {it:att_it} / {it:ame_it} for the {it:att_it}, {it:att_i}, and {it:att_itime} estimands (which share a single event-study horizon internally). With {cmd:subgroups()}, per-subgroup versions {it:e(b_es_{it:type}_s{it:sgcase})} / {it:e(V_es_{it:type}_s{it:sgcase})} are returned as well. After {it:ereturn post} these matrices are valid inputs for {it:lincom} and {it:test}. A {ul:direct positional handoff} to {it:honestdid} (Rambachan-Roth robustness) additionally requires ALL of:
{cmd:esfixedbaseperiod} with the default {cmd:esrelativeto(-1)}; no {cmd:timetrends}; at least one pre- and one post-treatment coefficient; horizons consecutive apart from the intended
jump over the omitted -1 reference (checkable signature on the colnames: {it:pre1} absent, rightmost pre column is {it:pre2}, {it:post0} present, no other gaps); and a covariance with no
structurally degenerate rows (see {it:e(es_ndegen_{it:type})} below; honestdid inference inherits whatever VCE you specified). {it:numpre} is purely positional and cannot communicate
gaps, non-(-1) references, or the (default) pooled cohort-specific reference — those geometries hand honestdid the wrong event-time geometry {ul:silently} and need a user-supplied
transformation. Check the regime via {it:e(estimate_type)} == "fixed_reference_period" and the numeric spacing via the {it:relyr} column of {it:e(escomp)}/{it:e(escomp_s{it:sgcase})}.
The matrices are returned whenever at least one horizon has an identified variance; horizons whose standard errors are not identified then carry zero (structural-zero) or missing
diagonal rows, flagged by a Caution and by {it:e(es_ndegen_{it:type})} (the count of degenerate rows, posted only when such rows exist — its absence marks a clean covariance). Fully
non-identified VCVs ({cmd:impvalues}, {cmd:fastnose}, typically {cmd:timetrends}) are not returned. Additionally, when the cell-containment caution fires (see the continuous-treatments
standard-error caution above), {it:e(V_es_ame*)} is withheld — the captured covariance is not the sampling covariance of the ame estimates — while {it:e(b_es_ame*)} still posts; att-side
event-study objects are unaffected. Separately, {it:contreatcontrols(byCohort)} requires each treated cohort's untreated/reference-period observations to carry variation in the
continuous treatment (that variation is all that identifies the cohort-specific baseline dose slope): wooldid refuses with r(459), naming the offending cohort, when a treated cohort has
none, and prints a computational-fragility caution (suggesting a larger {it:regtol()}) when that variation is below 1/1000 of the sample standard deviation of the treatment variable.{p_end}

{synopt:{it:e(numpre_{it:type})}} The number of pre-treatment relative periods (relative year < 0) in the corresponding event-study vector, provided for use as {it:honestdid}'s {it:numpre()} argument (subject to the geometry conditions in the {it:e(b_es)} entry above).{p_end}

{synopt:{it:e(gteffects)}} The cohort-by-relative-period building-block table produced by the {it:gteffects} option; see that option's description for the column layout.{p_end}

{synopt:{it:Bootstrap inference}} When {cmd:bootstrap()} is used: scalars {it:e(boot_reps)}, {it:e(boot_completed)}, {it:e(boot_failed)}, {it:e(boot_level)}, {it:e(boot_N_V)} (replications behind the covariance) and {it:e(boot_v_posted)}; macros {it:e(boot_method)}, {it:e(boot_stat)} (z, or t under {cmd:bootdf()}), {it:e(boot_seed)}, {it:e(boot_alpha)} and {it:e(boot_ties)}; matrices {it:e(V_boot)}, {it:e(b_boot)} (replication means), {it:e(bias_boot)}, {it:e(boot_navail)} (available replications per estimate), {it:e(boot_ci_normal)}, {it:e(boot_ci_percentile)}, {it:e(boot_ci_bc)}, {it:e(boot_z0)}, and - up to 10,000 completed replications - {it:e(bootdist)}, the replication matrix itself. The analytic covariances the bootstrap replaced are kept as hidden {it:e(V_analytic)} and {it:e(V_es_analytic_*)}. Under {cmd:bootdf()} the scalar {it:e(boot_df)} records the degrees of freedom in force. {it:e(df_r)} is aligned with the reference distribution actually reported, so a {cmd:test} or {cmd:lincom} run by hand after a bootstrap fit reproduces the p-values in the reported tables under every {cmd:bootdf()} setting: under {cmd:bootdf(clusters)} or an explicit integer it holds those degrees of freedom, and under {cmd:bootdf(normal)} it is set to missing so that {cmd:test} and {cmd:lincom} fall back to chi-squared and z, matching the reported statistics. The analytic residual degrees of freedom the bootstrap replaced are kept as hidden {it:e(df_r_analytic)}, alongside the analytic covariances. {p_end}

{synopt:{it:Permutation inference}} When {cmd:permute()} is used: scalars {it:e(perm_p_{it:est})}, {it:e(perm_lb_{it:est})}/{it:e(perm_ub_{it:est})} (the 5th-95th percentile null band),
{it:e(perm_sd_{it:est})}, and under {it:permtci} {it:e(perm_tlb_{it:est})}/{it:e(perm_tub_{it:est})}, for each requested estimand; run-level {it:e(perm_reps)}, {it:e(perm_completed)}, {it:e(perm_failed)},
{it:e(perm_npossible)} (also as the string {it:e(perm_npossible_str)} for very large counts); macros {it:e(perm_method)} and {it:e(perm_stat)} (b or t); the comparator matrix
{it:e(perm_comparator)} (1 x k in estimand order, on the {it:e(perm_stat)} scale: the same-path rep-0 comparison statistics actually used for the tail counts - computed on the fast
no-SE path unless {it:permtstat}, so they can differ from the displayed estimates in the last digits; a failed or missing rep-0 statistic is refused loudly rather than substituted -
with p-values exactly reconstructable as count/{it:e(perm_completed)} of draws with |draw| >= |comparator| for exact enumerations, and (1 + count)/({it:e(perm_completed)} + 1),
add-one in both numerator and denominator, for sampled draws); and the draw matrix {it:e(permdist)} (completed
replications by estimand, returned for up to 10,000 draws).{p_end}

{synopt:{it:e(escomp)}} Event-study composition diagnostics, one row per reported horizon (rows align one-for-one with the ES rows of {it:e(wdidmainresults)}; horizons without a surviving
coefficient get no row). Columns: {it:ncohorts} (distinct treatment-onset values contributing), {it:nunits} (distinct units of {it:i} — equivalently, the number of unit-by-period
coefficients averaged at that horizon; that identity is enforced for plain binary treated-cell coefficients under default safety by the hard dropped-coefficient error — it can fail under
{it:safety(off)} or {it:contreatwithin}, and in continuous/interacted designs the screen cautions rather than errors), {it:nclusters} (distinct clusters, when a single {cmd:cluster()} variable is used; missing otherwise), {it:nobs} (estimation-sample observations),
{it:sumw_att} (total {it:att}-style weight: sample weight times any {cmd:customweightmultiplier()}), {it:wshare_att} and {it:wshare_att_it} (the horizon's share of the reported-window
treated weight under the {it:att} and {it:att_it} weighting concepts respectively; each column sums to 1 across the reported horizons), and {it:relyr}. With {cmd:subgroups()},
per-subgroup versions {it:e(escomp_s{it:sgcase})} are returned as well. Composition shifts across horizons are a common cause of jumpy event-study profiles; note this is a descriptive counts-and-totals screen,
not a test — flat composition columns rule out changes in these {ul:reported} summaries (a cohort entering or leaving, a unit-count or weight-mass shift), but they do not rule out
turnover of {ul:which} units contribute at unchanged counts (identities are not stored), moving columns do not by themselves establish composition bias (that additionally requires effect
heterogeneity across the cohorts entering or leaving), and composition change {ul:within} a coarse {it:i} is visible only through {it:nobs}/{it:sumw_att}, not the count columns. The table
is displayed when {cmd:summarystats} is specified.{p_end}

{synopt:{it:e(escomp_refN)}, {it:e(escomp_refnunits)}, {it:e(escomp_refncoh)}, {it:e(escomp_refW)}} Observation count, distinct treated units, distinct onset cohorts, and total
{it:att}-style weight in the complete reference-normalization pool. Under the pooled regime that pool is the union of treated rows earlier than the requested pre-window
({it:relyr} < -{it:esprelength}) {ul:and} each treated unit's forced first in-sample observation (first two under {cmd:timetrends}) — the forced rows exist even when the requested
window covers every period. Under {cmd:esfixedbaseperiod} it is the {cmd:esrelativeto()} period plus, under {cmd:timetrends}, the immediately preceding trend-anchor period (two relative
periods, not the named reference alone). Available under both regimes.{p_end}

{synopt:{it:e(escomp_trtN)}, {it:e(escomp_trtW)}} Observation count and total {it:att}-style weight over the entire treated estimation sample, for assessing how much treated weight lies
outside the reported event-study window.{p_end}

{synopt:{it:e(anticipation)}} The anticipation shift applied to the treatment onset, returned when {cmd:anticipation()} is active.{p_end}

{synopt:{it:e(roth_power_{it:type})}} With {cmd:pretrendpower()}: one row per requested power level, columns {it:power}, {it:mdt_wald}, {it:mdt_nis} (the minimum detectable linear
pre-trend under the joint-Wald and any-individually-significant pre-tests respectively), and {it:bias_wald}, {it:bias_nis} (the bias such a trend would imply for that block's pooled estimand if it
continued post-treatment, i.e. mdt × (hbar − reference); for {it:att_i}/{it:att_itime} requests the att_it block stands in and its bias columns describe att_it weighting only). Companion scalars: {it:e(roth_slope_{it:type})}, {it:e(roth_slope_se_{it:type})}, {it:e(roth_slope_ll/ul_{it:type})} (the GLS
pre-trend slope and 95% CI), {it:e(roth_prep_{it:type})} (the joint pre-trend test p-value), {it:e(roth_hbar_{it:type})} (the estimand-weighted mean post-treatment relative period over
the full post support of the estimation sample), {it:e(roth_npre_{it:type})}, {it:e(roth_keff_{it:type})} (pre-period count and effective rank), {it:e(roth_W_{it:type})} (the joint
pre-trend test statistic), and {it:e(roth_nissize_{it:type})} (the actual size of the any-individually-significant pre-test). {it:e(roth_ok)} = 1 marks a completed diagnostics
pass: a type absent from the roth_* results was then legitimately skipped or unavailable, and roth_ok = 1 with no {it:e(roth_power_*)} matrix at all means the diagnostics were
requested but none were computable (that caution is replayed). A roth_* object present {ul:without} {it:e(roth_ok)} marks an interrupted pass; replay then prints a caution and
does not render the partial block.{p_end}

{pstd} In addition to these core results, the following materials are also returned in e(): {p_end}
{synopt: Joint Tests by Horizon} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented in the scalars e(joint_horiz)
(ATT), e(joint_horiz_i) (ATT_i), e(joint_Choriz) (AME; continuous treatments only), and e(joint_Choriz_i) (AME_i; continuous treatments only). {p_end}

{synopt: Joint Tests by Subgroup & Time Period} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented by
subgroup, when used, in the local macros: e(joint_horiz_s{it:sgcase}), e(joint_horiz_i_s{it:sgcase}), e(joint_Choriz_s{it:sgcase}), and e(joint_Choriz_i_s{it:sgcase}), where {it:sgcase}
is the value of the subgroup. {p_end}

{synopt: Adversarial Line Tests:} P-values from the composite-null tests that all (non-reference) relative-time period treatment effect estimates fall on the same line
 are found in the scalars {it:e(adversarial_pval_att)}, {it:e(adversarial_pval_ame)}, {it:e(adversarial_pval_att_s[SG])}, etc.; they are posted only under {it:esfixedbaseperiod}
 and only when the event-study covariance is not degenerate (a caution is printed otherwise). The display-oriented fitted line is always returned in matrices named
 {it:e(adversarial_line_att)}, etc., with columns {it:slope} and {it:intercept_rt0}; {it:intercept_rt0} is the fitted value at numeric relative time 0, regardless of wooldid's
 internal nonnegative time coding. {p_end}

{synopt: SG Effect Homogeneity Tests} For each estimate type (ATT, ATT_it, AME, etc.), the joint test that all subgroups have the same such post-treatment effect (namely, the
same effect as the lowest numbered SG) may be found in the scalar {it:e(joint_att_allsgssame)} and analogously named scalars for other treatment effect types. For pre-treatment effects,
the test can be found in objects with names like {it:e(joint_pre_att_allsgssame)}.{p_end}

{synopt:{it:e(N)}, {it:e(N_dropped)}, {it:e(Nclus)}} Estimation sample observation count, number of observations dropped from initially specified sample, and cluster count. {p_end}

{synopt:R^2 measures} For standard estimation, these include e(r2), e(r2adj), e(r2_within), e(r2_withinadj). For poisson estimation, these include e(r2pseud), a pseudo-R^2 measure. {p_end}

{synopt:{it:e(histogramestimates)}} The cohort i-specific treatment effects used to generate histograms, when requested. Row names carry internal temporary-variable identifiers,
not cohort values; the histogram itself is the intended output.{p_end}

{synopt:{it:e(summarystats)}} Summary statistics on the outcome variable and (if present) continuous treatment variable, when requested. Some of the most commonly referenced summary
statistics from this matrix are also returned directly as scalars into e(). {p_end}

{synopt: CFX Plot Info} Information used to generate the continuous treatment effects plots. The information used to produce each plot will be stored in a matrix with the following naming
pattern: e(cfx[item]_[type]), where [item] is either one of att, att_i, att_it, or att_itime (for overall average treatment effects) or one of ame, ame_i, ame_it, or ame_itime (for marginal
effects of the continuous treatment). Meanwhile, [type] is either cfxpost (estimated in the post-treatment period), cfxpre (estimated in the pre-treatment period), or relyrZ (estimated in
relative time period Z, where the letter n represents a negative sign). Matrix rows are keyed {it:lattice_#} (the lattice or {cmd:cfxlattice()} level, in {it:margins} order) with a null
row equation; the same keys apply to the one-column matrices produced under {cmd:cfxplotnose}/{cmd:fastnose}, and AME rows no longer carry the dose-variable equation name. {p_end}

{synopt:{it:e(sample)}} A function characterizing the estimation sample used; useful for identifying which observations are dropped by {it:wooldid} during its sample formation and sanitation
procedures.{p_end}

{synopt:{it:e(cmdline)}} The full {it:wooldid} command line for the successful estimation, beginning with the command name.{p_end}

{synopt:{it:e(estimate_type)}} A macro recording whether your treatment effect estimates use a pooled reference period or a fixed reference period. {p_end}

{synopt:{it:Other Items:}} {it:e(ll)}, {it:e(ll_0)}, {it:e(rmse)}, {it:e(rss)}, {it:e(F)}, and {it:e(chi2)} (poisson only) are passed through directly
from the underlying regression estimator, wherever available. {p_end}

{pstd} In addition to these objects returned in e(), if saveplots() is specified, a dataset containing the information used to produce event study plots along with .gph files containing
all graphs generated will be saved to a file with your specified name in STATA's working directory. {p_end}

{pstd} If the option {it:verbose} is specified, the stage 1 underlying regression that {it:wooldid} estimates will also be stored to a model called {it:___wooldidfullmodel}, which can be
accessed using standard {it:estimates restore} type syntax. Note that {it:___wooldidfullmodel} is a reserved scratch estimate name: {it:wooldid} overwrites and
removes a stored estimate carrying this name on every run, including runs that end in an error (under {it:verbose}, {it:___wooldidfullmodel} is retained on failed runs as well, for
inspection). Do not store your own estimates under this name. Under {it:verbose} with {it:poisson} (in the high-dimensional-fixed-effects estimation path), the {it:_ppmlhdfe_d} scaffolding
variable that {it:ppmlhdfe} creates is also retained in the dataset rather than dropped at exit. Note that coefficient-level replay of the stored model ({it:lincom}, {it:test},
{it:matrix list e(b)}) works, but data-touching replay ({it:predict}, {it:margins}) does not: the model's regressors are internally constructed temporary variables that no longer exist
after {it:wooldid} returns. {p_end}

{pstd} {ul:Tips for accessing output and remaking graphs:} Use the {it:matlist} command to print the contents of one of {it:wooldid}'s output matrices. To convert one of {it:wooldid}'s 
output matrices into a dataset that you can directly manipulate (e.g., to remake the automatically generated event study graphs), use the {it:svmat} command. For example, the command 
{it:svmat e(histogramestimates), names(col)} will create in memory a dataset containing the contents of {it:e(histogramestimates)} (i.e., the cohort i-specific treatment effects), 
using as variable names the names of the matrix's columns. The program {help event_plot} can also be used to automatically generate event study plots using the information stored by
{it:wooldid} in {it:e()}, though the user will need to point {it:event_plot} to the name of the results of interest (e.g., for ATT-type event study coefficients, by specifying 
{it:stub_lag(ES_att_post#)} and {it:stub_lead(ES_att_pre#)}). If taking this approach, be aware that the {it:event_plot} does not automatically insert 0s into the event study plot 
to represent the missing reference period. {p_end}

{pstd} {ul:Compatible Post-Estimation Commands:} {it:estimates save}, {it:estimates store}, {it:estimates restore}, {it:esttab}, and similar commands. {p_end}

{pstd} {ul:Replay:} Typing {it:wooldid} by itself after a successful estimation redisplays the full stored results from {it:e()} - the header, summary statistics, all results tables,
the joint, adversarial, and pre-trend power diagnostics, and the event-study composition, cohort-by-period, and randomization-inference blocks when present - byte-for-byte as originally shown. Replay works after {it:estimates restore} of a stored {it:wooldid} model (subject to the failed-call
refusal below), after {it:cleanmatrices}, and with no dataset in memory. It does not re-draw plots, does not repeat estimation-time Caution messages (which describe the estimation run, not the results),
and accepts no options: any non-empty command line - including {it:wooldid, update} - is processed as a normal call. If a {it:wooldid} call fails anywhere in its main pass - parsing, estimation, diagnostics, display, or plots - {it:e()} still holds the previous successful estimation
(standard Stata behavior), but replay then refuses with an informative error until a {it:wooldid} call completes successfully: redisplaying results that predate the failed call would
be misleading. The refusal is session-scoped (an internal {it:WOOLDID_LASTRC} marker; {it:WOOLDID_*} global names are reserved by {it:wooldid}): it also blocks replay of models
explicitly selected with {it:estimates restore} until some {it:wooldid} call completes. If a failure instead occurs in the separate post-estimation randomization-inference phase of a
{cmd:permute()} call, that call's own completed ordinary estimation is restored, valid, and remains replayable (the permutation request itself failed loudly). Replay requires results posted by a
{it:wooldid} version with replay support; internally, the display state travels in hidden {it:e(wdid_*)} macros, which are implementation details and not part of the documented
{it:e()} surface. {p_end}

{pstd} {ul:Incompatible Post-Estimation Commands:} {it:lincom}, {it:test}, and all other post-estimation commands that rely on receiving a valid variance-covariance matrix from
an estimator. While {it:wooldid} does return a coefficient vector ({it:e(b)}) and a variance-covariance matrix ({it:e(V)}), these are for display purposes only. The variance-covariance
matrix has all off-diagonal elements set to 0, meaning that while {it:esttab} and the like will present correct estimates using these inputs, efforts to combine or otherwise test different
estimates will not, in general, yield valid results. This limitation applies to the summary {it:e(b)}/{it:e(V)}; the event-study coefficients are additionally returned with a genuine full variance-covariance matrix in {it:e(b_es_{it:type})} / {it:e(V_es_{it:type})} (see Stored Results), and {it:lincom} and {it:test} {ul:are} valid when applied to those (inspect the diagonal first when {it:e(es_ndegen_{it:type})} is present: structurally degenerate rows carry zero or missing variance). {it:honestdid} is valid only under the positional-geometry conditions listed in the {it:e(b_es)} Stored Results entry. For {it:lincom}/{it:test}, first repost the event study with {it:ereturn post e(b_es_{it:type}) e(V_es_{it:type})}. {p_end}



{marker comparison}{...}
{title:Comparison to Related Estimators}

{pstd} {ul:Comparing Output to hdidregress twfe:}
The STATA 18 command {help hdidregress} has an option for implementing Wooldridge (2025) style estimates. Its implementation takes the Mundlak approach to estimation
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
The Borusyak, Jaravel, and Spiess (2023) difference-in-differences imputation estimator is implemented in the package {help did_imputation}. Wooldridge (2025)
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

{phang}7) Event study with a fixed reference period, treating treatment as anticipated one period before its recorded onset, with Roth (2022) pre-trend power diagnostics
at 50 and 80 percent power:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)}

{phang}8) Difference-in-differences with a treatment intensity variable (one that is 0 outside the treatment period), estimated using a stepwise-linear function and a lattice with 5 quantiles
(formed based on the treated group's treatment period), and generating continuous treatment effects plots for the post treatment period:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(intensity) contreatcontrols(untreatedZero) makecfxplot lattice(baseTreatedPost 5)}

{phang}9) Difference-in-differences with a continuous treatment variable that varies within all cohorts, with a lattice containing 5 quantiles (formed separately by treatment group and
by period), controlling for the continuous treatment variable by treatment group and by time period, with continuous treatment effects plots for the pre-treatment period (which
require an event-study window - here via {it:esprelength}):{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(contreatvar) contreatcontrols(byTreated byTime) espre(3) espost(3) makecfxplot lattice(byTreatedByTime 5) cfxplottypes(pre)}

{phang}10) Event study followed by a Rambachan & Roth (2023) HonestDiD pre-trends sensitivity analysis, using the full event-study variance-covariance matrix wooldid returns. (Requires the {stata ssc install honestdid:honestdid} package. A fixed reference period at the default {it:esrelativeto(-1)} is {ul:required} for this positional handoff — the example relies on the reference being the single omitted period at relative time -1; see the {it:e(b_es)} Stored Results entry for the full geometry conditions.):{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att espre(4) espost(4) esfixedbaseperiod}

{p 10 10}   {cmd:matrix b = e(b_es_att)}

{p 10 10}   {cmd:matrix V = e(V_es_att)}

{p 10 10}   {cmd:honestdid, numpre(`=e(numpre_att)') b(b) vcov(V) delta(rm) mvec(0.5(0.5)2)}

{phang}11) Fast exploratory event study (point estimates only), also returning the cohort-by-relative-period building-block table:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) espre(2) espost(2) fastnose gteffects}

{phang}12) Randomization inference: permute treatment-onset dates across clusters, holding onset-group counts fixed, with 500 replications (set the seed for reproducibility):{p_end}

{p 10 10}   {cmd:set seed 12345}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(500)}

{phang}13) Cluster bootstrap standard errors with 500 replications, reproducible, and joint pre-trend tests computed on the bootstrap covariance - the sanctioned route when clustering at or above the cohort level:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) jointtests bootstrap(standard) reps(500) rseed(42)}

{phang}14) Bayesian (Dirichlet reweighting) bootstrap, which leaves every cluster in every replication and so never loses the comparison group to an unlucky resample:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) bootstrap(bayes) reps(500) rseed(42)}



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

{p 10 10}   Under various circumstances, {it:wooldid} cannot estimate standard errors for all requested objects (for example, when an Ibragimov-Muller estimate at some event-study horizon
has only a single contributing cluster, or when per-cluster estimates are numerically identical on effectively noiseless data). If standard errors are identified for some but not all
reported estimates, {it:wooldid} still returns {it:e(V)}: the identified estimates carry their variances as usual, while each affected estimate carries a placeholder entry that displays
as a missing standard error (and that {it:test} drops and {it:lincom} propagates as missing). The placeholder is stored as a large negative number, which every Stata display and inference
surface renders as missing - but code that reads {it:e(V)} cells directly for its own matrix algebra should treat negative diagonal entries as not identified. If no reported estimate has
an identified standard error, {it:wooldid} returns point estimates without {it:e(V)}, which can cause the functionality of certain post-estimation tools to break. In these situations,
{it:wooldid} prints a caution; when the cause is unidentified standard errors, the caution states how many estimates were affected.

{pstd}If you encounter a bug - be it the program crashing or just the program delivering a weird result - please reach out to me about it. This program is still
under development, so assistance finding and correcting bugs is much appreciated. I am also happy to help and offer guidance on how to use this program.{p_end}

{pstd} Note: if you encounter unusual, hard to track down errors, try uninstalling and reinstalling {it:reghdfe} and {it:ftools} first. This can fix some problems. {p_end}



{marker citations}{...}

{title:Citations}

{phang} Borusyak, Kirill, Xavier Jaravel, and Jann Spiess. Revisiting Event Study Designs: Robust and Efficient Estimation. Working Paper, 2023. {p_end}

{phang} Rambachan, Ashesh and Jonathan Roth. A More Credible Approach to Parallel Trends. Review of Economic Studies, 2023, 90(5): 2555-2591. {p_end}

{phang} Roth, Jonathan. Pretest with Caution: Event-Study Estimates after Testing for Parallel Trends. American Economic Review: Insights, 2022, 4(3): 305-322. {p_end}

{phang} Correia, Sergio. reghdfe: Stata module for linear and instrumental-variable/gmm regression absorbing multiple levels of fixed effects. Statistical Software Components s457874, 2017.
Boston College Department of Economics. https://ideas.repec.org/c/boc/bocode/s457874.html {p_end}

{phang} Ibragimov, Rustam and Muller, Ulrich K. t-Statistic Based Correlation and Heterogeneity Robust Inference. Journal of Business & Economic Statistics, 2010. {p_end}

{phang} Wooldridge, Jeffrey M. Two-Way Fixed Effects, the Two-Way Mundlak Regression, and Difference-in-differences Estimators. Empirical Economics, 2025. {p_end}



{marker author}{...}
{title:Author & How to Cite}

{pstd} Author: Thomas A. Hegland (Agency for Healthcare Research and Quality)  {p_end}

{pstd} Contact: thomas.hegland@ahrq.hhs.gov, thomashegland.com, @thomas_hegland  {p_end}

{pstd} Citation: Hegland, Thomas A. wooldid: Estimation of Difference-in-Differences Treatment Effects with Staggered Treatment Onset Using Heterogeneity-Robust Two-Way
Fixed Effects Regressions. Statistical Software Components s459238, 2023.{p_end}



{marker disclaimer}{...}
{title:Disclaimer}
{pstd}This document reflects only the views of the author and does not necessarily represent the views of the Agency for Healthcare Research and Quality, the Department
of Health and Human Services, or the United States government. The program with which this help file is associated is the work of the author and does not come
with any endorsement by or warranty from the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government.  {p_end}
