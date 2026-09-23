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
{viewerjumpto "Bootstrap Standard Errors" "wooldid##bootstrap"}{...}
{viewerjumpto "Randomization Inference" "wooldid##permute"}{...}
{viewerjumpto "Supporting Features" "wooldid##supporting"}{...}
{viewerjumpto "Event Studies" "wooldid##eventstudy"}{...}
{viewerjumpto "Subgroup Effects" "wooldid##interactions"}{...}
{viewerjumpto "Continuous Treatments" "wooldid##continuous"}{...}
{viewerjumpto "Controls and Fixed Effects" "wooldid##controls"}{...}
{viewerjumpto "Technical Features" "wooldid##technical"}{...}
{viewerjumpto "Runtime" "wooldid##runtime"}{...}
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
across subgroups, event study style estimates and related parallel trends assumption tests, and average and marginal treatment effects in the presence of a
continuous treatment (or treatment intensity) variable [experimental]. It also produces event study plots, histograms of treatment effects by cohort, and plots of treatment
effects across the distribution of a specified continuous treatment variable. Options for inference include standard analytic clustered standard errors,
unconditional standard errors that account for sampling variation in covariates,
Ibragimov and Muller (2010) style cluster-robust inference, bootstrap standard errors, and randomization inference for settings with few treated clusters.
Other features include anticipation windows ({it:anticipation()}), a table of cohort-by-period effects ({it:gteffects}), pre-trend power diagnostics ({it:pretrendpower()}),
a point-estimates-only mode ({it:fastnose}), and replay of the stored results by typing {it:wooldid} alone after estimation. {p_end}

{pstd}{ul:Shortcuts:} {help wooldid##description:Description}, {help wooldid##options:Options}, {help wooldid##output:Working with Output},
{help wooldid##runtime:Runtime}, {help wooldid##comparison:Comparison with Related Estimators}, {help wooldid##examples:Examples},
{help wooldid##bugs:Bugs and Assistance}, {help wooldid##citations:Citations}, {help wooldid##author:Author Info and How to Cite}, {help wooldid##disclaimer:Disclaimer} {p_end}

{pstd}{ul:Dependencies:} {cmd:wooldid} requires Stata 16 or higher, {help reghdfe} version 5.7.3, {help ftools},
and {help gtools}. Poisson estimation also requires {help ppmlhdfe}. {p_end}

{pstd}{ul:Github Link:} {browse "https://github.com/thegland/wooldid/"} - the command {it:wooldid, update} updates this program from the materials at this link. A longer usage guide written
for AI agents and advanced users is at {browse "https://github.com/thegland/wooldid/blob/main/wooldidinstall/wooldid_agent_guide.md"}; it is also an ancillary file of the package, copied to
the working directory by {cmd:net get wooldid}. {p_end}

{marker syntax}{...}
{title:Syntax}

{phang}{cmd: wooldid} {it:y i t ttre} [if] [aw/pw=weights] [{cmd:,} {help wooldid##specifyestimand:estimands} {help wooldid##standarderrors:standard_errors}
{help wooldid##bootstrap:bootstrap_standard_errors} {help wooldid##permute:randomization_inference} {help wooldid##supporting:supporting_features} {help wooldid##eventstudy:event_study}
{help wooldid##interactions:subgroup_effects} {help wooldid##continuous:continuous_treatment} {help wooldid##controls:controls_and_FE} {help wooldid##technical:technical_features} ]{p_end}

{phang}{ul:The if qualifier and weights:} {cmd:wooldid} keeps the observations for which the {it:if} expression equals exactly 1, where Stata's own {it:if} keeps any nonzero value; write {cmd:if flag != 0} for Stata's behaviour.
Weights may be aweights or pweights; frequency weights (fw) are not accepted. {p_end}

{synoptset 22 tabbed}{...}
{synopt : {it:y}}Outcome variable. {p_end}
{synopt : {it:i}}Cohort identifier; must consist of integers greater than 0. Cohorts can be fine (i = individual units) or coarse (i = collections of individual units), but should be no coarser
than to consist of treatment onset cohorts. Coarser cohort definitions tend to enable much faster estimation, so do not default {it:i} to the finest level: coarsen it and specify unit-level
fixed effects in {opt fe()}, which a coarse {it:i} also needs when the panel is unbalanced within cohorts. {p_end}
{synopt : {it:t}}Time period identifier; must be numeric and on the same grid as {it:ttre}. Zero, negative and non-consecutive values are fine on the default path. With {opt uncond:itionalse},
{opt customvce()} or {opt inter:activecontrols()}, {it:t} enters the model as a Stata factor variable and must then consist of non-negative integers. {p_end}
{synopt : {it:ttre}}Variable containing cohort-specific dates of treatment onset; missing indicates never-treated cohorts, cohorts treated before the sample begins are dropped. {it:ttre} must be
constant within {it:i} and on the same grid as {it:t} ({it:t} minus {it:ttre} an integer), and some cohort must be never-treated or first treated after the estimation window; each is checked
with an informative error unless {opt safety(off)}. Store large onset codes, such as 8-digit dates, as double: float rounds them and silently merges adjacent cohorts. {p_end}

{phang} {ul:Options by Purpose:} {p_end}
{synoptset 1 tabbed }{...}
{synopt : {help wooldid##specifyestimand:Specifying the Estimand(s):}} {opt att} {opt att_it} {opt att_i} {opt att_itime} {p_end}

{synopt : {help wooldid##standarderrors:Standard Errors:}} {opt clus:ter(variable(s)[numerical])}  {opt rob:ust} {opt uncond:itionalse} {opt customvce(customstring)} {opt imp:values} {opt level(#)}
{opt fast:nose} {opt nofastcollapse} {opt simult:aneousbands} {p_end}

{synopt : {help wooldid##bootstrap:Bootstrap Standard Errors:}} {opt boot:strap(standard|bayes)} {opt reps(int)} {opt bootts:tat} {opt bootdisp:lay(spec)} {opt bootd:f(spec)} {opt bootalp:ha(spec)}
{opt ties} {opt bootsa:ve(name)} {opt rseed(#)} {opt dots(int)} {p_end}

{synopt : {help wooldid##permute:Randomization Inference:}} {opt perm:ute(shuffle|probabilities)} {opt permreps(int)} {opt permstr:ata(variable)} {opt permpr:obs(spec)} {opt permex:act}
{opt permts:tat} {opt permt:ci} {opt perman:alyticse} {opt permab:sp} {opt permsi:de(left|right)} {opt permraw:p} {opt permstep:down} {opt permonset:cohorts} {opt permsa:ve(name)} {p_end}

{synopt : {help wooldid##supporting:Supporting Features:}} {opt hist:ogramcohorteffects} {opt save:plots(filename_stem)} {opt replace} {opt sum:marystats}
{opt pois:son} {p_end}

{synopt : {help wooldid##eventstudy:Event Study Estimation:}}  {opt esf:ixedbaseperiod} {opt espre:length(int>=0)} {opt espost:length(int>=0)} {opt esr:elativeto(-int)}
{opt joint:tests}   {opt adv:ersarial} {opt altmethod_adversarial} {opt makep:lots} {opt oos:treatedcontrols(int>=0)} {opt antic:ipation(int>=0)} {opt pretrendpow:er(numlist)}  {p_end}

{synopt : {help wooldid##interactions:Subgroup Effects:}} {opt sub:groups(variable[>=0 int])} {opt suppress:sgprimaryeffects}  {p_end}

{synopt : {help wooldid##continuous:Continuous Treatments:}} {opt contreat(variable[numerical])} {opt contreatelas:ticitytype(dydx|dyex)}  {opt contreatw:ithin}
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

{pstd} The exact margins calculation performed in step 2 depends on the estimands
specified by the user. The overall average treatment effect would be obtained by computing fitted values from the underlying regression among all treated observations (those where
  TreatedFlag_it is equal to 1), subtracting the counterfactual fitted values for these observations that set TreatedFlag_it to 0, and then taking the average of these differences
  (applying any sample weights). To target estimands other than the overall average treatment effect, this procedure can be modified by either
  introducing a set of weights (or modifying the pre-existing sample weights) when averaging the computed differences in the margins calculation or by restricting the sample so that the
  average only contains observations from a subgroup or time period of interest.
For more on the different estimand-types available, see {help wooldid##specifyestimand:estimands}. {p_end}
{pstd} With respect to inference, there are several approaches to standard error estimation available. The default estimates standard errors in the underlying regression and then obtains
standard errors on the treatment effects produced by the margins calculation using the delta method. These errors are conditional standard errors in the sense that they are conditional
on the data: the difference in fitted values potentially depends on the value of various covariates, and this approach takes these values as fixed when performing the margins calculation.
P-values are then obtained from a t-distribution with the residual degrees of freedom of the underlying regression.
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
at varying points along the distribution of the continuous treatment variable. It can also report the effect of the continuous treatment variable with respect to a proportional change in the dose
rather than a unit change. For more detail on options related to continuous treatments, see {help wooldid##continuous:continuous_treatment}. {p_end}

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
and thus that these cohorts should be allowed to contribute to pre-treatment relative time period-specific estimates in the event study context. The option {it:anticipation(#)}
treats the effective onset as {it:ttre} - #, so that periods of anticipatory response count as treated rather than pre-treatment periods.
{p_end}

{marker options}{...}
{title:Options}

{marker specifyestimand}{...}
{dlgtab:Specifying the Estimand(s)}

{phang}{opt att}: Estimate average treatment on the treated (ATT) effects; this option is the default and may be combined with options att_it, att_i, and att_itime to estimate multiple different
effect types. ATT effects are computed by performing the margins calculation (see step 2 in {help wooldid##description:Description} for more on the margins calculation)
using whatever weights (if any) are used in the main regression command, multiplied by the {it:customweightmultiplier()} variable if one is specified; that multiplier enters the margins weights only, not the regression weights.{p_end}

{phang}{opt att_it}: Estimate simple averages of the i-t cell-specific treatment effects (ATT_it). The treatment effect in each cohort-period is calculated via the margins calculation
using any weights that may apply to the observations falling in each i-t cell. The final ATT_it estimate is the simple, unweighted mean of all of the cohort-period-specific estimates.
This is done by reweighting the data used in the margins calculation (after estimation of the underlying regression) by w/W_it where w is an observation's sample weight and W_it is the
total sample weight associated with all observations appearing in a given i-t cell. If no weights are specified, w is taken to have a value of 1 for all observations. If cohorts are
fine and thus each i-t cell contains only one observation, this is equivalent to ignoring any specified sample weights when performing margins calculations.{p_end}

{phang}{opt att_i}: Estimate the simple average of cohort-specific treatment effects (ATT_i). The final ATT_i estimate is the simple, unweighted mean of the cohort-specific treatment effects.
ATT_i effects are computed by performing the margins calculation, but reweighting observations by either w/W_ipre or w/W_ipost, where w is an observation's sample weight, W_ipre is the
total sample weight associated with cohort-i in the pre-treatment period (but outside the reference period), and W_ipost is the total sample weight associated with cohort-i in the
post-treatment period. Under {it:timetrends} the period immediately before the reference period is also normalized to zero and excluded from W_ipre. Each observation is then
reweighted based on whether it falls in the pre- or post-treatment period. If weights are not used, w is taken to have a value of 1 for all observations.{p_end}

{phang}{opt att_itime}: Estimate the simple average of cohort-specific simple average treatment effects (ATT_itime). Treatment effect calculation here begins like for {it:att_it}, obtaining
treatment effects for individual cohort-periods of interest using weights. Then, cohort-specific effects are obtained by taking the simple mean of the cohort-period-specific effects.
The overall presented ATT_itime estimate is then equal to the simple, unweighted mean of these cohort-specific effects. ATT_itime effects are calculated by performing the margins calculation,
but reweighting observations by w/(W_it * T_iprepost), where w is an observation's sample weight, W_it is the amount of sample weight falling in a given i-t cell, and T_iprepost contains
(depending on the time period) either the number of post-treatment time periods in which a cohort-i is observed or the number of pre-treatment, non-reference period time periods in which a
cohort i is observed. If weights are not used, w is taken to have a value of 1 for all observations.{p_end}

{pstd}Under {it:contreat()} each requested estimand is also reported as the average marginal effect of the continuous treatment, labelled {it:ame}, {it:ame_it}, {it:ame_i} and {it:ame_itime}.
Because {it:att_it}, {it:att_i} and {it:att_itime} weight cells or cohorts equally rather than observations, some options treat them differently. {it:unconditionalse} and
{it:bootstrap(bayes)} are refused for them unless {it:safety(off)} is specified. {it:makecfxplots} produces its plots and the {it:e(cfx*)} matrices for {it:att} only unless
{it:safety(off)} is specified. A {it:subgroups()} variable that splits a unit across subgroups over time is refused; under {it:safety(off)} such a unit enters each subgroup by its share of
sample weight, so {it:att_i} and {it:att_itime} then no longer weight units equally.{p_end}

{marker standarderrors}{...}
{dlgtab:Standard Errors}

{phang}{opt clus:ter(clusvars)}: Variable(s) specifying clustering level. Clustered standard errors are computed analytically. By default,
regression estimation using this option will rely on either {help reghdfe} or, for poisson, {help ppmlhdfe}. Several variables give multi-way clustering on
this default path; {it:unconditionalse}, {it:impvalues}, {it:bootstrap()} and {it:permute()} each require a single cluster variable. A cluster variable that leaves
fewer than two clusters in the estimation sample stops the run with an error; a finer cluster variable, {it:robust} or {it:fastnose} are the ways out. {p_end}

{phang}{opt rob:ust}: Heteroskedasticity robust standard errors; relevant only if cluster() is not specified. By default, regression estimation using this
option will rely on either {it:reghdfe} or, for poisson, {it:ppmlhdfe}. {p_end}

{phang}{opt uncond:itionalse}: Compute standard errors that account for sampling variation in covariates. This option also forces {it:wooldid} to abandon
high dimensional fixed effect estimation packages for the underlying two-way fixed effect regression in favor of the commands {it:regress} and {it:poisson}.
As such, this option should be slower to use on average. {it:fe()} and {it:ccc_absorb()} cannot be absorbed on this path: {it:wooldid} cautions and enters them as written, so put them in {it:controls()} with your own
{it:i.} prefixes. This option requires specifying either robust or clustered standard errors, with a single cluster variable; a robust or clustered vce supplied
through {it:customvce()}, or pweights, also satisfies this requirement. Unless {it:safety(off)} is specified, it is refused with {it:att_i}, {it:att_it} and
{it:att_itime} and cautions with aweights. It is refused with {it:impvalues}, {it:fastnose}, {it:bootstrap()} and {it:permute()}, and it turns the default cell
collapse off. {p_end}

{phang}{opt customvce()}: Forces {it:wooldid} to switch to using either the {help regress} or {help poisson} command for underlying regression estimation,
and then passes the user specified text without modification into the {it:vce()} option for these commands. This allows specifying a broader set of
standard error calculation methods than allowed by the default high dimensional fixed effects estimation packages used by {it:wooldid}. Not compatible with
{it:robust} or {it:cluster()}: specify the clustering or robustness desired inside the {it:vce()} text itself. Like {it:unconditionalse}, it cannot absorb {it:fe()}
or {it:ccc_absorb()}, which are entered as written under a caution. Refused with {it:impvalues} and {it:bootstrap()}, and it cannot be combined with {it:permute()}, which needs a {it:cluster()} variable to
randomize over. It turns the default cell collapse off. {p_end}

{phang}{opt imp:values}: Implements Ibragimov and Muller (2010) style cluster-robust inference for all estimates. The idea behind the Ibragimov and Muller approach is to compute
cluster-specific treatment effect estimates, present their average as one's final treatment effect estimate, and conduct inference by t-testing the cluster-specific estimates against 0.
Using IM-style inference requires specifying a single cluster variable at the level of {it:i} or coarser. A {it:cluster()} finer than {it:i}, or one that cuts across it, is
refused unless {it:safety(off)} is specified, since the cluster-specific estimates are then re-weightings of the same coefficients and the reported standard error is far too small.
IM-style inference changes {it:wooldid}'s target estimand: the cluster-specific estimates are averaged with equal weight, so what is reported is the cluster-equally-weighted mean
whatever estimand was requested. For example, if clustering at the cohort level, {it:att} will yield the same effect estimates as {it:att_i}. Each reported row uses a t
distribution with degrees of freedom one less than the number of clusters contributing to that row, and {it:e(df_r)} is set to the smallest of these. IM-style inference is not
allowed with the use of subgroups, or with the options {it:jointtests}, {it:makecfxplots}, {it:adversarial}, {it:customvce()}, {it:unconditionalse}, {it:fastnose},
{it:pretrendpower()}, {it:bootstrap()} and {it:permute()}. No event study covariance is returned under {it:impvalues}, so {it:e(V_es_*)} is absent and the {it:honestdid} handoff is unavailable. {p_end}

{phang}{opt fast:nose}: Fast exploratory mode, for specification search on large panels. Every {it:margins} call runs without variance estimation, so only point estimates
are produced: standard errors, test statistics, p-values, confidence intervals, joint tests and plot confidence bands are omitted, and the corresponding table cells display
as missing. {it:e(V)} is not posted, so {it:test} and {it:lincom} will refuse, while {it:esttab} prints an empty inference row. Refused with {it:unconditionalse},
{it:impvalues}, {it:pretrendpower()}, {it:bootstrap()}, {it:permtstat}, {it:permtci}, {it:adversarial} and {it:jointtests}; the last two are accepted under
{it:permute()}, where fast mode is already the default. {it:simultaneousbands} is accepted but no band is built. Re-run the chosen specification without {it:fastnose}
before reporting results. {p_end}

{phang}{opt nofastcollapse}: Turns off the default cell collapse, which exists to significantly speed up {it:wooldid}: the {it:margins} calls are evaluated on one row per
cohort-by-period-by-treatment cell rather than on every observation, while the regression itself still runs on the full data. Results agree with an uncollapsed run to
machine precision but are not bit-identical. It applies to binary and continuous treatments, to poisson, to subgroups and inside {it:bootstrap()} and {it:permute()}
replications. It is skipped, with a caution naming the option responsible, when the request is ineligible: {it:interactivecontrols()}, {it:contreatpoly()} above 1,
{it:lattice()}, {it:makecfxplots}, the {it:contreat()} variable also appearing inside {it:controls()}, {it:coarsecohortcontrols()}, {it:fe()} or {it:ccc_absorb()},
{it:unconditionalse}, and {it:customvce()}. Disable it yourself when you are using controls or fixed effects interacted with treatment in some way {it:wooldid} cannot
detect, such as a control built from a renamed copy of the treatment variable, since the collapse assumes the treatment effect is constant within a cell; when your
analysis is memory constrained, since the collapse preserves the data and so holds an extra copy of it; and when you need byte-for-byte reproduction of an uncollapsed
result. {p_end}

{phang}{opt simult:aneousbands}: Compute a simultaneous (sup-t) confidence band over each event study's horizons at {it:level()}, which covers every horizon of the
family at once with that probability. The printed table keeps its pointwise intervals; the band endpoints are stored in {it:e(esband_*)} and drawn as a dashed pair on
the event study plots. Without {it:bootstrap()} the band is simulated from the family's analytic event study covariance; with {it:bootstrap()} it is taken from the
replications. Refused with {it:permute()}, whose multiplicity tool is {it:permstepdown} instead. No band is built under {it:fastnose} or for a family with fewer than
two horizons with a usable variance, and a caution says so. {p_end}

{phang}{opt level(#)}: Confidence level for every reported confidence interval, including the plots, the {it:gteffects} cells and the {it:bootstrap()} and
{it:permute()} intervals. It must lie between 10 and 99.99 with at most two digits after the decimal point; the default is {it:c(level)}, normally 95. The interval
columns of the stored results are labelled {it:lb_#ci} and {it:ub_#ci} for the level in force, and the level is stored in {it:e(level)}. The pre-test sizes reported
by {it:pretrendpower()} are fixed at 5%. {p_end}

{pstd} {ul: Discussion of What Level to Cluster At:}
With {it:cluster()} at the level of {it:i} or coarser, each treated cohort-by-period cell lies wholly inside one cluster, and the saturated cell dummies leave the
cluster-robust covariance nothing within the cluster to identify. Several consequences are specific to {it:wooldid}: the analytic standard errors of the continuous treatment
({it:ame}) estimates are not identified and are reported as missing; the analytic joint pre-trend tests and the adversarial line test are not calibrated; pairwise
subgroup contrasts can be degenerate; and the {it:att} rows' analytic standard errors are understated without a caution, the {it:gteffects} cells more so (a {it:gteffects} run counts the
affected cells and cautions). The remedies are a {it:cluster()} strictly finer than the cohort-by-period cells, {it:bootstrap()}, or {it:permute()}. {p_end}

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
units and where time fixed effects are not shared across clusters. Note that {it:wooldid} reads each cluster's estimate over its treated post-treatment cells alone, so control units inside
a cluster do not give it its own control arm unless, for example, cluster-by-time effects are placed in {it:fe()}. In the standard diff-in-diff setting, {it:wooldid} will still present
IM-style estimates by computing cluster-specific treatment
effect estimates, despite the fact that the existence of a shared control group (and perhaps a shared reference period, if {it:esfixedbaseperiod} is not being used) likely causes mischief with the
independence assumption. {p_end}

{pstd} Replication-based alternatives to the analytic family are described under {help wooldid##bootstrap:Bootstrap Standard Errors} and {help wooldid##permute:Randomization Inference}. {p_end}

{marker bootstrap}{...}
{dlgtab:Bootstrap Standard Errors}

{phang}{opt boot:strap(mode)}: Bootstrap standard errors and inference: the reported standard errors, test statistics, p-values, confidence intervals, joint tests
and subgroup contrasts become bootstrap quantities. Point estimates are unchanged; the analytic covariance is kept in the hidden {it:e(V_analytic)}. Each
replication re-estimates the full model. {it:bootstrap(standard)} is the pairs cluster bootstrap: the {it:cluster()} groups, or
the observations when {it:cluster()} is not given, are resampled with replacement. {it:bootstrap(bayes)} is the Bayesian bootstrap: each replication draws one
Dirichlet weight per cluster and multiplies that cluster's observations by it, so no cluster ever leaves the sample. Requires {it:reps()}. Refused with {it:permute()}, {it:impvalues},
{it:unconditionalse}, {it:customvce()}, {it:fastnose}, {it:oldsyntax},
{it:suppresssgprimaryeffects} and multi-way {it:cluster()}. {it:bootstrap(bayes)} with {it:att_it}, {it:att_i} or {it:att_itime} omits the cell-sampling component
of the variance and is refused unless {it:safety(off)} is given, which makes it a caution. Plots, {it:gteffects} cells, {it:jointtests}
and {it:pretrendpower()} all follow the bootstrap. Only the {it:histogramcohorteffects} output is untouched: it draws the observed point estimates, which the bootstrap does not change. Under
{it:cfxplotnose} the continuous treatment cells receive no bootstrap inference, so drop that option if bands are wanted. {it:wooldid} cautions when fewer than 6
clusters are treated. Neither mode is a remedy for few treated clusters and no wild cluster bootstrap is offered; use {it:permute()} there, or
{help boottest} on the underlying regression outside {it:wooldid}. When {it:cluster()} does not nest {it:i} and some unit spans two clusters, a caution says the
standard errors of the cell-weighted estimands may be understated; the remedies are a finer {it:i} or a {it:cluster()} that nests {it:i}. {p_end}

{phang}{opt reps(#)}: Number of bootstrap replications. Required with {it:bootstrap()} and refused without it. No minimum is imposed, so a small value can be used to
time a specification. Replications that fail to estimate are excluded and counted in {it:e(boot_failed)}. The smallest p-value a run can print is set by the number
of completed replications {it:n}: 2/{it:n} for the percentile and bias-corrected coefficient p-values, which can print exactly 0, 2/({it:n}+1) under
{it:boottstat}, and 1/({it:n}+1) for the joint tests, and a caution names the limit whenever it exceeds .05. {p_end}

{phang}{opt bootalp:ha(spec)}: The Dirichlet convention for {it:bootstrap(bayes)}; refused with {it:bootstrap(standard)}. {it:bootalpha(equal)}, the default, draws
Dirichlet(1,...,1) across clusters and applies each cluster's draw as a multiplier on its observations, so a cluster's total mass stays proportional to its size.
{it:bootalpha(size)} instead gives each cluster a Gamma shape equal to its observation count and spreads the drawn mass equally over its rows, which matches
{help bayesboot} with {it:id()}. It is for comparison against that command only; on panel data it understates the standard errors badly. {p_end}

{phang}{opt bootd:f(spec)}: The reference distribution for the Wald family of bootstrap statistics; refused without {it:bootstrap()}. {it:bootdf(normal)} uses standard
normal critical values and reports z statistics. {it:bootdf(clusters)} uses the t distribution with {it:G}-1 degrees of freedom, {it:G} the number of
clusters in the estimation sample, and is refused without {it:cluster()}; a positive integer sets the degrees of freedom directly. The default is
{it:bootdf(clusters)} when {it:cluster()} is given and {it:bootdf(normal)} otherwise. Only the Wald family p-values and intervals change; the bootstrap standard
errors, the other families and the joint tests do not, and under {it:boottstat} nothing printed does. {it:e(df_r)} follows the same choice, so a {it:test} or
{it:lincom} run afterwards on {it:e(V)} reproduces the Wald family, which under the default {it:bootdisplay(percentile)} is not the printed p-value; the footer
says so. {p_end}

{phang}{opt bootdisp:lay(spec)}: Which p-value and confidence interval the printed tables carry; refused without {it:bootstrap()}. The four families are
{it:percentile}, the default, from the percentile ranks of the draws; {it:bc}, the bias-corrected percentile method; {it:bca}, the bias-corrected
and accelerated method, whose acceleration comes from a delete-one-cluster jackknife of {it:G} extra fits and which requires {it:cluster()}; and {it:wald}, a normal
or t approximation governed by {it:bootdf()}. Every family except {it:bca} is computed and stored in {it:e()} on every run; this option chooses only what is displayed
and never changes the reported standard error, which is the standard deviation of the bootstrap distribution. The percentile and bias-corrected intervals are not
symmetric about the estimate, and the statistic column beside them, labelled {it:t_wald} or {it:z_wald}, is kept for reference only. Those two p-values can print
exactly 0, meaning no completed replication fell beyond the null; a line under the table gives the resolution. Under
{it:bootdisplay(bca)} the {it:gteffects} and continuous treatment cells carry the {it:bc} family. With {it:boottstat} only {it:percentile} is available. {p_end}

{phang}{opt bootts:tat}: A studentized, or percentile-t, bootstrap; refused without {it:bootstrap()}. The object bootstrapped is each replication's t statistic on
its own analytic standard error, rather than the coefficient. The reported standard error is therefore the analytic one, labelled {it:se_analytic}, the statistic
is labelled {it:t_pt}, and the p-value and confidence interval come from the ranks of the studentized draws. The interval is not symmetric about the estimate.
Studentization covers every coefficient, the {it:gteffects} and continuous treatment cells and the subgroup contrasts, and the joint pre-trend, subgroup equality
and {it:adversarial} tests are bootstrapped as test statistics; without {it:boottstat} the adversarial tests are reported as absent. Refused with
{it:bootstrap(bayes)}, with {it:poisson}, and with {it:bootdisplay(wald)}, {it:bootdisplay(bc)} and {it:bootdisplay(bca)}. It is materially slower than a plain
bootstrap, since every replication estimates its own variance. With {it:cluster()} at or above the level of {it:i} the analytic pairwise
contrast standard errors can be degenerate, so the studentized contrasts and the subgroup equality test may be absent or unstable; use the plain bootstrap for them
there. {p_end}

{phang}{opt bootsa:ve(name)}: Also write the replication draws to {it:name}.dta in the current working directory; refused without {it:bootstrap()}. {it:name} is a
plain file name, and an existing file is overwritten only with {it:replace}. The file has one row per replication, {it:_rep} = 1 to {it:reps()}, with the
coefficients under their {it:e(b)} names, and a row {it:_rep} = 0 with the observed estimates; a failed replication keeps its row with every value missing. Under
{it:boottstat} the replication analytic standard errors are added as {it:se_*}, and with {it:cluster()} one {it:w_} column per cluster records how many times it was
drawn, or its Dirichlet weight under {it:bootstrap(bayes)}. Its characteristics record the engine, seed and random-number state. {p_end}

{phang}{opt ties}: Adjust the bias-corrected intervals' median-bias constant for ties between the replications and the observed statistic, as the {it:ties} option
of Stata's {help bootstrap} does; refused without {it:bootstrap()}. It has no effect on continuous statistics. {p_end}

{phang}{opt dots(#)}: Show a progress dot every {it:#} replications; {it:dots(0)} suppresses the display, and the default is 1. After the fifth replication a projected
runtime for the loop is printed on the dots line. Requires {it:bootstrap()} or {it:permute()} and is refused with neither. {p_end}

{phang}{opt rseed(#)}: Seed for the bootstrap or permutation draws, an integer between 0 and 2,147,483,647 passed to {it:set seed}.
Requires {it:bootstrap()} or {it:permute()} and is refused with neither. The seed is stored in {it:e(boot_seed)} or {it:e(perm_seed)}. The random-number state
before the first draw is always stored in {it:e(boot_rngstate)} or {it:e(perm_rngstate)}, so a run without {it:rseed()} can be replayed with {it:set rngstate}
before repeating the call. {p_end}

{marker permute}{...}
{dlgtab:Randomization (Permutation) Inference}

{phang}{opt perm:ute(method)}: Randomization inference over the assignment of treatment onset, intended for settings with few treated clusters. {it:wooldid} reassigns onset dates across
clusters, re-estimates every requested quantity under each reassignment, and ranks the observed estimates in the resulting distributions. {it:shuffle} permutes the observed onset labels
across clusters, holding the count at each onset fixed; {it:probabilities} draws each cluster's onset from the probabilities in {it:permprobs()}. Requires a single {it:cluster()} variable,
the randomization unit, that nests {it:i} unless {it:permonsetcohorts} is specified; every cluster must sit inside one treatment-onset group, with the never-treated their own group and coded
as plain missing. That screen runs on every row the {it:if} expression keeps, including rows the estimator later drops for a missing outcome or weight, so a specification that runs without
{it:permute()} can be refused once {it:permute()} is added; restrict the sample with {it:if}. Also requires {it:permreps()} or, for shuffle only, {it:permexact}. Not compatible with
{it:impvalues}, {it:unconditionalse}, {it:contreatwithin}, {it:contreatcontrols(untreatedZero)}, {it:oldsyntax}, {it:suppresssgprimaryeffects}, {it:pretrendpower()}, {it:simultaneousbands}
or {it:bootstrap()}. A design that admits only one possible assignment is refused. With {it:subgroups()} the subgroup screens apply over the whole sample, since any cluster can be treated in
a replication. Every p-value and interval the run reports is permutation-based, including the subgroup contrasts, the joint pre-trend, subgroup equality and adversarial tests, and the
continuous-treatment plot cells. P-values are two-sided equal-tailed, using the add-one form for sampled draws. Intervals are shift intervals at {it:level()}: the estimate minus the
equal-tailed quantiles of the permutation distribution. By default the run and the replications are in fast mode, without analytic standard errors, so the {it:se_analytic} column is blank
and {it:e(V)}, {it:e(b_es_*)} and {it:e(V_es_*)} are not posted; {it:permanalyticse} restores them. Quantities with no permutation distribution are left blank rather than filled
analytically: the {it:gteffects} cells, whose keys are the onsets the replications reassign, and the continuous-treatment plot cells of the {it:baseTreated}, {it:byTreated} and {it:byTtre}
lattice types. The {it:histogramcohorteffects} output stays analytic, with a caution. Replications that fail to estimate are excluded and counted in {it:e(perm_failed)}; the run is refused
when fewer than max(20, half) complete. A caution prints when fewer than 6 clusters are ever treated, since the set of distinct assignments is then small whatever {it:permreps()} is. {p_end}

{phang}{opt permreps(#)}: Number of sampled reassignments, drawn with replacement; at least 20. Required unless {it:permexact} is specified. A sampled-draw p-value cannot fall below
2/({it:permreps}+1) under the default convention, or 1/({it:permreps}+1) under {it:permabsp} or {it:permside()}, and a caution prints when that floor is above .05. A caution also prints when
{it:permreps()} exceeds 10 percent of the possible assignments. {p_end}

{phang}{opt permex:act}: Enumerate every distinct assignment instead of sampling. Only with {it:permute(shuffle)}, and not together with {it:permreps()}. A caution prints above 1,000
assignments and enumeration is refused above 65,534, which is the largest number of rows a Stata matrix can hold and so the point at which the replication
accumulators stop working. Above 10,000 enumerated assignments the draw matrices are not posted and a caution says so; {it:permsave()} writes them regardless. P-values are then plain shares over the {it:M} assignments, with a floor of 2/{it:M}, or 1/{it:M} under {it:permabsp}. {p_end}

{phang}{opt permstr:ata(varname)}: Confine the shuffle to within strata of a numeric variable, so each onset group keeps its composition across strata. Every cluster must sit inside a single
stratum, with no missing values. Only with {it:permute(shuffle)}; under {it:permute(probabilities)} strata are expressed through {it:permprobs()} instead. {p_end}

{phang}{opt permpr:obs(spec)}: The assignment-probability design for {it:permute(probabilities)}, written as {it:(onsetvar level level ...)} followed by one {it:(clustervalue p p ...)} group
per cluster. The first group names the timing variable and lists every in-sample onset level, with never-treated as {it:.}; each later group gives one cluster's probabilities over those
levels, nonnegative and summing to 1. Every cluster value must appear exactly once. Requires {it:permreps()}; {it:permexact} is refused. A caution prints if the observed assignment has zero
probability under the design. {p_end}

{phang}{opt permts:tat}: Randomize on t-ratios instead of coefficients: each replication's estimate over its own analytic standard error. The replications then compute standard errors, so
the run is slower, and {it:fastnose} is refused. Intervals become the estimate minus the equal-tailed quantiles of the permuted t-ratios times the analytic standard error. The joint,
equality and adversarial tests and the pairwise contrasts are unchanged. It protects against a treated cluster with unusually dispersed errors only where the cluster-robust standard error
sees that cluster's own variance, which at {it:cluster()} equal to {it:i} it does not; a caution says so. {p_end}

{phang}{opt permt:ci}: Also report, for the pooled estimands, the symmetric inverted-|t| interval: the estimate plus or minus the {it:level()} quantile of the permuted absolute t-ratios
times the observed standard error. Implies {it:permtstat}. {p_end}

{phang}{opt perman:alyticse}: Without {it:permtstat}, compute the analytic standard errors anyway and show them in the {it:se_analytic} column for reference; the p-values and intervals
remain permutation-based, and the posted {it:e(V)} is analytic. {p_end}

{phang}{opt permab:sp}: Report absolute-value permutation p-values instead of equal-tailed ones: the share of replications whose absolute statistic is at or above the observed absolute
statistic. Intervals become the estimate plus or minus the {it:level()} quantile of the permuted absolute statistics, so p-value and interval keep agreeing. The joint, equality and
adversarial tests are unchanged. At most one of {it:permabsp}, {it:permside()} and {it:permrawp} may be specified. {p_end}

{phang}{opt permsi:de(left|right)}: One-sided permutation p-values: {it:permside(left)} is the share of replications at or below the observed statistic and {it:permside(right)} the share at
or above it. The intervals stay the two-sided shift intervals. {p_end}

{phang}{opt permraw:p}: The plain absolute-value share with no add-one adjustment, which is the number {it:ritest} and Stata's {it:permute} prefix report. It can print exactly 0, meaning no
draw reached the observed statistic; intervals are those of {it:permabsp}. {p_end}

{phang}{opt permstep:down}: Report, for each family of two or more horizons, the Romano and Wolf max-|studentized| step-down computed from the same replications, which controls the
family-wise error rate across the horizons. Post-treatment and pre-treatment horizons are adjusted as separate families. Results are stored in {it:e(perm_stepdown*)} and summarised in the
footer; a family whose step-down cannot be formed is silently not posted. {p_end}

{phang}{opt permsa:ve(name)}: Write the permutation draws to {it:name}.dta in the working directory; an existing file needs {it:replace}. One row per replication plus a row 0 with the
observed estimates; one column per coefficient of {it:e(b)}, the replication standard errors under {it:permtstat}, and one {it:onset_clustervalue} column per cluster holding its onset in
each replication. {p_end}

{phang}{opt permonset:cohorts}: Assert that {it:i} is the treatment-onset grouping, one level of {it:i} per distinct onset with the never-treated a single level, and on that assertion accept
a {it:cluster()} finer than {it:i}; each replication then rebuilds {it:i} from the onsets it drew. The assertion is checked and the run stops if it fails; {it:safety(off)} does not lift the
check. Compatible with {it:permstrata()} and {it:permute(probabilities)}, not with {it:subgroups()}. {p_end}

{pstd} {ul: Discussion of Permutation Intervals and Pre-Treatment Tests:}
The permutation distribution of a coefficient is generated by relabelling outcomes that carry the true effect, so the shift interval widens as the true effect grows and is calibrated only
near the null; {it:permtstat} removes that growth for the post-treatment rows only, and the implied standard error in {it:e(perm_se_shift)} is not a sampling standard error. Every
pre-treatment p-value, including the joint pre-trend tests, tests the sharp null of no effect at any horizon. In staggered designs the reassigned onsets move post-treatment outcomes into the
draws' pre-treatment horizons, so the pre-period p-values rise toward 1 as the true effect grows; a large pre-period p-value is then not by itself evidence of no pre-trend. With a single
onset date this does not arise. {p_end}

{marker supporting}{...}
{dlgtab:Supporting Features}

{phang}{opt hist:ogramcohorteffects}: Produce a histogram of cohort i-specific treatment effects; return those effects in {it:e(histogramestimates)} for the user to manipulate. Under
{it:att_it} or {it:att_itime}, the column averages each cohort's cell effects. With {it:contreat()}, average marginal effects are histogrammed and appended to the matrix. Under
{it:bootstrap()} and {it:permute()}, the histogram estimates come from the initial analytic fit, not the replications. {p_end}

{phang}{opt save:plots(string)}: Save all graphs generated to .gph files beginning with the name specified in the option. Also saves .dta files beginning with the same
name containing the data used to generate any requested event study plots, one per subgroup, and, with {it:makecfxplots}, one {it:name}_cfx.dta of the plotted points. Graphs and datasets
generated have a suffix appended to the name specifying the exact graph/dataset
generated. Everything is saved to Stata's working directory; the argument is a name, not a path. {p_end}

{phang}{opt replace}: Overwrite any existing files with the names specified in {it:saveplots()}, {it:permsave()}, or {it:bootsave()}. Without it, an existing {it:permsave()} or
{it:bootsave()} file is an error, raised before the replications run. {p_end}

{phang}{opt sum:marystats}: Compute and present summary statistics of the outcome variable y and, if present, of the continuous treatment variable contreat.
Summary statistics are computed overall, within the treatment group by time period (post-treatment vs pre-treatment), and within the control group; {it:contreatwithin} omits the control
group rows, and {it:contreatcontrols(untreatedZero)} omits the pre-treatment and control group rows for contreat. Without a
fixed reference period, and with pre-treatment relative time period effects requested, the size of the reference group relative to the rest of the
pre-treatment period is also shown. When event study horizons are requested, the event study composition table {it:e(escomp)} is also displayed, overall and per subgroup. {p_end}

{phang}{opt pois:son}: Implement the specified {it:wooldid} command using poisson regression, achieved by switching from the regression command {it:reghdfe} to {it:ppmlhdfe}, or from
{it:regress} to Stata's {it:poisson} under {it:unconditionalse} or {it:customvce()}. Treatment effects
are presented in terms of differences in linear predictions (the poisson {it:predict(xb)} option). Poisson regression is useful when one would like to use
log(y) as an outcome but y contains natural 0s, or when treatment and control groups are only on parallel trends in proportions.
Poisson accepts pweights but not aweights, refuses a negative outcome, and may not be combined with {it:boottstat}. If {it:ppmlhdfe} fails with r(3200), {it:wooldid} retries once with the
{it:ccc_absorb()} and continuous treatment controls moved out of the absorbed effects. A treated cohort-period cell whose outcome is zero on every row is perfect separation rather than
multicollinearity; {it:wooldid} names the cell and exits. Coarsen {it:i}, drop the cohort, or specify {it:safety(off)}, which excludes the separated cells from the estimated effects. The
variable name {it:_ppmlhdfe_d} is reserved on poisson runs. {p_end}

{marker eventstudy}{...}
{dlgtab:Event Study Estimation}

{phang}{opt esf:ixedbaseperiod}: Use a fixed reference period for estimation, rather than a pooled reference period. Under this option the {it:pretreatment} row of the results
table averages every estimable pre-treatment period, not only the window named by {it:esprelength()}; the joint pre-trend test uses the requested window. {p_end}

{phang}{opt esr:elativeto(-integer)}: The relative time period used as the fixed reference period, if a fixed reference period is being used. The default is -1, the time period immediately
prior to treatment; the onset period or later is rejected with an error. Without {it:esfixedbaseperiod} this option has no effect and a caution is printed. {p_end}

{phang}{opt espre:length(integer)}: The number of pre-treatment relative time period-specific effects to be estimated. If any pre-treatment effects are requested and
{it:esfixedbaseperiod} is specified, all relative time period effects between the {it:esrelativeto()} period and treatment onset will be presented. Both this option and
{it:espostlength()} count consecutive integer relative periods rather than realized horizons, so on a panel whose relative periods have gaps (a biennial panel, for example) the
reported window can be shorter than the one requested; widen the request until it reaches the horizons you want. {p_end}

{phang}{opt espost:length(integer)}: The number of post-treatment relative time period-specific effects to be estimated. With {it:subgroups()} and {it:esprelength()} greater than 0,
a value of 0 is raised to 1 under a caution; {it:espostlength(-999)} keeps the post-treatment window empty; other negative values are an error. {p_end}

{phang}{opt joint:tests}: Calculate joint statistical significance tests of all requested pre-treatment relative time period-specific effects. At least two pre-treatment horizons
must fall in the estimation window (three under {it:esfixedbaseperiod} with {it:timetrends} and the default {it:esrelativeto(-1)}); with fewer, no test is produced. These tests are
unusually sensitive to how standard errors are clustered. With analytic standard errors clustered at the level of the cohort variable {it:i} or above, the saturated cohort-by-period
dummies leave the test's covariance degenerate and its p-value severely oversized, so {it:wooldid} computes and prints the test with a caution that its p-value is not calibrated. The remedies
are {it:bootstrap()}, {it:permute()}, or a {it:cluster()} strictly finer than the cohort-by-period cells that leaves many clusters inside each cell. Under a continuous treatment the
tests also require many observations per cohort-by-period cell. Refused with {it:impvalues} and with {it:fastnose}, except on the fast path {it:permute()} takes by default; under
{it:bootstrap()} and {it:permute()} the reported joint test is a bootstrap or permutation test. {p_end}

{phang}{opt adv:ersarial}: Tests the joint hypothesis that all estimated relative time period-specific effects fall on a single line through the reference period, where the slope of
the line is chosen to fit the estimates as well as possible, a null that is adversarial toward finding pre-trend violations. The p-value comes from a minimum-distance test on the
event-study coefficients and their covariance matrix. It requires
{it:esfixedbaseperiod}: a pooled reference period has no single omitted period, so under the default the fitted line is reported without a p-value, under a caution. Under
{it:timetrends} two periods are normalized to zero and the null reduces to every coefficient being zero. For display, {it:wooldid} also reports the line fit by the unweighted OLS
regression of {it:coef_rt} on {it:time_rt}, where {it:coef_rt} is the estimated coefficient in relative time period {it:rt} and includes a value of 0 for the reference period.
The test is cancelled, with a caution, when the requested window holds fewer than one pre-treatment and one post-treatment horizon (fewer than three coefficients in all under
{it:timetrends}), when the estimation sample fails to realize those horizons, or when the event-study covariance is degenerate; under {it:subgroups()} it is skipped for any
subgroup that lost a relative period as inestimable. Like {it:jointtests}, the test is sensitive to clustering: it is severely oversized unless each treated cohort-by-period cell
is spanned by many clusters, and clustering on {it:i} itself places every cell inside a single cluster. Where that condition is not met, use {it:permute()} or {it:bootstrap()}
with {it:boottstat}, under which the adversarial p-value is a permutation or bootstrap quantity (without {it:boottstat} it is absent under {it:bootstrap()}). Refused with
{it:impvalues} and with {it:fastnose}, except under {it:permute()}. {p_end}

{phang}{opt altmethod_adversarial}: Alternate method of implementing the test above, which projects a trend forward using the pre-treatment coefficients only (together with the
reference period), and then tests whether the post-treatment coefficients collectively fall in line with the projected values, accounting for the uncertainty in the projected line.
It modifies {it:adversarial} and is an error without it. It requires {it:esfixedbaseperiod}, and at least one pre-treatment coefficient beyond the reference period. {p_end}

{phang}{opt pretrendpow:er(numlist)}: Reports Roth (2022) style pre-trend power diagnostics for each requested power level, given in percent as integers strictly between 5 and 100
(for example {it:pretrendpower(50 80)}): a GLS estimate of the linear pre-trend slope with its standard error and confidence interval, the smallest linear pre-trend each pre-test
would detect with the requested power ({it:mdt_wald} for the joint pre-trend test and {it:mdt_nis} for the test that any pre-period coefficient is individually significant, both
at the 5% level), and the bias such a trend would imply for the pooled estimand if it continued into the post-treatment period. Requires {it:esfixedbaseperiod} and at least one
pre-treatment event-study coefficient. Refused with {it:impvalues}, {it:fastnose} and {it:permute()}, none of which produces the event-study covariance the diagnostics use; under
{it:bootstrap()} the diagnostics use the bootstrap covariance. With {it:subgroups()} only the overall event studies are analysed, and when only {it:att_i} or {it:att_itime} is
requested the reported block is the {it:att_it} event study. The diagnostics are skipped with a caution when the pre-period covariance block is unavailable (as under
{it:timetrends}) or ill-conditioned. The {it:mdt_nis} column is computed for up to 4 pre-periods; for larger designs, pass {it:e(b_es_type)}, {it:e(V_es_type)} and
{it:e(numpre_type)} to Roth's {help pretrends} package. Results are returned in {it:e(roth_power_type)} and companion scalars (see {help wooldid##output:Working with Output}). {p_end}

{phang}{opt makep:lots}: Produce event study plots. No plot is drawn unless at least one event-study horizon is requested. {it:wooldid} names its graphs (for example
{it:att_fullmodel}) and overwrites its own prior graphs of the same name on each run; it does not clear other graphs from memory. {p_end}

{phang}{opt oos:treatedcontrols(integer)}: By default, {it:wooldid} uses cohorts treated after the end of the sample the same way as never-treated control cohorts (i.e., such cohorts are
"out-of-sample treated controls"). When this option is specified with a value of T, cohorts treated within T periods of the end of sample are instead handled the same as other treated cohorts.
Such cohorts are potentially allowed to contribute to pre-treatment relative time period-specific effects when this option is specified,
but not otherwise. Note that if {it:esfixedbaseperiod} is specified, the out-of-sample-treated cohorts will be dropped if the specified reference period is not usable for these cohorts.
When {it:anticipation()} is also specified, the window applies to the shifted onset. Negative values are not allowed. {p_end}

{phang}{opt antic:ipation(integer)}: Treats the effective treatment onset as {it:ttre} minus #, so that periods in which units may respond in anticipation of a known upcoming
treatment count as treated periods rather than as clean pre-treatment periods. The shift is applied before relative time is constructed and is equivalent to re-running
{it:wooldid} with {it:ttre} replaced by {it:ttre} minus # (the convention of {it:did_imputation}'s {it:shift()} option): relative period 0 is # periods before the recorded
onset, so you may wish to raise {it:espostlength()} by # to keep the same number of genuine post-treatment horizons. The value is expressed in the units of {it:t}, so on a
biennial panel an odd value lands between observed periods and {it:wooldid} prints a caution. Units whose effective onset falls at or before their first usable period are dropped,
with a caution; under {it:esfixedbaseperiod} the same applies to units whose shifted reference period is unusable. Too small a value leaves a visible residual pre-trend, while too
large a value attenuates the estimated effects without any diagnostic, so prefer the smallest defensible value. With {it:subgroups(ttre)} the subgroup labels remain the recorded
onset values. The shift is recorded in {it:e(anticipation)}. Negative values are not allowed. {p_end}

{pstd}{ul: Two Approaches to Estimation -- Pooled Reference Period vs Fixed Reference Period:} Difference-in-differences treatment effect estimates involve comparing average outcomes during
treated cohorts' treatment periods to treated cohorts' average outcomes during a pre-treatment reference period. By default, {it:wooldid} does comparisons using a pooled reference period,
meaning that treated periods for a given cohort are compared to a reference period containing all pre-treatment observations for the cohort. If the user specifies the option
{it:esfixedbaseperiod}, {it:wooldid} will switch to employing a fixed reference period, meaning that all treatment effect estimates will be formed by comparing treated cohorts' average outcomes
in post-treatment periods to their outcomes in a particular pre-treatment reference period (specified in relative time period terms). By default, this pre-treatment reference period is the
time period immediately prior to treatment (relative period -1), though the user can specify alternative behavior. For more on how this is implemented, see {help wooldid##description:Description}.{p_end}

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
The SG variable must consist of non-negative integer values (or missing values) on the observations of treated units, and its levels must specify collections of cohorts i. For estimates to be valid, a given cohort i cannot
be split across multiple subgroups. Observations given missing values are not placed in any subgroup and are ignored when performing subgroup computations. A unit must be missing on all of
its rows to be left out: a treated unit labelled on some rows and missing on others is refused, and with {it:safety(off)} enters its subgroup weighted by the labelled share of its sample
weight, under a caution. Because subgroup values are
embedded in the names of stored results and plots, values 8 or more digits long can exceed Stata's limits on names; recode ID-style codes to small integers first (for example
{it:egen newsg = group(oldsg)}). With {it:esprelength()} greater than 0 and {it:espostlength()} 0, {it:wooldid} raises {it:espostlength()} to 1 under a caution, since subgroup
estimation requires a post-treatment horizon; {it:espostlength(-999)} keeps the post-treatment window empty. Not compatible with {it:impvalues}. {p_end}

{phang}{opt suppress:sgprimaryeffects}: Skips estimation of subgroup-specific treatment effects and computes only the differences in effects between SGs. This option should
speed estimation if these treatment effects are not otherwise of interest. It is refused by {it:bootstrap()} and {it:permute()}. {p_end}

{pstd} {ul: Description of Approach to Subgroup Effect Estimation:} Given a variable dividing cohorts i into different subgroups, {it:wooldid} will estimate treatment effects separately
for each subgroup (i.e., each level of the SG variable). This is achieved without modification to the underlying
regression specification in step 1 (see {help wooldid##description:here}). Rather, it is done by carrying out the margins calculations in step (2) within the subsamples of the data
specified by the subgroup variable.
Additionally, {it:wooldid} computes all possible pairwise comparisons of the overall
effects (and pre-treatment effects) estimates by SG (e.g., SG 1 vs 2, 1 vs 3, 2 vs 3, etc.), presenting for each comparison the difference in effect estimates and a p-value for each difference.
A joint test of the hypothesis that all SGs share the same effect estimate (normalized to be the effect estimate of the numerically lowest numbered SG) is also computed. Note: standard errors
for subgroup-specific estimates can fail when using clustered standard errors in conjunction with subgroups that contain few (or only one) cluster. Separately, with clustering at
the level of the cohort variable {it:i} or above, the analytic standard error of a pairwise contrast can be degenerate (its standard error is rounding residue rather than an estimated
variance); {it:wooldid} then reports that contrast, and the across-subgroup equality test covering it,
without inference, under a caution naming the count. With {it:safety(off)} they are reported as computed. The remedies are {it:permute()}, {it:bootstrap()} without
{it:boottstat}, or a {it:cluster()} strictly finer than the cohort-by-period cells. {p_end}

{pstd} Each treated unit of {it:i} should lie wholly within one subgroup: a subgroup variable that varies within a cohort-by-period cell, or over time within a unit, is refused
with an informative error. With {it:safety(off)} such a run proceeds, entering a split unit into each subgroup weighted by its share of sample weight there, which changes what
{it:att_i} and {it:att_itime} mean. The usual remedy is to define {it:i} at the finer level over which the subgroup varies, for example {it:egen newi = group(i sg), missing};
the {it:missing} option matters, since without it rows with a missing subgroup value would be dropped from the estimation sample. Note that splitting {it:i} this way makes
{it:att_i} and {it:att_itime} average over the split pseudo-units rather than over the original units. {p_end}

{pstd} A special case to consider is the difference-in-differences-in-differences (triple diff) case. In the triple diff case, all cohorts are categorized into two groups, with both groups
existing among treated cohorts and control cohorts. The target estimand may be thought of as the difference in the diff-in-diff estimates obtained separately for each of the two groups. To
implement a triple diff design in {it:wooldid}, one can specify the variable containing the two triple diff groups as a SG variable and then manually specify inclusion of SG by period
fixed effects (i.e., specifying {it:subgroups(sgvarname) fe(i.sgvarname##i.t)}). The difference in SG level estimates returned will be the triple diff estimate. {p_end}

{marker continuous}{...}
{dlgtab:Continuous Treatments}

{phang}{opt contreat(varname)}: The continuous treatment variable. Must be defined for all observations (set to 0 if not naturally defined within untreated cohort-periods);
observations on which it is missing are dropped from the estimation sample. If the continuous treatment carries no identifying variation for any treated cohort-period cell's
dose slope (for example, a dose constant within each cohort), the within-cell average marginal effects are not identified: under default safety {it:wooldid} names this situation, reports the ame
results as missing, and continues with the att results. Not available with {it:gteffects}. {p_end}

{phang}{opt contreatelas:ticitytype(dydx|dyex)}: Specifies which type of marginal effect of the continuous treatment to compute. {it:dydx} is the default and reports the
average marginal effect itself (change in y per unit change in the continuous treatment). {it:dyex} reports the change in y per proportional change in the continuous treatment.
Under {it:poisson} results are on the linear index, where {it:dydx} already gives the semi-elasticity of E[y] and {it:dyex} the elasticity. Only {it:dydx} and {it:dyex} are accepted; the
older {it:eyex} and {it:eydx} forms have been removed. {p_end}

{phang}{opt contreatw:ithin}: Specify this when you wish to use only within-treated cohort variation in the continuous treatment variable. This option relaxes the estimability
checks that would otherwise cause {it:wooldid} to stop estimation in the absence of any untreated cohorts (the data-integrity checks still run unless {it:safety(off)} is given),
and causes {it:wooldid} to skip sample cleaning routines that drop time periods for which a control cohort is not present. Requires {it:contreat()}. {p_end}

{phang}{opt contreatp:oly(integer)}: Includes polynomial terms in {it:contreat} of the specified degree, a positive integer (the default is 1). Be cautious using higher degree polynomials. {p_end}

{phang}{opt lattice(latticevar OR latticetype latticenum)}: Specifies a lattice variable, a variable specifying groups of observations across which the effect of {it:contreat} is
allowed to vary, even within cohort-periods and within the reference period. If just the name of a variable is passed to {it:lattice()}, this variable will be used as the lattice variable;
it must be numeric and contain only non-negative integers, and observations on which it is missing are dropped from the estimation sample.
If instead a lattice type and a lattice number are passed, this option will generate a lattice variable using {it:latticenum} quantiles (an integer greater than 1), with the quantiles being formed within groups
according to the {it:latticetype} (see below). Any specified sample weights will be used when constructing the lattice. If {it:lattice} is not specified, {it:wooldid} will proceed without
using a lattice variable (as though you specified a lattice variable that is always 1). A lattice is the most expensive option in runtime terms (see {help wooldid##runtime:Runtime}).
Requires {it:contreat()}. When a lattice is combined with {it:contreatcontrols(basic)}, the baseline dose level is partialled out within each lattice level. {p_end}

{phang} {it:latticetype:} Specifies groups in which to estimate {it:latticenum} quantiles of {it:contreat}; type names are case-sensitive. Choice of lattice type can be important for both interpretive reasons
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
variable or, if specified, the variable in {it:cfxlattice()}; one of the two must be present. These plots are only sure to be valid for the ATT estimand type, but may be valid for others if lattice groupings
do not crosscut cohorts (or, less stringently, cohort-periods, for ATT_it), so {it:wooldid} produces them for the {it:att} estimand only unless {it:safety(off)} is specified.
Requires {it:contreat()}; not compatible with {it:impvalues}. {p_end}

{phang}{opt cfxlattice(varname)}: A variable containing the bins in which effects are to be calculated for the continuous treatment effects plots. If not specified,
the lattice variable will be used. The variable must be numeric and contain only non-negative integers on treated observations; where it is missing on treated observations,
{it:margins} omits those observations from the plotted averages and {it:wooldid} prints a caution. Without {it:makecfxplots} the option has no effect and draws a caution. {p_end}

{phang}{opt cfxplott:ypes(pre post integers)}: A list of types of cfxplots to be produced. If {it:post} is among the types specified, plots will be produced giving
treatment effects across the continuous treatment variable within the post-treatment period. If {it:pre} is among the types specified, plots will be produced within the pre-treatment period.
If any integers are specified, plots will be produced just within the specified periods relative to treatment. If nothing is specified, the option post is assumed. {it:pre} requires
{it:esfixedbaseperiod} or {it:esprelength()} greater than 0. Without {it:esfixedbaseperiod} an integer must lie within the requested event-study window and entries outside it are
dropped without a message; with {it:esfixedbaseperiod} any integer is accepted except the reference period itself, which is dropped. An entry that is neither {it:post}, {it:pre} nor an integer is an error. {p_end}

{phang}{opt makecfxplotsbysg}: Produce continuous treatment effects plots for each group specified in the subgroup variable. Has no effect, and prints nothing, without both {it:makecfxplots} and {it:subgroups()}.
Each subgroup's plot data is returned in {it:e(cfxitem_type_sgN)}, where N is the subgroup value. {p_end}

{phang}{opt cfxplotnose}: Do not estimate standard errors for the cfxplots; speeds estimation. {p_end}

{phang}{opt contreatc:ontrols(controltype)}: Specifies what controls for {it:contreat} will be added: the specified controls are always equal to f(Contreat_jt) (the function used to capture
  the effect of contreat within treated i-t cells) interacted with a set of fixed effects-like indicator variables (time period, cohort, etc.). The specified controls will partial out
the effect of {it:contreat} in a particular way, affecting the nature of the reference effect of {it:contreat} to which estimated effects of {it:contreat} in treated cohort-periods are
compared. This option must be specified to proceed with a continuous treatment, and care must be chosen when specifying it: different options correspond with very different
substantive assumptions. See more on this topic in the description section below. {p_end}

{phang} {it:conttreatcontroltype}: untreatedZero (which cannot be paired with any other type); or one or more of the following: basic, byTreated, byTtre, byCohort, and byTime  {p_end}
{synopt : {it:basic}} Controls for the effect of contreat (i.e., for f(Contreat_jt)) are added without any interactions.  {p_end}
{synopt : {it:byTreated}} Controls for the effect of contreat are interacted with indicators for membership in the treatment group and in the control group. {p_end}
{synopt : {it:byTtre}} Controls for the effect of contreat are interacted with indicators for membership in each treatment onset group and in the control group.  {p_end}
{synopt : {it:byCohort}} Controls for the effect of contreat are interacted with cohort (i) membership indicators. {p_end}
{synopt : {it:byTime}} Controls for the effect of contreat are interacted with time period (t) indicator variables. {p_end}
{synopt : {it:untreatedZero}} Special option to formulate the regression to work when a continuous treatment variable is always equal to 0 outside of treated cohort-periods. Setting
the continuous treatment variable to always equal 0 outside of treated cohort-periods can be an appropriate way to handle cases where the variable is unknowable or undefined outside
of treated cohort-periods. If this option is specified when the always 0 condition is not met, {it:wooldid} will return poorly identified estimates. {p_end}

{pstd} Each control type requires variation in {it:contreat} across the reference-period observations of every group for which it fits a baseline dose slope: the whole
sample under {it:basic}, each treatment-status group under {it:byTreated}, each onset group under {it:byTtre}, each cohort under {it:byCohort} and each calendar period under
{it:byTime}. A group without it is named and the run stops; under {it:byTime} this happens whenever never-treated units carry a zero or constant dose in the periods after the
last onset. {p_end}

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

{pstd} Take care also not to over-saturate the specification with lattice terms. Each treated cohort-period cell must hold at least two distinct values of Contreat_jt within each
lattice bin it occupies, since that is the only variation identifying that piece's slope; pieces without it are dropped for collinearity, the reported average marginal effect then
averages only the pieces that survived, and the analytic standard errors can be far smaller than the estimates' sampling variability. A useful rule of thumb is that the number of
units contributing to a treated cell should be several times the number of lattice bins. {p_end}

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
be the integral of the estimated average marginal effects across the distribution of the continuous treatment variable. Note: the control type {it:byTime} is refused in this case,
since in the absence of an untreated group, the time period interactions will force {it:wooldid} to partial out genuine treatment effects for some group, threatening identification;
{it:safety(off)} downgrades the refusal to a caution. {p_end}

{pstd} {ul:Continuous Treatments and Cluster-Robust Standard Errors:} With {it:cluster()} at the level of {it:i} or coarser, every treated cohort-by-period cell lies within a single
cluster, and the analytic standard errors of the ame-family quantities are then not identified. {it:wooldid} prints a caution, reports those standard errors, the ame joint tests and
the ame adversarial p-values as missing, withholds {it:e(V_es_ame*)}, and draws the ame continuous treatment effects plots without confidence intervals; the att results are
unaffected. Valid ame inference under such designs comes from {it:bootstrap()}, {it:permute()}, {it:impvalues}, {it:robust} (when errors are independent across observations), or a
{it:cluster()} strictly finer than the cohort-by-period cells. {p_end}

{marker controls}{...}
{dlgtab:Controls and Fixed Effects}

{phang}{opt control:s()}: The contents of this option are passed directly to {it:reghdfe} as controls, without any modification (to {it:ppmlhdfe} under {it:poisson}, and to
{it:regress} or {it:poisson} under {it:customvce()} or {it:unconditionalse}). {p_end}

{phang}{opt fe()}: The contents of this option are passed directly to {it:reghdfe} as fixed effects (into its absorb statement) without modification. Variables given will be interpreted as
categorical variables, though i.X##c.X2 type notation is also permitted for passing heterogeneous slope controls and the like. Under {it:customvce()} or {it:unconditionalse},
estimation runs through {it:regress} rather than {it:reghdfe}, and a bare variable name is then read as a continuous variable; write the {it:i.} prefix explicitly there. {p_end}

{phang}{opt timetrend:s}: Control for cohort specific time trends in the regression (i.e., for i.i##c.t). Importantly, this option is required for heterogeneous time trend controls to work,
even when they are specified manually, as it also updates {it:wooldid}'s various estimability checks, drops cohorts without sufficient pre-treatment data to estimate a time trend,
and, in the fixed reference period case, ensures that the two reference periods that must be normalized to 0 to estimate the time trend are chronologically adjacent. This option does not
automatically include any time trends that are a function of any included continuous treatment variable. With {it:cluster()} on the cohort identifier {it:i}, the cohort trend
terms typically leave the covariance unidentified: {it:wooldid} says so and continues, reporting the estimates without standard errors; the same specification under
{it:robust} is unaffected. {p_end}

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
when specifying that variable as a continuous treatment variable and then separately including controls that interact the continuous treatment variable with time dummy variables.
The treatment-effect slopes exist only where a control varies within cohort-by-period cells; where it does not (for example {it:i} is the unit and the control is a unit-level
variable), the run is identical to {it:controls(i.t#c.x)} and {it:wooldid} prints a caution naming the control. This option also turns the default cell collapse off (see
{it:nofastcollapse}). {p_end}

{pstd} {ul:Caution Regarding Continuous Treatment Variables and Controls:} If seeking to manually specify additional controls as a function of a continuous treatment variable, be sure not
to include the continuous treatment variable itself among your controls, as this can adversely affect estimates produced in {it:wooldid}'s margins calculation stage. It is okay to include
a control with values identical to that of the continuous treatment variable: it just needs to have a different variable name. It is recommended that any custom such controls be specified
in conjunction with the option {it:contreatcontrols(basic)}, or any option other than {it:untreatedZero}. {p_end}

{marker technical}{...}
{dlgtab:Technical Features}

{phang}{opt update}: If specified, the command {it:wooldid,} {it:update} will update your local installation of {it:wooldid} based on the most recent edition of the program available on github.
A running Stata session keeps using the copy it has already loaded, so the update takes effect after you restart Stata or run {cmd:discard}. {p_end}

{phang}{opt customw:eightmultiplier(varname)}: When computing various treatment effects using the margins calculation (see step 2 {help wooldid##description:here}), this variable is used to
weight the observations (or if a weight would already be used, that weight is multiplied by this variable). This allows production of certain bespoke types of treatment effects.
Values must be non-negative (negative values are rejected with an error, bypassable with {it:safety(off)}); missing values drop the observation. {p_end}

{phang}{opt clean:matrices}: If specified, {it:wooldid} will clear out all matrices stored in Stata's main memory and in mata. Occasionally, {it:wooldid} enters a state where it repeatedly
crashes for difficult to determine reasons (relating to {it:reghdfe} crashing and leaving certain content in mata). Doing this can sometimes recover {it:wooldid} from this problem. {p_end}

{phang}{opt regtol(#)}: Specifies the convergence tolerance to be used by {it:reghdfe} and {it:ppmlhdfe}, as the negative base-ten exponent of the tolerance ({it:regtol(9)} means
a tolerance of 1e-9). By default, {it:reghdfe} sets a value of 8; {it:wooldid} has a default of 9.
You may wish to set this higher or lower; higher may be advisable in settings where you ask {it:reghdfe} to absorb interaction terms with many continuous variables -- as tends
to occur when using continuous treatment variables. It has no effect under {it:unconditionalse} or {it:customvce()}, which estimate through {it:regress} or {it:poisson}. {p_end}

{phang}{opt ccno:absorb}: Forces continuous treatment controls to be included in the main regression statement, rather than in the {it:absorb()} option, which is the default
whenever using the commands {it:reghdfe} or {it:ppmlhdfe} to estimate the underlying regression. Under certain conditions, this can improve estimation stability and avoid problems
with variables being inappropriately dropped from the regression. {p_end}

{phang}{opt ccc_absorb()}: An alternate version of coarsecohortcontrols() that admits only continuous variables and passes the controls to the {it:absorb()} statement when using
a high dimensional fixed effects estimator. This is often less efficient than the alternative approach in terms of computation time, but can be important in terms of memory
efficiency in the main regression or under certain other circumstances. Be extra careful to ensure that specifications featuring {it:ccc_absorb()} are identified; not all cases
of identification failure are detected. In particular, with an event-study window, variables that are constant within cohort-by-period cells (a unit-level variable with a
unit-level {it:i}, for example) absorb the treatment interactions themselves and the run stops with an error naming this option. Under {it:contreat()}, naming the treatment
variable itself here turns the default cell collapse off (see {it:nofastcollapse}). {p_end}

{phang}{opt ver:bose}: Prints the regression command {it:wooldid} uses and prints progress reports throughout the estimation process; the regression itself is run quietly.
This option also stores the underlying regression model used by {it:wooldid} to {it:___wooldidfullmodel}; this model can be accessed using {it:estimates restore}, though only for
coefficient-level work ({it:lincom}, {it:test}), since its regressors are temporary variables that no longer exist after {it:wooldid} returns. {p_end}

{phang}{opt gte:ffects}: Returns in {it:e(gteffects)} (and prints, for tables of up to 200 rows) the individual cohort-by-relative-period effect estimates underlying {it:wooldid}'s
aggregates, one row per estimated {it:i}-by-relative-period cell, read directly off the underlying regression. Columns: {it:i}, {it:onset}, {it:caltime}, {it:relyr}, {it:estimate},
{it:se}, {it:p}, the interval bounds (named for the run's {it:level()}), {it:inwindow} (1 if the cell falls inside the reported event-study window), {it:nobs}, and the cell's total
estimation weight under each requested estimand ({it:w_att}, {it:w_att_it}, {it:w_att_i}, {it:w_att_itime}). The weighted mean of {it:estimate} under one of these weight columns,
taken over the post-treatment rows ({it:relyr} at or above 0), the pre-treatment rows, or a single {it:relyr}, reproduces the corresponding aggregate exactly, so the table supports
custom aggregations. With {it:subgroups()} a single full-sample table is returned. Standard errors are the regression's coefficient-level
ones: with {it:cluster()} at the level of {it:i} or coarser a treated cell lies inside a single cluster and its standard error is understated, which a caution reports; a cell whose
variance is not identified is reported with a missing standard error. Not available with {it:contreat()}, where the cohort-by-period coefficients are dose-zero intercepts rather than
treatment effects; available under {it:fastnose}; under {it:bootstrap()} the cell inference is bootstrap-based, and under {it:permute()} the cells carry no inference. {p_end}

{phang}{opt emptycellsoverride}: By default, {it:wooldid} tells Stata to drop interaction terms that correspond with empty cells, which is achieved by the command
{it:set emptycells drop}, which improves estimation efficiency in this setting. This option overrides this behavior and roughly doubles runtime. This should only be of interest in conjunction
with the option {it:verbose}. {p_end}

{phang}{opt safety(on|off)}: When {it:safety(off)} is not explicitly specified, {it:wooldid} will caution the user or halt estimation when {it:wooldid} is likely to yield unreliable estimates.
Most important among these checks are tests for whether or not a variable intended to capture the treatment effect in a given i-t cell was dropped from the underlying regression for
multicolinearity, as may occur if one's fixed effects or continuous treatment control variables partial out one of the main effects. Identification failures of this sort will generally be detected,
though not necessarily when the problem concerns continuous treatments being analyzed with a lattice or with polynomial terms. Default safety also checks the data (that onset is
constant within {it:i}, that weights and the {it:customweightmultiplier()} variable are non-negative, that lattice variables are non-negative integers), cautions when {it:i} takes
more than 100 distinct values (see {help wooldid##runtime:Runtime}), and halts with an informative error before estimation when the specification would exceed the current
{help maxvar} setting. {it:safety(off)} skips these checks, and downgrades a number of refusals described elsewhere in this help to cautions. If a run changes your specification,
{it:wooldid} says so whatever {it:safety()} is set to. Only {it:on} and {it:off} are accepted. {p_end}

{phang}{opt oldsyntax}: In older versions, {it:wooldid} used the suffix {it:contemp} to flag the event study coefficient containing the period contemporaneous with treatment.
This turned out to make {it:wooldid} relatively incompatible with the popular event study plot production program {help event_plot}, so the {it:contemp} suffix has been
replaced with the suffix {it:post0}. This option reverts {it:wooldid} to using the older coefficient syntax, for compatibility with any older code that relies on it. It cannot be
combined with {it:permute()} or {it:bootstrap()}. {p_end}

{marker runtime}{...}
{title:Runtime}

{pstd} Runtime is driven by the number of cohort-by-period interaction cells {it:wooldid} builds, roughly the number of distinct values of {it:i} times the number of relative
periods in the event-study window, and it grows faster than linearly in that count, while sample size, subgroup count and event-study width each cost only linearly. Coarsening
{it:i} is therefore normally the single largest speedup available: define {it:i} as cohorts or treatment onset groups and recover unit-level variation with unit fixed effects
in {it:fe()}. Dropping observations buys little, which also makes subsampling a poor way to shorten a long {it:bootstrap()} or {it:permute()} run; shrinking the window or
coarsening {it:i} works there. {it:lattice()} is the most expensive option, so prefer the smallest {it:latticenum} that answers the question. {p_end}

{pstd} Nearly all of a run is spent inside {help margins}, most of it in the delta-method variance calculation. {it:fastnose} skips that calculation entirely and is worth
several times over when standard errors are not needed, and the default cell collapse (see {it:nofastcollapse}) shrinks the data each {it:margins} call sees; the two compose.
Runs that must be evaluated on the full data ({it:unconditionalse}, {it:customvce()}, {it:interactivecontrols()}, and the continuous-treatment forms listed under
{it:nofastcollapse}) pay the full per-row cost. {it:emptycellsoverride} roughly doubles runtime and is rarely wanted. {p_end}

{marker output}{...}
{title:Working With Output}

{pstd}
{it:wooldid} generates output primarily in the form of matrices and scalars stored in e(), though all main results can also be accessed using standard utilities like {it:esttab}.
The commands {it:est restore} and {it:est save} are compatible with {it:wooldid} for saving its output. The visible {it:wooldid} output in e() can be viewed using the command
{it:ereturn list}; a number of results are posted hidden (the preserved analytic covariances, the per-family bootstrap and permutation tables, and the {it:e(wdid_*)} state that
replay uses) and appear only under {it:ereturn list, all}, though each can be read by name either way. Typing {it:wooldid} by itself after a successful estimation redisplays
the stored results, including after {it:estimates restore}; it does not redraw plots or repeat estimation-time cautions, and it refuses with an error after a failed call until
some {it:wooldid} call completes successfully. {p_end}

{pstd} The coefficient vector {it:e(b)} and its covariance {it:e(V)} carry these column names, which are what {it:esttab}, {it:coefplot}, {it:lincom} and {it:test} address a
coefficient by: {it:main_att} (the headline average treatment effect on the treated), {it:pretreatment_att} (the pooled pre-treatment average), {it:ES_att_pre#} and
{it:ES_att_post#} (the event-study horizons, numbered outward from the reference period), and for the other estimand types the same names with {it:att} replaced by
{it:att_it}, {it:att_i}, {it:att_itime} or {it:ame}. Under {it:subgroups()} each subgroup's rows take the suffix {it:_s#}. {p_end}

{pstd} The primary results are stored in the following matrices: {p_end}
{synoptset 30 tabbed}{...}
{synopt:{it:e(wdidmainresults)}} All primary and event study results. The columns are {it:estimate}, {it:se}, {it:t} ({it:z} under {it:poisson}), {it:p}, the lower and upper interval
bounds (named {it:lb_95ci} and {it:ub_95ci}, or for the run's {it:level()}), and {it:relyr}. The pooled estimands sit under the equations {it:main} and {it:pretreatment}, and the
event-study rows under {it:ES_att}, {it:ES_att_it}, {it:ES_ame} and {it:ES_ame_it}, named {it:post#} and {it:pre#}.{p_end}

{synopt:{it:e(wdidsgresults_{it:sgcase})}} Primary and event study results for the subgroup {it:sgcase}, if subgroups are used and presentation of their
primary results has not been suppressed.{p_end}

{synopt:{it:e(pwsgs)}} Pairwise differences between the subgroup-specific treatment effects, if subgroups are used, with the same columns less {it:relyr}. A row whose analytic
contrast covariance is degenerate carries its point estimate with the inference columns missing (see {help wooldid##interactions:Subgroup Effects}). {p_end}

{synopt:{it:e(b_es_{it:type})}, {it:e(V_es_{it:type})}} The event study coefficient vector and its full variance-covariance matrix, labelled {it:pre#}/{it:post#} with the reference
period omitted, per estimand type ({it:att}; {it:ame} under {it:contreat()}; and {it:att_it} or {it:ame_it} for the cell-weighted estimands, which share one event-study
horizon block), with per-subgroup versions {it:e(b_es_{it:type}_s{it:sgcase})} and {it:e(V_es_{it:type}_s{it:sgcase})}. After {it:ereturn post} these matrices are valid inputs for
{it:lincom} and {it:test}. A direct positional handoff to {it:honestdid} additionally requires {it:esfixedbaseperiod} with the default {it:esrelativeto(-1)}, no {it:timetrends},
at least one pre- and one post-treatment coefficient, consecutive horizons, and a covariance with no degenerate rows; other geometries hand {it:honestdid} the wrong event-time
geometry silently. The covariance is withheld, under a caution, when no horizon has an identified variance ({it:fastnose}, typically {it:timetrends}) and for the ame families
under the cell-containment caution; the coefficient vector still posts. {it:e(es_ndegen_{it:type})} counts horizons whose variance is degenerate and is posted only when such rows
exist. {p_end}

{synopt:{it:e(numpre_{it:type})}} The number of pre-treatment relative periods in the corresponding event-study vector, provided for use as {it:honestdid}'s {it:numpre()} argument.{p_end}

{synopt:{it:e(gteffects)}} The cohort-by-relative-period table produced by the {it:gteffects} option; see that option for the column layout. {it:e(gt_ncontained)} counts the cells
lying wholly inside a single cluster, whose standard errors are understated, and is posted only when such cells exist; on the same rule {it:e(gt_ndegen)} counts cells whose variance is not
identified.{p_end}

{synopt:{it:e(escomp)}} Event-study composition diagnostics, one row per reported horizon, aligned with the event-study rows of {it:e(wdidmainresults)}: {it:ncohorts},
{it:nunits}, {it:nclusters}, {it:nobs}, {it:sumw_att}, {it:wshare_att}, {it:wshare_att_it} and {it:relyr}. Composition shifts across horizons are a common cause of jumpy event-study
profiles. Per-subgroup versions are {it:e(escomp_s{it:sgcase})}; the scalars {it:e(escomp_refN)}, {it:e(escomp_refnunits)}, {it:e(escomp_refncoh)} and {it:e(escomp_refW)} describe
the reference-normalization pool, and {it:e(escomp_trtN)} and {it:e(escomp_trtW)} the whole treated estimation sample. The table is displayed when {it:summarystats} is
specified.{p_end}

{synopt:{it:e(roth_power_{it:type})}} With {it:pretrendpower()}: one row per requested power level, columns {it:power}, {it:mdt_wald}, {it:mdt_nis}, {it:bias_wald} and
{it:bias_nis}. Companion scalars: {it:e(roth_slope_{it:type})}, {it:e(roth_slope_se_{it:type})}, {it:e(roth_slope_ll_{it:type})}, {it:e(roth_slope_ul_{it:type})} (the GLS
pre-trend slope and its interval), {it:e(roth_prep_{it:type})} (the joint pre-trend p-value), {it:e(roth_hbar_{it:type})} (the weighted mean post-treatment relative period),
{it:e(roth_npre_{it:type})}, {it:e(roth_keff_{it:type})}, {it:e(roth_W_{it:type})} and {it:e(roth_nissize_{it:type})}; {it:e(roth_ok)} = 1 marks a completed diagnostics pass, and the hidden
{it:e(roth_cellmax_{it:type})} records the largest treated cell size when the pre-period covariance was rank-deficient.{p_end}

{synopt:{it:e(esband_{it:type})}} Under {it:simultaneousbands}: the band endpoints per event-study family (rows {it:lb}, {it:ub}, {it:relyr}; {it:_s#} per subgroup), with
{it:e(esband_crit_{it:type})}, {it:e(esband_n_{it:type})}, {it:e(esband_engine)} ({it:analytic-simulated}, {it:bootstrap} or {it:bootstrap-studentized}), {it:e(esband_level)},
{it:e(esband_families)} and {it:e(esband_ndraw)}.{p_end}

{pstd} In addition to these core results, the following materials are also returned in e(): {p_end}
{synopt: Joint Tests by Horizon} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented in the scalars e(joint_horiz)
(ATT), e(joint_horiz_i) (the shared ATT_it / ATT_i / ATT_itime horizon block), e(joint_Choriz) (AME; continuous treatments only), and e(joint_Choriz_i) (AME_it; continuous treatments only).
On an analytic run the p-value is that of the {it:margins} joint test; under {it:bootstrap()} and {it:permute()} it is a rank p-value of the replication statistics, which are stored in
{it:e(boot_jointW)} and {it:e(boot_jointWobs)} or {it:e(perm_jointW)} and {it:e(perm_jointWobs)}. {p_end}

{synopt: Joint Tests by Subgroup & Time Period} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented by
subgroup, when used, in e(joint_horiz_s{it:sgcase}), e(joint_horiz_i_s{it:sgcase}), e(joint_Choriz_s{it:sgcase}), and e(joint_Choriz_i_s{it:sgcase}), where {it:sgcase}
is the value of the subgroup. {p_end}

{synopt: Adversarial Line Tests:} P-values from the tests that all relative-time period treatment effect estimates and the reference period fall on the same line
are found in the scalars {it:e(adversarial_pval_att)}, {it:e(adversarial_pval_ame)}, {it:e(adversarial_pval_att_s[SG])}, etc.; they are posted only under {it:esfixedbaseperiod}
and only when the event-study covariance is not degenerate. The fitted line is returned in matrices named {it:e(adversarial_line_att)}, {it:e(adversarial_line_att_it)},
{it:e(adversarial_line_ame)}, {it:e(adversarial_line_ame_it)} and their {it:_s{it:sgcase}} twins, with columns {it:slope} and {it:intercept_rt0}. {p_end}

{synopt: SG Effect Homogeneity Tests} For each estimate type (ATT, ATT_it, AME, etc.), the joint test that all subgroups have the same such post-treatment effect (namely, the
same effect as the lowest numbered SG) may be found in the scalar {it:e(joint_att_allsgssame)} and analogously named scalars for other treatment effect types. For pre-treatment effects,
the test can be found in objects with names like {it:e(joint_pre_att_allsgssame)}. A test whose pairwise covariance is degenerate is stored as missing.{p_end}

{synopt:{it:e(N)}, {it:e(N_dropped)}, {it:e(Nclus)}} Estimation sample observation count, number of observations dropped from initially specified sample, and cluster count.
{it:e(N_clust)}, {it:e(clustvar)}, {it:e(vce)} and {it:e(vcetype)} are the conventional variance-type entries, posted so that {it:estout} and {it:esttab} find what they expect;
{it:e(vce_stage1)} records the first-stage type under {it:unconditionalse}. {it:e(level)} holds the confidence level in force, and {it:e(df_r)} the residual degrees of freedom; under {it:impvalues} {it:e(df_r)} is the smallest of the
per-estimate values, which {it:e(imdf_min)} and {it:e(imdf_max)} bracket. {p_end}

{synopt:R^2 measures} For standard estimation, these include e(r2), e(r2adj), e(r2_within), e(r2_withinadj). For poisson estimation, these include e(r2pseud), a pseudo-R^2 measure. {p_end}

{synopt:{it:e(histogramestimates)}} The cohort i-specific treatment effects used to generate histograms, when requested.{p_end}

{synopt:{it:e(summarystats)}} Summary statistics on the outcome variable and (if present) continuous treatment variable, when requested. The columns are {it:mean}, {it:median},
{it:sd}, {it:iqr}, {it:p25}, {it:p75}, {it:min} and {it:max}; the rows are {it:y_insample}, {it:y_pretreat}, {it:y_posttreat} and {it:y_control}, and under {it:contreat()}
{it:c_insample}, {it:c_pretreat}, {it:c_posttreat} and {it:c_control}. Thirteen are also returned as scalars when {it:summarystats} is given: {it:e(ypre_mean)}, {it:e(ypre_med)},
{it:e(ypre_sd)} and {it:e(ypre_iqr)}, and under {it:contreat()} {it:e(contreat_mean)}, {it:e(contreat_med)} and {it:e(contreat_sd)} with {it:contreatpre_} and {it:contreatpost_}
versions of the same three. {p_end}

{synopt: CFX Plot Info} Information used to generate the continuous treatment effects plots. The information used to produce each plot will be stored in a matrix with the following naming
pattern: e(cfx[item]_[type]), where [item] is either one of att, att_i, att_it, or att_itime (for overall average treatment effects) or one of ame, ame_i, ame_it, or ame_itime (for marginal
effects of the continuous treatment). Meanwhile, [type] is either cfxpost (estimated in the post-treatment period), cfxpre (estimated in the pre-treatment period), or relyrZ (estimated in
relative time period Z, where the letter n represents a negative sign). Under {it:makecfxplotsbysg} each matrix carries the further suffix {it:_sg[sgcase]}. {it:e(wdid_cfxmats)}
lists every cfx matrix a run posted. {p_end}

{synopt:{it:e(sample)}} A function characterizing the estimation sample used; useful for identifying which observations are dropped by {it:wooldid} during its sample formation and sanitation
procedures.{p_end}

{synopt:{it:e(cmd)}, {it:e(depvar)}, {it:e(cmdline)}} The command name, the outcome variable, and the input you fed to {it:wooldid} on the command line.{p_end}

{synopt:{it:e(estimate_type)}} A macro recording whether your treatment effect estimates use a pooled reference period or a fixed reference period. {p_end}

{synopt:{it:e(anticipation)}} The anticipation shift applied to the treatment onset, when {it:anticipation()} is active.{p_end}

{synopt:{it:e(wdid_vcvstate)}} What {it:e(V)} is on this run, in one word: {it:analytic-diagonal}, {it:analytic-blocks} (the default on an analytic run with an event study, whose
within-family blocks {it:e(wdid_vblocks)} names and {it:e(wdid_vblockn)} counts), {it:bootstrap}, or {it:absent}. The cell-collapse results {it:e(wdid_fcol)}, {it:e(wdid_fcskip)}, {it:e(wdid_fcguard)} and
{it:e(wdid_fcncell)} record whether the default cell collapse applied (see {it:nofastcollapse}). These are hidden results; read them by name. {p_end}

{synopt:{it:Other Items:}} {it:e(ll)}, {it:e(ll_0)}, {it:e(rmse)}, {it:e(rss)}, {it:e(F)}, and {it:e(chi2)} (poisson only) are passed through directly
from the underlying regression estimator, wherever available. {p_end}

{synopt:{it:Bootstrap inference}} When {it:bootstrap()} is used, the reported tables carry the bootstrap inference, and e() additionally holds: run-level scalars and macros
({it:e(boot_reps)}, {it:e(boot_completed)}, {it:e(boot_failed)}, {it:e(boot_level)}, {it:e(boot_N_V)}, {it:e(boot_method)}, {it:e(boot_display)}, {it:e(boot_stat)},
{it:e(boot_df)}, {it:e(boot_seed)}, {it:e(boot_rngstate)}, {it:e(boot_alpha)}, {it:e(boot_ties)}, and {it:e(boot_nstraddle)} with {it:e(boot_reid)}, which record whether the replications re-identified duplicated
units); the replication covariance
{it:e(V_boot)} (also posted as {it:e(V)}, with the analytic one preserved as {it:e(V_analytic)}), the replication means {it:e(b_boot)} and bias {it:e(bias_boot)}; the draw
matrices {it:e(bootdist)} and, under {it:boottstat}, {it:e(bootdist_se)} (up to 10,000 replications); the interval and p-value families aligned to {it:e(b)},
{it:e(boot_ci_percentile)}, {it:e(boot_ci_normal)}, {it:e(boot_ci_bc)}, {it:e(boot_ci_bca)}, {it:e(boot_p_percentile)}, {it:e(boot_p_bc)}, {it:e(boot_p_bca)}, with
{it:e(boot_z0)} and {it:e(boot_accel)}; the per-family result tables {it:e(boot_pctl_mainresults)}, {it:e(boot_bc_mainresults)}, {it:e(boot_bca_mainresults)},
{it:e(boot_wald_mainresults)} and {it:e(boot_tpctl_mainresults)} with their {it:_pwsgs} and {it:_sgresults_{it:sgcase}} twins, and {it:e(boot_{it:family}_gteffects)} for the
cell table; the test statistics {it:e(boot_jointW)}, {it:e(boot_jointWobs)}, {it:e(boot_advW)}, {it:e(boot_advWobs)} and the studentized records {it:e(boot_tstud)} and
{it:e(boot_tstud_pwsgs)}; and the macros {it:e(boot_jt_family)}, {it:e(boot_adv_absent)}, {it:e(boot_bc_absent)} and {it:e(boot_gteffects_family)} recording which families were
formed and which were withheld. Further records: {it:e(boot_navail)} (replications available per estimate), {it:e(boot_novar)} and {it:e(boot_sesmall)} (replications
with no, or a very small, standard error under {it:boottstat}), {it:e(boot_tstat)}, {it:e(boot_jk_completed)}, {it:e(boot_v_posted)}, {it:e(boot_V_shortfall)} and
{it:e(boot_ves_failed)} (whether the covariance and the event-study covariances were rebuilt from the replications), {it:e(boot_jt_absent)}, {it:e(boot_sgjt_absent)} and
{it:e(boot_sgjt_reason)} (tests withheld and why), {it:e(boot_sgW)} and {it:e(boot_sgWobs)} (the subgroup-equality statistics), {it:e(bootdist_pwse)}, the hidden
{it:e(boot_ptails)}, {it:e(boot_comparator)}, {it:e(boot_ase_obs)}, {it:e(V_es_analytic_{it:type})} and {it:e(df_r_analytic)}, and the per-family cfx twins
{it:e(boot_{it:family}_{it:cfxname})}. {it:e(df_r)} is aligned to the reported reference distribution, and {it:e(boot_percentile_mainresults)} is an alias of
{it:e(boot_pctl_mainresults)}. {p_end}

{synopt:{it:Permutation inference}} When {it:permute()} is used, the reported tables carry the permutation inference (columns {it:se_analytic}, {it:t_perm}, {it:p_perm} and the
interval), the observed analytic tables are preserved hidden as {it:e(perm_analytic_mainresults)}, {it:e(perm_analytic_sgresults_{it:sgcase})} and {it:e(perm_analytic_pwsgs)},
and e() additionally holds, one named column per coefficient of {it:e(b)}: {it:e(perm_p)}, {it:e(perm_ci)} (rows {it:ll} and {it:ul}), {it:e(perm_se_shift)}, {it:e(perm_p_lo)},
{it:e(perm_p_hi)}, {it:e(perm_p_se)}, {it:e(perm_p_ci)}, {it:e(perm_sd)}, {it:e(perm_n)}, {it:e(perm_nullband)}, {it:e(perm_comparator_b)}, {it:e(perm_ase_obs)} and, under
{it:permstepdown}, {it:e(perm_stepdown)} and its family, pre-treatment and subgroup variants; the draw matrices {it:e(permdist_b)} and, under {it:permtstat}, {it:e(permdist_se)}
(up to 10,000 replications); the test statistics {it:e(perm_jointW)}, {it:e(perm_sgW)}, {it:e(perm_advW)} and their observed twins {it:e(perm_jointWobs)}, {it:e(perm_sgWobs)},
{it:e(perm_advWobs)}; run-level scalars and macros {it:e(perm_reps)}, {it:e(perm_completed)}, {it:e(perm_failed)}, {it:e(perm_npossible)}, {it:e(perm_level)},
{it:e(perm_exact)}, {it:e(perm_method)}, {it:e(perm_stat)}, {it:e(perm_pconv)}, {it:e(perm_seed)}, {it:e(perm_rngstate)}, {it:e(perm_analyticse)}, {it:e(perm_absp)},
{it:e(perm_adv)}, {it:e(perm_precaveat)}, {it:e(perm_irebuilt)}, {it:e(perm_blank)}, {it:e(perm_gt_blank)}, {it:e(perm_cfx_blank)}, {it:e(perm_offcentre)},
{it:e(perm_incomplete)} and {it:e(perm_addone)}; and the pooled-estimand scalars {it:e(perm_p_{it:est})}, {it:e(perm_bandlo_{it:est})}, {it:e(perm_bandhi_{it:est})},
{it:e(perm_sd_{it:est})} and, under {it:permtci}, {it:e(perm_tlb_{it:est})} and {it:e(perm_tub_{it:est})}, with {it:e(perm_comparator)} and {it:e(permdist)}. Under the fast-mode
default, {it:e(p)}, {it:e(ci_l)} and {it:e(ci_u)} are also posted for {it:esttab}. Further records: {it:e(perm_npossible_str)}, {it:e(perm_comparator_se)} (under
{it:permtstat}), the hidden {it:e(perm_kernel)} and {it:e(perm_pwsgs_kernel)}, and behind each cfx matrix the hidden {it:e(perm_analytic_{it:cfxname})},
{it:e(permdist_{it:cfxname})}, {it:e(permcmp_{it:cfxname})} and, under {it:permtstat}, {it:e(permdistse_{it:cfxname})} and {it:e(permcmpse_{it:cfxname})}. The ame
step-down tables are {it:e(perm_Cstepdown)} and {it:e(perm_Cstepdown_i)}. {it:e(perm_lb_{it:est})} and {it:e(perm_ub_{it:est})} are deprecated aliases of the band
scalars. {p_end}

{pstd} In addition to these objects returned in e(), if saveplots() is specified, a dataset containing the information used to produce event study plots along with .gph files containing
all graphs generated will be saved to a file with your specified name in Stata's working directory. {p_end}

{pstd} If the option {it:verbose} is specified, the stage 1 underlying regression that {it:wooldid} estimates will also be stored to a model called {it:___wooldidfullmodel}, which can be
accessed using standard {it:estimates restore} type syntax. This is a reserved name that {it:wooldid} overwrites on every run; do not store your own estimates under it. {p_end}

{pstd} {ul:Tips for accessing output and remaking graphs:} Use the {it:matlist} command to print the contents of one of {it:wooldid}'s output matrices. To convert one of {it:wooldid}'s
output matrices into a dataset that you can directly manipulate (e.g., to remake the automatically generated event study graphs), use the {it:svmat} command. For example, the command
{it:svmat e(histogramestimates), names(col)} will create in memory a dataset containing the contents of {it:e(histogramestimates)} (i.e., the cohort i-specific treatment effects),
using as variable names the names of the matrix's columns. The program {help event_plot} can also be used to automatically generate event study plots using the information stored by
{it:wooldid} in {it:e()}, though the user will need to point {it:event_plot} to the name of the results of interest (e.g., for ATT-type event study coefficients, by specifying
{it:stub_lag(ES_att_post#)} and {it:stub_lead(ES_att_pre#)}). If taking this approach, be aware that the {it:event_plot} does not automatically insert 0s into the event study plot
to represent the missing reference period. {p_end}

{pstd} {ul:Compatible Post-Estimation Commands:} {it:estimates save}, {it:estimates store}, {it:estimates restore}, {it:esttab}, and similar commands. {p_end}

{pstd} {ul:Incompatible Post-Estimation Commands:} {it:lincom}, {it:test}, and all other post-estimation commands that rely on receiving a valid variance-covariance matrix from
an estimator. While {it:wooldid} does return a coefficient vector ({it:e(b)}) and a variance-covariance matrix ({it:e(V)}), these are for display purposes only: on an analytic run
the off-diagonal elements are 0 except within a single event-study family, so while {it:esttab} and the like will present correct estimates using these inputs, efforts to combine or
otherwise test different estimates will not, in general, yield valid results ({it:e(wdid_vcvstate)} records what {it:e(V)} is on a given run). The event-study coefficients are
additionally returned with a genuine full variance-covariance matrix in {it:e(b_es_{it:type})} and {it:e(V_es_{it:type})}, on which {it:lincom} and {it:test} are valid once they are copied to matrices and posted
({it:matrix b = e(b_es_{it:type})}, {it:matrix V = e(V_es_{it:type})}, {it:ereturn post b V}). Under {it:bootstrap()} the posted {it:e(V)} is the replication covariance whenever one could be
formed, which {it:e(wdid_vcvstate)} confirms by reading
{it:bootstrap}. {p_end}

{marker comparison}{...}
{title:Comparison to Related Estimators}

{pstd} {ul:Comparing Output to hdidregress twfe:}
The Stata 18 command {help hdidregress} has an option for implementing Wooldridge (2025) style estimates. Its implementation takes the Mundlak approach to estimation
and uses the command {help regress}, in conjunction with {help margins}, to compute estimates. Estimates from {it:hdidregress twfe} with the option {it:notyet}
(i.e., allowing use of not-yet treated cohorts as controls) should be comparable to estimates from {it:wooldid} using the {it:att} estimand-type, a pooled
reference period (i.e., not using {it:esfixedbaseperiod}), unconditional standard errors, and sample weights which (if present) do not vary within units.
Estimates produced by the two will still differ, however, when certain types of panel imbalance are present. For example, when a treatment onset group contains
multiple cohorts and some of those cohorts are not observed in periods when others are observed, {it:wooldid} and {it:hdidregress twfe} will yield different results.
In such cases, {it:wooldid} approaches the situation by presenting the average treatment effect among only the observed cohort-periods. For example, suppose cohorts 1 and 2
are in the same treatment onset group, cohort 1 is observed for the first 3 periods after treatment onset, and cohort 2 is only observed for the first 2. In this case, {it:wooldid}
and {it:hdidregress} will report the same treatment effects for the treatment onset group for the first 2 periods after treatment. The two programs will disagree on the third
effect: {it:wooldid} will report the same 3rd period effect that it (and {it:hdidregress}) would report if cohort 2 were not included in the dataset at all; {it:hdidregress} does not.
Finally, in order to get estimates using controls to
match between the two methods, a {it:wooldid} user will have to manually adjust the controls as needed (e.g., given a control that is constant within cohorts, a
{it:wooldid} user will need to specify i.t#c.control to replicate the behavior of including that control in {it:hdidregress}).  {p_end}

{pstd} {ul:Comparing Output to did_imputation:}
The Borusyak, Jaravel, and Spiess (2023) difference-in-differences imputation estimator is implemented in the package {help did_imputation}. Wooldridge (2025)
proves some equivalence results between their estimator and the Wooldridge estimator. For binary treatments, {it:wooldid} will yield the same estimates as
{it:did_imputation} given that one uses the {it:att} estimand-type, uses a pooled reference period (i.e., not using {it:esfixedbaseperiod}), requests no
pre-treatment horizons ({it:esprelength(0)}), does not use unconditional standard errors, and does not use controls. The {it:esprelength(0)} condition matters:
the pooled reference period is formed from the horizons the run asks for, so requesting pre-treatment horizons removes them from the reference period and shifts every
post-treatment estimate. With controls, {it:wooldid} and {it:did_imputation} will still yield the same results provided
that each cohort-time period cell contains only 1 observation. In this 'fine cohorts' case, the inclusion of i-t cell indicators ensures that the included controls
have their coefficients determined solely from the same portion of the sample as is used to estimate the regression in the BJS estimator. However, outside of special
cases, results with controls will not match if cohorts are coarse (contain multiple observations per i-t cell). To assist in this case, the {it:wooldid} option
{it:coarsecohortcontrols()}
will add to the {it:wooldid} regression your controls plus interactions between your controls and the i-t cell indicators associated with periods for which treatment
effects are requested. While the {it:coarsecohortcontrols()} option will not yield the same estimates as {it:did_imputation}, it arguably reflects an approach closer to
the spirit of BJS in that it increases the degree to which the effect of the controls is independent between the control/reference portion of the sample and
the portion of the sample where treatment effects are of interest. {p_end}

{pstd} {ul:Comparing Output to csdid and xthdidregress:}
The {it:gteffects} option returns the cohort-by-relative-period effects underlying {it:wooldid}'s aggregates, which are comparable to the group-time tables of
{help csdid} and {help xthdidregress} when {it:i} indexes treatment cohorts. {p_end}

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

{phang}4) Event study with 3 pre-treatment and 3 post-treatment estimates using a pooled reference period and generating event study plots, as well as testing the hypothesis
that all pre-treatment effects are 0 and fitting the adversarial line (whose p-value requires {it:esfixedbaseperiod}, so here the line is reported without one): {p_end}

{p 10 10}   {cmd:wooldid y i t ttre, espre(3) espost(3)  makeplots jointtests adversarial}

{phang}5) Event study with 3 pre-treatment and 3 post-treatment event study estimates, calculating all estimates relative a fixed reference period (fixed to the
period 2 periods before treatment):{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, espre(3) espost(3) esfixedbaseperiod esrelativeto(-2)}

{phang}6) Diff-in-diff with separate treatment effects calculated by treatment onset cohort:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, subgroups(ttre)}

{phang}7) Event study with a fixed reference period, treating treatment as anticipated one period before its recorded onset, with Roth (2022) pre-trend power diagnostics
at 50 and 80 percent power:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)}

{phang}8) Difference-in-differences with a treatment intensity variable (one that is 0 outside the treatment period), estimated using a stepwise-linear function and a lattice with 5 quantiles
(formed based on the treated group's treatment period), and generating continuous treatment effects plots for the post treatment period:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(intensity) contreatcontrols(untreatedZero) makecfxplot lattice(baseTreatedPost 5)}

{phang}9) Difference-in-differences with a continuous treatment variable that varies within all cohorts, with a lattice containing 5 quantiles (formed separately by treatment group and
by period), controlling for the continuous treatment variable by treatment group and by time period, with continuous treatment effects plots for the pre-treatment period (which
require an event-study window, here via {it:espre()}):{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, contreat(contreatvar) contreatcontrols(byTreated byTime) espre(3) espost(3) makecfxplot lattice(byTreatedByTime 5) cfxplottypes(pre)}

{phang}10) Event study followed by a Rambachan and Roth (2023) HonestDiD sensitivity analysis, using the event-study variance-covariance matrix {it:wooldid} returns (requires the
{stata ssc install honestdid:honestdid} package, and a fixed reference period at the default {it:esrelativeto(-1)}; see {it:e(b_es)} under Working With Output):{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att espre(4) espost(4) esfixedbaseperiod}

{p 10 10}   {cmd:matrix b = e(b_es_att)}

{p 10 10}   {cmd:matrix V = e(V_es_att)}

{p 10 10}   {cmd:honestdid, numpre(`=e(numpre_att)') b(b) vcov(V) delta(rm) mvec(0.5(0.5)2) parallel(0)}

{p 8 8}{cmd:parallel(0)} keeps {it:honestdid} from starting extra Stata processes, which it does by default wherever the {it:parallel} package is installed. A robust confidence
set that comes back empty is a finding about the data (the pre-trend evidence rules out every value at that {it:M}), not a failure of the handoff; extend {it:mvec()} to larger {it:M}
to see where the set becomes non-empty.{p_end}

{phang}11) Fast exploratory event study (point estimates only), also returning the cohort-by-relative-period building-block table:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) espre(2) espost(2) fastnose gteffects}

{phang}12) Randomization inference: permute treatment-onset dates across clusters, holding onset-group counts fixed, with 500 replications; every reported p-value and interval is
then permutation-based, and the run is a fast (no standard error) run unless {it:permtstat} or {it:permanalyticse} is added:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(500) rseed(12345)}

{phang}13) Cluster bootstrap standard errors with 500 replications, reproducible, with joint pre-trend tests formed from the replications' own test statistics ({it:boottstat}),
the route to a calibrated joint test when clustering at or above the cohort level:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) jointtests bootstrap(standard) reps(500) rseed(42) boottstat}

{phang}14) Bayesian (Dirichlet reweighting) bootstrap, which leaves every cluster in every replication and so never loses the comparison group to an unlucky resample:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) bootstrap(bayes) reps(500) rseed(42)}

{marker bugs}{...}

{title:Bugs and Assistance}

{pstd}{ul:Known Issues:} {p_end}

{p 10 10}   With cohort-specific time trends, do not manually specify time-trend variables for inclusion in the fixed effects (i.e., in the {it:absorb()} statement when using
{it:reghdfe}); let {it:timetrends} build them. Estimates are better conditioned when {it:esfixedbaseperiod} is also specified.

{p 10 10}   Estimation with continuous treatment controls sometimes faces numerical stability issues. For more consistent results, specify either a higher {it:regtol()} or
            the option {it:ccnoabsorb}, which keeps the continuous treatment controls out of {it:reghdfe}'s {it:absorb()} statement. A distinct case can look similar: when a
            treated cohort carries almost no variation in the dose across its untreated or reference-period observations, that cohort's baseline dose slope is nearly
            unidentified and the reported average marginal effect is unreliable in a way that tightening {it:regtol()} does not repair. {it:wooldid} cautions when it sees
            this and refuses when the variation is smaller still; the remedy is to examine or drop the affected cohort(s).

{p 10 10}   The automatic tests for model identification that {it:wooldid} implements are not fully comprehensive: {it:wooldid} will almost always succeed in identifying
when a treatment effect has been partialled out by a fixed effect in the binary treatment case and will often succeed in doing so when a continuous treatment variable is specified
without a lattice variable or any polynomial terms. However, when a continuous treatment is specified with a lattice and/or a polynomial term, and especially in the presence of an
unbalanced panel, {it:wooldid} will be unlikely to detect problems related to important treatment effects being partialled out by controls or fixed effects. As such, one should
be very careful when specifying one's fixed effects and {it:contreatcontrols} when using a continuous treatment.

{p 10 10}   Under various circumstances, {it:wooldid} cannot estimate standard errors for all requested objects. If standard errors are identified for some but not all
reported estimates, {it:wooldid} still returns {it:e(V)}, with the affected estimates carrying a placeholder entry that displays as a missing standard error; code that reads
{it:e(V)} directly should treat negative diagonal entries as not identified. If no reported estimate has an identified standard error, {it:wooldid} returns point estimates
without {it:e(V)}, which can cause the functionality of certain post-estimation tools to break. In these situations, {it:wooldid} prints a caution.

{pstd}If you encounter a bug - be it the program crashing or just the program delivering a weird result - please reach out to me about it. This program is still
under development, so assistance finding and correcting bugs is much appreciated. I am also happy to help and offer guidance on how to use this program. A longer
guide to using {it:wooldid}, written for AI agents and advanced users, is at {browse "https://github.com/thegland/wooldid/blob/main/wooldidinstall/wooldid_agent_guide.md"}. {p_end}

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
