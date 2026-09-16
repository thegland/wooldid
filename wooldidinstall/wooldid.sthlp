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
{viewerjumpto "Stored Results" "wooldid##storedresults"}{...}
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

{phang}{ul:The {it:if} qualifier:} {cmd:wooldid} keeps the observations for which the {it:if} expression evaluates to exactly 1, where Stata's own {it:if} qualifier keeps every nonzero value, missing included. The two agree for comparisons such as {cmd:if x > 5}, which return 0 or 1, and differ for a bare flag variable holding other values: {cmd:if flag} excludes the rows where {cmd:flag} is 2 (or missing), which Stata would keep. The rule is deliberate - it is what stops missing and out-of-range values entering the sample - and it applies identically to the estimation sample, to the randomization design under {help wooldid##standarderrors:permute()} and to the bootstrap frame, so the three always agree. {cmd:wooldid} prints a Caution naming the count whenever the expression is nonzero but not exactly 1 somewhere, and those observations are included in {it:e(N_dropped)}. Write {cmd:if flag != 0} for Stata's behaviour. {p_end}

{synoptset 22 tabbed}{...}
{synopt : {it:y}}Outcome variable. {p_end}
{synopt : {it:i}}Cohort identifier; must consist of integers greater than 0. Cohorts can be fine (i = individual units) or coarse (i = collections of individual units), but should be no coarser
than to consist of treatment onset cohorts. Coarser cohort definitions tend to enable faster estimation. For a typical run, prefer treatment onset
 cohorts as {it:i} and, on the default analytic path, cluster below that level (see {help wooldid##standarderrors:What level to cluster at});
 {cmd:permute()} is the one engine that still ties the two choices together. It draws one onset per cluster, so it needs {cmd:cluster()} to nest {it:i} unless {cmd:permonsetcohorts} asserts that {it:i} {ul:is} the treatment onset grouping, in which case each replication rebuilds {it:i} from the drawn onsets and a finer {cmd:cluster()} is admissible.
 Choose {it:i} by the treatment heterogeneity the estimands should average over, not by the level at which you want fixed effects - those are
 specified separately in {opt fe()}, so a coarse {it:i} combines with unit-level ones - and go finer only where that heterogeneity, or a
 {cmd:subgroups()} variable varying inside a cohort, requires it. {p_end}
{synopt : {it:t}}Time period identifier; must be numeric and on the same grid as {it:ttre} ({it:t} minus {it:ttre} must be an integer wherever {it:ttre} is nonmissing - see the {it:ttre} entry). Integer values are recommended, but the value itself only labels periods on the default estimation path: zero, negative, and non-consecutive values (gaps) are all fine. With {opt uncond:itionalse}, {opt customvce()}, or {opt inter:activecontrols()}, {it:t} enters the model as a native Stata factor variable and must then consist of non-negative integers no larger than 2,147,483,647 (Stata's factor-variable cap): recode large calendar codes such as %tc datetimes to small integers first.{p_end}
{synopt : {it:ttre}}Variable containing cohort-specific dates of treatment onset; missing indicates never-treated cohorts, cohorts treated before the sample begins are dropped. Binary designs require a comparison group: some units must be never-treated or first treated after the estimation window (an informative error is raised otherwise; see {it:oostreatedcontrols}). Store {it:ttre} (and {it:t}) as {cmd:double} when onset codes are large: 8-digit calendar dates exceed float precision, so a float-stored onset variable is rounded before {it:wooldid} sees it, silently merging adjacent cohorts. {it:wooldid} prints a caution when {it:i}, {it:t} or {it:ttre} is stored as float with values above 16,777,216 (2^24, beyond which float cannot hold consecutive integers); the values must be regenerated as double from their source, since recasting after the fact does not recover them. {it:ttre} must take a single value (or be missing on every row) within each group defined by {it:i} over the estimation sample: {it:wooldid} builds relative time row by row from {it:ttre}, so an {it:i} group carrying different onsets on different rows - including an onset recorded only from the onset period onward, with missing values pre-treatment - would be silently split across cohorts. {it:wooldid} verifies this and exits with an informative error if it is violated (bypass with {it:safety(off)}). Onset times must also sit on the same grid as {it:t} - {it:t} minus {it:ttre} must be an integer wherever {it:ttre} is nonmissing (integer spacing is not required: biennial or five-year panels are fine for estimation, but see {it:esprelength()} - the event-study window is counted in relative-period slots, so such a panel realizes only every second or fifth horizon of a requested window). Fractional offsets previously survived to the estimation stage, where they either failed with opaque errors or silently mislabeled the onset period as pre-treatment; {it:wooldid} now rejects them at parse time (bypass with {it:safety(off)}). If every offending onset falls outside the usable estimation window (pre-sample, or beyond the last observed period), the rejection softens to a Caution, since such onsets are dropped or treated as never-treated downstream anyway. {p_end}

{phang} {ul:Options by Purpose:} {p_end}
{synoptset 1 tabbed }{...}
{synopt : {help wooldid##specifyestimand:Specifying the Estimand(s):}} {opt att} {opt att_it} {opt att_i} {opt att_itime} {opt semi:elasticity} {p_end}

{synopt : {help wooldid##standarderrors:Standard Errors:}} {opt clus:ter(variable(s)[numerical])}  {opt rob:ust} {opt uncond:itionalse} {opt customvce(customstring)} {opt imp:values} {opt fast:nose} {opt fastcollapse} {opt fullprecisionslow} {opt simult:aneousbands}{p_end}
{synopt : {space 4}{it:shared by both replication engines:}} {opt level(#)} {opt dots(int)} {opt rseed(#)}{p_end}
{synopt : {space 4}{it:permutation only:}} {opt perm:ute(method)} {opt permreps(int)} {opt permex:act} {opt permstr:ata(variable)} {opt permpr:obs(spec)} {opt permts:tat} {opt permt:ci} {opt permsa:ve(name)} {opt perman:alyticse} {opt permab:sp} {opt permsi:de(left|right)} {opt permraw:p} {opt permstep:down} {opt permonset:cohorts}{p_end}
{synopt : {space 4}{it:bootstrap only:}} {opt boot:strap(mode)} {opt reps(int)} {opt ties} {opt bootalp:ha(spec)} {opt bootd:f(spec)} {opt bootts:tat} {opt bootdisp:lay(spec)} {opt bootsa:ve(name)}

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
using whatever weights (if any) are used in the main regression command, multiplied by the {it:customweightmultiplier()} variable if one is specified - that multiplier enters the margins weights only and is deliberately not folded into the regression weight.{p_end}

{phang}{opt att_it}: Estimate simple averages of the i-t cell-specific treatment effects (ATT_it). The treatment effect in each cohort-period is calculated via the margins calculation
using any weights that may apply to the observations falling in each i-t cell. The final ATT_it estimate is the simple, unweighted mean of all of the cohort-period-specific estimates.
This is done by reweighting the data used in the margins calculation (after estimation of the underlying regression) by w/W_it where w is an observation's sample weight and W_it is the
total sample weight associated with all observations appearing in a given i-t cell. If no weights are specified, w is taken to have a value of 1 for all observations. If cohorts are
fine and thus each i-t cell contains only one observation, this is equivalent to ignoring any specified sample weights when performing margins calculations.   {p_end}

{phang}{opt att_i}: Estimate the simple average of cohort-specific treatment effects (ATT_i). The treatment effect for each cohort is calculated using the margins calculation and any sample
weights. The final ATT_i estimate is the simple, unweighted mean of the cohort-specific treatment effects. ATT_i effects are computed by performing the margins calculation, but reweighting
observations by either w/W_ipre or w/W_ipost, where w is an observation's sample weight, W_ipre is the total sample weight associated with cohort-i in the pre-treatment period
(but outside the reference period, and, under {it:timetrends}, outside the period immediately preceding it, which {it:timetrends} normalizes to zero as well), and W_ipost is the total sample weight associated with cohort-i in the post-treatment period. Each observation is then reweighted based on whether it
falls in the pre- or post-treatment period. If weights are not used, w is taken to have a value of 1 for all observations.{p_end}

{phang}{opt att_itime}: Estimate the simple average of cohort-specific simple average treatment effects (ATT_itime). Treatment effect calculation here begins like for {it:att_it}, obtaining
treatment effects for individual cohort-periods of interest using weights. Then, cohort-specific effects are obtained by taking the simple mean of the cohort-period-specific effects.
The overall presented ATT_itime estimate is then equal to the simple, unweighted mean of these cohort-specific effects. ATT_itime effects are calculated by performing the margins calculation,
but reweighting observations by w/(W_it * T_iprepost), where w is an observation's sample weight, W_it is the amount of sample weight falling in a given i-t cell, and T_iprepost contains
(depending on the time period) either the number of post-treatment time periods in which a cohort-i is observed or the number of pre-treatment, non-reference period time periods in which a
cohort i is observed. If weights are not used, w is taken to have a value of 1 for all observations.{p_end}

{pmore}Each of these four options names a weighting concept rather than a single number. Under {it:contreat()} every requested estimand is therefore reported twice on that same weighting: once as the effect of the binary onset ({it:att}, {it:att_it}, {it:att_i}, {it:att_itime}), and once as the average marginal effect of the dose ({it:ame}, {it:ame_it}, {it:ame_i}, {it:ame_itime}). Both families are displayed, and both are stored under the same names with {it:att} replaced by {it:ame}. {p_end}

{pmore}{ul:Options that treat the cell-weighted estimands differently:} because {it:att_it}, {it:att_i} and {it:att_itime} give each treatment cell equal weight rather than each observation, several options restrict or reinterpret them. {it:unconditionalse} is refused with an error under default safety ({it:safety(off)} proceeds). {it:bootstrap(standard)} re-identifies the duplicated copies of a resampled cluster for these estimands, so that a twice-drawn cluster contributes its cells twice rather than pooling them into one. Where {it:cluster()} nests {it:i} the copies become distinct units; where it is finer than {it:i}, or cuts across it, the copies stay inside their cohort and the replication's own standard error clusters on the drawn cluster identity. With no {it:cluster()} at all the resampling is observation-level. {it:bootstrap(bayes)} is refused for them under default safety, because cluster reweighting cannot resample treatment cells and the standard errors would omit the cell-sampling variance component ({it:safety(off)} proceeds with a caution). A {it:subgroups()} variable that splits a unit across subgroups over time is refused whatever the estimand, and under {it:safety(off)} such a unit enters each subgroup weighted by its share of sample weight there, which is not the equal unit weighting {it:att_i} and {it:att_itime} otherwise imply. Finally, {it:makecfxplots} draws its continuous-treatment plots, and returns the {it:e(cfxatt_*)} and {it:e(cfxame_*)} matrices, for the {it:att} weighting only: for {it:att_it}, {it:att_i} and {it:att_itime} it prints a caution and produces neither the plots nor the matrices unless {it:safety(off)} is specified ({it:contreatwithin} with {it:contreat()} sets safety off internally and so has the same effect). {p_end}

{phang}{opt semi:elasticity}: Compute treatment effect estimates in terms of semi-elasticities (% change in y / unit change in x). This is achieved by adjusting the margins calculation
to present the average percentage change, rather than the average level difference, between the two sets of predicted outcomes. This is implemented using the
{it:eydx} option of the {help margins} command. Note that for estimation to succeed, this option tends to require large samples and a y-variable that is not sparsely distributed.
{ul:Under {it:poisson} without {it:poisexpresults}} the reporting scale is already the linear index {it:xb} = ln E[y|x], on which a plain marginal effect {ul:is} the semi-elasticity of E[y];
{it:margins}' {it:eydx} there would take a further log of the linear index and report a number that depends on the units of {it:y}. On that combination {it:wooldid} computes {it:dydx} on the
{it:xb} scale instead and prints a caution saying so, so the option changes the table header rather than the estimates - which are exactly a plain {it:poisson} run's. Use {it:poisexpresults}
if you want the semi-elasticity of the exponentiated prediction. Note also that {it:eydx} divides by the prediction, which makes a run ineligible for the cell collapse: it is evaluated on the full data, silently unless {it:verbose} is on, and an explicitly typed {it:fastcollapse} is refused with an error (see {help wooldid##standarderrors:fastcollapse}). The remapped {it:poisson} case just described is the exception, since what it computes there is {it:dydx}. {p_end}


{marker standarderrors}{...}
{dlgtab:Standard Errors}

{phang2}{ul:What level to cluster at.} The intended practice is to cluster {ul:below} the level of {it:i} - at a level whose units are nested inside cohorts, so that a treated cohort-by-period cell contains several clusters. At or above the level of {it:i} the cluster-robust covariance stops seeing the treated cells' own variation, and two things follow that the table does not announce: the pointwise standard errors on the {it:att} rows are understated (measured at roughly 0.66 to 0.75 of the reference on the diagnostic fixtures), and the joint pre-trend tests are oversized, which is why those carry a caution of their own (see {cmd:jointtests}). There is no plug-in correction - at that level the analytic covariance is the wrong object rather than a mis-scaled one - so the remedies are a finer {cmd:cluster()}, or {cmd:bootstrap()} or {cmd:permute()}, which build their inference from replications instead. Those two remedies now compose: a {cmd:cluster()} finer than {it:i} is accepted under {cmd:bootstrap(standard)} with the cell-weighted estimands and under {cmd:boottstat}, and under {cmd:permute()} with {cmd:permonsetcohorts}. {p_end}

{phang}{opt clus:ter(clusvars)}: Variable(s) specifying clustering level. Several variables give multi-way clustering on the default estimation path; {it:unconditionalse}, {it:impvalues}, {cmd:permute()} and {cmd:bootstrap()} each require a single cluster variable and refuse a multi-way {it:cluster()} with an error. Clustered standard errors are computed analytically on the default estimation path; under {cmd:bootstrap()} the reported standard errors are replication quantities, under {cmd:permute()} without {cmd:permtstat} or {cmd:permanalyticse} no analytic standard error is computed at all and the reported inference is the permutation distribution's, and under {it:impvalues} the reported standard error is the cross-cluster Ibragimov-Muller one. A cluster variable that leaves fewer than two clusters in the estimation sample stops the run with an error, since cluster-robust inference needs at least two of them; the count is taken after {it:wooldid} has formed its sample, so a variable that looks adequate in the raw data can still fail here. The message names the three ways out: a finer cluster variable, {it:robust}, or {it:fastnose}, which asks {cmd:margins} for no variance at all and is exempt from the check. The same gate stops a run that is left with no residual degrees of freedom whatever. By 
default,
regression estimation using this option will underlying rely on either {help reghdfe} or, for poisson, {help ppmlhdfe}. {p_end}

{phang}{opt rob:ust}: Heteroskedasticity robust standard errors; relevant only if cluster() is not specified. By default, regression estimation using this
option will underlyingly rely on either {it:reghdfe} or, for poisson, {it:ppmlhdfe}. {p_end}

{pmore}Note on degrees of freedom: the default estimation path and the {cmd:regress}-based path used by {it:customvce()} and {it:unconditionalse} apply different rules for fixed effects nested within the clustering variable, so the same specification can return slightly different standard errors depending on which path it takes. The divergence is conservative on the {cmd:regress} path rather than anti-conservative, and it does not generalise to the non-cluster variance strings, which is why wooldid documents the difference rather than rescaling one path onto the other. {p_end}

{phang}{opt uncond:itionalse}: Compute standard errors that account for sampling variation in covariates. With {cmd:cluster()}, a single cluster variable is supported: the variance comes from {it:margins, vce(unconditional)}, which accepts one. This option also forces {it:wooldid} to abandon
high dimensional fixed effect estimation packages for the underlying two-way fixed effect regression in favor of the commands {it:regress} and {it:poisson}.
As such, this option should be slower to use on average. This option requires specifying either robust or clustered standard errors (a robust or clustered vce supplied via {it:customvce()} also satisfies this requirement, as do pweights, which imply robust standard errors). Under the default {it:safety()} setting it is refused with {it:att_i}, {it:att_it} and {it:att_itime}, whose unconditional standard errors are unreliable ({it:safety(off)} proceeds anyway), and it cautions with {it:aweight}s and with {it:fe()} or {it:ccc_absorb()}, which the {cmd:regress} path cannot absorb - specify those as {it:controls()} with your own {it:i.} prefixes. It is refused outright with {it:impvalues}, {it:fastnose}, {cmd:permute()} and {cmd:bootstrap()}, and it turns the default cell collapse off, since {it:margins} checks the collapsed row count against {it:e(N)} when building the unconditional variance (an explicit {opt fastcollapse} alongside it is then an error). {p_end}

{phang}{opt customvce()}: Forces {it:wooldid} to switch to using either the {help regress} or {help poisson} command for underlying regression estimation,
and then passes the user specified text without modification into the {it:vce()} option for these commands. This allows specifying a broader set of
standard error calculation methods than allowed by the default high dimensional fixed effects estimation packages used by {it:wooldid}. Not compatible with {it:robust} or
{it:cluster()}: specify the clustering or robustness desired inside the {it:vce()} text itself. Also refused with {it:impvalues} and with {cmd:bootstrap()}, whose reported standard errors would be replication quantities rather than the custom VCE's. It cannot be combined with {cmd:permute()} either, though by construction rather than by a screen of its own: {cmd:permute()} requires a single {it:cluster()} variable to randomize onsets over, and {it:customvce()} refuses {it:cluster()}, so the pair cannot be written down. Specifying the clustering inside the {it:vce()} text does not get you there, because it is the {it:cluster()} variable itself that the permutation engine needs. It turns the default cell collapse off pending measurement, so an explicit {opt fastcollapse} alongside it is refused. With {it:fe()} or {it:ccc_absorb()} {it:wooldid} cautions, since the {cmd:regress} path it switches to cannot absorb them - specify those as {it:controls()} with your own {it:i.} prefixes. {p_end}

{phang}{opt imp:values}: Implements Ibragimov and Muller (2010) style cluster-robust inference for all estimates. The idea behind the Ibragimov and Muller approach is to compute
cluster-specific treatment effect estimates, present their average as one's final treatment effect estimate, and conduct inference by t-testing the cluster-specific estimates against 0.
Using IM-style inference requires specifying a single cluster variable and is not compatible with the {it:unconditionalse} option. IM-style inference changes {it:wooldid}'s
target estimand, and it does so in every run rather than in special cases: the procedure averages the cluster-specific estimates with {ul:equal} weight, so what it reports is the
cluster-equally-weighted mean, whatever estimand was requested and whatever weights that estimand ordinarily carries. Clustering at the cohort level is where this becomes
visible - {it:att} then yields the same effect estimates as {it:att_i} - and {it:cluster()} at the level of {it:i} with {it:att_i} is the pairing that is coherent by
construction, since the cluster-equally-weighted mean is what IM computes in any case. The equal weighting is not an implementation choice that could be reweighted away: it is
what earns the {it:q}-1 degrees of freedom below, so trading it for the requested estimand's weights would cost the inference that motivates the method. Note that {it:wooldid} does not adjust the
weights used for computing different estimands when IM-style inference is requested: in general, IM-style inference will yield incorrect results for non-att estimands if clusters cross-cut cohorts
or (worse) cohort-time period cells. IM-style inference is not currently allowed with the use of subgroups, the option {it:jointtests}, or the option {it:makecfxplots}, nor with {it:adversarial}, {it:customvce()}, {it:fastnose}, {it:pretrendpower()}, {cmd:permute()} or {cmd:bootstrap()}: each of these combinations is refused with an error. Each reported row's test statistic, p-value and confidence interval are on the t distribution with {it:q}-1 degrees of freedom, where {it:q} is the number of clusters that contributed a non-missing estimate to that row, so {it:q} can differ from row to row. {it:e(df_r)} is set to the SMALLEST of those degrees of freedom, so a {cmd:test}, {cmd:lincom} or {cmd:esttab} run by hand on {it:e(V)} is conservative relative to the printed table rather than narrower than it; the range is stored in {it:e(imdf_min)} and {it:e(imdf_max)} whether or not the degrees of freedom are common, and when they are not common the run additionally says so under the table. The stage-1 regression's residual degrees of freedom - which the {it:gteffects} cells still use - are kept as hidden {it:e(df_r_analytic)}. Every IM {it:margins} call runs with {it:nose}, so the event-study b/V capture never fires: {it:e(b_es_)}, {it:e(V_es_)} and {it:e(numpre_)} are not returned and the {it:honestdid} and {it:pretrends} handoffs have no inputs. {it:e(wdid_vcvstate)} therefore reads {it:analytic-diagonal} rather than {it:analytic-blocks} - the within-family event-study covariance blocks an ordinary run splices into {it:e(V)} are absent here, so every off-diagonal of {it:e(V)} is zero. {p_end}

{phang}{opt fast:nose}: Fast exploratory mode. Every {it:margins} call runs with {it:nose}, so no standard errors, test statistics, p-values, or confidence intervals are computed
for any estimand: only point estimates are produced, the corresponding cells of the results tables display as missing (.), and the output carries the line "SEs and related materials
are not estimated in fast mode." Intended for specification search on large panels, where the delta-method variance calculation inside {it:margins} dominates run time (measured
speedups of roughly 2-4x on typical specifications, and far larger on poisson event-study runs, where fast mode completes in seconds on configurations that otherwise run for tens of
minutes); re-run the chosen specification without {it:fastnose} before reporting results. Under {it:fastnose}: {it:e(V)} is not posted ({it:e(properties)} is {it:b}, and the disclosure marker {it:e(wdid_vcvstate)} reads {it:absent}), so {it:test},
{it:lincom}, {it:_se[]}, {it:estat vce} and anything else needing a variance-covariance matrix will refuse. {it:esttab} and {it:estout} do {ul:not} refuse, and under a bare {it:fastnose} they print a coefficient with an empty inference row beneath their usual footer, so read the table rather than trusting the footer. Under the {cmd:permute()} fast-mode default they do better than that: the run posts {it:e(p)}, {it:e(ci_l)} and {it:e(ci_u)} under the names {it:estout} reads {ul:ahead} of anything derived from {it:e(V)}, so {cmd:esttab, p} and {cmd:esttab, ci} carry the run's real permutation p-values and shift intervals, with truthful footers, and a plain {cmd:esttab} tables the coefficients. No {it:e(se)} is posted there and none is invented: a permutation run has no standard error to report; {it:e(V_es_)} is not returned, so
{it:honestdid} and Roth-style diagnostics have no covariance to work from, while {it:e(b_es_)} and {it:e(numpre_)} {ul:are} returned - the estimates are meaningful even where their uncertainty is unavailable, and a Caution per family says so (under the {cmd:permute()} fast-mode default none of the three is returned); no joint tests of any kind are computed (including the subgroup homogeneity tests); the one exception is the subgroup-by-horizon margins call, which keeps its variance so that the estimability screen described under {it:subgroups()} still drops a non-estimable horizon (fast mode therefore reports the same event-study rows as the full run, at a small cost per subgroup); the continuous-treatment-effects plot
matrices carry only the point-estimate column, exactly as under {it:cfxplotnose}; and event-study plots are drawn without confidence bands. {cmd:simultaneousbands} is accepted rather than refused, but no band is built - the analytic engine needs the event-study covariance fast mode does not compute - and a caution says so; the pointwise interval columns are unaffected and no {it:e(esband_*)} results are returned. Not compatible with
{it:unconditionalse}, {it:impvalues}, {it:pretrendpower()} or {cmd:bootstrap()} (each is refused with an error), nor with {it:permtstat} or {it:permtci}, which randomize on t statistics fast mode does not compute. {it:adversarial} and {it:jointtests} are likewise refused, {ul:except} under the {cmd:permute()} fast-mode default described next (no {cmd:permtstat}, no {cmd:permanalyticse}), where fast mode is already in force and naming it is accepted rather than an error. Under {cmd:permute()} without {cmd:permtstat} the observed run takes this fast path by default - the permutation inference does not consume the analytic standard errors - and {cmd:permanalyticse} restores them; the adversarial and joint-test surfaces are then permutation quantities built after the replications, not the analytic ones fast mode lacks. Point estimates can differ from a full run in
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
applied, 0 when it was not), beside the scalar {it:e(wdid_fcncell)} - the collapsed cell count - and the macro {it:e(wdid_fcskip)}, which names the reason when an eligible run did not collapse after all; {it:verbose} reports the decision as the run proceeds. {p_end}

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
{opt fullprecisionslow} - which restores the summation order, but not a computation that has itself changed since; see its entry. {p_end}

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
not linear in the coefficients. ({it:poisson} with {it:contreatelasticitytype(dyex)} left this list on 2026-09-06: {it:dyex} now uses the analytic Jacobian under {it:poisson} as {it:dydx} does, and the collapse reproduces the full-data {it:ame} rows there.)
{ul:The read splits a cell on something the cell key does not contain:}
{it:lattice()} and {it:makecfxplots} (with the bins of {it:lattice()} or {it:cfxlattice()}), which bin the dose; {it:cfxlattice()} on its own is ignored and does not affect the collapse. ({it:impvalues}, which computes effects over the clusters, left this list on 2026-09-06: the cluster variable joins the cell key under {it:impvalues}, which makes the per-cluster reads exact - measured at 1e-12 both for a cluster nesting {it:i} and for one crossed with it - at the cost of more cells only when the cluster is crossed with {it:i}.)
{ul:{it:margins} objects to the collapsed row count itself:} {it:unconditionalse}, whose variance is built from observation-level influence functions and whose row count
{it:margins} checks against {it:e(N)}. {ul:Unverified rather than known-broken:} {it:customvce()}, refused pending measurement. {p_end}

{phang2}Everything else is eligible, including several combinations that are easy to assume are not. {it:contreat()} is eligible on the default {it:dydx} scale and under
{it:contreatelasticitytype(dyex)}, with every flavour of {it:contreatcontrols()} and with {it:ccnoabsorb}. Poisson estimation on the default (log) reporting scale is eligible,
including with {it:contreat()} and with {it:contreatelasticitytype(dyex)} - and, since 2026-09-07, with {it:semielasticity} and {it:contreatelasticitytype(eyex)}/{it:(eydx)}, which on the
linear-index scale are computed as {it:dydx}/{it:dyex} and so are linear in the coefficients (they remain ineligible off the poisson path and under {it:poisexpresults}). So are {it:fastnose}, {it:jointtests}, {it:adversarial} and {it:altmethod_adversarial}, {it:esfixedbaseperiod} and {it:esrelativeto()},
{it:pretrendpower()}, {it:histogramcohorteffects}, {it:summarystats}, {it:gteffects} (which {it:wooldid} refuses alongside {it:contreat()} for a separate reason), {it:makeplots}
and {it:saveplots()}, and the replications of {cmd:permute()} and {cmd:bootstrap()}, which are evaluated collapsed like any other run. {p_end}

{phang2}{ul:{bf:subgroups()} is eligible, but screened, for a subtler reason worth understanding.} When {it:margins} judges a subgroup-by-horizon effect not estimable it does not
error - it computes the number anyway and flags it only by marking the name it returns, and that flag rests on a rank decision made at the last digit. wooldid's subgroup code is
the one place that reads the flag, and it does so three separate ways: to drop the affected horizon with a caution, to skip the adversarial pre-trend test, and inside the subgroup
joint tests. A fourth exact-zero branch is the adversarial line test's own rank truncation, at 1e-12 of the largest eigenvalue of the event-study covariance, which runs on every event study (subgroup or not) and is three orders tighter than the collapse screen's tolerance - so a residue-level direction is dropped identically collapsed or not. So it is the one place where the collapse's last-digit change in summation order can change {bf:what} is reported rather than its low-order digits: a horizon that
should be dropped can come back as a meaningless near-zero, or a reported joint test can turn to missing. This can only happen where the variance matrix has a direction that is
singular in truth - a subgroup effect, or a contrast between subgroups, whose true variance is zero, as when a fixed effect absorbs every cell ({it:fe()} on a cohort-by-period
variable), the fit is exact, or the clustered variance is rank-deficient - so that the number {it:margins} computes there is cancellation residue rather than an estimate. wooldid
therefore screens each binary-treatment subgroup {it:margins} call after the fact (the continuous-treatment subgroup reads under {cmd:contreat()} are not screened; see the design record), inspecting the returned variance matrix for a not-estimable marker, a subgroup effect or a between-subgroup
contrast whose variance is below 1e-9 of the largest, or a contrast direction singular at that tolerance (the smallest genuine subgroup contrast measured across the test inventory
is 3.5e-5 and the residue is around 1e-16, so the band has ample room on both sides). The screen is a few matrix operations per call. If it finds nothing, the collapsed result
stands. If it finds a singular direction, then under the default safety setting wooldid restores the full data at that point and evaluates that call and every remaining subgroup
call on it, so everything downstream of the flag is exactly what {opt fullprecisionslow} would have produced (the main results, and any subgroup horizons already read before the
flag, were computed collapsed and agree with their full-data values to the last digit as described above); {it:verbose} reports when this happens. Under {it:safety(off)} wooldid
instead keeps the collapsed rows, as requested, and prints a caution naming what it found, since in that situation which subgroup horizons and joint tests are reported can depend on
the summation order. The outcome is recorded in the hidden macro {it:e(wdid_fcguard)}: empty when nothing was flagged, {it:fallback}, or {it:caution}. An explicit
{it:fastcollapse} is a request the screen may set aside in the same way; and a fallback carries into the replications of both {cmd:permute()} and {cmd:bootstrap()}, which are otherwise evaluated collapsed. {p_end}

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
{it:margins} result name, and 123 non-subgroup comparisons on aggregating data showed no structural difference at all. (The adversarial rank truncation just described applies to non-subgroup event studies too, but it is a fixed rule of the test itself, not a read of a {it:margins} result name.) {p_end}

{phang2}{ul:Memory.} The full data is held aside while the collapsed window runs, so peak memory includes an extra copy of your dataset. This is immaterial for typical panels but
worth knowing if you are already near your memory limit, particularly if you have issued your own {help preserve} before calling {it:wooldid}; {opt fullprecisionslow} removes
the extra copy. {p_end}

{phang}{opt fullprecisionslow}: Turn the default cell collapse off, so that every {it:margins} call is evaluated on the full estimation sample as it was before the collapse
existed. Results are then bit-identical to those of earlier {it:wooldid} versions given the same data and options wherever the estimator itself is unchanged - the option restores the summation order, not a computation that has since moved: under {it:poisson} without {it:poisexpresults}, {it:semielasticity} and {it:contreatelasticitytype(eyex)}/{it:(eydx)} are now computed as {it:dydx}/{it:dyex} on the linear-index scale (see {opt semi:elasticity}), and no setting of this option recovers the {it:eydx}/{it:eyex} number an earlier version reported, at the cost of the run time the collapse would have saved
(see {opt fastcollapse} for the magnitudes on both sides). Use it when you need byte-for-byte reproduction of an archived result, when you are comparing output across runs at
the last digit, or when memory is tight and the extra copy of the data held during the collapsed window is unwelcome. Forwarded to {cmd:permute()} and {cmd:bootstrap()}
replications. Cannot be combined with {it:fastcollapse} (the pair is refused as contradictory). {p_end}

{phang}{opt permute(method)}: Randomization (permutation) inference over the treatment assignment, for settings - most importantly few treated clusters - where cluster-asymptotic standard errors are unreliable.

{phang2}{ul:What the test assumes.} The reference distribution is generated by reassigning treatment onset across clusters (within {cmd:permstrata()} strata), so the p-value tests the sharp null under the assumption that onset was assigned {ul:exchangeably} across clusters: any cluster could have received any onset. It is not robust to the treated clusters being atypical in size or in error variance - relabelling then mixes clusters of different variance, and the reference distribution is narrower than the sampling distribution of the observed statistic (MacKinnon and Webb 2020). Measured on a 12-cluster null design whose 3 treated clusters have three times the control error standard deviation, the coefficient statistic rejected at the 5% level in 6 of 30 datasets. The default statistic is the raw coefficient (fast: the replications run without variance estimation); {cmd:permtstat} randomizes the studentized statistic, MacKinnon and Webb's recommendation when clusters may be heterogeneous, but it protects only when the cluster-robust standard error is informative about the treated clusters - many observations per treated cohort-period cell, with {cmd:cluster()} coarser than the cell. With {cmd:cluster()} at {it:i} and one observation per {it:i}-period the treated cells are absorbed by their own cohort-period terms, the standard error cannot see their variance, and neither statistic is protected (the same design rejected in 7 of 30 datasets under {cmd:permtstat}); a caution says so whenever {cmd:cluster()} is at {it:i}, whichever statistic is being randomized: under {cmd:permtstat} it says the studentization cannot see the treated clusters' own variance there, and on the default coefficient path it says the coefficient statistic carries no protection against a treated cluster whose errors are unusually dispersed, pointing at {cmd:permtstat} with that same proviso. The remedy the literature recommends for that case - a studentization by a standard error that sees the treated clusters, or the wild cluster bootstrap-t - is not offered here.
{it:wooldid} re-randomizes which clusters carry which treatment-onset dates, re-estimates every requested quantity under each counterfactual assignment, and ranks the observed estimates among the
resulting null distributions. {ul:Every p-value and interval reported under permute() is permutation-based}: the main, pre-treatment and event-study rows, the subgroup tables and their pairwise contrasts,
the joint pre-trend tests, the across-subgroup equality tests, the adversarial line tests and the continuous-treatment plot cells all carry permutation p-values and permutation intervals, and the event-study
and continuous-treatment plots - drawn after the replications from the rewired tables - carry the same intervals. No analytic p-value is ever shown under {cmd:permute()}; a quantity with no usable permutation
distribution is left blank rather than filled from another source; only the {it:histogramcohorteffects} output is computed from the stage-1 analytic covariance before the replications run and remains analytic, and a caution says so. Requires {cmd:cluster()} with a single variable - the randomization
unit (the randomization design consists of the clusters having at least one row where the {it:if} expression evaluates to exactly 1 with a nonmissing cluster value; rows {it:wooldid}
later removes by its deterministic markout, e.g. missing outcomes or zero weights, still count toward the design and its possible-assignment count) - and every cluster must sit inside exactly one treatment-onset group (onset constant within cluster; never-treated is its own group; the cluster must also nest {it:i}, i.e. cluster at or above the unit level, since each replication assigns one onset per cluster, unless {cmd:permonsetcohorts} asserts that {it:i} is exactly the treatment onset grouping, in which case the replication rebuilds {it:i} from the onsets it drew and a {cmd:cluster()} finer than {it:i} is admissible; never-treated onsets must be coded as plain missing {cmd:.} - tagged missing codes {cmd:.a}-{cmd:.z} are refused, recode them to {cmd:.} first). Two methods: {cmd:permute(shuffle)} permutes the observed onset labels
across clusters, holding the number of clusters at each onset fixed ({cmd:permstrata(varname)} confines the shuffling to within strata, so each onset group keeps its composition with respect to the strata; the stratum variable must be numeric, a string one being rejected by Stata's option parser before {it:wooldid} can say anything about it, and it is accepted only with {cmd:permute(shuffle)}; with {cmd:permute(probabilities)} it is an error, because there the strata design is expressed through the per-cluster rows of {cmd:permprobs()} instead;
every cluster must sit inside a single stratum, with no missing stratum values); {cmd:permute(probabilities)} draws each cluster's onset from user-supplied probabilities via {cmd:permprobs()}. Specify {cmd:permreps(#)} (draws with replacement) or,
for shuffle only, {cmd:permexact} (enumerate every distinct assignment; a caution prints above 1,000 assignments and enumeration is refused above 10,000). The number of possible assignments is always computed and
returned, with a caution when {cmd:permreps()} exceeds 10 percent of it; a design admitting only one possible assignment (e.g. every stratum onset-homogeneous, or every {cmd:permprobs()} row deterministic) is refused - there is nothing to randomize.
{ul:What is reported.} The permuted statistic is the coefficient itself, or its t-ratio under {cmd:permtstat}. p-values are two-sided {ul:equal-tailed} permutation p-values: with {it:P_lo} and {it:P_hi} the shares of
replications at or below and at or above the observed statistic, {it:p} = min(1, 2 min({it:P_lo}, {it:P_hi})); sampled draws use the add-one Monte-Carlo form (1 + count)/(replications + 1), and exact enumerations plain
shares (the observed assignment is a member of its own reference set, so an exact two-sided p has a floor of 2/{it:M}; an element that some enumerated assignments could not estimate keeps the exact share over the assignments that completed - its floor is 2/{it:n} on its own count {it:n},
and a caution names it. Only an element whose reference set does not contain the observed statistic at all - the observed assignment's own replication failed to carry it -
falls back to the add-one form, under its own caution). {cmd:permabsp} selects the absolute-value convention instead: {it:p} is the share of replications whose absolute statistic is at or above the observed absolute
statistic, under the same add-one and exact rules. Intervals are {ul:shift intervals}: the estimate minus the equal-tailed quantiles of the permutation distribution (times the analytic standard error under {cmd:permtstat}), at
{cmd:level()}; under {cmd:permabsp} they are the symmetric duals of the absolute-value p-value, the estimate plus or minus the {cmd:level()} quantile of the permuted absolute statistics, so that in either
convention the interval excludes zero exactly when the p-value falls below the level (up to the add-one adjustment). A shift interval inverts the permutation test under the fixed-zero-null pivotality approximation - the estimate's sampling error is taken to have the null permutation distribution - and is
{ul:not} a sharp-null randomization-test inversion (which would re-permute adjusted outcomes at every candidate value). {bf:Its width grows with the true effect}: the permutation distribution of the coefficient is generated by relabelling the observed outcomes, which carry the effect, so its spread scales with the effect (measured on a 12-cluster staggered fixture with the analytic standard error held at .38: 95% shift-interval width 1.9 at a true effect of 0, 4.7 at 5 and 18.3 at 20), and the interval - like the implied standard error in {it:e(perm_se_shift)}, which is not a sampling standard error away from the null - is calibrated only near the null. {cmd:permtstat} removes that growth for the post-treatment rows (widths 2.0, 1.4 and 1.6 on the same fixture) but not for the pre-treatment rows or the pairwise contrasts, whose permuted statistics still carry the post-period effect. This is not the Conley-Taber construction, which draws its reference from placebo assignments among the never-treated units only. The standard-error column is the analytic standard
error when one is computed (labelled {it:se_analytic}; blank otherwise), the statistic column {it:t_perm} is the permuted t-ratio under {cmd:permtstat} and blank otherwise, and the interval-implied standard errors are
in {it:e(perm_se_shift)}. The pooled-estimand summaries in the footer (p, the {cmd:level()} null band - the equal-tailed quantiles of the null draws themselves, in either convention - and the SD of the permutation distribution) are the same numbers the tables carry.
{ul:Fast-mode default.} Without {cmd:permtstat} the observed run and the replications all run in fast (no-SE) mode: nothing the permutation inference reports needs the analytic standard errors, so computing them
is wasted {it:margins} work. {it:e(V)}, {it:e(b_es_*)} and {it:e(V_es_*)} are then not posted and the {it:honestdid}/{it:pretrends} handoff is unavailable. {ul:So that the standard table tools still work}, the run posts its permutation inference under the names {it:estout} looks for before it looks at any covariance: {it:e(p)} (the same numbers as {it:e(perm_p)}) and {it:e(ci_l)}/{it:e(ci_u)} (the shift interval, the same as {it:e(perm_ci)}). {cmd:esttab, p} and {cmd:esttab, ci} therefore carry this run's real permutation figures under truthful footers, and a plain {cmd:esttab} tables the coefficients. There is deliberately no {it:e(se)}: a permutation run has no standard error, and the shift standard error in {it:e(perm_se_shift)} is a different object, so naming it {it:se} would put a number in a column that has none. Nothing else changes - with no {it:e(V)} posted, {it:test}, {it:lincom}, {it:_se[]} and {it:estat vce} go on refusing, which is correct. {cmd:permanalyticse} restores the analytic standard errors
(shown for reference only; {it:e(V)} and {it:e(V_es_*)} are then analytic, and {it:test} or {it:lincom} run by hand on them are analytic quantities). The cell collapse applies to the replications as to the observed
run ({it:fullprecisionslow} opts out); its degenerate-subgroup fallback reads a variance, which under the fast-mode default only the subgroup-by-horizon margins call carries, so elsewhere it engages only under {cmd:permtstat} or {cmd:permanalyticse}.
{ul:Joint, equality and adversarial tests} are one-sided permutation tests of a quadratic-form statistic in the {ul:permutation-covariance metric} - the covariance of the replications' coefficients over the
relevant elements (replications plus the observed comparator, listwise complete; under a complete {cmd:permexact} enumeration the enumerated assignments alone, the observed one being among them) - so they exist in both modes and consume no per-replication covariance: the joint pre-trend tests ({cmd:jointtests}) test the
pre-treatment horizons jointly, the across-subgroup equality tests use the chain contrast of the subgroup levels, and the adversarial line test ({cmd:adversarial} with {cmd:esfixedbaseperiod}) is {it:wooldid}'s
line-fit statistic evaluated on that metric. Their p-values, like every pre-treatment row and contrast, test the {ul:sharp null of no effect at any horizon}: in staggered designs the reassigned onsets move actual
post-treatment outcomes into the draws' pre-treatment horizons, so the null distribution widens with the true post-treatment effect while the observed pre-period statistics do
not, and the pre-period p-values rise toward 1 as the true effect grows. A large pre-period p-value is therefore not by itself evidence of no pre-trend (the tests lose power
against pre-trends as the effect grows; with a single onset date the reference distribution is unaffected), and the adversarial p is not the composite on-a-line test of the
analytic run. The replication statistics are in {it:e(perm_jointW)},
{it:e(perm_sgW)} and {it:e(perm_advW)}. Pairwise subgroup contrasts are permuted as coefficient differences in both modes. {cmd:jointtests} under {cmd:permute()} is exempt from the analytic joint-test screen
described under {cmd:jointtests} - randomization inference is what that screen recommends.
{ul:Subgroups} enter the replications. Because any cluster can be assigned treatment in a replication, the subgroup variable must satisfy the subgroup screens (non-negative integers; constant within
{it:i}-by-{it:t} cells and within units) over the whole sample, not only on the observed treated units, and the run is refused otherwise; a subgroup level absent from some replications contributes fewer draws
({it:e(perm_n)}). Under {cmd:permute(shuffle)} without {cmd:permstrata()} the shuffle also re-randomizes which subgroups the treated clusters fall in; {cmd:permstrata({it:sgvar})} holds each subgroup's onset
composition fixed where membership is cluster-level. With {cmd:contreat()} the dose stays with its unit and only onset timing is re-randomized (the null: timing as-if random given dose); the {it:ame} rows are
permuted like the {it:att} rows, with replication cells that lack dose variation (exact structural zeros) excluded and counted. The {it:gteffects} cells are keyed by cohort onset - exactly what a replication
reassigns - so they have no permutation distribution and their inference columns are blank. The continuous-treatment plot cells are permuted for assignment-invariant lattice types ({it:all}, {it:byCohort},
{it:byTime}, {it:byCohortByTime}, or a user {cmd:cfxlattice()} variable) and blanked, with a caution, for the {it:baseTreated*}, {it:byTreated*} and {it:byTtre*} types, whose bins move with the assignment. These matrices are the one permutation surface whose column names are left exactly as {cmd:margins} wrote them. Column 1 is the observed estimate on every type. Column 2 ships under the plain name {it:se} and holds the analytic standard error, not a permutation quantity. It is not renamed {it:se_analytic} as it is on the main and subgroup tables, because these names are consumed by {it:wooldid}'s own plotting code and by {cmd:marginsplot}-shaped code outside it, so the provenance is disclosed in the footer instead. On the permuted types, columns 4 to 6 are the permutation p-value and the permutation shift interval, and column 3 is the estimate over that analytic standard error under {cmd:permtstat} and missing otherwise. On the assignment-dependent types, "blanked" means columns 4 to 6 are set missing always and column 3 missing unless {cmd:permtstat} is in force; column 2 still holds the analytic standard error. So a reader exporting these matrices after a {cmd:permute()} run must not take the {it:se} column for a permutation standard error, and must not read a blank interval as an absent estimate. The shift construction presumes the null distribution is centred at zero, and {it:wooldid} checks that presumption rather than assuming it: when the central {cmd:level()}% band of an element's null draws (the signed equal-tailed quantiles, in whichever p-value convention is in force) excludes zero, the presumption has failed and the interval it would produce would not even contain the estimate. That element's shift interval and the implied standard error in {it:e(perm_se_shift)} are left blank, a Caution names the elements at run time, and they are listed in {it:e(perm_offcentre)}; the permutation p-value itself stands. It fires on quantities that are not zero under reassigned onsets, and under a {cmd:permute(probabilities)} design far from the observed assignment. Reassigned
onsets can change which periods are estimable in a given draw: each draw is the estimate {it:wooldid} would compute under that counterfactual assignment, estimability screens included. Replications that
fail to estimate are excluded and counted ({it:e(perm_failed)}), with a caution above 10 percent and a refusal when fewer than max(20, half) succeed (exact enumerations are exempt from that floor, but any
exact-mode failure draws a caution). When fewer than 6 of the resampling clusters are ever treated a caution names the count: few ever-treated clusters means few distinct onset assignments, so sampled draws repeat and the reference set (and with it the set of p-values the run can produce) stays coarse however large {cmd:permreps()} is. That is a limit of the design rather than of the replication count, and {cmd:permexact} answers it by enumerating the assignments and reporting the exact floor instead of sampling them. It is not the same condition as {cmd:bootstrap()}'s caution at fewer than 6 treated clusters, which is about resampling variance rather than about how many assignments exist. {cmd:permreps()} must be at least 20, and sets the resolution: a sampled-draw p-value cannot fall below 2/({cmd:permreps()} + 1) under the default equal-tailed convention (1/({cmd:permreps()} + 1) under {cmd:permabsp}), so 50 replications resolve to .039, 100 to .020, 199 to .010 and 999 to .002; below 39 a 5% rejection is unreachable under the default convention and a caution says so. Exact enumeration has its own floor of 2/{it:M} (1/{it:M} under {cmd:permabsp}). A degenerate permutation distribution (every completed draw equal) draws a caution; a missing or unreadable observed comparison
statistic is refused loudly with r(459). Draws consume the session random-number stream: give {cmd:rseed(#)}, or {cmd:set seed} beforehand, for reproducibility; {cmd:dots()} controls the progress display and {cmd:permsave()} writes the
draws. Not compatible with {it:impvalues}, {it:unconditionalse}, {it:contreatwithin}, {it:contreatcontrols(untreatedZero)} (replications reassign onset while the dose path stays fixed, contradicting
untreatedZero's dose-zero-outside-treatment premise), {it:oldsyntax}, {it:suppresssgprimaryeffects} (the suppressed subgroup entries are what the replications permute), {it:pretrendpower()} (the Roth diagnostics
need an event-study covariance, and randomization inference produces none), {cmd:simultaneousbands} (a confidence-band object; the permutation engine's multiplicity tool is the {cmd:permstepdown} step-down), or {cmd:bootstrap()}. The conventions the reported p-values follow, and the floor each one carries, are mapped at {help wooldid##pvalues:Which p-value is which}; why this engine and {cmd:bootstrap()} differ where they do is at {help wooldid##enginediff:Why bootstrap() and permute() do not behave identically}. {p_end}

{phang}{opt permprobs(spec)}: The assignment-probability design for {cmd:permute(probabilities)}, structured as {cmd:(onsetvar level level ...)} followed by one {cmd:(clustervalue p p ...)} group per cluster.
The first group names the actual timing variable and lists every in-sample onset level (never-treated is the level {cmd:.}) in the column order the probabilities use; each later group gives one cluster's
probabilities over those levels (nonnegative, summing to 1; rows within the acceptance tolerance of 1e-6 are rescaled to sum exactly 1 before drawing). Every in-sample cluster value must appear exactly once, and the {cmd:cluster()} variable must be numeric, since the grammar names clusters by numeric value. {cmd:permreps()} is required and {cmd:permexact} is refused. The possible-assignment count reported is the product of the
per-cluster support sizes; if the observed assignment has zero probability under the design, a caution notes that the p-value is a tail probability under that design distribution rather than a randomization
test. {p_end}

{phang}{opt permtstat}: Randomize on t-ratios instead of raw estimates: the permuted statistic is each replication's estimate over its own analytic standard error, ranked against the observed
{it:t_perm} = estimate/{it:se_analytic}. The observed run and the replications then go through the full inference path (slower: fast mode cannot be exploited, and requesting {it:fastnose} alongside is refused);
an element whose replication or observed standard error is unusable is left blank rather than filled from the coefficient draws, and a caution counts them. The shift intervals are then the estimate minus the
equal-tailed quantiles of the permuted t-ratios times {it:se_analytic}. The joint, equality and adversarial tests and the pairwise contrasts are unchanged by {cmd:permtstat} (they are formed from the coefficient
draws on the permutation-covariance metric). The equal-tailed exact-test floor is 2/{it:M} in both modes (1/{it:M} under {cmd:permabsp} when no other assignment ties the observed absolute statistic). {cmd:permtstat} is the MacKinnon-Webb (2020) recommendation when clusters may be heterogeneous, with the proviso described under {cmd:permute()}: it protects only when the cluster-robust standard error is informative about the treated clusters, and at {cmd:cluster()} = {it:i} with one observation per {it:i}-period it is not (a caution says so). {p_end}

{phang}{opt permtci}: Retained from the earlier permutation surface: additionally report, for the pooled estimands, the symmetric inverted-|t| {cmd:level()} interval - the observed estimate plus or minus
the {cmd:level()} quantile of the permuted absolute t-ratios times the observed standard error ({it:e(perm_tlb_*)}/{it:e(perm_tub_*)}). It rests on the same pivotality approximation as the shift intervals every row
now carries and is not a sharp-null randomization-test inversion; under {cmd:permabsp} it coincides with the interval the table carries for those estimands. Implies {cmd:permtstat}. Unavailable (with a caution) when the observed run has no usable standard error. {p_end}

{phang}{opt permsa:ve(name)}: Also write the permutation draws to {it:name}{cmd:.dta} in the current working directory (a plain file name; an existing file needs {cmd:replace}). One row per replication
({cmd:_rep} = 1 to {cmd:permreps()}, or 1 to {it:M} under {cmd:permexact}; a failed replication keeps its row, all missing) plus a {cmd:_rep} = 0 row carrying the observed estimates (the rep-0 comparator the p-values rank against equals this row except under {cmd:permanalyticse}, where the comparator is a fast-mode refit; the continuous-treatment plot cells are not written - their draws live in the hidden {it:e(permdist_*)} twins), one column per
coefficient of {it:e(b)} under its {it:e(b)} name; under {cmd:permtstat} also the replication and observed analytic standard errors ({cmd:se_}{it:name}), so the permuted t-ratio is the coefficient over
its {cmd:se_}; and one {cmd:onset_}{it:clustervalue} column per cluster holding the onset that cluster carried in each replication (row 0: the observed onset). A trailing {cmd:.dta} on the argument is stripped, so {cmd:permsave(mydraws.dta)} and {cmd:permsave(mydraws)} are the same request; it is still a plain file name, with no directory part. The assignment block can come up short in two ways, each of which prints a Caution at the time: it is dropped in full when the cluster-value and observed-onset lists do not conform to the assignment matrix, leaving a file that carries the replication draws and no {cmd:onset_} columns at all, and an individual cluster's column is omitted when its {cmd:onset_}{it:clustervalue} name is already taken, either by a coefficient column or by two cluster values that map to the same Stata name. Because the columns are named rather than positional, nothing inside the saved file marks either omission, so check the {cmd:onset_} column count against the cluster count before relying on the block. {cmd:bootsave()} strips a trailing {cmd:.dta} in the same way and skips its own {cmd:w_}{it:clustervalue} columns on the same kind of name collision. Every reported permutation quantity is otherwise
reconstructable from the file. {p_end}

{phang}{opt perman:alyticse}: Under {cmd:permute()} without {cmd:permtstat}, compute the analytic standard errors anyway (the observed run then takes the full inference path) and show them in the
{it:se_analytic} column for reference; the p-values and intervals remain permutation-based, and the {it:t_perm} column stays blank because the permuted statistic is the estimate itself. This also
restores {it:e(V)}, {it:e(b_es_*)} and {it:e(V_es_*)} (analytic) and lets the cell-collapse fallback screen engage in the replications. {p_end}

{phang}{opt permab:sp}: Report absolute-value permutation p-values instead of equal-tailed ones: each p-value is the share of replications whose absolute statistic (the coefficient, or the t-ratio under
{cmd:permtstat}; the coefficient difference for a subgroup contrast) is at or above the observed absolute statistic, with the same add-one form for sampled draws and plain shares under {cmd:permexact}.
The intervals then become the symmetric duals of that p-value - the estimate plus or minus the {cmd:level()} quantile of the permuted absolute statistics (times the analytic standard error under
{cmd:permtstat}) - so the table's p-value and interval keep agreeing. The joint, equality and adversarial tests are one-sided quadratic-form tests and are unchanged. The null band, the SD and
the off-centre check are convention-free. The convention in force is recorded in {it:e(perm_pconv)} ({it:equaltailed} or {it:absolute}; {it:raw}, {it:left} or {it:right} under the two options below). The three are alternative conventions: at most one of {cmd:permabsp}, {cmd:permside()} and {cmd:permrawp} may be specified, and the combination is refused. {p_end}

{phang}{opt permsi:de(left|right)}: One-sided permutation p-values: {cmd:permside(left)} reports the share of replications whose statistic is at or below the observed one and {cmd:permside(right)} the share at or above it, under the same add-one / exact-share rule as the default; the two shares are always stored as {it:e(perm_p_lo)} and {it:e(perm_p_hi)}, so the default two-sided p is 2 min of them. The intervals stay the two-sided equal-tailed shift intervals at {cmd:level()} (a one-sided test has no two-sided dual). {p_end}

{phang}{opt permraw:p}: The plain absolute-value share {it:c}/{it:n} with no add-one adjustment: {it:c} is the number of replications whose absolute statistic is at or above the observed absolute statistic and {it:n} the completed replications. This is the number {helpb ritest} and Stata's own {helpb permute} prefix report, offered so a referee's cross-check can be reproduced; its intervals are the symmetric duals used under {cmd:permabsp}. {p_end}

{phang2}{ul:Which convention, and the floors.} wooldid's default is Rosenbaum's equal-tailed two-sided p-value, 2 min({it:P_lo}, {it:P_hi}), with the add-one Monte-Carlo form (1 + count)/({it:n} + 1) for sampled draws (Davison and Hinkley 1997; Phipson and Smyth 2010) - the convention {helpb boottest}'s ptype(equaltail) uses - because it is the exact dual of the equal-tailed shift interval and spends alpha/2 in each tail, which matters for asymmetric nulls. {helpb ritest}, Stata's {helpb permute} and Young (2019) instead report the absolute-value share {it:c}/{it:n} without the add-one adjustment ({cmd:permrawp} here; {cmd:permabsp} is the same statistic with the add-one form), and MacKinnon and Webb (2020) use the absolute-value convention too. The conventions differ in their floors: with {it:n} sampled draws the smallest attainable p is 2/({it:n} + 1) under the default, 1/({it:n} + 1) under {cmd:permabsp} and 0 under {cmd:permrawp}; under {cmd:permexact} with {it:M} assignments it is 2/{it:M} under the default and 1/{it:M} under the absolute-value conventions (2/{it:n} and 1/{it:n} for an element that some assignments could not estimate, with {it:n} its count in {it:e(perm_n)}) - with one treated cluster among 12 the exact default p cannot fall below .167 and the absolute-value p below .083. Every convention also comes with its Monte-Carlo precision: {it:e(perm_p_se)} is the standard error of the reported p-value's underlying share, sqrt({it:p}(1 - {it:p})/{it:n}) (doubled under the equal-tailed rule), and {it:e(perm_p_ci)} its {cmd:level()} Clopper-Pearson interval, as {helpb permute} and {helpb ritest} print. Under {it:permexact} the shares are exact and carry no Monte-Carlo error, so both are left missing; the exception is an element that fell back to the add-one form, whose share is over a set that excludes the observed statistic and does carry both. {p_end}

{phang2}{ul:Multiplicity across the horizons.} {cmd:permstepdown} reports, for each family and subgroup with two or more horizons, the max-|studentized| step-down of Westfall and Young / Romano and Wolf from the same replications: each horizon's draws are studentized by the standard deviation of its null draws (the fixed metric the joint tests use), the horizons are ordered by the observed |t|, and the j-th adjusted p-value is the share of replications whose largest |t| over the horizons not yet rejected reaches the observed one, made monotone. It is exact under the sharp null of no effect at any horizon of the family and controls the family-wise error rate across those horizons, which the per-horizon p-values in the table do not. The post-treatment and pre-treatment horizons are run as {ul:separate} families, adjusted among themselves and reported and stored apart: pooling them would change each post-treatment adjusted p-value by enlarging its family, and the two answer different questions. The pre-treatment family inherits the pre-treatment caveat below - in staggered designs the reassigned onsets widen the null for pre-period statistics - so its adjusted p-values are conservative in a way the post-treatment family's are not. The procedure is permutation-only; under {cmd:bootstrap()} the corresponding object is a simultaneous band ({cmd:simultaneousbands}), which is a uniform confidence band rather than a family-wise-error test. The stored kernel's {it:p_raw} row (and the footer's first column) is the two-sided |t| share, add-one, in the step-down's own studentized metric - not the table's convention - so compare it with the step-down p rather than with the table. The results are stored in {it:e(perm_stepdown)} ({it:_i} for {it:att_it}, {it:perm_Cstepdown} for the continuous-treatment family, {it:_pre} for the pre-treatment family, {it:_s#} per subgroup; rows {it:t_obs}, {it:p_raw}, {it:p_singlestep}, {it:p_stepdown}, one named column per horizon) and summarised in the footer. The step-down runs on the replications that are listwise complete over {ul:every} horizon of the family rather than horizon by horizon, so its reference set can be considerably smaller than the per-horizon counts in {it:e(perm_n)}, which is exactly what happens on a design where a reassigned onset makes some horizon inestimable in a draw. A family's step-down can also fail to be formed at all, when fewer than two replications are complete over the family, when the comparator row is missing a horizon, or when some horizon's null draws have no dispersion; the matrix is then simply not posted and nothing is printed, so test for {it:e(perm_stepdown*)} by name rather than reading the footer's silence as evidence that multiplicity is not a problem. A family whose stored name would run past Stata's 32-character limit is skipped the same way. The footer summarizes the four overall families and their pre-treatment twins; the per-subgroup matrices are returned in {it:e()} and are not summarized there. Before 2026-09-07 this procedure ran on every {cmd:permute()} run with two or more post-treatment horizons; it is now opt-in, and {cmd:permstepdown} reproduces exactly what it used to print for the post-treatment family. {p_end}

{phang}{opt permstepdown}: Compute and report the max-|studentized| step-down described above. Requires {cmd:permute()}; off by default. {p_end}

{phang}{opt permonset:cohorts}: Assert that {it:i} {ul:is} the treatment-onset grouping, with exactly one level of {it:i} for each distinct onset and the never-treated forming a single level of their own, and on that assertion accept a {cmd:cluster()} finer than {it:i}. The two choices answer different questions. {it:i} says which cohorts the estimands average over; {cmd:cluster()} says which observations are dependent and what the design randomized. An {it:i} set at the treatment-onset level can therefore sit above the randomization unit, which is the pairing {help wooldid##standarderrors:What level to cluster at} recommends and the one this engine refused. The refusal was mechanical rather than methodological: each replication assigns one onset per cluster and writes it to every row of that cluster, so a unit straddling two clusters would carry two onsets at once, which the estimator has no reading for. Under this option each replication instead rebuilds {it:i} from the onsets it drew, since the units sharing a drawn onset {ul:are} that draw's cohort, so the onset is single-valued within the rebuilt {it:i} by construction. {p_end}

{phang2}{ul:Why it has to be declared.} {it:wooldid} cannot tell from the data whether {it:i} is an onset grouping or a substantive unit definition that happens to coincide with one on this sample, and rebuilding the wrong one would silently redefine the estimand. The declaration is then checked, on the estimation sample and with the never-treated included: an {it:i} carrying two onsets, or an onset spread across two levels of {it:i}, stops the run with r(459), and {it:safety(off)} does not lift it, since without the one-to-one correspondence the rebuild is not a defined map rather than an unwise one. Where {cmd:cluster()} already sits at {it:i} the rebuild has nothing to do and is skipped, so the option is then a declaration that moves no number, and {it:e(perm_irebuilt)} records which of the two happened. {p_end}

{phang2}{ul:What it composes with.} {cmd:permstrata()} composes, since the strata constrain which onsets a cluster may draw and the rebuild sits downstream of that; a rebuilt cohort can then pool clusters from different strata, which is the cohort definition doing its job, a cohort being an onset and not a stratum. {cmd:permute(probabilities)} composes as well, and there an onset can go undrawn, so a replication's rebuilt {it:i} can hold fewer cohorts than the observed sample and the cohort-averaged estimands are then averages over the cohorts that draw produced, which is the ordinary reading of a functional of the replication's own data. {cmd:subgroups()} is refused: a rebuilt cohort-by-period cell pools the clusters that drew one onset, those clusters come from different observed cohorts, and the subgroup variable then varies inside the cell, which the estimator refuses in every replication. A finer {cmd:cluster()} is also what makes {cmd:permtstat} protective, since the studentized statistic guards against a treated cluster atypical in variance only where the cluster-robust standard error can see that cluster's own variance, which at {cmd:cluster()} equal to {it:i} it cannot. {p_end}

{phang}{opt simult:aneousbands}: Compute a {ul:simultaneous} (sup-t, or uniform) confidence band over each event study's horizons, at {cmd:level()}. Off by default; no shipped figure or table changes unless it is asked for. The pointwise interval printed on each table row covers that
horizon with probability {cmd:level()}; a simultaneous band covers {ul:every} horizon of the family at once with that probability, which is what a reader scanning a whole event study for "no effect anywhere before onset" actually needs. It is the estimate plus or minus one common critical value times each horizon's own standard error, not the printed pointwise interval scaled up, and it is symmetric about the
estimate. So where the printed pointwise interval is a rank interval and not symmetric - under {cmd:bootstrap()} with the default {cmd:bootdisplay(percentile)}, or with {cmd:bootdisplay(bc)} or {cmd:bootdisplay(bca)} - the band does not nest it endpoint for endpoint: one of the band's two ends can sit inside the printed interval, even where the band is the wider of the two. The band is a property of the family rather than of a row, so the printed table keeps its pointwise intervals: the endpoints are stored in {it:e(esband_{it:type})} (rows {it:lb}, {it:ub}, {it:relyr}; one column per horizon, per subgroup as
{it:e(esband_{it:type}_s{it:sgcase})}) with the critical value in {it:e(esband_crit_{it:type})}, and the event-study plot draws them as a dashed pair over the pointwise shading. Two engines, chosen by the run: with {cmd:bootstrap()} the band is the {cmd:level()} quantile, over
replications, of the largest studentized deviation from the observed estimates, computed on the listwise-complete replications in one pass over the accumulator that already feeds the joint test ({it:e(esband_n_{it:type})} records how many). Under {cmd:boottstat} that deviation is studentized by each replication's {ul:own} standard error and the endpoints use the observed analytic ones - a genuine percentile-{it:t}, on the same scale as the pointwise intervals printed beside it, so the band can never come out narrower than the interval it accompanies; {it:e(esband_engine)} reads {it:bootstrap-studentized} there. Without {cmd:boottstat} the denominator is the standard deviation of the coefficient draws, which is the only quantity available. Without {cmd:bootstrap()} the band is simulated from the
family's own analytic covariance {it:e(V_es_*)} - 10,000 draws from the implied correlation matrix, on a fixed internal seed with the caller's random-number state saved and restored, so the band is reproducible and the user's stream is undisturbed. {it:e(esband_engine)} records which.
Not available with {cmd:permute()}, whose multiplicity object is the {cmd:permstepdown} step-down instead, and unavailable under {cmd:fastnose} or where a family has fewer than two horizons with a usable variance. The two engines treat an unusable horizon differently. The analytic engine drops a horizon whose variance is missing or not positive, whose coefficient is missing, or whose name does not parse as a relative period, and bands the rest: {it:e(esband_{it:type})} can therefore carry fewer columns than {it:e(b_es_{it:type})}, the critical value is a maximum over the retained horizons only, and the band is simultaneous over those horizons rather than over the family as the table prints it. The bootstrap engine refuses the whole family for a single unusable horizon (a missing observed coefficient, a missing or non-positive observed analytic standard error under {cmd:boottstat}, or a draw standard deviation of zero), because its critical value is a maximum taken replication by replication, so a horizon that cannot supply a value at every retained replication changes what that maximum is over; it also needs at least two replications that are complete over every horizon of the family. A family with five good horizons and one bad one therefore gets a band from the analytic engine and none from the bootstrap engine. A Caution names every family that was attempted and declined, with its reason, both when no family at all could be banded and when some were banded and others were not, and under the bootstrap engine the reason names the horizon responsible; the figures for those families show their pointwise intervals only. Only families that were attempted are named: one whose inputs are absent from the run altogether ({it:e(V_es_ame*)} withheld by the cell-containment screen, say) is passed over before there is a reason to give. {it:e(esband_families)} lists the families that were banded. {p_end}

{phang}{opt boot:strap(mode)}: Bootstrap standard errors and inference. When specified, the reported standard errors, z statistics, p-values, confidence intervals, joint tests and subgroup contrasts are bootstrap quantities, replacing the analytic ones; point estimates are unchanged, and the analytic covariance is preserved in the hidden {it:e(V_analytic)} (and {it:e(V_es_analytic_*)}). Each replication re-estimates the full model, so runtime scales with {cmd:reps()}. Two modes. {cmd:bootstrap(standard)} is the cluster bootstrap: clusters - the {cmd:cluster()} variable, or observations when no {cmd:cluster()} is given - are resampled with replacement. For the observation-weighted estimands ({it:att}, {it:ame}) replications are computed by cluster frequency reweighting, which is the same statistic as physical resampling and matches Stata's own {helpb bootstrap} prefix replication for replication under a common seed. When a cell-weighted estimand ({it:att_it}, {it:att_i}, {it:att_itime}) is requested, replications instead physically resample, and what becomes of the unit identities on that frame depends on where {cmd:cluster()} sits relative to {it:i}. Where the cluster nests {it:i}, each drawn copy of a cluster becomes a fresh set of units, which is what makes a twice-drawn cluster contribute its treatment cells twice. Where the cluster is finer than {it:i}, or cuts across it, the units are kept exactly as they are and the drawn copies arrive as duplicate rows inside their own cells: re-identifying there would split a unit straddling two clusters into two units in {ul:every} replication whether or not anything had been duplicated, turning the unit fixed effects into effects per cluster per unit, so the rule that is right under nesting is silently wrong under a finer cluster and the leg carries both. With no {cmd:cluster()} the observations themselves are resampled and the units are again kept whole. None of the three is a compromise for the estimand: {it:att_i} is a functional of whatever data a replication holds, and each replication rebuilds its weights from that data. The same physical leg is taken by {cmd:boottstat} whenever {cmd:cluster()} is given, for a different reason: a studentized bootstrap needs each replication's own standard error to be the same functional of the resampled data as the observed standard error is of the original, and only a physical expansion lets the replication cluster on the {ul:resampled} cluster identity - under frequency reweighting a twice-drawn cluster's copies were pooled into one score and the replication standard errors were inflated by roughly the square root of two (their mean square relative to the observed one is 2 - 1/{it:G} in expectation). So {cmd:boottstat} with {cmd:cluster()} takes that same physical leg under either regime, since what is given a fresh identity is always the drawn {ul:cluster}, which is what the replication's own standard error has to cluster on, and that is independent of whether the unit identities move with it. It costs the re-identification overhead on top of the replications running with full variance estimation; on a 12-cluster fixture the corrected replication standard errors sit at about 0.92 (plus or minus 0.04) of the observed one in mean square. Without {cmd:cluster()}, {cmd:boottstat} replications carry the resample multiplicity as a frequency weight, which is again the same functional as the observed standard error; when a user weight ({cmd:aweight} or {cmd:pweight}) is also present that frequency form is unavailable, so the replications physically duplicate the drawn observations, carrying the user's weight onto every copy and keeping the original {it:i} (no re-identification is needed, because the resampling unit is the observation). Re-identification roughly doubles the cost of each replication, because the resampled data carries distinct cohort identities for duplicated clusters and so presents a larger cell count to the regression. What is re-identified is the cell and aggregation identity only. A term inside {it:fe()}, {it:controls()}, {it:interactivecontrols()}, {it:coarsecohortcontrols()} or {it:ccc_absorb()} that names the {it:i} variable itself keeps its original codes across the copies, so a twice-drawn cluster is two cohorts to the aggregation and one unit to those terms, the same convention Stata's own {helpb bootstrap} prefix follows under {cmd:idcluster()}. {it:wooldid} cautions when it sees the bare {it:i} variable named in one of those options; it is a caution rather than a refusal, because the built-in prefix behaves the same way. The check compares bare tokens, so {it:fe(state)} is noticed where {it:fe(i.state)} is not. That cost is inherent to making those estimands' bootstrap correct rather than an inefficiency, but it is worth budgeting for when choosing {cmd:reps()}. {cmd:bootstrap(bayes)} is the Bayesian bootstrap: rather than resampling, each replication draws one Dirichlet(1,...,1) weight per cluster and multiplies every observation of that cluster by it (on top of any existing weights), so a cluster's total mass stays proportional to its size in every replication - exactly as under the pairs cluster bootstrap, whose first two moments it shares - no cluster ever leaves the sample, and replication estimability is stable; under {cmd:bootalpha(size)} the Dirichlet shape is instead the cluster's observation count and the drawn mass is spread equally over the cluster's rows, so on unweighted data the cluster's replicate total is the draw itself. Because reweighting cannot resample cells, {cmd:bootstrap(bayes)} combined with a cell-weighted estimand omits the cell-sampling variance component and can understate those standard errors when treatment effects vary across clusters: it is refused under default safety and proceeds with a caution under {cmd:safety(off)}, while {cmd:bootstrap(standard)} carries no such limitation. {bf:Neither mode is a remedy for few treated clusters.} What is implemented is the pairs cluster bootstrap (and its Bayesian counterpart, which shares its first two moments): replications that draw no treated cluster fail and are excluded, so with few treated clusters the retained distribution is conditional on drawing at least one, and the few-cluster literature (MacKinnon and Webb 2017, 2018; Roodman, MacKinnon, Nielsen and Webb 2019) finds the pairs cluster bootstrap badly sized there and recommends the wild cluster bootstrap-t instead. wooldid does not offer a wild cluster bootstrap: a self-contained one of the margins functional would be a new engine with a full refit per replication, and the route through {helpb boottest} on the stored model was assessed and ruled infeasible (2026-09), so the tool for that setting is {helpb boottest} run on the stage-1 regression outside wooldid, or {cmd:permute()}. A caution fires when fewer than 6 treated clusters are available for resampling (beside the existing caution at fewer than 10 clusters in total), and {cmd:boottstat} is not a cure for it. Below two resampling units it is a refusal rather than a caution: {cmd:bootstrap()} requires at least two clusters in the estimation sample, or two observations on the unclustered leg, since resampling or reweighting a single unit produces no variation at all. That one fires when the options are screened, before the estimation runs. In a 24-cluster staggered simulation (review_process/bootstrap_report_2026-08-30.md) the Wald-family interval covered the true {it:att} 89.4% of the time at a 90% nominal level under {cmd:bootstrap(standard)} (standard error 0.967 of the truth) and 85.1% under {cmd:bootstrap(bayes)} (0.877; a balanced design, on which the per-cluster multiplier and the pre-2026-09-06 total-mass split coincide up to scale), against 89.8% (0.869) for the analytic clustered interval; the percentile and bias-corrected families were not measured there. Requires {cmd:reps()}. Incompatible with {cmd:permute()} (one engine at a time), {cmd:impvalues} and {cmd:unconditionalse} (each replaces the inference this option produces), {cmd:fastnose} (fast mode reports point estimates only, so there is no inference surface left for this option to replace, and the engine already runs every replication in fast mode itself whenever {cmd:boottstat} is not requested, {cmd:boottstat} being what turns it off, since each replication must then estimate its own variance), {cmd:customvce()} (the engine has to know what a replication's variance is; an arbitrary vce string is not inspectable), {cmd:suppresssgprimaryeffects} (the suppressed subgroup entries are what the replications resample), {cmd:oldsyntax}, and multi-way {cmd:cluster()} (the engine resamples one unit, so there must be exactly one). The same requirements in the same shape, beside {cmd:permute()}'s, are in {help wooldid##enginediff:the engine-compatibility table}. {cmd:pretrendpower()} {bf:is now supported}: the Roth diagnostics are computed after the bootstrap rewiring, so Roth's own formulas are fed the bootstrap event-study covariance rather than the analytic one, and the reported power, detectable trend and slope figures are bootstrap-based. The method is unchanged - only the covariance it consumes - and the existing rank guard applies to the bootstrap block as it does to the analytic one, skipping with a caution if it is not invertible. With {cmd:jointtests} the joint pre-trend tests (and, with {cmd:subgroups()}, the across-subgroup equality tests) are bootstrap quantities in both forms: under {cmd:boottstat} each is the bootstrap of its test statistic (each replication's own covariance), and without {cmd:boottstat} it is the fixed-metric rank test the permutation engine also uses - the replication draws are recentred at the observed estimates, W is the quadratic form of the tested contrast on the listwise bootstrap covariance of those columns, and the p-value is the add-one share of replications whose W is at or beyond the observed one; {it:e(boot_jt_family)} records which ({it:studentized} or {it:fixedmetric}), the replication statistics are in {it:e(boot_jointW)} and {it:e(boot_sgW)} with the observed statistics in {it:e(boot_jointWobs)} and {it:e(boot_sgWobs)}, and a test with no usable replication statistic is absent with a caution rather than analytic. Neither is the Wald test on {it:e(V)} that a hand {cmd:test} reproduces; that one is a normal-theory test on the bootstrap covariance and is not printed. That combination is exempt from the analytic joint-test screen described under {cmd:jointtests} - bootstrap inference is what that screen recommends. The event-study and continuous-treatment plots are drawn after the rewiring from the bootstrapped tables and matrices, so they and their saved plot data carry the displayed bootstrap family's intervals, and the {it:gteffects} cells and {it:e(cfx*)} matrices are rewired too, except under {it:cfxplotnose}, which truncates the {it:e(cfx*)} matrices to their point-estimate column before the replications run, so there is nothing for the rewiring to fill and those cells come back with no bootstrap standard error, p-value or interval and their plot rows are drawn as points, without a message. {cmd:permute()} does not share that behaviour: it carves {it:cfxplotnose} out of the truncation precisely so its cells can be rewired. If bands are wanted on the continuous-treatment figures under {cmd:bootstrap()}, drop {it:cfxplotnose}. Only the {it:histogramcohorteffects} output is computed from the analytic covariance before the replications run and remains analytic, and a caution says so. Stata's own {helpb bootstrap} prefix can also be used on an unweighted wooldid call ({cmd:bootstrap _b, reps(#) seed(#) cluster(c) idcluster(nc): wooldid y i t ttre, ... fastnose}), but it rejects weighted commands, cannot deliver the Bayesian flavor, and discards replications whose coefficient vector lost a horizon - the built-in {cmd:bootstrap()} option avoids all three. The conventions the reported p-values follow, and the floor each one carries, are mapped at {help wooldid##pvalues:Which p-value is which}; why this engine and {cmd:permute()} differ where they do is at {help wooldid##enginediff:Why bootstrap() and permute() do not behave identically}. {p_end}

{phang}{opt reps(#)}: Number of bootstrap replications; required with {cmd:bootstrap()}, and refused with an error on a run that does not use it. No minimum is imposed, so small values can be used to time a specification - but bootstrap inference needs many replications: hundreds for standard errors, and a thousand or more before the percentile and bias-corrected intervals are worth reading. Replications that fail to estimate - for example a resample retaining no comparison group - are excluded and counted in {it:e(boot_failed)}, the convention Stata's own {helpb bootstrap} follows.
{ul:Resolution.} A rank-based p-value cannot be finer than its own grid, and {cmd:wooldid} prints a Caution naming that grid whenever it is coarser than .05 - the counterpart of the one {cmd:permreps()} prints, but without {cmd:permreps()}'s hard floor of 20. Which number binds depends on the family the table
carries: the {cmd:bootdisplay(percentile|bc|bca)} coefficient p-values are plain two-sided shares on a grid of 2/{it:completed}, whose only value below .05 is exactly 0 (read as "below 2/{it:completed}", not as strong evidence); {cmd:boottstat}'s studentized p-values are add-one ranks with both grid and floor
2/({it:completed} + 1), so they can never print 0; and under {cmd:bootdisplay(wald)} without {cmd:boottstat} the coefficient p-values are analytic and unlimited, leaving only the joint, equality and adversarial family - one-sided add-one shares with floor 1/({it:completed} + 1). About 40 completed replications
are needed before a 5% rejection is reachable, 100 for .02, 200 for .01 and 1000 for .002. {p_end}

{phang}{opt level(#)}: Confidence level for every reported confidence interval - the analytic tables, the event-study and continuous-treatment plots, the {it:gteffects} cells and the Roth slope interval, as well as the {cmd:bootstrap()} and {cmd:permute()} intervals - between 10 and 99.99; the default is {cmd:c(level)}, normally 95. The interval columns are labelled {it:lb_#ci} / {it:ub_#ci} for the level in force (a non-integer level writes its decimal point as {it:p}: {it:lb_92p5ci}), and that level is stored in {it:e(level)}. The analytic intervals honour {cmd:level()} since 2026-09-06; before that it was accepted only with {cmd:bootstrap()} or {cmd:permute()}. The pre-test sizes behind {cmd:pretrendpower()}'s {it:mdt_wald} and {it:mdt_nis} columns are fixed at 5% and do not follow it; within that entry {cmd:level()} reaches the slope interval only. {p_end}

{phang}{opt rseed(#)}: Seed for the bootstrap or permutation draws, passed to {cmd:set seed} before the first replication; as with the built-in prefixes this advances the session's random-number state. Without it the draws consume the current state and the run is not reproducible. Applies to {cmd:bootstrap()} and {cmd:permute()}, is refused with an error on a run with neither, and is recorded in {it:e(boot_seed)} or {it:e(perm_seed)}. Whether or not a seed is given, the random-number state in force before the first draw is stored in {it:e(boot_rngstate)} or {it:e(perm_rngstate)} (as Stata's own {helpb bootstrap} stores {it:e(rngstate)}), so a run without {cmd:rseed()} can still be replayed exactly by {cmd:set rngstate} to that value before repeating the call; {cmd:bootsave()} and {cmd:permsave()} files carry the engine, the seed and that state as dataset characteristics {it:_dta[wooldid_engine]}, {it:_dta[wooldid_seed]} and {it:_dta[wooldid_rngstate]}. {p_end}

{phang}{opt dots(#)}: {bf:Requires} {cmd:bootstrap()} or {cmd:permute()}; on a run with neither it is refused with an error rather than quietly ignored, so {cmd:dots(0)} on a plain analytic call is an error and not a way of asking for silence. Show a progress dot every {it:#} replications, with {cmd:dots(0)} suppressing the display; the default is 1. After the fifth replication the engines print, on the dots line, a projected runtime for the whole loop from the time the first five took (Stata timer 100 is used for those five replications). {p_end}

{phang}{opt ties}: {bf:Requires} {cmd:bootstrap()} - it is refused with an error on its own. Adjust the bias-corrected confidence intervals' median-bias constant for ties between the replications and the observed statistic, exactly as the built-in {helpb bootstrap}'s {cmd:ties} option does. It has no effect on continuous statistics, which do not tie, and under {cmd:boottstat} it touches only the stored coefficient-draw bias-corrected tables, never anything displayed. {p_end}

{phang}{opt bootalp:ha(spec)}: {bf:Requires} {cmd:bootstrap(bayes)}; specifying it with {cmd:bootstrap(standard)} is an error rather than a setting that is quietly ignored. The Dirichlet convention for {cmd:bootstrap(bayes)}. {cmd:bootalpha(equal)}, the default, draws Dirichlet(1,...,1) across clusters - the Rubin convention - and applies each cluster's draw as a multiplier on its observations, so a cluster's total mass stays proportional to its size and the replication distribution centres where the pairs cluster bootstrap's does when cluster sizes differ (before 2026-09-06 the draw was split as a total, which centred unbalanced-cluster replications on the equal-cluster-weighted estimand instead). {cmd:bootalpha(size)} instead gives cluster {it:c} the Gamma shape {it:n_c}, its observation count, and divides the drawn mass by that same count so that every row of the cluster carries {it:p_c}/{it:n_c}, reproducing what {helpb bayesboot} does when passed {cmd:id(}{it:clustervar}{cmd:)}. That parity is an unweighted statement, because {helpb bayesboot} refuses user weights outright; the divisor is the Gamma shape rather than the cluster's weight total, which is what keeps a weighted run's replications centred on the weighted estimand the point estimate targets. It is available for exact comparison against that command, but it is not a sound default for panel data, because tying the Dirichlet's concentration to the observation count rather than the cluster count shrinks the cluster-weight dispersion by roughly the square root of the mean cluster size and the resulting standard errors are badly understated. {p_end}

{phang}{opt bootts:tat}: {bf:Requires} {cmd:bootstrap()}, and is refused with an error without it. A studentized (percentile-t) bootstrap: the object bootstrapped is the t-statistic, {it:t(r)} = ({it:b(r)} - {it:b}) / {it:se(r)}, with each replication's own analytic standard error in the denominator, rather than the coefficient. The reported standard error is the {ul:analytic} one, labelled {it:se_analytic} in the tables - it is the studentization's own denominator, and a studentized distribution has no standard deviation in coefficient units. The reported statistic, labelled {it:t_pt}, is {it:b} / {it:se_analytic}. The p-value and confidence interval come from the ranks of the studentized draws: by default the equal-tailed percentile-t interval [{it:b} - {it:q_hi} {it:se}, {it:b} - {it:q_lo} {it:se}] from the signed quantiles of {it:t(r)}, which is not symmetric about {it:b}, and the two-sided rank of the observed statistic among the draws, 2 (1 + min(#{it:t(r)} at or below {it:t}, #{it:t(r)} at or above {it:t})) / (1 + {it:completed}), under the add-one Monte-Carlo-null convention: the studentized draws are recentred at the observed estimate, so they form a reference set drawn under the null and the observed statistic is a member of it. That is the same convention every bootstrapped joint statistic in the package uses. The interval endpoints stay plain order statistics, so the p and the interval are near-duals rather than exact ones - at the boundary an interval can exclude the null while the p sits just above the level. {cmd:bootdisplay(wald)} is refused with {cmd:boottstat}, since there is no bootstrap standard deviation of a studentized statistic for a normal or t approximation to rest on, and so are {cmd:bootdisplay(bc)} and {cmd:bootdisplay(bca)}: the percentile-t interval is already second-order accurate, and applying the bias-corrected shift to the studentized quantiles and then reflecting them moves the interval against the median bias, so there is no bias-corrected or accelerated studentized family (a bias-corrected one was offered before 2026-09-06 and has been withdrawn). Studentization applies to every coefficient in {it:e(b)} - the main, pretreatment and event-study rows and the subgroup tables - and to the {it:gteffects} cells and the {it:e(cfx*)} matrices, each by its own replication standard error. The pairwise subgroup contrasts are studentized the same way from each replication's own pairwise standard error, and the joint pre-trend tests, the across-subgroup equality tests and the adversarial line tests are bootstrapped as test statistics: each replication's statistic is recentred at the observed estimates and evaluated on that replication's own covariance (the joint tests on its {it:e(V_es_*)}, the equality tests on a covariance rebuilt from its level and pairwise standard errors), and the p-value is the add-one share of replications at or beyond the observed statistic. Without {cmd:boottstat} the joint pre-trend and across-subgroup equality tests are the fixed-metric rank tests described under {cmd:bootstrap()}; the adversarial line tests are not bootstrapped at all, since bootstrapping them means bootstrapping the statistic itself. They are posted as missing with {it:e(boot_adv_absent)} = 1 and a caution, never as analytic p-values under a bootstrap header. {cmd:boottstat} is refused with {cmd:poisson} (whose replications run on sampling weights) and with {cmd:bootstrap(bayes)} (whose continuous Dirichlet mass is not a duplication count), because neither can produce a replication standard error of the required form; the percentile, bias-corrected and Wald families remain available there. What the help does {ul:not} claim: the replication standard error is the same functional as the observed one, which is what a percentile-t construction requires; whether the resulting intervals attain nominal coverage has not been measured. {p_end}

{phang2}{ul:Costs and limits.} The replications are no longer run under {it:fastnose}, since each one has to estimate a variance to supply its own standard error, so a {cmd:boottstat} run is materially slower than a plain bootstrap of the same {cmd:reps()}, and with {cmd:cluster()} it also pays the physical re-identification leg (see {cmd:bootstrap()}). The p-value is a multiple of 2 / (1 + {it:completed}) with a floor of 2 / (1 + {it:completed}), so it is never exactly 0 and {cmd:reps()} sets the resolution. A replication whose standard error for a given quantity is missing or non-positive is dropped from that quantity's studentized distribution; the row {it:n} of {it:e(boot_tstud)} reports how many were used per coefficient, and a caution names any coefficient that used fewer than half of the completed replications. {bf:Under clustering at or above the} {it:i} {bf:level the analytic pairwise-contrast standard errors can be degenerate} (missing, or numerically zero with a p-value of 0 - the same limitation named in the caution {cmd:jointtests} prints under analytic clustered standard errors), and {cmd:boottstat} studentizes each contrast by exactly that standard error in every replication, so the studentized contrasts and the across-subgroup equality test can be absent or numerically unstable in that configuration; the non-studentized bootstrap forms the contrasts from the coefficient draws and is the one to use for them there. {p_end}

{phang}{opt bootd:f(spec)}: {bf:Requires} {cmd:bootstrap()}, and is refused with an error without it. The reference distribution for the reported bootstrap statistics. {cmd:bootdf(normal)} uses standard normal critical values and reports z statistics - the convention Stata's own {helpb bootstrap} prefix follows. {cmd:bootdf(clusters)} instead uses the t distribution with {it:G}-1 degrees of freedom, where {it:G} is the number of {cmd:cluster()} groups, and a positive integer sets the degrees of freedom directly; {cmd:bootdf(clusters)} without {cmd:cluster()} is refused rather than resolved to t({it:N}-1). {bf:The default depends on whether you clustered}: {cmd:bootdf(clusters)} when {cmd:cluster()} was given, and {cmd:bootdf(normal)} when it was not. The bootstrap standard error is estimated rather than known, so where a cluster count exists the t({it:G}-1) correction is the appropriate one; with no {cmd:cluster()} the resampling unit is the observation, and a {it:G}-1 rule would be {it:N}-1 and indistinguishable from the normal. An explicit {cmd:bootdf()} always overrides the default. This changes only the reference distribution, never the bootstrap standard errors themselves: the Wald-family p-values and normal-theory confidence intervals widen; the percentile and bias-corrected intervals are quantile-based and are unaffected, as are the joint tests, which are formed from the replication statistics (studentized under {cmd:boottstat}, on the fixed bootstrap-covariance metric otherwise). The t option is worth considering when the number of clusters is modest: simulations on 24-cluster designs found coverage running below nominal even where the bootstrap standard error itself was calibrated, because the standard error is estimated rather than known, and a t({it:G}-1) critical value is the standard correction. Under every setting {it:e(df_r)} is aligned with the reference distribution actually reported - the same degrees of freedom under {cmd:bootdf(clusters)} or an explicit integer, and missing under {cmd:bootdf(normal)} so that the fallback is chi-squared and z - so a {cmd:test} or {cmd:lincom} run by hand afterwards reproduces the reported p-values rather than diverging from them - {bf:this holds under} {cmd:bootdisplay(wald)}{bf:, which is what} {it:e(V)} {bf:and} {it:e(df_r)} {bf:always carry}. Under the default {cmd:bootdisplay(percentile)} and under {cmd:bootdisplay(bc)} the printed standard error is still the bootstrap standard deviation that {it:e(V)} holds, but the printed p-value and interval are rank-based, so a {cmd:test} or {cmd:lincom} on {it:e(V)} reproduces the Wald family rather than the printed p; the run's footer says so explicitly. Under {cmd:boottstat} this option has no effect on any reported p-value or interval - every one of them, for coefficients, contrasts, joint and adversarial statistics alike, is a rank among the studentized draws - and it only aligns {it:e(df_r)} as described. The analytic value is preserved as hidden {it:e(df_r_analytic)}. {p_end}
{phang}{opt bootdisp:lay(spec)}: {bf:Requires} {cmd:bootstrap()}, and is refused with an error without it. Which p-value and confidence interval the printed tables carry. Every family is computed and stored in {it:e()} on every bootstrap run regardless; this option chooses only what is displayed, and it does {ul:not} change the reported standard error, which is always the standard deviation of the bootstrap distribution. {cmd:bootdisplay(percentile)}, {bf:the default}, takes the p-value and interval from the percentile ranks of the bootstrap distribution. {cmd:bootdisplay(bc)} uses the bias-corrected percentile method, shifting the interval's percentiles by the median-bias correction {it:z0} that {it:e(boot_z0)} reports and reporting the p-value dual to that shifted interval. {cmd:bootdisplay(bca)} uses the bias-corrected and accelerated method (Efron 1987; the family Stata's own {helpb bootstrap} prefix also reports): the same {it:z0} plus an acceleration {it:a} per coefficient, estimated from a delete-one-cluster jackknife - {it:G} extra fast fits, one per resampling cluster, run after the replications - stored in {it:e(boot_accel)}, with the intervals in {it:e(boot_ci_bca)} and the dual p-values in {it:e(boot_p_bca)}; because of that jackknife cost the BCa family is computed only when it is requested, unlike the other three. {cmd:bootdisplay(bca)} requires {cmd:cluster()}: the jackknife deletes one resampling unit per fit, which without {cmd:cluster()} would be one fit per observation. The pairwise subgroup contrasts get a BCa twin from the jackknife contrasts; the {it:gteffects} and {it:e(cfx*)} cells carry the bias-corrected family under {cmd:bootdisplay(bca)} (a caution says so), since the acceleration would need a jackknife of every cell. {cmd:bootdisplay(wald)} instead reports the p-value and interval from a normal or t approximation, the choice between the two governed by {cmd:bootdf()}. {bf:Note that a percentile or bias-corrected interval is not symmetric about the estimate}, so the estimate plus or minus a critical value times the reported standard error will not in general land on the reported interval under those two settings; that is the ordinary reading of a bootstrap table rather than an inconsistency, and it is only under {cmd:bootdisplay(wald)} that the standard error, statistic, p-value and interval all reproduce one another. Under the percentile and bias-corrected displays the statistic column is therefore labelled {it:t_wald} ({it:z_wald} under {cmd:bootdf(normal)}): it is the estimate over the bootstrap standard deviation, a Wald object whose tail probability is not the rank-based p printed beside it, and it is kept for reference rather than as the test statistic of that row. The percentile and BC p-values use plain order-statistic counts with no add-one adjustment, which keeps each one dual to its own interval - {cmd:p} at or below {it:alpha} exactly when that interval excludes zero, up to a boundary case at even replication counts where the interval rule averages two adjacent order statistics and no count-based p can match an averaged endpoint (measured: 6 disagreements in 35,646 constructed cases, against 36 for the add-one alternative). Both can print as exactly 0, meaning no completed replication fell on the far side of the null, with resolution 2/{it:n}, so {cmd:reps()} sets how small a p is reachable - when it happens a line under the table and the run's footer say that the 0 means p below 2/{it:n}, while the cell and {it:e()} keep the 0; neither moves with {cmd:level()}, while the intervals do. A neighbouring edge case withholds the bias-corrected family for a coefficient outright: when every completed replication of that coefficient falls on one side of the rep-0 comparator, the median-bias correction {it:z0} is {it:invnormal(0)} or {it:invnormal(1)} and is undefined, so {it:e(boot_z0)}, {it:e(boot_ci_bc)} and {it:e(boot_p_bc)} are missing in that column - and the BCa twin with them, which starts from the same {it:z0}. The percentile and Wald families are unaffected. A Caution names the affected coefficients when it happens and {it:e(boot_bc_absent)} lists them afterwards; under {cmd:bootdisplay(bc)} those rows print a blank p-value and interval rather than one borrowed from another family. A low {cmd:reps()} makes a one-sided draw set likelier, but it is not the only cause and raising {cmd:reps()} is not always the remedy: a replication distribution that sits away from the estimate rather than scattered about it produces a one-sided draw set at any replication count. {bf:With} {cmd:boottstat} {bf:only the default} {cmd:bootdisplay(percentile)} {bf:is available}: the studentized draws are summarised by their equal-tailed ranks. {cmd:bootdisplay(wald)} is refused because there is no bootstrap standard deviation of a studentized statistic for a normal or t approximation to rest on, and {cmd:bootdisplay(bc)} because the percentile-t interval is already second-order accurate and a bias-corrected shift of its studentized quantiles moves the interval against the median bias. {cmd:boottstat} chooses {ul:what} is bootstrapped; {cmd:bootdisplay()} chooses how the coefficient draws are read when nothing is studentized. {p_end}

{phang}{opt bootsa:ve(name)}: {bf:Requires} {cmd:bootstrap()}, and is refused with an error without it. Also write the bootstrap replication draws to {it:name}{cmd:.dta} in the current working directory (a plain file name; an existing file is overwritten only with {cmd:replace}; both rules are checked before the replications run). The file has one row per replication, {it:_rep} = 1 to {cmd:reps()}, with the replication coefficients under their {it:e(b)} names (a replication that failed keeps its row with every value missing, so {it:_rep} indexes the draws exactly as {it:e(bootdist)} does after its missing rows are dropped), plus a row {it:_rep} = 0 carrying the observed estimates, so the draws can be recentred from the file alone. Under {cmd:boottstat} it also carries each replication's analytic standard errors as {it:se_*} (the observed analytic standard errors on the {it:_rep} = 0 row) and, with {cmd:subgroups()}, the replication pairwise-contrast standard errors as {it:pwse}{it:#}{it:_*}, numbered by column of {it:e(bootdist_pwse)} (missing on the observed row). With {cmd:cluster()} it also carries the resample itself: one {it:w_}{it:value} column per cluster (map-frame order, named by the cluster variable's value) holding that cluster's multiplicity in the replication under {cmd:bootstrap(standard)} - the number of times it was drawn, 0 when it was not, summing to the number of clusters - or its Dirichlet multiplier under {cmd:bootstrap(bayes)} (the multipliers sum to 1 per replication); the {it:_rep} = 0 row holds 1 under {cmd:bootstrap(standard)} and 1/{it:G} under {cmd:bootstrap(bayes)}, the observed run's weights on the same scale. Under {cmd:bootalpha(size)} the column holds the per-observation multiplier actually applied - the Dirichlet draw divided by the cluster's observation count {it:n_c} - so it is the rows, sum over c of {it:n_c} times {it:w_c}, that sum to 1, and the {it:_rep} = 0 row holds 1/{it:N}, {it:N} the number of observations in the estimation sample. Without {cmd:cluster()} the resampling unit is the observation and no such columns are written. Under {cmd:boottstat} the replication joint-test statistics are added as {it:jointW_}{it:name} and the adversarial statistic as {it:advW}, with the observed statistics on the {it:_rep} = 0 row. The file's characteristics record the engine, the seed and the random-number state (see {cmd:rseed()}). It is written from the accumulators, so it is not subject to the 10,000-replication cap on the posted {it:e(bootdist)}. {p_end}


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


{marker pvalues}{...}
{dlgtab:Which p-value is which, and the floor it carries}

{pstd}
{it:wooldid} prints p-values from four different constructions, and they do not all obey the same rule. The
split is deliberate: the add-one Monte-Carlo-null form is used wherever the observed statistic genuinely
belongs to its own reference set, and the plain share is used where it does not. This table is the map. In
it, {it:n} is the number of {ul:completed} replications - not the number requested - and {it:s(r)} is the
{it:r}th draw of whatever statistic the row reports. {p_end}

{synoptset 40 tabbed}{...}
{synopt:{bf:surface}}{bf:convention, and the smallest p it can print}{p_end}
{synopt:{it:analytic rows}} A normal or {it:t} tail on the reported standard error. No floor and no grid; a p can underflow to 0 but is not bounded away from it by a draw count. {p_end}
{synopt:{cmd:impvalues}} A {it:t} tail on {it:G}-1 degrees of freedom (Ibragimov-Muller). Same: no draw-count floor. {p_end}
{synopt:{cmd:permute()}, coefficient rows} {bf:Add-one, equal-tailed}: min(1, 2 min({it:P_lo}, {it:P_hi})) with {it:P_lo} = (1 + #{c -(}{it:s(r)} {ul:<} {it:s_obs}{c )-})/({it:n} + 1) and {it:P_hi} its mirror. {bf:Floor 2/({it:n} + 1)}; it cannot print 0. {p_end}
{synopt:{cmd:permute()} + {cmd:permabsp}} The same add-one rule on absolute values: (1 + #{c -(}|{it:s(r)}| {ul:>} |{it:s_obs}|{c )-})/({it:n} + 1). {bf:Floor 1/({it:n} + 1)}; it cannot print 0. {p_end}
{synopt:{cmd:permute()} + {cmd:permside()}} The one-sided add-one tail, {it:P_lo} or {it:P_hi} alone. {bf:Floor 1/({it:n} + 1)}; it cannot print 0. {p_end}
{synopt:{cmd:permute()} + {cmd:permrawp}} {bf:The plain share}, #{c -(}.{c )-}/{it:n}, which is what {it:ritest} and {helpb permute} report. It {ul:can} print exactly 0, and 0 there means "no draw was this extreme", not "p is zero". {p_end}
{synopt:{cmd:permute()} + {cmd:permexact}} Plain shares over the {it:M} enumerated assignments, because under a complete enumeration the observed assignment {ul:is} one of the {it:M} - the add-one is already there. {bf:Floor 1/{it:M}}. An element absent from the observed assignment's own replication falls back to the add-one form, and the run says which. {p_end}
{synopt:{cmd:permute()}, joint and equality tests} One-sided upper tail, add-one: (1 + #{c -(}{it:W(r)} {ul:>} {it:W_obs}{c )-})/({it:n} + 1) on a fixed metric. {bf:Floor 1/({it:n} + 1)}. {p_end}
{synopt:{cmd:bootstrap()}, {cmd:bootdisplay(percentile)} and {cmd:(bc)}/{cmd:(bca)} coefficient rows} {bf:The plain share}, min(1, 2 min({it:m0}, {it:m1})/{it:n}) - no add-one. {bf:Floor 2/{it:n}}, and it {ul:can} print exactly 0. {p_end}
{synopt:{cmd:bootstrap()} + {cmd:bootdisplay(wald)} coefficient rows} A normal or {it:t} tail on the bootstrap standard deviation. No draw-count floor. {p_end}
{synopt:{cmd:bootstrap()} + {cmd:boottstat} coefficient rows} {bf:Add-one}, 2(1 + min({it:m0}, {it:m1}))/({it:n} + 1) on the studentized draws. {bf:Floor 2/({it:n} + 1)}; it cannot print 0. {p_end}
{synopt:{cmd:bootstrap()}, joint, equality and adversarial tests} One-sided upper tail, add-one, under either the fixed-metric or the studentized form: (1 + #{c -(}{it:W(r)} {ul:>} {it:W_obs}{c )-})/({it:n} + 1). {bf:Floor 1/({it:n} + 1)}. {p_end}

{pstd}
{bf:Why the split, and why it is not an inconsistency.} Under a permutation null the observed assignment is
exchangeable with the permuted ones, so it belongs in the reference set and counts itself. That is what
makes the p exactly valid at {ul:any} draw count, where the plain share over-rejects - measured at six
percent against a nominal five with ninety-nine draws on the studentized tail. The same argument carries to
every statistic whose draws are null draws: the fixed-metric joint tests under both engines, and the
studentized bootstrap tail, whose draws are recentred at the observed estimate and so form a null reference
set the observed statistic is a member of. It does {ul:not} carry to the bootstrap percentile and
bias-corrected coefficient rows: that distribution is centred at the {ul:estimate}, not at the null, so
adding the observed value to it would be a conservatism choice rather than an exactness result. There the
plain share stands, which is also what {helpb bootstrap}'s own {it:estat bootstrap} reports. {p_end}

{pstd}
{bf:Reading a p of exactly 0.} Only the two plain-share families can produce one: {cmd:permrawp} and the
bootstrap percentile and bias-corrected rows. It means no draw reached the observed statistic, which at
{it:n} draws is a statement about resolution, not about the evidence - the honest reading is "p {ul:<}
1/{it:n}" (or 2/{it:n} two-sided). Every add-one family reports that resolution as its floor instead, which
is why those families never print 0. {p_end}

{pstd}
{bf:The floor is a property of the run, not of the data. Raising the draw count is the only thing that lowers it, and the two engines differ on what they will accept: {cmd:bootstrap()} imposes no minimum at all ({cmd:reps()} need only be positive, and a run is refused only when no replication completes), while {cmd:permute()} requires {cmd:permreps()} to be at least 20 and additionally refuses a run in which fewer than {it:max(20, half the requested draws)} complete. A {cmd:permexact} enumeration is exempt from that completion floor, because the shares stay exact over whatever completed and what incompleteness costs there is resolution.} {cmd:reps()} and {cmd:permreps()} both print a
caution when the floor they imply sits above .05, naming the number, because at that setting no result the
run can produce would be conventionally significant. Raising the draw count is the only fix; neither engine
has a hard minimum, since a handful of draws is a legitimate smoke test. {p_end}

{marker enginediff}{...}
{dlgtab:Why {cmd:bootstrap()} and {cmd:permute()} do not behave identically}

{pstd}
The two engines differ in what they refuse, what they default to, what they store and what they can even be
asked. Most of those differences are not accidents of implementation: they follow from what a bootstrap and
a randomization test {ul:are}, and harmonizing them would make one of the two wrong. Three mechanical facts
generate almost all of them. {p_end}

{pstd}
{bf:1. Recentring.} Bootstrap draws are sampling draws centred at the observed estimate, so a null reference
distribution has to be {ul:manufactured} by recentring - an asymptotic device. Permutation draws are null
draws by construction, nothing is recentred, and the observed statistic is a member of its own reference
set. This is why the p-value conventions differ (see {help wooldid##pvalues:the p-value map}), why the
permutation interval is a shift interval rather than a quantile interval, and why an exact enumeration under
{cmd:permexact} needs no add-one while sampled draws do. {p_end}

{pstd}
{bf:2. Reweighting.} A bootstrap replication carries duplication or Dirichlet weights, so a replication's
standard error is a {ul:different functional} of the data than the observed standard error is of the
original sample - which is why {cmd:boottstat} has to re-identify duplicated clusters to be a genuine
percentile-{it:t}, and why {cmd:bootstrap(bayes)} cannot studentize at all. A permutation replication
relabels onset over the same rows, so it is the {ul:same} functional, and {cmd:permtstat} costs nothing
structural. {p_end}

{pstd}
{bf:3. The assignment is the randomized object.} Anything keyed on treatment status or onset has no
permutation distribution at all: relabelling onset changes what the quantity {ul:means}, not just its value.
This is why the continuous-treatment lattice types whose bins are functions of the assignment have their
inference columns blanked under {cmd:permute()} and not under {cmd:bootstrap()}, why the randomization unit
must nest the cohort variable {it:i} unless {cmd:permonsetcohorts} lets each replication rebuild it, and why {cmd:pretrendpower()} is refused under {cmd:permute()}: Roth's
diagnostic consumes an event-study {ul:covariance}, and the permutation engine does not produce a covariance
of the estimates - it produces a null distribution of a statistic. That is a difference of kind, not a gap
to be filled. {p_end}

{pstd}
{bf:What this means in practice.} The two engines' exchangeability assumptions are {ul:orthogonal}, not one
strong and one weak. The pairs cluster bootstrap assumes clusters are exchangeable as {ul:sampling} draws
and fails with few or atypical treated clusters. Randomization inference assumes clusters are exchangeable
in {ul:assignment} and fails when onset is systematically related to the outcome path - which is the
identifying assumption the design already needs. Neither is a fallback for the other, and the few-cluster
referral from {cmd:bootstrap()} to {cmd:permute()} is a referral to a {ul:different} assumption, not a
weaker one. {p_end}

{pstd}
Where a difference does {ul:not} trace to one of the three facts - a threshold number, a message wording,
which engine happened to receive a stored object first - it is an accident, and those are being removed
rather than defended. The engine-compatibility table below is the one place both engines' requirements are
stated in the same shape, with the reason attached. {p_end}

{marker supporting}{...}
{dlgtab:Supporting Features}

{phang}{opt hist:ogramcohorteffects}: Produce a histogram of the cohort-specific treatment effects behind the requested estimand and return those effects in {it:e(histogramestimates)} for the user to manipulate. The effects are one per treated cohort {it:i} under every estimand - both margins calls average over {it:i} - and a run requesting both gets both sets of columns: the {it:att_it}/{it:att_itime} column is the cohort-level average of that cohort's period-by-period cell effects, taken under the cell weights, not a distribution over cohort-by-period cells. With {cmd:contreat()} a matching histogram of the average marginal effects is drawn as well and its column ({it:ame_i}, {it:ame_itime}) is appended to the returned matrix. {p_end}

{pmore}Under {cmd:bootstrap()} and {cmd:permute()} the histogram outputs are the one reported surface still taken from the stage-1 analytic covariance rather than from the replications; a caution saying so prints with the results. {p_end}

{phang}{opt save:plots(string)}: Save all graphs generated to .gph files beginning with the name specified in the option. Also saves .dta files beginning with the same
name containing the data used to generate any requested event study plots (one file per subgroup value under {cmd:subgroups()}) and, with {it:makecfxplots}, one {it:name}{cmd:_cfx.dta} holding every plotted continuous-treatment point, written on an ordinary analytic run as well as under {cmd:bootstrap()} or {cmd:permute()}, wherever there is plot data to write; on the analytic path it is written after the results are posted and draws no figure of its own. If that write fails, typically because a file of the name already exists and {opt replace} was not given, the run is not lost: a caution says the {cmd:.dta} is missing, and the figures, the {cmd:.gph} files, the estimates and every stored result are unaffected. Graphs and datasets generated have a suffix appended to the name specifying the exact graph/dataset
generated. Everything is saved to Stata's working directory; the argument must be a plain name stem (paths are rejected by the option's syntax).{p_end}

{phang}{opt replace}: Overwrite existing files when {it:wooldid} writes materials: the .gph and .dta files named by {it:saveplots()}, and the replication datasets named by {cmd:permsave()} and {cmd:bootsave()}. For those two the name is checked before the replications run and an existing file is a hard error (r(602)) rather than a silent skip; {it:replace} is what turns that error off. {p_end}

{phang}{opt sum:marystats}: Compute and present summary statistics of the outcome variable y and, if present, of the continuous treatment variable contreat. Also displays the event-study composition matrix {it:e(escomp)} itself - and, with {cmd:subgroups()}, each per-subgroup {it:e(escomp_s{it:sgcase})} - together with its coverage scalars: the reference-normalization pool's N, treated-unit and onset-cohort counts and att-style weight mass, and, for the entire treated estimation sample, the N and att-style weight mass only (see {it:e(escomp_*)} under Stored Results). None of this composition material is computed, displayed or returned unless event-study horizons are requested ({it:esprelength()} + {it:espostlength()} > 0).
Summary statistics are computed overall, within the treatment group by time period (post-treatment vs pre-treatment), and within the control group; under {it:contreatwithin} the control-group rows are omitted, and under {it:contreatcontrols(untreatedZero)} the continuous treatment's pre-treatment and control rows are omitted as well. If a
fixed reference period is not being used and pre-treatment relative time period effects are requested, statistics on the size of the reference group relative to the rest of the
pre-treatment period will also be shown. {p_end}

{phang}{opt pois:son}: Implement the specified {it:wooldid} command using poisson regression, achieved by switching from the regression command {it:reghdfe} to {it:ppmlhdfe} (or, where {it:unconditionalse} or {it:customvce()} already routes estimation through {it:regress} on explicit dummies, from {it:regress} to Stata's own {it:poisson}). Treatment effects
are presented in terms of differences in linear predictions (the poisson {it:predict(xb)} option). Poisson regression is especially useful in situations where one would like to use
log(y) as an outcome, but where y contains natural 0s, as well as when treatment and control groups are only on parallel trends in proportions (e.g., one group is always 2x the other).
Poisson estimation accepts pweights but not aweights (refused with an error). Two further refusals apply: an outcome containing negative values is rejected with r(411), and {cmd:boottstat} may not be combined with poisson (r(198)), because the poisson replications are estimated on sampling weights, so a replication's own standard error is not the same functional of the resampled data as the observed one and the studentized bootstrap would be anti-conservative; the percentile, bias-corrected and Wald bootstrap families remain available. If {it:ppmlhdfe} fails on the initial specification with r(3200), {it:wooldid} retries once with the {it:ccc_absorb()}
and continuous-treatment control terms moved out of the absorbed-effects list into the regression itself (a caution prints); results are equivalent up to numerical
implementation. The retry is keyed on r(3200) alone; any other {it:ppmlhdfe} failure is re-raised unchanged, with no retry. The post-estimation collinearity screen refuses with r(322) when a treated cohort-period dummy has been dropped, and under poisson a dropped dummy whose cell has an outcome of zero on every row is perfect separation rather than multicollinearity: {it:wooldid} says so by name, giving the {it:i} value and the relative-time horizon, and exits r(322). The remedies are to coarsen {it:i}, to drop the cohort, or to specify {it:safety(off)}, which in this case does not merely silence the message: the run proceeds with the separated cohort-period cells excluded from the estimated effects, which is a different realized sample rather than the same one reported more loudly. A pre-existing variable named {it:_ppmlhdfe_d} is overwritten by {it:ppmlhdfe}'s own scaffolding on poisson runs and then removed by {it:wooldid}'s cleanup, except under {it:verbose}, which deliberately retains it alongside the retained full model - treat
that name as reserved in your data. {p_end}

{phang}{opt poisexp:results}: When implementing poisson, present treatment effects in terms of differences in exponentiated linear predictions (the poisson "number of events" or {it:predict(n)}
option). Results of this sort will thus not have a log-type interpretation. The {it:e(gteffects)} cohort-by-relative-period cells are exempt: being stage-1 regression coefficients they stay on the linear-index (xb) scale whether or not {it:poisexpresults} is specified, and a note above that table says so. Specifying {it:poisexpresults} without {it:poisson} is an error, since the option has no effect outside of poisson estimation. {p_end}


{marker eventstudy}{...}
{dlgtab:Event Study Estimation}

{phang}{opt esf:ixedbaseperiod}: Use a fixed reference period for estimation, rather than a pooled reference period.

{pmore}One reading note about the output under this option. The {it:pretreatment} row of the results table averages EVERY estimable pre-treatment period, not the window named by {it:esprelength()}: it is a whole-pre-history placebo, and it does not narrow when you shorten the requested pre-window. The joint pre-trend test on the same run does use the requested window, so the two quantities can cover different sets of periods in one table. Read the row as a placebo over all available pre-treatment time, and the joint test as a statement about the horizons you asked for. {p_end}

{phang}{opt esr:elativeto(-integer)}: The relative time period used as the fixed reference period, if a fixed reference period is being used. The default is -1, the time period immediately
prior to treatment. Under {it:esfixedbaseperiod}, values above -1 (the onset period or later) are rejected with an error: the fixed reference period must lie strictly before onset. Without {it:esfixedbaseperiod}, this option has no effect and a Caution is printed. {p_end}

{phang}{opt espre:length(integer)}: The number of pre-treatment relative time period-specific effects to be estimated. If any pre-treatment effects are requested and
{it:esfixedbaseperiod} is specified, all relative time period effects between the {it:esrelativeto()} period and treatment onset will be presented.

{pmore}Both this option and {it:espostlength()} count CONTIGUOUS INTEGER relative-period slots, not realized horizons. On a panel whose relative-period axis has gaps - a biennial panel, or any design whose onsets sit off the observed time grid - the requested block contains slots no cell realizes, each of which consumes one horizon from the far end of the window. A biennial panel asked for two pre and three post horizons therefore returns one and two, and a realized horizon that the stage-1 regression estimates and that enters the reported {it:att} - which pools every post-onset cell, not only the reported window - can fall outside both the reported event-study table and {it:e(b)}. Widen the request until the window reaches the horizons you want.{p_end}

{phang}{opt espost:length(integer)}: The number of post-treatment relative time period-specific effects to be estimated. With {it:subgroups()} and {it:esprelength()} > 0, a value of
0 is raised to 1: subgroup estimation requires at least one post-treatment horizon. The caution announcing the raise prints either way, since the run is reporting a post-treatment horizon you did not ask for and saying so is a matter of disclosure rather than of safety; what {it:safety(off)} suppresses is only the two lines of accompanying advice. Specify {it:espostlength(-999)} to override the raise and keep zero post-treatment horizons; every other negative value is an error. {p_end}

{phang}{opt joint:tests}: Calculate joint statistical significance tests of all requested pre-treatment relative time period-specific effects. At least two pre-treatment horizons must fall in the estimation window; with only one the option is accepted, no joint test is produced, and nothing is printed to say so. Under {it:esfixedbaseperiod} with {it:timetrends} and the default {it:esrelativeto(-1)} the requirement is one period deeper: the estimation sample must realize at least three pre-onset periods rather than two, so two requested pre-treatment horizons are not enough, and a window that is short by this rule is dropped in the same silent way. An {it:esrelativeto()} deeper than -1 lifts the extra requirement. Caution: simulations suggest these
tests tend to be particularly sensitive to whether you have a sufficient number of clusters in your sample for standard error estimation. Two further reliability conditions apply.
First, {it:jointtests} is unreliable under wooldid's analytic clustered SEs when clustering is at the level of the cohort variable {it:i} or above: the joint test p-values are
severely oversized in that configuration (the cluster-robust covariance cannot capture the within-cluster cross-horizon variation the tests depend on, a problem that does not
diminish with more cohorts or larger cells), so {it:wooldid} prints a caution and computes the test anyway - under both safety settings - reporting the p-value as-is; it should
not be read as a calibrated test. If joint-type inference is needed in such settings, use bootstrap-based clustered standard errors or randomization inference
(see {it:permute()}), or cluster strictly finer than the cohort-by-period cells. Second, in the continuous treatment case the joint tests are only correctly sized when the i-by-t cells each contain many observations: the tests aggregate
within-cell dose-slope variances that carry a downward small-sample bias which compounds across horizons, so designs with only a handful of observations per cohort-period cell
yield oversized joint tests even under robust or finely clustered SEs (note that {it:lattice()} refines cells smaller and {it:contreatpoly} adds within-cell parameters, both of
which sharpen this concern). Two options refuse {it:jointtests} outright: {cmd:impvalues} and {cmd:fastnose} (fast mode reports point estimates only; that refusal stands down on the fast path {cmd:permute()} takes by default). Under {cmd:bootstrap()} and {cmd:permute()} the analytic joint tests are not computed at all: the reported joint pre-trend p-value is a bootstrap or permutation test instead, and where one cannot be formed it is reported as absent rather than as an analytic test. {p_end}

{phang}{opt adv:ersarial}: Tests the composite null hypothesis that all estimated relative time period-specific effects fall on a single line together with the reference period - some line through the omitted base period, with the slope
unknown - a null that is adversarial toward finding pre-trend violations in that the slope is chosen to fit the estimates as well as possible. Under {it:esfixedbaseperiod} each reported coefficient is a difference from the reference period, so an underlying path that is linear implies estimates proportional to (relative period - reference period): one free parameter, not two. Fitting a free intercept as well - which {it:wooldid} did before 2026-09-08 - makes the statistic exactly insensitive to a uniform level shift of every reported coefficient, which is what an anomalous or anticipatory reference period produces. The p-value comes from a calibrated
minimum-distance test on the event-study coefficients and their full variance-covariance matrix (the same objects returned in {it:e(b_es_type)} and {it:e(V_es_type)}): the nuisance
line is estimated by generalized least squares under a generalized inverse of that covariance, and the minimized quadratic distance is referred to an F distribution whose numerator
degrees of freedom equal the residual rank (the rank of the event-study covariance minus one, the single line parameter identified), with the regression's residual degrees of freedom
in the denominator (a chi-squared reference is used when no residual degrees of freedom are available). The p-value requires {it:esfixedbaseperiod}: under the pooled reference regime
no single omitted period exists, so a caution is printed and the fitted line is reported without a p-value. Under {it:timetrends} two periods are normalized to zero (the reference period and the one immediately before it), so the only line passing through both of them is the zero line, and the primary null collapses from "the coefficients lie on some line through the reference period" to "every coefficient is zero". There is then no free line parameter to estimate: the statistic is the whole quadratic form on the event-study covariance and its numerator degrees of freedom are that covariance's rank, one more than the rank-minus-one rule that governs every other run. A p-value is produced whenever the covariance has positive rank, which is the ordinary case under {it:robust}; it is {it:cluster()} at the level of {it:i} that typically leaves it entirely zero. The test under {it:timetrends} is therefore a strictly stronger one than this option's description, and a rejection there is not by itself evidence that the effects depart from a line. If the event-study covariance is degenerate (as is typical under
{it:timetrends}), a caution is printed and no p-value is computed. For display, wooldid also reports a line fit by the unweighted OLS regression of {it:coef_rt} on {it:time_rt} (over the pre-treatment horizons only under {it:altmethod_adversarial}),
where {it:coef_rt} is the estimated coefficient in relative time period {it:rt} and {it:time_rt} is the relative time period; the variable {it:coef_rt} includes values of 0 for the
reference period (or periods, if {it:timetrends} and {it:esfixedbaseperiod} are specified). These tests will not be conducted if an insufficient number of relative time
period-specific treatment effects are requested (either overall, or within specific subgroups). The conditions are these, and each of them leaves the rest of the run intact and prints a caution rather than stopping it. The requested window must hold at least one pre-treatment and one post-treatment horizon, and under {it:timetrends} at least three event-study coefficients in all, the reference period counting as one of them. The pre-and-post rule is then applied a second time to the horizons the estimation sample actually realizes, so a window that was adequate when the options were parsed can still lose the test after the sample is formed. With {it:altmethod_adversarial}, {it:esfixedbaseperiod} and the default {it:esrelativeto(-1)}, a single pre-treatment horizon leaves nothing but the reference row from which to fit the pre-period line, and the test is cancelled. Under {it:subgroups()} the test is skipped for any subgroup that lost a relative period as inestimable (its coefficients no longer span a complete relative-time grid) or that is left with too few of them.
 Under {cmd:permute()} the adversarial p-value is a permutation p-value of the line-fit statistic on the permutation-covariance metric and tests the sharp null of no effect at any horizon (see {cmd:permute()}). Under {cmd:bootstrap()} the adversarial p-value is produced only with {cmd:boottstat}, which bootstraps the test statistic itself; without it the p-value is absent rather than reported as an analytic test, and {it:e(boot_adv_absent)} is set. {it:adversarial} is refused with an error by {cmd:impvalues} and by {cmd:fastnose} (whose refusal stands down on the fast path {cmd:permute()} takes by default). {p_end}

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
design. Magnitudes, from 200 replications per cell under a true null on one cohort-by-period design (30 cohorts, 12 periods, the statistic in force since
2026-09-08): with six clusters spanning each treated cell the test rejected 4.5 percent of the time at a nominal 5 percent, and the permutation p on the same draws
rejected 3.5 percent; with two clusters per cell it rejected 18.0 percent, and with each cell inside a single cluster 72.5 percent. The threefold degradation described
above is therefore nearer fourfold at two clusters, and the single-cluster case is worse than an order of magnitude. Under a continuous treatment on the same design the
{it:att} line test rejected 5.0 percent, while the {it:ame} line test - the one the cell-size condition concerns - rejected 40.2 percent with six observations per cell
and 11.0 percent with fifty, its reported standard errors running about a fifth below the sampling standard deviation at the smaller size. Fifty observations per cell
improves the {it:ame} line test greatly without making it exact. {p_end}

{phang}{opt altmethod_adversarial}: Modifies {it:adversarial} and does not stand alone: specifying it without {it:adversarial} is an error. It is an alternate method of implementing the test above, which projects a trend forward using the pre-treatment coefficients only (plus the planted
zero-valued reference row or rows), and then tests whether the post-treatment coefficients collectively fall in line with the projected values. The test propagates the full
forecast covariance: writing the post-treatment coefficients as {it:b_Q}, the pre-treatment fit set as {it:b_P}, and the linear operator taking {it:b_P} to the forecasts as {it:A},
the forecast error {it:b_Q - A b_P} is tested against its complete covariance {it:V_QQ + A V_PP A' - A V_PQ - V_QP A'}, so uncertainty in the fitted pre-period line and the
pre/post covariance both count against the null. Like the primary implementation, the p-value requires {it:esfixedbaseperiod}. On a relative-period axis with gaps it additionally requires at least one real pre-treatment coefficient beyond the planted reference row or rows; where none remains, the adversarial test is cancelled with a caution (for that subgroup, under {it:subgroups()}). {p_end}

{phang}{opt pretrendpow:er(numlist)}: Reports Roth (2022)-style pre-trend power diagnostics for each requested power level (specified in percent as integers strictly between 5 and 100,
e.g. {it:pretrendpower(50 80 90 95)}): a GLS estimate of the linear pre-trend slope with its standard error and confidence interval (fit through the reference period on the pre-treatment
event-study coefficients using their full variance-covariance matrix from {it:e(V_es_{it:type})} — a more efficient and validly-inferenced line than the display-oriented
{it:e(adversarial_line_{it:type})}, which is an unweighted fit over all horizons, or over the pre-period horizons only under {it:altmethod_adversarial}); the minimum linear pre-trend each pre-test would detect with the requested power ({it:mdt_wald} powers
the joint pre-trend test, an F test with the regression's residual degrees of freedom; {it:mdt_nis} powers the test that any pre-period coefficient is individually significant — the
convention of Roth's {it:pretrends} package, whose size is also reported since it exceeds its nominal level when there are several pre-periods — both pre-tests are conducted at the 5% level and, unlike the slope interval, neither follows {cmd:level()}); and the bias each such trend would imply
for that block's pooled estimand {ul:if it continued into the post-treatment period} (the minimum detectable trend itself depends only on the design's precision, never on the estimated
coefficients). When only {it:att_i} or {it:att_itime} is requested, the reported block is the shared {it:att_it} event-study curve ({it:ame_it} under {cmd:contreat()}) standing in for them: its
hbar and bias columns describe {it:att_it}-family post-period weighting only and do {ul:not} quantify the bias such a trend would imply for the requested {it:att_i}/{it:att_itime}
(or {it:ame_i}/{it:ame_itime}) estimands, whose post-period weights differ. With {cmd:subgroups()}, diagnostics are computed for the overall event-study types only ({it:att}, {it:att_it},
{it:ame}, {it:ame_it}); subgroup event-study objects ({it:e(b_es_*_s#)}/{it:e(V_es_*_s#)}) are not analyzed and no {it:e(roth_*_s#)} results are returned. Requires {cmd:esfixedbaseperiod} (under the pooled reference the single omitted reference period the linear-trend normalization needs is not defined) and at least one
pre-treatment event-study coefficient; not compatible with {cmd:impvalues}, {cmd:fastnose}, or {cmd:permute()} - each is refused with an error, since none of the three produces the event-study variance-covariance matrix the diagnostics are computed on. Under {cmd:bootstrap()} the diagnostics do run, and that covariance is then the bootstrap one rather than the analytic one; a family whose bootstrap covariance is unavailable is skipped with a caution rather than computed on the analytic matrix. Results are returned in {it:e(roth_power_{it:type})} and companion scalars (see Stored Results) and skipped
with a caution when the pre-period variance-covariance block is unavailable (e.g. under {cmd:timetrends}) or carries missing or non-positive variance entries; rank-deficient blocks are
handled with a generalized inverse and reduced degrees of freedom, with a caution (note that on an exactly rank-deficient block every reported statistic reflects only the linearly
independent subset of pre-period coefficients the generalized inverse retains - not a summary of all pre-period evidence - and depends on which columns are retained, so treat those
numbers as one representative of an underdetermined family and interpret with care). A requested power level at or below the NIS test's actual size is achievable by a zero trend, so {it:mdt_nis} and {it:bias_nis} are exactly 0 there (a missing cell always means the
computation was unavailable, never a boundary answer). The exact-quadrature {it:mdt_nis} column is computed only for up to 4 pre-periods (its cost grows
steeply with dimension); for larger designs, feed {it:e(b_es_{it:type})}/{it:e(V_es_{it:type})}/{it:e(numpre_{it:type})} to Roth's {it:pretrends} package (subject to the same positional-geometry conditions as the honestdid handoff — see the {it:e(b_es)} Stored Results entry), which implements the same
convention with a normal critical value. When the underlying estimator posts no residual degrees of freedom (e.g. the
ppmlhdfe path), the joint test and minimum detectable trends use chi-squared rather than F distributions and the slope CI uses a normal critical value. Under {cmd:poisson} without
{cmd:poisexpresults}, the diagnostics are on the linear-index (xb) scale like the event-study coefficients themselves. {p_end}

{phang}{opt makep:lots}: Produce event study plots. No plot is drawn unless at least one event-study horizon is requested via {it:esprelength()} or {it:espostlength()}. Note on graph management (applies to all of {it:wooldid}'s plotting options): {it:wooldid} names the graphs it draws
(e.g. att_fullmodel, hist_att_i, cfxatt_cfxpost) and creates them with the graph-name {it:replace} suboption, so re-running overwrites {it:wooldid}'s own prior graphs of the
same name; it does not clear other graphs from memory. Two consequences: a user graph whose name exactly matches a {it:wooldid} graph name will be overwritten when that plot
is recreated, and {it:wooldid} graphs from an earlier, differently-configured run linger in memory until that configuration is run again (or the user drops them). {p_end}

{phang}{opt oos:treatedcontrols(integer)}: By default, {it:wooldid} uses cohorts treated after the end of the sample the same way as never-treated control cohorts (i.e., such cohorts are
"out-of-sample treated controls"). When this option is specified with a value of T, cohorts treated within T periods of the end of sample are instead handled the same as other treated cohorts.
Such cohorts are potentially allowed to contribute to pre-treatment relative time period-specific effects when this option is specified,
but not otherwise. Note that if {it:esfixedbaseperiod} is specified, the out-of-sample-treated cohorts will be dropped if the specified reference period is not usable for these cohorts
(not observed within the sample, or removed because no comparison unit exists in that calendar period). When {it:anticipation()} is also specified, the out-of-sample window applies to the {ul:effective} (shifted) onset, so a cohort whose recorded onset
lies just beyond the window can become usable once its effective onset falls within it. Negative values are not allowed. {p_end}

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
The SG variable must contain only non-negative values (or missing values) on in-sample observations of treated units - both pre- and post-onset rows - and only integers (or missing values) on those rows at or after the unit's onset, which are the rows whose values become the subgroup levels; values on never-treated observations and on rows wooldid's deterministic screens discard (e.g. missing outcomes) are ignored by the check, since no subgroup computation uses them and its levels must specify collections of cohorts i. Under {cmd:permute()} these screens, and the containment screens described below, instead run over every in-sample row rather than over treated units' rows only, because any cluster can be assigned treatment in a replication (see the {cmd:permute()} entry). A given cohort i cannot be split across multiple
subgroups: {it:wooldid} exits with an informative error if it is. Under {it:safety(off)} such a run proceeds using share weighting (see the subgroup discussion below). Observations given missing values are not placed in any subgroup and are ignored when performing subgroup computations. When {it:esprelength()} is greater than zero and {it:espostlength()} is 0 - its default, and equally when 0 is requested explicitly - {it:wooldid} raises {it:espostlength()} to 1, because subgroups for which pre-treatment effects are not estimable can otherwise crash the run; the change prints a caution under default safety and is made silently under {it:safety(off)}. Specify {it:espostlength(-999)} to keep the post-treatment window empty. Because subgroup values are
embedded in the names of stored results and plots, values 8 or more digits long - 7, with {cmd:makecfxplotsbysg} and 3-digit relative-year plot requests - can hit Stata's 31/32-character name limits and crash some option combinations with r(198); recode
ID-style subgroup codes to small integers first (e.g., {cmd:egen newsg = group(oldsg)}). {it:impvalues} is refused outright alongside {it:subgroups()}: Ibragimov and Muller inference over subgroup estimates is not implemented, and the combination exits with r(198) rather than dropping either request.   {p_end}

{phang}{opt suppress:sgprimaryeffects}: Skips estimation of subgroup-specific treatment effects and computes only the differences in effects between SGs. This option should
speed estimation if these treatment effects are not otherwise of interest. {cmd:bootstrap()} and {cmd:permute()} both refuse it, because the suppressed subgroup entries are exactly what their replication vectors need in order to resample the subgroup statistics. It also switches off the rank-deficiency screen described below, which reads the subgroup level standard errors this option never computes: degenerate contrast standard errors are then reported unscreened and with no caution. {p_end}

{pstd} {ul: Description of Approach to Subgroup Effect Estimation:} Given a variable dividing cohorts i into different subgroups, {it:wooldid} will estimate treatment effects separately
for each subgroup (i.e., each level of the SG variable). This is achieved without modification to the underlying
regression specification in step 1 (see {help wooldid##description:here}). Rather, it is done by carrying out the margins calculations in step (2) within the subsamples of the data
specified by the subgroup variable.
Additionally, {it:wooldid} computes all possible pairwise comparisons of the overall
effects (and pre-treatment effects) estimates by SG (e.g., SG 1 vs 2, 1 vs 3, 2 vs 3, etc.), presenting for each comparison the difference in effect estimates and a p-value for each difference.
A joint test of the hypothesis that all SGs share the same effect estimate (normalized to be the effect estimate of the numerically lowest numbered SG) is also computed, except under {cmd:fastnose} (no such test is posted at all) and for any family the rank-deficiency screen below withholds (posted as missing). Note: standard errors
for subgroup-specific estimates can fail when using clustered standard errors in conjunction with subgroups that contain few (or only one) cluster.
{marker sgdegen}{...}
Separately, when the clustering is at the level of the cohort variable {it:i} or above, the analytic covariance of two subgroups' effects can be numerically rank-deficient in the contrast
direction - the two level standard errors agree to several digits and the contrast standard error is the cancellation residue between them, not an estimated variance (the same degeneracy named
in the {cmd:jointtests} caution). {it:wooldid} screens each pairwise row on the implied correlation between the two subgroups' estimates and, above 1 - 1e-2 (or where the contrast standard error is itself missing, the extreme of the same failure), reports that row's standard error,
statistic, p-value and interval and the across-subgroup equality test covering it as {ul:absent}, with a caution naming the count; the point estimates are unaffected. {cmd:safety(off)} reports
them with the same caution instead of withholding them, and the screen is inert under {cmd:suppresssgprimaryeffects}, where the subgroup level standard errors it needs are never computed, under
{cmd:fastnose}, which computes no standard errors at all, and under {cmd:bootstrap()} and {cmd:permute()}, whose replication shells rewrite these cells with their own inference. Use {cmd:permute()}, {cmd:bootstrap()} without {cmd:boottstat} (which studentizes each contrast by exactly that standard error),
or a clustering strictly finer than the cohort-by-period cells. {p_end}

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

{phang}{opt contreat(varname)}: The continuous treatment variable. Must be defined for all observations (set to 0 if not naturally defined within untreated cohort-periods); observations on which it is missing are dropped from the estimation sample, appearing only in the header's dropped-observation count.
With continuous treatments, reliable analytic inference additionally requires that cohort-by-period cells contain a large number of observations relative to the number of
within-cell parameters estimated. If the continuous treatment carries no identifying variation for any treated cohort-period cell's dose slope (for example, a dose constant
within each cohort - every treated-cell slope term is then dropped for multicollinearity), the within-cell average marginal effects are not identified: under default safety
{it:wooldid} names this situation, reports the ame results as missing rather than as the regression's structural zeros, and continues with the unaffected att results
(expected zero doses in untreated periods, as under {it:contreatcontrols(untreatedZero)}, do not trigger this). Only the raw {it:e(b)} coefficient vector retains the structural zeros, since Stata cannot post missing coefficients, and the check is skipped under {it:fastnose}. Not available with {opt gte:ffects}: under a continuous treatment the cohort-by-period interaction coefficients are dose-zero intercepts rather than treatment effects, so that combination is rejected with an error; aggregate {it:e(b_es_ame)} or the cfx materials instead.{p_end}

{phang}{opt contreatelas:ticitytype(dydx/eyex/eydx/dyex)}: Specifies which type of treatment effect to compute. {it:eyex} corresponds with computing the elasticity
(% change in y/% change in contreat), {it:eydx} corresponds with semi-elasticities comparable to what can be produced for the main effect of treatment
(% change in y/unit change in continuous treatment), {it:dyex} gives the other semi-elasticity (change in y/% change in continuous treatment), and {it:dydx} is the default
(not an elasticity). A supplied value outside this list is rejected with an error, whether or not {opt contreat(varname)} is specified; without {opt contreat(varname)} the option is otherwise inert, its value read by nothing and no message printed.
{ul:Under {it:poisson} without {it:poisexpresults}} results are on the linear index {it:xb} = ln E[y|x], where {it:dydx} already gives the semi-elasticity of E[y] and {it:dyex} already gives the
elasticity. {it:wooldid} therefore computes {it:eydx} as {it:dydx} and {it:eyex} as {it:dyex} on that scale, with a caution: the option changes the table header rather than the estimates. Use
{it:poisexpresults} to obtain the same quantities from the exponentiated prediction. Note that successful estimation of these elasticities tends to require a larger sample, as well as y and contreat variables that are not sparsely distributed. {p_end}

{phang}{opt contreatw:ithin}: Specify this when you wish to use only within-treated cohort variation in the continuous treatment variable. This option disables some of {it:safety}'s features (the estimability-relaxation checks; the data-integrity screens, such as onset constancy within units, still run unless {it:safety(off)} is given explicitly), turns off
estimability checks that would otherwise cause {it:wooldid} to stop estimation in the absence of any untreated cohorts, and causes {it:wooldid} to skip sample cleaning routines that drop time
periods for which a control cohort is not present. Requires {opt contreat(varname)}: because this option disables those checks, requesting it without a continuous treatment is rejected with an error rather than silently ignored.   {p_end}

{phang}{opt contreatp:oly(integer)}: Includes polynomial terms in {it:contreat} of the specified degree. The value is validated whether or not {it:contreat()} is given, but without {it:contreat()} the option is otherwise inert: nothing reads it and nothing is said about it, unlike {it:contreatwithin}, {it:lattice()} and {it:contreatcontrols()}, each of which is an error on its own. The degree must be a positive integer (default 1; 0 and negative values are rejected). Be cautious using higher degree polynomials. {p_end}

{phang}{opt lattice(latticevar OR latticetype latticenum)}: Specifies a lattice variable, a variable specifying groups of observations across which the effect of {it:contreat} is
allowed to vary, even within cohort-periods and within the reference period. If just the name of a variable is passed to {it:lattice()}, this variable will be used as the lattice variable.
If instead a lattice type and a lattice number are passed, this option will generate a lattice variable using {it:latticenum} quantiles, where {it:latticenum} must be an integer greater than 1, with the quantiles being formed within groups
according to the {it:latticetype} (see below). Any specified sample weights will be used when constructing the lattice, unless the option {it:latticeignoreweights} is specified (which
constructs the lattice quantiles without using weights). 
If {it:lattice} is not specified, {it:wooldid} will proceed without
using a lattice variable (as though you specified a lattice variable that is always 1); {it:lattice()} without {opt contreat(varname)} is rejected with an error rather than silently ignored. A user-supplied lattice variable must be numeric and contain only non-negative
integers (it is used as a factor variable; the value check can be bypassed with {it:safety(off)}); a lattice type passed without a quantile count is rejected with a
specific error at parse time. Observations on which a user-supplied lattice variable is missing are dropped from the estimation sample, surfacing only in the header's dropped-observation count; a generated lattice instead places null values in bin 0. Be aware that a lattice is the most expensive scaling axis in {it:wooldid}: runtime grows at roughly the 1.9th power of the
cohort-by-period cell count once a lattice is present (versus about the 1.3rd to 1.7th power without one), and the absorbing regression rather than
{it:margins} becomes the dominant cost. Each additional lattice point multiplies the terms the regression must carry, so prefer the smallest
{it:latticenum} that answers the question. See {help wooldid##runtime:Runtime and what drives it}.
When a lattice is combined with {it:contreatcontrols(basic)}, the baseline dose level is always partialled out, within each lattice level: the lattice-level dose control is in the
model, and because it spans the pooled one the default estimation path absorbs that one alone, while the paths carrying these terms as explicit regressors ({it:ccnoabsorb},
{it:unconditionalse}, {it:customvce()}) carry both, the spelling there governing which column the collinearity drop takes. Before 2026-09-05, plain {it:basic} with a lattice inserted no baseline dose control on the default estimation path, so every ame-family
row (the main and pre-treatment ame and each event-study ame horizon) was off by the baseline dose slope; ame results under {it:lattice()} with {it:contreatcontrols(basic)} changed on
that date (the att results did not). {p_end}

{phang} {it:latticetype:} Specifies groups in which to estimate {it:latticenum} quantiles of {it:contreat}. Type names are matched case-sensitively and must be typed exactly as spelled below. Choice of lattice type can be important for both interpretive reasons
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
variable or, if specified, the variable in {it:cfxlattice()}. One of {it:lattice()} or {it:cfxlattice()} must be present: requesting the plots with neither is rejected with an error. These plots are only sure to be valid for the ATT estimand type, but may be valid for others if lattice groupings
do not crosscut cohorts (or, less stringently, cohort-periods, for ATT_it). For that reason {it:wooldid} produces them only for the {it:att} estimand unless {it:safety(off)} is in force; {it:contreatwithin} relaxes safety internally and has the same effect. Requires a continuous treatment ({opt contreat(varname)}), and is not compatible with {opt imp:values}. {p_end}

{phang}{opt cfxlattice(varname)}: A variable containing the bins in which effects are to be calculated for the continuous treatment effects plots. If not specified,
the lattice variable will be used. The variable must be numeric and, on observations with a nonmissing treatment onset within the sample restriction, contain only
non-negative integers ({it:margins} accepts only non-negative integers in {it:over()}; the value check can be bypassed with {it:safety(off)}). Values on never-treated
observations are ignored: no continuous treatment effects plot uses those rows. If the variable has missing values on observations with a nonmissing onset, {it:wooldid}
prints a Caution: {it:margins} silently omits those observations from the plotted bin averages (this Caution, like the value check above, is suppressed by {it:safety(off)}). These screens apply only when {it:makecfxplots} is actually requested; {it:cfxlattice()} without {it:makecfxplots} has no effect and draws a Caution.
Note that {it:makecfxplots} requires a continuous treatment ({it:contreat()}); previously, requesting it with {it:cfxlattice()} alone was silently ignored. {p_end}

{phang}{opt cfxplott:ypes(pre post integers)}: A list of types of cfxplots to be produced. If {it:post} is among the types specified, plots will be produced giving
treatment effects across the continuous treatment variable within the post-treatment period. If {it:pre} is among the types specified, plots will be produced within the pre-treatment period.
If any integers are specified, plots will be produced just within the specified periods relative to treatment. If nothing is specified, the option post is assumed. {it:pre} and
out-of-window relative-period integers are dropped from the request: {it:pre} requires {it:esfixedbaseperiod} or {it:esprelength()} > 0, a positive integer must be no larger than {it:espostlength()}, and a negative integer's magnitude must be no larger than {it:esprelength()}. {it:esfixedbaseperiod} lifts both integer limits but discards the {it:esrelativeto()} reference period itself. Dropping and rejection are different things: a token that is neither {it:post} nor {it:pre} must be a number, and one that is neither ({it:postpost}, {it:2post2}, a variable name) stops the run with an error naming the entry, as does a numeric token that is missing or not an integer. A leading plus sign is accepted and stripped, so {cmd:cfxplottypes(+2)} and {cmd:cfxplottypes(2)} are the same request and the repetition is removed; repeated relative periods are collapsed to one, first occurrence kept. Only the window screens described above drop silently, so a request in which nothing survives completes normally with no continuous treatment effects plot and no message. {p_end}

{phang}{opt makecfxplotsbysg}: Produce continuous treatment effects plots for each group specified in the subgroup variable. It needs both {it:makecfxplots} and {it:subgroups()} and is silently inert without either: no error and no caution, unlike {it:cfxlattice()} without {it:makecfxplots}. Being part of the {it:makecfxplots} output, it carries that option's estimand restriction; see {opt makec:fxplots} above. Each subgroup's plot data comes back in a matrix named {it:e(cfx{it:item}_{it:type}_sg{it:sgcase})}, for example {it:e(cfxatt_cfxpost_sg1)} or {it:e(cfxame_relyr2_sg3)}, where the suffix is {cmd:_sg} followed by the subgroup value. Note that this is the one family of stored results using {cmd:_sg}: the subgroup twins elsewhere use {cmd:_s} ({it:e(escomp_s1)}, {it:e(b_es_att_s1)}) or a bare underscore ({it:e(wdidsgresults_1)}). {it:e(wdid_cfxmats)} lists every cfx matrix name a run actually posted.   {p_end}

{phang}{opt cfxplotnose}: Do not estimate standard errors for the cfxplots; speeds estimation. Without {it:makecfxplots} the option is inert: nothing reads it, and no message is printed, unlike {it:cfxlattice()}, which draws a caution when it is specified on its own. Under {cmd:permute()} the truncation does not apply: the six-column shape of the returned {it:e(cfx*)} matrices is kept, missing-padded because {cmd:margins} returns a short {it:r(table)} when it computes no standard errors, precisely so the permutation rewiring can fill the p-value and interval columns from the replications. A {cmd:permute()} run that specifies {it:cfxplotnose} therefore still gets permutation inference on the plot cells wherever the lattice type is assignment-invariant, and the blanking caution and {it:e(perm_cfx_blank)} apply to it as to any other {cmd:permute()} run. {cmd:bootstrap()} does not carve it out, so there the cells stay analytic. {p_end}

{phang}{opt contreatc:ontrols(controltype)}: Specifies what controls for {it:contreat} will be added: the specified controls are always equal to f(Contreat_jt) (the function used to capture
  the effect of contreat within treated i-t cells) interacted with a set of fixed effects-like indicator variables (time period, cohort, etc.). The specified controls will partial out
the effect of {it:contreat} in a particular way, affecting the nature of the reference effect of {it:contreat} to which estimated effects of {it:contreat} in treated cohort-periods are
compared. This option must be specified to proceed with a continuous treatment (and is rejected with an error if specified without one), and care must be chosen when specifying it: different options correspond with very different
substantive assumptions. See more on this topic in the description section below. {p_end}

{phang} {it:conttreatcontroltype}: untreatedZero, which cannot be paired with any other type; or one or more of the following: basic, byTreated, byTtre, byCohort, and byTime  {p_end}
{synopt : {it:basic}} Controls for the effect of contreat (i.e., for f(Contreat_jt)) are added without any interactions. The baseline dose level is always partialled out, including when a
lattice is present (since 2026-09-05; see {it:lattice()}).  {p_end}
{synopt : {it:byTreated}} Controls for the effect of contreat are interacted with indicators for membership in the treatment group and in the control group. {p_end}
{synopt : {it:byTtre}} Controls for the effect of contreat are interacted with indicators for membership in each treatment onset group and in the control group.  {p_end}
{synopt : {it:byCohort}} Controls for the effect of contreat are interacted with cohort (i) membership indicators. Each treated cohort must therefore have variation in {it:contreat} across its untreated/reference-period observations; a cohort that does not is named and the run stops, with no {it:safety(off)} escape.  {p_end}
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
be the integral of the estimated average marginal effects across the distribution of the continuous treatment variable. Note: the control type {it:byTime} is refused in this case,
since in the absence of an untreated group, the time period interactions will force {it:wooldid} to partial out genuine treatment effects for some group, threatening identification. Specify {it:safety(off)} before the run to proceed anyway, which downgrades the refusal to a caution. {p_end}


{marker controls}{...}
{dlgtab:Controls and Fixed Effects}

{phang}{opt control:s()}: The contents of this option are passed directly to {it:reghdfe} as controls, without any modification (to {it:ppmlhdfe} under {it:poisson}, and to {cmd:regress} or {cmd:poisson} under {it:customvce()} or {it:unconditionalse}). {p_end}

{phang}{opt fe()}: The contents of this option are passed directly to {it:reghdfe} - or to {it:ppmlhdfe} under {it:poisson} - as fixed effects (into its absorb statement) without modification. Variables given will be interpreted as
categorical variables, though i.X##c.X2 type notation is also permitted for passing heterogeneous slope controls and the like.

{pmore}This holds on the default estimation path. Under {it:customvce()} or {it:unconditionalse}, estimation runs through {cmd:regress} - or {cmd:poisson} under {it:poisson} - rather than {it:reghdfe}, and a BARE token there is read as a CONTINUOUS variable: {it:fe(state)} fits a linear term in the state code rather than a set of state indicators. Write the prefix explicitly - {it:fe(i.state)} - whenever those options are in use, or on any specification you may later run both ways.{p_end}

{phang}{opt timetrend:s}: Control for cohort specific time trends in the regression (i.e., for i.i##c.t). Importantly, this option is required for heterogeneous time trend controls to work,
even when they are specified manually, as it also updates {it:wooldid}'s various estimability checks, drops cohorts without sufficient pre-treatment data to estimate a time trend,
and, in the fixed reference period case, ensures that the two reference periods that must be normalized to 0 to estimate the time trend are chronologically adjacent. This option does not
automatically include any time trends that are a function of any included continuous treatment variable.  {p_end}

{pmore}With {cmd:cluster()} on the cohort identifier {it:i}, the cohort trend terms typically leave the covariance non-identified: the estimation command returns an all-zero covariance, {it:wooldid} says so and continues, the reported estimates stand, and every standard error, test statistic, p-value and interval in the run is missing. {it:e(V_es_*)} is then not posted and {it:e(wdid_vcvstate)} reads {it:absent}, so {it:honestdid}, {it:pretrendpower()}, the joint tests and {it:simultaneousbands} have no input. The same specification under {it:robust} is unaffected. {p_end}

{phang}{opt coarse:cohortcontrols()}: Includes the specified controls as well as interactions between them with all Tz_it dummies included in the model, meaning one interaction with each
cohort-time period for which a treatment effect is estimated (if {it:esfixedbaseperiod} is specified, that includes every estimable treated cohort by relative time period effect). Note that
reported treatment effects will not vary as a direct function of these controls; the effect of these controls are partialled out from the reported treatment effects.
This option differs from simply specifying cohort
by time period interactions with a given control in that it pools the effect of the controls in the reference period and the control group. Note that all specified coarse cohort controls must be
specified using factor notation, so that {it:wooldid} knows if they are to be interacted as a continuous variable (prefixed with "c.") or as a categorical variable (prefixed with "i."). {p_end}

{phang}{opt inter:activecontrols()}: Includes the specified variables in the regression in two capacities. The cell slopes exist only where a control varies within cohort-period cells; if it does not (for example {it:i} is the unit and the control is a unit-level variable), every slope term drops out and the run is identical to {it:controls(i.t#c.x)} - {it:wooldid} then prints a caution naming the control (that check does not run under {it:fastnose} or {it:safety(off)}, where the reduction is silent). First, the variables are included as controls interacted with time period-specific
indicator variables. Second, the variables are de-meaned by cohort (i) and then treatment effects are permitted to vary linearly with the de-meaned variables. This is achieved by including
interactions between the de-meaned controls and all Tz_it dummies included in the model. The specification used is similar to as with {it:coarsecohortcontrols}, but differs in that these terms
are marked so that variation associated with them is reported as part of the treatment effects (as opposed to being partialled out from the treatment effects). This specification most closely
mimics the approach to including covariates that appears in Wooldridge (2025). When using only one interactive control, ATTs obtained using this option should be similar to those obtained
when specifying that variable as a continuous treatment variable and then separately including controls that interact the continuous treatment variable with time dummy variables.
The treatment-effect interactions are built for the treated cohort-by-period cells only: the pooled reference cell (untreated observations) carries no interaction slope, so its
de-meaned-control variation is partialled out through the time-period interactions and the reference-period effect is a structural zero. Before 2026-09-05 the reference cell also received a
slope column that was exactly collinear with the time-period interactions plus the treated-cell slopes, so the estimation engine's collinearity drop decided the reported horizons: the
horizon estimates, their standard errors, the pre-treatment joint tests and {it:e(b_es_*)}/{it:e(V_es_*)} were solver- and regressor-order-dependent wherever the control varies within
cohort over time, under weights, or on an unbalanced panel. Results under this option changed on that date. {p_end}

{pmore}This option also puts the run outside the default cell collapse: the {it:margins} stage is evaluated on the full estimation sample, and an explicit {it:fastcollapse} is refused with an error naming the option (see {it:fastcollapse}). {p_end}


{pstd} {ul:Caution Regarding Continuous Treatment Variables and Controls:} If seeking to manually specify additional controls as a function of a continuous treatment variable, be sure not
to include the continuous treatment variable itself among your controls, as this can adversely affect estimates produced in {it:wooldid}'s margins calculation stage. It is okay to include
a control with values identical to that of the continuous treatment variable: it just needs to have a different variable name. It is recommended that any custom such controls be specified
in conjunction with the option {it:contreatcontrols(basic)}, or any option other than {it:untreatedZero}. {p_end}

{pstd} {ul:Caution Regarding Continuous Treatments and Cluster-Robust Standard Errors:} Valid clustered standard errors cannot be obtained analytically for the ame-family
quantities when clustering at the cohort ({it:i}) level or above, necessitating use of bootstrapping or other methods (see {help wooldid##standarderrors:standard_errors} above).
With {cmd:contreat()}, the regression gives every treated cohort-by-period cell its own
treatment-intensity slope. When any such cell lies entirely within a single cluster of the specified {cmd:cluster()} (which is automatic when clustering on the cohort identifier
{it:i} itself), the cluster-robust variance estimator is structurally blind to within-cell variation in the continuous treatment: the "standard errors" it would deliver for the ame-family
quantities do not measure their sampling variability (in the {it:contreatcontrols(basic)}/{it:byTime} cases they collapse to the variance of the absorbed baseline dose control, one number
repeated across every ame quantity — including configurations where the covariance appears full rank). {it:wooldid} detects this design, prints an explanatory caution, sets the ame-family
standard errors, tests, and joint diagnostics (adversarial ame p-values, ame joint pre-trend tests, subgroup ame comparisons, ame pre-trend power) to missing, and withholds
{it:e(V_es_ame*)} - and with it any simultaneous band ({it:simultaneousbands}) for the ame families, which is dropped without a further caution whenever some other family still gets one. The same screen blanks the inference columns of every {it:e(cfx}{it:ame item}{it:_}{it:type}{it:)} matrix and its per-subgroup twins, keeping the point-estimate column alone, and so draws the ame continuous-treatment figures and their saved {cmd:.gph} files with no confidence interval: {cmd:marginsplot} is called with {it:noci} rather than the usual shaded band. The att-side cfx matrices and figures are untouched. Keep this case apart from the dose-slope refusal described under {it:e(b_es_}{it:type}{it:)}, where the point-estimate column goes too and the ame cfx matrices come back entirely missing. The run still returns the unaffected point estimates and {it:e(b_es_ame*)}. The screen runs only on the default analytic path: it is skipped under {it:fastnose}, which asks {it:margins} for no variance at all; under {it:impvalues}, whose cross-cluster dispersion inference stays valid here; and under {it:unconditionalse}, which is a different variance estimand and reports ame standard errors on such a design without this caution. Valid ame inference under such designs: bootstrap standard errors, IM p-values
({it:impvalues} — cross-cluster dispersion inference, which remains valid), randomization inference, {it:robust} (or a robust {it:customvce()}) when errors are independent across observations (not appropriate
under within-unit serial dependence), or clustering on a variable strictly finer than the cohort-by-period cells. Note that binary-treatment (att) standard errors under cohort-level clustering are also somewhat understated on such saturated designs (a milder, documented
phenomenon); wooldid does not suppress them. On some of these designs the estimation command returns an ALL-ZERO covariance outright, in which case no standard error in the run is identified: {it:wooldid} says so and continues, the reported estimates stand, and every standard error, test statistic, p-value and interval is missing. Builds from 2026-08-24 to 2026-09-08 re-estimated the specification through {cmd:regress} with the fixed effects written out as explicit dummies, to salvage the standard errors that survived; that rescue was removed on 2026-09-08, because it re-solved a large collinear system to recover part of an inference the containment had already invalidated, and it delayed the point at which the specification's problem became visible. {p_end}


{marker technical}{...}
{dlgtab:Technical Features}

{phang}{opt update}: If specified, the command {it:wooldid,} {it:update} will update your local installation of {it:wooldid} based on the most recent edition of the program available on github. The new files are written to disk, but a running Stata session keeps using the copy of {it:wooldid} it has already loaded, so the update takes effect only after you restart Stata or run {cmd:discard}. {p_end}

{phang}{opt customw:eightmultiplier(varname)}: When computing various treatment effects using the margins calculation (see step 2 {help wooldid##description:here}), the weights normally 
used to produce the requested type of treatment effect (ATT, ATT_i, etc.) will be multiplied by this variable. When requesting an ATT without specifying sample weights, this is 
equivalent to estimating all regressions without weights, but then using these weights when averaging treatment effects to obtain the ATT. This allows production of certain 
bespoke types of treatment effects. Values must be non-negative on in-sample observations (negative values would silently corrupt the weight masses and estimates; they are rejected
with an error, bypassable with {it:safety(off)}); missing values drop the observation. Under {cmd:impvalues} the reported event effects are equal-weight averages of per-cluster
margins, so the {it:e(escomp)} weight-share columns describe sample weight composition, not the Ibragimov-Muller cross-cluster aggregation. {p_end}

{phang}{opt clean:matrices}: If specified, {it:wooldid} will clear out all matrices storied in STATA's main memory and in mata. Occasionally, {it:wooldid} enters a state where it repeatedly
crashes for difficult to determine reasons (relating to {it:reghdfe} crashing and leaving certain content in mata). Doing this can sometimes recover {it:wooldid} from this problem. {p_end}

{phang}{opt regtol(#)}: Specifies the convergence tolerance to be used by {it:reghdfe} and {it:ppmlhdfe}. The argument is the negative base-ten exponent of that tolerance rather than the tolerance itself ({it:regtol(9)} means {it:tol(1e-9)}), so a larger number is a tighter tolerance and a slower, more precise fit; the {it:wooldid} default of 9 is one digit tighter than {it:reghdfe}'s own 1e-8. It reaches the {it:reghdfe}/{it:ppmlhdfe} path only. Under {it:unconditionalse} or {it:customvce()} the estimation runs through {cmd:regress} or {cmd:poisson} on explicit dummies, neither of which takes a convergence tolerance, so {it:regtol()} is accepted there and has no effect. The value is not range-checked, so a negative argument produces a malformed tolerance and an opaque failure inside the estimation command. By default, {it:reghdfe} sets a value of 8; {it:wooldid} has a default of 9.
You may wish to set this higher or lower; higher may be advisable in settings where you ask {it:reghdfe} to absorb interaction terms with many continuous variables -- as tends
to occur when using continuous treatment variables. {p_end}

{phang}{opt ccno:absorb}: Forces continuous treatment controls to be included in the main regression statement, rather than in the {it:absorb()} option, which is the default. Since the controls it relocates exist only under {it:contreat()} with {it:contreatcontrols()}, the option is inert and silent on a run with no continuous treatment
whenever using the commands {it:reghdfe} or {it:ppmlhdfe} to estimate the underlying regression. Under certain conditions, this can improve estimation stability and avoid problems
with variables being inappropriately dropped from the regression. {p_end}

{phang}{opt ccc_absorb()}: An alternate version of coarsecohortcontrols() that admits only continuous variables and passes the controls to the {it:absorb()} statement when using
a high dimensional fixed effects estimator. This is often less efficient than the alternative approach in terms of computation time, but can be important in terms of memory
efficiency in the main regression or under certain other circumstances. Note: be extra careful to ensure that specifications featuring ccc_absorb() are identified; not all cases
of identification failure in this case are bound to be automatically detected. Under {cmd:contreat()}, naming the treatment variable itself inside {it:ccc_absorb()} puts the run
on the ineligible list for the default cell collapse (see {it:fastcollapse}), which silently turns the collapse off and makes an explicit {it:fastcollapse} an error. It is not refused with an event-study window, but it can collide with one: whenever the {it:ccc_absorb()} variables are constant within a cohort-by-period cell - a unit-level variable with a unit-level {it:i} is the usual case - the absorbed set takes in the cohort-by-period treatment interactions themselves, the fitted model then carries no coefficient at any relative period, and the run stops with an error naming {it:ccc_absorb()} as a likely cause. Where those variables vary within the cell the combination estimates normally. {p_end}

{phang}{opt ver:bose}: Prints the regression command {it:wooldid} uses and prints progress reports throughout the estimation process. The regression itself is run quietly on every estimation path, so its coefficient table is never displayed.
This option also stores the underlying regression model used by {it:wooldid} to {it:___wooldidfullmodel}; this model can be accessed using {it:estimates restore}. {p_end}

{phang}{opt gte:ffects}: Returns (and, for tables of up to 200 rows, prints) the individual cohort-by-relative-period effect estimates underlying {it:wooldid}'s aggregates
(the exact building blocks for the plain binary designs the identity notes below describe; under transformations, time-varying interactive controls, or {cmd:impvalues} the aggregates
are not simple reweightings of these rows), in {it:e(gteffects)} — one row per estimated {it:i}-by-relative-period cell, read directly off the underlying regression's coefficient vector at essentially zero additional
computation. Note that a coefficient dropped by the regression (collinearity omission) simply has no row: the table does not distinguish a partially identified horizon from a complete
one on its own - under default safety the dropped-coefficient screen errors first, but under {it:safety(off)} compare the row count against the expected cohort-by-period cells. Columns: {it:i} (the cohort/unit identifier value), {it:onset} (that unit's effective treatment onset, after any {cmd:anticipation()} shift), {it:caltime} (the calendar
period the cell covers), {it:relyr} (period relative to onset), {it:estimate}, {it:se}, {it:p}, {it:lb_95ci}, {it:ub_95ci} (both interval columns are named for the run's {cmd:level()}: {cmd:level(90)} gives {it:lb_90ci} and {it:ub_90ci}, and a fractional level writes the point as a p, {cmd:level(97.5)} giving {it:lb_97p5ci} and {it:ub_97p5ci}), {it:inwindow} (1 if the cell falls inside the reported
event-study window), {it:nobs} (estimation-sample observations in the cell), and the cell's total estimation weight under each requested estimand ({it:w_att}, {it:w_att_it},
{it:w_att_i}, {it:w_att_itime}; missing for estimands not requested). The weighted mean of {it:estimate} under any of these weight columns, taken over the rows that make up the aggregate in question, reproduces the
corresponding aggregate exactly: {it:relyr} >= 0 for a main row, {it:relyr} < 0 for a pretreatment row, and a single {it:relyr} for an event-study horizon. The post-treatment
filter is on {it:relyr}, not on {it:inwindow} - the main aggregate averages every estimable post-treatment period, so out-of-window rows belong in it - and the weight columns are
positive on pre-treatment rows as well, so a weighted mean taken over every row of the table reproduces no reported aggregate. With that filter the table supports arbitrary custom aggregations and is directly comparable to the group-time tables of {it:csdid} or {it:xthdidregress}
(the cells are ATT(g,t)-style group-time effects when {it:i} indexes cohorts; with a fine unit-level {it:i} they are unit-by-period effects and the table can be large). Notes:
standard errors are the regression's coefficient-level (conditional) VCE — {it:unconditionalse} and {it:impvalues} do not alter them; under {it:poisson} the estimates are on the
linear-index (xb) scale ({it:poisexpresults} and {it:semielasticity} are not applied); under {it:timetrends} the coefficient-level variances are typically not identified and the
{it:se} column is missing (degenerate zero variances are always reported as missing, never zero); with {cmd:subgroups()} a single full-sample table is returned; with time-varying
{cmd:interactivecontrols()} the interaction coefficients absorb control variation and the cell-level reading is unreliable. Not available with {cmd:contreat()}, where the
cohort-by-period interaction coefficients are dose-zero intercepts rather than treatment effects. Fully available under {it:fastnose} (the regression's own VCE is intact even
though the {it:margins} stage skips variances), making {it:gteffects} the one inference surface that survives fast mode. Under {cmd:bootstrap()} the cell standard errors, p-values and intervals are bootstrap quantities as well: each replication's cells are harvested by cohort and relative period, and the columns follow {cmd:bootdisplay()} and {cmd:boottstat} exactly as the main tables do (the standard error is the bootstrap standard deviation of the cell, or its analytic standard error under {cmd:boottstat}). The same holds for the {it:e(cfx*)} matrices behind the continuous-treatment plots, and under {cmd:bootstrap()} those plots are redrawn from the bootstrapped matrices after the rewiring (as {cmd:twoway} renderings of the matrices rather than {cmd:marginsplot} output, so they can differ cosmetically from an analytic run's plots, and {cmd:saveplots()} additionally writes one {it:name}{cmd:_cfx.dta} with every plotted point); only the {it:histogramcohorteffects} outputs remain analytic under {cmd:bootstrap()}. Under {cmd:permute()} the {it:gteffects} cells are keyed by cohort onset, which the replications reassign, so they have no permutation distribution and their inference columns are blank (the points stay); the {it:e(cfx*)} matrices are permuted for assignment-invariant lattice types and blanked otherwise, and the plots are redrawn from them as under {cmd:bootstrap()}. If no cohort-by-period coefficients can be read off the
fit (heavily collinear designs), no table is returned and a caution prints. Two reading notes. First, under {cmd:interactivecontrols()} the CELLS are approximate: the interaction coefficients absorb control variation, as the caveat above says, so a cell is not the clean cohort-by-period treatment effect it is elsewhere, and the weighted-mean identity that aggregates them inherits that approximation. (An earlier version of this passage asserted the identity remained exact there, contradicting the caveat fifteen lines above it; the assertion was wrong and has been withdrawn.) Recovering exact cells under that option would need an estimation-weighted within-cell mean per control and would turn each cell standard error from a diagonal read into a quadratic form; wooldid does not do this. Second, when the reported event-study window has no post-treatment horizons (e.g. {cmd:espostlength(0)}), {it:inwindow} is 0 on every post-treatment row by construction. {p_end}

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
{help maxvar} setting, {it:wooldid} halts with an informative error before estimation rather than failing mid-run with r(907). {it:wooldid} also cautions, regardless of the {it:safety()} setting, when {it:jointtests} is combined
with analytic clustered SEs clustered at the level of the cohort variable {it:i} or above, where the joint test p-values are severely oversized (see {it:jointtests}); the tests are
computed and printed either way, and the caution names the remedies ({it:permute()}, bootstrap-based clustered SEs, or clustering strictly finer than the cohort-by-period cells).
Under {cmd:bootstrap()} or {cmd:permute()} the screen does not run at all, since the joint tests are then bootstrap or permutation quantities rather than analytic ones.   {p_end}

{marker refusepolicy}{...}
{pstd}
{ul:The refuse / caution / skip policy.} {it:wooldid} has one vocabulary for what it does when a
specification is questionable, and it is worth stating once because the choice between the three is not
arbitrary. {bf:Refused} means the run stops with a named error: the combination would produce a number that
is wrong, or one whose meaning is undefined. {bf:Reported without inference} means the estimate stands and
the standard error, statistic, p-value and interval are blank: the point estimate is identified and the
uncertainty around it is not, and blanking is how that is said - never a zero. {bf:Cautioned} means
everything is computed and printed, with a line naming the risk and the remedy. No accepted option is a
silent no-op, and no check downgrades itself: {it:safety(off)} is the single documented way to move a
refusal to a caution or to skip a check, and it moves nothing else. It has never suppressed a disclosure
about what was estimated - if the run changed your specification, it says so whatever {it:safety()} is set
to. {p_end}

{marker safetysites}{...}
{pstd}
{ul:What {it:safety()} actually touches.} The option is not a single switch on a single check; it gates
nearly sixty sites across the shell, the worker, the validator and a dozen subprograms. This is the published
inventory, grouped by what each site DOES rather than listed one line at a time - the form that is usable, and
the form that stays true as the code moves. It is also the prerequisite for ever narrowing the option. {p_end}

{pstd}
{bf:First, three values are in play, and they are not the same.} What you typed; a {ul:working} value, which
{it:contreatwithin} together with {it:contreat()} silently rewrites to {it:off}, because that combination
makes the base treatment effects irrelevant by design; and a preserved copy of your own value. Roughly half
the sites gate on your own value and are therefore untouched by that rewrite - every data-integrity assert
(onset constancy within units, weight sign, the lattice and cfxlattice domains, onset-period
commensurability, the {it:customweightmultiplier} domain), the subgroup-containment rules, and both engines'
replication settings, which inherit {ul:your} value rather than the rewritten one. The rest gate on the
working value and are relaxed by that rewrite whether or not you asked for it. {p_end}

{pstd}
{bf:What default safety buys you, by what it does when it fires.} {p_end}

{synoptset 32 tabbed}{...}
{synopt:{bf:refuses, and {it:safety(off)} proceeds}} {it:unconditionalse} with a cell-weighted estimand; the cohort identifier not being a positive integer; negative weights; a non-integer or negative {it:lattice()} or {it:cfxlattice()} value; onset not commensurable with the time variable; and the variable-count preflight that would otherwise fail mid-run at r(907). {p_end}
{synopt:{bf:refuses, and {it:safety(off)} downgrades to a caution}} {it:contreatwithin} with the {it:byTime} continuous-treatment control type; {it:bootstrap(bayes)} with a cell-weighted estimand; and, under {cmd:permute()} with {cmd:subgroups()}, the two design-population screens - a subgroup variable that varies within a cohort-by-period cell, and one that splits a unit across subgroups. {p_end}
{synopt:{bf:skips a check entirely}} The whole post-estimation collinearity family, including the screen that reports an unidentified dose slope as {ul:missing} rather than as the regression's structural zero; the realized relative-time range asserts; and the control-availability count. {bf:Note} that this family gates on the {ul:working} value, so {it:contreatwithin} with {it:contreat()} skips it too. The visible consequence is on the dose-slope screen: where the continuous treatment carries no identifying variation, the {it:ame} rows print the regression's structural zeros instead of missing values, and without the caution that names them. This is the same output an explicit {it:safety(off)} produces, and that caution says so in as many words - a column of zeros with no standard errors beside them is not a set of estimates. {it:contreatwithin} declares that the base effects are irrelevant to the design, which is why it carries the relaxation. {p_end}
{synopt:{bf:suppresses a caution, changing nothing else}} The granular-{it:i} runtime warning; unconditional standard errors with aweights; periods with no in-sample control; the first-period and time-trend forcing notes; and the announcement that {it:espostlength} was raised. {p_end}
{synopt:{bf:changes what is computed}} Continuous-treatment effect plots for the non-ATT estimands are produced only under {it:safety(off)} - the default refuses to draw them and says so. The subgroup-contrast degeneracy site belongs here rather than among the caution suppressions: under default safety the flagged pairwise rows have their standard error, statistic, p-value and interval blanked in {it:e(pwsgs)} and the across-subgroup equality tests that cover them are posted as missing, while {it:safety(off)} reports all of them, with a caution of its own naming them as cancellation residue. A caution prints either way, in two different wordings, so what {it:safety(off)} changes here is the numbers rather than the volume. This site reads the value the user typed, not the working value {it:contreatwithin} rewrites. {p_end}

{pstd}
{bf:What {it:safety(off)} does {ul:not} buy.} About forty-five refusals stand whatever it is set to, and they
are the ones that matter most: every option-combination and domain screen, the missing-package check, the
{it:poisson} non-negativity requirement, every {it:fastnose} incompatibility, and the two cell-collapse
refusals whose own message says it in as many words:
{it:"No safety(off) escape: these are measured wrong answers, not risks."}
That is the line the option draws. {it:safety(off)} means {it:"I accept a risk"}; it has never meant
{it:"give me numbers I have been shown are wrong"}, and no future narrowing of the option
should move a check across that line. {p_end}

{pstd}
{bf:The engine-compatibility table.} One row per option the two engines treat differently, in one shape,
with the reason attached. {it:refused} means an error at parse; {it:allowed} means it runs and does what its
own entry says. {p_end}

{synoptset 30 tabbed}{...}
{synopt:{bf:option}}{bf:under} {cmd:bootstrap()} {bf:/ under} {cmd:permute()}{p_end}
{synopt:{cmd:impvalues}, {cmd:unconditionalse}, {cmd:oldsyntax}} refused / refused. Each replaces or reshapes the inference the engine is there to produce. {p_end}
{synopt:{cmd:customvce()}} refused / refused. The engine has to know what the replication's variance is; an arbitrary vce string is not inspectable. {p_end}
{synopt:multi-way {cmd:cluster()}} refused / refused. Both engines resample or relabel {ul:one} unit, so there must be exactly one. {p_end}
{synopt:{cmd:suppresssgprimaryeffects}} refused / refused. Under {cmd:permute()} the suppressed subgroup entries are precisely what the replications permute. {p_end}
{synopt:{cmd:fastnose}} refused / it {ul:is} the default. Fast mode is what makes the permutation affordable; {cmd:permanalyticse} opts back in. {p_end}
{synopt:{cmd:contreatwithin}, {cmd:contreatcontrols(untreatedZero)}} allowed / refused. Replications reassign onset while the dose path stays fixed, which contradicts the dose-zero-outside-treatment premise. {p_end}
{synopt:{cmd:pretrendpower()}} supported / refused. Roth's diagnostic consumes an event-study {ul:covariance}. The bootstrap has one and is fed it; randomization inference produces a null distribution of a statistic, not a covariance of the estimates. A difference of kind. {p_end}
{synopt:{cmd:simultaneousbands}} supported / refused. A confidence-band object; the permutation engine's multiplicity device is {cmd:permstepdown}. {p_end}
{synopt:{cmd:permstepdown}} not available / opt-in. The step-down stays permutation-only by ruling; the bands give the bootstrap the same family-wise device. {p_end}
{synopt:the studentizing switch} {cmd:boottstat} / {cmd:permtstat}. Both studentize, for different reasons - see fact 2 above. {p_end}
{synopt:{cmd:cluster()} nesting {it:i}} not required / required unless {cmd:permonsetcohorts}. The bootstrap re-identifies on the drawn cluster, so a cluster finer than {it:i} simply resamples inside the cohort. The permutation draws one onset per cluster and writes it to every row, so a unit spanning two clusters carries two onsets and is undefined. {cmd:permonsetcohorts} resolves that by rebuilding {it:i} from the drawn onsets. {p_end}
{synopt:{cmd:jointtests}} fixed-metric rank test, or the bootstrap of the statistic under {cmd:boottstat} / fixed-metric rank test. Neither is a Wald test on {it:e(V)}. {p_end}
{synopt:{cmd:adversarial}} bootstrapped only under {cmd:boottstat}; otherwise absent with a caution and {it:e(boot_adv_absent)} set, never analytic under a bootstrap header / a fixed-metric rank test. {p_end}
{synopt:{cmd:histogramcohorteffects}} analytic under both, with a caution: it is computed from the stage-1 covariance before the replications run. {p_end}
{synopt:draw count} {cmd:reps()} required / {cmd:permreps()} required, and {cmd:permexact} enumerates instead. {p_end}



{phang}{opt oldsyntax}: In older versions, {it:wooldid} used the suffix {it:contemp} to flag the event study coefficient containing the period contemporaneous with treatment. 
This turned out to make {it:wooldid} relatively incompatible with the popular event study plot production program {help event_plot}, so the {it:contemp} suffix has been 
replaced with the suffix {it:post0}. This option reverts {it:wooldid} to using the older coefficient syntax, for compatibility with any older code that relies on it. It cannot be combined with {cmd:permute()} or {cmd:bootstrap()}: each refuses {it:oldsyntax} with an error. {p_end}




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
The {cmd:bootstrap()} re-identification leg - required by {it:att_it}, {it:att_i} and {it:att_itime}, and taken by {cmd:boottstat} whenever {cmd:cluster()} is given - roughly doubles per-replication cost where the cluster nests {it:i}, because the re-identified copies present a larger cell count to the regression; where the cluster is finer than {it:i} or cross-cuts it the units are kept whole, the cell count does not grow, and the overhead is the physical frame alone. On the other side,
{it:fastnose} skips the variance calculation that accounts for most of a {it:margins} call and is correspondingly worth several times over on runs where
standard errors are not needed. The default cell collapse (see {it:fastcollapse}) attacks the same cost from
the other side - it shrinks the data each {it:margins} call sees rather than skipping the variance - and the two compose. The collapse is already on unless
{it:fullprecisionslow} or an ineligible option turns it off, so a run that must be evaluated on the full data ({it:unconditionalse}, {it:interactivecontrols()}, {it:customvce()},
{it:poisexpresults}, {it:semielasticity}, {it:contreatpoly()} above 1, or one of the other {cmd:contreat()} reads on the ineligible list under {it:fastcollapse}) is the one paying
the full per-row {it:margins} cost described here. {cmd:subgroups()} is {ul:not} one of them: it collapses like any other run, and only sets the collapse aside part-way through
if the estimability screen described under {it:fastcollapse} fires. {p_end}

{pstd}{bf:Where the time goes.} Nearly all of a typical run is spent inside Stata's {it:margins}, and most of each {it:margins} call is the delta-method
variance calculation rather than the point estimates. Very little time is spent in {it:wooldid}'s own code, which is why the levers that matter are the ones
that change the size of the problem handed to {it:margins} - K above all - rather than anything about how the data are stored or passed in. {p_end}

{marker output}{...}
{marker storedresults}{...}
{title:Working With Output}

{pstd}
{it:wooldid} generates output primarily in the form of matrices and scalars stored in e(), though all main results can also be accessed using standard utilities like {it:esttab}.
Wherever a stored name is composed from a subgroup value or from a continuous-treatment plot name, a hidden twin or companion whose name would run past Stata's 32-character limit is
skipped silently; the visible table is never affected. That covers the per-family bootstrap subgroup tables {it:e(boot_*_sgresults_*)} and the studentized {it:e(boot_tpctl_sgresults_*)};
the per-family bootstrap cfx twins {it:e(boot_*_}{it:cfxname}{it:)}, which are lost once the cfx name passes 22 characters, reachable under {it:makecfxplotsbysg} without any long subgroup
value at all; the per-subgroup {cmd:permstepdown} tables; and the {cmd:simultaneousbands} companions, where the two lengths behave differently: a very
long family key drops the band and its companions together, so the family simply never appears in {it:e(esband_families)}, while an intermediate one posts {it:e(esband_{it:type})} without
{it:e(esband_crit_{it:type})} and {it:e(esband_n_{it:type})}.
The commands {it:est restore} and {it:est save} are compatible with {it:wooldid} for saving its output. The visible {it:wooldid} output in e() can be viewed using the command
{it:ereturn list}; a number of the results described below are posted {it:hidden} - the preserved analytic covariances and degrees of freedom, the per-family bootstrap and permutation tables, the permutation kernel, the tail and comparator records, and {it:e(wdid_vcvstate)} and its companions - and those appear only under {it:ereturn list, all}, though each can be read by name ({it:di e(wdid_vcvstate)}) either way. {p_end}

{pstd} The coefficient vector {it:e(b)} and its covariance {it:e(V)} carry these column names, which are what
{it:esttab}, {it:coefplot}, {cmd:lincom}, {cmd:test} and the {cmd:permtstat}/{cmd:bootsave()} tables address a coefficient by:
{it:main_att} (the headline average treatment effect on the treated), {it:pretreatment_att} (the pooled pre-treatment
average), {it:ES_att_pre#} and {it:ES_att_post#} (the event-study horizons, numbered outward from the reference
period), and for the other estimand types the same names with {it:att} replaced by {it:att_it}, {it:att_i}, {it:att_itime}
or {it:ame}. Under {cmd:subgroups()} each subgroup's rows take the suffix {it:_s#}. {p_end}

{pstd} The primary results are stored in the following matrices, scalars, and macros: {p_end}
{synoptset 30 tabbed}{...}
{synopt:{it:e(wdidmainresults)}} All primary and event study results. The columns are {it:estimate}, {it:se}, {it:t} ({it:z} under {it:poisson}), {it:p}, {it:lb_95ci}, {it:ub_95ci} and {it:relyr}, the two interval columns being named for the run's {cmd:level()}: {cmd:level(90)} gives {it:lb_90ci} and {it:ub_90ci}, and a fractional level writes its point as a {it:p}, so {cmd:level(97.5)} gives {it:lb_97p5ci}, the same rule {it:e(gteffects)} follows. Rows are equation-and-name pairs: the pooled estimands sit under the equations {it:main} and {it:pretreatment} ({it:main:att}, {it:pretreatment:att}, {it:main:ame}, and so on) with {it:relyr} missing, and the event-study rows sit under {it:ES_att}, {it:ES_att_it}, {it:ES_ame} and {it:ES_ame_it}, named {it:post#} and {it:pre#} and carrying their horizon in {it:relyr}.{p_end}

{synopt:{it:e(wdidsgresults_{it:sgcase})}} Primary and event study results for the subgroup {it:sgcase}, if subgroups are used and presentation of their
primary results has not been suppressed.{p_end}

{synopt:{it:e(pwsgs)}} Pairwise differences between the subgroup-specific treatment effects, if subgroups are used. The columns are {it:estimate}, {it:se}, {it:t} ({it:z} under {it:poisson}), {it:p}, {it:lb_95ci} and {it:ub_95ci}, named for the run's {cmd:level()} as elsewhere; there is no {it:relyr} column here, unlike {it:e(wdidmainresults)}. A row whose analytic contrast covariance is rank-deficient carries its point estimate with every inference column missing - see {help wooldid##sgdegen:the note on subgroup contrasts}. {p_end}

{synopt:{it:e(b_es_{it:type})}, {it:e(V_es_{it:type})}} The event study coefficient vector and its {ul:full} variance-covariance matrix (carrying genuine off-diagonal covariances; since 2026-09-07 the same blocks are spliced into {it:e(V)} - see {it:e(wdid_vcvstate)}), relative-period-labeled as {it:pre#}/{it:post#} with the reference period omitted. Returned per estimand type: {it:att} (binary) and, with {cmd:contreat()}, {it:ame}, plus {it:att_it} / {it:ame_it} for the {it:att_it}, {it:att_i}, and {it:att_itime} estimands (which share a single event-study horizon internally). With {cmd:subgroups()}, per-subgroup versions {it:e(b_es_{it:type}_s{it:sgcase})} / {it:e(V_es_{it:type}_s{it:sgcase})} are returned as well. After {it:ereturn post} these matrices are valid inputs for {it:lincom} and {it:test}. A {ul:direct positional handoff} to {it:honestdid} (Rambachan-Roth robustness) additionally requires ALL of:
{cmd:esfixedbaseperiod} with the default {cmd:esrelativeto(-1)}; no {cmd:timetrends}; at least one pre- and one post-treatment coefficient; horizons consecutive apart from the intended
jump over the omitted -1 reference (checkable signature on the colnames: {it:pre1} absent, rightmost pre column is {it:pre2}, {it:post0} present, no other gaps); and a covariance with no
structurally degenerate rows (see {it:e(es_ndegen_{it:type})} below; honestdid inference inherits whatever VCE you specified). {it:numpre} is purely positional and cannot communicate
gaps, non-(-1) references, or the (default) pooled cohort-specific reference — those geometries hand honestdid the wrong event-time geometry {ul:silently} and need a user-supplied
transformation. Check the regime via {it:e(estimate_type)} == "fixed_reference_period" and the numeric spacing via the {it:relyr} column of {it:e(escomp)}/{it:e(escomp_s{it:sgcase})}.
Where the dose-slope screen voids the ame block entirely (every treated cohort-period dose-slope term dropped for multicollinearity, as when the dose is constant within cohort), both {it:e(b_es_ame*)} and {it:e(V_es_ame*)} are withheld, because the captured ame points are the regression's structural zeros rather than estimates; the Caution that normally announces a coefficient vector posted without its covariance is deliberately suppressed on those families, so nothing in the output reports the absence either. A {it:honestdid} or {it:pretrends} handoff therefore finds no ame input at all on such a design, the att-side event-study objects being unaffected and the structural zeros surviving only in the raw {it:e(b)}. Under cell containment {it:e(b_es_ame*)} does still post. Otherwise the matrices are returned whenever at least one horizon has an identified variance; horizons whose standard errors are not identified then carry zero (structural-zero) or missing
diagonal rows, flagged by a Caution and by {it:e(es_ndegen_{it:type})} (the count of degenerate rows, posted only when such rows exist — its absence marks a clean covariance). When NO horizon has an identified variance ({cmd:fastnose}, typically {cmd:timetrends}, and any run whose estimation returns an all-zero covariance), the covariance is not returned but the coefficient vector still is, under a Caution naming what has no input on that run: the estimates remain meaningful even where their uncertainty is unavailable. Additionally, when the cell-containment caution fires (see the continuous-treatments
standard-error caution above), {it:e(V_es_ame*)} is withheld — the captured covariance is not the sampling covariance of the ame estimates — while {it:e(b_es_ame*)} still posts; att-side
event-study objects are unaffected. Separately, {it:contreatcontrols(byCohort)} requires each treated cohort's untreated/reference-period observations to carry variation in the
continuous treatment (that variation is all that identifies the cohort-specific baseline dose slope): wooldid refuses with r(459), naming the offending cohort, when a treated cohort has
none, and prints a computational-fragility caution (suggesting a larger {it:regtol()}) when that variation is below 1/1000 of the sample standard deviation of the treatment variable.{p_end}

{synopt:{bf:Design conditions for every consumer of} {it:e(V_es_*)}} The event-study covariance is a delta-method object carried out of the {cmd:margins} call, so it inherits the design's {ul:identification}, not only its arithmetic - and neither condition below is visible in the matrix itself.
{bf:One:} the treated cells must not be absorbed by their own terms. Where a treated ({it:i}, relative-period) cell lies wholly inside one cluster of {cmd:cluster()}, the cluster-robust covariance cannot see that cell's own variation; {it:wooldid} cautions on exactly that predicate (see {cmd:jointtests}), and the same predicate governs how much the matrix says about the treated horizons at all.
{bf:Two - singleton-cell designs:} with a unit-level {it:i}, each treated cohort-by-period cell is a single observation, the treated residuals vanish, and the covariance carries control-cluster noise only. Measured on that geometry, the event-study standard errors sit at 0.84 to 1.07 of the pooled event study's while the point estimates agree to 1e-15 - the numbers look entirely ordinary and the object is not measuring what it appears to.
{it:honestdid} and {cmd:pretrendpower()} consume this matrix directly and inherit both conditions; a coarser {it:i} (cohorts, or treatment-onset groups) with unit fixed effects supplied through {cmd:fe()} is the design that satisfies them. {p_end}

{synopt:{it:e(numpre_{it:type})}} The number of pre-treatment relative periods (relative year < 0) in the corresponding event-study vector, provided for use as {it:honestdid}'s {it:numpre()} argument (subject to the geometry conditions in the {it:e(b_es)} entry above).{p_end}

{synopt:{it:e(gteffects)}} The cohort-by-relative-period building-block table produced by the {it:gteffects} option; see that option's description for the column layout.{p_end}

{synopt:{it:Bootstrap inference}} When {cmd:bootstrap()} is used: scalars {it:e(boot_reps)}, {it:e(boot_completed)}, {it:e(boot_failed)}, {it:e(boot_level)}, {it:e(boot_N_V)} (replications behind the covariance) and {it:e(boot_v_posted)}; macros {it:e(boot_method)}, {it:e(boot_stat)} ({it:z} or {it:t} under {cmd:bootdisplay(wald)}, {it:z_wald} or {it:t_wald} under the percentile and bias-corrected displays, {it:t_pt} under {cmd:boottstat}), {it:e(boot_seed)} (under {cmd:rseed()}), {it:e(boot_rngstate)} (the random-number state before the first draw), {it:e(boot_alpha)} (under {cmd:bootstrap(bayes)} only) and {it:e(boot_ties)} (under {cmd:ties} only); matrices {it:e(V_boot)}, {it:e(b_boot)} (replication means), {it:e(bias_boot)}, {it:e(boot_navail)} (available replications per estimate), {it:e(boot_ci_normal)}, {it:e(boot_ci_percentile)}, {it:e(boot_ci_bc)}, {it:e(boot_z0)}, under {cmd:bootdisplay(bca)} also {it:e(boot_ci_bca)}, {it:e(boot_p_bca)}, {it:e(boot_accel)} and the scalar {it:e(boot_jk_completed)}, and - up to 10,000 completed replications - {it:e(bootdist)}, the replication matrix itself. The analytic covariances the bootstrap replaced are kept as hidden {it:e(V_analytic)} and {it:e(V_es_analytic_*)}. {it:e(boot_ves_failed)} lists the event-study families ({it:att}, {it:att_it}, {it:ame}, {it:ame_it} and their {it:_s#} subgroup twins) whose {it:e(V_es_*)} could not be rebuilt from the replications and therefore remains the analytic matrix under a bootstrap header; the equivalent test is the absence of that family's hidden {it:e(V_es_analytic_*)}, which is created only on the success branch. The consumers of that matrix, {cmd:pretrendpower()} in particular and the {it:honestdid} handoff with it, are skipped for the listed families rather than run against an analytic covariance. The macro is what survives into a replayed or restored estimate once the run-time caution has scrolled away; its absence means every family was rebuilt. Three inference families are stored on every bootstrap run whichever one is displayed, and a fourth, BCa, whenever {cmd:bootdisplay(bca)} is requested, the acceleration being computed only on request; the BCa tables are {it:e(boot_bca_mainresults)}, with {it:e(boot_bca_sgresults_{it:sgcase})} and {it:e(boot_bca_pwsgs)} alongside, the {it:_sgresults_*} and {it:_pwsgs} counterparts existing for all four families. Of the three standing ones: {it:e(boot_p_percentile)} and {it:e(boot_p_bc)} hold the percentile and bias-corrected p-values, and the hidden {it:e(boot_pctl_mainresults)} (also posted under the derivable name {it:e(boot_percentile_mainresults)}, so that {it:e(boot_}{cmd:`e(boot_display)'}{it:_mainresults)} resolves for every display including the default), {it:e(boot_bc_mainresults)} and {it:e(boot_wald_mainresults)} - with their {it:_sgresults_*} and {it:_pwsgs} counterparts - hold the complete tables under each. The reported standard error is the bootstrap standard deviation under all three, so it does not vary with {cmd:bootdisplay()}; under {cmd:boottstat} it is the analytic standard error instead, the displayed family is the studentized one (the hidden {it:e(boot_tpctl_*)} tables), and the coefficient-draw Wald, percentile and bias-corrected tables remain stored for reference but are not displayable. A {cmd:boottstat} run also stores {it:e(boot_tstud)} (rows {it:n}, {it:t_obs}, {it:p_pctl}, {it:qlo_pctl}, {it:qhi_pctl}, {it:live}, one named column per coefficient), {it:e(bootdist_se)} (the replication standard errors, on the same row mask as {it:e(bootdist)}), the hidden {it:e(boot_ase_obs)} (the observed analytic standard errors by name), {it:e(boot_tstud_pwsgs)} and {it:e(bootdist_pwse)} for the contrasts, {it:e(boot_jointW)} (the replication joint statistics, one named column per test) and {it:e(boot_advW)} (the replication adversarial statistics), so every reported studentized number is reconstructable from {it:e()} by name. The {it:gteffects} and {it:e(cfx*)} cell tables are stored per family as the hidden {it:e(boot_pctl_gteffects)}, {it:e(boot_bc_gteffects)}, {it:e(boot_wald_gteffects)} and {it:e(boot_pctl_}{it:cfxname}{it:)} etc., with {it:e(boot_gteffects_family)} recording the displayed one; the adversarial p-value's absence without {cmd:boottstat} is recorded in {it:e(boot_adv_absent)}. The macro {it:e(boot_display)} records which family was displayed. The bootstrap analogue of the permutation tail record is the hidden {it:e(boot_ptails)} (rows {it:p_lo}, {it:p_hi}, {it:se_p}, {it:cp_lo}, {it:cp_hi}: the one-sided percentile shares about 0, the Monte-Carlo standard error of the percentile p-value and its {cmd:level()} Clopper-Pearson interval, one named column per coefficient), and the rep-0 comparator vector is the hidden {it:e(boot_comparator)} - the counterparts of {it:e(perm_p_lo)}/{it:e(perm_p_hi)}, {it:e(perm_p_se)}/{it:e(perm_p_ci)} and {it:e(perm_comparator_b)}. {it:e(boot_jt_family)} records how the joint and equality tests were formed ({it:studentized} under {cmd:boottstat}, {it:fixedmetric} otherwise), with the replication statistics in {it:e(boot_jointW)} / {it:e(boot_sgW)} and the observed ones in {it:e(boot_jointWobs)} / {it:e(boot_sgWobs)} (the latter two without {cmd:boottstat}); when some joint or equality test could not be computed, {it:e(boot_jt_absent)} or {it:e(boot_sgjt_absent)} is set, so the disclosure survives replay. Those two flags say {ul:that} a test was withheld; {it:e(boot_sgjt_reason)} says {ul:why}, carrying one or more of {it:nostatistic} (attempted on the fixed metric and no statistic survived - fewer than two complete replications carry every subgroup level, or the bootstrap covariance of those levels is rank-deficient), {it:nostudentizedinputs} (under {cmd:boottstat}, the observed or per-replication pairwise standard errors the studentization needs were unavailable) and {it:unmapped} (the tested columns could not be matched to replication columns). The run says the same thing in a caution at the time; the macro is what survives into a replayed estimate. Absence of the flag means nothing was withheld, not that nothing was attempted. {it:e(boot_bc_absent)} is the same kind of record for the bias-corrected family, one level down: it lists the COEFFICIENTS whose {it:z0} was undefined because every completed draw fell on one side of the rep-0 comparator, and whose {it:e(boot_z0)}, {it:e(boot_ci_bc)}, {it:e(boot_p_bc)} and BCa twin are missing for that reason while their percentile and Wald entries stand; its absence means the family was withheld nowhere. The scalar {it:e(boot_df)} records the degrees of freedom in force whenever the reference distribution is a t - under {cmd:bootdf()} and equally under its default, which is {cmd:bootdf(clusters)} when {cmd:cluster()} was given and {cmd:bootdf(normal)}, which posts no {it:e(boot_df)} at all, otherwise. Under {cmd:boottstat} the scalar {it:e(boot_tstat)} is 1 and {it:e(boot_stat)} reads {it:t_pt}; the studentized kernel is {it:e(boot_tstud)} (described above), {it:e(bootdist_se)} holds the replications' own standard errors beside {it:e(bootdist)}, and the coefficient-draw tables stay stored but undisplayed. {it:e(df_r)} is aligned with the reference distribution actually reported, so a {cmd:test} or {cmd:lincom} run by hand after a bootstrap fit reproduces the {cmd:bootdisplay(wald)} p-values under every {cmd:bootdf()} setting - and note that {it:e(V)} and {it:e(df_r)} carry the standard-deviation family whichever family is displayed, so under the default percentile display they reproduce the stored {it:e(boot_wald_mainresults)} numbers rather than the columns printed in the tables: under {cmd:bootdf(clusters)} or an explicit integer it holds those degrees of freedom, and under {cmd:bootdf(normal)} it is set to missing so that {cmd:test} and {cmd:lincom} fall back to chi-squared and z, matching the reported statistics. The analytic residual degrees of freedom the bootstrap replaced are kept as hidden {it:e(df_r_analytic)}, alongside the analytic covariances. {p_end}

{synopt:{it:Permutation inference}} When {cmd:permute()} is used, the reported tables themselves carry the permutation inference (columns {it:se_analytic}, {it:t_perm}, {it:p_perm} and the
{cmd:level()} interval), and the hidden {it:e(perm_analytic_mainresults)}, {it:e(perm_analytic_sgresults_*)}, {it:e(perm_analytic_pwsgs)} and {it:e(perm_analytic_}{it:cfxname}{it:)} keep the tables as the
observed run posted them (a hidden twin whose name would exceed Stata's 32-character limit - a subgroup value of 9 or more digits, or a long cfx name - is skipped; the visible table is unaffected. For the cell-draw twins each pair is guarded separately: a cfx name of 22-23 characters keeps {it:e(permdist_*)}/{it:e(permcmp_*)} and loses only the {cmd:permtstat} se twins). Stripe-aligned surfaces, one named column per coefficient of {it:e(b)}: {it:e(perm_p)} (the equal-tailed permutation p, or the absolute-value p under {cmd:permabsp}), {it:e(perm_ci)} (rows {it:ll}/{it:ul}, the shift interval),
{it:e(perm_se_shift)} (the interval-implied standard error), {it:e(perm_p_lo)} and {it:e(perm_p_hi)} (the one-sided shares), {it:e(perm_stepdown)} and its family/pre/subgroup variants (the max-|studentized| step-down over a family's horizons, under {cmd:permstepdown}), {it:e(perm_p_se)} and {it:e(perm_p_ci)} (the Monte-Carlo standard error of the p-value's share and its {cmd:level()} Clopper-Pearson interval), {it:e(perm_sd)} (SD of the permutation distribution), {it:e(perm_n)} (replications carrying the element), {it:e(perm_nullband)} (rows {it:lo}/{it:hi}: the
equal-tailed {cmd:level()} quantiles of the null draws),
{it:e(perm_comparator_b)} (the rep-0 comparator vector), the hidden {it:e(perm_kernel)} (every row of the per-element kernel, by name) and {it:e(perm_ase_obs)} (the observed analytic standard errors), and -
up to 10,000 completed replications - the draw matrices {it:e(permdist_b)} (and {it:e(permdist_se)} under {cmd:permtstat}). Test statistics: {it:e(perm_jointW)}, {it:e(perm_sgW)} and {it:e(perm_advW)} (the
replication statistics behind the joint, equality and adversarial p-values, one named column per test), their observed twins {it:e(perm_jointWobs)}, {it:e(perm_sgWobs)} and {it:e(perm_advWobs)} (the observed statistic each p ranks; the p is reconstructable from the pair), the scalar {it:e(perm_adv)} (= 1 when the adversarial line tests were formed from the permutations) and the matrix {it:e(perm_pwsgs_kernel)} (the pairwise subgroup contrasts' kernel); {it:e(perm_comparator_se)} (the rep-0 analytic standard errors under {cmd:permtstat}); the hidden cell-draw twins {it:e(permdist_}{it:cfxname}{it:)}, {it:e(permcmp_}{it:cfxname}{it:)} and, under {cmd:permtstat}, {it:e(permdistse_}{it:cfxname}{it:)}/{it:e(permcmpse_}{it:cfxname}{it:)} behind the continuous-treatment plot cells; the macro {it:e(perm_offcentre)} (elements whose shift interval was blanked); the tests themselves post under the analytic names
({it:e(joint_horiz)} etc., {it:e(joint_*_allsgssame)}, {it:e(adversarial_pval_*)}). Run-level: {it:e(perm_reps)}, {it:e(perm_completed)}, {it:e(perm_failed)}, {it:e(perm_npossible)} (also as the string
{it:e(perm_npossible_str)}), {it:e(perm_level)}, {it:e(perm_exact)}, {it:e(perm_blank)} (structural-zero ame cells excluded), {it:e(perm_precaveat)} (1 when the permuted coefficient set contains any pre-treatment row (a pretreatment-equation estimate, or any horizon at a negative relative year) and 0 otherwise; when it is 1 the staggered-design sharp-null caveat described under {cmd:permute()} in {help wooldid##standarderrors:Standard Errors} applies, and this scalar is what replays that printed caution into a restored estimate), {it:e(perm_gt_blank)}, {it:e(perm_cfx_blank)},
{it:e(perm_analyticse)}, {it:e(perm_absp)}, {it:e(perm_irebuilt)} (1 when the replications rebuilt {it:i} from the onsets they drew under {cmd:permonsetcohorts}, 0 when that option was given with {cmd:cluster()} already at {it:i}, where the rebuild is a no-op; not posted without the option); macros {it:e(perm_method)}, {it:e(perm_stat)} ({it:b} or {it:t}), {it:e(perm_pconv)} ({it:equaltailed}, the default; {it:absolute} under {cmd:permabsp}; {it:raw} under {cmd:permrawp}; {it:left} or {it:right} under {cmd:permside()}, the same five conventions mapped at {help wooldid##pvalues:Which p-value is which}, of which at most one of the three options may be given), {it:e(perm_seed)} (under {cmd:rseed()}), {it:e(perm_rngstate)} (the random-number state before the first draw) and {it:e(perm_incomplete)}, {it:e(perm_addone)} (the elements whose p-values took the add-one form). The legacy pooled-estimand scalars {it:e(perm_p_{it:est})} (now the same p
the table carries, in the convention in force), {it:e(perm_bandlo_{it:est})}/{it:e(perm_bandhi_{it:est})} (the {cmd:level()} null band - the region the statistic occupies under the null, NOT an interval for the estimand; the older names {it:e(perm_lb_{it:est})}/{it:e(perm_ub_{it:est})} hold the same two numbers and are deprecated, since lb/ub elsewhere in wooldid means a genuine interval bound), {it:e(perm_sd_{it:est})} and, under {cmd:permtci}, {it:e(perm_tlb_{it:est})}/{it:e(perm_tub_{it:est})}, the
estimand-keyed comparator {it:e(perm_comparator)} and the estimand-keyed draw matrix {it:e(permdist)} (on the {it:e(perm_stat)} scale) are kept. Note that {it:e(p)}, {it:e(ci_l)} and {it:e(ci_u)} are also posted under the {cmd:permute()} fast-mode default, under the names {it:estout} looks for ahead of anything derived from a covariance; see the {cmd:permute()} entry for what they hold and why no {it:e(se)} accompanies them. {it:e(permdist)} is {ul:not} the counterpart of {it:e(bootdist)} despite the name: it carries one column per requested estimand on the statistic scale, where {it:e(bootdist)} is the full coefficient-draw matrix aligned to {it:e(b)}. The counterpart is {it:e(permdist_b)}, and {it:e(permdist)} is deprecated in favour of it.{p_end}

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
pre-trend slope and its {cmd:level()} confidence interval (the run's own level, {cmd:c(level)} by default and so normally 95), built from a {it:t} critical value where residual degrees of freedom exist and a normal one otherwise), {it:e(roth_prep_{it:type})} (the joint pre-trend test p-value), {it:e(roth_hbar_{it:type})} (the estimand-weighted mean post-treatment relative period over
the full post support of the estimation sample), {it:e(roth_npre_{it:type})}, {it:e(roth_keff_{it:type})} (pre-period count and effective rank), {it:e(roth_W_{it:type})} (the joint
pre-trend test statistic), and {it:e(roth_nissize_{it:type})} (the actual size of the any-individually-significant pre-test). {it:e(roth_ok)} = 1 marks a completed diagnostics
pass: a type absent from the roth_* results was then legitimately skipped or unavailable, and roth_ok = 1 with no {it:e(roth_power_*)} matrix at all means the diagnostics were
requested but none were computable (that caution is replayed). A roth_* object present {ul:without} {it:e(roth_ok)} marks an interrupted pass; replay then prints a caution and
does not render the partial block.{p_end}

{pstd} In addition to these core results, the following materials are also returned in e(): {p_end}
{synopt: Joint Tests by Horizon} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented in the scalars e(joint_horiz)
(ATT), e(joint_horiz_i) (the shared ATT_it / ATT_i / ATT_itime horizon block), e(joint_Choriz) (AME; continuous treatments only), and e(joint_Choriz_i) (AME_it; continuous treatments only). What the p-value is depends on the inference engine: on an analytic run it is the p of {it:margins}' {cmd:contrast(overjoint)} test over the pre-treatment horizons; under {cmd:bootstrap()} it is an add-one rank p of the replication statistics, of the family named by {it:e(boot_jt_family)} ({it:studentized} or {it:fixedmetric}), the statistics themselves in {it:e(boot_jointW)}/{it:e(boot_jointWobs)}; under {cmd:permute()} it is a rank p of the permutation statistics in {it:e(perm_jointW)}/{it:e(perm_jointWobs)}, in the add-one form unless the observed assignment is itself one of the draws or the enumeration is exact. The three are not comparable to one another. {p_end}

{synopt: Joint Tests by Subgroup & Time Period} Joint tests for the statistical significance of pre-treatment relative time period-specific event study effects are presented by
subgroup, when used, in e(joint_horiz_s{it:sgcase}), e(joint_horiz_i_s{it:sgcase}), e(joint_Choriz_s{it:sgcase}), and e(joint_Choriz_i_s{it:sgcase}), where {it:sgcase}
is the value of the subgroup. These are e() {it:scalars} on an analytic run and e() {it:local macros} under {cmd:bootstrap()} and {cmd:permute()}; an {it:e(}{it:name}{it:)} reference reads either, but code that keys on the type - {cmd:confirm scalar}, or filtering {it:ereturn list} - must expect both. The engine-dependent p-value conventions described just above apply here too. {p_end}

{synopt: Adversarial Line Tests:} P-values from the composite-null tests that the relative-time period treatment effect estimates and the reference period fall on the same line
 are found in the scalars {it:e(adversarial_pval_att)}, {it:e(adversarial_pval_ame)}, {it:e(adversarial_pval_att_s[SG])}, etc.; they are posted only under {it:esfixedbaseperiod}
 and only when the event-study covariance is not degenerate (a caution is printed otherwise). The display-oriented fitted line is always returned in matrices named
 {it:e(adversarial_line_att)} and {it:e(adversarial_line_att_it)}, joined under {cmd:contreat()} by {it:e(adversarial_line_ame)} and {it:e(adversarial_line_ame_it)}, with per-subgroup twins named {it:e(adversarial_line_{it:type}_s{it:sgcase})} exactly as the p-value scalars are, and with columns {it:slope} and {it:intercept_rt0}; the line matrix and the p-value scalar are gated independently, so the line can be present where the p-value is not; {it:intercept_rt0} is the fitted value at numeric relative time 0, regardless of wooldid's
 internal nonnegative time coding. {p_end}

{synopt: SG Effect Homogeneity Tests} For each estimate type (ATT, ATT_it, AME, etc.), the joint test that all subgroups have the same such post-treatment effect (namely, the
same effect as the lowest numbered SG) may be found in the scalar {it:e(joint_att_allsgssame)} and analogously named scalars for other treatment effect types. For pre-treatment effects,
the test can be found in objects with names like {it:e(joint_pre_att_allsgssame)}. A test whose pairwise covariance is rank-deficient is stored as missing rather than computed from that covariance - see {help wooldid##sgdegen:the note on subgroup contrasts}.{p_end}

{synopt:{it:e(esband_{it:type})}, {it:e(esband_crit_{it:type})}, {it:e(esband_engine)}, {it:e(esband_level)}, {it:e(esband_families)}, {it:e(esband_ndraw)}, {it:e(esband_n_{it:type})}} Under {cmd:simultaneousbands}: the simultaneous band endpoints per event-study family (rows {it:lb}, {it:ub}, {it:relyr}; one column per horizon, {it:_s#} per subgroup), the common critical value, which engine produced them, the level, the families that got one, and, per engine, the simulation draw count or the listwise-complete replication count ({it:e(esband_ndraw)} exists on the simulated engine only). {it:e(esband_engine)} takes three values: {it:analytic-simulated} when {cmd:bootstrap()} was not used, drawing from the family's own analytic correlation matrix; {it:bootstrap} under {cmd:bootstrap()} without {cmd:boottstat}, where the deviations are scaled by the coefficient-draw standard deviation; and {it:bootstrap-studentized} under {cmd:bootstrap()} with {cmd:boottstat}, where each replication's own standard error is the denominator and the observed analytic ones set the endpoints. Code keying on the exact string must accept all three, since {it:bootstrap-studentized} is still a bootstrap band. Note also that the column count of {it:e(esband_{it:type})} is the number of horizons actually banded, which on the simulated engine can be fewer than the family has; see {opt simult:aneousbands}. Absent when the option was not given. {p_end}
{synopt:{it:e(wdid_vblocks)}, {it:e(wdid_vblockn)}} The event-study families whose genuine covariance blocks were spliced into {it:e(V)}, and the number of off-diagonal pairs written. {bf:What this means for} {it:e(V)}: it is a {ul:mixed} matrix. Within one event-study family - the relative-period rows of {it:att}, {it:att_it}, {it:ame} or {it:ame_it}, per subgroup - the off-diagonals are the real joint covariance of those estimates, because a family's horizons all come from a single {cmd:margins} call whose {it:r(V)} is exactly that covariance, so a {cmd:test} or {cmd:lincom} {ul:across the horizons of one family} is valid. Everywhere else the off-diagonals are still zero and still mean "not estimated", not "zero covariance": a test spanning two families, or a family and a pooled estimate, is not valid on this matrix and will overstate its precision. The general cross-family matrix is deliberately not built - three modes ({cmd:impvalues}, {cmd:unconditionalse}, {cmd:fastnose}) have no common delta-method object, and a pooled effect is a linear combination of its own horizons, so the full matrix would be rank-deficient by construction. A family is spliced only when every horizon carrying a variance in {it:e(V)} matches a row of {it:e(V_es_*)} by name and agrees with it numerically; rows whose variance was never identified keep a zero row and column.
{ul:Under} {cmd:bootstrap()} these two results describe the {ul:preserved analytic} matrix, not the reported one: the splice runs before the replications are wired in, so the blocks are in the hidden {it:e(V_analytic)} while {it:e(V)} holds the full replication covariance ({it:e(wdid_vcvstate)} reads {it:bootstrap} there and is the field to check). {p_end}
{synopt:{it:e(wdid_vcvstate)}} What {it:e(V)} is on this run, in one word, so downstream code can check rather than infer: {it:analytic-diagonal} (variances only, every off-diagonal zero - a {cmd:test} or {cmd:lincom} spanning two estimates is NOT valid on it), {it:analytic-blocks} (as diagonal, plus the genuine within-family event-study blocks described under {it:e(wdid_vblocks)} - the default on any analytic run that produced an event study), {it:bootstrap} (the full replication covariance, with {it:e(df_r)} aligned so a hand-run {cmd:test} reproduces the {cmd:bootdisplay(wald)} table), or {it:absent} (no {it:e(V)} posted at all - the {cmd:permute()} fast-mode default, every {cmd:fastnose} run, and any run whose variance could not be posted, which the run says at the time). Under {cmd:permute()} with {cmd:permanalyticse} or {cmd:permtstat} the posted matrix is the ANALYTIC one, shown for reference: it is not a permutation quantity. {p_end}
{synopt:{it:e(N)}, {it:e(N_dropped)}, {it:e(Nclus)}, {it:e(N_clust)}, {it:e(clustvar)}, {it:e(vce)}, {it:e(vcetype)}} Estimation sample observation count, number of observations dropped from initially specified sample (including rows removed because the {it:if} expression was nonzero but not exactly 1 - see the {it:if} note under Syntax; rows where it was missing are excluded from the sample and, by that same rule, are not counted here), and cluster count. {it:e(N_clust)} is Stata's conventional name for the same count {it:e(Nclus)} carries, and {it:e(clustvar)}, {it:e(vce)} and {it:e(vcetype)} are the conventional variance-type entries, posted so that {cmd:estout}, {cmd:esttab} and anything else keying on them finds what it expects. {it:e(vce)} is {it:cluster}, {it:robust}, {it:unconditional}, {it:bootstrap}, the first word of {cmd:customvce()}, or {it:ols} ({it:oim} under {cmd:poisson}) when none of those applies; {it:e(vcetype)} is {it:Robust} or {it:Bootstrap} where one of those fits and is absent otherwise. Any run that posts no {it:e(V)} - see {it:e(wdid_vcvstate)} - leaves {it:e(vce)} and {it:e(vcetype)} empty, since they name what that matrix is; {it:e(clustvar)} names the design and stands either way. {p_end}

{synopt:R^2 measures} For standard estimation, these include e(r2), e(r2adj), e(r2_within), e(r2_withinadj). For poisson estimation, these include e(r2pseud), a pseudo-R^2 measure. {p_end}

{synopt:{it:e(histogramestimates)}} The cohort i-specific treatment effects used to generate histograms, when requested. Row names carry internal temporary-variable identifiers,
not cohort values; the histogram itself is the intended output.{p_end}

{synopt:{it:e(summarystats)}} Summary statistics on the outcome variable and (if present) continuous treatment variable, when requested. The columns are {it:mean}, {it:median}, {it:sd}, {it:iqr}, {it:p25}, {it:p75}, {it:min} and {it:max}. The rows are fixed keys rather than variable names: {it:y_insample}, {it:y_pretreat}, {it:y_posttreat} and {it:y_control} for the outcome, and {it:c_insample}, {it:c_pretreat}, {it:c_posttreat} and {it:c_control} for the continuous treatment, the {it:c} block only under {cmd:contreat()}, the {it:_control} rows omitted under {it:contreatwithin}, and {it:c_pretreat} and {it:c_control} omitted under {it:contreatcontrols(untreatedZero)}. Thirteen of these are also returned directly as scalars, and only when {opt summarystats} is given: {it:e(ypre_mean)}, {it:e(ypre_med)}, {it:e(ypre_sd)} and {it:e(ypre_iqr)} are the mean, median, standard deviation and interquartile range of the outcome over treated units in the periods before their own onset (the {it:y_pretreat} row) taken with the regression weights; and under {cmd:contreat()} the same three moments of the continuous treatment on three samples, {it:e(contreat_mean)}/{it:e(contreat_med)}/{it:e(contreat_sd)} over the whole estimation sample, {it:e(contreatpost_mean)}/{it:e(contreatpost_med)}/{it:e(contreatpost_sd)} over treated units at or after onset, and {it:e(contreatpre_mean)}/{it:e(contreatpre_med)}/{it:e(contreatpre_sd)} over treated units before it. The {it:contreatpre} trio is not posted under {it:contreatcontrols(untreatedZero)}, following the matrix row that is omitted there for the same reason. {p_end}

{synopt: CFX Plot Info} Information used to generate the continuous treatment effects plots. The information used to produce each plot will be stored in a matrix with the following naming
pattern: e(cfx[item]_[type]), where [item] is either one of att, att_i, att_it, or att_itime (for overall average treatment effects) or one of ame, ame_i, ame_it, or ame_itime (for marginal
effects of the continuous treatment). Meanwhile, [type] is either cfxpost (estimated in the post-treatment period), cfxpre (estimated in the pre-treatment period), or relyrZ (estimated in
relative time period Z, where the letter n represents a negative sign). Under {it:makecfxplotsbysg} each plot's matrix carries a third component, {it:e(cfx[item]_[type]_sg[sgcase])}, where [sgcase] is the subgroup value ({it:e(cfxatt_cfxpost_sg1)}, {it:e(cfxame_relyr2_sg3)}), with the same row keys and column shape as the twins above. That suffix is {cmd:_sg}, and this is the one family of stored results that uses it: the subgroup twins elsewhere use {cmd:_s} ({it:e(escomp_s1)}, {it:e(b_es_att_s1)}) or a bare underscore ({it:e(wdidsgresults_1)}). {it:e(wdid_cfxmats)} lists every cfx matrix name a run actually posted. Matrix rows are keyed {it:lattice_#} (the lattice or {cmd:cfxlattice()} level, in {it:margins} order) with a null
row equation; the same keys apply to the one-column matrices produced under {cmd:cfxplotnose}/{cmd:fastnose}, and AME rows no longer carry the dose-variable equation name. {p_end}

{synopt:{it:e(sample)}} A function characterizing the estimation sample used; useful for identifying which observations are dropped by {it:wooldid} during its sample formation and sanitation
procedures.{p_end}

{synopt:{it:e(cmd)}, {it:e(depvar)}, {it:e(cmdline)}} {it:e(cmd)} holds {it:wooldid} on every successful estimation, which is what a bare {cmd:wooldid} replay and any post-estimation wrapper dispatching on {it:e(cmd)} keys on; {it:e(depvar)} holds the outcome variable, the first variable of the varlist; and {it:e(cmdline)} holds the full {it:wooldid} command line for the successful estimation, beginning with the command name. These are posted, like the conventional variance entries above, so that {cmd:estimates store} and {cmd:estimates restore}, {cmd:estout}, {cmd:esttab} and anything else keying on the conventional eclass identifiers finds what it expects.{p_end}

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
and accepts no options: any non-empty command line - including {it:wooldid, update} - is processed as a normal call. If a {it:wooldid} call fails at the parse or estimation stage, {it:e()} still holds the previous successful estimation (standard Stata behavior); if it fails later - diagnostics, display, or plots, after the results were posted - {it:e()} holds the failed call's own partial post. Either way replay refuses with an informative error until a {it:wooldid} call completes successfully: redisplaying results that predate the failed call would
be misleading. The refusal is session-scoped (an internal {it:WOOLDID_LASTRC} marker; {it:WOOLDID_*} global names are reserved by {it:wooldid}): it also blocks replay of models
explicitly selected with {it:estimates restore} until some {it:wooldid} call completes. A failure inside the replication loop of a {cmd:permute()} or {cmd:bootstrap()} call is the one case that preserves replayability: that call's own completed ordinary estimation is restored, valid, and remains replayable, the marker is cleared, and the call exits with the loop's return code (the engine request itself failed loudly). Everything after the loop behaves like any other post-estimation failure (the rewiring of {it:e()} to the engine's inference, the deferred display, and the deferred continuous-treatment and event-study plots, {cmd:pretrendpower()} and the {cmd:bootsave()}/{cmd:permsave()} write), because {it:e()} is then a half-rewired mixture of observed and replication results, and replay is refused until some {it:wooldid} call completes successfully. A failure during the rewiring reports a message naming that state; a failure in one of the deferred plots is caught and reports the ordinary form. Either way, refit before replaying. Replay requires results posted by a
{it:wooldid} version with replay support; internally, the display state travels in hidden {it:e(wdid_*)} results, which are implementation details and not part of the documented
{it:e()} surface, with the exception of those that carry their own description somewhere in this help, whether in Stored Results or in an option's entry, which {ul:are} documented and supported. Today those are the macros {it:e(wdid_vcvstate)} and {it:e(wdid_vblocks)} and the scalar {it:e(wdid_vblockn)} (see Stored Results), and the cell-collapse results {it:e(wdid_fcol)}, {it:e(wdid_fcskip)} and {it:e(wdid_fcguard)} with the scalar {it:e(wdid_fcncell)} (see {opt fastcollapse}). Being hidden, none of the {it:e(wdid_*)} results appear under a plain {it:ereturn list}; read them by name, or use {it:ereturn list, all}. {p_end}

{pstd} {ul:Incompatible Post-Estimation Commands:} {it:lincom}, {it:test}, and all other post-estimation commands that rely on receiving a valid variance-covariance matrix from
an estimator. While {it:wooldid} does return a coefficient vector ({it:e(b)}) and a variance-covariance matrix ({it:e(V)}), these are for display purposes only. On an analytic run that produced an event study the variance-covariance
matrix is a {ul:mixed} one: within a single event-study family its off-diagonals are the genuine joint covariance of those horizons, so a {it:test} or {it:lincom} across the horizons of one family {ul:is} valid on it; everywhere else - across two families, or between a family and a pooled estimate - they are still 0 and still mean "not estimated". {it:e(wdid_vcvstate)} records what {it:e(V)} is on any given run, and {it:e(wdid_vblocks)} names the families whose blocks were spliced. So while {it:esttab} and the like will present correct estimates using these inputs, efforts to combine or otherwise test different
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
{it:did_imputation} given that one uses the {it:att} estimand-type, uses a pooled reference period (i.e., not using {it:esfixedbaseperiod}), requests no
pre-treatment horizons ({it:esprelength(0)}), does not use unconditional standard errors, and does not use controls. The {it:esprelength(0)} condition is easy
to miss and it matters: the pooled reference is the set of horizons the run asks for, so requesting pre-periods removes them from that pool and shifts every
post horizon by a constant (measured 4.8e-02 on a common-onset panel) while both commands still run cleanly. {it:test/suite/newfeat/bjs_identity.do} pins both
directions - the match at {it:esprelength(0)} and the divergence without it. With controls, {it:wooldid} and {it:did_imputation} will still yield the same results provided
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

{phang}4) Event study with 3 pre-treatment and 3 post-treatment estimates using a pooled reference period and generating event study plots, as well as testing the hypothesis
that all pre-treatment effects are 0 (the {cmd:adversarial} composite test that the estimates and the reference period fall on the same line requires {cmd:esfixedbaseperiod} - the pooled reference has no single omitted period - so here it returns the fitted line without a p-value, under a caution): {p_end}

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

{p 10 10}   {cmd:honestdid, numpre(`=e(numpre_att)') b(b) vcov(V) delta(rm) mvec(0.5(0.5)2) parallel(0)}

{p 8 8}The {cmd:parallel(0)} is not optional housekeeping. {cmd:honestdid} defaults to {cmd:parallel(4)} on any machine
where the {cmd:parallel} package is installed, so each {cmd:delta()} call silently starts four additional Stata
processes - multiplicative inside a loop over cases or over {cmd:mvec} values.{p_end}

{p 8 8}A Rambachan-Roth robust confidence set can legitimately come back {ul:empty} on data with strong pre-trends: under the C-LF construction no parameter value survives the assumed bound on post-treatment violation. That is a finding about the data - the pre-trend evidence rules out every value at that {it:M} - and not a failure of the handoff or a sign that the wrong matrix was passed. Observed on real fixtures; widen {it:mvec()} to see where the set becomes non-empty.

{phang}11) Fast exploratory event study (point estimates only), also returning the cohort-by-relative-period building-block table:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) espre(2) espost(2) fastnose gteffects}

{phang}12) Randomization inference: permute treatment-onset dates across clusters, holding onset-group counts fixed, with 500 replications (set the seed for reproducibility); every reported p-value and interval is then permutation-based, and the run is a fast (no-SE) run unless {cmd:permtstat} or {cmd:permanalyticse} is added:{p_end}

{p 10 10}   {cmd:set seed 12345}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(500)}

{phang}13) Cluster bootstrap standard errors with 500 replications, reproducible, and joint pre-trend tests formed from the replications' own test statistics ({cmd:boottstat}; without it they are instead the fixed-metric rank test of the replication draws, as {it:e(boot_jt_family)} records) - the sanctioned route when clustering at or above the cohort level:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) jointtests bootstrap(standard) reps(500) rseed(42) boottstat}

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
