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

{pstd} {bf:wooldid} -- Suite for difference-in-differences analyses with staggered treatment onset, using the two-way fixed effects approach of Wooldridge (2025)
and the high dimensional fixed effects estimators of Correia (2017). It estimates several types of average treatment effect, compares effects across subgroups,
produces event study estimates and parallel trends tests, and estimates average and marginal effects of a continuous treatment or treatment intensity variable [experimental].
It draws event study plots, histograms of cohort effects, and plots of effects across the distribution of a continuous treatment. Inference options are analytic
clustered standard errors, unconditional standard errors that account for sampling variation in covariates, Ibragimov and Muller (2010) cluster-robust inference,
bootstrap standard errors, and randomization inference for settings with few treated clusters. Other features include {it:anticipation()} windows, a {it:gteffects}
table of cohort-by-period effects, {it:pretrendpower()} diagnostics, a point-estimates-only {it:fastnose} mode, and replay of the stored results by typing {it:wooldid} alone. {p_end}

{pstd}{ul:Shortcuts:} {help wooldid##description:Description}, {help wooldid##options:Options}, {help wooldid##output:Working with Output},
{help wooldid##runtime:Runtime}, {help wooldid##comparison:Comparison with Related Estimators}, {help wooldid##examples:Examples},
{help wooldid##bugs:Bugs and Assistance}, {help wooldid##citations:Citations}, {help wooldid##author:Author Info and How to Cite}, {help wooldid##disclaimer:Disclaimer} {p_end}

{pstd}{ul:Dependencies:} {cmd:wooldid} requires Stata 16 or higher, {help reghdfe} version 5.7.3, {help ftools},
and {help gtools}. Poisson estimation also requires {help ppmlhdfe}. {p_end}

{pstd}{ul:Github Link:} {browse "https://github.com/thegland/wooldid/"}. The command {it:wooldid, update} updates this program from that page. A longer guide
for AI agents and advanced users is at {browse "https://github.com/thegland/wooldid/blob/main/wooldidinstall/wooldid_agent_guide.md"}. It is also an ancillary file of the package,
copied to the working directory by {cmd:net get wooldid}. {p_end}

{marker syntax}{...}
{title:Syntax}

{phang}{cmd: wooldid} {it:y i t ttre} [if] [aw/pw=weights] [{cmd:,} {help wooldid##specifyestimand:estimands} {help wooldid##standarderrors:standard_errors}
{help wooldid##bootstrap:bootstrap_standard_errors} {help wooldid##permute:randomization_inference} {help wooldid##supporting:supporting_features} {help wooldid##eventstudy:event_study}
{help wooldid##interactions:subgroup_effects} {help wooldid##continuous:continuous_treatment} {help wooldid##controls:controls_and_FE} {help wooldid##technical:technical_features} ]{p_end}

{phang}{ul:The if qualifier and weights:} {cmd:wooldid} keeps the observations for which the {it:if} expression equals exactly 1, where Stata's own {it:if} keeps any nonzero value; write {cmd:if flag != 0} for Stata's behaviour.
Weights may be aweights or pweights; frequency weights (fw) are not accepted. {p_end}

{synoptset 22 tabbed}{...}
{synopt : {it:y}}Outcome variable. {p_end}
{synopt : {it:i}}Cohort identifier; must consist of integers greater than 0, and at most 2,147,483,647 on units with a treatment onset. Cohorts can be fine (i = individual units) or coarse (i = collections of individual units), but should be no coarser
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

{synopt : {help wooldid##bootstrap:Bootstrap Standard Errors:}} {opt boot:strap(standard|bayes)} {opt reps(int)} {opt bootts:tat} {opt bootdisp:lay(spec)} {opt bootd:f(spec)}
{opt bootsa:ve(name)} {opt rseed(#)} {opt dots(int)} {p_end}

{synopt : {help wooldid##permute:Randomization Inference:}} {opt perm:ute(shuffle|probabilities)} {opt permreps(int)} {opt permstr:ata(variable)} {opt permpr:obs(spec)} {opt permex:act}
{opt permts:tat} {opt permt:ci} {opt perman:alyticse} {opt permab:sp} {opt permsi:de(left|right)} {opt permstep:down} {opt permonset:cohorts} {opt permsa:ve(name)} {p_end}

{synopt : {help wooldid##supporting:Supporting Features:}} {opt hist:ogramcohorteffects} {opt save:plots(name)} {opt replace} {opt sum:marystats}
{opt pois:son} {p_end}

{synopt : {help wooldid##eventstudy:Event Study Estimation:}}  {opt esf:ixedbaseperiod} {opt espre:length(int>=0)} {opt espost:length(int>=0)} {opt esr:elativeto(-int)}
{opt joint:tests}   {opt adv:ersarial} {opt altmethod_adversarial} {opt makep:lots} {opt oos:treatedcontrols(int>=0)} {opt antic:ipation(int>=0)} {opt pretrendpow:er(numlist)}  {p_end}

{synopt : {help wooldid##interactions:Subgroup Effects:}} {opt sub:groups(variable[>=0 int])}  {p_end}

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
{opt update} {opt customw:eightmultiplier(variable[>=0 numerical])} {opt regtol(+int)} {opt ccno:absorb} {opt ccc_absorb(varlist[numerical])}
{opt ver:bose} {opt gte:ffects} {opt safety(on|off)} {p_end}

{marker description}{...}
{title:Description of Estimation Approach}

{pstd}
{bf:wooldid} estimates difference-in-differences treatment effects in the vein of Wooldridge (2025), with extensions for continuous treatments and for effects by
subgroup or treatment arm. Treatment must be irreversible: once treated, cohorts must remain treated.{p_end}

{pstd}
{it:wooldid} works in two steps. Step 1 estimates the underlying two-way fixed effects regression. At baseline this regresses the outcome on cohort fixed effects,
time fixed effects, and one indicator for each cohort-period (i-t) cell whose treatment effect is to be estimated. Step 2 uses {help margins} to turn the coefficients
on those indicators into the treatment effects of interest. For every treated observation it takes the difference between the fitted value and the counterfactual
fitted value without treatment, and then averages these differences. The underlying regression is: {p_end}

{p2col 5 8 8 0 : Step 1 (Underlying Regression):} Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it) {p_end}

{pstd} where i indexes cohorts, t time periods and j individual units, which may coincide with cohorts or be more numerous. Z is the set of (i,t) pairs where i is a
treated cohort and t is one of its treated periods. Tz_it is 1 in the cell z=(iz,tz) of Z and 0 otherwise, and TreatedFlag_it is 1 for any observation in a cell of Z.
fe_i and fe_t are cohort and time fixed effects, and Y_jt is the outcome. {p_end}

{pstd} The margins calculation in step 2 depends on the estimands requested. The overall average treatment effect averages, over all observations with TreatedFlag_it
equal to 1, the fitted value less the counterfactual fitted value with TreatedFlag_it set to 0, applying any sample weights. Other estimands reweight this average,
or restrict it to a subgroup or time period. See {help wooldid##specifyestimand:estimands}. {p_end}
{pstd} By default, standard errors come from the underlying regression and pass through the margins calculation by the delta method. They are conditional
standard errors: they take the observed covariate values as fixed. P-values use a t distribution with the degrees of freedom stored in {it:e(df_r)}, described under
{help wooldid##output:Working with Output}. Under {it:poisson} they use the normal distribution. {it:unconditionalse} instead accounts for sampling variation in
covariates. See {help wooldid##standarderrors:standard_errors}. {p_end}

{pstd} The approach extends with little change. {it:poisson} swaps the regression command without changing the equation or the margins step. Controls and fixed
effects enter the regression without changing the margins step. Subgroup effects come from running the margins step within each subgroup, and {it:wooldid} can test
whether subgroup effects differ, for example across treatment arms, as interaction terms would in a traditional difference-in-differences. See
{help wooldid##interactions:subgroups}. {p_end}

{pstd} Event study estimates run the margins step on the observations in each period relative to treatment. Pre-treatment effects also need extra indicators,
for the cohort-periods of treated cohorts in the pre-treatment horizons requested. The equation is unchanged: only the set Z grows. The joint pre-trend tests and the
adversarial line test are described under {help wooldid##eventstudy:event_study}.{p_end}

{pstd} Event studies use a pooled reference period by default. It pulls each treated cohort's onset earlier by just enough periods to estimate the requested
horizons, adds no more (i,t) pairs to Z than needed, and so may pool several periods into each cohort's reference period. Its event study is normalized to 0 in the
earliest period reported. Under {it:esfixedbaseperiod}, Z holds every (i,t) pair of the treated cohorts except the reference period, by default the period just before
treatment. The regression is then saturated with i-t indicators, and every effect is measured against that one period. See {help wooldid##eventstudy:event_study}. {p_end}

{pstd} With more change, the approach handles continuous treatments and treatment intensity. When cohorts are coarse, so that cells hold several units, and the
dose varies within i-t cells, the effect of the dose can be estimated within each treated cell, using a flexible function of the dose where needed. This extension is
experimental and is not covered in Wooldridge (2025), though Wooldridge has
{browse "https://twitter.com/jmwooldridge/status/1695115782739435922?s=61&t=pXC6R4euc7LaN6MMB5VApQ":informally suggested}
interactions between treated cell indicators and the continuous treatment. Reading these within-cell effects as the effect of a larger dose needs more than parallel
trends. See {help wooldid##continuous:continuous_treatment}. {p_end}

{pstd} The continuous treatment regression is: {p_end}
{p2col 5 8 8 0 : Step 1 (Continuous):} Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it * f(Contreat_jt)) + g(f(Contreat_jt)) {p_end}

{pstd} where Contreat_jt is the continuous treatment, f() is a flexible function of it, and g() interacts f() with fixed-effect-like indicators. Average treatment
effects use the same margins calculation as before. The average marginal effect shifts Contreat_jt rather than TreatedFlag_it. {it:wooldid} can also report average
treatment effects along the distribution of the dose, and the effect of a proportional rather than a unit change in the dose. See {help wooldid##continuous:continuous_treatment}. {p_end}

{pstd} {it:wooldid} also offers Ibragimov and Muller (2010) cluster-robust inference. The regression is unchanged, but the margins step computes each effect,
average treatment effect or effect of a continuous treatment, separately within each cluster. The estimate is the simple average of the cluster estimates, and
inference t-tests them. See {help wooldid##standarderrors:standard_errors}. {p_end}

{pstd}
{ul: Comment on Sample Formation, Cohorts Treated Prior to Sample Start:} The estimation sample can be smaller than your sample. {it:wooldid} drops observations with
missing data or zero weights, treated cohorts not observed before treatment, or not observed in the reference period under {it:esfixedbaseperiod}, and time periods
in which no untreated or not yet treated unit is observed. Under {it:timetrends} it also drops cohorts with too few pre-treatment periods to estimate a trend. {p_end}

{marker options}{...}
{title:Options}

{marker specifyestimand}{...}
{dlgtab:Specifying the Estimand(s)}

{phang}{opt att}: Estimate average treatment effects on the treated (ATT). This is the default, and it may be combined with {it:att_it}, {it:att_i} and {it:att_itime}.
The margins calculation, step 2 in {help wooldid##description:Description}, weights observations by the regression's sample weights, if any, times any
{it:customweightmultiplier()} variable. That multiplier enters the margins weights only, not the regression weights.{p_end}

{phang}{opt att_it}: Estimate the simple, unweighted mean of the cohort-period (i-t) cell treatment effects (ATT_it). Each cell's effect uses the sample weights
of its own observations. With fine cohorts, where each cell holds one observation, this ignores the sample weights.{p_end}

{phang}{opt att_i}: Estimate the simple, unweighted mean of the cohort-specific treatment effects (ATT_i). Each observation is weighted by its share of its cohort's
sample weight in the post-treatment period, or in the pre-treatment period outside the reference period. Under {it:timetrends} the period immediately before the
reference period is also normalized to zero and excluded from the pre-treatment weight.{p_end}

{phang}{opt att_itime}: Estimate the simple average across cohorts of each cohort's simple average cell effect (ATT_itime). Cell effects use sample weights,
as for {it:att_it}. Each cohort's effect is the unweighted mean of its cell effects in the post-treatment periods in which it is observed, or in its pre-treatment
periods outside the reference period. ATT_itime is the unweighted mean of these cohort effects.{p_end}

{pstd}Under {it:contreat()} each requested estimand is also reported as the average marginal effect of the continuous treatment, labelled {it:ame}, {it:ame_it}, {it:ame_i} and {it:ame_itime}.
Because {it:att_it}, {it:att_i} and {it:att_itime} weight cells or cohorts equally rather than observations, some options treat them differently. {it:unconditionalse} is
refused for them unless {it:safety(off)} is specified, and {it:bootstrap(bayes)} is refused for them. {it:makecfxplots} produces its plots and the {it:e(cfx*)}
matrices for {it:att} only unless {it:safety(off)} or {it:contreatwithin} is specified. A {it:subgroups()} variable that splits a unit across subgroups over time is refused.
Except under {it:permute()}, {it:safety(off)} admits such a unit by its weight share in each subgroup, and {it:att_i} and {it:att_itime} then no longer weight units equally.{p_end}

{marker standarderrors}{...}
{dlgtab:Standard Errors}

{phang}{opt clus:ter(clusvars)}: Variable(s) specifying the clustering level. Clustered standard errors are computed analytically. Several variables give multi-way
clustering on the default path. {it:unconditionalse}, {it:impvalues}, {it:bootstrap()} and {it:permute()} each require a single cluster variable. A cluster variable that
leaves fewer than two clusters in the estimation sample stops the run with an error. A finer cluster variable, {it:robust} or {it:fastnose} are the ways out. Cluster at
the level at which treatment is assigned, or coarser, and never below the unit whose rows are correlated over time. The discussion of what level to cluster at below
adds a second condition specific to {it:wooldid}. {p_end}

{phang}{opt rob:ust}: Heteroskedasticity robust standard errors, relevant only without {it:cluster()}. Without {it:cluster()}, {it:robust} or {it:customvce()}, the
linear regression reports conventional standard errors, which assume independent, homoskedastic errors. With pweights, or under {it:poisson}, the default is robust instead.
None of these defaults allows for correlation within a unit over time, so on panel data they usually overstate precision. Specify {it:cluster()} at the level at which
treatment is assigned. {p_end}

{phang}{opt uncond:itionalse}: Compute standard errors that account for sampling variation in covariates. The regression then runs through {it:regress} or
{it:poisson} rather than the high dimensional fixed effects packages, which is usually slower. {it:fe()} and {it:ccc_absorb()} are not absorbed on this path.
{it:wooldid} cautions and enters them as written, so put them in {it:controls()} with your own {it:i.} prefixes. Requires {it:robust} or clustered standard errors
with a single cluster variable. A robust or clustered {it:customvce()}, or pweights, also qualify. Unless {it:safety(off)} is specified, it is refused with
{it:att_i}, {it:att_it} and {it:att_itime} and cautions with aweights. It is refused with {it:impvalues}, {it:fastnose}, {it:bootstrap()} and {it:permute()}. It
turns the default cell collapse off. {p_end}

{phang}{opt customvce()}: Estimate the underlying regression with {help regress} or {help poisson} and pass the text unchanged into their {it:vce()} option.
This allows standard error methods that the default high dimensional fixed effects packages do not offer. Not compatible with {it:robust} or {it:cluster()}: specify
the clustering or robustness inside the {it:vce()} text itself. Like {it:unconditionalse}, it cannot absorb {it:fe()} or {it:ccc_absorb()}, which are entered as
written under a caution. Refused with {it:impvalues} and {it:bootstrap()}. It cannot be combined with {it:permute()}, which needs a {it:cluster()} variable to
randomize over. {p_end}

{phang}{opt imp:values}: Ibragimov and Muller (2010) style cluster-robust inference for all estimates. Each effect is estimated within each cluster, the
cluster estimates are averaged, and inference t-tests them against 0. Requires a single cluster variable at the level of {it:i} or coarser. A {it:cluster()} finer
than {it:i}, or one that cuts across it, is refused, since the cluster estimates are then re-weightings of the same coefficients and the reported standard error is
far too small. The target estimand changes: the cluster estimates are averaged with equal weight, so what is reported is the cluster-equally-weighted mean whatever
estimand was requested. Clustering at the cohort level, for example, makes {it:att} equal to {it:att_i}. Each reported row uses a t distribution with degrees of
freedom one less than the number of clusters contributing to that row, and {it:e(df_r)} is the smallest of these. Not allowed with subgroups or with {it:jointtests},
{it:makecfxplots}, {it:adversarial}, {it:customvce()}, {it:unconditionalse}, {it:pretrendpower()}, {it:bootstrap()} and {it:permute()}. {it:fastnose} is accepted with
a caution and has no effect. No event study vector or covariance is returned: {it:e(b_es_*)}, {it:e(V_es_*)} and {it:e(numpre_*)} are absent, and the {it:honestdid}
handoff is unavailable. {p_end}

{phang}{opt fast:nose}: Fast exploratory mode, for specification search on large panels. Every {it:margins} call runs without variance estimation, so only point
estimates are produced. Standard errors, test statistics, p-values, confidence intervals, joint tests and plot confidence bands are omitted and display as missing.
{it:e(V)} is not posted, so {it:test} and {it:lincom} refuse, and {it:esttab} prints an empty inference row. Refused with {it:unconditionalse}, {it:pretrendpower()},
{it:boottstat}, {it:permtstat}, {it:permtci}, {it:permanalyticse}, {it:adversarial} and {it:jointtests}. The last two are accepted under {it:permute()}, where fast mode
is already the default. With {it:impvalues} or {it:bootstrap()} it is accepted with a caution and has no effect, since {it:impvalues} already runs {it:margins} without
variance estimation and the bootstrap replications already run in fast mode. {it:simultaneousbands} is accepted but no band is built. Re-run the chosen specification
without {it:fastnose} before reporting results. {p_end}

{phang}{opt nofastcollapse}: Turns off the default cell collapse, which speeds up {it:wooldid} considerably by evaluating the {it:margins} calls on one row per
cohort-by-period-by-treatment cell rather than on every observation. The regression itself still runs on the full data. Results agree with an uncollapsed run to
{it:margins}' own accuracy, typically machine precision, but up to about 1e-6 relative on standard errors with sampling weights, controls and an event-study window,
where neither run is exact. The collapse applies to binary and continuous treatments, to poisson, to subgroups and inside {it:bootstrap()} and {it:permute()}
replications. It is skipped, with a caution naming the option responsible, under {it:interactivecontrols()}, {it:contreatpoly()} above 1, the {it:contreat()} variable
also appearing inside {it:controls()}, {it:coarsecohortcontrols()}, {it:fe()} or {it:ccc_absorb()}, and {it:unconditionalse}. Disable it yourself when controls or fixed
effects interact with treatment in a way {it:wooldid} cannot detect, such as a control built from a renamed copy of the treatment variable, since the collapse assumes
the treatment effect is constant within a cell. Disable it also when memory is tight, since the collapse holds an extra copy of the data, and when you need byte-for-byte
reproduction of an uncollapsed result. {p_end}

{phang}{opt simult:aneousbands}: Compute a simultaneous (sup-t) confidence band over each event study's horizons at {it:level()}. The printed table keeps its
pointwise intervals. The band endpoints are stored in {it:e(esband_*)} and drawn as a dashed pair on the event study plots. Without {it:bootstrap()} the band is
simulated from the family's analytic event study covariance, on the same t distribution as the pointwise intervals where the regression reports degrees of freedom.
With {it:bootstrap()} it is taken from the replications. With {it:cluster()} at or above the level of {it:i} the analytic band is too narrow, and a caution says so.
Use {it:bootstrap()} or a {it:cluster()} strictly finer than the cohort-by-period cells there. Refused with {it:permute()}, which offers {it:permstepdown} instead. No
band is built under {it:fastnose} or for a family with fewer than two horizons with a usable variance, and a caution says so. {p_end}

{phang}{opt level(#)}: Confidence level for every reported confidence interval, including the plots, the {it:gteffects} cells and the {it:bootstrap()} and
{it:permute()} intervals. It must lie between 10 and 99.99 with at most two digits after the decimal point. The default is {it:c(level)}, normally 95. The interval
columns of the stored results are labelled {it:lb_#ci} and {it:ub_#ci} for the level in force, which is stored in {it:e(level)}. The pre-test sizes reported by
{it:pretrendpower()} are fixed at 5%. {p_end}

{pstd} {ul: Discussion of What Level to Cluster At:}
With {it:cluster()} at the level of {it:i} or coarser, each treated cohort-by-period cell lies wholly inside one cluster, and the saturated cell dummies leave the
cluster-robust covariance nothing within the cluster to identify. The analytic standard errors of the {it:ame} estimates are then not identified and are reported as
missing. The analytic joint pre-trend tests and the adversarial line test are not calibrated. Pairwise subgroup contrasts can be degenerate. The {it:att} rows'
analytic standard errors are understated without a caution, and the {it:gteffects} cells more so, which a {it:gteffects} run counts under a caution. The analytic
{it:simultaneousbands} band is scaled by the understated event study standard errors, so it is too narrow. The remedies are a {it:cluster()} strictly finer than the
cohort-by-period cells, {it:bootstrap()}, or {it:permute()}. {it:permute()} refuses the band, and for the {it:gteffects} cells only the finer {it:cluster()} is a
remedy. A finer {it:cluster()} is a remedy only when it still sits at or above the level at which treatment is assigned. When no such level is finer than the cells,
use {it:bootstrap()} or {it:permute()}. {p_end}

{pstd} {ul: Discussion of Unconditional Standard Errors:}
The effects {it:wooldid} reports depend on the regression coefficients and also on the observed values of any covariates or continuous treatment. The default
delta-method standard errors take those values as given. Unconditional standard errors also account for their sampling variation, as in the discussion of
{it:vce(unconditional)} in {help margins}. Which is appropriate likely varies from context to context. {p_end}

{pstd} Wooldridge (2025) expects conditional and unconditional standard errors often to be similar. They are similar when treatment effects vary little
across cohort-time period cells, covariate values and treatment intensities. Where effects are heterogeneous the unconditional standard errors are larger, and the ratio
between them does not shrink as the sample grows. The conditional standard errors describe the effect averaged over the treated observations in the sample. The
unconditional standard errors describe the population average. They also matter more when sample weights make particular observations very influential. With
unconditional standard errors, subgroup and event study estimates can be sensitive to observations outside the focal subgroup. {p_end}

{pstd} {ul: Caution Regarding Unconditional Standard Errors When Used with Certain Weights and Estimands:} With aweights, or with estimands other than the ATT,
tests using unconditional standard errors are sometimes badly oversized, rejecting the null hypothesis at very high rates. The aweights problem does not always occur, but where
it does it affects {help hdidregress} as well. Unless {it:safety(off)} is specified, {it:wooldid} therefore cautions against unconditional standard errors with aweights
and stops estimation if they are used with a non-ATT estimand. Unconditional standard errors tend to be highly sensitive to the choice of weights. {p_end}

{pstd} {ul: Discussion of Ibragimov and Muller-Style Inference:}
The approach can improve on analytic clustered standard errors when clusters are few and very heterogeneous in size. It assumes the cluster-specific estimates are
independent. That is most plausible when a continuous treatment's effect is estimated from within-cluster variation only, or in a bespoke design where clusters contain
both treated and control units and do not share time fixed effects. {it:wooldid} reads each cluster's estimate over its treated post-treatment cells alone, so control
units inside a cluster do not give it its own control arm unless, for example, cluster-by-time effects are placed in {it:fe()}. In the standard design {it:wooldid}
still reports IM-style estimates, though the shared control group, and a shared reference period without {it:esfixedbaseperiod}, likely violate the independence
assumption. {p_end}

{pstd} Replication-based alternatives to the analytic family are described under {help wooldid##bootstrap:Bootstrap Standard Errors} and {help wooldid##permute:Randomization Inference}. {p_end}

{marker bootstrap}{...}
{dlgtab:Bootstrap Standard Errors}

{phang}{opt boot:strap(mode)}: Bootstrap inference. Every reported standard error, test statistic, p-value and confidence interval becomes a bootstrap quantity,
including the joint tests, subgroup contrasts, plots, {it:gteffects} cells and {it:pretrendpower()}. Point estimates are unchanged, and the analytic covariance is kept
in the hidden {it:e(V_analytic)}. Each replication re-estimates the full model. {it:bootstrap(standard)} is the pairs cluster bootstrap, resampling the {it:cluster()}
groups, or the observations without {it:cluster()}. {it:bootstrap(bayes)} is the Bayesian bootstrap. It multiplies each cluster's observations by a Dirichlet weight, so
no cluster ever leaves the sample. Requires {it:reps()}. Refused with {it:permute()}, {it:impvalues}, {it:unconditionalse}, {it:customvce()} and multi-way
{it:cluster()}. {it:fastnose} is accepted with a caution and has no effect, except with {it:boottstat}, where it is refused. {it:bootstrap(bayes)} is refused with
{it:att_it}, {it:att_i} and {it:att_itime}, since it omits the cell-sampling component of their variance. {it:wooldid} cautions when fewer than 6 clusters are treated.
Neither mode is a remedy for few treated clusters, and no wild cluster bootstrap is offered. Use {it:permute()} there, or {help boottest} on the underlying regression
outside {it:wooldid}. When {it:cluster()} does not nest {it:i}, some unit spans two clusters and a cell-weighted estimand is requested, a caution says its standard
errors may be understated. The remedies are a finer {it:i} or a {it:cluster()} that nests {it:i}. When {it:cluster()} nests {it:i} and the run carries {it:att_it},
{it:att_i}, {it:att_itime} or {it:boottstat}, each drawn copy of a cluster gets fresh unit codes, also where {it:i} is named in {it:fe()}, {it:controls()},
{it:interactivecontrols()}, {it:coarsecohortcontrols()} or {it:ccc_absorb()}. A variable that encodes {it:i} under another name, or was merged in on {it:i}, keeps its
original codes. {p_end}

{phang}{opt reps(#)}: Number of bootstrap replications. Required with {it:bootstrap()} and refused without it. There is no minimum, so a small value can time
a specification. Replications that fail to estimate are excluded and counted in {it:e(boot_failed)}. Above 10,000 completed replications the draw matrices are not posted,
and a caution says so. {it:bootsave()} writes them regardless. With {it:n} completed replications the smallest printable p-value is 2/{it:n} for the percentile coefficient
p-values, 2/({it:n}+1) under {it:boottstat}, and 1/({it:n}+1) for the joint tests. A caution names the limit whenever it exceeds .05. The {it:bc} and {it:bca} p-values
invert their intervals and have no such limit. They can fall below 2/{it:n}, and they and the percentile p-values can print exactly 0. {p_end}

{phang}{opt bootd:f(spec)}: The reference distribution for the Wald family of bootstrap statistics. Refused without {it:bootstrap()}. {it:bootdf(normal)} uses
standard normal critical values and reports z statistics. {it:bootdf(clusters)} uses the t distribution with {it:G}-1 degrees of freedom, {it:G} the number of clusters in
the estimation sample, and is refused without {it:cluster()}. A positive integer sets the degrees of freedom directly. The default is {it:bootdf(clusters)} with
{it:cluster()} and {it:bootdf(normal)} otherwise. Only the Wald family p-values and intervals change, and under {it:boottstat} nothing printed does. {it:e(df_r)} follows
the same choice, so {it:test} or {it:lincom} run afterwards on {it:e(V)} reproduces the Wald family. Under the default {it:bootdisplay(percentile)} that is not the printed
p-value, and the footer says so. {p_end}

{phang}{opt bootdisp:lay(spec)}: Which p-value and confidence interval the printed tables carry. Refused without {it:bootstrap()}. {it:percentile} is the default.
{it:bc} is the bias-corrected percentile method. {it:bca} is bias-corrected and accelerated, requires {it:cluster()}, and costs {it:G} extra fits for its delete-one-cluster
jackknife. {it:wald} is a normal or t approximation set by {it:bootdf()}. Every family except {it:bca} is computed and stored in {it:e()} on every run. This option only
chooses what is displayed. It never changes the reported standard error, the standard deviation of the bootstrap distribution. The percentile and bias-corrected intervals
are not symmetric about the estimate, and the {it:t_wald} or {it:z_wald} column beside them is for reference only. The percentile and bias-corrected p-values can print exactly 0, meaning no completed replication fell beyond the null. Under {it:percentile} a line under the table
gives the resolution. Under {it:bootdisplay(bca)} the {it:gteffects} and continuous treatment
cells carry the {it:bc} family. With {it:boottstat} only {it:percentile} is available. {p_end}

{phang}{opt bootts:tat}: A studentized, or percentile-t, bootstrap of each replication's t statistic on its own analytic standard error. Refused without
{it:bootstrap()}. The reported standard error is the analytic one, labelled {it:se_analytic}, and the statistic is labelled {it:t_pt}. The p-value and the interval come
from the ranks of the studentized draws, and the interval is not symmetric about the estimate. Studentization covers every coefficient, the {it:gteffects} and continuous
treatment cells and the subgroup contrasts. The joint pre-trend, subgroup equality and {it:adversarial} tests are bootstrapped as test statistics. Without {it:boottstat}
the adversarial tests are reported as absent. Refused with {it:bootstrap(bayes)}, {it:poisson}, {it:bootdisplay(wald)}, {it:bootdisplay(bc)} and {it:bootdisplay(bca)}.
It is materially slower, since every replication estimates its own variance. With {it:cluster()} at or above the level of {it:i} the studentized contrasts and the subgroup
equality test may be absent or unstable. Use the plain bootstrap for them there. {p_end}

{phang}{opt bootsa:ve(name)}: Also write the replication draws to {it:name}.dta in the working directory. Refused without {it:bootstrap()}. {it:name} is a plain
file name, and an existing file is overwritten only with {it:replace}. The file has a row {it:_rep} = 0 with the observed estimates and one row per replication, {it:_rep} =
1 to {it:reps()}, with the coefficients under their {it:e(b)} names. A failed replication's row is all missing. Under {it:boottstat} the replication analytic standard
errors are added as {it:se_*}. With {it:cluster()}, one {it:w_} column per cluster records how many times it was drawn, or its Dirichlet weight under {it:bootstrap(bayes)}.
The file's characteristics record the engine, seed and random-number state. If the file cannot be written, the run still completes with a caution and posts the return code
in {it:e(boot_save_rc)}. {p_end}

{phang}{opt dots(#)}: Show a progress dot every {it:#} replications. The default is 1, and {it:dots(0)} suppresses the display. After the fifth replication a
projected runtime for the loop is printed on the dots line. Requires {it:bootstrap()} or {it:permute()}. {p_end}

{phang}{opt rseed(#)}: Seed for the bootstrap or permutation draws, an integer between 0 and 2,147,483,647 passed to {it:set seed}. Requires {it:bootstrap()}
or {it:permute()}. The seed is stored in {it:e(boot_seed)} or {it:e(perm_seed)}. The random-number state before the first draw is always stored in {it:e(boot_rngstate)} or
{it:e(perm_rngstate)}, so a run without {it:rseed()} can be replayed with {it:set rngstate} before repeating the call. {p_end}

{marker permute}{...}
{dlgtab:Randomization (Permutation) Inference}

{phang}{opt perm:ute(method)}: Randomization inference over the assignment of treatment onset, for settings with few treated clusters. {it:shuffle} permutes the
observed onset labels across clusters, holding the count at each onset fixed. {it:probabilities} draws each cluster's onset from the probabilities in {it:permprobs()}.
Requires {it:permreps()} or, for shuffle only, {it:permexact}. Requires a single {it:cluster()} variable, the randomization unit, that nests {it:i} unless
{it:permonsetcohorts} is specified. Every cluster must sit inside one treatment-onset group. The never-treated are their own group, coded as plain missing. This screen covers
every row the {it:if} expression keeps, even rows later dropped for a missing outcome or weight. A specification that runs without {it:permute()} can therefore be refused
with it. Restrict the sample with {it:if}. An onset more than {it:oostreatedcontrols()} periods after the last period, net of {it:anticipation()}, is permuted as
never-treated, as the estimation treats it, and a caution counts the clusters affected. Refused with {it:impvalues}, {it:unconditionalse}, {it:contreatwithin},
{it:contreatcontrols(untreatedZero)}, {it:pretrendpower()}, {it:simultaneousbands} and {it:bootstrap()}, and when the design admits only one possible assignment. With
{it:subgroups()} the subgroup screens apply over the whole sample. A unit split across subgroups, or labelled on some rows and missing on others, is refused whatever
{it:safety()} is. Every p-value and interval is permutation-based, including the subgroup contrasts, the
joint pre-trend, subgroup equality and adversarial tests, and the continuous-treatment plot cells. The subgroup contrasts and equality tests test the sharp null of no effect
at any horizon, not equality of effects. P-values are two-sided and equal-tailed, with the add-one form for sampled draws. Intervals are shift intervals at {it:level()}:
the estimate minus the equal-tailed quantiles of the permutation distribution. By default the run and the replications are in fast mode. The {it:se_analytic} column is
then blank, and {it:e(V)}, {it:e(b_es_*)} and {it:e(V_es_*)} are not posted. {it:permanalyticse} restores them. The {it:gteffects} cells, whose onsets the replications reassign, and the
continuous-treatment plot cells of the {it:baseTreated}, {it:byTreated} and {it:byTtre} lattice types have no permutation distribution. They are left blank, not filled
analytically. Replications that fail to estimate are excluded and counted in {it:e(perm_failed)}. The run is refused when fewer than max(20,
half) complete. A caution prints when fewer than 6 clusters are ever treated. The test's size then rests on their being exchangeable in assignment with the others, and a
few atypical treated clusters can mis-size it. Few treated clusters do not by themselves make the possible assignments few. The number of possible assignments is in {it:e(perm_npossible)}. Under {it:permute(shuffle)} without {it:permstrata()} it is
small only when there are few clusters in all or only one or two are treated. Strata make it smaller. {p_end}

{phang}{opt permreps(#)}: Number of sampled reassignments, drawn with replacement. At least 20. Required unless {it:permexact} is specified. A sampled-draw
p-value cannot fall below 2/({it:permreps}+1) under the default convention, or 1/({it:permreps}+1) under {it:permabsp} or {it:permside()}. A caution prints when that floor
is above .05, and another when {it:permreps()} exceeds 10 percent of the possible assignments. {p_end}

{phang}{opt permex:act}: Enumerate every distinct assignment instead of sampling. Only with {it:permute(shuffle)}, and not together with {it:permreps()}. A
caution prints above 1,000 assignments, and enumeration is refused above 65,534. Above 10,000 enumerated assignments the draw matrices are not posted, and a caution says so.
{it:permsave()} writes them regardless. P-values are then plain shares over the {it:M} assignments, with a floor of 2/{it:M}, or 1/{it:M} under {it:permabsp}. {p_end}

{phang}{opt permstr:ata(varname)}: Confine the shuffle to within strata of a numeric variable, so each onset group keeps its composition across strata. Every
cluster must sit inside a single stratum, with no missing values. Only with {it:permute(shuffle)}. Under {it:permute(probabilities)} express strata through {it:permprobs()}
instead. {p_end}

{phang}{opt permpr:obs(spec)}: The assignment-probability design for {it:permute(probabilities)}. Write it as {it:(onsetvar level level ...)} followed by one
{it:(clustervalue p p ...)} group per cluster. The first group names the timing variable and lists every in-sample onset level. The never-treated, including onsets after
the last period, are {it:.}. Each later group gives one cluster's probabilities over those levels, nonnegative and summing to 1. Every cluster value must appear exactly
once. Requires {it:permreps()}. {it:permexact} is refused. A caution prints if the observed assignment has zero probability under the design. {p_end}

{phang}{opt permts:tat}: Randomize on t-ratios instead of coefficients: each replication's estimate over its own analytic standard error. The run is slower, and
{it:fastnose} is refused. Intervals become the estimate minus the equal-tailed quantiles of the permuted t-ratios times
the analytic standard error. The joint, equality and adversarial tests and the pairwise contrasts are unchanged. It protects against a treated cluster with unusually
dispersed errors only where each treated cohort-period cell contains several clusters. Under {it:permute()} that needs {it:permonsetcohorts} with a {it:cluster()} finer
than {it:i}. Whenever {it:cluster()} nests {it:i} it does not protect, and at {it:cluster()} equal to {it:i} a caution says so. {p_end}

{phang}{opt permt:ci}: Also report, for the pooled estimands, the symmetric inverted-|t| interval: the estimate plus or minus the {it:level()} quantile of the permuted absolute t-ratios
times the observed standard error. Implies {it:permtstat}. {p_end}

{phang}{opt perman:alyticse}: Without {it:permtstat}, compute the analytic standard errors anyway and show them in the {it:se_analytic} column for reference.
The p-values and intervals remain permutation-based, and the posted {it:e(V)} is analytic. {it:fastnose} is refused. {p_end}

{phang}{opt permab:sp}: Report absolute-value permutation p-values instead of equal-tailed ones: the share of replications whose absolute statistic is at or
above the observed one. Intervals become the estimate plus or minus the {it:level()} quantile of the permuted absolute statistics, so p-value and interval agree. The joint,
equality and adversarial tests are unchanged. {it:permabsp} and {it:permside()} may not be combined. {p_end}

{phang}{opt permsi:de(left|right)}: One-sided permutation p-values: {it:permside(left)} is the share of replications at or below the observed statistic and {it:permside(right)} the share at
or above it. The intervals stay the two-sided shift intervals. {p_end}

{phang}{opt permstep:down}: Report, for each family of two or more horizons, the Romano and Wolf max-|studentized| step-down computed from the same
replications. It controls the family-wise error rate across the horizons under the sharp null of no effect at any horizon. Post-treatment and pre-treatment horizons are
separate families. Results are stored in {it:e(perm_stepdown*)} and summarised in the footer. A family whose step-down cannot be formed is silently not posted. {p_end}

{phang}{opt permsa:ve(name)}: Write the permutation draws to {it:name}.dta in the working directory. {it:name} is a plain file name, and an existing file needs
{it:replace}. Rows: one per replication plus a row 0 with the observed estimates. Columns: one per coefficient of {it:e(b)}, the replication standard errors under
{it:permtstat}, and one {it:onset_clustervalue} per cluster holding its onset in each replication. If the file cannot be written, the run still completes with a caution
and posts the return code in {it:e(perm_save_rc)}. {p_end}

{phang}{opt permonset:cohorts}: Assert that {it:i} is the treatment-onset grouping, one level of {it:i} per distinct onset with the never-treated a single level.
On that assertion a {it:cluster()} finer than {it:i} is accepted, and each replication rebuilds {it:i} from the onsets it drew. Onsets after the last period belong to the
never-treated level under a {it:cluster()} finer than {it:i}. At {it:cluster()} equal to {it:i} they keep a level of their own, and nothing is rebuilt. The assertion is
checked and the run stops if it fails. {it:safety(off)} does not lift the check. Compatible with {it:permstrata()} and {it:permute(probabilities)}, not with
{it:subgroups()}. {p_end}

{pstd} {ul: Discussion of Permutation Intervals and Pre-Treatment Tests:} The shift interval widens as the true effect grows and is calibrated only near the
null. {it:permtstat} removes that growth for the post-treatment rows only. The implied standard error in {it:e(perm_se_shift)} is not a sampling standard error. Every
pre-treatment p-value, including the joint pre-trend tests, tests the sharp null of no effect at any horizon. In staggered designs the reassigned onsets move post-treatment
outcomes into the draws' pre-treatment horizons, so the pre-period p-values rise toward 1 as the true effect grows. A large pre-period p-value is then not by itself
evidence of no pre-trend. With a single onset date this does not arise. {p_end}

{marker supporting}{...}
{dlgtab:Supporting Features}

{phang}{opt hist:ogramcohorteffects}: Produce a histogram of cohort i-specific treatment effects, and return those effects in {it:e(histogramestimates)}. Under
{it:att_it} or {it:att_itime}, the column averages each cohort's cell effects. With {it:contreat()}, average marginal effects are histogrammed and appended to the matrix.
Under {it:bootstrap()} and {it:permute()}, the histogram draws the observed point estimates, which neither engine changes, and uses no covariance. {p_end}

{phang}{opt save:plots(name)}: Save every graph as a .gph file, and the data behind any event study plots as .dta files, one per subgroup, all under names
beginning with {it:name}. With {it:makecfxplots} the plotted points are saved as {it:name}_cfx.dta. A suffix on each file name identifies the graph or dataset. {it:name}
must be a Stata name, so it takes no directory part and at most 32 characters. Everything is saved to Stata's working directory. {p_end}

{phang}{opt replace}: Overwrite any existing files named by {it:saveplots()}, {it:permsave()} or {it:bootsave()}. Without it an existing file is an error, raised
under {it:permute()} and {it:bootstrap()} before the replications run. The exception is {it:name}_cfx.dta on a run without {it:bootstrap()} or {it:permute()}, which gets a
caution. {p_end}

{phang}{opt sum:marystats}: Present summary statistics of the outcome and any {it:contreat()} variable: overall, in the treatment group before and after
treatment, and in the control group. {it:contreatwithin} omits the control group rows, and {it:contreatcontrols(untreatedZero)} omits the pre-treatment and control group
rows for the dose. Without a fixed reference period, and with pre-treatment horizons requested, the size of the reference group relative to the rest of the pre-treatment
period is also shown. With event study horizons, the event study composition table {it:e(escomp)} is also displayed, overall and per subgroup. {p_end}

{phang}{opt pois:son}: Estimate by poisson regression, with {it:ppmlhdfe} in place of {it:reghdfe}, or Stata's {it:poisson} in place of {it:regress} under
{it:unconditionalse} or {it:customvce()}. Treatment effects are differences in linear predictions, as from {it:predict(xb)}. It suits an outcome with natural zeros that
one would otherwise log, and groups that share parallel trends only in proportions. Accepts pweights but not aweights. Refused with a negative outcome and with
{it:boottstat}. {it:ppmlhdfe} cannot absorb a slope on a one-level factor or a factor with two or more slopes. So under plain {it:contreatcontrols(basic)}, with no
{it:lattice()} and no other control type, the pooled dose control and, without {it:ccnoabsorb}, any {it:ccc_absorb()} controls enter the regression explicitly. So do the
continuous treatment controls under {it:contreatpoly()} above 1 with a {it:lattice()} or another control type. If {it:ppmlhdfe} still fails with r(3200), as a one-level term in {it:fe()} or a one-level
{it:lattice()} can make it, {it:wooldid} retries once with the {it:ccc_absorb()} and continuous treatment controls moved out of the absorbed effects. {it:fe()} stays
absorbed, so a term there that made the first fit fail stops the run with r(3200). {it:ppmlhdfe} can drop an explicit pooled dose control as collinear with the fixed effects
when the dose's mean is large relative to its spread. {it:wooldid} then prints a caution. Subtracting the mean from the dose avoids it, though under
{it:contreatelasticitytype(dyex)} that also changes the elasticity, which scales by the dose. A treated cohort-period cell whose outcome is zero on every row is perfect
separation. {it:wooldid} names the cell and exits. Coarsen {it:i}, drop the cohort, or specify {it:safety(off)}, which excludes the separated cells from the estimated
effects. The variable name {it:_ppmlhdfe_d} is reserved on poisson runs. {p_end}

{marker eventstudy}{...}
{dlgtab:Event Study Estimation}

{phang}{opt esf:ixedbaseperiod}: Use a fixed reference period rather than a pooled one. The {it:pretreatment} row of the results table then averages every estimable
pre-treatment period, not only the window named by {it:esprelength()}. The joint pre-trend test uses the requested window. {p_end}

{phang}{opt esr:elativeto(-integer)}: The relative time period used as the fixed reference period. The default is -1, the period immediately before treatment. The onset
period or later is an error. Without {it:esfixedbaseperiod} this option has no effect and a caution is printed. {p_end}

{phang}{opt espre:length(integer)}: The number of pre-treatment relative time period-specific effects to be estimated. Under {it:esfixedbaseperiod}, any pre-treatment
request presents every relative period between the {it:esrelativeto()} period and treatment onset. This option and {it:espostlength()} count consecutive integer relative
periods, not realized horizons. On a panel whose relative periods have gaps, such as a biennial panel, the reported window can then be shorter than requested. Widen the
request until it reaches the horizons you want. {p_end}

{phang}{opt espost:length(integer)}: The number of post-treatment relative time period-specific effects to be estimated. With {it:subgroups()} and {it:esprelength()}
greater than 0, a value of 0 is raised to 1 under a caution, since subgroup estimation needs a post-treatment horizon. {it:espostlength(-999)} keeps the post-treatment
window empty. Other negative values are an error. {p_end}

{phang}{opt joint:tests}: Calculate joint significance tests of all requested pre-treatment relative time period-specific effects. At least two pre-treatment horizons
must fall in the estimation window, or three under {it:esfixedbaseperiod} with {it:timetrends} and the default {it:esrelativeto(-1)}. With fewer, no test is produced.
With analytic standard errors clustered at the level of {it:i} or above, the test's p-value is severely oversized, and {it:wooldid} prints the test with a caution that
its p-value is not calibrated (see {help wooldid##standarderrors:Standard Errors}). The remedies are {it:bootstrap()}, {it:permute()}, or a {it:cluster()} strictly finer
than the cohort-by-period cells that leaves many clusters inside each cell. Under a continuous treatment the tests also require many observations per cohort-by-period
cell. Refused with {it:impvalues}, and with {it:fastnose} except on the fast path {it:permute()} takes by default. Under {it:bootstrap()} and {it:permute()} the reported
joint test is a bootstrap or permutation test. {p_end}

{phang}{opt adv:ersarial}: Tests the joint hypothesis that all estimated relative time period-specific effects fall on a single line through the reference period, with
the slope that fits them best. This null is adversarial toward finding pre-trend violations. The p-value comes from a minimum-distance test on the event-study
coefficients and their covariance matrix. Requires {it:esfixedbaseperiod}. The default pooled reference period has no single omitted period, so there the test fails with
a caution and nothing is returned. Under {it:timetrends} two periods are normalized to zero and the null reduces to every coefficient being zero. Beside each p-value,
{it:wooldid} returns the unweighted OLS line through the coefficients, with 0 at the reference period. The test is cancelled with a caution when the window holds fewer
than one pre-treatment and one post-treatment horizon, or fewer than three coefficients in all under {it:timetrends}, when the sample fails to realize those horizons, or
when the event-study covariance is degenerate. Under {it:subgroups()} the {it:att} and {it:att_it} tests are skipped for any subgroup that lost a relative period as
inestimable. A rank-deficient but not degenerate covariance gives a test at reduced rank, and without {it:bootstrap()} or {it:permute()} a caution gives the rank and the
degrees of freedom used. The numerator degrees of freedom are stored in {it:e(adversarial_df_att)} and its twins. Like {it:jointtests}, the test is severely oversized
unless many clusters span each treated cohort-by-period cell, which clustering on {it:i} never gives. There, use {it:permute()}, or {it:bootstrap()} with {it:boottstat}.
Without {it:boottstat} the test is absent under {it:bootstrap()}. Refused with {it:impvalues}, and with {it:fastnose} except under {it:permute()}. {p_end}

{phang}{opt altmethod_adversarial}: An alternative form of the test above. It projects a trend from the pre-treatment coefficients and the reference period, then tests
whether the post-treatment coefficients jointly fall in line with it, allowing for the uncertainty in the projection. It modifies {it:adversarial} and is an error
without it. It requires {it:esfixedbaseperiod} and at least one pre-treatment coefficient beyond the reference period. {p_end}

{phang}{opt pretrendpow:er(numlist)}: Reports Roth (2022) pre-trend power diagnostics for each requested power level, in percent as integers strictly between 5 and 100,
as in {it:pretrendpower(50 80)}. The diagnostics are a GLS estimate of the linear pre-trend slope with its standard error and confidence interval, the smallest linear
pre-trend each pre-test would detect with that power, and the bias such a trend would imply for the pooled estimand if it continued after onset. {it:mdt_wald} is for the
joint pre-trend test and {it:mdt_nis} for the test that any pre-period coefficient is individually significant under the normal reference, both at the 5% level. Requires
{it:esfixedbaseperiod} and at least one pre-treatment event-study coefficient. Refused with {it:impvalues}, {it:fastnose} and {it:permute()}. Under {it:bootstrap()} the
diagnostics use the bootstrap covariance. With {it:subgroups()} only the overall event studies are analysed. When only {it:att_i} or {it:att_itime} is requested the
reported block is the {it:att_it} event study. The diagnostics are skipped with a caution when the pre-period covariance block is unavailable, as under {it:timetrends},
or ill-conditioned. A block of partial rank gives diagnostics at reduced rank under a caution, usually from a singleton-cell design (see condition Two under
{it:e(V_es_type)} in {help wooldid##output:Working with Output}). The {it:mdt_nis} column is computed for up to 4 pre-periods. For larger designs, pass
{it:e(b_es_type)}, {it:e(V_es_type)} and {it:e(numpre_type)} to Roth's {help pretrends} package. Results are returned in {it:e(roth_power_type)} and companion scalars.
{p_end}

{phang}{opt makep:lots}: Produce event study plots. No plot is drawn unless at least one event-study horizon is requested. {it:wooldid} names its graphs, for example
{it:att_fullmodel}, and overwrites its own prior graphs of the same name on each run. It does not clear other graphs from memory. {p_end}

{phang}{opt oos:treatedcontrols(integer)}: By default, cohorts treated after the end of the sample serve as never-treated controls. With a value of T, cohorts treated
within T periods of the end of the sample are handled like other treated cohorts, and can then contribute to pre-treatment relative period effects. Under
{it:esfixedbaseperiod} such cohorts are dropped if their reference period is not usable. With {it:anticipation()} the window applies to the shifted onset. Negative
values are not allowed. {p_end}

{phang}{opt antic:ipation(integer)}: Treats the effective treatment onset as {it:ttre} minus #, so periods of anticipatory response count as treated rather than as clean
pre-treatment periods. This equals re-running {it:wooldid} with {it:ttre} minus # in place of {it:ttre}, as {it:did_imputation}'s {it:shift()} option does. Relative
period 0 is then # periods before the recorded onset, so you may wish to raise {it:espostlength()} by # to keep the same number of genuine post-treatment horizons. The
value is in the units of {it:t}. On a biennial panel an odd value lands between observed periods and {it:wooldid} prints a caution. Units whose effective onset falls at
or before their first usable period are dropped with a caution, as are, under {it:esfixedbaseperiod}, units whose shifted reference period is unusable. Too small a value
leaves anticipation in the reference period, which shows in the pre-treatment horizons. Too large a value keeps the pre-trends clean but attenuates the pooled estimates,
and the first post-treatment horizons read near zero. A run without an event study shows neither, so request one and prefer the smallest defensible value. With
{it:subgroups(ttre)} the subgroup labels remain the recorded onset values. The shift is recorded in {it:e(anticipation)}. Negative values are not allowed. {p_end}

{pstd}{ul: Two Approaches to Estimation, Pooled versus Fixed Reference Period:} By default, {it:wooldid} compares each cohort's treated periods to a pooled reference
period holding all of that cohort's pre-treatment observations. Event study estimates on this default are normalized to 0 in the earliest period reported, the pooled
reference period. Under {it:esfixedbaseperiod}, every estimate compares treated cohorts' post-treatment outcomes to their outcome in one pre-treatment reference period,
relative period -1 unless {it:esrelativeto()} sets another. {p_end}

{pstd} The choice also changes the comparison group. Under the pooled reference period every untreated observation serves as a control: never-treated cohorts, cohorts
treated after the end of the sample, and the pre-treatment periods of cohorts not yet treated, apart from any pre-treatment horizons the run asks for. Parallel trends
and no anticipation must then hold across all of those periods. Under {it:esfixedbaseperiod} each treated cohort's pre-treatment periods other than the reference period
carry their own indicators, and the reference period only anchors the cohort's own level. With {it:timetrends}, the period just before the reference period carries no
indicator either, and the two anchor the cohort's level and trend. The period effects then come from never-treated cohorts and cohorts treated after the end of the
sample alone. Cohorts not yet treated do not serve as controls. Parallel trends is needed only between each treated cohort and those cohorts, from the reference period
on. This is the never-treated comparison of Callaway and Sant'Anna (2021) and Sun and Abraham (2021). The pooled default matches the imputation estimator of Borusyak,
Jaravel, and Spiess (2023) under the conditions given in {help wooldid##comparison:Comparison with Related Estimators}. {p_end}

{pstd} A fixed reference period is transparent: every cohort is compared to the same relative period. A pooled one can be more precise, since it averages over more
periods, but different cohorts' reference periods then hold different numbers of periods. A pooled period is likely a poor choice when treated cohorts' outcomes evolve
differently from the comparison cohorts' before treatment. Trends shared by all cohorts are absorbed by the time fixed effects. It also tends to be a poor choice when
requesting more pre-treatment horizons than every cohort can supply. A cohort with fewer observed pre-treatment periods then takes its earliest observed period as its
reference period, even where other cohorts report an effect for that period. With differing outcome trends, the event study can then show counterintuitive time paths
that hide even severe parallel trends violations. {p_end}

{pstd}{ul: Invalid Estimands:} For relative-period estimates the estimand types {it:att_i}, {it:att_itime} and {it:att_it} coincide. Each weights every cohort's
contribution to a relative period equally, and they differ only in cells outside that period. Event studies therefore report only {it:att} and {it:att_it}, and a request
for {it:att_i} or {it:att_itime} shows the {it:att_it} results. {p_end}

{marker interactions}{...}
{dlgtab:Subgroup Effects}

{phang}{opt sub:groups()}: Specifies the subgroup (SG) variable, whose levels define the subgroups for which {it:wooldid} computes separate estimates. The SG variable
must hold non-negative integers or missing values on the observations of treated units, and its levels must be collections of cohorts i. Observations with missing values
belong to no subgroup and are ignored in subgroup computations. A unit is left out only if it is missing on all of its rows. A treated unit split across subgroups, or
labelled on some rows and missing on others, is refused. With {it:safety(off)}, except under {it:permute()}, such a unit enters each of its subgroups weighted by its
labelled share of sample weight there, under a caution. A split unit then no longer carries equal weight in {it:att_i} and {it:att_itime}. Subgroup values are embedded
in the names of stored results and plots, so values 8 or more digits long can exceed Stata's limits on names. Recode ID-style codes to small integers first, for example
with {it:egen newsg = group(oldsg)}. With {it:esprelength()} greater than 0, an {it:espostlength()} of 0 is raised to 1 under a caution (see {it:espostlength()}). Not
compatible with {it:impvalues}. {p_end}

{pstd} {ul: Description of Approach to Subgroup Effect Estimation:} {it:wooldid} estimates treatment effects for each level of the SG variable by running the margins
step within that subgroup. The regression itself is unchanged. It also reports every pairwise difference between the subgroups' overall and pre-treatment effects, each
with a p-value, and a joint test that all SGs share the effect of the lowest-numbered SG. Clustered standard errors for a subgroup estimate can fail when the subgroup
holds few clusters, or only one. With clustering at the level of {it:i} or above, the analytic standard error of a pairwise contrast can be degenerate. {it:wooldid} then
reports that contrast, and the across-subgroup equality test covering it, without inference, under a caution naming the count. With {it:safety(off)} they are reported as
computed. Under {it:permute()} the contrasts are permuted as coefficient differences and stand, and the caution counts the rows whose {it:se_analytic} is degenerate.
Under {it:boottstat} each contrast is studentized by that standard error, and a caution says so when it is degenerate. The remedies are {it:permute()}, {it:bootstrap()}
without {it:boottstat}, or a {it:cluster()} strictly finer than the cohort-by-period cells. {p_end}

{pstd} A subgroup variable that varies within a cohort-by-period cell, or over time within a unit, is refused with an informative error. Where it varies within a cell
and {it:safety(off)} lets the run proceed, every subgroup estimate only reweights the same pooled cell effects. Contrasts then cannot show within-cell differences, and
can come back identical when the true effects differ. The usual remedy is to define {it:i} at the finer level over which the subgroup varies, for example
{it:egen newi = group(i sg), missing}. Without {it:missing}, rows with a missing subgroup value would be dropped from the estimation sample. Splitting {it:i} this way
makes {it:att_i} and {it:att_itime} average over the split pseudo-units rather than over the original units. {p_end}

{pstd} A triple difference is a special case. Its target is the difference between the diff-in-diff estimates of two groups, each present among both treated and control
cohorts. Define {it:i} as the treatment onset cohort by group, for example {it:egen newi = group(ttre sgvarname), missing}. Specify the group variable as the SG variable
and add SG by period fixed effects, as in {it:subgroups(sgvarname) fe(i.sgvarname##i.t)}. The contrast removes shocks that hit both groups in a treated cohort and period
only when the two SG estimates weight the same onset-cohort-by-period cells in the same way. Use {it:att_it}, with both groups present in every treated
onset-cohort-by-period cell. Its contrast is then the simple average of the cell-level triple differences. Under {it:att} each group weights its cells by the sample
weight of its own observations, so that contrast keeps those shocks whenever the group shares of that weight differ across cells. {p_end}

{marker continuous}{...}
{dlgtab:Continuous Treatments}

{phang}{opt contreat(varname)}: The continuous treatment variable. It must be defined for all observations, set to 0 where it has no natural value in untreated
cohort-periods. Observations on which it is missing are dropped from the estimation sample. If no treated cohort-period cell has dose variation that identifies its dose
slope, for example with a dose constant within each cohort, the within-cell average marginal effects are not identified. Under default safety {it:wooldid} says so,
reports the ame results as missing and continues with the att results. Under {it:contreatcontrols(untreatedZero)} that design is refused instead. If only some treated
cells lack within-cell dose variation, their dose slopes are not identified either and would enter the ame as zeros. {it:wooldid} drops those cells from the estimation
sample with a caution, so the att and the ame are estimated on the remaining cells. The cells checked are those that carry a dose slope, which include the pre-treatment
periods under {it:esfixedbaseperiod}. The number dropped is stored in {it:e(contreat_ndropcells)}, and their observations are counted in {it:e(N_dropped)}.
{it:safety(off)} silences the caution but keeps the drop. A cell whose dose is zero throughout is never dropped. No cell is dropped when no cell on its side of treatment
onset has dose variation. Under {it:contreatcontrols(untreatedZero)} the check does not run. Not available with {it:gteffects}. {p_end}

{phang}{opt contreatelas:ticitytype(dydx|dyex)}: The type of marginal effect of the continuous treatment. {it:dydx}, the default, reports the change in y per unit change
in the dose. {it:dyex} reports the change in y per proportional change in the dose. Under {it:poisson} results are on the linear index, where {it:dydx} already gives the
semi-elasticity of E[y] and {it:dyex} the elasticity, and the tables and plots are titled that way. Only {it:dydx} and {it:dyex} are accepted. {p_end}

{phang}{opt contreatw:ithin}: Use only within-treated-cohort variation in the continuous treatment. It lifts the checks that stop estimation when there are no untreated
cohorts, and stops {it:wooldid} objecting to i-t fixed effects in {it:fe()}. It also skips the sample cleaning that drops time periods with no control cohort. The
data-integrity checks and the dose-slope check under {it:contreat()} still run unless {it:safety(off)} is given. Requires {it:contreat()}. {p_end}

{phang}{opt contreatp:oly(integer)}: Includes polynomial terms in {it:contreat} up to the given degree, a positive integer. The default is 1. Be cautious using higher
degree polynomials. {p_end}

{phang}{opt lattice(latticevar OR latticetype latticenum)}: Specifies a lattice variable, whose groups of observations let the effect of {it:contreat} vary, even within
cohort-periods and within the reference period. Given a variable name, that variable is the lattice. It must be numeric and contain only non-negative integers, and
observations on which it is missing are dropped from the estimation sample. Given a lattice type and a lattice number, {it:wooldid} builds the lattice from
{it:latticenum} quantiles, an integer greater than 1, formed within the groups the {it:latticetype} names (see below). Sample weights are used in building it. Without
{it:lattice()} the whole sample is a single lattice bin. A lattice is the most expensive option in runtime terms (see {help wooldid##runtime:Runtime}). Requires
{it:contreat()}. When a lattice is combined with {it:contreatcontrols(basic)}, the baseline dose level is partialled out within each lattice level. {p_end}

{phang} {it:latticetype:} The groups within which the {it:latticenum} quantiles of {it:contreat} are formed. Type names are case-sensitive. The choice matters for
interpretation, as of {it:makecfxplots} output, and for feasibility: clustered standard errors often fail unless multiple clusters are represented at each lattice point.
The types are: {p_end}
{synoptset 20 tabbed}{...}
{synopt : {it:all}} Quantiles are formed over the whole within-sample distribution of the continuous treatment variable, untreated observations included. {p_end}
{synopt : {it:byTime}} Quantiles are formed within each time period t separately. An observation is in the lowest bin if its {it:contreat} is in the lowest
1/{it:latticenum} of its period's distribution, even if its value is high in the sample overall. {p_end}
{synopt : {it:byTreated}} Quantiles are formed separately in the treated and untreated groups. {p_end}
{synopt : {it:byTtre}} Quantiles are formed separately for the control group and for each treatment onset cohort. {p_end}
{synopt : {it:byCohort}} Quantiles are formed separately within each cohort i. {p_end}
{synopt : {it:byTreatedByTime}} Quantiles are formed by time period t, separately for the treated and untreated groups.  {p_end}
{synopt : {it:byTtreByTime}} Quantiles are formed by time period t, separately for the control group and each treatment onset cohort. {p_end}
{synopt : {it:byCohortByTime}} Quantiles are formed separately by cohort and time period t. {p_end}
{synopt : {it:baseTreated}} Quantiles are formed within the treated group, and their cut points define the lattice points for the rest of the sample. {p_end}
{synopt : {it:baseTreatedRef}} As {it:baseTreated}, with the quantiles formed within the treated group's reference period. {p_end}
{synopt : {it:baseTreatedPost}} As {it:baseTreated}, with the quantiles formed within the treatment group's treated period. {p_end}
{synopt : {it:baseTreatedPre}} As {it:baseTreated}, with the quantiles formed within the treatment group's full pre-treatment period. {p_end}

{phang}{opt latticeig:noreweights}: Ignore sample weights when forming the lattice quantiles. Under {it:bootstrap()} each replication's resampling multiplicity is kept,
so its lattice is formed on the resample. {p_end}

{phang}{opt makec:fxplots}: Plot the average marginal effect of {it:contreat} and the average treatment effect within each bin of the lattice variable or, if given, the
{it:cfxlattice()} variable. One of the two must be present. The plots are produced for {it:att} only unless {it:safety(off)} or {it:contreatwithin} is specified. They
are sure to be valid only for {it:att}. For other estimands they may be valid if the bins do not cut across cohorts, or, for {it:att_it}, across cohort-periods. Requires
{it:contreat()}. Not compatible with {it:impvalues}. {p_end}

{phang}{opt cfxlattice(varname)}: The variable whose bins the continuous treatment effects plots use. The default is the lattice variable. It must be numeric and contain
only non-negative integers on treated observations. Where it is missing on treated observations, {it:margins} omits those observations from the plotted averages and
{it:wooldid} prints a caution. Without {it:makecfxplots} the option has no effect and draws a caution. {p_end}

{phang}{opt cfxplott:ypes(pre post integers)}: The cfxplots to produce. {it:post} plots effects across the continuous treatment within the post-treatment period,
{it:pre} within the pre-treatment period, and an integer within that relative period. The default is {it:post}. {it:pre} requires {it:esfixedbaseperiod} or
{it:esprelength()} greater than 0. Without {it:esfixedbaseperiod} an integer must lie within the requested event-study window, and entries outside it are dropped without
a message. With {it:esfixedbaseperiod} any integer is accepted except the reference period itself, which is dropped. Any other entry is an error. {p_end}

{phang}{opt makecfxplotsbysg}: Produce the continuous treatment effects plots for each subgroup. Without both {it:makecfxplots} and {it:subgroups()} it has no effect and
prints nothing. Each subgroup's plot data is returned in {it:e(cfxitem_type_sgN)}, where N is the subgroup value. {p_end}

{phang}{opt cfxplotnose}: Do not estimate standard errors for the cfxplots, which speeds estimation. Dropped, with a caution, under {it:permute()} with {it:permtstat} or
{it:permtci} and under {it:bootstrap()} with {it:boottstat}, which need them. Otherwise the plot cells keep their replication inference under both engines. {p_end}

{phang}{opt contreatc:ontrols(controltype)}: The controls for {it:contreat}. Each is f(Contreat_jt), the function that captures the dose effect within treated i-t cells,
interacted with a set of fixed-effect-like indicators such as time period or cohort. They set the reference dose effect against which the effect in treated
cohort-periods is measured. This option is required with a continuous treatment. Choose it with care, since the types make very different substantive assumptions (see
below). {p_end}

{phang} {it:conttreatcontroltype}: {it:untreatedZero} alone, or one or more of {it:basic}, {it:byTreated}, {it:byTtre}, {it:byCohort} and {it:byTime}. {p_end}
{synopt : {it:basic}} f(Contreat_jt) is added without interactions. {p_end}
{synopt : {it:byTreated}} Controls for the effect of contreat are interacted with indicators for the treatment group and the control group. {p_end}
{synopt : {it:byTtre}} Controls for the effect of contreat are interacted with indicators for each treatment onset group and the control group. {p_end}
{synopt : {it:byCohort}} Controls for the effect of contreat are interacted with cohort (i) indicators. {p_end}
{synopt : {it:byTime}} Controls for the effect of contreat are interacted with time period (t) indicators. {p_end}
{synopt : {it:untreatedZero}} For a continuous treatment that is always 0 outside treated cohort-periods, which suits a dose that is unknowable or undefined there. If
the dose is not always 0 outside treated cohort-periods, the baseline effect of {it:contreat} is left out of the model, and {it:wooldid} will return biased estimates,
however precise they appear. {p_end}

{pstd} Each control type needs variation in {it:contreat} across the reference-period observations of every group for which it fits a baseline dose slope. That is the
whole sample under {it:basic}, each treatment-status group under {it:byTreated}, each onset group under {it:byTtre}, each cohort under {it:byCohort} and each calendar
period under {it:byTime}. If any group lacks it the run stops. Under {it:byTime} this happens whenever never-treated units carry a zero or constant dose in the periods
after the last onset. A group whose variation is below 1/100000 of the sample standard deviation of {it:contreat} is refused as well. One between that and 1/1000 of it
draws a caution instead, since its baseline dose slope is nearly unidentified and the estimates can be wrong rather than noisy. {p_end}

{pstd} {ul:On Specifying contreatcontrols() and lattice:} With a continuous treatment the regression is:
Y_jt = fe_i + fe_t + sum_over_it_in_Z(beta_it * TreatedFlag_it * Tz_it + TreatedFlag_it * Tz_it * f(Contreat_jt)) + g(f(Contreat_jt))
Each treated i-t cell gets its own treatment effect intercept and its own effect of the continuous treatment through f(), which can be made flexible where that effect may be non-linear within a cell. The options {it:contreatpoly()} and {it:lattice()} set f(), in general a stepwise function of contreat. Without a lattice, f(Contreat_jt) = beta_it1*Contreat_jt^1 + beta_it2*Contreat_jt^2 + ... + beta_itp*Contreat_jt^p, where p is set by {it:contreatpoly(p)}. With a lattice, f(Contreat_jt) = sum_over_l_in_L[L_l + L_l * Poly(Contreat_jt)], where L_l is 1 for observations in lattice bin l and 0 otherwise, and Poly(Contreat_jt) holds the requested polynomial terms of Contreat_jt. {p_end}

{pstd} A lattice lets the dose slope vary within i-t cells, and it also adds treatment effect intercepts within cells, much as a finer {it:i} would. Average marginal
effects do not capture those intercept jumps, nor the cell intercepts themselves. With many lattice groups the intercepts can absorb nearly all of the dose effect and
threaten identification of the ame. The ame and the att can also disagree, even in sign. In such cases use a coarser lattice or rely on the att, which does capture the
intercepts, above all in the continuous treatment effects plots. {p_end}

{pstd} Do not over-saturate the specification with lattice terms. Each treated cohort-period cell must hold at least two distinct values of Contreat_jt within each
lattice bin it occupies. That variation is all that identifies the piece's slope. Without it, one slope term in that cell drops for collinearity and the piece takes the
slope of another bin in the same cell. The reported average marginal effect and the plotted average marginal effects by bin include that borrowed slope. A caution gives
the number of such pieces, and {it:verbose} lists the dropped terms. Separately, a piece with no more observations than parameters is fitted exactly, and the analytic
standard errors can then be far smaller than the estimates' sampling variability. As a rule of thumb, the units contributing to a treated cell should number several
times the lattice bins. {p_end}

{pstd} The control function g(f(Contreat_jt)) plays for the dose the role the fixed effects play for a binary treatment. It largely decides which variation identifies
the dose effect, since it sets the reference dose effect against which the effects in treated cells are measured. Each control type adds f(Contreat_jt) interacted with
the indicators listed under {it:contreatcontrols()}. For example, {it:basic} gives g(f(Contreat_jt)) = f(Contreat_jt), and {it:byCohort} gives g(f(Contreat_jt)) =
sum_over_i[fe_i * f(Contreat_jt)]. {p_end}

{pstd} {it:byTime} partials out period-specific effects of the dose, much as the second difference does for a binary treatment. It can be specified alone or with another
control type. With {it:byCohort} and {it:byTime}, g() = sum_over_i[fe_i*f(Contreat_jt)] + sum_over_t[fe_t*f(Contreat_jt)]. {p_end}

{pstd} With {it:untreatedZero}, g() = 0, so the dose slopes are identified only from variation within treated cohort-periods. Use it only for a continuous treatment that
is not meaningfully defined outside treated cohort-periods. {p_end}

{pstd} To identify the average marginal effects from dose variation within treated cohorts, without a strictly untreated control group, specify {it:contreatwithin}.
{it:wooldid} still reports att estimates, but read them with care. In a period where some unit is not yet treated, a treated cell's att includes the level effect of
treatment, measured against that period's untreated units. In a period with no untreated reference rows, the treated cells' level effects are not identified.
{it:wooldid} then sets one cell's level effect in that period to zero, measures the other cells against it, and names those periods in a caution. The att there depends
on which cell was set to zero. It can change when the cohorts are relabelled, and it can take the wrong sign. Rows moved into an {it:esprelength()} window also leave the
reference pool, so the pretreatment att can be affected the same way. The att and its event study can therefore shift sharply at the last treatment onset for this reason
alone. The att is the integral of the estimated average marginal effects across the dose distribution only when {it:fe()} absorbs every cohort-period level effect, for
example with i-t fixed effects. Even then, it extrapolates the estimated slope down to a dose of zero. Rely on the ame results under {it:contreatwithin}. The control
type {it:byTime} is refused here, since without an untreated group the time period interactions would partial out genuine treatment effects for some group.
{it:safety(off)} downgrades the refusal to a caution. {p_end}

{pstd} {ul:What the ame Measures:} Within a treated cell, the ame compares units that received different doses. The choice of g() sets the parallel trends assumption:
without treatment, units at different doses would have followed the same outcome path, apart from the dose effects g() allows for. That identifies the att, and the att
within each bin of the continuous treatment effects plots. The ame is the effect of a marginal increase in the dose only under a stronger condition: units at different
doses would also have gained the same from treatment at any given dose. This is the strong parallel trends assumption of Callaway, Goodman-Bacon and Sant'Anna (2024).
Without it, the ame also carries selection bias, since the units that took higher doses may simply be those that gain more from treatment. The same holds for the ame
within each bin and for differences in att across bins. Where the att under {it:contreatwithin} is the integral of the ame, it needs the same stronger condition. {p_end}

{pstd} {ul:Continuous Treatments and Cluster-Robust Standard Errors:} With {it:cluster()} at the level of {it:i} or coarser, the analytic standard errors of the
ame-family quantities are not identified. {it:wooldid} prints a caution, reports those standard errors, the ame joint tests and the ame adversarial p-values as missing,
withholds {it:e(V_es_ame*)}, and draws the ame continuous treatment effects plots without confidence intervals. The att rows are still reported, but their analytic
standard errors are understated in this design as well (see {help wooldid##standarderrors:What Level to Cluster At}). A cell holding a single observation does not
trigger this. Valid ame inference then comes from {it:bootstrap()}, {it:permute()}, {it:impvalues}, {it:robust} when errors are independent across observations, or a
{it:cluster()} strictly finer than the cohort-by-period cells. {p_end}

{marker controls}{...}
{dlgtab:Controls and Fixed Effects}

{phang}{opt control:s()}: Passed unmodified as controls to the regression command: {it:reghdfe}, {it:ppmlhdfe} under {it:poisson}, and {it:regress} or
{it:poisson} under {it:customvce()} or {it:unconditionalse}. {p_end}

{phang}{opt fe()}: Passed unmodified to {it:reghdfe}'s absorb statement as fixed effects. Variables are read as categorical, and i.X##c.X2 notation is
permitted for heterogeneous slopes and the like. Under {it:customvce()} or {it:unconditionalse} estimation runs through {it:regress}, where a bare variable name is read as
continuous. Write the {it:i.} prefix explicitly there. {p_end}

{phang}{opt timetrend:s}: Control for cohort-specific time trends, i.i##c.t. This option is required for heterogeneous time trend controls to work. It also
updates {it:wooldid}'s estimability checks, drops cohorts without enough pre-treatment data to estimate a trend, and, with a fixed reference period, makes the two
reference periods normalized to 0 chronologically adjacent. It adds no trends in a continuous treatment variable. With {it:cluster()} on {it:i}, the cohort trend terms
typically leave the covariance unidentified. {it:wooldid} says so and reports the estimates without standard errors. The same specification under {it:robust} is
unaffected. {p_end}

{phang}{opt coarse:cohortcontrols()}: Includes the specified controls and their interactions with every Tz_it dummy in the model, one for each cohort-time
period with an estimated treatment effect. With {it:esfixedbaseperiod} that is every estimable treated cohort by relative time period. The control slopes are pooled across
the reference period and the control group, and each treated cohort-period gets its own slope. Each reported effect is the cell's average outcome less its predicted
untreated outcome at the pooled slope, so it does not depend on where the control is centred or on its units. Cohort by time period interactions you specify yourself
differ, because they do not pool the reference period and the control group. Every coarse cohort control must use factor notation, {it:c.} for continuous and {it:i.} for
categorical. {p_end}

{phang}{opt inter:activecontrols()}: Includes the specified variables in two capacities. They enter as controls interacted with time period indicators. They
are also de-meaned by cohort {it:i} and interacted with every Tz_it dummy, so treatment effects vary linearly with them. Unlike {it:coarsecohortcontrols()}, that variation
is reported as part of the treatment effects rather than partialled out. This most closely mimics the approach to covariates in Wooldridge (2025). The treatment-effect
slopes exist only where a control varies within cohort-by-period cells. Where it does not, for example a unit-level control with a unit-level {it:i}, the run is identical
to {it:controls(i.t#c.x)} and {it:wooldid} prints a caution naming the control. This option also turns the default cell collapse off. See {it:nofastcollapse}. {p_end}

{pstd} {ul:Caution Regarding Continuous Treatment Variables and Controls:} Controls you build from a continuous treatment variable must not include that variable
itself, which can distort {it:wooldid}'s margins calculation. A control with identical values under a different name is fine. Pair such controls with
{it:contreatcontrols(basic)}, or any control type other than {it:untreatedZero}. {p_end}

{marker technical}{...}
{dlgtab:Technical Features}

{phang}{opt update}: If specified, the command {it:wooldid,} {it:update} will update your local installation of {it:wooldid} based on the most recent edition of the program available on github.
A running Stata session keeps using the copy it has already loaded, so the update takes effect after you restart Stata or run {cmd:discard}. {p_end}

{phang}{opt customw:eightmultiplier(varname)}: Weights the observations in the margins calculation, multiplying any sample weight, to produce bespoke
treatment effects. It scales the margins weights only, not the regression weights. See step 2 under {help wooldid##description:Description}. Negative values are an error
unless {it:safety(off)} is specified. Missing values drop the observation. {p_end}

{phang}{opt regtol(#)}: The convergence tolerance for {it:reghdfe} and {it:ppmlhdfe}, as the negative base-ten exponent: {it:regtol(9)} means a tolerance of
1e-9. The {it:wooldid} default is 9, where {it:reghdfe}'s own is 8. A higher value may help when many interaction terms with continuous variables are absorbed, as with
continuous treatments. It has no effect under {it:unconditionalse} or {it:customvce()}, which estimate through {it:regress} or {it:poisson}. {p_end}

{phang}{opt ccno:absorb}: Enter the continuous treatment controls in the main regression rather than in the {it:absorb()} option of {it:reghdfe} or
{it:ppmlhdfe}, where they go by default. Some {it:poisson} runs already enter them explicitly. See {it:poisson}. Under certain conditions this improves estimation stability
and keeps variables from being dropped wrongly. {p_end}

{phang}{opt ccc_absorb()}: A version of {it:coarsecohortcontrols()} for continuous variables. The pooled control slopes enter the regression. The treated
cells' own slopes, on the control centred at each cell's mean, are passed to the {it:absorb()} statement. The reported effects are those of {it:coarsecohortcontrols()}.
It is often slower but can save memory in the main regression. Be extra careful that specifications with {it:ccc_absorb()} are identified, since not all identification
failures are detected. Under {it:contreat()}, naming the treatment variable itself here turns the default cell collapse off. See {it:nofastcollapse}. On some
{it:poisson} runs these slopes enter the regression explicitly. See {it:poisson}. {p_end}

{phang}{opt ver:bose}: Print the regression command {it:wooldid} uses and progress reports throughout. The regression itself runs quietly. The underlying
regression is also stored as {it:___wooldidfullmodel}. {it:estimates restore} can recover it for coefficient-level work only, such as {it:lincom} and {it:test}. {p_end}

{phang}{opt gte:ffects}: Returns in {it:e(gteffects)} the cohort-by-relative-period effect estimates underlying {it:wooldid}'s aggregates, one row per estimated
{it:i}-by-relative-period cell, read directly off the underlying regression. The table is printed when it has at most 200 rows. Columns: {it:i}, {it:onset}, {it:caltime},
{it:relyr}, {it:estimate}, {it:se}, {it:p}, the interval bounds named for the run's {it:level()}, {it:inwindow}, which is 1 inside the reported event-study window,
{it:nobs}, and the cell's total estimation weight under each requested estimand: {it:w_att}, {it:w_att_it}, {it:w_att_i} and {it:w_att_itime}. The table supports custom
aggregations. The weighted mean of {it:estimate} under a weight column, over the post-treatment rows ({it:relyr} at or above 0) or the pre-treatment rows, reproduces the
corresponding main or pretreatment row exactly. Over a single {it:relyr} it reproduces the event-study row under {it:w_att} or {it:w_att_it} only, since {it:att_i} and {it:att_itime}
report the {it:att_it} event study. With {it:subgroups()} one full-sample table is returned. Standard errors are the regression's coefficient-level ones. A cell's standard
error is understated, and a caution reports it, when the cell lies inside one cluster, as with {it:cluster()} at the level of {it:i} or coarser. The same holds for a cell
with a single observation under {it:robust}, pweights, a robust {it:customvce()} or {it:poisson} on {it:ppmlhdfe}. A cell whose variance is not identified has a missing
standard error. Not available with {it:contreat()}, where the cell coefficients are dose-zero intercepts rather than treatment effects. Available under {it:fastnose}. Under
{it:bootstrap()} the cell inference is bootstrap-based, but a cell the caution counts is resampled whole, so its standard error stays understated. Under {it:permute()} the
cells carry no inference. {p_end}

{phang}{opt safety(on|off)}: By default {it:wooldid} cautions or halts when it is likely to yield unreliable estimates. The main check is whether an i-t cell's treatment
variable was dropped for multicollinearity, as when fixed effects or continuous treatment controls partial out a main effect. Such failures are generally detected, but not always for a continuous treatment with a lattice or polynomial terms. Default safety also checks that onset is
constant within {it:i}, that weights and the {it:customweightmultiplier()} variable are non-negative, and that lattice variables are non-negative integers. It cautions when
{it:i} takes more than 100 distinct values. See {help wooldid##runtime:Runtime}. It halts with an informative error before estimation if the specification would exceed
the current {help maxvar} setting. {it:safety(off)} skips these checks, and downgrades a number of refusals described elsewhere in this help to cautions. If a run changes
your specification, {it:wooldid} says so whatever {it:safety()} is set to. Only {it:on} and {it:off} are accepted. {p_end}

{marker runtime}{...}
{title:Runtime}

{pstd} Runtime is driven by the number of cohort-by-period interaction cells {it:wooldid} builds, roughly the number of distinct values of {it:i} times the number
of relative periods in the event-study window. It grows faster than linearly in that count, while sample size, subgroup count and event-study width each cost only
linearly. Coarsening {it:i} is therefore normally the single largest speedup: define {it:i} as cohorts or treatment onset groups and put unit fixed effects in {it:fe()}.
Dropping observations buys little, so subsampling is a poor way to shorten a long {it:bootstrap()} or {it:permute()} run. Shrinking the window or coarsening {it:i} works
there. With {it:lattice()}, prefer the smallest {it:latticenum} that answers the question. {p_end}

{pstd} Nearly all of a run is spent inside {help margins}, most of it in the delta-method variance calculation. {it:fastnose} skips that calculation and saves
several times over when standard errors are not needed. The default cell collapse shrinks the data each {it:margins} call sees, and the two compose. Runs the collapse
skips, listed under {it:nofastcollapse}, pay the full per-row cost. {p_end}

{marker output}{...}
{title:Working With Output}

{pstd}
{it:wooldid} stores its results as matrices and scalars in {it:e()}, where {it:esttab} and similar tools find them. {it:ereturn list} shows the visible results.
Hidden results, such as the preserved analytic covariances, the per-family bootstrap and permutation tables and the {it:e(wdid_*)} records, appear only under
{it:ereturn list, all}. Each can be read by name either way. A test the output reports as absent or withheld may be missing from {it:e()} or stored there as missing,
for example the adversarial test under {it:bootstrap()} without {it:boottstat}. Check it with {it:missing()}, which is true either way, not with {it:confirm scalar}.
Typing {it:wooldid} alone after estimation, or after {it:estimates restore}, redisplays the stored results. It does not redraw plots or repeat estimation-time cautions.
A call that fails, at any stage, leaves {it:e()} as it was, so typing {it:wooldid} alone then redisplays the earlier results. {p_end}

{pstd} Stata's {help bootstrap} and {help jackknife} prefixes work with {it:_b}, as in {it:bootstrap _b, reps(200) cluster(i): wooldid y i t ttre, att cluster(i)}.
The point estimates are those of the plain run. A {it:cluster()} given to {it:wooldid} must name the prefix's {it:cluster()} variable. The prefix refuses any other,
including its own {it:idcluster()} variable. The {it:bootstrap} prefix does not allow weights, and {it:svy} is not supported. After a prefix, typing {it:wooldid} alone
shows the prefix's table. The {it:bootstrap()} option is the supported route to bootstrap inference. {p_end}

{pstd} {it:esttab}, {it:coefplot}, {it:lincom} and {it:test} address a coefficient of {it:e(b)} and {it:e(V)} by these column names.
{it:main_att} is the headline average treatment effect on the treated, and {it:pretreatment_att} the pooled pre-treatment average. {it:ES_att_pre#} and {it:ES_att_post#}
are the event-study horizons, named by relative period: {it:pre#} is relative period -# and {it:post#} is relative period #, so {it:post0} is the onset period.
The reference period has no column. The other estimand types replace {it:att} with {it:att_it}, {it:att_i}, {it:att_itime} or {it:ame}. Under {it:subgroups()} each
subgroup's rows carry the infix {it:_SG{it:sgcase}_}, where {it:sgcase} is the subgroup value, as in {it:main_SG1_att}, {it:pretreatment_SG1_att} and
{it:ES_att_SG1_post0}. The suffix {it:_s{it:sgcase}} names per-subgroup {it:e()} matrices such as {it:e(b_es_att_s1)}, not these columns. {p_end}

{pstd} The primary results are stored in the following matrices: {p_end}
{synoptset 30 tabbed}{...}
{synopt:{it:e(wdidmainresults)}} All primary and event study results. The columns are {it:estimate}, {it:se}, {it:t} ({it:z} under {it:poisson}), {it:p}, the lower and upper interval
bounds (named {it:lb_95ci} and {it:ub_95ci}, or for the run's {it:level()}), and {it:relyr}. The pooled estimands sit under the equations {it:main} and {it:pretreatment}, and the
event-study rows under {it:ES_att}, {it:ES_att_it}, {it:ES_ame} and {it:ES_ame_it}, named {it:post#} and {it:pre#}.{p_end}

{synopt:{it:e(wdidsgresults_{it:sgcase})}} Primary and event study results for the subgroup {it:sgcase}, if subgroups are used.{p_end}

{synopt:{it:e(pwsgs)}} Pairwise differences between the subgroup-specific treatment effects, if subgroups are used, with the same columns less {it:relyr}. A row whose analytic
contrast covariance is degenerate carries its point estimate with the inference columns missing (see {help wooldid##interactions:Subgroup Effects}). {p_end}

{synopt:{it:e(b_es_{it:type})}, {it:e(V_es_{it:type})}} The event study coefficient vector and its full covariance matrix, labelled {it:pre#} and {it:post#}
with the reference period omitted, for each estimand type: {it:att}, {it:ame} under {it:contreat()}, and {it:att_it} or {it:ame_it} for the cell-weighted estimands, which
share one horizon block. The per-subgroup versions are {it:e(b_es_{it:type}_s{it:sgcase})} and {it:e(V_es_{it:type}_s{it:sgcase})}. They are valid for {it:lincom} and
{it:test} once posted, as shown under Incompatible Post-Estimation Commands below. A direct positional handoff to {it:honestdid} also requires {it:esfixedbaseperiod} with
the default {it:esrelativeto(-1)}, no {it:timetrends}, at least one pre- and one post-treatment coefficient, consecutive horizons, and a covariance with no degenerate rows.
Any other layout hands {it:honestdid} the wrong event times without warning. The covariance is withheld, under a caution, when no horizon has an identified variance,
typically under {it:timetrends}, and for the ame families under the cell-containment caution. The vector still posts. Under {it:fastnose} the covariance is never
computed, so it is not posted and no caution prints. Under {it:impvalues} neither matrix is posted. {it:e(es_ndegen_{it:type})} counts horizons whose variance is
degenerate and is posted only when there are any. Two design conditions shape the analytic matrix without showing in it, and {it:honestdid} and {it:pretrendpower()} inherit both. One:
a treated cohort-by-period cell wholly inside one cluster adds none of its own variation, see {help wooldid##standarderrors:Standard Errors}. A {it:cluster()} strictly finer than the cells avoids this. Two: with a unit-level
{it:i} each treated cell is a single observation with a zero residual, so the covariance carries control-side variation only, under {it:robust} as well. A coarser {it:i}
with unit fixed effects in {it:fe()} avoids it. {p_end}

{synopt:{it:e(numpre_{it:type})}} The number of pre-treatment relative periods in the corresponding event-study vector, provided for use as {it:honestdid}'s {it:numpre()} argument. Not posted under {it:impvalues}.{p_end}

{synopt:{it:e(gteffects)}} The cohort-by-relative-period table of the {it:gteffects} option, whose entry gives the columns. {it:e(gt_ncontained)} counts the
cells whose standard errors are understated: cells lying wholly inside a single cluster, or holding a single observation under {it:robust}, pweights, a {it:customvce()}
of robust, hc2 or hc3, or {it:poisson} on {it:ppmlhdfe}. {it:e(gt_ndegen)} counts cells whose variance is not identified. Each is posted only when such cells exist.{p_end}

{synopt:{it:e(escomp)}} Event-study composition diagnostics, one row per reported horizon, aligned with the event-study rows of {it:e(wdidmainresults)}: {it:ncohorts},
{it:nunits}, {it:nclusters}, {it:nobs}, {it:sumw_att}, {it:wshare_att}, {it:wshare_att_it} and {it:relyr}. Composition shifts across horizons are a common cause of jumpy event-study
profiles. Per-subgroup versions are {it:e(escomp_s{it:sgcase})}. The scalars {it:e(escomp_refN)}, {it:e(escomp_refnunits)}, {it:e(escomp_refncoh)} and {it:e(escomp_refW)} describe
the reference-normalization pool, and {it:e(escomp_trtN)} and {it:e(escomp_trtW)} the whole treated estimation sample. {it:summarystats} displays the table.{p_end}

{synopt:{it:e(roth_power_{it:type})}} Under {it:pretrendpower()}: one row per requested power level, with columns {it:power}, {it:mdt_wald}, {it:mdt_nis},
{it:bias_wald} and {it:bias_nis}. {it:e(roth_slope_{it:type})}, {it:e(roth_slope_se_{it:type})}, {it:e(roth_slope_ll_{it:type})} and {it:e(roth_slope_ul_{it:type})} hold the
GLS pre-trend slope and its interval, {it:e(roth_prep_{it:type})} the joint pre-trend p-value, and {it:e(roth_hbar_{it:type})} the weighted mean post-treatment relative
period. The other scalars are {it:e(roth_npre_{it:type})}, {it:e(roth_keff_{it:type})}, {it:e(roth_W_{it:type})} and {it:e(roth_nissize_{it:type})}. {it:e(roth_ok)} = 1 marks
a completed pass. The hidden {it:e(roth_cellmax_{it:type})} records the largest treated cell size when the pre-period covariance was rank-deficient.{p_end}

{synopt:{it:e(esband_{it:type})}} Under {it:simultaneousbands}: the band endpoints per event-study family (rows {it:lb}, {it:ub}, {it:relyr}; {it:_s#} per subgroup), with
{it:e(esband_crit_{it:type})}, {it:e(esband_n_{it:type})}, {it:e(esband_engine)} ({it:analytic-simulated}, {it:bootstrap} or {it:bootstrap-studentized}), {it:e(esband_level)},
{it:e(esband_families)} and {it:e(esband_ndraw)}.{p_end}

{pstd} Other results in {it:e()}: {p_end}
{synopt: Joint Tests by Horizon} The joint tests of the pre-treatment event-study effects are in the scalars e(joint_horiz) for the ATT, e(joint_horiz_i) for the
shared ATT_it, ATT_i and ATT_itime horizon block, e(joint_Choriz) for the AME and e(joint_Choriz_i) for the AME_it, the last two under continuous treatments only. On an
analytic run the p-value is that of the {it:margins} joint test. Under {it:bootstrap()} and {it:permute()} it is a rank p-value of the replication statistics, stored in
{it:e(boot_jointW)} and {it:e(boot_jointWobs)} or {it:e(perm_jointW)} and {it:e(perm_jointWobs)}. {p_end}

{synopt: Joint Tests by Subgroup & Time Period} With subgroups, the same tests for each subgroup are in e(joint_horiz_s{it:sgcase}), e(joint_horiz_i_s{it:sgcase}),
e(joint_Choriz_s{it:sgcase}) and e(joint_Choriz_i_s{it:sgcase}), where {it:sgcase} is the subgroup value. {p_end}

{synopt: Adversarial Line Tests:} The p-values of the tests that all relative-period effects and the reference period fall on one line are in the scalars
{it:e(adversarial_pval_att)}, {it:e(adversarial_pval_ame)}, {it:e(adversarial_pval_att_s[SG])} and so on. They are posted only under {it:esfixedbaseperiod} and only when
the event-study covariance is not degenerate. The fitted line, with columns {it:slope} and {it:intercept_rt0}, is in the matrices {it:e(adversarial_line_att)},
{it:e(adversarial_line_att_it)}, {it:e(adversarial_line_ame)}, {it:e(adversarial_line_ame_it)} and their {it:_s{it:sgcase}} twins. A test that fails returns neither the
p-value nor the line. Without {it:bootstrap()} or {it:permute()}, the scalars {it:e(adversarial_df_att)}, {it:e(adversarial_df_att_it)}, {it:e(adversarial_df_ame)},
{it:e(adversarial_df_ame_it)} and their {it:_s{it:sgcase}} twins hold the numerator degrees of freedom each test used. {p_end}

{synopt: SG Effect Homogeneity Tests} For each estimate type, the joint test that all subgroups share the post-treatment effect of the lowest numbered SG is in
the scalar {it:e(joint_att_allsgssame)} and its analogues for the other types. The pre-treatment tests have names like {it:e(joint_pre_att_allsgssame)}. A test whose
pairwise covariance is degenerate is stored as missing.{p_end}

{synopt:{it:e(contreat_ndropcells)}} Under {it:contreat()}, the number of treated cohort-period cells dropped from the estimation sample because the dose does not vary within
them (see {help wooldid##continuous:contreat()}). Posted only when such cells exist.{p_end}

{synopt:{it:e(N)}, {it:e(N_dropped)}, {it:e(Nclus)}} The estimation sample size, the number of observations dropped from the specified sample, and the cluster count.
{it:e(N_clust)}, {it:e(clustvar)}, {it:e(vce)} and {it:e(vcetype)} are the conventional entries that {it:estout} and {it:esttab} expect. Under {it:poisson} or pweights
without {it:cluster()}, {it:robust}, {it:customvce()} or {it:unconditionalse} they read {it:robust} and {it:Robust}, the variance the regression estimates there.
{it:e(vce_stage1)} records the first-stage type under {it:unconditionalse}. {it:e(level)} holds the confidence level in force. {it:e(df_r)} holds the degrees of freedom of
the reported t statistics. Under {it:cluster()}, or a {it:customvce()} of cluster, that is one less than the number of clusters, and under multi-way {it:cluster()} one less
than the smallest cluster count. Otherwise it is the residual degrees of freedom. Under {it:poisson} it is missing, since the statistics are z. Under {it:bootstrap()} it
follows {it:bootdf()}. Under {it:impvalues} it is the smallest of the per-estimate values, which {it:e(imdf_min)} and {it:e(imdf_max)} bracket. {p_end}

{synopt:R^2 measures} For standard estimation, these include e(r2), e(r2adj), e(r2_within), e(r2_withinadj). For poisson estimation, these include e(r2pseud), a pseudo-R^2 measure. {p_end}

{synopt:{it:e(histogramestimates)}} The cohort i-specific treatment effects used to generate histograms, when requested.{p_end}

{synopt:{it:e(summarystats)}} Summary statistics on the outcome and any continuous treatment, under {it:summarystats}. The columns are {it:mean}, {it:median}, {it:sd},
{it:iqr}, {it:p25}, {it:p75}, {it:min} and {it:max}. The rows are {it:y_insample}, {it:y_pretreat}, {it:y_posttreat} and {it:y_control}, and under {it:contreat()}
{it:c_insample}, {it:c_pretreat}, {it:c_posttreat} and {it:c_control}. Some are also returned as scalars: {it:e(ypre_mean)}, {it:e(ypre_med)}, {it:e(ypre_sd)} and
{it:e(ypre_iqr)}, and under {it:contreat()} {it:e(contreat_mean)}, {it:e(contreat_med)} and {it:e(contreat_sd)} with {it:contreatpre_} and {it:contreatpost_} versions of
the same three. Under {it:contreatcontrols(untreatedZero)} the {it:c_pretreat} and {it:c_control} rows and the {it:contreatpre_} scalars are not returned. Under
{it:contreatwithin} neither control row is. {p_end}

{synopt: CFX Plot Info} The data behind each continuous treatment effects plot is a matrix named e(cfx[item]_[type]). [item] is att, att_i, att_it or att_itime for
average treatment effects, or ame, ame_i, ame_it or ame_itime for marginal effects of the continuous treatment. [type] is cfxpost for the post-treatment period, cfxpre for
the pre-treatment period, or relyrZ for relative period Z, with the letter n for a minus sign. Under {it:makecfxplotsbysg} each matrix carries the further suffix
{it:_sg[sgcase]}. {it:e(wdid_cfxmats)} lists every cfx matrix a run posted. When the plot standard errors are estimated, a bin whose dose-slope terms were all dropped for
collinearity has missing entries. If {it:saveplots()} cannot write {it:name}_cfx.dta on a run without {it:bootstrap()} or {it:permute()}, the run completes with a caution
and the hidden {it:e(wdid_cfxsave_rc)} holds the return code. {p_end}

{synopt:{it:e(sample)}} Marks the estimation sample, which shows the observations {it:wooldid} dropped during sample formation.{p_end}

{synopt:{it:e(cmd)}, {it:e(depvar)}, {it:e(cmdline)}} The command name, the outcome variable, and the input you fed to {it:wooldid} on the command line.{p_end}

{synopt:{it:e(estimate_type)}} A macro recording whether your treatment effect estimates use a pooled reference period or a fixed reference period. {p_end}

{synopt:{it:e(anticipation)}} The anticipation shift applied to the treatment onset, when {it:anticipation()} is active.{p_end}

{synopt:{it:e(wdid_vcvstate)}} What {it:e(V)} is on this run, in one word: {it:analytic-diagonal}, {it:analytic-blocks}, {it:bootstrap} or {it:absent}.
{it:analytic-blocks} is the default on an analytic run with an event study. Its within-family blocks are named in {it:e(wdid_vblocks)} and counted in {it:e(wdid_vblockn)}.
{it:e(wdid_fcol)}, {it:e(wdid_fcskip)}, {it:e(wdid_fcguard)} and {it:e(wdid_fcncell)} record whether the default cell collapse applied, see {it:nofastcollapse}. All are
hidden results, read by name. {p_end}

{synopt:{it:Other Items:}} {it:e(ll)}, {it:e(ll_0)}, {it:e(rmse)}, {it:e(rss)}, {it:e(F)}, and {it:e(chi2)} (poisson only) are passed through directly
from the underlying regression estimator, wherever available. {p_end}

{synopt:{it:Bootstrap inference}} Under {it:bootstrap()}, {it:e()} also holds the following. Run records:
{it:e(boot_reps)}, {it:e(boot_completed)}, {it:e(boot_failed)}, {it:e(boot_jk_completed)}, {it:e(boot_level)}, {it:e(boot_N_V)}, {it:e(boot_method)}, {it:e(boot_display)},
{it:e(boot_stat)}, {it:e(boot_tstat)}, {it:e(boot_seed)}, {it:e(boot_rngstate)}, and {it:e(boot_df)} when the reference distribution is a t. {it:e(boot_save_rc)} is posted
when a {it:bootsave()} write failed. {it:e(boot_nstraddle)} and {it:e(boot_reid)} record whether the replications re-identified duplicated units, and {it:e(boot_covsub)},
present only when that happened, counts the covariate terms rewritten to the fresh unit code. Covariance: {it:e(V_boot)}, also posted as {it:e(V)}, with the analytic one
kept as {it:e(V_analytic)}, and the replication means {it:e(b_boot)} and bias {it:e(bias_boot)}. {it:e(boot_v_posted)}, {it:e(boot_V_shortfall)} and {it:e(boot_ves_failed)}
record whether the covariance and the event-study covariances were rebuilt from the replications. Draws, up to 10,000 replications: {it:e(bootdist)}, and under {it:boottstat} {it:e(bootdist_se)} and,
for the subgroup contrasts, {it:e(bootdist_pwse)}. Families aligned to {it:e(b)}: {it:e(boot_ci_percentile)}, {it:e(boot_ci_normal)}, {it:e(boot_ci_bc)}, {it:e(boot_ci_bca)},
{it:e(boot_p_percentile)}, {it:e(boot_p_bc)} and {it:e(boot_p_bca)}, with {it:e(boot_z0)} and {it:e(boot_accel)}. Per-family tables: {it:e(boot_pctl_mainresults)}, also
under the alias {it:e(boot_percentile_mainresults)}, {it:e(boot_bc_mainresults)}, {it:e(boot_bca_mainresults)}, {it:e(boot_wald_mainresults)} and
{it:e(boot_tpctl_mainresults)}, with their {it:_pwsgs} and {it:_sgresults_{it:sgcase}} twins, {it:e(boot_{it:family}_gteffects)} for the cell table and
{it:e(boot_{it:family}_{it:cfxname})} for the cfx matrices. Test statistics: {it:e(boot_jointW)}, {it:e(boot_jointWobs)}, {it:e(boot_advW)}, {it:e(boot_advWobs)},
{it:e(boot_sgW)} and {it:e(boot_sgWobs)}, and the studentized records {it:e(boot_tstud)} and {it:e(boot_tstud_pwsgs)}. Families formed and tests withheld, with the reason:
{it:e(boot_jt_family)}, {it:e(boot_gteffects_family)}, {it:e(boot_adv_absent)}, {it:e(boot_bc_absent)}, {it:e(boot_jt_absent)}, {it:e(boot_sgjt_absent)} and
{it:e(boot_sgjt_reason)}. Per estimate: {it:e(boot_navail)} counts the replications available, and {it:e(boot_novar)} and {it:e(boot_sesmall)} those with no, or a very
small, standard error under {it:boottstat}. Hidden: {it:e(boot_ptails)}, {it:e(boot_comparator)}, {it:e(boot_ase_obs)}, {it:e(V_es_analytic_{it:type})} and
{it:e(df_r_analytic)}. {it:e(p)}, {it:e(ci_l)} and {it:e(ci_u)} hold the displayed family, so {it:esttab} reports the printed p-values and intervals. The standard errors
and statistics it reports come from {it:e(V)}. {p_end}

{synopt:{it:Permutation inference}} Under {it:permute()} the reported tables carry the permutation inference, in the columns {it:se_analytic}, {it:t_perm},
{it:p_perm} and the interval. The observed analytic tables are kept hidden as {it:e(perm_analytic_mainresults)}, {it:e(perm_analytic_sgresults_{it:sgcase})} and
{it:e(perm_analytic_pwsgs)}. {it:e()} also holds the following. One named column per coefficient of {it:e(b)}: {it:e(perm_p)}, {it:e(perm_ci)} with rows {it:ll} and
{it:ul}, {it:e(perm_se_shift)}, {it:e(perm_p_lo)}, {it:e(perm_p_hi)}, {it:e(perm_p_se)}, {it:e(perm_p_ci)}, {it:e(perm_sd)}, {it:e(perm_n)}, {it:e(perm_nullband)},
{it:e(perm_comparator_b)}, {it:e(perm_ase_obs)}, under {it:permtstat} {it:e(perm_comparator_se)}, and under {it:permstepdown} {it:e(perm_stepdown)} with its family, pre-treatment and subgroup variants. The ame
step-down tables are {it:e(perm_Cstepdown)} and {it:e(perm_Cstepdown_i)}. Draws, up to 10,000 replications: {it:e(permdist_b)}, and under {it:permtstat}
{it:e(permdist_se)}. Test statistics: {it:e(perm_jointW)}, {it:e(perm_sgW)} and {it:e(perm_advW)}, with their observed twins {it:e(perm_jointWobs)}, {it:e(perm_sgWobs)}
and {it:e(perm_advWobs)}. Run records: {it:e(perm_reps)}, {it:e(perm_completed)}, {it:e(perm_failed)}, {it:e(perm_npossible)}, {it:e(perm_npossible_str)},
{it:e(perm_level)}, {it:e(perm_exact)}, {it:e(perm_method)}, {it:e(perm_stat)}, {it:e(perm_pconv)}, {it:e(perm_seed)}, {it:e(perm_rngstate)}, {it:e(perm_analyticse)},
{it:e(perm_absp)}, {it:e(perm_adv)}, {it:e(perm_precaveat)}, {it:e(perm_irebuilt)}, {it:e(perm_blank)}, {it:e(perm_gt_blank)}, {it:e(perm_cfx_blank)},
{it:e(perm_offcentre)}, {it:e(perm_incomplete)} and {it:e(perm_addone)}, and {it:e(perm_save_rc)} when a {it:permsave()} write failed. Pooled-estimand scalars:
{it:e(perm_p_{it:est})}, {it:e(perm_bandlo_{it:est})}, {it:e(perm_bandhi_{it:est})}, {it:e(perm_sd_{it:est})}, and under {it:permtci} {it:e(perm_tlb_{it:est})} and
{it:e(perm_tub_{it:est})}, with {it:e(perm_comparator)} and {it:e(permdist)}. Hidden: {it:e(perm_kernel)},
{it:e(perm_pwsgs_kernel)}, and behind each cfx matrix {it:e(perm_analytic_{it:cfxname})}, {it:e(permdist_{it:cfxname})}, {it:e(permcmp_{it:cfxname})}, and under
{it:permtstat} {it:e(permdistse_{it:cfxname})} and {it:e(permcmpse_{it:cfxname})}. {it:e(p)}, {it:e(ci_l)} and {it:e(ci_u)} are posted, so {it:esttab} reports the
permutation p-values and intervals. Where {it:e(V)} is posted, the standard errors and t statistics it reports are the analytic ones. {p_end}

{pstd} Under {it:verbose} the underlying regression is stored as {it:___wooldidfullmodel} for {it:estimates restore}. {it:wooldid} overwrites this reserved
name on every run, so do not store your own estimates under it. {p_end}

{pstd} {ul:Tips for accessing output and remaking graphs:} {it:matlist} prints one of {it:wooldid}'s output matrices. {it:svmat} turns one into a dataset, for
example to remake the event study graphs. {it:svmat e(histogramestimates), names(col)} loads the cohort-specific effects, with the matrix's column names as variable names.
{help event_plot} can draw event study plots from {it:e()} when pointed at the coefficient names, for example {it:stub_lag(ES_att_post#)} and {it:stub_lead(ES_att_pre#)}
for the ATT. It does not insert the 0 for the omitted reference period. {p_end}

{pstd} {ul:Compatible Post-Estimation Commands:} {it:estimates save}, {it:estimates store}, {it:estimates restore}, {it:esttab}, and similar commands. {p_end}

{pstd} {ul:Incompatible Post-Estimation Commands:} {it:lincom}, {it:test}, and other commands that need a valid covariance matrix from the estimator.
{it:e(b)} and {it:e(V)} are for display. On an analytic run the off-diagonal elements of {it:e(V)} are 0 except within a single event-study family. {it:esttab} and the
like present correct estimates from them, but combining or testing different estimates is not, in general, valid. Under {it:bootstrap()} the posted {it:e(V)} is the replication covariance whenever one could be
formed, and {it:e(wdid_vcvstate)} then reads {it:bootstrap}.
The event-study coefficients also come with a full covariance matrix in {it:e(b_es_{it:type})} and {it:e(V_es_{it:type})}. {it:lincom} and {it:test} are valid on them
once they are copied and posted: {it:matrix b = e(b_es_{it:type})}, {it:matrix V = e(V_es_{it:type})}, {it:ereturn post b V, dof(`=e(df_r)')}. The {it:dof()} option
carries {it:e(df_r)}, so on an analytic run they use the degrees of freedom of the printed event study. Without it they use the normal and the chi-squared, which are
narrower than the table when clusters are few. Where {it:e(df_r)} is missing, {it:dof()} changes nothing. {p_end}

{marker comparison}{...}
{title:Comparison to Related Estimators}

{pstd} {ul:Comparing Output to hdidregress twfe:}
The Stata 18 command {help hdidregress} has an option for Wooldridge (2025) style estimates. {it:hdidregress twfe} with the option
{it:notyet}, which lets not-yet-treated cohorts serve as controls, should be comparable to {it:wooldid} with {it:att}, a pooled reference period without
{it:esfixedbaseperiod}, unconditional standard errors, and sample weights, if any, that do not vary within units. The two still differ under some kinds of panel imbalance, for example when a
treatment onset group contains several cohorts and some are not observed in periods when others are. {it:wooldid} then averages over the observed cohort-periods only. If cohorts 1 and 2 share
an onset and cohort 2 is observed for 2 post-treatment periods to cohort 1's 3, the two agree on the first 2 periods. For the third, {it:wooldid} reports the effect
both would report without cohort 2, and {it:hdidregress} does not. To match estimates with controls, adjust the {it:wooldid}
controls by hand. For example, a control that is constant within cohorts must be entered as {it:i.t#c.control} to replicate including it in {it:hdidregress}. {p_end}

{pstd} {ul:Comparing Output to did_imputation:}
The Borusyak, Jaravel, and Spiess (2023) imputation estimator is implemented in {help did_imputation}. For binary treatments, {it:wooldid} gives the same estimates
with {it:att}, a pooled reference period without {it:esfixedbaseperiod}, no pre-treatment horizons, that is {it:esprelength(0)}, no unconditional standard errors and no
controls. The pooled reference period is formed from the horizons the run asks for, so requesting pre-treatment horizons removes them from it and shifts every
post-treatment estimate. The level of {it:i} matters as well. {it:did_imputation} uses unit fixed effects, so
the estimates match only when {it:i} identifies individual units or every unit within an {it:i} is observed in the same periods. Otherwise they differ, and adding the
unit to {it:fe()} does not restore the match. With controls the estimates still
match when each cohort-period cell holds a single observation. With coarse cohorts, several observations per cell, they do not match outside special cases. {it:coarsecohortcontrols()}
does not restore the match there, but it arguably comes closer to the spirit of {it:did_imputation}. It makes the controls' effect more independent between the
treated cells on one side and the control and reference observations on the other. {p_end}

{pstd} {ul:Comparing Output to csdid and xthdidregress:}
The {it:gteffects} option returns the cohort-by-relative-period effects underlying {it:wooldid}'s aggregates. When {it:i} indexes treatment cohorts these are group-time effects, like the tables of
{help csdid} and {help xthdidregress}, but by default they are not the same numbers. Under the default pooled reference period each cell uses the not-yet-treated observations as controls
and the cohort's pre-treatment periods outside the requested event-study window as its baseline. With {it:esfixedbaseperiod} at the default {it:esrelativeto(-1)}, a balanced panel and no controls, the cells equal the group-time table
of {help csdid} with its {it:long2} option. {p_end}

{marker examples}{...}
{title:Examples}

{pstd}Examples 1 to 10 omit {it:cluster()} for brevity. On panel data, add it at the level at which treatment is assigned.{p_end}

{phang}0) Conventional diff-in-diff, targeting the ATT, with {it:i} the treatment onset cohort, unit fixed effects, and standard errors clustered by unit:{p_end}

{p 10 10} {cmd:wooldid y i t ttre, cluster(unit) fe(unit) att}

{phang}1) Conventional diff-in-diff, targeting both ATT and ATT_it, with production of histogram showing cohort-i specific effects: {p_end}

{p 10 10}   {cmd:wooldid y i t ttre, att att_it histogramcohorteffects}

{phang}2) Conventional diff-in-diff, with added controls and fixed effects:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, controls(controlvar1) fe(i.fe_bonus)}

{phang}3) Conventional diff-in-diff, using coarse cohort controls and including cohort-specific time trends:{p_end}

{p 10 10}   {cmd:wooldid y i t ttre, coarsecohortcontrols(c.controlvar1 i.categoricalcontrolvar) timetrends }

{phang}4) Event study with 3 pre-treatment and 3 post-treatment estimates using a pooled reference period and generating event study plots, as well as testing the hypothesis
that all pre-treatment effects are 0 and requesting the adversarial test (which requires {it:esfixedbaseperiod}, so here it fails with a caution): {p_end}

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

{p 10 10}   With cohort-specific time trends, let {it:timetrends} build the trend terms rather than placing them in {it:fe()} yourself. Estimates are better
conditioned when {it:esfixedbaseperiod} is also specified.

{p 10 10}   Continuous treatment controls can be numerically unstable. A higher {it:regtol()} or {it:ccnoabsorb} gives more consistent results. A cohort with
            almost no variation in the dose across its untreated or reference-period observations can look similar. Its baseline dose slope is then nearly unidentified,
            and the average marginal effect is unreliable in a way that {it:regtol()} does not repair. {it:wooldid} cautions when it sees this and refuses when the
            variation is smaller still. Examine or drop the affected cohorts.

{p 10 10}   The identification checks are not complete. With a binary treatment {it:wooldid} almost always detects a treatment effect partialled out by a
fixed effect, and it often does with a continuous treatment that has no lattice and no polynomial terms. With a lattice or polynomial terms, and especially in an
unbalanced panel, it is unlikely to detect important treatment effects partialled out by controls or fixed effects. Specify fixed effects and {it:contreatcontrols}
with care under a continuous treatment.

{p 10 10}   {it:wooldid} cannot always estimate standard errors for every requested object. If some reported estimates have identified standard errors and
others do not, {it:e(V)} is still returned, and the others carry a placeholder entry that displays as a missing standard error. Code that reads {it:e(V)} directly should
treat negative diagonal entries as not identified. If no reported estimate has an identified standard error, {it:wooldid} returns point estimates without {it:e(V)},
which can break some post-estimation tools. A caution prints in both cases.

{pstd}If you encounter a bug, whether a crash or a strange result, please reach out to me about it. Help finding and correcting bugs is much appreciated.
I am also happy to offer guidance on using this program. {p_end}

{pstd} If you encounter unusual errors that are hard to track down, try uninstalling and reinstalling {it:reghdfe} and {it:ftools} first. If {it:wooldid} keeps failing after a
{it:reghdfe} crash, type {cmd:clear mata} and {cmd:matrix drop _all} before the next call. These also remove your own matrices and Mata objects. {p_end}


{marker citations}{...}

{title:Citations}

{phang} Borusyak, Kirill, Xavier Jaravel, and Jann Spiess. Revisiting Event Study Designs: Robust and Efficient Estimation. Working Paper, 2023. {p_end}

{phang} Callaway, Brantly, Andrew Goodman-Bacon, and Pedro H. C. Sant'Anna. Difference-in-Differences with a Continuous Treatment. NBER Working Paper 32117, 2024. {p_end}

{phang} Callaway, Brantly and Pedro H. C. Sant'Anna. Difference-in-Differences with Multiple Time Periods. Journal of Econometrics, 2021, 225(2): 200-230. {p_end}

{phang} Sun, Liyang and Sarah Abraham. Estimating Dynamic Treatment Effects in Event Studies with Heterogeneous Treatment Effects. Journal of Econometrics, 2021, 225(2): 175-199. {p_end}

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
