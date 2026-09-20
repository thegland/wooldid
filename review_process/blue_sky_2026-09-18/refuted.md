# Findings that did not survive review

**Regenerated 2026-09-19 from the corrected dataset.** These are the leads the adversarial rounds
killed or reduced to nothing. They are recorded so nobody re-files them, and because a refuted lead
is evidence about the program too: each one is a place wooldid behaves correctly under a plausible
suspicion.

12 findings, of 223 total.

### `SKEP-N3` — permstrata() strata are built on the same unscreened cluster set

- **Verdicts:** REFUTED x4 (the prober who filed it, 4 later-round reviewers)
- **Was filed as:** medium; final severity **none**.
- **Claim:** The stratified permutation design counts clusters that contribute no rows to e(sample), so the reported number of possible assignments under permstrata() is overstated in the same way and by a larger factor, because the inflation compounds within strata.

### `EC-8` — After any failed wooldid, `estimates restore` gives back a valid result set that wooldid then refuses to replay, with a false explanation

- **Verdicts:** REFUTED x2 (the prober who filed it, 2 later-round reviewers)
- **Was filed as:** medium; final severity **none**.
- **Claim:** The replay-staleness refusal is keyed on a session global WOOLDID_LASTRC that no `estimates restore` can clear, so after a failed call a user who restores a stored good wooldid result (e(cmd)=wooldid, valid b/V) is refused with 'the stored results predate that failed call' - which is untrue of the restored results - even though sthlp:1299 states est restore/est save are compatible with wooldid.

### `MCBP-JOINTSIZE` — The analytic joint horizon test rejects a true sharp null 75.5% of the time at nominal 5% — far outside the 11-36% oversizing the project documents

- **Verdicts:** REFUTED x2 (the prober who filed it, 2 later-round reviewers)
- **Was filed as:** medium; final severity **none**.
- **Claim:** On a G=10 event-study spec with `jointtests`, the ANALYTIC joint test over event-study horizons rejects a TRUE SHARP NULL at .755 [.619,.854] at nominal 5%. The project's own record characterises the analytic joint/adversarial tests as oversized at 11-36%; .755 is 2-7x that band and its Wilson lower bound .619 clears .36 decisively. The adversarial test at .286 [.178,.424] does sit inside the documented band.

### `SK-01` — semielasticity averages (dp/dx)/p over a subpopulation containing negative predictions and reports the result with se, t, p and a 95% interval, unscreened

- **Verdicts:** PARTIAL x2 REFUTED x2 (the prober who filed it, 4 later-round reviewers)
- **Was filed as:** medium; final severity **none**.
- **Claim:** wooldid's semielasticity asks margins for eydx, which divides by the prediction, and wooldid never checks the sign of the divisor. On a panel where 14 of 720 fitted values are non-positive -- three of them inside the treated-and-post subpopulation the att is averaged over -- the run completes, prints 'Treatment Effects Presented as Semi-elasticities', and reports a point estimate with a standard error, t = 3.05, p = 0.003 and a 95% interval, with no caution of any kind. A '% change in y per unit change in x' is not defined at those rows, so the headline figure is contaminated to an extent the…

### `SKEP-03` — Where a subgroup's continuous event study is 'skipped entirely', e(b_es_ame_s#) is still posted as a 1x1 matrix of MISSING under the junk colname c1, with e(numpre_ame_s#) missing

- **Verdicts:** REFUTED x2 (the prober who filed it, 2 later-round reviewers)
- **Was filed as:** low; final severity **none**.
- **Claim:** In an ordinary ANALYTIC contreat() + subgroups() + ES run, a subgroup for which wooldid announces 'no continuous-treatment relative period effect is estimable ... so its continuous-treatment event study is skipped entirely' nonetheless gets e(b_es_ame_s#) posted as a 1x1 matrix holding a missing value under the column name 'c1', and e(numpre_ame_s#) is missing rather than 0 or absent.

### `fp-contreat-D` — fe() with customvce()/unconditionalse: wooldid cautions, then enters the fe() variables as CONTINUOUS regressors and completes with rc 0

- **Verdicts:** REFUTED x2 (the prober who filed it, 1 independent refuter writing its own repro, 1 third-lens reviewer)
- **Was filed as:** low; final severity **none**.
- **Claim:** On the plain-regress paths wooldid prints a caution that absorbed fixed effects are not supported and that the user must specify them as controls with i. prefixes -- and then proceeds to paste the bare fe() varlist into the regression with no i. prefix, fitting a linear trend in the identifier instead of a fixed effect. The run completes rc 0 and reports fit statistics and an average marginal effect from that model.

### `skept-roworder-01` — The same rows in a different physical order change an event-study standard error by 3.7e-08 relative, on the reghdfe path with continuous weights only

- **Verdicts:** REFUTED x2 (the prober who filed it, 2 later-round reviewers)
- **Was filed as:** low; final severity **none**.
- **Claim:** With a weight variable whose values are not exactly representable in binary (any ordinary survey/administrative weight), wooldid's reghdfe path returns a different event-study standard error depending on the caller's physical row order. On a 30-state x 10-year panel with `att esprelength(3) espostlength(3) [aw=wt] cluster(st)`, sorting the identical rows differently moves one horizon's se/t/p/lb/ub by 3.732e-08 relative - 370x the project's stated 1e-10 indifference threshold. The point estimates do not move. The root cause is reghdfe's weighted clustered VCE, which I isolated directly: bare r…

### `FASTNOSE-02` — fast mode prints "SEs and related materials are not estimated in fast mode" directly above a full gteffects table of standard errors, p-values and 95% confidence intervals

- **Verdicts:** REFUTED (the prober who filed it, 1 independent refuter writing its own repro)
- **Was filed as:** low; final severity **none**.
- **Claim:** Under fastnose with gteffects, the header says standard errors were not computed and every Main Results inference cell is missing, while the gteffects table in the same output carries complete se/p/lb/ub columns that are byte-identical to the analytic run's. sthlp:286 says fast mode computes 'no standard errors, test statistics, p-values, or confidence intervals ... for any estimand', and the gteffects paragraph never mentions fastnose.

### `SKEP-N5` — NEGATIVE RESULT: an event-study window does not drop rows, so one of the three triggers ado:9367 names does not exist

- **Verdicts:** PARTIAL REFUTED (the prober who filed it, 2 later-round reviewers)
- **Was filed as:** low; final severity **none**.
- **Claim:** ado:9367's comment names three things that make a replication drop rows the bootstrap frame still holds - 'missing y, an event-study window, singletons'. Missing y and singletons reproduce; the event-study window does not change e(sample) at all.

### `hv-lat-04` — The granular-i caution tells the user runtime, not the variable limit, is what will bite under lattice() -- on both designs measured the variable limit bit first, at unremarkable runtimes

- **Verdicts:** REFUTED (the prober who filed it, 1 independent refuter writing its own repro)
- **Was filed as:** low; final severity **none**.
- **Claim:** ado:1371 and sthlp:1163 assert that under lattice() 'runtime grows ... about its 1.9th power ... so a granular i can dominate runtime long before it approaches Stata's variable limits'. At the cell counts where the variable ceiling is first breached under lattice(), runtime is still tens to low hundreds of seconds, so a user reading that sentence has no reason to suspect maxvar. The caution also does not print at all at i = 50 distinct values, which is where the failure I measured happened.

### `state-leak-03-mata-namespace-leftovers` — A plain run leaves three globals in the user's Mata workspace, one of them named `k` - attribution to wooldid NOT established

- **Verdicts:** REFUTED (the prober who filed it, 1 independent refuter writing its own repro)
- **Was filed as:** low; final severity **none**.
- **Claim:** One plain `wooldid y unit year ttre, cluster(state) att` on a session whose Mata workspace was empty leaves it holding three real scalars - __00000L, __00000M and `k` - and leaves two scalars in the Stata scalar namespace (__gtools_gfile_ghdfeabs = 65, __gtools_gfile_gregvcov = 65). A user with their own Mata variable named `k`, an utterly ordinary name, would have it silently overwritten by running wooldid. I could not establish whether wooldid or one of its hard dependencies creates them.

### `SKEP-N4` — NEGATIVE RESULT: impvalues is clean on this axis

- **Verdicts:** CONFIRMED x2 (the prober who filed it, 2 later-round reviewers)
- **Was filed as:** low; final severity **none**.
- **Claim:** The Ibragimov-Muller path counts only clusters with rows in the estimation sample, so wholly dead clusters do not reach its cluster count or its degrees of freedom.
