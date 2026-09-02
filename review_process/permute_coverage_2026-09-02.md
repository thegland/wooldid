# Every p-value under permute(): the coverage scan (2026-09-02)

Assessment requested by the author: "I recently found out that permute only focused on
permutation testing the high-level statistics, like the ATT, but not the event study
horizons or anything else. We want to make sure that when permute is enabled, *every*
p-value that wooldid might produce is produced via a permutation test. That includes
inference on event study estimates, adversarial test p-values, joint test p-values, etc.
-- everything! Please do a scan to document all parameters that can in principle be
tested via permutation testing but which are not right now. If any would be difficult or
problematic to approach with permutation testing, please let me know and make the case
for this - I am open to some exclusions (e.g., ruling that permute cannot be combined
with adversarial or joint, or with subgroups, etc. - make the case and let me know).
Please surface these possible exclusion cases to me." Assessment only; nothing in the ado
changes here. Companion document for `bootstrap()` is drafted separately.

**Headline.** Today `permute()` produces one permutation p-value per requested estimand
(`att att_it att_i att_itime`), stored as `e(perm_p_<est>)` with its lb/ub/sd companions,
read from `_b[main_<est>]` at `wooldid.ado:5943`/`5983` and printed only in the separate
randomization-inference block at `7252-7271`; the main table's own p column stays
analytic and unlabelled (chk_a.log 448-449 prints p 1.66e-17 while `e(perm_p_att)` =
.0323 at 493). **Eighteen reported p-values are not permutation-based under `permute()`
today**, and every one of them prints - or, for the cfx lattice matrices, is posted to `e()` -
without a provenance marker (with two reachability riders: the four joint pre-trend items reach a
permute run only on a mixed-cluster-granularity design or under `safety(off)`, per ruling 1, and
the three adversarial items only under `esfixedbaseperiod`): the pretreatment and event-study horizon rows, all four joint
pre-trend tests plus their per-subgroup twins, the adversarial line tests (pooled and
per-subgroup), the three subgroup tables and their `ame` twins, the pairwise subgroup contrasts,
the equality-across-subgroups tests, every `ame` row under `contreat()`, the `gteffects` cell
p's, the cfx lattice p's and `pretrendpower()`'s joint pre-trend p. `wooldid_permctl` refuses
five option combinations outright (section 1). The recommended scope is:
ship the missing analytic-leftovers caution first (permute has no analogue of the
bootstrap's `wooldid.ado:256`), then cover the free and unambiguous families - `ame`
under `contreat()`, event-study post horizons, a family-level `gteffects` max statistic,
and the V-free adversarial statistic - on top of one shared per-element harvest
refactor, and exclude the cfx lattice, plot bands, per-cell `gteffects` p's and the
V-only `pretrendpower` quantities. **Six rulings govern the rest** (section 5 and section 6 list the same six): (0) label every
un-permuted p, or blank it - the literal reading of the request; (1) subgroups in replications,
one ruling governing four families; (2) pre-period statistics, one cross-cutting ruling governing
seven items; (3) the `jointtests` exit - refuse the combination or carve it in with the analytic
tests suppressed; (4) the display and provenance contract; (5) the `e()` naming scheme, which
must be fixed before any family ships its own. Section 4 carries seventeen numbered rulings
below those, including the `ame` blanked-draw threshold and the `gteffects` coverage form.

**A caveat on evidence provenance.** The judges who produced sections 2-4 held the Stata
slots; this drafter had none. Line citations marked *(verified this session)* were read
directly out of `wooldid.ado` while writing; everything else is carried forward from the
judges and the advocate unaltered, including their `~` approximations. Where two judges
conflicted, both positions are stated and the advocate's resolution is marked as such.

A separate verification pass afterwards re-read 145 of the load-bearing citations, re-checked
every row of the section-2 matrix against the code (all 62 statuses held; both counting claims
are arithmetically correct), and re-ran the defect claims in Stata. It corrected about a dozen
pointers in place and **refuted one structural claim**: permute + `jointtests` is not unreachable
on every dataset, only on designs whose `cluster()` granularity is uniform - ruling 1 now carries
the counter-example. Logs: `pb/VER_defect/`.

---

## 1. What `permute()` covers today

**The engine.** The permute shell occupies `wooldid.ado:101-177` (verified this session:
the block opens `if "\`wdp_on'" == "1" {` at 101 and closes after
`wooldid_display , part(perm)` at 175). It runs *after* the observed full estimation is
complete and stored; the observed `e()` is restored at ~141 and the permutation summaries
are appended on top of it. Nothing in the observed surface is overwritten - measured:
`mreldif(e(wdidmainresults), analytic)` == 0 (chk_d.log 315).

**What the rep string keeps and strips** (verified this session, `wooldid.ado:105-133`).
Kept: the estimand flags, `esfixedbaseperiod`, `timetrends`, `poisson`/`poisexpresults`,
`robust`, `ccnoabsorb`, `latticeignoreweights`, `emptycellsoverride`, `cluster()`, `customvce()`,
`fe()`, the control families,
`customweightmultiplier`, `fastcollapse`/`fullprecisionslow`, the full `contreat()`
family (`contreat`, `contreatpoly`, `contreatcontrols`, `semielasticity`,
`contreatelasticitytype`), `lattice()`, `esprelength`/`espostlength`/`esrelativeto`,
`safety`, `regtol`, `oostreatedcontrols`, `anticipation`. Reps run under `fastnose`
unless `permtstat` is in force (110). Stripped by explicit v1 scope choice, with the
comment "subgroups are stripped: v1 permutation statistics cover the pooled estimands
only": `subgroups()`, `jointtests`, `adversarial`, `gteffects`,
`histogramcohorteffects`, `pretrendpower`, and every plot option.

**The statistic.** `wooldid_permrun` (`~5828-6092`) allocates
`matrix res = J(wdp_R, : word count estimands, .)` at 5910 (verified) and reads exactly
one scalar per estimand per rep: `_b[main_<est>]`, or `_b[main_<est>]/_se[main_<est>]`
under `permtstat` (5983/5991, verified). The rep-0 comparator lives in four fixed
tempname scalars, `tempname wdpo1 wdpo2 wdpo3 wdpo4` (5833, verified) - a hard cap of
four elements that no judge noticed and that any per-element harvest must replace with a
comparator vector as its *first* change, not as an optimization. The two-sided p is
formed at 6058-6062 (verified): `(1 + nex)/(nuse + 1)` sampled, `nex/nuse` under
`permexact`.

**The `e()` surface added.** `e(perm_p_<est>)`, `e(perm_lb/ub/sd_<est>)`,
`e(perm_tlb/tub_<est>)` under `permtci`, `e(permdist)` (posted only when completed reps
<= 10000, shell 169), `e(perm_comparator)`, `e(perm_method)`, `e(perm_stat)`,
`e(perm_reps/completed/failed/npossible)`.

**The display.** `part(perm)` at `7252-7271` prints one line per estimand plus the
5th-95th null band, correctly labelled "a null band, not a parameter CI" (7255), and the
`permtci` pivotality hedge at 7266-7268. Everything else in the run - `part(1)`,
`part(2)`, the Roth block, the plots - prints inside the worker before the shell block
ever runs: `part(1)` and `part(2)` are gated at `4237`/`4241` on
`if "\`bootstrap'" == ""` with **no permute gate**, and `wooldid_ptp` prints between them
at `4240`. That ordering fact governs the whole display contract in §5: a display-side
provenance tag keyed on `e(perm_method)` would fire only on replay, never live.

**What is missing relative to bootstrap.** The bootstrap shell prints, at `wooldid.ado:256` -
option-gated on `makeplots`/`makecfxplots`/`saveplots`/`gteffects`/`histogramcohorteffects`/
`adversarial`/`altmethod_adversarial`, so it fires only when one of those is requested
(verified this session), "Caution: plots, saved plot data, gteffects cell standard errors,
histogram outputs and adversarial line-test p-values are produced from the ANALYTIC
covariance during estimation; the reported tables and e() inference above are
bootstrap-based." The permute block has no analogue of that caution. The help file is not
silent about permute's scope, though: `sthlp` 378-380 already documents two carve-outs - "Under
`subgroups()` the permutation statistics cover the pooled estimands only (rep runs skip the
subgroup machinery)" and "Permutation statistics cover the att-family estimands (the
continuous-treatment ame quantities are not permuted in this version)". What is neither
documented nor printed is that the analytic p's for those surfaces keep printing unmarked inside
the same run.

**Refusals** (`wooldid_permctl`, 5232-5391). Five option combinations are refused outright:
a `permute()` argument other than `shuffle`/`probabilities` (5238), `impvalues` (5242),
`unconditionalse` (5246), `contreatwithin` (5250) and `contreatcontrols(untreatedZero)` (5260,
whose comment records the adjudication that `contreat()` in general stays allowed); on the design side permute also requires a single `cluster()` variable (5312-5315) whose
clusters nest inside one onset group (5343-5347), and refuses a design in which every stratum is
onset-homogeneous so the number of distinguishable reassignments is 1 (5897-5900). Those design
refusals are as much "exclusion cases" as the option ones, and ruling 1 turns on them.

---

## 2. The coverage matrix

Every inference item wooldid can report, with where it is computed, what reference
distribution it uses today, and whether `permute()` touches it. Status values: `covered`
= a permutation p exists; `not_covered` = an analytic quantity prints under `permute()`;
`na` = no p-value exists to permute; `refused` = the option is refused with `permute()`.

| # | item | options required | computed at (`wooldid.ado`) | analytic reference | permute status | evidence / note |
|---|---|---|---|---|---|---|
| 1 | `main_att` - pooled `main:<est>` rows: se, t, p, CI | none (`att` default) | 1985-2020, `r(table)` at 1991 | `margins r(table)`, `ttail(e(df_r))` | **covered** (p only, in a separate block) | `e(perm_p_<est>)` appended 145-168 from `_b[main_<est>]` (5943/5983), summarised 6026-6068, printed 7252-7271. The main-table p column itself stays ANALYTIC and unlabelled: chk_a.log 448-449 p 1.66e-17 vs `e(perm_p_att)` .0323 at 493 |
| 2 | `perm_p_summary` - `e(perm_p/lb/ub/sd_*)`, `e(perm_tlb/tub_*)`, `e(permdist)`, `e(perm_comparator)`, counts | `permute()` + `permreps()`; `permtstat`/`permtci` | 6026-6068 (`_pctile` 6042/6055); shell 145-174 | add-one two-sided perm p; 5/95 pctile band; SD | **covered** | covers only `wdp_est` = att att_it att_i att_itime (134-135); `perm_tlb/tub` mix an analytic observed se with a permuted \|t\| quantile (159-165), hedged at 7268 |
| 3 | `pre_att` - `pretreatment:<est>` rows | none (`esprelength`>0 default) | 2037-2070, `r(table)` 2043 | `margins r(table)`, analytic VCE | **not_covered** | column present in every rep's `e(b)` but discarded; `e(perm_p_pretreatment_att)` absent (chk_a.log 499); printed p analytic .0479 (chk_a.log 452), no label |
| 4 | `es_horizon_att` - `ES_<est>_pre#/post#` rows | none | 2100-2120, `r(table)` 2106, V capture 2116 | `margins over(_t) r(table)` | **not_covered** | colnames present in every rep (chk_a.log 428-438); analytic p printed unlabelled (chk_a.log 456-459) |
| 5 | `es_plots_ci` - event-study plot bands and saved plot data | `makeplots` / `saveplots` | CI columns from 2100-2120 consumed by the worker plot block | analytic `r(table)` lb/ub | **not_covered** (not a p-value) | no permutation band; plots use analytic CIs; permute prints no analogue of the bootstrap caution at 256 |
| 6 | `ves_att` - `e(b_es_*)`/`e(V_es_*)` blocks incl. `_s#` | none; `subgroups()` adds `_s#` | 3990-4030 (captures 2116/2570/2782/3234) | analytic `e(V)` subblock | **na** | permute produces no covariance; blocks stay analytic and unlabelled |
| 7 | `ves_withheld` - `e(V_es_*)` withheld/downgraded | `contreat()` with coarse cluster; degenerate horizons | 4012-4030; bootwire 6484-6485 | analytic screen decides existence | **na** | degenerate-row caution 4026 |
| 8 | `joint_pre_att` - `e(joint_horiz)`, `e(joint_horiz_i)` | `jointtests` (refused with coarse cluster: SI-A1 1232-1263) | 2205-2220, `r(p)` 2218 | `test` after `margins over(_t)` | **not_covered** | `jointtests` stripped from reps (109); analytic p prints unlabelled under `safety(off)` (chk_b.log 380-383; `e(joint_horiz)` 1.2085e-24 at 464). SI-A1's own error text (1257) recommends `permute()` for this test |
| 9 | `joint_pre_ame` - `e(joint_Choriz)`, `e(joint_Choriz_i)` | `jointtests contreat()` | 2650-2666, `r(p)` 2664 | `test` after `margins dydx over(_t)` | **not_covered** | as above; `ame` columns never read by permrun |
| 10 | `joint_pre_sg` - `e(joint_horiz_s#)` etc. (locals 4137-4140) | `jointtests subgroups()` | 2965-2998 (`r(p)` 2996) att; 3275-3296 (3294) ame | `test` after `margins subpop() over(_t)` | **not_covered** | no refusal exists in permctl (5231-5391); subgroups stripped (107-109) and the analytic p prints (part(2) 7182-7196) |
| 11 | `adv_line_att` - `e(adversarial_pval_att/att_it)`, `e(adversarial_line_*)` | `adversarial` (or `altmethod_adversarial`) + `esfixedbaseperiod` | 2140-2160, `r(p)` 2153 via `wooldid_advx` 8388-8389; posted 4110-4135 | Mata GLS line fit + Wald on analytic V | **not_covered** | `adversarial` stripped (109); analytic p prints unlabelled (chk_b.log 141-144 "single line: 0.0000"). On `cluster(cohort)` the composite-null covariance is degenerate and no p is posted at all (chk_b.log 162-165) |
| 12 | `adv_line_ame` | `adversarial esfixedbaseperiod contreat()` | 2595-2602, `r(p)` 2600 | `wooldid_advx` on analytic V | **not_covered** | as `adv_line_att` |
| 13 | `adv_line_sg` - `e(adversarial_pval_<est>_s#)` | `adversarial esfixedbaseperiod subgroups()` | 2951 (att), 3260 (ame); posted 4142-4155 | `wooldid_advx` on analytic V | **not_covered** | printed part(2) 7200-7206; inherits the SI-7 dropped-horizon refusal (2938-2946) |
| 14 | `sg_main_att` - subgroup tables `main_SG#_*`, `pretreatment_SG#_*`, `ES_*_SG#_*` | `subgroups()` | 2700-2795 (`r(table)` 2707/2731/2772/2793) | `margins subpop() r(table)` | **not_covered** | three sg tables print analytic p's with no label (chk_a.log ~200-300); no `e(perm_p_*SG*)` |
| 15 | `sg_ame` - subgroup ame rows | `subgroups() contreat()` | 3165-3220 (`r(table)` 3170/3189/3216) | `margins dydx subpop()` | **not_covered** | analytic |
| 16 | `pwsgs` - pairwise subgroup contrasts, `e(pwsgs)` cols 2-6 | `subgroups()` with >1 level | 3039/3044/3114/3119 (att); 3324/3329/3368/3373 (ame) | `margins pwcompare r(table_vs)` | **not_covered** | printed analytic under permute; degenerate on the probe data - se 2.33e-10, p 6.53e-87, one row missing (chk_a.log 510-554) |
| 17 | `allsgssame_att` - `e(joint_<est>_allsgssame)`, `e(joint_pre_<est>_allsgssame)` | `subgroups()` (>1 level) | 3060-3080 (`test` 3065, contrast fallback 3070-3076), 3137-3155; posted 4160-4163 | `test` on analytic VCE chain | **not_covered** | analytic; all '.' on the probe data (chk_a.log 306-309), 1.5e-81 elsewhere (chk_d.log 1031-1034) |
| 18 | `allsgssame_ame` | `subgroups() contreat()` | 3343-3354 / 3389-3400; posted 4164-4169 | `test` on analytic VCE | **not_covered** | analytic |
| 19 | `allsgssame_ame_i` - `ame_i`/`ame_itime` equality tests | `subgroups() contreat()` + `att_i`/`att_itime` | 4164-4169 (`subinstr` att->ame) | `test` on analytic VCE | **not_covered** | analytic. Covering it by name-driven list closes a bootstrap-side gap - see §4 ruling 17 |
| 20 | `ame_main` - `main:ame`, `pretreatment:ame`, `ES_ame_*` | `contreat()` | 2245-2310 (`r(table)` 2249/2302), 2552 (ES; V 2570) | `margins dydx r(table)` | **not_covered** | `contreat` kept in reps (121-122) but `wdp_est` is the att family only (134-135); `e(perm_p_ame)` absent (chk_d.log 305), `wdidmainresults` identical to analytic (315), RI block lists att only (233) |
| 21 | `cfx_lattice` - lattice cell estimate/se/p/CI, saved plot data, `e(cfx*)` | `contreat() makecfxplots` | 2405-2520 (`r(table)` 2411/2416/2440/2445/2490/2516); posted 4191-4194 | `margins r(table)` | **not_covered** | `makecfxplots` stripped from reps; analytic |
| 22 | `gteffects_cell_se` - `e(gteffects)` cell se, t/z, p, CI | `gteffects` | `wooldid_gtx` 7818-7823; posted 4034-4037; refused with contreat() at 4422-4424 | `2*ttail(df,\|z\|)`, `invnormal` crit on analytic V | **not_covered** | `gteffects` stripped from reps; analytic. The cells are the summands of `main_att`, which does get a permutation p |
| 23 | `hist_cohort` - `e(histogramestimates)`, histogram graphs | `histogramcohorteffects` | 1857-1920 (`margins ... nose`); posted 4175 | none (`nose`: no VCE) | **na** | no inference exists on this surface at all |
| 24 | `roth_prep` - joint pre-trend F/chi2 p, `e(roth_prep_<est>)` | `pretrendpower() esfixedbaseperiod` | `wooldid_ptp` 6724-6725 | Mata GLS on analytic `e(b_es_*)`/`e(V_es_*)` | **not_covered** | `pretrendpower` stripped from reps; analytic p prints unlabelled under permute (chk_b.log 390/411 p 0.0000; `e(roth_prep_att)` 1.2e-24 at 470). permctl has no ptp gate at all (4609 omits it; 4623 passes it to bootctl) |
| 25 | `roth_slope` - `e(roth_slope_se/ll/ul_*)`, `e(roth_nissize_*)` | `pretrendpower() esfixedbaseperiod` | 6610-6611 and surrounding GLS | analytic `e(V_es_*)` | **not_covered** (not a p-value) | analytic; se .0013738 printed under permute (chk_b.log 389). The se and nissize are functions of S and m alone - no permutation analogue exists |
| 26 | `roth_power` - `e(roth_power_<est>)`, `e(roth_W/hbar/keff/npre)` | `pretrendpower() esfixedbaseperiod` | 6593-6870 | noncentral quadrature on analytic `e(V_es_*)` | **na** | power diagnostic; no null contrasted |
| 27 | `im_rows` - Ibragimov-Mueller cross-cluster t rows | `impvalues` | 2015-2017, 2064-2066, 2177-2179, 2274-2276, 2327-2329, 2623-2625 | `ttest r(se)/r(p)` over cluster-level estimates | **refused** | permctl 5242 |
| 28 | `boot_ci_flavors` - `e(boot_ci_*)`, `e(bootdist)`, `e(V_boot)` etc. | `bootstrap() reps()` | 5808-5819; `wooldid_bootci` 8713-8741 | replication matrix | **na** | mutually exclusive with `permute()` |
| 29 | `e_V` - posted `e(V)` and `e(df_r)` | none | worker post ~3960-3980; bootwire 6326-6500 | analytic VCE | **na** | RI produces no covariance, so `e(V)` necessarily stays analytic |
| 30 | `sescan_caution` - identical-SE screen, placeholder-entry caution, AMESE containment caution | none (data-triggered) | 7424-7540; 3960-3980 | analytic `e(V)` diagonal | **na** | appropriate under permute, except that the 7443 text recommends "randomization inference" *inside* a run that already has permute on (chk_d.log 156-160) |
| 31 | `summarystats_desc` - `e(summarystats)`, `e(escomp*)`, N/Nclus/R2 | `summarystats` (escomp always) | worker | descriptive counts | **na** | nothing to permute |

**Totals.** 31 items. Two carry a permutation p today (#1, #2), and #1's permutation p
lives in a separate block rather than in the table. **Eighteen** items are reported
p-values that are analytic under `permute()` (#3, 4, 8-22, 24). One is refused (#27).
The remaining ten carry no p-value to permute.

---

## 3. Family-by-family feasibility

### 3.1 Event study and pretreatment rows (`pre_att`, `es_horizon_att`, `es_plots_ci`)

**Statistic.** Per-element two-sided randomization p, one per named `e(b)` column in the
att-family ES block, formed exactly as the existing pooled p (6038-6041, 6061-6062) with
the statistic taken from the rep's named column instead of `_b[main_<est>]`. Under
`permtstat`, `_b/_se` per horizon. Output per horizon: p, 5th-95th null band, permuted
SD, per-element completed count.

**Null and validity.** The null is the *same for every element*: reassignment of onset
labels across clusters (5965-5974) makes any statistic computed from the reassigned data
an exact test of the single sharp null that timing affects no unit in any period.
Choosing `ES_att_post3` over `main_att` changes the direction of power, not the null - so
a per-horizon p cannot be read as "the effect at horizon 3 is zero given the other
horizons". Post horizons are clean. **Pre horizons are contaminated**: under a real effect
a reassignment drags truly-post periods into the rep's pre-window, inflating the reference
distribution, so the pre-row p is *conservative exactly when the effect is large* - the
wrong direction for an alarm. `es_plots_ci` is not a p-value at all: RI yields a null band
centred at zero, never a parameter CI, and a real RI confidence band needs test inversion
at R x G rep runs per horizon.

**Engineering.** Zero rep-string change is needed: `esprelength`/`espostlength`/
`esrelativeto`/`esfixedbaseperiod` are already emitted (109, 126-129), which is why every
rep's `e(b)` already carries the ES columns. Replace the fixed R x #estimands matrix
(5910) with an R x K accumulator over the rep-0 column set and share `wooldid_bootgrab`
(8437-8464) verbatim - it is already name-aligned via an asarray of the rep colstripe
(8448-8453) and already leaves a cell missing when the rep lacks that column
(8451-8459). Comparator: `matrix compname = e(b)` after rep-0, bootrun's pattern at 6185,
which preserves the NFPERM-12 bitwise-reconstruction property element by element. The
real work is **availability accounting**: today a rep failing on any estimand is blanked
whole-row (6002-6010, verified) so `nuse == ndone` for every column, and permdist packing
keys on column 1 (6072). At K ~ 10-40 that rule is catastrophic - one missing horizon
would discard the rep for `main_att` too.

**Cost.** Zero additional replications and essentially zero per-rep time; the ES columns
are already computed and discarded (5983). Marginal cost is one mata call per rep, a
K-times-longer summary loop, and an R x K permdist.

**Recommendation: cover_with_caveat.** Post horizons: **cover** - the cleanest per-element
RI target in the package. Pre rows: cover only with a printed sentence stating that the
permutation p for a pre row is a sharp-null p, not a parallel-trends test, and that a true
effect widens its null distribution; suppression is preferable to shipping it unlabelled.
`es_plots_ci`: **exclude**, matching the bootstrap's own precedent (256), but store the
per-horizon null band in `e()` so a user can plot it by hand. Optional and nearly free: a
Westfall-Young max-\|studentized\| step-down over the *post* block, exactly valid because
the joint null *is* the sharp null.

**Effort.** ~300-450 changed lines across 3-4 gated commits plus a battery re-pin.
Prerequisite commit: fix shell:126 (see §4 ruling 15).

### 3.2 Joint pre-trend tests (`joint_pre_att`, `joint_pre_ame`, `joint_pre_sg`)

**Statistic.** A per-rep *scalar built from coefficients only*, never from a per-rep VCE:
the fixed-metric Wald `Q_r = b_r,pre' Sinv b_r,pre`, with `Sinv` held fixed across every
rep and the pre-block read by name (`ES_<fam>_pre#`, `ES_<fam>_SG<lvl>_pre#`). One-sided
tail count (Q is a nonnegative quadratic form), so the two-sided `abs()` of 6038 does not
apply. Preferred metric: the covariance of the completed reps' own pre-block - it needs no
analytic covariance at all, so it survives the SI-A1 configuration where 1257 declares the
cluster-robust covariance unfit, and it survives `amescrub`, where `e(V_es_ame*)` is
withheld (4016-4021) and the analytic ame joint test is skipped outright (2656). A
permutation test's *exactness* does not depend on `Sinv` being a correct sampling
covariance - only its power does.

**Null and validity.** This is the sharpest null mismatch in the scan. The label promises
"conditional on real post effects, are the pre coefficients zero?"; RI delivers "is there
any effect anywhere?". Contamination direction: the permutation joint pre-trend p **rises
with the true effect size**, so the pre-trend test becomes easier to pass exactly when
there is an effect to identify. On the one reachable probe the scales make this dominant,
not marginal: the analytic joint p is 1.2085e-24 (chk_b.log 464) while `e(perm_p_att)` is
.0323 (chk_b.log 457) on 12 clusters. A middle path worth recording so it is not
re-proposed: conditioning the reference set on "assigned onset <= true onset for every
cluster" removes the contamination but the conditioning event is a function of the
realized assignment and the set is not invariant under permutation within itself - the
conditional-randomization argument fails.

**Engineering.** Reachability is the first blocker (see §4 ruling 1). `jointtests` must
*not* be re-added to the rep string: `fastnose` + `jointtests` is refused at parse time
(4969-4972), so re-adding it forces non-fastnose reps and a full VCE each. The pooled
items are therefore free. Do not extend permrun's whole-row blanking (6002-6008) to the
new columns; keep `badrep` keyed on the estimand columns so the shipped `e(perm_p_att)`
stays bit-identical, and handle per-element completion listwise at the statistic (the
bootjw rule at 8622/8638/8642). **Trap:** computing Q on whatever pre-block columns a rep
happens to have makes Q a *different function per rep*, which is not a permutation
statistic at all - listwise-complete reps only. Degenerate columns need the calibrated
scale-free guard from `wooldid_bootjw` (8655-8672).

**Cost.** Zero incremental per-rep estimation for the pooled items; one k x k covariance
and R quadratic forms post-loop. `joint_pre_sg` is the only real cost and it belongs to
the subgroups ruling.

**Recommendation: cover_with_caveat**, with the mandatory caveat printed in the display
and not only in the help, and with the reachability precondition. `joint_pre_sg` falls
out cleanly if subgroups-in-reps is declined. Preferred `e()` design: suppress the
analytic joint tests under permute (add `permute` to the four gates that already exclude
`bootstrap` at 2208, 2656, 2964, 3275) and reuse the existing names, so the display at
7043-7057 and 7180-7196 needs no change - which also dissolves the SI-A1 circularity
outright. Deferred, not rejected: reps on the not-yet-treated sample, the only variant
that tests the null the label promises.

**Effort.** Pooled items ~200-240 net lines across 4 commits; `joint_pre_sg` +60-80.

### 3.3 Adversarial line tests (`adv_line_att`, `adv_line_ame`, `adv_line_sg`)

**Statistic.** Three candidates; only (C) is recommended. (A) fixed-metric Wald with V
frozen at the observed covariance - unavailable exactly where RI is most wanted, because
`advx` aborts on a rank-deficient V (8341-8342) and this is **measured** to fire on the
flagship permute config (chk_b.log B3, `cluster(cohort)`: "the composite-null covariance
is degenerate; no p-value computed", both types). (B) the rep's own analytic p as the
statistic - mechanically available (the only `WOOLDID_REPMODE` gate in the file is the F8
composition block at 3717-3721), but 4965-4966 hard-refuses `fastnose` + `adversarial`,
so reps must pay a full VCE each. **(C) recommended: the V-free residual SSR around the
published display line.** `wooldid_advstat` already fits that line by unweighted OLS in a
temp frame (7354-7365) over exactly the rows the user sees, with the altmethod cut at
7361; `local advrss = e(rss)` beside the existing `matrix advB = e(b)` at 7363 yields the
statistic. Always computable, no rep-string change, no VCE.

**Null and validity.** The analytic null is composite - there *exists* a line fitting the
horizon coefficients (8287-8296) - and the sharp null is a **strict subset** of it. So an
RI p here rejects when "not every effect is zero", which in any study that requested an
adversarial pre-trend test is not in doubt. Two structural consequences, both analytical
and unprobed: the distortion direction is indeterminate (misalignment inflates W;
smoothing a step over misaligned onsets deflates it toward a line), and the statistic has
**zero power against an exactly-linear-in-event-time effect** - precisely the alternative
the analytic test is built to treat as a null. W is also non-monotone in effect size by
construction, so the two-sided convention at 6038 is meaningless for it.

**Scope refinement (advocate, not in any judge's write-up).** The analytic adversarial p
exists *only* under `esfixedbaseperiod`: `wooldid_advstat` gates its entire inference
block on `if "\`esf'" != ""` (~7371-7392), returning a fitted line and status `nop` with
no p on pooled-reference runs. So the "unlabelled analytic p under permute" defect is
confined to `esf` runs, and roughly half of adversarial usage has nothing to cover or
mislabel.

**Recommendation: cover_restricted_form.** Mandatory regardless of the ruling: the
shell-side caution, which must be shell-side rather than a display tag because part(1) and
part(2) print at 4237/4241 before `e(perm_method)` is appended at 145. Cover
`adv_line_att`/`adv_line_ame` with statistic (C), printed in `part(perm)` under its own
sharp-null heading, never reusing the wording at 7064/7203. Exclude `adv_line_sg` pending
the subgroups ruling - it inherits the SI-7 dropped-horizon refusal (2938-2946) at a rate
that makes the missingness effect-dependent. Do not build impose-and-permute in this pass;
it is the only route to the advertised null but imposes a unit-level *homogeneous* linear
effect strictly stronger than the estimand-level null.

**Effort.** P0 caution ~3-5 lines; P2 (this family) ~60-90 lines / 1 commit on top of the
shared prerequisite. Impose-and-permute is ~150-250 lines and research-grade.

### 3.4 Subgroups (`sg_main_att`, `sg_ame`, `pwsgs`, `allsgssame_*`)

**Statistic.** Levels: per-element `|b|` tail count on the `main_SG#_*` /
`pretreatment_SG#_*` / `ES_*_SG#_*` columns, already present by name in every rep's
`e(b)` (48 columns in the author's own 2-estimand / 3-subgroup / 4-horizon probe).
Pairwise: per-rep difference of the two named sg columns, tail-counted against the rep-0
comparator difference - the exact analogue of `wooldid_bootpw` (8681-8695) with a tail
count replacing the SE. Equality-across-subgroups: a Mahalanobis
`W_r = d_r' inv(S_perm) d_r` on the chain contrast (`wooldid_bootjw` chain mode,
8624-8632), with `S_perm` the *permutation* covariance, degrading to
`max_j |d_j| / sd_perm(d_j)` when rank-deficient - which is exactly the case where the
analytic path posts '.' or 1.5e-81 today.

**Null and validity.** Subgroup membership is genuinely fixed under reassignment - wooldid
already refuses a unit split across subgroups (866-877) and sg varying within an i-by-t
cell (846-859) - so each rep is a well-defined re-randomization. Levels: exact test of the
global sharp null with a subgroup-focused statistic; the ordinary Fisher-vs-Neyman gap,
honest with one line of label. **Contrasts are the material substitution**: under real but
*common* effects the observed contrast is ~0 while each rep's subgroups are contaminated
differently (they differ in cohort composition), so the permutation distribution of the
contrast is inflated and the test is **conservative and underpowered for heterogeneity**.
It will not manufacture heterogeneity findings; it will miss real ones.

**Engineering.** Add `subgroups()` to the rep string mirroring the bootstrap at 214.
Three new rep-failure modes appear, all consequences of reassignment changing which units
are treated, and all silent today: the negative-sg assert over `ifstatement & ttre < .`
(818-822, r411) can fire on a unit that is never-treated observed but treated in a rep -
whole-rep loss *correlated with the assignment*; `sglist` is `glevelsof` over
`ifstatement & t >= ttre` (823), so a subgroup can vanish from a rep (tolerated by
bootgrab) or all can vanish (hard r459 at 825-827); the integer screen (834-838) and F13
screens have `ttre`-dependent domains. Two hard prerequisites: port the bootstrap's
cell-collapse fallback rule (201-208, verified) into the permute shell, since `fcdeg`
fires specifically on subgroup margins calls (2761, 3031, 3086) and the permute shell
passes `fastcollapse` unconditionally at 119 (verified); and add a permctl-time screen for
negative/non-integer sg over the *full* in-sample population. Also mirror bootctl's
`suppresssgprimaryeffects` refusal - under it the sg columns never enter `e(b)`
(3902-3909), so nothing is harvestable.

**Cost.** Per-rep margins count goes from ~3 to ~(5 + G) per estimand - roughly
(1 + 0.6G)x, i.e. ~2.5x at G=3, doubled again under `contreat`. The decisive existence
proof: the bootstrap already pays exactly this multiplier and ships at `reps(1000)` on the
same machinery.

**Recommendation: cover_with_caveat**, in three tiers. Tier 1 **cover** the level rows
outright. Tier 2 **cover** `pwsgs`/`allsgssame` *only* with the printed null-stating line;
without it, exclude. Tier 3 sg pre rows follow the cross-cutting pre-period ruling
verbatim. Do not gate the equality tests on the analytic scalar existing the way bootjnt
does at 5708 - the degenerate-analytic case is precisely where RI adds value.

**Effort.** Shared prerequisite ~200 ado + 40 mata over 2 commits; sg-family marginal
~470 ado + 95 mata over 4 commits; ~650-800 lines total from HEAD.

### 3.5 Continuous treatment / `ame` (`ame_main`)

**Statistic.** Six-to-many scalar elements already present by name in every replication's
`e(b)` and currently discarded: `main_ame`, `main_ame_it/_i/_itime`,
`pretreatment_ame(+flavors)`, `ES_ame_pre#/post#`, `ES_ame_it_*`. Confirmed live: the D1
probe's `e(b)` is 1 x 12 with `e(b_es_ame)` 1 x 4 (chk_d.log 291-297). Same tail count as
the pooled att p, on a different named column.

**Null and validity.** Already adjudicated by the author, in code and in the help:
permctl 5257-5258 records "contreat in general stays ALLOWED (adjudicated: the null is
timing-as-if-random given dose)", and the help states that null verbatim to users. The
remaining help sentence - "the continuous-treatment ame quantities are not permuted in
this version" - is a scope note, not a statistical position. The pre rows carry the same
contamination as everywhere else, in the conservative direction.

**The one real asymmetry versus att.** The `ame` is a within-cell dose *slope*, and
whether it is identified depends on the assignment: a reassignment leaving a cell with no
within-cell dose variation drops the slope terms for collinearity, `amevoid` fires
(1607-1613) and the entries re-enter `e(b)` as exact **structural zeros** (3912-3919,
because `ereturn post` cannot carry missing coefficients). Harvested naively, `|0| <
|observed|` always, so every voided rep counts as a non-exceedance and the p is
**anti-conservative**. The remedy exists and is value-based: `wooldid_bootgrab` blanks
exactly this pattern at 8455-8458. Reusing bootgrab is a correctness requirement, not an
optimization. Blanking is however a *missing-not-at-random* deletion - the voiding
reassignments are the ones concentrating low-dose-variance clusters into the treated
window - so exactness becomes conditional and must be disclosed.

**The `permtstat` hazard.** `permtstat` drops `fastnose` (110), `amescrub`'s gate fires
when `fastnose` is empty (1333-1344), and the canonical permute design (cluster inside one
onset group: a single cluster variable, 5312-5315, that must nest inside one onset group, 5343-5347) makes cluster == i, so every treated i-t cell is contained and the
ame se is scrubbed (the scrub column floor at 1495; the scrub itself at 2256-2260, 2309-2313, 2558-2564). Measured: "72 treated i-t
cell(s) lie completely within a single cluster" (chk_d.log 155-158). Adding `ame` to the
t-set naively makes rep-0 abort at 5945-5952 and blanks every rep at 5987-6004 - turning
working runs into `r(459)` in exactly the configuration RI is recommended for.

**Cost.** Marginal per-rep cost is **zero**: `contreat` and its modifiers are already in
the rep string (121-125, verified) and the ame margins passes are not gated on `fastnose`
(1961 only downgrades to `nose`; 2250-2254 pads `r(table)`).

**Recommendation: cover.** The strongest case in the scan, and the only family where
permute produces p-values the package *cannot* produce analytically: when containment
kills the analytic ame SEs, wooldid prints "Consider switching to use of bootstrap SEs, IM
p-values, or randomization inference" (7443, echoed 7546) - and the measured permute run
prints that caution, then shows all six ame rows with se/t/p/lb/ub entirely blank
(chk_d.log 180-226) while `e(perm_p_ame)` does not exist (305). Three non-negotiable
conditions: (a) the bootgrab ame structural-zero blank; (b) per-element statistic
selection so ame falls back to `b` under `permtstat`; (c) disclose blanked-draw counts
with a withholding threshold (author's call, §4 ruling 4).

**Effort.** ~110-155 net ado lines plus tests; ~35-55 if a sibling family lands the shared
accumulator commit first.

### 3.6 Diagnostics (`gteffects_cell_se`, `cfx_lattice`)

**`gteffects` - statistic.** Family-level scalar: `gtmax_post` = max over treated cells
with `relyr >= 0` of `|estimate/se|`, with a **one-sided** add-one p (a max-of-absolutes
is already one-tailed, so 6038's two-sided form must not be reused), plus `gtmax_pre` over
`relyr < 0` carrying the pre-period caveat. Restrict to `permute(shuffle)`.

**Why it is uniquely cheap.** The block is `gen` + five `gegen` totals + a tag + `mkmat`
+ a mata regex pass (1653-1676) - no regression, no margins; the help calls it "at
essentially zero additional computation" (849). And its SEs come from
`vecdiag(e(V))` off the **stage-1 fit** (1668, with a `J(1,K,.)` fallback at 1669), so they
survive `fastnose` (help 860-862). This is the one place in the permute family where
t-scale randomization is free, in contrast to `permtstat`, which forces the whole loop off
`fastnose` (110).

**Why per-cell p's are refused.** A per-cell RI p tests the *global* sharp null with that
cell as the statistic, while a reader of a group-time cell reads "this cohort had no
effect in this period", and the two diverge in both directions. Cells exist only where
`targetflag == 1` (7797), so observed cell (i,t) appears in a rep only when i's reassigned
onset is <= t - every per-cell p is a *conditional* randomization test on a different
reference set with small heterogeneous counts. `e(gteffects)` is one row per treated
(i,relyr) cell, up to N x T (the display refuses to print above 200 rows, 7242), so
marginal per-cell p's are indefensible without a Romano-Wolf step-down. And
`e(gteffects)` is a **matrix**, not `e(b)` columns, so it cannot reuse `wooldid_bootgrab` -
it needs a key-based `(i, caltime)` merge.

**`cfx` - why it is excluded.** The estimand *drifts across reps* for exactly the
documented usage: `baseTreated*`, `byTreated*` and `byTtre*` lattices (950-970) re-derive
their bins from the permuted onset each rep, so `lattice_3` denotes a different dose range
every rep - and those are the types in the shipped help examples (1158, 1164). It is the
most expensive per-rep item in the family (2 margins passes per cfxplottype x (1 + #sg),
the code's own comment calling it "a pretty taxing request" at 2407). It cannot run in
reps at all until a no-graph internal mode exists: `marginsplot` is unconditional at 2417,
2456, 2496, 2532. And the artifact users consume is a plot whose band is a *parameter CI*,
which an RI null band cannot replace.

**Recommendation: cover_restricted_form for `gteffects`, exclude `cfx`.** Ship the caution
for both. Add `gtmax_post`/`gtmax_pre`, shuffle-only, refused or noted under `timetrends`
(where the se column is missing by design). Surface the **coverage-by-subtraction**
alternative to the author: under permute, *blank* the analytic se/p/lb/ub columns of
`e(gteffects)` rather than computing permutation ones (~10-20 lines, in-repo precedent in
`amescrub`'s cfx blanking at 2446-2453 and `fastnose`'s cfx truncation at 2360) - it
satisfies "every p under permute is a permutation p" trivially, at the cost of destroying
the one inference surface that survives fast mode.

**Effort.** Caution ~12-18 lines; `gtmax` ~90-150 lines over 2-3 commits; per-cell p's
~300-450 lines; cfx ~250-400 lines and highest risk.

### 3.7 Pre-trend power (`roth_prep`, `roth_slope`)

**The organizing cut.** The family splits on what depends on `b` versus on `V` alone, and
only the b-dependent half admits a permutation statistic. b-dependent:
`roth_W = b_pre' Sinv b_pre` (6722-6723) and its p (6724-6725), and the GLS slope
(6717-6718). **V-only, with no permutation analogue under any design**:
`roth_slope_se = 1/sqrt(m' Sinv m)` (6719) does not involve `b` at all;
`roth_nissize` (6743) is a null rejection rate of a pre-test, a function of S alone;
`mdt_wald` (6763) and `mdt_nis` (6790-6803) are minimum-detectable-trend power
quantities. A permutation distribution of a power calculation contrasts no null.

**Consequence for the two inventory rows.** `roth_prep` is the family's only p-value and
the only coverable object. `roth_slope` (is_pvalue false) has no coverable p: its se and
nissize are V-only, and an RI p for the slope *point* would be a **new** reported
statistic - a new reported number rather than a permutation version of an existing one, which is
the reason to leave it out of scope (the plan's R4 list at 279-280 names three other items, not
this one).

**Statistic if covered.** Fixed-metric Wald: hold `Sinv` at the observed value built once
by `wooldid_ptp`'s own route (pre-index 6655-6686, symmetrize 6702, `invsym` + `keff`
6704-6710), permute only `b_pre`, one-sided tail count. Needs no per-rep VCE, so reps stay
under `fastnose`. `W_obs` **must** be recomputed from the rep-0 fastnose pre columns with
the same `Sinv` - not from `e(roth_W_<type>)` (6860, full-run provenance) - or permexact's
observed-is-a-member property breaks; the discipline and its measured failure mode ("b >=
local-copy-of-b evaluates FALSE") are recorded at 5915-5921.

**Cost.** Zero marginal per-rep estimation for the fixed-V form. What is *not* proposed,
for scale: a per-rep refit of S would require dropping `fastnose` and re-adding
`pretrendpower` per rep with its NIS quadrature - the in-code budget note at 6731-6733
records ~2.9s per evaluation at k=6 with ~85 evaluations per power level per type, plus an
80-iteration bisection at 6790. Hours added at R=1000, to calibrate against the very
analytic covariance permute exists to avoid.

**Recommendation: exclude_pending_ruling**, in two tiers. Tier 0, unconditional: permute +
`pretrendpower` runs *today* with no gate at all - permctl (5232-5391) has no ptp check
and validate omits `ptp()` from the permctl call at 4609 while passing it to bootctl at
4623 - so an analytic F p prints inside a run headlined "Randomization Inference" (7254)
with no marker anywhere. Do **not** mirror bootctl's refusal: allow the combination and
add the Roth block to the Tier-0 caution. Unlike bootstrap, permute does not invalidate
the Roth block, and the V-only mdt/hbar diagnostics are exactly what a few-cluster user
wants. Tier 1: `roth_prep` **is** the same joint test as `e(joint_horiz)` up to covariance
source and g-inverse handling (compare 6659-6665/6722-6723 with bootjnt's 5653 and
`wooldid_bootjw` 8613-8677) and must be ruled with it, sharing the code path.

**Effort.** Tier 0 10-30 lines / 1 commit. Tier 1 150-250 lines, of which the great
majority is shared with the joint family; 90-130 lines / 1 commit if the joint family
lands first.

---

## 4. Exclusion cases for the author's ruling

Seventeen rulings. Each states the case for exclusion, the case against, and the
advocate's verdict. Rulings 1, 2 and 11 govern the others and should be decided first.

**1. `jointtests` + `permute()`.** *For exclusion:* every `jointtests` p in this family is
a pre-trend test (display text 7047) and the permutation null is the sharp null, so the p
*rises* with the true effect - the alarm goes quiet exactly when there is an effect to
identify. Every other misspecification in this package errs toward over-rejection with a
visible caution; this one errs toward silent reassurance on an identification check. And a
refusal costs nothing that works today. *Against:* the package promises this in its own
error text, twice - 1257 refuses `jointtests` under coarse clustering and 1261 repeats it
as a caution, both saying "Consider wooldid's built-in randomization inference (the
permute() option) ... for joint-type inference here", and `permute()` provides no joint
test. Zero incremental per-rep cost. And a permutation test is exact for any *fixed*
statistic, so the covariance defect that makes the analytic test unusable costs the
permutation version power, not validity. *Verdict:* **author must rule; the status quo is
indefensible either way.** New structural finding upgrading every judge's "measured on
this config only": on any design whose `cluster()` granularity is **uniform**, permute +
`jointtests` is unreachable under default safety, by three independent legs - (a) every i inside
one cluster triggers SI-A1 at 1232, r(198); (b) every i split across clusters means shuffled reps
assign different onsets within one i, so the f13 onset-constant-within-i screen (447-463, error
459) kills every rep, which is exactly the measured "only 0 of 30 permutation replications
estimated successfully"; (c) the apparent escape `permstrata(i)` makes every stratum
onset-homogeneous, so M == 1 and 5897-5900 refuses. Leg (b) holds only under default safety,
since f13 is gated on `safetyinput != "off"` (447). All three legs were reproduced directly
(`pb/VER_defect/d1.log` L1/L2/L3 with the f13 mechanism proven in `d2.log`).

**The universal version of that claim is false, and the counter-example matters.** A *mixed*
granularity design escapes all three: `pb/VER_defect/route4.dta` (18 units, 24 clusters; six
never-treated units split across two clusters each inside one onset-homogeneous stratum, twelve
treated units one cluster each in a heterogeneous stratum) runs
`wooldid y unit t ttre, att att_it cluster(clv) esprelength(2) espostlength(2) esfixedbaseperiod
jointtests adversarial permute(shuffle) permreps(30) permstrata(...)` to completion at **default**
safety - rc 0, 30 of 30 reps estimated, the joint tests computed, printed and posted
(`d2.log` R4A). SI-A1 never fires because its `jtnic == jtni` test at 1240-1250 finds 18 units
against 24 clusters. Two consequences for the ruling: exit (A) is a **breaking change** for such
a run, not a no-op; and on that run the analytic joint p prints directly above the randomization
-inference headline at default safety, which is the single strongest argument for the
gate-suppression half of exit (B). Two coherent exits: **(A)**
refuse in permctl (~6 lines, threading `jt()` into the 4609 call which today omits it)
*and* rewrite 1257/1261 to stop recommending a capability that does not exist - the
smallest honest move, and the advocate's lean; **(B)** carve `permute` out of 1232 as
`bootstrap` already is, *and* suppress the analytic joint tests under permute (add
`permute` to the gates at 2208/2656/2964/3275), *and* post fixed-metric-Wald permutation
p's under the existing `e(joint_*)` names. **Never ship (B) without the gate suppression**
- carving 1232 alone leaves the disowned analytic p printing under a permute headline.

**2. `subgroups()` in replications - one ruling, four families.** *For exclusion:* a
recorded, deliberate v1 scope choice in two places (`f6_permutation_plan.md` §2 item 7 and
R4 ~279-280) *and* documented to users in the help. Real cost on the slowest path,
~2.5x at G=3 before `contreat` doubles it, paid by every permute run that merely happens
to have `subgroups()` in the command - **and that multiplier is an operation count, not a
measurement**: no timing run was made for this assessment, so the ruling that turns on it should
be preceded by one (the repo's own house style is to quote a measured ratio). New assignment-correlated rep-kill risk nobody has
characterised on real data (818-822, 834-838, 846-877), so the surviving reference
distribution can be *biased* rather than merely thinned. And the contrast items substitute
the null. *Against:* the current state is the worst of both worlds - the sg tables,
`e(pwsgs)` and the allsgssame scalars print analytic p's with no label, and part(1)/part(2)
are gated only on `bootstrap` (4237, 4241) so they print before the shell block runs. The
analytic numbers are frequently degenerate in exactly the few-cluster regime permute exists
for (pwsgs se 2.33e-10 / p 6.53e-87; allsgssame '.'), while the bootstrap's rewiring of the
same rows produced usable numbers. The level items are unimpeachable and free. And the
bootstrap already carries `subgroups()` in its rep string (214) at `reps(1000)`. *Verdict:*
**author must rule, once, for all four families** - this is the highest-leverage decision
in the scan. Recommended: Tier 1 cover the level rows; Tier 2 cover `pwsgs`/`allsgssame`
only with a printed null-stating line, else exclude Tier 2; Tier 3 sg pre rows follow
ruling 11. Two non-negotiable prerequisites if covered: port the cell-collapse fallback
rule (201-208 -> after 119) and add the full-population sg screen in permctl. If declined,
the only deliverable is the Tier-0 caution and a help sentence naming the sg tables as
analytic.

**3. `adversarial` / `altmethod_adversarial` + `permute()`.** (Both judgments of this family
recommend the same primary statistic; one also offered a second, **under `permtstat` only**: use
the rep's own analytic adversarial p as the statistic, which is mechanically available since the
only `WOOLDID_REPMODE` gate in the file is the F8 composition block at 3717-3721. It costs the
fastnose saving, so it belongs with ruling 12's t-statistic decision.) *For exclusion:* the null
mismatch is structural and uniquely bad - the sharp null is a strict *subset* of the
composite adversarial null, so an RI p rejects when "not every effect is zero", which is
not in doubt in any study that requested this test. Zero power against an exactly-linear
effect, the alternative the analytic test is built to tolerate (analytical, unprobed). The
statistic is non-monotone in effect size by construction. *Against:* today's state is a
pure labelling defect the author's request names directly, and coverage via the V-free
route is nearly free and is *available exactly where the analytic test fails* - measured
degenerate on the flagship config. *Verdict:* **cover_restricted_form**, with a scope
refinement no judge stated: the analytic adversarial p exists only under
`esfixedbaseperiod` (`wooldid_advstat` gates its whole inference block on `esf`,
~7371-7392), so the defect and the deliverable are confined to `esf` runs. Ship the
shell-side caution regardless; cover the pooled items with the residual-SSR statistic under
its own sharp-null heading; exclude `adv_line_sg` pending ruling 2; do not build
impose-and-permute in this pass.

**4. `contreat()` / `ame` + `permute()`.** *For exclusion:* one argument that the first judgment
of this family raised and the second dropped - `contreatcontrols(byCohort)` reps can fail at
higher rates through the by-cohort screen (7449-7503), thinning the reference distribution in a
way that is correlated with the assignment. Plus the one real asymmetry versus
att - the ame is a dose slope whose *identification* depends on the assignment, so the
reference set is a set of reassignments on which the estimand is not uniformly defined,
and blanking voided draws is a missing-not-at-random deletion that makes exactness
conditional rather than exact. Plus scope discipline: the exclusion is recorded and stated
to users, and two working remedies already exist. *Against:* the package prescribes the
remedy and then withholds it - the containment caution says "Consider ... randomization
inference" (7443) and the measured permute run shows all six ame rows blank with
`e(perm_p_ame)` absent. It is the *modal* case, not a corner (72 contained cells measured).
The null is already adjudicated in code (5253-5258). And the marginal per-rep cost is zero.
*Verdict:* **cover** - the strongest inclusion case in the scan, with the three
non-negotiable conditions in §3.5. **Author must rule** on the blanked-draw threshold
(suggested: caution always, withhold the ame p above ~10% blanked).

**5. `pretrendpower()` + `permute()`.** *For exclusion:* most of the printed Roth block has
no permutation analogue at all - the slope se, nissize and mdt columns are functions of S
and m only. Covering one number out of five leaves a block needing a legend. And bootstrap
refuses the pair outright, so refusal is the house precedent. *Against:* `roth_prep` is the
analytic p *least* defensible in permute's home regime, and the block's own code concedes
the calibration stakes ("a chi-squared shortcut is optimistic by ~25 percent at 10
clusters"). Unlike bootstrap, permute does not invalidate the block. And the combination
runs silently today with no permctl gate. *Verdict:* **author must rule between two coherent
exits.** (A) Mirror `bootctl`'s refusal (5578-5583) in `wooldid_permctl`: the parity move, it
closes the unlabelled-analytic-p gap by removing the combination, and it is what one of the two
independent judgments of this family recommended as its primary remedy. (B) Cover the label and
exclude the statistic: allow the combination, add the Roth block to the Tier-0 caution, and leave
`roth_prep` analytic-but-marked - which keeps a diagnostic the user asked for and which permute
(unlike bootstrap) does not invalidate. The two judgments of this family split on exactly this,
so it is presented as a fork rather than a recommendation. Rule `roth_prep` jointly with ruling 1
- they are the same statistic. Note the ordering consequence: because `jointtests` is structurally
unreachable on uniform-granularity designs while `pretrendpower` + permute runs freely,
`roth_prep` is on those designs the
*only reachable* pre-trend joint p under permute, which makes its labelling more urgent
than its coverage. Exclude `roth_slope` outright.

**6. `gteffects` + `permute()`.** *For exclusion:* per-cell p's test the global sharp null
while the table invites a cell-specific reading; cells do not exist independently of the
assignment; multiplicity up to N x T; and the v1 plan listed gteffects as observed-run
only. *Against:* the cells are literally the summands of a statistic that already gets an
RI p; coverage is almost free and uniquely free of the fastnose tradeoff; and the
family-level max has no null mismatch at all, since "is any cell nonzero" *is* the sharp
null. *Verdict:* **author must rule between three shapes** - the two independent judgments of this
family disagreed on which, so all three are on the table. (i) **`gtmax` only**: one family-level
max-|t| statistic with a one-sided tail, shuffle-only, refused or noted under `timetrends`; no
per-cell numbers. This is the lean. (ii) **Per-cell RI p's**: one permutation p per treated cell
over completed reps, with a per-cell count column and an `n_k < 20` blank rule, posted as
`e(perm_gteffects)` plus an optional raw-draw matrix `e(permdist_gt)` - the form that makes the
table genuinely useful for user-built aggregates, at ~150-220 lines. Both judgments agree that
marginal per-cell p's are indefensible **without** a Romano-Wolf max-|t| step-down over the
`inwindow==1` rows, so (ii) means (ii)-plus-step-down or not at all. (iii)
**Coverage-by-subtraction**: blank the analytic se/p/CI columns under permute (~10-20 lines) -
ruling 0's answer applied here. In every shape the caution is mandatory and the note at 7241 must
say the columns are analytic. Note the
new cross-cutting dependency: `gteffects`' window is `gthstart = 1 + min_t - espre` (1670, with the window end at 1672 and the pooled-regime override at 1671),
so it inherits the `espre` drift of ruling 15.

**7. `makecfxplots` / cfx lattice + `permute()`.** *For exclusion:* the estimand drifts
across reps for exactly the documented lattice types; highest per-rep cost in the family;
cannot run in reps at all until a no-graph mode exists; and the deliverable is a plot whose
band is a parameter CI that an RI null band cannot replace. *Against:* under `amescrub`
the cfx ame inference columns are already blanked, so a permutation p would be the only
inference available there; and the fastnose harvest target is clean. *Verdict:* **exclude**
- the clearest exclusion in the scan. The prerequisites alone exceed the whole cost of the
ame family, which delivers strictly more. Strong optional add: refuse
`makecfxplots` + `permute()` under default safety, or at minimum name it in the caution.

**8. `histogramcohorteffects` + `permute()`.** *For exclusion:* definitional, not a
judgment call - the histogram harvests `r(b)` from margins calls that already carry `nose`
(1857-1885), so there are no SEs, no t-ratios and no p-values anywhere on that surface.
The same one-liner applies to `summarystats` and `e(escomp*)`. *Against:* none on the
statistics; only that the bootstrap caution at 256 names "histogram outputs", so omitting
it permute-side would make the two engines' cautions disagree. *Verdict:* **exclude,
definitionally, no author ruling needed.** The entire deliverable is naming it in the
permute-side caution for parity. Recorded explicitly only so it is not re-litigated: the
author listed it among the combinations to rule on, and the correct answer is "there is no
p-value here", not "it would be hard".

**9. `makeplots` / `saveplots` - plot confidence bands.** *For exclusion:* not a p-value
and not convertible. RI yields null bands, never parameter CIs; substituting one would be
read as a CI - strictly worse than the analytic graphic - and the display text already
draws that distinction (7255). Plots are produced inside the worker before the shell block,
so attaching a permutation band needs a display-order refactor, precisely why the bootstrap
left plots analytic. *Against:* saved plot data and the drawn band are the artifacts users
circulate, and under permute they are analytic and unlabelled. *Verdict:* **exclude the
band swap; ship the label.** Two cheap partials: join the Tier-0 caution, and if
per-horizon coverage lands, store the null band in `e()` so a user can plot it by hand.

**10. `permute(probabilities)` versus `permute(shuffle)` for any per-element family.** *For
exclusion:* shuffle preserves the onset multiset, so relative-year support is structurally
stable; probabilities does not, so which horizons, cells and subgroup levels are estimable
genuinely varies, and a per-element p over only the reps where the element was estimable is
conditional on a data-dependent event and no longer exact. For a *max-over-elements*
statistic it is worse than conditioning: the number of elements maximised over varies, so
the statistic is not the same function in every rep, which voids the permutation argument
outright. *Against:* the conditioning is disclosable for single-element statistics, and
refusing probabilities would fragment the option surface. *Verdict:* **author must rule,
with a sharp split by statistic type.** Recommended: single-element statistics allowed
under both methods with mandatory per-element counts and a probabilities-specific caution;
**max-over-elements statistics refused under probabilities**; fixed-metric Wald statistics
require listwise-complete reps and post missing rather than computing Q on a rep-specific
column subset.

**11. Pre-period statistics - one cross-cutting ruling.** Spans ES pre horizons,
`pretreatment_*`, sg pre rows, `joint_horiz`, `roth_prep`, `gtmax_pre` and the altmethod
forecast statistic. *For exclusion:* the distortion runs the dangerous way, and this is the
single most important statistical fact in the scan. wooldid defines the pre-window relative
to the *assigned* onset, so a cluster whose true onset is t0 and whose simulated onset is
t0+3 contributes genuinely post-treatment periods at relative times -3,-2,-1. The reference
distribution of every pre-period statistic is therefore inflated while the observed
statistic is not, so the p **rises with the true effect size**. The reference distribution
is a mixture - roughly half the reassignments move onset earlier (clean) and half later
(contaminated) - so the upper tail that sets the p is dominated by contaminated draws.
Shipping a conservative-in-the-wrong-direction number beside an honest analytic pre-trend
p invites users to prefer the permutation number *because* permute was requested.
*Against:* the pre rows are exact tests of the sharp null and are free; the joint null is
strictly *stronger* than "no effect", so a rejection is unambiguously informative; the
failure direction is conservative, never anti-conservative; and a principled restricted
form exists. *Verdict:* **author must rule once, for all families** - a build that permutes
pooled pre rows one way, sg pre rows another and `roth_prep` a third is incoherent. Three
honest options: (A) cover all pre rows with a single printed null-stating sentence
(recommended, since the alternative today is an unlabelled analytic number); (B) cover post
rows only, leaving pre rows analytic-with-a-label (the strictly safer choice); (C) gate pre
p's behind an acknowledgment. **Defer the not-yet-treated placebo loop in every case**: it
changes the estimation sample, so k, S and the estimand all move, the observed statistic
would have to be recomputed on a sample the observed run never used (breaking
observed-is-a-member), and `espostlength(-999)` at shell:127 trips only the *reported* post
window, not the estimation sample. It is a second engine, roughly the size of the original
permute build.

**12. `permtstat` + per-element coverage.** *For exclusion:* `permtstat` already drops
`fastnose` (110), so every rep pays a full VCE - compounding it with `subgroups()` produces
a blow-up with no warning. There is **no per-rep contrast SE to studentize with**: `e(V)`
is built as `diag(hadamard(col2, col2))` (3931) with the explicit comment "without giving
any information about covariance (which we did not estimate anyway)", and per-rep contrast
SEs live only in `r(table_vs)`, which is not harvested. Under `contreat` it is actively
destructive (see §3.5). And standardizing a *single* element by its own permutation SD is a
monotone transform yielding identical tail counts, so it buys nothing where available.
*Against:* for level elements it works and is already paid for, and studentization
genuinely matters *across* elements, since a raw max would be dominated by the widest
horizon. *Verdict:* **cover_restricted_form for the mechanics, but one piece needs the author.**
Settled by the code: level elements may be studentized; contrast and joint statistics stay
coefficient-based in *both* modes, and the help must say so, or `e(perm_stat)` - a single
run-level label at shell:152 - would silently mean two things in one run. **Not settled - the two
independent judgments of this family chose different defaults**: what `permtstat` should do with
`ame` elements, whose SE is scrubbed in the canonical permute design. The options are (a) the
`ame` elements fall back silently to the `b` statistic with one caution while the att elements
stay on `t`, producing a mixed-statistic run under a single `e(perm_stat)` label; (b) refuse
`permtstat` together with `ame` coverage; or (c) drop `ame` from the t-set and say so in the
display. (a) is this document's lean, but it is a ruling, not a code consequence. Studentize by the element's own *permutation* SD for any max statistic.
If `permtstat` ever goes per-element, the cell-collapse fallback rule must land first.

**13. `permexact` + per-element coverage.** *For exclusion:* the plain-share p (6061) is
valid only because the observed assignment is a member of a *complete* enumeration; with
per-element availability that becomes element-specific, and the existing incompleteness
caution (6023) would understate the problem by naming one run-level count while individual
elements rest on far smaller reference sets. Reporting an "exact" p that is not exact is
worse than reporting an approximate one. *Against:* the accounting is disclosable, and the
machinery exists on the bootstrap side (`wooldid_bootfin` per-column counts and
low-availability disclosure, 8491-8497). *Verdict:* **cover with mandatory per-element
accounting.** Element k uses the plain share only when `n_k == R`; the 6023 caution must
name *which* elements lost assignments; the 2/M floor argument holds per element.
**Judge correction embedded here:** five judges cited the `max(20, R/2)` refusal at 6016 as
a constraint on exact runs; it is gated `& "\`permexact'" == ""` (verified this session),
so exact enumerations are **already exempt**, exactly as the help states. Any per-element
floor design must preserve that exemption or it will newly break working `permexact` runs.

**14. `safety(off)` + `permute()` with a cluster finer than i.** A silent-wrongness path,
not a coverage question, and no judge surfaced it. *For exclusion (i.e. for adding a
refusal):* under default safety a sub-i cluster is self-limiting - the f13
onset-constant-within-i screen fires in every rep and the run dies loudly at the 6016
floor. But that screen is gated `if "\`safetyinput'" != "off"` (447), and the permute shell
forwards the post-validate working `safety` into the rep string at 130 (verified), so reps
inherit `safety(off)` too. Every replication then runs with cohorts silently split across
onsets and the reference distribution is built from mis-specified fits: silently *wrong*
randomization inference, not a loud failure. *Against:* `safety(off)` is the documented
escape hatch, and the configuration may be rare since permctl's nesting requirement points
users toward cluster >= i. *Verdict:* **add an upfront refusal** - permctl should require
the permute `cluster()` to nest i, refusing with a fixable message. "Rare" is not a defense
for a path producing a plausible-looking p from mis-specified fits. **Related latent parity
item:** shell:130 passes the post-validate working `safety` while shell:221-223 deliberately
passes `safetyinput` for bootstrap, with the comment that the working copy may have been
rewritten to `off` by `contreatwithin`. Harmless today only because `contreatwithin` is
refused with permute (5249-5251); bring 130 to parity while the fix is free.

**15. `esfixedbaseperiod` versus the pooled regime - the shell:126 `esprelength` defect.**
Not an exclusion candidate but a prerequisite and a scope boundary, recorded because two
judges disagreed about its severity. The defect: shell:126 emits the post-validate
*normalized* `esprelength` where the bootstrap shell deliberately emits `esprelengthraw`
(217, with the A8 comment that validate's esf normalization re-fires on every nested
re-parse). Verified this session: 4728 takes the raw copy, 126 emits the mutated one, 217
emits the raw one, and `esprelengthraw` is in validate's c_local blanket at 5036, so it is
in scope at 126. *Advocate's resolution of the judge conflict:* the ES judge called this "a
live defect - today's pooled permute p is already computed against mis-specified reps"; the
ptp judge said the opposite. **The ptp judge is right.** Under `esf` the worker generates
`_t` at 649 as `t - ttre` with no `esprelength` masking, so the rep's *regression
specification* is unchanged; the pooled branch *does* mask by `esprelength` (543-546) but
is never incremented. What actually drifts is `espre = min(min_t, esprelength)` (757) and
`horizonstarter` (772, pinned to 1 in the pooled regime by 773) - the *reported* pre-horizon set - and transitively the `gteffects`
window (1670/1672) and the espre-dependent adversarial cancellation checks (791, 795-800).
Corrected severity: **a reporting-window/alignment defect, not a live inference defect in
`e(perm_p_att)`** - but a hard prerequisite for any per-element ES harvest and for `gtmax`,
and it can still move pins *indirectly* if a wider reported window renders a horizon
inestimable and aborts a rep, changing `ndone` (plausible, unverified - no Stata slot).
Under `esf` + `timetrends` it compounds, since the worker applies a second increment at
645. *Verdict:* fix 126 for parity, treat as prerequisite plus cleanup, probe before and
after, and expect a pristine-twice per-config battery re-pin.

**16. Few clusters and the exact-p floor.** *For exclusion:* `permexact`'s p has a floor -
in symmetric designs the mirror-image assignment ties the observed estimate in absolute
value, so exact p's cannot go below 2/M (1/M under `permtstat`). With few clusters, no
permutation p can reach conventional significance regardless of effect size, so covering
*additional* statistics is pointless there - every one inherits the same floor. *Against:*
the floor is a property of the design, not of any statistic, and it is already computed,
documented and reported (`e(perm_npossible)`, the 10-percent caution, the single-assignment
refusal at 5897-5900). A disclosed floor is a correct answer to a limited design.
*Verdict:* **not an exclusion - record it as a coverage limit with one help sentence.** Two
operational consequences: the floor is per-element, so adding K elements does not add K
independent chances to reject; and it *strengthens* the case for a max-\|studentized\|
step-down, which pools power across horizons under exactly the null RI supports.

**17. `bootstrap(bayes)` + cell-weighted estimands - out of scope.** Named in the author's
combination list but it is a bootstrap-engine ruling: `bootstrap()` and `permute()` are
mutually exclusive, so no permute run can reach it, and the relevant screen already exists
and is adjudicated (refused under default safety, caution under `safety(off)`). *Verdict:*
**out of scope, routed to the bootstrap document.** Recorded here only so the ruling is
confirmed considered rather than dropped. Same routing for one bootstrap-side finding the
sg judge raised: `wooldid_bootjnt`'s allsgssame loop is
`foreach it in att att_it att_i att_itime ame ame_it` (5695), which omits `ame_i` and
`ame_itime` - but each item is additionally gated on the analytic scalar already existing
(5708). Since drafting, that gap has been **measured**: the analytic `ame` allsgssame scalars do
post on `continuous_control_dense` once `cluster()` is finer than the cohort-by-period cell
(`pb/VER_defect/d6.log` A1), and on the `bootstrap(bayes)` route the four `ame_i`/`ame_itime`
scalars then come back bit-identical to analytic beside rewired neighbours (A3). It is a
confirmed defect on that route; `bootstrap(standard)` is structurally barred from reaching it
(the nesting requirement at 5546-5551). The derived
recommendation is sound regardless and should be adopted for the permute build: derive the
item list from the `e(b)` column names rather than a literal list, and do not gate on the
analytic scalar existing.

---

## 5. Recommended scope and build order

**Governing principle - and this one is the author's to rule, not the assessment's.** The
recommendation below assumes that *labelling* is sufficient: permutation p's where they are
well-founded, an explicit provenance label everywhere else. The author asked for something
stronger - "*every* p-value that wooldid might produce is produced via a permutation test" - and
there is a route that meets it literally: **coverage by subtraction**, i.e. blank (or suppress)
every analytic p that permute does not cover, so that no un-permuted p is printed at all. It is
cheap (~10-20 lines per surface) and it is offered below for `gteffects` alone; generalised, it
would apply to the subgroup tables, the adversarial p's, the Roth block and the cfx bands. The
case against generalising it is that it deletes numbers users rely on today in exchange for
uniformity, and that some of those surfaces (plot bands, the V-only Roth columns) have no
permutation analogue to replace them with, so the table would simply lose columns. **That
trade - label everything, or blank what cannot be permuted - is ruling 0, and it governs all
seventeen below.**

What is *not* achievable under either choice is a literal reading in which every printed number
becomes a permutation quantity: `histogramcohorteffects` and `summarystats` carry no p-values at
all (margins run with `nose`, 1857-1885); most of the Roth block is V-only quantities
contrasting no null; plot bands are parameter CIs and RI produces only null bands; and
**permute produces no covariance**, so a row's se/t/CI can never become permutation
quantities the way bootstrap rewires them (6476-6493) - a p-only overwrite would leave
every row mixed-provenance *within* the row. The coherent contract is: permutation p's
where they are well-founded, an explicit provenance label everywhere else, and one caution
naming every analytic-only surface. What is not acceptable is the status quo, where
analytic p's print unlabelled under a "Randomization Inference" headline.

### Build order

**Tier 0 - ship regardless, no ruling needed (~15-25 lines, 1 commit, no numeric
movement).** One cross-engine note before the detail: the caution at 256 and the `sthlp` entry it
mirrors are a **shared edit surface**. The bootstrap companion proposes narrowing 256's clause
list; this document proposes a permute-side twin with a wider trigger and a three-category
structure. They are one edit with one owner - do not let the two documents produce conflicting
partial rewrites. The missing permute-side analytic-leftovers caution, **shell-side** inside the
permute block before `wooldid_display , part(perm)` at 175 - not a display-side tag,
because part(1)/part(2)/`wooldid_ptp` all print at 4237-4241 before `e(perm_method)` is
appended at 145. Name plots, saved plot data, `gteffects` cell SEs, histogram outputs,
adversarial line-test p's (`esf` runs only), the Roth block, the subgroup tables/contrasts/
equality tests, and any analytic joint tests. **Distinguish three categories in the
wording**, which the bootstrap's single line at 256 conflates, and give it a **wider trigger**
than 256's - the permute-side list must add `subgroups`, `pretrendpower` and `jointtests`, none
of which 256 gates on: (a) surfaces with analytic
p's that *could* be permuted; (b) analytic CIs that RI cannot replace in kind; (c) surfaces
with no inference at all. One caution saying "these are analytic" across all three invites
the reader to think a permutation version is pending for all of them. Extend the `sthlp`
`permute()` entry to match.

**Tier 1 - parity and safety fixes, own commit, before any harvest (~15-30 lines; forces a
battery re-pin, so isolate it).** (a) shell:126 emit `esprelengthraw` (ruling 15). (b) Port
the cell-collapse fallback rule (201-208, keyed on `e(wdid_fcguard)`) replacing the
unconditional pass-through at 119 - mandatory if subgroups-in-reps lands. (c) Bring
shell:130 to parity with shell:223 (ruling 14). (d) The permctl nesting refusal (ruling 14).
Probe pristine-twice with per-config isolation before re-pinning, per the adjudication
doctrine.

**Shared prerequisite for everything in Tiers 2-3, costed once (~200-250 lines, 2
commits).** Replace permrun's narrow accumulator (5910) **and the four-scalar comparator
(5833)** with a name-aligned R x K accumulator and a comparator *vector*, sharing
`wooldid_bootgrab` (8437-8464) verbatim plus bootrun's rep-0 name-set check (6186-6197);
convert whole-row blanking (6002-6010) to per-element completed counts. Regression test for
the commit: `e(perm_p_att)` must stay **bit-identical** for existing configs.

**Tier 2 - cover, subject to one ruling each.** (a) `ame`/`contreat` - the strongest case;
zero marginal per-rep cost; null already adjudicated. (b) ES **post** horizons and the sg
**level** rows (if ruling 2 goes cover). (c) `gteffects` `gtmax_post` - free, no fastnose
sacrifice, one-sided, shuffle-only. (d) `adv_line_att`/`adv_line_ame` via the V-free
residual SSR, `esf` runs only, under its own sharp-null heading.

**Tier 3 - the two decisions that govern everything else; rule these first.** Ruling 2
(subgroups in reps) and ruling 11 (pre-period statistics).

**Tier 4 - exclude.** cfx/lattice; `histogramcohorteffects` and `summarystats`; plot band
swaps; per-cell `gteffects` p's; `roth_slope`; per-rep VCE refits for any joint statistic.
`bootstrap(bayes)` + cell-weighted routes to the bootstrap document.

**The `jointtests` decision sits outside the tiers** because the status quo is neither
covered nor refused: unreachable on uniform-granularity designs under default safety, reachable
on a mixed-granularity one, and recommended twice
by the package's own error text. Pick exit A or B per ruling 1, and rule `roth_prep` with
it.

### Display contract

1. **Leave part(1) and part(2) numerically untouched.** RI supplies no covariance, so
   overwriting only the p column leaves se/t/CI analytic in the same row.
2. **Print one provenance line under part(1)** when permute is active: the table's
   se/t/p/CI are analytic; the permutation p-values are in the randomization block below.
3. **Extend `part(perm)` (7252-7271)** from its one-line-per-estimand loop into a real
   table: one row per element with name, estimate, permutation p, null band, SD, completed
   count - plus the sharp-null sentence and, for pre rows, the contamination sentence.
4. If the author instead wants the p column overwritten, the overwritten cells **must** be
   flagged in-table and the analytic table kept hidden, mirroring `e(V_es_analytic_*)`
   (6490). The alternative structural route - extend the deferral gate at 4237/4241 to
   `"\`bootstrap'\`permute'" == ""` and replay in the shell, the bootstrap pattern at
   253-254 - is cleaner but changes *when* results appear for every permute run.
5. **Everything analytic gets the Tier-0 caution**, in the three categories above.

### Accounting rule for elements missing in some reps

The current invariant is stated in the code: whole-row blanking exists "so ... nuse ==
ndone for every column" (5984-5985, verified), and permdist packing keys on column 1
(6072). It is load-bearing today and **unachievable for any ragged family** - at K ~ 10-50
one inestimable horizon or one vanished subgroup level would discard the rep for
`main_att` too, silently moving the shipped `e(perm_p_att)`. The replacement:

1. Keep a **row-level** "rep estimated at all" counter for the `max(20, R/2)` refusal
   (6016) and the `permexact` incompleteness caution (6023) - **preserving 6016's existing
   `permexact` exemption** (verified this session; five judges missed it).
2. Add **per-element counts `n_k`**, computed the way `wooldid_bootfin` already does
   (8491-8500); the existing summary frame at 6033-6041 already drops missings and counts
   per column, so the blanking is the only thing forcing `nuse == ndone`.
3. **Per-element floors**: an element with `n_k < max(20, R/2)` posts missing with a named
   caution rather than a p; the run-level refusal survives only if *no* element clears.
4. **Pack `permdist` on the any-non-missing rule** (`wooldid_bootfin`'s convention at
   8487), not on column 1.
5. **Disclose `n_k`** in `e()` and caution when an element's completion rate falls
   materially below the run's; under `permexact`, name which elements lost assignments.
6. **Joint/Wald statistics take reps listwise-complete only** - computing Q on a
   rep-specific column subset makes Q a different function per rep.
7. **Regression test**: `e(perm_p_att)` bit-identical on every existing config.

### `e()` naming proposal

The judges split. The ES judge proposed per-family matrices `e(perm_es_<est>)` with
relyr-anchored rownames; the sg and ame judges converged on **1 x K row vectors
column-named with the `e(b)` stripe**. Recommend the stripe form:

- `e(perm_p)`, `e(perm_lb)`, `e(perm_ub)`, `e(perm_sd)`, `e(perm_n)` - each 1 x K with
  `colnames` equal to `e(b)`'s, so every element is addressable without inventing a name
  per element and the vectors merge against `e(b)`/`e(b_es_*)` by name (the key the
  bootstrap's `V_es` rewire already uses at 6478).
- Keep the existing `e(perm_p/lb/ub/sd_<est>)` scalars untouched as documented aliases, and
  add `e(perm_p_ame)` etc. as their exact parallels.
- `e(permdist)` becomes R x K with the same stripe; the shell:169 guard should become
  **K-aware** (cap on R x K, not R alone).
- Family-specific additions: `e(perm_joint_stat)` / `e(perm_joint_n_<name>)`;
  `e(perm_p_adv_<typetext>)`, `e(perm_adv_stat)`, `e(perm_adv_regime)`;
  `e(perm_p_gtmax_post/_pre)`, `e(perm_gtstat)`.

Record the per-family matrix form as the alternative: it is more legible in a printed
table, at the cost of a second naming convention.

### Multiplicity

Offer the max-\|studentized\| step-down **once, across the ES and joint families**, not per
family. Because the joint null *is* the sharp null, a Westfall-Young / Romano-Wolf
step-down is exactly valid with no independence assumption and is nearly free from the same
accumulator - and it delivers the joint p as its first step, which is why those two
families should be scoped together rather than sequentially. Two constraints: studentize by
each element's own permutation SD, or the widest horizon dominates; and a max over **pre**
horizons inherits the contamination, so only a max over post horizons is clean.

### Cost compound

`subgroups` x `permtstat` x R is the blow-up nobody prices as a product: `permtstat` drops
`fastnose` so every rep pays a full VCE; `subgroups()` adds ~(1 + 0.6G) margins passes per
estimand; `contreat` doubles the margins work again; `permexact` multiplies R by the
enumeration size (369,600 possible assignments on the probe design). Any two are tolerable;
three compound multiplicatively with no warning and no progress indicator - `permrun` has
no `dots` option, unlike `bootrun`. Measure the per-rep subgroups multiplier before
defaulting it on, and consider a permctl-time cost caution when these co-occur.

---

## 6. Open questions

**Governing rulings (decide first; each governs multiple families).**

1. **Subgroups in replications** (ruling 2): cover, cover-levels-only, or decline? One
   ruling for four families, gating two more. If covered, are the two prerequisites (the
   cell-collapse fallback rule and the full-population sg screen) approved?
2. **Pre-period statistics** (ruling 11): cover-all-with-caption, post-only, or gate behind
   an acknowledgment? Governs seven items across five families.
3. **`jointtests`** (ruling 1): exit A (refuse and rewrite 1257/1261) or exit B (carve-out
   plus gate suppression plus permutation p's under the existing names)? Rule `roth_prep`
   with it.
4. **Display contract** (§5): append-only with a provenance line (recommended), or overwrite
   the p column with the analytic values kept hidden, or defer part(1)/part(2) and replay
   in the shell?
5. **`e()` naming** (§5): the 1 x K stripe-aligned vectors (recommended) or per-family
   matrices?

**Per-family thresholds and choices.**

6. `ame`: what blanked-draw threshold triggers withholding rather than cautioning
   (suggested: caution always, withhold above ~10 percent)?
7. `gteffects`: coverage by computation (`gtmax`) or coverage by subtraction (blank the
   analytic se/p/CI columns, ~10-20 lines)?
8. Adversarial: should the printed "on a single line" p (7064/7203) stay analytic with a
   caution, be replaced, or should both appear? Replacement is defensible only via
   impose-and-permute, which is a separate build. Is the unit-level homogeneous linear
   effect it imposes an acceptable assumption?
9. Adversarial: should `altmethod_adversarial`'s forecast statistic get an RI p while the
   primary composite statistic does not? `||d||` is monotone in effect magnitude; the
   primary W is not.
10. Joint: preferred metric - the permutation covariance of the reps' own pre-block
    (recommended; needs no analytic V and works under `amescrub`) or
    `invsym(e(V_es_<fam>))`? And on guard failure, degrade to max-|t| or post missing?
11. Joint / sg: if per-subgroup joint tests are covered, should `permstrata(<sg var>)` be
    recommended or **required**? Without it the reference distribution randomizes the
    subgroup-onset composition itself.
12. Multiplicity: is the max-\|studentized\| step-down in scope now, or listed for later?
13. `permute(probabilities)` (ruling 10): confirm the split - single-element statistics
    allowed with disclosure, max-over-elements statistics refused.
14. Should the restricted not-yet-treated placebo loop be scheduled at all as a future
    `permplacebo`-style option? It doubles the rep budget and is a second engine.

**Prerequisites and cleanups needing approval.**

15. shell:126 `esprelengthraw` (ruling 15) and the battery re-pin it forces. Do you want it
    probed before the change?
16. The permctl nesting refusal and the shell:130 `safetyinput` parity fix (ruling 14).
17. Should `makecfxplots` + `permute()` be refused under default safety, or merely named in
    the caution?
18. Formal graduation: `f6_permutation_plan.md` §2 item 7 and R4 (279-280) record the v1
    exclusions, and `wooldid.sthlp` states them to users. Please record the graduation
    ruling explicitly so the plan document and the `permute()` help entry can be updated in
    the same wave. **Note what R4's list actually names**, because three of this document's own
    Tier-2 recommendations are on it verbatim: "sg-specific perm stats, one-sided p's,
    per-horizon perm bands are all OUT of v1". The one-sided tail counts proposed for the joint
    statistic (§3.2), for `gtmax` (§3.6) and for the Roth prep test (§3.7), and the per-horizon
    null bands proposed in §3.1, are therefore graduation *requests* against the author's own
    recorded ruling, not free additions.

---

## 7. Artifacts

The scratch do-files and logs the reconciler cited (session scratch, not committed; all
under `pb/X/`):

- **(a) `chk_a.do/.log`**: `att att_it cluster(cohort) esprelength(2) espostlength(2)
  subgroups(sg_three) permute(shuffle) permreps(30)`, rc 0. Only
  `e(perm_p_att)` = `e(perm_p_att_it)` = .032258 (log 246-253, RI block 226-229) are
  permutation p's; the main table (448-465), three subgroup tables, `e(pwsgs)` (510-554,
  degenerate se 2.33e-10 p 6.5e-87) and four allsgssame p's (306-309, '.') all print
  analytic p's with no caution or label; `e(perm_p_pretreatment_att)` absent (499);
  `e(permdist)` 30x2, `e(perm_comparator)` 1x2; `e(b)` 1x48 with every SG/ES column named
  (428-438).
- **(b1) `chk_b.log` B1**: adding `esfixedbaseperiod jointtests adversarial` to (a)'s
  `cluster(cohort)` under default safety gives the SI-A1 refusal r(198) (60-73), whose text
  recommends "wooldid's built-in randomization inference (the permute() option)" for
  joint-type inference.
- **(b2) `chk_b.log` B2**: same with `cluster(person_id)` - analytic joint p's 0.0000 and
  adversarial "single line: 0.0000" print unlabelled (136-144), then 0 of 30 reps estimable
  -> r(459) (147-153; the worker requires onset constant within i, so person-level
  reassignment always fails).
- **(b3) `chk_b.log` B3**: `adversarial esfixedbaseperiod cluster(cohort) permute` (not
  SI-A1 gated), rc 0; adversarial composite-null covariance degenerate on this data so no p
  (162-165); `e(adversarial_line_att/att_it)` posted analytic (296-299); RI block att/att_it
  only (221-229).
- **(b4) `chk_b.log` B4**: `jointtests adversarial pretrendpower(80) safety(off)
  permute(shuffle) permreps(30)`, rc 0: `e(joint_horiz)` = `e(joint_horiz_i)` = 1.2e-24
  printed as 0.0000 (380-383, 464-465); Roth block prints analytic slope se .0013738 and
  joint p 0.0000 (386-428, `e(roth_prep_att)` 470); nothing is labelled analytic; the
  `safety(off)` caution (315-322) again recommends `permute()` while permute is on.
- **(c) `chk_c.do/.log`**: the bootstrap spec `att att_it cluster(cohort) esprelength(2)
  espostlength(2) subgroups(sg_three) jointtests adversarial esfixedbaseperiod gteffects
  histogramcohorteffects bootstrap(standard) reps(30) rseed(1)` - reps complete then r(198)
  "invalid syntax" (1085-1089); partial `e()` at 1096-1205.
- **(c-root) `chk_c2.do/.log`**: `local x = +0.012688b70e62bX-3ff` (=
  `strofreal(1e-310,"%21x")`) gives r(198) invalid syntax, while zero and normal-range hex
  literals parse (62-89); reproduced in wooldid at 121-134 (the s3 p is denormal); the same
  call without `subgroups()` completes (143-221, bootstrap block 204-213).
- **(c') `chk_c3.do/.log`**: (c) without `esfixedbaseperiod`, analytic run with
  `safety(off)` as comparator, rc 0. `mreldif`: `e(V)` vs analytic .0207,
  `e(V_analytic)` 0; `e(V_es_att)` .0043; `e(wdidmainresults)` .90, `e(wdidsgresults_1)`
  .95 (rewired); `e(pwsgs)` rewired (se .165 p .78 vs analytic 2e-10); `e(gteffects)` 0 and
  `e(histogramestimates)` 0 (analytic/unchanged); `e(adversarial_line_att)` 0 (analytic);
  `e(joint_att_allsgssame)` .691 / pre .762 vs analytic '.'; caution 256 printed
  (1361-1364).
- **(d1) `chk_d.log` D1**: `att cluster(cohort) esprelength(2) espostlength(2)
  contreat(dose_vary) contreatcontrols(basic) permute(shuffle) permreps(30)`, rc 0:
  `e(perm_p_att)` = .0323, `e(perm_p_ame)` **absent** (305), RI block lists att only
  (226-233), `e(wdidmainresults)` identical to analytic (`mreldif` 0, 315) - ame rows
  analytic and unlabelled; the AMESE caution recommends "randomization inference" inside
  the permute run (156-160).
- **(d2) `chk_d.log` D2**: same with `bootstrap(standard) reps(30) rseed(1)`, rc 0:
  `main:ame` se .199 p .138, `pretreatment:ame` and `ES_ame` rows all rewired (337-397);
  `e(V_es_ame)` **absent** in the bootstrap `e()` (413-475) - the analytic block was
  withheld under cell containment (325-334) so bootwire 6484-6485 skips it.
- **(d3) `chk_d.log` D3**: `att att_i contreat(dose_vary) contreatcontrols(basic)
  subgroups(sg_three)`, analytic vs bootstrap: allsgssame att .691 (analytic 1.5e-81),
  att_i .699 (analytic '.'), pre_att .762; on that `cluster(cohort)` run `e(joint_ame*_allsgssame)`
did not post
  analytically on this data, so the 5695 `ame_i`/`ame_itime` omission is code-confirmed
  only; pwsgs `ame_i` rows rewired (1505-1518).

Design documents referenced: `review_process/f6_permutation_plan.md` (§2 item 7 and R4
~279-280 record the v1 scope choices), `review_process/bootstrap_plan.md` and
`review_process/bootstrap_report_2026-08-30.md`. Help-file entries:
`wooldidinstall/wooldid.sthlp` ~363-412 (`permute()` and its modifiers) and ~956-962
(stored results).
