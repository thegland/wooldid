# bootstrap(): what is covered, what is not, and what should be (2026-09-02)

Assessment requested by the author: "I recently found out that permute only focused on
permutation testing the high-level statistics, like the ATT, but not the event study
horizons or anything else. We want to make sure that when permute is enabled, *every*
p-value that wooldid might produce is produced via a permutation test. That includes
inference on event study estimates, adversarial test p-values, joint test p-values, etc.
-- everything! Please do a scan to document all parameters that can in principle be
tested via permutation testing but which are not right now. If any would be difficult or
problematic to approach with permutation testing, please let me know and make the case
for this - I am open to some exclusions ... Please surface these possible exclusion cases
to me. Please also do the same thing for bootstrapping." This is the **bootstrap** half.
Assessment only; nothing in the ado changes here. Companion to the permute scan.
`bootstrap_plan.md` and `bootstrap_report_2026-08-30.md` are the build record this
amends.

**Headline.** `bootstrap()` already covers far more than `permute()` does - the main,
pretreatment, event-study, subgroup, ame and pairwise-contrast tables, `e(V)`,
`e(V_es_*)`, the joint pre-trend tests and the subgroup-equality chains are all rewired
from the replication accumulator - but of the **21 reported p-value surfaces** in the
inventory, **12 are bootstrap-based, 6 still ship an ANALYTIC p inside a completed
bootstrap run** (the three adversarial line tests, `gteffects` cell p's, the cfx lattice
bands - all behind the caution at `wooldid.ado:256` - plus `e(joint_ame_i_allsgssame)` /
`e(joint_ame_itime_allsgssame)`, which are silent), 2 are refused at parse
(`pretrendpower()`, `impvalues`) and 1 is unreachable (`permute()` is mutually refused).
Ahead of every coverage question sits a **defect**: on the default `bootdf(normal)` path
a joint-test p in the denormal band `(4.9e-324, 2.2e-308)` is handed from Mata to Stata
as `%21x` (8676-8677), Stata's parser rejects the zero-leading-digit hex literal with
r(198), and the user's whole command aborts with `e(pwsgs)` left carrying **analytic**
se/t/p under a bootstrap header - reproduced on the suite's own
`continuous_control_dense` with ordinary options (`chk_c.log 1085-1089`). Recommended
scope: fix six defects, then cover four cheap items (event-study plot bands, the
adversarial line tests, `e(V_es_ame*)`, `ame_i`/`ame_itime` chains) at zero
per-replication cost; permanently exclude `gteffects` cell SEs, the auto-lattice cfx
bands and `histogramcohorteffects`; make `pretrendpower()` conditional on three
preconditions. **Nine rulings are needed from the author**: whether to label every
un-bootstrapped p or blank it (section 5.5 - this one governs the rest); whether to lift the
standing `pretrendpower()` refusal (entry 19); mirror-vs-structural gating; the reference
distribution (entry 18, the only one with no recommended default); the Roth slope doctrine; the
`bayes` + cell-weighted marker; a `jointtests` cluster floor; the `timetrends` fill
(probe-needed); and fixed-lattice cfx coverage.

**A caveat on evidence provenance.** Sections 2-4 were produced by family judges holding the
Stata slots; this document's drafter had none, so line citations were read from the ado while
writing and the empirical claims are carried forward from the judges. A separate verification
pass afterwards re-read every load-bearing citation and re-ran the defect claims in Stata: the
denormal crash, the withheld `e(V_es_ame)`, the `ame_i`/`ame_itime` omission and the unmarked
analytic p-values all reproduced (logs under `pb/VER_defect/`), and roughly twenty citation
pointers were corrected in place. What is still **unverified in Stata** is listed in section 6.2;
four of those items sit under recommendations ranked MUST or CHEAP, so read them before building.

---

## 1. What `bootstrap()` covers today

The engine is a shell block at `wooldid.ado:186-258`, inside `prog def wooldid` (19).

**Rep option string** (`wdb_ro`, 193-227). KEPT: `att att_it att_itime att_i`,
`esfixedbaseperiod`, `timetrends`, `poisson`, `poisexpresults`, `robust`, `ccnoabsorb`,
`latticeignoreweights`, `emptycellsoverride`, `contreatwithin`, `fastnose`, plus
`cluster()`, `fe()`, `coarsecohortcontrols()`, `ccc_absorb()`, `controls()`,
`interactivecontrols()`, `customweightmultiplier()`, `fastcollapse`/`fullprecisionslow`,
`contreat()`/`contreatpoly()`, `contreatcontrols()`, `lattice()`, `semielasticity`,
`contreatelasticitytype()`, **`subgroups()`** (214), `esprelength()`/`espostlength()`
(the raw user value, A8), `esrelativeto()`, `safety()`, `regtol()`,
`oostreatedcontrols()`, `anticipation()`. STRIPPED: everything display- or
diagnostic-shaped - `makeplots`, `saveplots`, `makecfxplots`, `gteffects`,
`histogramcohorteffects`, `jointtests`, `adversarial`, `pretrendpower()`, `summarystats`
- and the bootstrap options themselves. Reps run under `fastnose`, so only `e(b)` is
consumed.

**Ordering.** `est sto` the observed run and capture `e(b)` (228-229) → `wdb_reid`
resolved (230: `bootstrap(standard)` **and** any of `att_it`/`att_i`/`att_itime`) →
`cap noi wooldid_bootrun` (231, the only wrapped call) → on a nonzero rc the shell
`exit`s (232-236) → `wooldid_bootwire` (248, **bare**) → `wooldid_bootjnt` (250,
**bare**) → `wooldid_display , part(1)` / `part(2)` / `part(boot)` (253-255) → the
analytic-provenance caution (256). The observed full run with analytic SEs completes
*before* the loop; the display is deferred because the worker skips `part(1)`/`part(2)`
under `bootstrap` (4237, 4241) and the shell re-calls them after the rewiring, so the
existing display code prints bootstrap numbers with no formatting change (amendment A14,
5594-5596).

**Statistic.** `wooldid_bootrun` (6094-6324) accumulates the full name-aligned `e(b)` per
replication via `wooldid_bootgrab` (8437-8463). `wooldid_bootwire` (6326-6549) forms
listwise complete-case moments and reposts `e(V)` (analytic hidden as `e(V_analytic)`,
6462), cols 2-6 of `e(wdidmainresults)` (6396-6412) and `e(wdidsgresults_<sg>)`
(6428-6452), and `e(V_es_<fam>)` via `wooldid_bootves` (6478-6493; analytic hidden as
`e(V_es_analytic_<fam>)`, 6490). `wooldid_bootjnt` (5588-5826) rewires the joint
pre-trend tests, the subgroup-equality chains, the pairwise contrasts `e(pwsgs)` and the
CI flavors, every inverted quadratic form going through `wooldid_bootjw` (8613-8677) with
a listwise covariance and a scale-free PSD guard.

**`e()` surface added.** `e(V_boot)`, `e(b_boot)`, `e(bias_boot)`, `e(bootdist)`,
`e(boot_ci_normal)`, `e(boot_ci_percentile)`, `e(boot_ci_bc)`, `e(boot_z0)`,
`e(boot_navail)`, `e(boot_N_V)`, `e(boot_v_posted)`, `e(boot_method)`, `e(boot_reps)`,
`e(boot_completed)`, `e(boot_failed)`, `e(boot_level)`, `e(boot_stat)`, `e(boot_seed)`,
`e(boot_df)`, hidden `e(V_analytic)`, `e(V_es_analytic_*)`, `e(wdid_vcespec)`. `e(df_r)`
is overwritten **only** when `bootdf()` resolves finite (6511-6516).

**Refusals** (`wooldid_bootctl`, 5399-5585): `permute()`, `impvalues`, `fastnose`,
`unconditionalse`, `customvce()`, `oldsyntax`, `suppresssgprimaryeffects`, multiway
`cluster()`, and `pretrendpower()` (5578-5583). `reps()` has no floor by explicit author
ruling (5430-5437). The `cleanmatrices` refusal was deliberately removed (5569-5577).

**What the replication distribution is, and is not, used for** (verified 2026-09-02). The
accumulator is `matrix accname = J(reps, K)` over the observed column set - **coefficients only**,
one row per replication. From it wooldid builds the bootstrap covariance and, separately, three
confidence-interval flavours in Mata `wooldid_bootci` (5890): `e(boot_ci_normal)` (obs +/- z*se),
`e(boot_ci_percentile)` (order statistics of the replication distribution) and `e(boot_ci_bc)`
(bias-corrected order statistics, with `e(boot_z0)` and the `ties` adjustment). Two things follow
that are easy to assume otherwise. First, **there is no percentile-rank p-value anywhere**, and no
standard error obtained by inverting a percentile interval: the reported table's se is the
bootstrap standard deviation and its p comes from the Wald/`bootdf()` reference (6480-6486), so the
percentile and BC objects are side outputs that nothing in the printed table or in `e(V_boot)`
consumes. (`e(bootdist)` is posted when reps <= 10000, so a user can compute a percentile-rank p by
hand.) Second, **there is no studentized (percentile-t) bootstrap**: the rep option string ends
with a hardcoded `fastnose` (193, with the comment "reps run under fastnose (only e(b) is
consumed)"), so no per-replication standard error exists to studentize with. The permutation engine
*does* have this switch - `permtstat` changes the statistic to `_b/_se` and drops `fastnose` from
its rep string for exactly that reason - so `permtstat` has no bootstrap twin. Both gaps are
statistic-choice gaps rather than coverage gaps, which is why the 21-row inventory below does not
show them.

**Declared exceptions.** The caution at 256 - option-gated on
`makeplots`/`makecfxplots`/`saveplots`/`gteffects`/`histogramcohorteffects`/`adversarial`/
`altmethod_adversarial`, so it fires only when one of those is requested - names "plots, saved
plot data, gteffects cell standard errors, histogram outputs and adversarial line-test p-values"
as ANALYTIC. One of those clauses is **factually wrong**: all four histogram margins calls run
`nose` (1866, 1871, 1896, 1903) and `e(histogramestimates)` (4175) carries point estimates only -
no covariance of any flavor enters it, so histogram output does not belong in a list of things
produced from the analytic covariance. The adversarial clause is correct as written, with one
nuance worth keeping: `e(adversarial_line_*)` is itself an unweighted OLS of coefficients on
relative time (7354-7368) and so is V-free; it is the adversarial *p-values* that are analytic.

---

## 2. Coverage matrix

Status key: **rewired** = bootstrap-based today; **analytic (cautioned)** = analytic p
shipped inside a completed bootstrap run, named at 256; **analytic (silent)** = analytic
value survives with no caution and no marker; **refused** = combination errors at parse;
**n/a** = no inference object.

| Item | Options required | Computed at | Analytic reference | Bootstrap status | Evidence |
|---|---|---|---|---|---|
| `main_att` — pooled `main:<est>` rows (att att_it att_i att_itime): se/t/p/CI | none (`att` default) | 1985-2020 (`r(table)` 1991) | `margins r(table)`, `ttail(e(df_r))` | **rewired** | bootwire rewires cols 2-6 of `e(wdidmainresults)` + `e(V)`; `chk_c3.log 1525` mreldif .90, analytic hidden at 1505 |
| `perm_p_summary` — `e(perm_p/lb/ub/sd_*)`, `e(permdist)` | `permute()` `permreps()` | 6026-6068; shell 145-174 | add-one two-sided perm p | **n/a** | `bootstrap()`+`permute()` mutually refused (5488-5491) |
| `pre_att` — `pretreatment:<est>` rows | none | 2037-2070 (`r(table)` 2043) | `margins r(table)` analytic | **rewired** | `chk_c3.log`: pretreatment se .0237 p .166 (boot) vs analytic .0148/.0479 |
| `es_horizon_att` — `ES_<est>_pre#/post#` rows | none | 2100-2120 (`r(table)` 2106; V 2116) | `margins over(_t) r(table)` | **rewired** | rows at `chk_c3.log 1525`; `e(V_es_att)` mreldif .0043 (1509), analytic hidden 1513 |
| `es_plots_ci` — ES plot bands + saved `_ES.dta` | `makeplots` / `saveplots` | plotting block consumes the CI cols from 2100-2120 | analytic `r(table)` lb/ub | **analytic (cautioned)** | caution 256 names "plots, saved plot data" (`chk_c3.log 1361-1364`) |
| `ves_att` — `e(V_es_*)`/`e(b_es_*)` blocks (honestdid entry), incl. `_s#` | none; `subgroups()` adds `_s#` | 3990-4030 (captures 2116/2570/2782/3234) | analytic `e(V)` subblock | **rewired** | 6476-6493 via `wooldid_bootves`; `chk_c3.log 1509-1521`, `chk_d.log 553` |
| `ves_withheld` — `e(V_es_*)` where the analytic block was withheld/downgraded | `contreat()` with cluster at/above cell, or degenerate horizons | 4012-4030; bootwire 6484-6485 `cap confirm` / `continue` | analytic screen decides existence | **analytic (silent)** | re-verified in `VER_defect/d4.log` V2 (`e(V_es_att)` present, `e(V_es_ame)` **absent**, `e(b_es_ame)` present, and no `e(V_es_analytic_ame)` either - the object has no hidden analytic twin and no caution names it). `chk_d.log 413-475`: D2 posts `e(b_es_ame)` but **no** `e(V_es_ame)` though `e(boot_navail)=30` on every `ES_ame` column; no caution |
| `joint_pre_att` — `e(joint_horiz)`, `e(joint_horiz_i)` | `jointtests` (SI-A1 gated, 1232-1263) | 2205-2220 (`r(p)` 2218) | `test` after `margins over(_t)` | **rewired** | bootjnt 5648-5667; `chk_c3.log 1555` 3.06e-20 vs analytic 2.95e-14. **DEFECT**: denormal `%21x` → r(198) aborts the run (`chk_c.log 1085-1205`) |
| `joint_pre_ame` — `e(joint_Choriz)`, `e(joint_Choriz_i)` | `jointtests` `contreat()` | 2650-2666 (`r(p)` 2664) | `test` after `margins dydx over(_t)` | **rewired** | bootjnt 5649 loops `ame ame_it`; same denormal exposure |
| `joint_pre_sg` — `e(joint_horiz_s#)` etc. (ereturn **locals**, 4137-4140) | `jointtests` `subgroups()` | 2965-2998 (2996) att; 3275-3296 (3294) ame | `test` after `margins subpop over(_t)` | **rewired** | bootjnt 5669-5687; `chk_c3.log 1563` 5.57e-06 vs analytic 1.35e-14. Denormal crash reproduced here: `chk_c.log 1140-1141` s1=`0`, s2=`1.4175550282e-156`, s3 missing, abort |
| `adv_line_att` — `e(adversarial_pval_att/att_it)`, `e(adversarial_line_*)` | `adversarial` (or `altmethod_adversarial`) + `esfixedbaseperiod` | 2140-2160 (2153) via `wooldid_advx` 8388-8389; posted 4110-4135 | Mata GLS line fit + Wald on analytic V | **analytic (cautioned)** | caution 256; `chk_c3.log 1544` line identical to analytic (correct - the line is V-free) |
| `adv_line_ame` — ame/ame_it twin | `adversarial esfixedbaseperiod contreat()` | 2595-2602 (2600) | `wooldid_advx` on analytic V | **analytic (cautioned)** | caution 256 |
| `adv_line_sg` — `e(adversarial_pval_<est>_s#)` (**scalars**, 4150) | `adversarial esfixedbaseperiod subgroups()` | 2951 (att), 3260 (ame) | `wooldid_advx` on analytic V | **analytic (cautioned)** | caution 256; `chk_c3.log 1552` `e(adversarial_pval_att_s1)` unchanged |
| `sg_main_att` — subgroup tables `main_SG#_*`, `ES_*_SG#_*`, `e(wdidsgresults_#)` | `subgroups()` | 2700-2795 (2707/2731/2772/2793; V 2782) | `margins subpop r(table)` | **rewired** | `subgroups` kept in `wdb_ro` (214); `chk_c3.log 1529` mreldif .95; `e(V_es_att_s1)` 1517 |
| `sg_ame` — subgroup ame rows | `subgroups() contreat()` | 3165-3220 (3170/3189/3216) | `margins dydx subpop r(table)` | **rewired** | `chk_d.log` D3 sg tables rewired |
| `pwsgs` — pairwise subgroup contrasts, cols 2-6 | `subgroups()` >1 level | 3039/3044/3114/3119 (att); 3324/3329/3368/3373 (ame) | `margins pwcompare r(table_vs)` | **rewired** | bootjnt 5724-5802 (5758-5759); `chk_c3.log 1532` mreldif .22, se .165 p .78 vs analytic degenerate |
| `allsgssame_att` — `e(joint_<est>_allsgssame)`, `e(joint_pre_<est>_allsgssame)` | `subgroups()` >1 level | 3060-3080, 3137-3155; posted 4160-4163 | `test` on analytic VCE chain | **rewired** | bootjnt 5694-5720, gated on the analytic scalar (5708); `chk_c3.log 1571/1575` .691/.762 vs analytic `.`; `chk_d.log 1543` att_i .699 |
| `allsgssame_ame` — ame / ame_it twins | `subgroups() contreat()` | 3343-3354 / 3389-3400; posted 4164-4169 | `test` on analytic VCE | **rewired** | 5695 list includes `ame ame_it`; the analytic ame scalars post on `continuous_control_dense` once `cluster()` is finer than the cohort-by-period cell (`cluster(person_id)`: `VER_defect/d6.log` A1 gives `joint_ame_allsgssame` .4109968435418874), and the bootstrap then rewires that item (A3: .4313975428385769). Under `cluster(cohort)` the AMESE containment screen scrubs the chain and nothing posts (`chk_d.log 1031-1035`) |
| `allsgssame_ame_i` — `ame_i` / `ame_itime` twins | `subgroups() contreat()` with `att_i`/`att_itime` | 4164-4169 (`subinstr("att","ame")` over att att_it att_i att_itime) | `test` on analytic VCE | **analytic (silent)** | bootjnt 5695 omits `ame_i ame_itime`; the analytic scalar survives beside rewired neighbours with no caution. Code-confirmed only (`chk_d.log` D3) |
| `ame_main` — `main:ame`, `pretreatment:ame`, `ES_ame_*` | `contreat()` | 2245-2310 (2249/2302), 2552 (V 2570) | `margins dydx r(table)` | **rewired** | `chk_d.log 337-397`: `main:ame` se .199 p .138; but see `ves_withheld` |
| `cfx_lattice` — lattice cell estimate/se/p/CI, `e(cfx*)` (4191-4194) | `contreat() makecfxplots` | 2405-2520 (2411/2416/2440/2445/2490/2516) | `margins r(table)` analytic | **analytic (cautioned)** | caution 256 ("plots, saved plot data") |
| `gteffects_cell_se` — `e(gteffects)` cols 6-9 | `gteffects` | `wooldid_gtx` 7818-7823; posted 4034-4037; refused with `contreat()` at 4422-4424 | analytic `e(V)` diagonal, **stage-1** `df_r` (read at 7805; the stage-1 window is 1646/1668/1673) | **analytic (cautioned)** | caution 256; `chk_c3.log 1536` mreldif 0 vs analytic |
| `hist_cohort` — `e(histogramestimates)`, histogram graphs | `histogramcohorteffects` | 1857-1920 (`margins ... nose`, 1866/1871/1896/1903); posted 4175 | **none** (no VCE) | **n/a** | caution 256 lists it, wrongly; `chk_c3.log 1540` identical |
| `roth_prep` — `e(roth_prep_<est>)` joint pre-trend F/chi2 | `pretrendpower() esfixedbaseperiod` | `wooldid_ptp` 6724-6725 | Mata GLS on analytic `e(b_es_*)`/`e(V_es_*)` | **refused** | 5582 "bootstrap() may not be combined with pretrendpower()" |
| `roth_slope` — slope se, CI, `e(roth_nissize_*)` | same | 6610-6611 + surrounding GLS | analytic `e(V_es_*)` | **refused** | 5582 |
| `roth_power` — `e(roth_power_*)` table | same | 6593-6870 | analytic `e(V_es_*)`, noncentral quadrature | **refused** | 5582 |
| `im_rows` — Ibragimov-Mueller cross-cluster t rows | `impvalues` | 2015-2017, 2064-2066, 2177-2179, 2274-2276, 2327-2329, 2623-2625 | `ttest r(se)`/`r(p)` | **refused** | 5493 |
| `boot_ci_flavors` — `e(boot_ci_normal/percentile/bc)`, `e(boot_z0)`, `e(bootdist)` | `bootstrap()` `reps()` | 5808-5819; `wooldid_bootci` 8702-8751 | replication matrix | **rewired** | by construction; **not posted** when the run aborts in bootjnt (`chk_c.log 1151-1205`) |
| `e_V` — `e(V)`, `e(df_r)` | none | worker post ~3960-3980; bootwire 6326-6500 | analytic VCE | **rewired** | listwise boot V in `e(V)`, analytic hidden (`chk_c3.log 1501-1505`); `e(df_r)=11` stays posted while `e(boot_stat)='z'` (`chk_d.log 419-450`) |
| `sescan_caution` — identical-SE screen (7531), placeholder disclosure (3969/3975), containment caution (7443) | none (data-triggered) | 7424-7540, 3960-3980 | analytic `e(V)` diagonal | **analytic (silent)** | printed *inside* bootstrap runs about SEs that are then replaced: `chk_c3.log 920-925`, `chk_d.log 325-334` precede a rewired table with no reconciliation |
| `summarystats_desc` — `e(escomp*)`, `e(N)`, `e(Nclus)` | `summarystats` | worker | descriptive counts | **n/a** | — |

**Counting rule for the headline.** 31 inventory items; 21 carry `is_pvalue`. Of those:
12 rewired (`main_att`, `pre_att`, `es_horizon_att`, `joint_pre_att`, `joint_pre_ame`,
`joint_pre_sg`, `sg_main_att`, `sg_ame`, `pwsgs`, `allsgssame_att`, `allsgssame_ame`,
`ame_main`); 6 analytic inside a completed bootstrap run (`adv_line_att`, `adv_line_ame`,
`adv_line_sg`, `gteffects_cell_se`, `cfx_lattice`, `allsgssame_ame_i`); 2 refused at
parse (`roth_prep` via `pretrendpower()`, `im_rows` via `impvalues`); 1 unreachable
(`perm_p_summary`).

---

## 3. Family-by-family feasibility

Mechanics only. The case-for / case-against / verdict for each ruling is in §4.

### 3.1 `adv` — adversarial line tests (`adv_line_att`, `adv_line_ame`, `adv_line_sg`)

**Statistic.** A plug-in recomputation of the *same* statistic with the bootstrap
covariance substituted - not a percentile interval (the object is a composite-null p) and
not a Wald on `e(V_boot)` (the restriction is not linear in `b`: the GLS weighting depends
on V). Primary (`alt==0`, 8324-8346): `W = min_g (b - Xg)' V^+ (b - Xg)`, `X = [1,
relyr]`, in the whitened eigenspace, `q = rank(V) - rank(M)`. Altmethod (`alt==1`,
8347-8386): `d = b_Q - A b_P` with `A` design-only (8360-8364), `S = V_QQ + A V_PP A' - A
V_PQ - V_QP A'`, `W = d' S^+ d`. Because `A` depends only on `relyr`, `S` under the
bootstrap covariance *equals* the bootstrap covariance of the replicate forecast errors
`d*` exactly - the two constructions coincide, no separate accumulation. Reference:
`advx` already forks `Ftail(q, dfr, W/q)` vs `chi2tail(q, W)` at 8388-8389; feed the
resolved bootstrap df computed as bootjnt does at 5599-5606, **never** `e(df_r)`.
`e(adversarial_line_*)` needs no recomputation and must not change (V-free OLS,
7354-7368).

Optional variant: a Monte-Carlo replicate reference. With the weight matrix fixed at
`V_boot^+` the min-distance residual is a fixed linear map `P`, so `W = ||P b||^2` and
`W*_r = ||P (b*_r - mean(b*))||^2` is computable per replication from the accumulator;
`p = (1 + #{W*_r >= W_obs})/(R_eff + 1)`. The residual projection annihilates any line, so
recentering imposes the composite null for free. `O(R k^2)`.

**Null and validity.** Null unchanged by the engine swap: primary = "all (non-reference)
horizon effects lie on *some* line" (sthlp 496-497; ado 8287-8295); altmethod = "the post
coefficients equal the forecast of the pre-period line" (sthlp 521-525; ado 8297-8305).
Only the covariance changes. **First-order valid, no refinement** - a genuine refinement
needs per-rep V (studentized), and reps run under `fastnose` with no `r(V)`. On the
documented size problems: the few-clusters-per-cell channel (sthlp 509-513, "roughly
threefold when cells were spanned by only two clusters, and by an order of magnitude when
each cell lay entirely within a single cluster") *is* the channel a cluster bootstrap
addresses in mechanism, and the help itself prescribes the remedy (sthlp 517-518) - but
this is mechanism plus the help's own recommendation, with **no size probe in this repo**.
The cohort-by-period common shock (sthlp 515-517, "no analytic variance estimator can
capture it") is **not** fixed: a shock common to every cluster in a cell has zero
cross-cluster dispersion, so cluster resampling nested inside the cell reproduces the same
realization in every replicate. `permute()` remains the only remedy there. Sharpest new
hazard: `advx` sets its own df from `rank(V)` by eigenvalue tolerance (8332-8341,
8379-8381); a listwise bootstrap covariance has rank `<= R_listwise - 1`, so thin reps
set the numerator df from the replication budget - the residual-rank p that A11 forbids
(5592-5594). `reps()` has no floor (5430-5435). `bootstrap(bayes)` + a cell-weighted
estimand understates the ES covariance, hence **oversized** `att_it`/`ame_it` adversarial
p's; reachable only under `safety(off)`.

**Engineering.** Computable from what the loop already produces - zero new per-rep
harvesting; every ES column is already a named accumulator column (`bootveskey`
8588-8594), and bootwire already reposts the bootstrap ES covariance at 6491.
**Correction to the recorded plan**: `bootstrap_plan.md:311` says "re-run
wooldid_advstat with the boot V_es post-rewiring". Not implementable - `wooldid_advstat`
contains, at 7367 (about two-thirds through; the program runs to 7396), a
`qui est restore ___wooldidfullmodel`, and that estimate is dropped at
4273 unless `verbose`; post-rewire the call either errors or (under `verbose`) restores
the pre-rewire analytic model over the rewired `e()` the deferred display is about to
read. Call the Mata kernel `wooldid_advx` (8311-8393) directly: it touches only
`st_matrix`/`st_numscalar`/`st_local` and plants the altmethod reference rows itself.
New eclass subprogram `wooldid_bootadv` between bootjnt (250) and `part(1)` (253) - the
display is `e()`-driven (7061-7067, 7200-7206), so ordering alone makes the printed p's
bootstrap-based with zero display change. Gate on the p's the analytic path actually
posted (the 5708 precedent), which inherits `finaladvtest` cancellation (2140-2146,
2590-2594), the B16 gapskip (7326-7332), `amescrub` (2595), the degenerate-covariance
skip (7386-7392) and the per-sg skips for free. Build the covariance **from the
accumulator, not from `e(V_es_<fk>)`**: `bootves` declines silently when a column is
excluded or fewer than 2 listwise rows exist (8560-8575), leaving the analytic block in
place. Guards: every label maps to a column; `R_listwise >= k + 2` (the 8638 precedent);
a **rank-aware** condition floor applied *after* screening structurally-zero rows -
`bootjw`'s `min(sc) <= 0` (8654) would wrongly kill families carrying legitimate
zero-variance planted reference rows (8283-8284, 8293-8295). Any failure posts MISSING
plus a named caution. `e()`: replace `e(adversarial_pval_*)` in place, preserve
`e(adversarial_pval_analytic_*)` hidden.

**Cost.** Zero per-replication. One-off per family: one listwise selection (`R x k`) plus
one `symeigensystem` on `k x k` - `O(R k^2 + k^3)`, k single digits to low teens.
Negligible against one replication.

**Recommendation.** Cover, with caveats: (1) hard rank/condition guard failing to MISSING
plus a caution; (2) the help's size paragraph amended to say the bootstrap addresses the
few-clusters-per-cell channel but **not** the cohort-by-period common shock; (3) a
`bayes` + cell-weighted caution naming the `att_it`/`ame_it` p's with the direction
(oversized).

**Effort.** One focused build session. ~100-150 lines: `wooldid_bootadv` ~80-110, one
shell line, optional ~30-40 line Mata helper only for the Monte-Carlo variant. Docs:
sthlp 398, 509-519, 1008-1011, plus the 256 string. Tests: pooled att+att_it; contreat
ame; subgroups; a thin-reps run that must post MISSING; and `e(adversarial_line_*)`
byte-identical to analytic.

### 3.2 `gteffects_cell_se` — `e(gteffects)` cols 6-9

**Statistic.** No post-rewire route exists. The cells are raw stage-1 interaction
coefficients `<g>._i#<t>._t#1.targetflag` read from the active stage-1 `e(b)`/`e(V)`
(1644-1673; parse 7793-7801), **not** columns of the final `e(b)` (3930-3960); the
aggregates are weighted means of the cells (sthlp 849-856), so a ~20x20 aggregate
covariance cannot be inverted back into 72+ cell variances (`bootstrap_plan.md:57-58`
already records this). The only candidate is per-rep harvesting into a second accumulator
keyed `(i, relyr)` - the two stable columns of `e(gteffects)` (col 1 = the user's own `i`,
628-629/752-753; col 4 = `relyr = tt - min_t - 1`, recomputed per run at 7829 so it
survives a rep's `min_t` shift). Then per-cell deviation-from-rep-mean SE (mirroring
`wooldid_bootfin` 8489-8498) and cols 6-9 rewritten with the 7813-7827 formula.

**Null and validity.** The decisive constraint is **containment**. A cell is keyed
(`_i`,`_t`) and `_i` *is* the user's `i`, so every cell is a subset of one `i`'s rows.
Whenever `cluster()` nests `i` - the norm, and **mandatory** on the cell-weighted leg
(5548-5559: "requires the cluster() variable to nest i") - every cell lies wholly inside
one cluster, and no cluster-resampling design perturbs within-cell data: `standard`
duplicates the cluster's rows verbatim, `bayes` scales them by one Dirichlet draw
(6155-6158, 6272-6278). The replication distribution then carries only comparison-side /
FE / reference-structure variation and omits the **own-cell sampling component** - the
dominant term. This is not a different valid variance estimator; it is a systematically
and silently downward-biased one. Under `bayes` every cell is present in every rep
(Dirichlet draws a.s. positive), the variance is near-degenerate but strictly positive,
so it clears `bootfin`'s exact-zero sentinel `V[k,k] == 0` (8517-8520) and the
`cnt < 0.5*ndone` disclosure (8497) - output is `|z|` astronomical, `p ~ 0` on essentially
every row. Under `standard` a cluster is undrawn with probability `(1-1/C)^C ~ 0.37`, so a
cell is available in ~63% of reps (conditional on a zero-weight cluster's rows being dropped so
the cell's coefficient disappears from `e(b)` - unverified in Stata; the exclusion verdict rests
on the containment argument, not on this figure) - comfortably **above** the 0.5 disclosure
threshold, so
a selected, comparison-only SE ships silently. The bootstrap does not fix the documented
few-cluster size problems for cells; it degenerates them. The one sound configuration is
`cluster()` strictly **finer** than `i` (i = cohort, cluster = person_id) - legal today on
the att/ame-only path, a minority usage, detected nowhere, and "finer" can still mean 2
clusters per cell. The canonical probe is the pathological case: `pb/X/chk_c3.do:38` runs
`cluster(cohort)` with `i = cohort`, so all 72 cells (`chk_c3.log:1469`) sit inside single
clusters and `e(gteffects)` is byte-identical to analytic (`chk_c3.log:1536`, mreldif 0).

**Engineering.** Reachability is favourable: `gteffects` is not in `wdb_ro`, the block is
not repmode-gated (the only worker gate is F8 escomp at 3721), and it survives `fastnose`
(sthlp 861-863). Key alignment holds on the reweighting legs (the rep call at 6279 uses
`repcall`, built at 6164 from the original `i`) but is **destroyed on the re-id leg**:
`qui gegen double __wdbnewi = group(__wdbns \`i')` (6243)
(6236-6252) gives draw-specific codes, and a crosswalk buys nothing statistically because
duplicated clusters are identical rows. A repmode-lean branch is mandatory, not an
optimization: 1654-1666 does 1 gen + 4 `gegen` totals + a tag + an `mkmat` over the full
sample purely for columns a rep never reads. Needs two new Mata functions
(`wooldid_bootgtgrab`/`wooldid_bootgtfin`), an `obsgt()`/`accgt()` thread through
`bootrun`, a bootwire section, a `level()` column rename (7851 hardcodes
`lb_95ci`/`ub_95ci` and ignores `level()`), and a **near-degeneracy guard** the existing
sentinel will not provide. Display needs no change (part(2) is deferred).

**Cost.** Per-rep, not one-off. With the lean gate: a `st_matrixcolstripe` plus one
`regexm` per coefficient name - milliseconds. Without it: 4 `gegen` totals + a tag + an
`mkmat` per rep, plausibly a 10-30% per-rep surcharge. Memory `R x ncell` in Mata: ~80MB
at reps(1000) with 10k cells.

**Recommendation.** **Exclude**, permanently, with a cheap honesty mitigation: post
`ereturn local gteffects_inference "analytic"` in bootwire, extend the caution at 256 with
the *reason* ("each cohort-by-period cell lies inside a single cluster, so cluster
resampling generates no within-cell variation to bootstrap"), and add the same sentence at
sthlp 398 and 859. If overruled, the acceptable shape is restricted (reweighting leg AND
`cluster()` not nesting `i` AND min per-cell distinct-cluster count `>= 10`) plus a
degeneracy guard - never an unguarded rewire.

**Effort.** Recommended path: ~5-8 lines across 3 files, under an hour. If overruled:
~90-120 net lines plus two new screens and a probe matrix covering bayes vs standard,
reid vs reweighting, nesting vs finer, `timetrends`, a vanishing cell, `level() != 95`,
`bootdf(clusters)` and a fine unit-level `i`. 2-4 focused days.

### 3.3 `J_B_plots` — `es_plots_ci` and `cfx_lattice`

**Split verdict; the two items are not the same problem.**

**`es_plots_ci` — statistic.** The band is a per-horizon **marginal** interval, so the
right quantity is what bootwire already computes and prints: `b +/- z_level * se_boot`,
i.e. cols 5/6 of the rewired `e(wdidmainresults)` (6396-6412) / `e(wdidsgresults_<sg>)`
(6428-6442). Do **not** default to a percentile band from `e(boot_ci_percentile)` - the
plotted band would then disagree with the printed table row for the same estimate.
`bootdf()` is already wired (`zc` at 6396 is `invt` under `bootdf(clusters)`/integer,
`invnormal` otherwise), so width tracks `bootdf()` with zero item work. The joint-test
oversizing screens (sthlp 487, 881) do not apply to a marginal band. The cohort-period
common shock is **partially** fixed for `ES_att_it`/`ES_ame_it` under
`bootstrap(standard)`, whose re-id leg (230, 6230-6272) contributes a twice-drawn
cluster's cells twice so cell-level dispersion enters the horizon average; **not** fixed
under `bayes`, where the omitted cell-sampling caution propagates verbatim to those bands
(the observation-weighted `ES_att`/`ES_ame` rows are unaffected).

**Validity.** Not a test; the band inherits whatever validity the printed table row has,
by construction. Today the plot and the table disagree about the covariance for the *same*
estimate - that is the defect. Two behavior changes need blessing: (1) the AME-SE scrub
(`amescrublo=2`, 1495; loops 2256-2258, 2309-2312) nulls the ame ES rows' se/t/p/CI while
keeping the estimate, and bootwire keys exclusion off a **missing estimate** (6363-6366),
so those rows would gain a bootstrap band they have never had - statistically correct,
visible; (2) where `amevoid` fires (1605-1610) the estimate itself is nulled and the band
correctly stays absent. Rows bootwire cannot serve get missing cols 2-6 (6407-6411), so
`rarea` breaks while `connected estimate` still draws - graceful degradation.

**Engineering.** Post-rewire regeneration, **no per-rep harvesting**. `wooldid_esplots`
(7562-7775) reads only the replay gate `e(wdid_replayv)` (7570), `e(wdidmainresults)`
(7574/7594), `e(wdidsgresults_<sg>)` (7695) and display macros (7583-7592) - all posted by
the worker at 4184 and 4196-4226, i.e. **before** the shell block. It consumes only cols
1/5/6/7, and 5/6 are exactly what bootwire overwrites. Three edit sites: gate the worker
call at 4254-4258 on `bootstrap` **only** (mirroring 4237/4241 - `permute` must keep the
worker-side call); add the call in the shell after `part(boot)` (255) and before 256; and
supply the `espre'+`espl' > 0` gate the shell lacks (a hidden `e(wdid_hases)` posted
beside 4216-4220, or a scan of col 7). **One real breakage**: bootwire renames the CI
columns at `level() != 95` (`lb_95ci -> lb_90ci`, 6415-6422 main, 6446-6449 sg; bootjnt
does the same to `pwsgs` at 5791-5795) while esplots hardcodes `lb_95ci`/`ub_95ci` at
7611, 7622-7623, 7642-7648, 7683, 7705, 7716-7717, 7736-7741 and 7764 - so
`bootstrap(standard) level(90) makeplots` would die r(111) on the new path. Fix: resolve
cols 5/6 **positionally** after each `svmat` (7602, 7696) before the `keep`. Saved plot
data comes free (7688, 7768). No new `e()` names.

**`cfx_lattice` — statistic and hazards.** Would need a per-lattice-cell SE from the
replication distribution; **cannot** come from the rewired covariance, because `e(V_boot)`
and `e(V_es_*)` are indexed by `e(b)` column names and the lattice cells live in separate
`e(cfx*)` matrices (4191-4194). Strictly per-rep harvesting, with two statistical hazards.
(a) **First-order alignment**: with `lattice(N type)` the lattice variable is a per-run
quantile tempvar built by `gquantiles` under the run's weights (947-987) and `wdb_ro`
passes the raw spec (211), so every replication **re-cuts the quantiles** -
`lattice_3` in rep 17 is a different dose interval than in the observed run, and the
accumulated dispersion mixes cell-boundary noise with effect noise. Only
`cfxlattice(fixedvar)` or `lattice(existingvar)` (991) gives an honest fixed partition.
(b) **Cell support**: cfx cells are conditional-on-dose subpopulations whose documented
failure mode is few clusters per cell (sthlp 664); a cluster bootstrap does not
manufacture clusters, so coverage converts today's honest blank into a
computable-but-weak band. `bayes`: the cfx block runs for `item=="att"` unless
`safety(off)` (2353) and margins is weighted by `w_att`, i.e. observation-weighted, so
the omitted-cell-sampling caution does not bite in the default configuration.

**cfx engineering.** `wdb_ro` must gain `makecfxplots`, `cfxlattice()`,
`cfxplottypes()`, `makecfxplotsbysg`, `cfxplotnose`. A **repmode gate is mandatory**: the
block calls `marginsplot` and `graph save` unconditionally, so 500 reps would draw 500x
the graphs and overwrite the user's `.gph` files. Accumulators: one `R x (#lattice rows)`
matrix per cfx key aligned on the `lattice_#` stripe (6888-6904) - up to 40 matrices with
`cfxplottypes(post pre -1 0 1)` and 3 subgroups. **Redraw**: `marginsplot` reads the live
`r()` left by the preceding `margins` call and cannot be replayed against a rewired
matrix, so the plots must be hand-rolled (~80-120 lines modeled on 7677-7685) with a
visible style change.

*Judge conflict, advocate-resolved.* The J_B_plots sketch names the repmode-gate sites as
"2417-2418, 2456-2457, **2516-2517** and 2532-2533". The advocate verified by grep that
`marginsplot` occurs at 2417, 2456, **2496** and 2532 and `graph save` at 2418, 2457,
**2497** and 2533, with nothing at 2516-2517; following the sketch would leave the
by-subgroup att-type `marginsplot`/`graph save` un-gated, so every replication redraws
that graph and overwrites the user's `saveplots` `.gph` - exactly the failure the gate
exists to prevent. Use 2496-2497.

**Cost.** `es_plots_ci`: zero extra per-rep work; net runtime delta ~0. `cfx_lattice`:
the dominant per-rep cost of any run requesting it - two `margins`-over-lattice calls per
`cfxplottype` per rep, times `(1 + |sglist|)` under `makecfxplotsbysg`; the worker's own
comment calls a single such call "a pretty taxing request" (2407). With
`cfxplottypes(post pre -1 0 1)` that is 10 margins calls per rep before subgroups, 40 with
three - it can dominate estimation and make reps(500) impractical.

**Recommendation.** `es_plots_ci`: **cover, unrestricted** - the cheapest genuine win in
the inventory. `cfx_lattice`: **cover restricted form at best** (fixed lattice only, with
a per-cell availability disclosure in the shape of `e(boot_navail)`), and **exclude** is
defensible; land `es_plots_ci` on its own and put cfx behind an author ruling.

**Effort.** `es_plots_ci`: ~20-25 lines across three sites plus the caution and two sthlp
edits; half a day, with the regression risk concentrated in the `level() != 95` path,
which needs one explicit probe. `cfx_lattice` restricted: `wdb_ro` additions, four repmode
gates, a bootctl fixed-lattice gate, an accumulator per key through `bootrun`, a new Mata
harvest, a post-loop rewire and a new redraw subprogram - 3-5 days plus timing work and
its own probe matrix.

### 3.4 `ptp` — `pretrendpower()` (`roth_prep`, `roth_slope`, `roth_power`)

**Statistic.** Not a percentile interval and not per-rep harvesting: all three are a
**post-rewire recomputation** of the existing analytic formula with the bootstrap
covariance substituted. `wooldid_ptp` already consumes exactly the two objects bootwire
replaces - `e(b_es_<type>)` at 6649 and `e(V_es_<type>)` at 6650 - and after the `hbar`
block ends at 6639 it touches **no data at all** (verified: no `qui`/`gen`/`gegen`/`sum`/
`e(sample)` over 6640-6876). With `S` = the bootstrap pre-period subblock, `m` = the
trend vector (6660-6672), `b_p` = observed pre coefficients: `roth_prep` is `W = b_p'
S^-1 b_p` (6721-6722) against `Ftail(keff, dfden, W/keff)` or `chi2tail(keff, W)`
(6724-6725); `roth_slope` is `m'S^-1 b_p / m'S^-1 m` (6718) with `se = 1/sqrt(m'S^-1 m)`
(6719) - **the only item whose point estimate moves**, because the slope is GLS with
weights `S^-1`; `roth_power`'s `mdt_wald = sqrt(lambda/c2)` with `c2 = m'S^-1 m`
(6763, with `c2 = m'S^-1 m` at 6711-6712) and `mdt_nis`/`roth_nissize` from `mvnormalcvqp` (`ptpsd` 6738-6741, the call at 6743, the bisection 6766-6806) are
pure functionals of `S`, while `hbar` and the bias multiplier are covariance-free.

**Null and validity.** `roth_prep` tests the same null bootjnt already tests under
bootstrap. The SI-A1 screen's failure is a **rank** pathology - a CRVE has rank `<= G`, so
a k-horizon quadratic form with k near G inverts a near-singular block (1256) - and the
bootstrap covariance is not rank-limited by G, so that degeneracy is removed and the help
already exempts bootstrap+jointtests on exactly this ground (sthlp 398). Not fixed: this
is a bootstrap-**standard-error** Wald, not a bootstrap-t; with few clusters the
covariance is still downward-biased and the test still over-rejects. `bootdf(clusters)` is
a partial small-sample patch, not a refinement. `roth_power` is not a test; under `S_boot`
it becomes the power of the *bootstrap* pre-trend test, which is the internally consistent
object once `roth_prep` is bootstrap-based - today's numbers power a test wooldid does not
report under bootstrap. `bootstrap(bayes)` + cell-weighted propagates in the **dangerous
direction**: `S` understated → `c2` overstated → `mdt_wald` **understated** → the
diagnostic claims the design can detect smaller pre-trends than it can, i.e. it
manufactures false reassurance. Exception: `roth_nissize` is exactly invariant to a
uniform rescaling of `S` (the bounds are `+/- crit*sd_j`, so only the correlation matrix
enters 6738-6741) - do not claim the bayes understatement biases it.

**Engineering.** Computable entirely from what bootwire posts. The **call site** is the
whole problem, and A13 says so: "wooldid_ptp computes/displays/posts inside the worker
from analytic V_es before any rewiring; 'automatic bootstrap Sigma_pre' was unfulfillable
in this ordering" (`bootstrap_plan.md:354`, amending A6 at :310 which had promised exactly
this). Line 4240 calls `cap noi wooldid_ptp` with no `bootstrap` gate (only
`if "\`pretrendpower'" != ""`), while 4237/4241 do gate the
displays, so merely deleting the 5582 refusal is a live mis-output path - ptp would run on
the analytic covariance and print *before* the rewired tables. Blocker for a naive move:
ptp's data block (6618-6639) needs `ttre` and `w`, both worker **tempvars** (423/437;
334-336). Fix: split on the covariance boundary - the worker computes the two
covariance-free `hbar` values and `c_local`s them out; the shell-sited ptp then has zero
data dependence. New site: between `part(1)` (252) and `part(2)` (253), keeping `cap noi`;
`wooldid_ptp` is eclass (6593) and the shell is eclass (19), so its ereturns stick.
**Provenance guard (mandatory)**: bootwire 6483-6492 posts the bootstrap subblock only
when `wdb_vesok=="1"` (6489) and there is **no else** - on any of `bootves`'s three early
returns (8552-8562, 8563-8568, 8569-8574) `e(V_es_<fk>)` silently stays ANALYTIC with no
caution and no hidden twin, so a shell-sited ptp would print unlabelled analytic Roth
diagnostics inside a bootstrap run. bootwire must publish the rewired-family list. **df
plumbing**: under the default `bootdf(normal)` `e(df_r)` is still the analytic residual df
(6511-6516), so the live run would print `chi2(keff)` while `part(roth)` on replay reads
`e(df_r)` at 7097 and prints `F(keff, analytic df)` beside the same stored p - resolve
from `e(boot_stat)`/`e(boot_df)` in **one** place. Not touched: `wooldid_ptpcheck`
(6561-6588), whose three requirements are all bootstrap-independent.

**Cost.** Per-rep zero. One-off: one extra `wooldid_ptp` pass. The cheap part is
closed-form matrix algebra; the expensive part (`mdt_nis` by bisection on `mvnormalcvqp`,
6763-6806) is unchanged and unrelated to the covariance source, already capped at
`k_pre <= 4` (6734) with a documented budget note (6731-6733). Invisible next to `reps()` refits.

**Recommendation.** Cover all three - but only as a **three-precondition package**:
(i) `roth_prep` routed through `wooldid_bootjw` so it agrees with `e(joint_horiz)` by
construction; (ii) the slope doctrine settled explicitly; (iii) bootwire publishing which
families it actually rewired. Adopt `bootjw`'s degeneracy policy (refuse and post MISSING,
8637/8646-8663) rather than ptp's g-inverse-and-report (6704-6712), and append
`pretrendpower` to the 5533 bayes caution list with the direction ("overstates power").

**Effort.** Medium-low: ~6 edit sites, no new Mata, no per-rep code; 1-2 days including
certification. Worker gate + `hbar` move ~30 lines; bootctl ~5; shell ~10; `wooldid_ptp`
signature and guards ~40; bootwire family list ~6; replay df + provenance token ~10; help
~15. Certification is the real cost: reps > k_pre at both `bootdf` settings with
live-vs-replay byte identity; reps <= k_pre must refuse; a `wdb_vesok==0` family must
caution and skip; `roth_prep` vs `e(joint_horiz)` agreement; the bayes caution; the
contreat leg.

### 3.5 `J_B_silent` — `ves_withheld`, `allsgssame_ame_i`, `sescan_caution`

**Unifying diagnosis.** The bootstrap surface is gated on **analytic availability**
throughout - bootwire `cap confirm matrix e(V_es_`fk')` / `if _rc != 0 continue`
(6484-6485) and bootjnt `if `ok' & "`e(joint_`pre'`it'_allsgssame)'" != ""` (5708). Every
analytic-specific withholding therefore propagates into a bootstrap run even when the
accumulator already holds everything needed.

**`ves_withheld` — statistic.** The listwise complete-case replication covariance of the
`ES_<fam>_<label>` accumulator columns - exactly what `wooldid_bootves` already computes
(8570-8578) - posted as `e(V_es_ame)`. **Not** a percentile object: `e(V_es_*)` is
consumed by honestdid/`lincom`/`test`, which need a covariance. Everything needed is
present: `chk_d.log 413-475` posts `e(b_es_ame)` 1x4 and `e(V_es_att)` 4x4 but **no**
`e(V_es_ame)`, with `e(boot_navail)=30` on all 12 columns, and `chk_d.log 480-503` shows
the `ES_ame` rows already reporting bootstrap se/z/p (.152246, .2155856, .2482963,
.3047624) because `amescrub` blanks only cols 2-6 so col 1 survives and bootwire's `excl`
stays 0. Same family, opposite direction and worse: when `bootves` fails, 6489 has no
`else`, so the **analytic** block stays posted with no caution, no hidden twin and no
marker - one dropped or degenerate horizon triggers it, and honestdid then silently
ingests an analytic covariance. Metadata: `e(es_ndegen_<fk>)` (4027) is computed on the
analytic block and never recomputed, while the help (~944-946) says its absence "marks a
clean covariance"; `e(numpre_<fk>)` is positional and stays correct.

**Validity.** The analytic failure (1326-1332, caution 7443) is that a treated i-t cell
lies wholly inside one cluster, so the sandwich has no independent replication of that
cell's dose slope. Under `bootstrap(standard)` with `cluster()` the whole cell enters
0,1,2,... times per resample, so the between-cluster component *is* captured - the
bootstrap is precisely the remedy 7443 names. Three honest limits: (i) `bayes` never
removes a cluster, so for cell-weighted estimands the resulting block is systematically
too small and would be consumed by honestdid, which never sees the 5533 caution; (ii) the
common cohort-period shock (sthlp ~515-517) is **not** fixed by a cluster bootstrap finer
than the cell; (iii) with G=12 and reps=30 the 4x4 is a 30-draw, 12-cluster object with no
df channel of its own, and under the default `bootdf(normal)` a hand-run `test`/`lincom`
would use the **analytic** df.

**`allsgssame_ame_i` — statistic and engineering.** The same bootstrap chain Wald bootjnt
already runs for its six covered items (`wooldid_bootjw` "chain", `df = #sg - 1`). The fix
is **one word**: bootjnt 5695 `foreach it in att att_it att_i att_itime ame ame_it` gains
`ame_i ame_itime`. The row names exist (the worker's "cont" round builds
`rownamespool` at 3441-3458 with the `Catt -> ame` rename at 3450), the 5708 gate makes
the change self-limiting, `suppresssgprimaryeffects` is refused with bootstrap so the
positional walk is safe, and the display formats these identically at 7223-7224.
**Second-order gap the one-liner does not close**: under `amescrub` the analytic ame chain
tests are never computed at all (3318 empties `sgpost`; 3338/3357 and the pre twin
3361-3404 then leave the scalar empty), so the 5708 mirror-gate suppresses the bootstrap
version in exactly the containment run where the bootstrap is the recommended remedy.

**`sescan_caution` — statistic.** Not a p-value; there is nothing to bootstrap. The
bootstrap-side object is the **same screen re-run on the post-rewire surface**: (a) the
identical-SE scan, directly consumable after 6424 (`wooldid_amescreen mode(sescan)` reads
col1=b, col2=se, 7516-7517, exactly `e(wdidmainresults)`'s layout) but with different text
(the bootstrap pathology is a degenerate *replication* distribution); (b) the placeholder
disclosure, whose bootstrap analogue half-exists as `bootfin`'s `-1e100` sentinels
(8516-8522) surfaced at 6547-6548, missing only the "e(V) contains placeholder entries"
clause; (c) the containment caution 7443, which should keep firing (it is the *reason* to
bootstrap) but must stop recommending the flavor already in force. The validity claim is a
display-truth claim: `chk_d.log 325-334` prints "standard errors were not identified for 6
of 12 reported estimates. e(V) contains placeholder entries" and `chk_d.log 480-503` then
shows all 12 rows with bootstrap SEs over a reposted bootstrap `e(V)`; `chk_c3.log
920-925` prints "7 reported standard errors appear to be identical under a cluster-type
VCE ... consider using a different approach to SE estimation" immediately above a rewired
table - i.e. advice to do what the user already did.

**Engineering (sescan).** Three sites: gate 3926-3927 on `"`bootstrap'" == ""` and re-run
the retitled screen from the shell after 250; gate the 3969/3975 display and stash the
count for a reconciled post-rewire line; and pass an "inference already in force" flag
into `mode(detect)` at 1342 so 7443 branches its last sentence. Site (c) must be reworded
**in place** rather than deferred - it fires at 1342, long before the loop, and the shell
`exit`s on a failed loop (232-236), so a deferred version would vanish exactly when the
user most needs the diagnosis.

**Cost.** All three one-off, zero per-replication. `ves_withheld`: one extra `bootves`
call per withheld family, `O(R m^2)`. `allsgssame_ame_i`: at most 4 more `bootjw` calls,
`O(R G^2 + G^3)`. `sescan`: one `amescreen` pass, `O(n^2)` with n < 60. The real cost of
the third is battery expected-output churn from moved and reworded cautions.

**Recommendation.** `allsgssame_ame_i`: **cover outright** (one word; ships regardless of
the other rulings). `ves_withheld`: **cover restricted form in two parts** - Part A builds
the block when `e(V_es_<fk>)` is absent but `e(b_es_<fk>)` is present, synthesizing the
stripe carrier from the `b_es` colnames (no Mata change needed: `bootves` uses `oldname`
only for its stripes, 8547/8581); Part B adds the missing `else` at 6489 and is
recommended **unconditionally**, being a strictly-worse instance triggered by a single
excluded horizon. `sescan_caution`: **cover restricted form** - gate, re-run, reword; do
not invent a p-value.

**Effort.** Small; roughly half a day plus battery reconciliation. `allsgssame_ame_i` ~15
minutes plus the work of *constructing* a probe dataset that reaches subgroups + contreat
+ `att_i`/`att_itime` + >1 realized subgroup + `amescrub==0` (D3 could not reach it). Part
B ~30 minutes / ~8 lines. Part A ~2-3 hours / ~30-40 lines confined to bootwire, plus the
`es_ndegen` recomputation and two help edits; +1 hour for the `f2ok==0` carrier form.
`sescan` ~2 hours of code, dominated by a full battery diff adjudication pass.

### 3.6 `J_B_jointcrash` — `joint_pre_att`, `joint_pre_ame`, `joint_pre_sg`, `allsgssame_att`, `allsgssame_ame`

**All five are already implemented.** The statistic is a Wald quadratic form on the
observed coefficients using the listwise complete-case bootstrap covariance of the
selected `e(b)` columns (`wooldid_bootjw`, 8613-8677): columns selected by a positional
walk of the results-matrix row names (5613-5645), `ok` = listwise complete rows (8622),
`Cm = I(k)` with `df = k` for the pre-trend tests and the successive-difference chain with
`df = k-1` for the equality tests (8624-8637), `S = Cm cov(A_complete) Cm'` (8642-8644), a
correlation-normalised eigenvalue guard with a 1e-8 relative floor (8653-8663),
`W = d'S^-1 d` (8669), `p = chi2tail(df, W)` by default (8677) or
`Ftail(df, dfden, W/df)` under `bootdf()` (8676). Under bootstrap the **analytic** joint
pre-trend tests are never computed at all, so the bootstrap value is the only one there
is.

*Judge conflict, advocate-resolved.* The family reports the exemption sites as "2208, 2656
and 2964 all gate on `"`suppressjointtests'`bootstrap'"==""`". The advocate found a
**fourth** site at 3275 (the per-subgroup continuous-treatment joint tests), carrying the
same gate plus `& `amescrub' == 0`. This broadens rather than weakens the claim, but
anyone editing the exemption from the three-site list would miss a site.

**Null and validity.** Standard Wald-with-bootstrap-VCE; validity asymptotic in G, under
exactly the conditions the bootstrap SEs already need. There is **no asymptotic
refinement** - the replication distribution of the Wald statistic is never formed. On the
documented problems: SI-A1 (1256) is genuinely repaired in mechanism, because the
bootstrap estimates cross-horizon covariation directly from the resampling distribution
rather than from within-cluster score products. **Few clusters is not repaired**, and this
is the substantive caveat: on the suite's own dense dataset with G=12, R=30 the bootstrap
and analytic joint p's differ by six to eight orders of magnitude **in both directions** -
pooled `e(joint_horiz)` 3.06e-20 vs analytic 2.95e-14, but subgroup `e(joint_horiz_s1)`
5.57e-06 vs analytic 1.35e-14 (`chk_c3.log ~1555`, `~1563`). That is instability, not a
uniform bias, and no rep count removes it. `bayes` + cell-weighted propagates to the
`_it`/`_i`/`_itime` columns only; the 5533 caution already says so and warns "p-values can
be extreme" - precisely the condition that pushes p into the crashing denormal band.

**Engineering: the defect.** Nothing needs per-rep harvesting; all five are post-rewire
recomputations. The Mata→Stata handoff `st_local("wdbjp", strofreal(p, "%21x"))`
(8676-8677) emits a hex literal with a **zero leading digit** for any p in
`(4.9e-324, 2.2250738585072014e-308)` - e.g. `+0.012688b70e62bX-3ff` - which Stata's
expression parser rejects with r(198) "invalid syntax". Verified by direct parse test
(`chk_c2.log 52-89`): `%21x` of exact zero is `+0.0000000000000X-3ff` and parses (rc 0,
value 0); a normal value `+1.b758d848fac55X-207` parses; the denormal returns 198 for both
`local x = ` and `scalar s = `. The **first failing evaluation** is the `if `wdbjp' != .`
*condition* at 5662/5679/5711, not the assignment. Because 250 is bare and the shell body
sits inside `prog def wooldid`, the whole command aborts. Reachability is the **default
path**: `Ftail` has a polynomial tail and effectively cannot reach the band, but
`chi2tail` (the default `bootdf(normal)`) reaches it at `W ~ 1417` for df=2. Blast radius
from `chk_c.log 1099-1205`: posted are `e(V)`, `e(V_boot)`, `e(V_es_*)` x8, the rewired
tables, `e(b_boot)`, `e(bias_boot)`, `e(bootdist)`, `e(boot_navail)`,
`e(boot_v_posted)=1`, `e(joint_horiz)=0`, `e(joint_horiz_s1)="0"`,
`e(joint_horiz_s2)="1.4175550282e-156"`; **missing** are `e(joint_horiz_s3)`,
`e(joint_horiz_i)` and every `_i`/`_ame` variant, the CI flavors, `e(boot_z0)`, the
allsgssame overwrite and the 256 caution; **worst**, `e(pwsgs)` is still 12x6 carrying
**analytic** se/t/p because 5724-5802 never ran - exactly the mixed surface the blanking
invariant at 5770-5772 exists to prevent. Only three `%21x` sites exist in the file: 8676,
8677, 8694.

**A suspected second defect, verified and cleared.** bootwire's caution at 6546 says
ame-family structural zeros "were excluded from those columns' bootstrap moments", raising
the question whether `bootjw` (reading the raw accumulator at 8620) sees contaminated ame
columns. It does not: `wooldid_bootgrab` writes an actual missing into the **accumulator
itself** at accumulation time (8455-8457), so `bootjw`'s listwise `ok` drops those rows.
The real consequence is availability - heavy blanking thins the listwise set until 8638
fires and the ame test posts missing with a caution while its att twin survives.

**Two real secondary gaps.** (1) `excl` is not threaded into bootjnt (syntax 5597 has no
`excl()`, unlike bootwire 6336): columns whose reported estimate is missing re-enter
`e(b)` as **structural zeros** (6353-6355) and 8455 explicitly does not blank them (its
`ob[1,k] != 0` test fails), so a horizon bootwire blanked out of the printed table can
still enter `d` and `S` as a zero restriction. (2) The allsgssame gate at 5708 requires the
**analytic** scalar, which the analytic chain test does not post when its rank guard fails
and no contrast fallback succeeds (3064-3077, 3140-3154) - exactly the rank-deficient case
where a well-conditioned bootstrap covariance would still deliver a test. The `ok` flag
already built at 5699-5707 is the natural replacement gate; the counterweight is the A14
mirror invariant.

**Cost.** Zero extra per-rep. One-off, already paid: the K-column walk plus at most 16
joint-test inversions and up to 12 chain inversions on matrices of order <= 4. The only
avoidable cost is the analytic chain test computed under bootstrap purely to serve as the
5708 gate. Fix cost: two lines of Mata, no new `e()` names, no display change, no battery
movement except on runs that currently abort.

**Recommendation.** Cover - already covered; the work is to make it survive. Ship, in
order: (1) **must**, the denormal clamp at 8676-8677 plus the twin at 8694; (2) **should**,
per-post `cap` at 5662/5679/5711 (wrapping the *condition*, not just the ereturn) routing
into the existing `wdbncaut` channel; (3) **should**, a jointtests-specific few-cluster
caveat, since the SI-A1 error routes users here precisely when G is small and the two
documented warnings currently point in opposite directions with nothing telling the user
so; (4) **consider**, the `excl` thread, the 5708 re-gate, and the reference distribution.

**Effort.** Crash fix ~1-2 hours including a suite pin asserting rc 0, `e(joint_horiz)`
finite, `e(joint_horiz_s3)` non-empty, `e(boot_ci_normal)` present and `e(pwsgs)` rewired;
battery risk essentially nil. `cap` hardening +0.5 hour. Few-cluster caution ~1 hour.
`excl` passthrough ~2-3 hours plus a constructed test. allsgssame re-gate ~1-2 hours plus
a full battery diff under the pristine-twice protocol. Reference-distribution change ~2-4
hours of code but it **moves every reported bootstrap joint p** and needs its own
adjudication and an author ruling first - do not bundle it with the crash fix.

---

## 4. Exclusion cases and rulings for the author

Nineteen entries. "Verdict" is the advocate's, except where it says AUTHOR RULING.

**1. `bootstrap()` + `gteffects` (cell se/p/CI, `e(gteffects)` cols 6-9).**
*For exclusion*: structural, not a tolerance question - a cell is keyed (`_i`,`_t`) and
`_i` **is** the user's `i` (628-629, 1655-1666, 7797-7801), and bootctl 5548-5559
*requires* `cluster()` to nest `i` on the cell-weighted leg, so every cell lies inside one
cluster and cluster resampling never perturbs within-cell data. `bayes` then leaves every
cell present in every rep with a near-degenerate but strictly positive variance that
clears the exact-zero sentinel (8517-8520) and the `0.5*ndone` disclosure (8497) → `p ~ 0`
across the whole table. The re-id leg destroys the `(i, relyr)` key outright, and under
`timetrends` the coefficient-level variances are typically not identified, so the analytic se
column is missing (sthlp 858-859) and a bootstrap se
would manufacture inference on an unidentified cell. `bootstrap_plan.md:57-58` already
classified cell SEs out of "everything".
*Against*: it is the only default-printed p-value in a bootstrap run with no bootstrap
twin, printed for tables up to 200 rows (7243-7247) inside the same output block as fully
rewired tables. The engineering is unusually cheap (one token into `wdb_ro`, no repmode
gate needed, survives `fastnose`, stable key), and `cluster()` strictly finer than `i` is
a textbook DiD arrangement where cell bootstrapping is genuinely right and is never
forbidden.
*Verdict*: **EXCLUDE, permanently.** The restricted form needs an inverted nesting probe +
a per-cell cluster count + a re-id carve-out + a near-degeneracy guard, and cannot fire on
any cell-weighted run by construction. Ship the mitigation:
`ereturn local gteffects_inference "analytic"` plus a reason-bearing clause in caution 256
and sthlp 398/859. Independently, this exposes two pre-existing inconsistencies worth
fixing: 7851 hardcodes `lb_95ci`/`ub_95ci` and ignores `level()`, and 7818-7819 uses the
**stage-1** `df_r` rather than `e(df_r)`.

**2. `bootstrap()` + `makecfxplots` / `cfxlattice()` / `lattice(N type)`.**
*For exclusion*: under the default auto lattice the quantile cut points are rebuilt per
replication (947-987; raw spec passed at 211), so `lattice_3` is a different dose interval
in every rep and the accumulated dispersion mixes bin-boundary noise with effect noise -
not the sampling distribution of the reported quantity. The cells are exactly where sthlp
490-493's few-clusters-per-cell problem lives, so coverage converts an honest blank into a
computable-but-weak band. `marginsplot` cannot be replayed against a rewired matrix, so
coverage forces a hand-rolled redraw and a visible style change, and the per-rep cost is
the **only** per-rep blow-up in the whole inventory.
*Against*: a real user-visible CI surface; under a fixed lattice the accumulation is a
textbook bootstrap of a subpopulation mean with alignment machinery that already exists
(8443-8463, 6888-6904); the ame cfx bands currently suppressed by the CRVE refusal would
become drawable from a covariance that never touches the refused CRVE; and it is the item
whose removal would let caution 256 disappear entirely rather than merely shrink.
*Verdict*: **EXCLUDE the auto-lattice case outright** (a scope call, not a deferral).
Fixed-lattice coverage is defensible but is an **author ruling**, not a recommendation.
Correction for whoever builds it: the repmode gates belong on 2417-2418, 2456-2457,
**2496-2497** and 2532-2533 (see §3.3).

**3. `bootstrap()` + `makeplots` / `saveplots` (ES band, saved `_ES.dta`).**
*For exclusion*: the band is a display object; the numbers are already in the rewired
`e(wdidmainresults)`. Moving the call moves graph creation out of the worker's error scope
(65-93), so a graph failure lands after an expensive loop. It silently changes what
archived `saveplots` `_ES.dta` files contain under identical file and column names - a
reproducibility trap. And esplots was never written against a rewired matrix: it hardcodes
`lb_95ci`/`ub_95ci` at nine sites while bootwire renames at `level() != 95`, so the move
is r(111) on any `level(90)` run.
*Against*: a pure **move** with zero per-rep cost - esplots reads only `e()`, and every
object it consumes is rewired by shell line 250. It fixes the loudest internal
contradiction in the package (table row and plotted band for the same estimate from
different covariances), and the subgroup bands come free.
*Verdict*: **COVER, with one mandatory companion and two things to settle first.** Cheapest
genuine win (~20-25 lines). Mandatory companion: resolve cols 5/6 **positionally** after each
`svmat` (7602, 7696) before the `keep`, or the `level() != 95` path dies. To settle: (a) it
silently changes what archived `saveplots` `_ES.dta` files contain under identical file and
column names - decide whether the file should record which covariance produced its bands; and
(b) bless in advance that ame ES rows blanked by the AMESE se-scrub will gain bands (entry 12).
On the gate, this document's lean is to gate the worker-side skip on `bootstrap` **only** -
`permute` runs its own loop at shell 101-177 and its plots are observed-run plots, which is
correct for that engine - but the permute companion proposes the same code site (4237/4241/
4254-4258) be written so both engines can defer. **It is one edit; reconcile the two documents
before touching it.**

**4. `bootstrap()` + `pretrendpower()` (currently a hard refusal, 5578-5583).**
*For exclusion*: covering creates a **second** bootstrap p on the identical null from the
identical covariance as `e(joint_horiz)`, with a different degeneracy policy (ptp
g-inverses and reports at reduced `keff`; `bootjw` refuses and posts missing) and a
different df - two authoritative-looking pre-trend p's that disagree exactly in the
thin-cluster case users care about. The Roth slope is GLS-weighted by `S^-1`, so
substituting `S_boot` **moves a reported point estimate**, contradicting bootwire's stated
invariant (6329-6330). Under `bayes` + cell-weighted the power columns move in the
reassuring direction. And the external route (feed the already-rewired `e(b_es_*)` /
`e(V_es_*)` to Roth's `pretrends`) is the reference implementation and is what the refusal
already names.
*Against*: A13 records the refusal as an **ordering artifact** - "unfulfillable in this
ordering" - not a statistical objection, and A6 originally specified the covered behavior.
`wooldid_ptp` does zero data work after the `hbar` block, so the move is a clean
covariance-boundary split. The `cap noi wooldid_ptp` at 4240 is **unconditional**, so
deleting the refusal alone is a live mis-output path. And `mdt_wald` is documented as
powering the joint pre-trend test, which under bootstrap is bootjnt's - today's numbers
would power a test wooldid does not report.
*Verdict*: **AUTHOR RULING; lean COVER**, but only as a three-precondition package, else
keep the refusal: (i) `roth_prep` routed through `wooldid_bootjw`; (ii) the slope doctrine
settled - recommended is the **fixed-weight variant** (analytic weights
`c = S_a^-1 m / (m'S_a^-1 m)`; slope byte-identical, `se = sqrt(c' S_boot c)`), which
preserves 6329-6330 and is still a genuine bootstrap SE of a well-defined statistic;
(iii) bootwire publishes which families it actually rewired. Without (iii), a shell-sited
ptp prints unlabelled **analytic** Roth diagnostics whenever `bootves` declines - strictly
worse than today's refusal.

**5. `bootstrap()` + `adversarial` / `altmethod_adversarial`.**
*For exclusion*: `pretrendpower` is refused on exactly this reasoning and is the nearer
analogue (both are GLS-on-`V_es` functionals, not linear Walds on `b`), so treating the
two identically is cheaper and more defensible. The primary statistic spans **all**
horizons, so k exceeds any joint test's and the thin-R rank hazard bites hardest
(`reps()` has no floor). The cohort-by-period common shock is not fixed by any cluster
bootstrap nested inside the cell, so the recomputed p would be equally oversized while
losing the explicit analytic-provenance caution. And there is **no in-repo size
measurement** of a bootstrapped adversarial test: sthlp 509-513 simulates the *analytic*
test, and the only bootstrap simulation the help cites is for att standard errors at 24
clusters.
*Against*: the help's own size paragraph closes with "prefer randomization inference (see
permute()) or bootstrap-based standard errors for pre-trend inference" (sthlp 517-518) -
the option does not apply the remedy to the test it is recommended for. The display
co-locates a bootstrap joint p (7046-7057) and an analytic adversarial p (7063-7064) four
lines apart in identical phrasing. `jointtests` + bootstrap already received exactly this
exemption. Every input already exists name-aligned with zero per-rep cost, and A14 means
the display needs no change at all.
*Verdict*: **COVER**, with the same guard package as ptp. The engineering correction is
load-bearing and verified: `wooldid_advstat` carries at 7367 a `qui est restore
___wooldidfullmodel` (7367) against an estimate dropped at 4273 unless `verbose`, so
`bootstrap_plan.md:311`'s recorded sketch is **not implementable** - call the Mata kernel
`wooldid_advx` directly. Do **not** adopt the Monte-Carlo replicate reference as the
default: that is a reference-distribution decision for the whole surface, not one family
(entry 17 / §5). `e(adversarial_line_*)` must stay untouched - verified V-free.

**6. `bootstrap()` + `histogramcohorteffects`.**
*For exclusion*: nothing to cover - all four histogram `margins` calls run `nose` (1866, 1871, 1896, 1903)
and `e(histogramestimates)` (4175) carries point estimates only, so no covariance of any
flavor enters it.
*Against*: none on inference grounds.
*Verdict*: **NO-OP for coverage**, but the documentation is **factually wrong** and must
be fixed: caution 256 and sthlp 398 both say histogram output is produced from the
analytic covariance, and it is produced from no covariance at all. The same sentence needs
a second precision fix: the adversarial clause is true of the p-*values* but false of
`e(adversarial_line_*)`. Since five families propose narrowing this one sentence, it needs
a single owner.

**7. `bootstrap()` + `jointtests` (few clusters, a G floor, the SI-A1 exemption).**
*For exclusion*: a joint test's content is the **off-diagonal** block, and it is G - not
reps - that bounds the information. At G=12, R=30 the bootstrap and analytic joint p's
differ by six to eight orders of magnitude in both directions (3.06e-20 vs 2.95e-14;
5.57e-06 vs 1.35e-14). The SI-A1 screen routes users to the bootstrap precisely when they
cluster coarsely, i.e. when G is small, and the exemption removes the only screen slowing
them down; under bootstrap the analytic joint tests are never computed (2208, 2656, 2964
**and 3275**), so there is no second number to disagree with. `bootjw`'s floor admits a
2-df Wald inverted from four complete replications (8638).
*Against*: the package promises this three times (1257, 1261, sthlp 487-491), and
`permute()` delivers no joint test at all - its rep option string strips `jointtests` and
`subgroups` - so excluding leaves the SI-A1 message naming two remedies and delivering
neither. The mechanism is genuinely right for the specific defect SI-A1 names, and
`bootjw`'s guards already follow A11.
*Verdict*: **KEEP COVERED; do not restrict.** Ship a jointtests-specific few-cluster
caution - this is the one place where two documented reliability warnings point in
opposite directions with nothing telling the user so. Requiring `bootdf()` to be
non-default for `jointtests` is a defensible alternative and is an **author ruling**.
**Must-fix regardless**: the denormal crash in this family aborts the user's whole command
on the **default** path, leaving a half-rewired `e()` whose `e(pwsgs)` still carries
analytic se/t/p under a bootstrap header.

**8. `bootstrap()` + `subgroups()` (per-sg joint tests, chains, `pwsgs`, sg ES bands).**
*For exclusion*: every per-subgroup object is the thinnest listwise object in the run - a
subgroup can lose a horizon in a replicate, so `bootjw`'s `df+2` floor bites there first -
and the observed instability is worst there. The sg joint tests are posted as ereturn
**locals** (5679) while their pooled twins are **scalars** (5662), so a missing sg test is
invisible to a scripted consumer checking for a scalar. Cost scales as `(1+|sglist|)`.
*Against*: `subgroups` is already in the rep option string (214) and is the one place
where the bootstrap surface is materially more complete than permute's, which strips
subgroups entirely. `suppresssgprimaryeffects` is refused with bootstrap, so the positional
walk cannot desynchronize. And `pwsgs` already carries the strictest blanking invariant in
the file (5770-5772).
*Verdict*: **KEEP COVERED - no exclusion.** Two defects to fix: (a) a per-sg joint test
that fails posts nothing at all, so absence is indistinguishable from "never requested"
except via the caution channel; (b) when bootjnt aborts mid-loop, `e(pwsgs)` is left
carrying analytic se/t/p under a bootstrap header (`chk_c.log 1189`, 12x6 un-blanked) -
which is why a bare `cap noi` on the bootjnt call at 250 would make things **worse**.

**9. `bootstrap()` + `contreat()` / the ame family.**
*For exclusion*: the ame columns are the only ones carrying a structural-zero blanking
rule (8455-8457), so their listwise R is systematically thinner than their att twins' and
an ame object can post MISSING while its att twin posts a p - an availability asymmetry
the user cannot predict from the option line. The ame families are also where the AMESE
containment screen fires, i.e. exactly where the analytic covariance was refused for cause
and the help documents the withholding as a contract.
*Against*: the run **already** reports bootstrap ame inference - `chk_d.log 480-503` shows
all four `ES_ame` rows carrying bootstrap se/z/p (.152, .216, .248, .305) while
`e(V_es_ame)` is absent - so the withholding is now internally inconsistent, not
conservative. The cluster bootstrap is precisely the remedy the containment caution
recommends in that same run (`chk_d.log 322-327`).
*Verdict*: **KEEP COVERED** and close two verified gaps. (a) bootjnt's chain loop at 5695
omits `ame_i ame_itime` while the analytic path posts both (4165-4169) and the display
formats them identically (7223-7224) - four scalars keep analytic values beside six
rewired neighbours. One-word fix, and now measured rather than code-confirmed
(`VER_defect/d6.log` A3) - though only on the `bootstrap(bayes)` + `safety(off)` +
`cluster(person_id)` route: `bootstrap(standard)` is structurally barred from reaching it by the
nesting requirement at 5546-5551, and when the cluster nests `i` the analytic scalars are absent
so the 5708 gate never opens (A2, A4). (b) Build `e(V_es_ame*)` from the accumulator when the
analytic block was withheld for containment: verified feasible with **no Mata change**,
because `bootves` uses `oldname` only for its stripes (8547, 8581).

**10. `bootstrap(bayes)` + cell-weighted estimands and everything derived from that
covariance.**
*For exclusion*: the understatement is **directional**, not merely imprecise - joint tests
and (if covered) adversarial p's become anti-conservative, and `pretrendpower`'s
`mdt_wald` would understate, a reassurance diagnostic biased toward reassurance. The 5533
caution already predicts "p-values can be extreme", which is exactly the condition that
pushes `chi2tail` into the crashing denormal band. And a newly-posted `e(V_es_*)` block
escapes the caution entirely, because honestdid/`lincom` consumers never see console text.
*Against*: it is already a `safety(off)`-only path with an explicit, accurate caution;
`bootstrap(standard)` is unaffected and one token away; and refusing derived surfaces
inside an opt-in escape hatch is paternalistic relative to the rest of the package's
`safety()` convention.
*Verdict*: **KEEP the `safety(off)` escape.** Two **author rulings**: (i) extend the 5533
enumeration to name the adversarial p's and, if covered, `pretrendpower` - with the
**direction** stated ("oversized p-values" / "overstated power"), not just "inherits the
understatement"; (ii) decide whether a newly-built `e(V_es_*)` posts at all on that path,
since that surface has no console-caution channel - a machine-readable marker is the
minimum. One exemption to preserve: `roth_nissize` depends only on the **correlation**
matrix, so a uniform rescaling of `S` leaves it unchanged.

**11. `bootstrap()` + `timetrends`.**
*For exclusion*: bootwire's exclusion rule keys on a **missing estimate**, not a missing SE
(6360-6376; the rewire at 6395-6412 then overwrites cols 2-6 for any row whose estimate is
present and whose bootstrap se is positive). Under `timetrends` the `f2ok` screen withholds
`e(b_es_*)`/`e(V_es_*)` entirely with the explicit comment "skip a fully non-identified VCV
(e.g. timetrends)" (3998-4010) and the adversarial p is refused with a caution - yet the ES
**rows** of `e(wdidmainresults)` would still receive bootstrap se/z/p/CI, creating
inference on coefficients the analytic path declared unidentified, under a bootstrap
label, with no caution. `timetrends` is not refused with bootstrap.
*Against*: a bootstrap SE of a coefficient that is **stably normalized** across
replications is a legitimate sampling quantity even where the analytic sandwich is
rank-deficient; refusing would remove numbers from a combination that is not currently
refused; and the identical estimate-keyed fill is the **correct** behavior in the
`amescrub` case, so a blanket rule against it would break that.
*Verdict*: **AUTHOR RULING, probe-needed** - no judge family covered this. The decisive
unknown is whether the `timetrends` ES coefficients carry an identical normalization in
every replication; if the per-rep normalization drifts, the accumulated dispersion is
normalization noise rather than sampling noise, and name-aligned accumulation cannot tell
the difference. No Stata slot was held for this assessment, so the mechanism is stated and
not the outcome. Cheapest safe interim: when the analytic run withheld `e(V_es_<fam>)` for
`f2ok==0`, caution that those ES rows' bootstrap inference has no analytic twin.

**12. `bootstrap()` + the AMESE se-scrub rows (`amescrublo==2`).**
*For exclusion*: a visible behavior change that reads as a regression to anyone who knows
why those cells are blank - the scrub is a deliberate refusal (1495; 2256-2258, 2309-2312)
and bootwire's estimate-keyed exclusion silently undoes it, including in the plotted band
once esplots is moved.
*Against*: verified in the probe that this is already happening and is **correct**:
`chk_d.log 322-327` prints the containment caution recommending bootstrap SEs, and
`chk_d.log 480-503` shows the ame rows carrying bootstrap SEs in that same run. The
bootstrap never touches the refused CRVE, so the objection to the analytic number does not
transfer.
*Verdict*: **BLESS**, with **positive labelling** rather than silence. The verified defect
in the same run is the stale advice, not the fill: "standard errors were not identified
for 6 of 12 reported estimates. e(V) contains placeholder entries" (`chk_d.log 328-333`)
prints **above** a table where all 12 rows carry finite bootstrap SEs over a reposted
bootstrap `e(V)`, and `chk_c3.log 920-925` advises "consider using a different approach to
SE estimation" above a rewired table. Reword 7443 **in place** (it fires at 1342, before
the loop, and the shell exits on a failed loop at 232-236, so deferring loses it exactly
when it is most needed); defer only the 3969/3975 count.

**13. `bootstrap()` + a pooled reference (no `esfixedbaseperiod`).**
*For exclusion*: nothing to decide - under the pooled reference no single omitted period
exists, so the adversarial p is not computed under **either** engine (advstat's inference
block is gated `if "`esf'" != ""` at 7374 (comment 7372); sthlp 501-502) and `pretrendpower` requires
`esfixedbaseperiod` at parse (`wooldid_ptpcheck`, 6561-6588).
*Against*: none. Manufacturing a pooled-reference adversarial p from a bootstrap
covariance would be a new statistic, not an engine swap, and the geometry the composite
null needs is absent.
*Verdict*: **NO NEW RULE** - mirror the analytic gating exactly. This matters as a **scope
statement**: any bootstrap adversarial or ptp coverage inherits the `esf` gate for free,
so a claim that "bootstrap covers every p" is true only over the p-values the
pooled-reference regime actually produces. Say that in the help rather than implying
broader coverage.

**14. `bootstrap()` with no `cluster()` (observation-level resampling on panel data).**
*For exclusion*: `bootrun` defines cluster := observation when `cluster()` is absent
(`qui gen double `obsclus' = _n`, 6118-6123). On panel data that is an iid row bootstrap
that ignores within-unit serial dependence entirely, inside a DiD package whose subject
matter is that dependence. bootctl's few-cluster caution and its 2-cluster floor are both
gated on `cluster()` being present, so this path receives no unit-count screen and no
caution, while the header prints "bootstrap, observation-level, reps(30)", which reads as
an upgrade.
*Against*: with neither `robust` nor `cluster()` the **analytic** default is
homoskedastic OLS (1266-1268), so the observation-level bootstrap is a lateral move that
additionally buys heteroskedasticity robustness - not a downgrade in flavor - and the help
documents the behavior explicitly (sthlp 398, "or observations when no cluster() is
given").
*Verdict*: **DO NOT EXCLUDE.** One narrow defect to fix, raised by no judge:
`bootdf(clusters)` with no `cluster()` resolves G to the **observation** count, so the
small-sample correction the user asked for becomes `t(N-1)`, numerically indistinguishable
from normal, silently. Either refuse `bootdf(clusters)` without `cluster()` or report what
it resolved to.

**15. `reps()` floor and the rank/availability guard policy across all consumers.**
*For exclusion (of the status quo)*: **three different floors coexist on the same
accumulator** - `bootves` posts `e(V_es_*)` from as few as 2 listwise rows (8572),
`bootjw` refuses below `df+2` (8638), and `bootfin` has no floor and merely discloses
below `0.5*ndone` (8497). So a `reps(3)` run can hand honestdid/`lincom` an `m x m`
event-study covariance of rank <= 2 while **refusing** the joint test built from the very
same columns. If adversarial and/or ptp are covered they add two more policies: `advx`
derives its numerator df from `rank(V)` by eigenvalue tolerance (8332-8341, 8379-8381),
which is the *definition* of the analytic test but becomes the A11-forbidden residual-rank
p when V is a thin replication covariance.
*Against*: the no-reps-floor ruling is explicit and reasoned (5430-5437: "no floor -
low-rep runs are legitimate ... see the help file on how many replications inference
needs"), and every guard failure converts a number the user sees today into a missing
value plus a caution - a visible completeness regression some users will read as a bug.
*Verdict*: **KEEP the no-reps-floor ruling**, but make the **per-object guards uniform**
under A11: every quadratic-form and covariance consumer requires `R_listwise >= k+2` and
fails to MISSING plus a named caution, never a rank-derived p. Raise `bootves`'s floor from
2 to `k+2` - the sharpest inconsistency found and one no judge flagged. **Ship-together
dependency**: the floor raise must land **with** the 6489 `else`/marker, because on its own
it makes the silent-analytic `e(V_es_*)` state **more** common. Any `bootadv`/`bootptp`
must screen structurally-zero ES rows (8283-8284, 8293-8295) before applying a condition
floor; `bootjw`'s `min(sc) <= 0` refusal would wrongly kill those families.

**16. The mirror-vs-structural gating principle** (bootwire 6484-6485; bootjnt 5708; any
`bootadv` family enumeration).
*For mirroring*: it inherits every upstream gate for free - `finaladvtest` cancellations,
gapskip, `amescrub`, degenerate-covariance skips, per-sg skips - with no duplicated logic,
and it guarantees the bootstrap `e()` surface is a same-shaped twin of the analytic one,
which is exactly what lets the existing display print it unchanged (A14, 5594-5596).
Computing structurally means a bootstrap run can carry a test with no analytic
counterpart, which breaks that invariant by design.
*Against mirroring*: it is wrong precisely where the analytic path withheld something
because **its own covariance was invalid** - the case the bootstrap exists to serve. Three
verified instances: bootwire skips `e(V_es_ame)` because the analytic block was withheld
for containment; the 5708 gate suppresses the ame chain because 3318/3338/3357 never
computed the analytic one under `amescrub`, i.e. in exactly the containment run where the
bootstrap is the remedy; and 5708 also suppresses a chain test whenever the analytic rank
guard failed.
*Verdict*: **AUTHOR RULING - the single ruling that resolves the most items at once** (the
`V_es_ame` build, the ame chain second-order gap, the adversarial family enumeration and
the `es_ndegen` question, together). Recommendation: **mirror for structural absences**
(estimand absent, columns missing, `esf`/gapskip/`finaladvtest` cancellations),
**compute structurally for analytic-specific withholdings** (`amescrub`, cell containment,
analytic rank failure), and attach a per-object marker wherever the structural route
produces something the analytic surface never had.

**17. Standing refusals to keep** (`permute()`, `impvalues`, `fastnose`,
`unconditionalse`, `customvce()`, `oldsyntax`, `suppresssgprimaryeffects`, multiway
`cluster()`).
*For keeping*: each is a genuine incompatibility rather than a scope gap, verified at
5488-5522 - `permute` and `bootstrap` each orchestrate their own replication loop;
`impvalues` changes the estimand; `fastnose` requests points only; `unconditionalse` and
`customvce()` directly contradict a bootstrap-SE request; `suppresssgprimaryeffects`
removes exactly the columns the subgroup replication vectors need; multiway clustering has
no single resampling unit.
*Against*: none material. Worth recording that the `cleanmatrices` refusal was correctly
**removed** (5569-5577): `matrix drop _all` at 4268-4271 does not touch `e()` matrices, so
the deferred display and any post-loop consumer are unaffected - the refusal had been
blocking a combination that works for adopath-installed users.
*Verdict*: **KEEP ALL.** `permtstat` and `permexact` are unreachable with `bootstrap()` by
the mutual refusal at 5488-5491, so on this engine they are N/A and belong entirely to the
permute half - no bootstrap-side ruling is needed or possible.

**18. The reference distribution, ruled once for the whole surface.** Today: `chi2(k)` by default
and `F(k, G-1)` under `bootdf(clusters)` - a cluster-based denominator - for every joint statistic
the bootstrap posts.
*For keeping the status quo*: it is what the current build documents and pins, every existing
test and battery number is calibrated to it, and changing it moves **every** reported bootstrap
joint p at once - the widest-blast-radius change in this document.
*Against*: the derivation the bootstrap covariance actually supports is a Monte-Carlo replicate
reference (compare the observed Wald against the distribution of the same Wald across
replications), which needs no distributional assumption and is the natural object once the
accumulator exists; the `chi2` shortcut is the one the package's own text elsewhere calls
optimistic at small cluster counts.
*Cost*: no per-rep cost - the replicate Walds are computable from the accumulator already held -
but a battery re-pin of every bootstrap joint p, so it must be its own commit.
*Verdict*: **RULED BY THE AUTHOR, 2026-09-02: make the default conditional on clustering.**
`bootdf()` keeps its three explicit values, but when the user does not specify it the default
becomes `normal` if no `cluster()` was given and `clusters` - t with G-1 df, F(k, G-1) for the
joint tests - if one was. Today the default is `normal` unconditionally and `clusters` is opt-in,
so this changes every *clustered* bootstrap run's reported p-values and confidence intervals and
needs its own commit and a battery re-pin; unclustered runs are unaffected. It is a small edit
(resolve the default in `wooldid_bootctl` / the shell where `wdbdf` is set at 5682-5686, using
the resampling-cluster count `wdb_r_nclus` that `bootjnt` and `bootwire` already receive), but
note two riders. With no `cluster()` the resampling unit is the observation, so a G-1 rule there
would be N-1 and indistinguishable from normal anyway - which is why the conditional default is
the right shape. And the ruling **closes** the replay mismatch rather than widening it, which is
the opposite of what an earlier reading of `sthlp` 956 suggested; measured 2026-09-02
(`ct/Q/qG2_dfr.do|log`, `cluster(cohort)`, 12 clusters, `bootstrap(standard) reps(40)`):

| setting | `e(df_r)` after the run | wooldid's reported p | hand-run `test` p |
|---|---|---|---|
| analytic, clustered | 11 (= G-1 already) | 0 to 14 d.p. | same |
| `bootdf(normal)` (today's default) | 11 | 0 to 14 d.p. (z) | 1.781e-11 - **mismatch** |
| `bootdf(clusters)` (the new default) | 11 | 1.781e-11 | 1.781e-11 - **match** |
| `bootdf(7)` | 7 (rewritten) | 2.2097e-08 | 2.2097e-08 - **match** |

So `e(df_r)` is already `G-1` on any clustered run, and wooldid rewrites it when `bootdf()` sets a
df explicitly; the only self-inconsistent setting is `bootdf(normal)`, which reports a z statistic
while leaving a finite `e(df_r)` behind for `test` to pick up. A small companion fix suggests
itself: under `bootdf(normal)` drop `e(df_r)` from the posted results (keeping the analytic value
as a hidden twin, the way `e(V_analytic)` already preserves the analytic covariance), so that a
hand-run `test` falls back to chi-squared and matches the reported z. The replicate-reference
option discussed above is not part of this ruling and remains available as a future `bootdf()`
value.

**19. Lifting the `pretrendpower()` refusal (5578-5583).** Entry 4 argues the case in full; it is
listed separately here because it is a *scope* ruling in its own right and the headline's list of
seven does not name it.
*Verdict*: **AUTHOR RULING; lean COVER, but only as the three-precondition package of entry 4** -
a shell-sited `wooldid_ptp`, a bootstrap-aware `bootves` decline path, and a labelling rule.
Without the third precondition a shell-sited ptp prints unlabelled analytic diagnostics whenever
`bootves` declines, which is worse than today's refusal. Keeping the refusal is a coherent answer.

---

## 5. Recommended scope and build order

### 5.1 MUST — defect fixes, independent of every coverage decision, ship first

1. **Denormal clamp** in `wooldid_bootjw` 8676-8677 plus the `wooldid_bootpw` twin at
   8694: `if (p < 2.2250738585072014e-308) p = 0` before `strofreal`. Highest-priority
   item in the inventory: a default-path abort of the user's whole command. (Principled
   alternative: `st_numscalar` and read `scalar(...)` at the three call sites - exact and
   parse-free, six lines instead of two.)
2. **Per-post `cap`** at 5662/5679/5711 routing into the existing `wdbncaut` channel -
   wrapping the **condition**, not just the `ereturn` - so a formatting surprise costs one
   joint test, not the `pwsgs` rewire and the entire display.
3. **The missing `else` at bootwire 6489**: when `bootves` declines, the analytic
   `e(V_es_*)` currently survives unlabelled with no hidden twin and honestdid ingests it
   silently.
4. **Truth-fix caution 256 and sthlp 398**: histogram output uses no covariance at all, so it
   does not belong in a list of things "produced from the ANALYTIC covariance"; and note that
   `e(adversarial_line_*)` is V-free while the adversarial **p-values** are the analytic part
   (that clause of 256 is correct as written - only the histogram clause is wrong).
Two further items are often grouped with these but are **not** defect fixes in the same sense -
each changes what a completed run reports, so ship them as their own commits with a battery
re-pin, not inside the crash-fix commit:
5. **Recompute `e(es_ndegen_*)`** on the posted block (it currently describes the analytic
   block, not the one the user gets).
6. **Raise `bootves`'s listwise floor** from 2 to `k+2` - but **only together with (3)**, and
   note entry 15's case against: the no-floor policy is an explicit, reasoned author ruling
   (5430-5437), and every guard failure converts a number the user sees today into a missing
   value plus a caution.

### 5.2 CHEAP COVERS — zero per-rep cost, no new statistics, no display work

7. Add `ame_i ame_itime` to bootjnt's `foreach` at 5695 (one word).
8. Move the `wooldid_esplots` call from the worker (4254) into the shell between 250 and
   253, gated on `bootstrap` only, with cols 5/6 resolved **positionally** after each
   `svmat`.
9. Build `e(V_es_ame*)` from the accumulator when the analytic block was withheld for cell
   containment (no Mata change needed).
10. `wooldid_bootadv` calling `wooldid_advx` directly (never `wooldid_advstat`), placed
    between bootjnt and `part(1)`, replacing `e(adversarial_pval_*)` in place with a hidden
    analytic twin.
11. Reword the 7443 containment caution **in place**; defer only the 3969/3975 no-SE count.

### 5.3 MEDIUM, CONDITIONAL

12. Move `wooldid_ptp` to the shell - only as a package with all three preconditions
    (entry 4). Without the third, a shell-sited ptp prints unlabelled analytic diagnostics
    whenever `bootves` declines: **worse than today's refusal**.

### 5.4 EXCLUDE — permanent scope calls, with disclosure rather than coverage

13. `gteffects` cell SEs/p/CIs — ship `e(gteffects_inference)` plus a reason-bearing
    caution.
14. cfx **auto-lattice** bands — the partition is re-cut per replication; fixed-lattice
    coverage is an author call and is the only per-rep cost blow-up in the inventory.
15. `histogramcohorteffects` — no inference exists to cover; fix the text only.

### 5.5 The display and labelling contract

**Labelling versus blanking - this is the author's ruling, not the assessment's.** The scope
below assumes *labelling* is sufficient. The literal reading of the request ("every p-value ...
produced via a permutation test", and its bootstrap twin) has an alternative route: **blank the
analytic p columns** that the bootstrap does not rewire, so that no un-bootstrapped p is printed
at all. It is cheap on the surfaces where it applies and is recommended below for `gteffects`
alone; generalised, it would strip the adversarial p's, the cfx bands and the Roth block. The
case against generalising it is that it removes numbers users have today and leaves columns that
no bootstrap quantity can refill. Rule this first - it governs the seventeen entries in section 4
and the same fork appears in the permute companion.

Under the labelling assumption, the concrete failure today is not mixture per se but **unmarked**
mixture in identical phrasing: `part(1)`
prints the bootstrapped joint test at 7046-7057 and the analytic adversarial p four lines
later at 7063-7064 with nothing distinguishing them but a trailing global caution. An
all-or-nothing rule would cost the joint tests - the surface the package most strongly
recommends the bootstrap for - to buy uniformity. The requirement instead:

- Every remaining analytic quantity carries a **machine-readable flavor marker** AND a
  truthful **named** caution. Today the only machine-readable signal for "this number is
  analytic inside a bootstrap run" is the *absence* of a hidden twin, which is unqueryable
  in the failure case.
- **Caution 256 / sthlp 398 is a shared edit surface with no owner.** Five families each
  propose narrowing the same sentence, and one of its clauses is factually wrong (histogram
  output uses no covariance at all; the adversarial clause is correct as written).
  Assign one owner and rewrite it once with the correct clause list, or the families will
  produce conflicting partial edits.
- The shipped coverage table must be **per-engine**: `permute` and `bootstrap` cover
  materially different surfaces (permute strips `subgroups` and covers the pooled att
  family only; bootstrap keeps `subgroups` and `contreat` and rewires the ame surface),
  and the help's size paragraphs must keep recommending `permute` for the cohort-by-period
  common shock even after bootstrap coverage widens.
- **Multiplicity stance is unchanged and should be said so explicitly.** Making every p
  bootstrap-based multiplies the number of p-values in a run but adjusts none of them; the
  joint tests **are** the multiplicity instrument the package offers. Do not invent a
  family-wise adjustment as part of this work - it would change shipped analytic numbers
  too. One sentence in the help.
- **`e(df_r)` semantics under the default `bootdf(normal)`**: `e(df_r)` is still the
  analytic residual df (6511-6516). Any new consumer that reads it for a reference
  distribution (a shell-sited ptp, `bootadv`, or `part(roth)` at 7097) silently produces an
  F reference with an analytic denominator. Resolve from `e(boot_stat)`/`e(boot_df)` in
  **one** place.
- **Scalar-vs-macro asymmetry**: pooled joint tests post as `e()` **scalars** (5662) while
  their per-subgroup twins post as ereturn **locals** (5679), mirroring the analytic
  convention. A consumer scripting `if e(joint_horiz_s1) < .` cannot detect a per-sg test at
  all. Any new per-sg surface (`bootadv`'s `adversarial_pval_*_s<sg>`, a **scalar** at
  4150) must state which convention it follows.

### 5.6 The accounting rule for elements missing in some reps

Every object uses a **listwise complete-case** rule over a **different column set**, so
`R_effective` differs per object within one run: `e(V)` over all non-excluded columns,
`e(V_es_<fam>)` over that family's horizons, each joint test over its own pre-block, each
`pwsgs` contrast over its own pair - and the ame columns are additionally thinned by the
structural-zero blanking that writes true missings into the accumulator (8455-8457).
Therefore **`e(boot_completed)` is not the R behind any particular statistic.** Only
`e(boot_navail)` (per column) and `e(boot_N_V)` are published; there is no per-joint-test
or per-`V_es` effective R. Publishing one per posted quadratic form is the honest fix and
is a precondition for any "how many reps stood behind this p" question. The uniform guard
policy is entry 15: `R_listwise >= k+2` at every consumer, failing to MISSING plus a named
caution.

### 5.7 `e()` naming proposal

- **Replaced in place** (bootstrap-sourced values under the analytic names, so the display
  prints them unchanged): `e(adversarial_pval_<t>[_s<sg>])`; under a covered ptp,
  `e(roth_prep_*)`, `e(roth_slope*)`, `e(roth_W_*)`, `e(roth_keff_*)`,
  `e(roth_nissize_*)`, `e(roth_power_*)`.
- **Hidden analytic twins** (extending the `e(V_analytic)` / `e(V_es_analytic_*)`
  convention at 6462/6490): `e(adversarial_pval_analytic_<t>[_s<sg>])`, and under a
  covered ptp optionally `e(roth_slope_se_boot_*)` if the fixed-weight fork is chosen.
- **New provenance macros**: `e(boot_ves_ok)` / `e(boot_ves_built)` (families whose
  `V_es` was rewired, and those *created* from replications and therefore having no
  analytic twin), `e(gteffects_inference)` = `"analytic"`, and one per-surface flavor local
  for any surface left analytic.
- **Unchanged**: `e(adversarial_line_*)`, `e(roth_hbar_*)`, `e(roth_npre_*)`,
  `e(roth_ok)`, `e(numpre_*)`.
- **Recomputed on the posted object**: `e(es_ndegen_*)`.
- Define the scheme **once** before any family ships its own.

### 5.8 Build order

(1) §5.1 items 1-2 with the suite pin. (2) §5.1 items 3-6 as one commit (the ship-together
dependency). (3) §5.2 items 7-9 and 11 - independent, cheap, small battery movement.
(4) The reference-distribution and mirror-vs-structural rulings. (5) §5.2 item 10
(`bootadv`). (6) §5.3 item 12 (ptp) if the three preconditions are accepted. (7) The
documentation rewrite under a single owner.

---

## 6. Open questions

### 6.1 Rulings requested (ordered by how many items each resolves)

1. **Mirror-vs-structural gating** (entry 16). Does the bootstrap surface mirror analytic
   availability, or compute structurally where the replication columns exist? Decides the
   `V_es_ame` build, the `amescrub` chain gap, the adversarial family enumeration and
   `es_ndegen` at once.
2. **Reference distribution, ruled once for the whole surface.** Today: `chi2(k)` by
   default and `F(k, G-1)` under `bootdf(clusters)` - a **cluster**-based denominator; the
   derivation-backed analogue for a Wald on a covariance estimated from R replications is
   Hotelling's `W*(m-k+1)/(m*k) ~ F(k, m-k+1)` with `m = R-1` - a **rep**-based
   denominator; and the adversarial family proposes a Monte-Carlo replicate reference.
   Adopting any of these for one family and not the others is an internal inconsistency;
   adopting Hotelling everywhere **moves every reported bootstrap joint p** and needs its
   own battery adjudication. Decide first, build second.
3. **The Roth slope doctrine fork** if ptp is covered: accept the movement (what the
   refusal message itself tells users to do by hand with `pretrends`), or fix the weights
   at the analytic `c` so the slope is byte-identical and only `se`/CI become bootstrap.
4. **`bayes` + cell-weighted**: extend the 5533 enumeration to name the adversarial p's and
   `pretrendpower` **with direction**, and decide whether a newly-built `e(V_es_*)` posts
   at all on that path. `roth_nissize` is exempt (correlation-only).
5. **`jointtests` + few clusters**: a caution, or require `bootdf()` to be non-default?
   The six-to-eight-orders-of-magnitude divergence at G=12 is the evidence, and the SI-A1
   screen routes users here precisely when G is small.
6. **`timetrends` fill** — probe-needed (does the per-rep normalization stay stable?).
   Interim: caution that those ES rows' bootstrap inference has no analytic twin.
7. **Fixed-lattice cfx coverage** yes or no.

### 6.2 Secondary questions

- **Guard-failure fallback for `bootadv`**: post MISSING plus a caution (A11 discipline),
  or fall back to the analytic p with documented provenance (the fallback recorded at
  `bootstrap_plan.md:311`)? The first is consistent with the rest of the surface; the
  second never removes a number a user sees today.
- **Availability asymmetry for `bootadv`**: mirror the analytic gating exactly, or also
  post where the analytic covariance was degenerate but the replication covariance is full
  rank? (Subsumed by ruling 1, but the adversarial case is where it is most visible.)
- **A reps floor for the adversarial test only**, given that its primary statistic spans
  all horizons and is the highest-dimensional quadratic form in the surface - against the
  standing no-floor ruling.
- **Size evidence**: is a simulation of the bootstrapped adversarial line test wanted
  before the help's size paragraph (509-519) is amended? The few-clusters-per-cell claim is
  currently mechanism plus the help's own recommendation, with no measurement in this repo.
- **`roth_prep` vs `e(joint_horiz)`**: route through `wooldid_bootjw` so the two agree by
  construction, or keep the ptp GLS form with its own `keff` convention?
- **Degeneracy policy for ptp**: `bootjw`'s refuse-and-say-so, or ptp's
  report-with-caution? For a bootstrap covariance, rank deficiency usually means "too few
  reps"; recommendation is refuse.
- **Subgroup Roth diagnostics**: `wooldid_ptp` never handles `_s<sg>` families (the enumeration at 6643; documented at sthlp 537)
  although bootwire rewires those blocks. Status-quo parity, or a separate wanted item?
- **The `hbar` split**: two new blanket-listed worker→shell locals, or a hidden `ereturn`
  the shell-sited ptp reads back (which posts a `roth_*` name before the pass completes,
  interacting with the `e(roth_ok)` marker at 6866-6871)?
- **`excl` passthrough into bootjnt**: fix, or accept and document? Columns whose reported
  estimate is missing re-enter `e(b)` as structural zeros and can enter `d`/`S` as a zero
  restriction.
- **The allsgssame gate**: move off the analytic scalar to the already-computed `ok` flag
  (5699-5707)? And is the analytic chain test worth computing under bootstrap at all when
  5711 overwrites it?
- **A caution on a bootstrap joint p of exactly 0** after the clamp: at G=12/R=30 a p at
  the underflow floor is arguably evidence the covariance is too thin rather than evidence
  against the null.
- **Wrapping the bare bootjnt call at 250**: a plain `cap noi` makes things **worse** (it
  leaves analytic `pwsgs` under a bootstrap header); any wrapper needs a cleanup path.
- **ES band flavor**: normal-theory from the rewired lb/ub (agrees with the printed table)
  or percentile from `e(boot_ci_percentile)`? Recommendation is the former; is a
  `plotci(percentile|bc)` opt-in wanted?
- **Saved `_ES.dta` provenance**: should the file record which covariance the band came
  from (a constant column or a dataset note carrying `e(wdid_vcespec)`)? Today the file and
  column names are identical under both engines.
- **Plots drawn after the loop**: acceptable that a user watching a long bootstrap sees no
  plots until the end, and that a graph failure lands outside the worker's `cap noi` scope?
- **The `esplots` move gate**: written so both engines can defer (`bootstrap``permute`), or
  bootstrap-only for now?
- **cfx**: is restricting to a fixed lattice acceptable, or should the observed auto-lattice
  assignment be frozen into a real variable and threaded through `wdb_ro`? Are you willing
  to lose `marginsplot`'s styling? Should cfx harvesting be opt-in rather than automatic
  whenever `makecfxplots` is present?
- **`gteffects` near-degeneracy threshold** if it is ever covered: author-set, or a
  ratio-based caution reporting the bootstrap/analytic SE ratio per cell?
- **Independent of any bootstrap decision**: `gteffects` ignores `level()` and hardcodes
  `lb_95ci`/`ub_95ci` (7851), and uses the stage-1 `df_r` rather than `e(df_r)` (7805,
  7818-7819). Fix separately?
- ~~**`allsgssame_ame` is unexercised**~~ - **closed by measurement.** On
  `continuous_control_dense` with `cluster(person_id)` the analytic ame scalars do post and the
  bootstrap rewires them (`VER_defect/d6.log` A1/A3); the earlier "never post" reading was
  specific to `cluster(cohort)`, where the containment screen scrubs the chain. What remains
  open is only whether the ame subgroup contrasts are identified across more than one subgroup.
- **`bootdf()` guidance**: should the help say `bootdf(clusters)` is the right choice when
  `e(V_es_*)` will be handed to honestdid/`lincom`, given that under the default `e(df_r)`
  keeps the analytic df?
- **Unverified in Stata** (no slot for this assessment): that a shell-sited `wooldid_ptp`'s
  ereturns stick; that `e(sample)` and the `e(b_es_*)`/`e(V_es_*)` surface survive the
  `est sto`/`restore` round trip intact; that `wdb_vesok==0` is reachable on a realistic
  dataset; and, for `gteffects`, whether a zero-weight cluster's rows are dropped so the
  cell's coefficient disappears from `e(b)` (the ~63% availability figure rests on this).

---

## 7. Artifacts

Session scratch under `pb/X` (not committed):

- **`chk_a.do/.log`** — `att att_it cluster(cohort) esprelength(2) espostlength(2)
  subgroups(sg_three) permute(shuffle) permreps(30)`, rc 0. Only `e(perm_p_att)` =
  `e(perm_p_att_it)` = .032258 (log 246-253, RI block 226-229) are permutation p's; the
  main table (448-465), three subgroup tables, `e(pwsgs)` (510-554, degenerate se 2.33e-10
  p 6.5e-87) and four allsgssame p's (306-309, `.`) all print analytic p's with no caution
  or label; `e(perm_p_pretreatment_att)` absent (499); `e(permdist)` 30x2,
  `e(perm_comparator)` 1x2; `e(b)` 1x48 with every SG/ES column named (428-438).
- **`chk_b.do/.log` B1** — adding `esfixedbaseperiod jointtests adversarial` under default
  safety → SI-A1 refusal r(198) (60-73), whose text recommends "wooldid's built-in
  randomization inference (the permute() option)" for joint-type inference.
- **B2** — same with `cluster(person_id)`: analytic joint p's 0.0000 and adversarial
  "single line: 0.0000" print unlabelled (136-144), then 0 of 30 reps estimable → r(459)
  (147-153).
- **B3** — `adversarial esfixedbaseperiod cluster(cohort) permute`, rc 0; the adversarial
  composite-null covariance is degenerate on this data so no p (162-165);
  `e(adversarial_line_att/att_it)` posted analytic (296-299); RI block att/att_it only.
- **B4** — `jointtests adversarial pretrendpower(80) safety(off) permute(shuffle)
  permreps(30)`, rc 0: `e(joint_horiz)` = `e(joint_horiz_i)` = 1.2e-24 printed as 0.0000
  (380-383, 464-465); the Roth block prints analytic slope se .0013738 and joint p 0.0000
  (386-428, `e(roth_prep_att)` 470); nothing labelled analytic; the `safety(off)` caution
  (315-322) recommends `permute()` while permute is on.
- **`chk_c.do/.log`** — the brief's exact bootstrap spec (`... jointtests adversarial
  esfixedbaseperiod gteffects histogramcohorteffects bootstrap(standard) reps(30)
  rseed(1)`): reps complete, then r(198) "invalid syntax" (1085-1089); partial `e()` at
  1096-1205 (mreldif `e(V)` vs analytic .0104, `e(V_analytic)` 1e-20, `e(V_es_att)` .0022;
  `e(joint_horiz)=0`, `e(joint_horiz_s1)="0"`, `e(joint_horiz_s2)="1.4175550282e-156"`;
  `e(joint_horiz_i)`, `joint_horiz_s3`, all att_it sg tests, the `pwsgs` rewire, the CI
  flavors and the 256 caution all absent).
- **`chk_c2.do/.log`** — root cause: `local x = +0.012688b70e62bX-3ff`
  (= `strofreal(1e-310,"%21x")`) → r(198), while zero and normal-range hex literals parse
  (62-89); reproduced in wooldid at 96-108 / 121-134 (the s3 p is denormal); the same call
  **without** `subgroups()` completes (143-221, bootstrap block 204-213).
- **`chk_c3.do/.log`** — the same spec without `esfixedbaseperiod`, with an analytic
  `safety(off)` comparator, rc 0. mreldif: `e(V)` .0207 (`e(V_analytic)` 0); `e(V_es_att)`
  .0043 (analytic hidden 0); `e(V_es_att_s1)` .0184; `e(wdidmainresults)` .90;
  `e(wdidsgresults_1)` .95; `e(pwsgs)` rewired (se .165 p .78 vs analytic 2e-10);
  `e(gteffects)` 0 and `e(histogramestimates)` 0 (unchanged); `e(adversarial_line_att)` 0;
  `e(joint_horiz)` 3.06e-20 vs analytic 2.95e-14; `e(joint_horiz_s1)` 5.57e-06 vs
  1.35e-14; `e(joint_att_allsgssame)` .691 / pre .762 vs analytic `.`; caution 256 printed
  (1361-1364); the identical-SE sescan caution printed **inside** the bootstrap run
  (920-925).
- **`chk_d.do/.log` D1** — `contreat(dose_vary) ... permute(shuffle) permreps(30)`, rc 0:
  `e(perm_p_att)` = .0323, `e(perm_p_ame)` **absent** (305), RI block lists att only
  (226-233), `e(wdidmainresults)` identical to analytic (mreldif 0, 315) - ame rows
  analytic and unlabelled; the amescreen caution recommends randomization inference inside
  the permute run (156-160).
- **D2** — same with `bootstrap(standard) reps(30) rseed(1)`, rc 0: `main:ame` se .199
  p .138, `pretreatment:ame` and `ES_ame` rows all rewired (337-397); `e(V)` rewired
  (mreldif 1); `e(V_es_att)` rewired (.0043); **`e(V_es_ame)` absent** (413-475, with
  `boot_navail` 30 for every `ES_ame` column) because the analytic block was withheld for
  cell containment (caution 325-334); "e(V) contains placeholder entries" printed inside
  the bootstrap run (330-334).
- **D3** — `att att_i contreat(dose_vary) subgroups(sg_three)`, analytic vs bootstrap:
  allsgssame att .691 (analytic 1.5e-81), att_i .699 (analytic `.`), pre_att .762;
  `e(joint_ame*_allsgssame)` did post analytically once `cluster(person_id)` was used
  (`VER_defect/d6.log` A1), and the four `ame_i`/`ame_itime` scalars were then measured
  bit-identical to analytic beside rewired neighbours on the `bootstrap(bayes)` route (A3); the
  5695
  `ame_i`/`ame_itime` omission is therefore measured, not code-confirmed only; `pwsgs` ame_i rows
  rewired (1505-1518).

Design documents amended by this assessment: `review_process/bootstrap_plan.md` (A6 at
:310-311, A13 at :354-356, the "not functions of e(b)" list at :57-58) and
`review_process/bootstrap_report_2026-08-30.md` (:248-252).
