# VERIFIER 1 (adversarial attack) report — AMESE R1/R3/R4 implementation

Target: commits 7054878 + ea69544, loaded via wooldidinstall/wooldid.ado.
Scripts/logs/data: this directory (v1a-v1j.do/.log, cdat*.dta). Builders: test/suite/mc/amese_cases.doh (d1).
~30 completed wooldid runs + 4 replays across 9 scripts (budget ~10-14; overrun spent attributing
three live findings — advisor-endorsed). Scale: OK / BREAK / DESIGN-QUESTION. All runs rc 0 unless stated.

## 1. False-fire attacks (v1a.log)

**P1** contreat + cluster(unitid) + espre/espost(4) + esf + anticipation(1) + jointtests, d1:
no fire — ame se .0848, V_es_ame posted, joint_Choriz posted, zero cautions. **OK**

**P2** cluster(region), region=mod(i,3), COARSER than i (cells nest inside region clusters):
fires ("165 treated cell(s)"), ame rows point-only, att rows + att joint test intact, V_es_ame
withheld, b_es_ame posts, plus the "not identified for 10 of 20" plumbing caution. Firing here is
CORRECT per the mechanism (cell scores still aggregate inside one cluster) — semantics check. **OK**

**P3** contreat + robust + [pw=w]: no fire, real ame SEs. Statically provable too: amescrubclus
empty without cluster()/customvce-cluster (ado 1139-1141); sescan cluster-gated (3394). **OK**

**P4** dose zeroed outside treated cells + contreatcontrols(untreatedZero) + safety(off) +
cluster(unitid): no fire, ame se .0628, V_es_ame posted. **OK**

**P2b** two-way cluster(i unitid): fires on the contained i-margin, rc 0, message names both vars.
The two-way CRVE's V_i term is blind while V_u/V_iu are not, so ANY-margin firing is a
conservative refusal, not exact. Defensible; surfacing. **DESIGN-QUESTION (minor)**

## 2. Scrub-leak hunt on firing configs (v1b.log, v1d.log)

**P5** kitchen sink (d1, basic + cluster(i) + lattice(all 3) + makecfxplots + cfxplottypes(post
pre) + summarystats + histogramcohorteffects + makeplots + saveplots + jointtests + adversarial +
pretrendpower(80)); predicate + vzero + refit all fire; rc 0:

- **BREAK 1 (headline): e(cfxame_cfxpre)/e(cfxame_cfxpost) carry full se/t/pvalue/ll/ul columns**
  from margins dydx(contreat) on the invalid cluster V. `cfxseopt` is empty unless cfxplotnose
  (ado 2078); r(table) rows 1..6 are captured (~2150-2157) and ereturned unfiltered (3655-3658);
  grep of all 24 amescrub sites shows none inside the cfx block 2065-2226, and an enumeration of
  all 21 ame-margins call sites confirms the two cfx sites (2157, 2223) are the only unguarded
  SE-producing ones. Observed values are the pathology itself: t = -27.8/-68.4/+73.5,
  p ~ 1e-22..1e-34 with 15 clusters, and the se column is IDENTICAL between cfxpre and cfxpost
  per lattice point (.034379/.04156/.0077391) — the collapsed-variance signature. The same
  numbers are drawn as 95% CI bands (rarea) in the saved cfx .gph files. Violates §3(a) "every
  ame-family SE/t/p/CI" and the new sthlp paragraph's own promise. cfxatt_* (semielasticity wrt
  targetflag) also posts SEs but is att-side — noted, not ruled. cfxlattice() routes through the
  same call, same leak (static; not separately run). **P5b**: subgroups + makecfxplotsbysg — e(cfxame_cfxpost) and
  e(cfxame_cfxpost_sg1) leak identically. Fix shape: force nose or scrub rows 2..6 at the cfx
  captures under amescrub.
- Everything else on the fired e() surface is clean: wdidmainresults/wdidsgresults_1 ame rows
  point-only; b_es_ame(+_s1) post without V_es_ame(+_s1); V_es_att posts; only att adversarial/
  joint/roth keys; no ame pretrendpower keys (self-gate works); sescan correctly silent;
  histogramestimates and summarystats are point/moment-only. **OK**
- v1sink_ES.dta: ES_ame rows carry missing lb/ub; ES_att intact. **OK**
- Bare replay: rc 0, ame cells still missing, att-only joint lines. est sto / clobber with
  regress / est restore / replay: rc 0, still scrubbed. (Replays do not re-print the cautions —
  table only; cosmetic observation.) **OK**

**P5c** (v1d.log) subgroups with treated cohorts in BOTH sgs, firing config: e(pwsgs) ame/pre_ame
rows scrubbed (estimate kept); joint_ame_allsgssame, joint_pre_ame_allsgssame, joint_Choriz_s#
all suppressed; att-side sg scalars survive per scope. **OK.** Out-of-scope observation: the ATT
pairwise rows show se = 8.7e-09, t = 14,242,721, p = 8.3e-188 on noise data — att-side
sg-comparison variance fully collapsed, displayed uncautioned (sescan scans only the main stack
and is off when amescrub fires). Corroborates the filed allsgssame observation (e53671f);
att-side is the author's scope ruling, so no verdict.

**gteffects**: statically unreachable with contreat (hard r(198) at ado 3880-3881; verified).
The mandate's gteffects-leak question is **moot**. Binary gteffects untouched.

## 3. R4 refit attacks (v1c/e/f/g/h/i/j.log)

**Weights/controls/est-state (R4a/c/e)** byCohort + cluster(i): [aw=w] refit att b matches the
same-data cluster(unitid) twin at 3.8e-15 (weight carried); fullmodel after refit is
cmd=regress, vce=cluster, clustvar=i, wtype=aweight (verbose retention + est restore); replay
after refit rc 0; [pw=w] pair likewise (3.8e-15, se recovered); controls(c.x) pair 5.8e-15. **OK**

- **BREAK 2 (headline): on byCohort configs the refit silently changes the reported ame POINT
  estimates by O(1).** P6a refit ame = 1.1431 vs P6b same-data reghdfe ame = 0.16778 (reldif
  .835) with att matching to 1e-15; controls pair reldif .686; clean unweighted attribution pair
  E2: refit ame = 0.99943 vs reghdfe-retained ame = 0.18763 (impvalues twin keeps reghdfe since
  amescrub=0; att reldif 2.6e-16). Mechanism: regress and reghdfe drop different members of the
  collinear {byCohort cohort-specific dose controls} x {treated-cell dose interactions} set;
  margins noestimcheck treats omitted as zero, so dydx(contreat) is not invariant to the drop.
  Contradicts the printed promise "ame point estimates are unaffected", the sthlp, and the
  battery characterization "~1e-12 b drift" — runbattery.do contains NO byCohort config (grep),
  and cert A5 (byCohort+refit) never pins ame b. Scope refinement: basic+lattice refit drift is
  benign 6.2e-12 (I pair), and configs whose vzero never fires are exactly identical (G, H pairs
  at reldif 0) — the O(1) drift is specific to cohort-collinear dose-control families
  (byCohort; byTreated/byTtre untested, same structure). e(b_es_ame*) consumers on refit
  byCohort runs receive these shifted points.

- **BREAK 3: bare fe() factors are silently linearized by the refit.** E1a fe(region2): refit
  fullmodel contains a single linear `region2` coefficient, R2 0.363, df_m 11; E1b fe(i.region2)
  twin keeps factor dummies (0b/1/2/3.region2), R2 0.733, df_m 13 — same data. reghdfe absorb()
  treats a bare varname as a categorical FE; the refit demotes it to a linear covariate with no
  warning, though wooldid's own reg-punt path prints exactly this warning for unconditionalse/
  customvce (ado 3829-3832 "manually specifying any needed i. prefixes"). In this balanced
  orthogonal design att b happened to be invariant (8.9e-15) and att se moved little, but the
  reported ame differs 10% between E1a/E1b and every residual-based quantity sits on a gross
  misfit. ccc_absorb() shares the exposure (same punt line, same missing guard; not separately run).

**R4d count logic + skip branch**: code reading — amevzi/amevzt are post-esamp glevelsof counts
of i and t; correct for what the refit newly materializes (the cell-interaction terms already
exist in the reghdfe model), EXCEPT fe()/ccc_absorb() levels are uncounted (with fe(i.bigfactor)
the refit would materialize an unguarded dummy set; with bare fe() it linearizes instead —
BREAK 3's other face). Empirical skip test J1 (30 i-levels + 440 realized t-levels = 470 > 450;
treated cohorts span only 12 periods so the model stays small): rc 0, containment + vzero
cautions print, NO refit caution, att b present with att se missing (b-only), ame scrubbed,
V_es_att withheld. Skip branch behaves exactly per contract. **OK** (with the fe-count gap noted).
[Practical note: my first skip probe (30x440 fully-balanced) was killed after 80 CPU-min —
treated cells are per-period beyond the ES window, so t-heavy treated panels make the internal
regression enormous; cost observation only.]

**DESIGN-QUESTION (R4 narrowing collateral, exhibit I1b)**: basic + lattice(all 3) +
cluster(unitid) vzeroes WITHOUT containment (amescrub=0), so the narrowed refit refuses to run:
the user gets a full 4-of-4 blackout (att and ame SEs all missing) plus advice to "re-run with
robust or cluster() finer than the cells" — but they already cluster finer than the cells. The
ea69544 narrowing (protecting designed b-only fallbacks) also withdrew the repair from
legitimately-broken non-containment contreat configs.

## 4. permtstat interaction (v1d.log)

**P11** contreat + cluster(i) + permute(shuffle) permreps(20) permtstat permtci: rc 0, 20/20
reps, 0 failed, perm_stat=t, perm_p_att=.762, perm_lb/ub and perm_tlb/tub posted; observed run
att se intact (.0928), ame scrubbed; EXACTLY one containment caution for the run; 59.5s total.
**OK behaviorally.** **DESIGN-QUESTION**: the design note asserts "permute rep runs are fastnose —
the predicate code must be unreachable there"; under permtstat the shell adds fastnose only when
permtstat is absent (line 101), so the predicate + scrub run quietly in every rep (cap qui at
4790/4831 swallows the caution). Behavior is clean — perm t-stats are att-scale and att se is
untouched — but the note's invariant as written is violated.

## 5. sescan false positive (v1d.log)

**P12** binary, balanced, SINGLE treated cohort, att+att_it+att_i+att_itime, cluster(i): the four
estimands legitimately coincide (identical b and se to 12+ digits, main and pretreatment rows)
and the G-b caution fires ("4 reported standard errors are numerically identical ... treat ...
with suspicion") — a false positive telling users to distrust correct SEs on a valid,
ordinary design. Caution-only, no wrong numbers. **DESIGN-QUESTION** (the note fixture-verified
only the unclustered balanced case; the clustered gate does not exclude coinciding estimands).

## 6. Caution-count preregistration (source-echo lines excluded)

containment expected/observed: A 2/2, B 2/2, C 3/3 (pre-kill), D 2/2, E 3/3 (E2b impvalues
correctly not evaluated), F 1/1, G 1/1, H 1/1, I 1/1 (only each pair's cluster(i) leg fires;
I1b's non-containment vzero prints no containment caution), J 1/1. vzero/refit: B 2/2, E 4/3
(E2b vzero prints, refit correctly gated off), I 2 vzero/1 refit (I1b gated off), J 1 vzero/0
refit, G,H 0/0. No double-prints, no swallowed cautions, no sescan misfires outside P12.

## Verdict summary

FAIL overall (implementation breaks its own contract in three places), with the predicate itself
(fire/no-fire semantics), the e()/display/replay scrub outside cfx, R3 suppression, permute
behavior, and the R4 skip/hygiene surfaces all holding exactly as designed.

- BREAK 1: e(cfxame_*)/e(cfxame_*_sg#) + saved cfx .gph CI bands carry invalid-V se/t/p/ci on
  firing configs (makecfxplots/makecfxplotsbysg; lattice or cfxlattice route).
- BREAK 2: R4 refit changes reported ame points by O(1) on byCohort refit configs (1.143 vs
  0.168; 0.999 vs 0.188) — contradicts the printed "points unaffected"; invisible to the battery
  (no byCohort configs) and to cert A5 (ame b unpinned).
- BREAK 3: refit silently linearizes bare fe() (and by code-identity ccc_absorb()) — the punt
  path's own i.-prefix warning is absent here; R2 0.363 vs 0.733 on the same data.
- DESIGN-QUESTIONS: multi-way cluster conservative fire; permtstat reps run the predicate
  contrary to the note's fastnose invariant; sescan FP on coinciding estimands; narrowed refit
  strands non-containment vzero blackouts with misdirected advice (I1b); feasibility count
  omits fe()/ccc_absorb() levels; att-side pwsgs se 8.7e-09/p 8e-188 uncautioned (out of scope,
  corroborates filed allsgssame observation).
