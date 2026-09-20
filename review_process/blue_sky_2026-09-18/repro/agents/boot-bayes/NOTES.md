# boot-bayes probe notes (adversarial, Dirichlet engine)

Target: `bootstrap(bayes)` — Dirichlet cluster reweighting. Method: REIMPLEMENTATION of
the weight generation + distributional testing.

## Code map (pinned build)
- `wooldid_bootrun` ado:8885-9552 — map frame, draw matrix, rep loop.
- Mata `wooldid_bootdw()` ado:14269-14318 — the Dirichlet generator.
  `E[k,r] = rgamma(1,1,a,1)`, a = n_c under size, 1 under equal; `P = E :/ colsum(E)`;
  `if size: P = P :/ N`.
- weight application ado:9444-9447: `repw = basew * drawv` (MULTIPLIES the user weight).
- rep call: `wooldid y i t ttre [aweight=repw], <repopts>` with repopts built ado:321-364.
- bootsave W columns: `wooldid_bootsavemk` ado:14905-14915 (obs row = wobs, ado:437).

## STATUS
- [DONE p01] build asserted (`which wooldid` -> pinned path), baseline bayes run OK.
- [DONE p02] **BITWISE REIMPLEMENTATION SUCCEEDED** for BOTH bootalpha(equal) and
  bootalpha(size): my own independent Mata implementation of the documented process
  (per-rep outer, per-cluster-ascending inner Gamma(alpha,1) from the seeded stream,
  normalize by colsum, divide by n_c under size) reproduces the multipliers wooldid
  actually applied, EXACTLY (`WF == MP` -> 1, maxabsdiff = 0).
  => I own the engine. RNG consumption order is as declared.
  Documented claims confirmed on this run: equal rowsum = 1; equal obs row = 1/G;
  size sum_c n_c*w_c = 1 exactly; size obs row = 1/N (N = 437 = e(N) here).

- [DONE p03] bootalpha() validation surface CLEAN. `equal`/`size` accepted; `EQUAL`,
  `Size`, `siz`, `eq`, `1`, `equal2` ALL refused r(198) with the full message — no
  silent default, no case-folding, no abbreviation. `bootalpha()` with
  `bootstrap(standard)` -> r(198); with no bootstrap() -> stray-companion r(198).
  boottstat+bayes -> r(198). bayes+att_i -> r(198) default safety, caution under
  safety(off). bayes unclustered, bayes+bootdisplay(bca), bayes+poisson all rc=0.
- [DONE p04] Dirichlet multiplier MULTIPLIES the user weight: replaying rep r by hand
  with `manw = bigpop * w_c` reproduces `e(bootdist)[r,1]` BITWISE in both modes;
  the negative control (aw=w_c alone) differs in the 1st significant digit.
- [DONE p07] No Stata-matrix ceiling: st_matrix() from Mata builds 15000 x 10 fine
  (c(max_matsize)=11000 does not bind on Mata-created matrices). Unclustered bayes
  (C := N) has no dimension limit. CLEAN.
- [DONE p05/p05b + gof.py] equal-mode GOF at R=1000, C=12: rowsum(p)=1 to 4.4e-16;
  cluster means vs Dirichlet(1..1) max|z|=1.84; per-cluster KS vs Beta(1,11) min
  p=0.113 (Bonferroni 0.0042); offdiag corr -0.09078 vs theory -0.09091; lag-1/2/3
  serial corr across reps <= 0.094; exchangeability permutation p=0.494; pooled PIT
  mean .50111 var .08316 (theory .5/.08333). CLEAN.
- [DONE p06] **EXTERNAL ORACLE FOUND AND PASSED.** Stata 19 ships `rwgen`/`bayesboot`
  (/Applications/Stata/ado/base/r/rwgen.ado v1.6.5). At a common seed:
    bootalpha(size)  == `rwgen bayes S, reps(R) id(st)`                  max rel diff 1.3e-14
    bootalpha(equal) == `rwgen bayes E, reps(R) id(st) priorpowers(-n_c)` * n_c, rel 9.5e-14
  (the residual is my CSV round-trip, not wooldid: `export delimited` truncates).
  So the code comment "matches bayesboot BITWISE on unweighted data" (ado:8899) and
  sthlp:548's bayesboot-parity claim are CONFIRMED against the real command.
- [DONE p05b size arm] size GOF at R=1000: rowsum(n_c*w)=1 to 4.4e-16; KS vs
  Beta(n_c, N-n_c) min p=0.0239 (Bonferroni 0.0042, passes); mean|z| ok; obs row =
  1/N = 1/108. CLEAN.
- [DONE p09] **LEAD FINDING.** bayes + cell-weighted estimands under safety(off) on a
  design with cluster-level effect heterogeneity (p09_panel):
     main_att     se_standard 0.315955  se_bayes 0.282868  ratio 0.895
     main_att_it  se_standard 0.315955  se_bayes 0.059040  ratio 0.187
     main_att_i   se_standard 0.332717  se_bayes 0.059647  ratio 0.179
     main_att_itime se_std    0.332717  se_bayes 0.059647  ratio 0.179
  AND on this data att and att_it are the SAME NUMBER (analytic table: both
  1.4669328, both se .07839043) - bootstrap(standard) reports them identically,
  bootstrap(bayes) reports SEs 4.8x apart. p15_attvsattit.do is the clean repro.
- [DONE p11] my "bayes silently keeps near-degenerate replications" hypothesis is
  NEGATIVE. 20 clusters, 2 never-treated: standard failed 61/500 reps (caution
  printed), bayes 0/500. bayes draws kurtosis 2.88 (NOT heavy-tailed), sd/IQR-sd
  0.957, corr(|dev|, control-cluster Dirichlet mass) = -0.057, mean|dev| 0.090 for
  control mass<0.02 vs 0.074 otherwise. No tail inflation. What IS true: bayes's
  interval is 14% NARROWER than standard's there (sd .0923 vs .1075) and reports no
  failures, while standard discloses that 12.2% of its replications were discarded.
- [DEAD] p08 died on an invalid macro name (`sd_unw_standard_.`); p08b is the fix.
- [DONE p14] **BITWISE against the official oracle**: bootalpha(size) == `rwgen bayes
  S, reps(20) id(st)` at seed 31415, max|diff| = exactly 0. bootalpha(equal) ==
  `rwgen bayes E, ... priorpowers(-n_c)` * n_c to 1.93e-16 relative (one ulp from the
  multiply). The ado's "matches bayesboot BITWISE" comment is TRUE.
- [DONE p13] doc-vs-code: bootsave size obs row = 1/437 (data rows) not 1/375 (e(N)).
  sthlp:557 says estimation sample. -> FINDING F2 (low).
- [DONE p08b] documented quantitative claims:
    unweighted  sd(bayes equal)/sd(standard)   = 0.9782  (theory sqrt(C/(C+1))=0.9837) OK
    unweighted  sd(bayes size)/sd(bayes equal) = 0.2406  (claim sqrt(C/N)=0.2620)      OK
    weighted*   sd(bayes equal)/sd(standard)   = 0.6362  <- but * = p04_panel has ONE
                cluster at 500x the weight of the rest, a deliberately pathological
                design. p17 repeats it on lognormal population weights.
- [DONE p16] reproducibility bitwise at a common rseed; reps() prefix property holds;
  centring within +/-0.20 boot-sd on 26 columns across event study, subgroups,
  jointtests; F001 weighted-centring fix CONFIRMED BUILT (-0.019 size, +0.086 equal).
  contreat() needs contreatcontrols() - my input error, informative refusal (p18 redo).
- [DONE p17] realistic weights: [none] sd(eq)/sd(std)=0.9258, [pop aw, max/min 8.3]
  0.9331; sd(size)/sd(equal) 0.3003 / 0.2959 vs sqrt(C/N)=0.2887. Both documented
  quantitative claims hold. p08b's 0.6362 was the 500x-weight-spike design only.
- [RUNNING] p12b coverage MC (45 draws x {bayes att+att_i safety(off), standard att,
  analytic}), nominal 95%, cluster()==i.
- [QUEUED] p19_homogeneous (does the cell-weighted understatement need heterogeneous
  effects, as the refusal text says?), p18_contreat_etc.
- [DONE p19] **F1 sharpened and it is decisive.** One run, two arms:
    het=0    standard se(att_it)=0.123593   bayes se(att_it)=0.074588  (bayes/std .604)
    het=1.2  standard se(att_it)=0.343761   bayes se(att_it)=0.074588  (bayes/std .217)
  standard's att_it SE rises 2.8x with cluster-level effect heterogeneity; the bayes
  one is IDENTICAL to the printed 6 decimals. bayes se(att) tracks standard (.921/.956).
  The analytic clustered SE is also blind (0.099348 both) - the cohort x relyr dummies
  absorb a cluster-constant treated-cell shift, so the residuals do not move.

## FILE MAP (for a replacement picking this up)
  p01_smoke.do            build assertion + baseline bayes run
  p02_reimpl.do           THE core result: independent reimplementation, bitwise match
  p03_valid.do            bootalpha()/gate validation surface
  p04_wreplay.do          multiplier multiplies the user weight (bitwise replay)
  p05_gof.do / p05b_analyze.do / gof.py   R=1000 distributional GOF (analysis in Python)
  p06_oracle_and_size.do  rwgen oracle (CSV) + size GOF at R=1000
  p07_bigN.do             Stata matrix ceiling (none)
  p08b_sescale.do         documented SE-scaling claims + F001 centring
  p09_cellwtd.do          cell-weighted SE magnitudes
  p11_fewcontrols.do      few control clusters (my hypothesis: NEGATIVE)
  p12b_coverage.do        coverage MC at cluster()==i
  p13b_docsample_standalone.do  FINDING F2 repro (self-contained)
  p14b_oracle_standalone.do     FINDING F3 repro (self-contained)
  p15_attvsattit.do       FINDING F1 repro, balanced design (self-contained)
  p16_surfaces.do         reproducibility, prefix property, ES/subgroups/jointtests
  p17_seratio_realistic.do  SE ratios on realistic weights
  p18_contreat_etc.do     if-sample, boundaries, BCa
  p19_homogeneous.do      FINDING F1 SHARPEST repro (self-contained)
- [DONE p12b/p20/p21] coverage at cluster()==i, 30 panels, nominal 95%, true att=1:
      analytic att   26/30 wald, se/sd .706
      analytic att_i 25/30 wald, se/sd .819
      bayes    att   23/30 pctl, 28/30 wald, se/sd .809
      bayes    att_i 24/30 pctl, 24/30 wald, se/sd .621   <- worst of the five
      standard att   28/30 pctl, 29/30 wald, se/sd .848
  All five reported SEs are below the true sampling sd. The bayes-vs-standard att
  percentile gap (23 vs 28) is paired and 5/5 same-direction; widths within 2%, so
  the bayes interval is SHIFTED, not narrower. Marked UNCERTAIN (n=30).
  NB my p12b postfile used str8 for the estimand name, which truncated "main_att_i"
  onto "main_att" and pooled them; p20 recovers them from postfile ORDER, and
  p22_coverage_selfcontained.do is the fixed single-file repro.
- [DONE p22] p22_coverage_selfcontained.do reproduces the p12b/p20 coverage table
  EXACTLY from a single self-contained file. That is the repro cited for F4.

## FINAL STATE: 4 findings (F1 high, F2 low, F4 medium; the oracle discovery is a
## method note, not a finding). FINDINGS.md is complete and self-describing.
