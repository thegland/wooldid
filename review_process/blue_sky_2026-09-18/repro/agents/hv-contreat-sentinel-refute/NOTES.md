# SKEPTIC notes — refuting the hv-contreat-sentinel findings

Build asserted pinned in every log (`which wooldid`).

## My own DGP (deliberately unlike the prober's)
`r1_own.do` builds `r1_panel.dta`: hospital-year administrative panel.
i = STATE (12), 30 hospitals per state, T = 2008..2019 (4320 rows).
Expansion years 2012 (4 states) / 2014 (3) / 2016 (2); 3 never-expanding.
`dose = expo_h * 1{year >= expyr_s}` — hospital-level exposure switched on at the
state's own onset (the canonical exposure x post continuous-DiD dose).
True dose coefficient 0.8.

## Runs
- `r1_own.do` / `.log` — sentinel + V_es dump, 3 VCEs.  (died at an uncaptured
  `dis` after the mata call failed; the sentinel evidence landed first.)
- `r2_adv.do` / `.log` — splice integrity + restore behaviour.
- `r2b_matatest.do` / `.log` — **harness fact**: wooldid's own Mata functions are
  NOT in Mata memory after a run (`mata describe` shows only my own), so a caller
  cannot reuse `wooldid_advW`. Transcribe instead.
- `r3_honest.do` / `.log` — adversarial kernel (my transcription), honestdid,
  pretrends, and the decisive synthetic control.
- `r4_boot.do` / `.log` — bootstrap neighbour + basic/time-invariant controls +
  ame point-estimate accuracy.

## RESULTS SO FAR

### hv-F5-verify — reproduces exactly on my data (IID)
4 sentinels in e(V) (`pretreatment_ame`, `ES_ame_pre4/pre3/pre2`), all b = 0.
e(V_es_ame) diag = 0,0,0,1.737,1.737,1.737. e(es_ndegen_ame)=3.
THREE separate announcements fire, not two:
 1. "Caution: a treatment period coefficient interacted with something ... was
    dropped for multicollinearity" (+ "triggered 41 times total")
 2. "Caution: standard errors were not identified for 4 of 16 reported estimates"
 3. "Caution: e(V_es_ame) contains 3 of 6 horizons with zero or missing variance
    (structurally degenerate rows); inspect the diagonal before joint inference -
    **honestdid**/lincom on those rows is not meaningful."
Documented at sthlp:903 (untreatedZero's expected zero doses deliberately do NOT
trigger the dose-slope screen), sthlp:1325 (V_es degenerate rows carry zero or
missing), sthlp:1616-1617 (e(V) carries the negative placeholder).
=> the two representations disagree BY DOCUMENTED DESIGN. Verdict PARTIAL,
severity low -> none as a defect.

### hv-adv-ame-rank — mechanism reproduces, bitwise
My transcription of wooldid's own `wooldid_advW` reproduces the PRINTED
adversarial p to <1e-12 for both families (MATCHPRINTED=1, att and ame, both VCEs).
- ame: W(full 6x6) == W(post-only 3x3) BITWISE, q = 2 both. IID and cluster(hid).
- att comparator on the same runs: W full != W post, q = 5 vs 2, p 0.917 vs 0.696.
- e(V_es_ame)'s entire 3 pre ROWS (diagonal AND the off-diagonals with the post
  block) are exactly 0, so the eigen-truncation provably cannot keep any
  pre-period direction. Verdict CONFIRMED; severity is disclosure-only.

### hv-honestdid-Minvariant — reproduces, but the CONTROL refutes the attribution
ame delta(rm), IID: Original [-0.3519, 4.8148]; M=0/0.5/1 all
[-0.3430373248, 4.828910034] — identical to every printed digit.
cluster(hid): all three M [-0.2242453055, 4.659319125]. att widens normally
([-0.227,0.502] / [-0.338,0.554] / [-0.510,0.621]).
**CONTROL (r3_honest.do, "SYNTHETIC CONTROL (no wooldid)"):** a hand-built
b = (0,0,0,1.2,1.0,0.8) with an exactly-zero 3x3 pre block and an ordinary post
block, fed straight to honestdid — NO wooldid anywhere — is M-invariant too
(SYN_RM_INVARIANT=1: [0.6186186186, 1.783783784] at M=0, 0.5 and 1). The same b
with a non-degenerate pre block widens (SYN2_RM_INVARIANT=0).
=> The M-invariance is honestdid's behaviour on ANY degenerate pre block. wooldid
ships inputs that are correct on its own documented convention, names honestdid
by name in a caution at estimation, posts e(es_ndegen_ame)=3, and lists
"no structurally degenerate rows" as a REQUIRED handoff condition (sthlp:1321).
Verdict PARTIAL, severity high -> medium.

### hv-noVsentinel-marker
`ereturn list` after the run and after `estimates use`: e(es_ndegen_ame)=3 present,
4 sentinel diagonals present, no scalar naming the 4. Verdict CONFIRMED, low.

### NEW (mine)
- `pretrends power 0.5` on the SAME degenerate ame inputs REFUSES loudly
  (r198: "Power level .5% too high ... PreTrendsBisect(): Invalid starting
  points"), on both VCEs, while it works on att (slope for 50% power = .10215).
  So the two advertised third-party handoffs disagree: pretrends fails loudly,
  honestdid fails silently.
- honestdid `delta(sd)` dies r(3498) on the degenerate block — and dies
  identically on my wooldid-free synthetic, so that is also honestdid's.
- e(V)'s spliced ES_ame post block vs e(V_es_ame)[4..6,4..6]: mreldif
  2.6e-11 (IID) / 2.9e-11 (cluster) — NOT bitwise equal, but below the project's
  1e-10 indifference threshold. Not filed.

## FINAL STATE (all runs complete)
Verdicts and the new findings are in FINDINGS.md. Repro files that a skeptic can
run with nothing else:
 - repro_seeded_claims.do  -> all four seeded claims + the wooldid-free control
 - finding_A_untreatedzero_regresspath.do -> the CRITICAL new finding
Supporting: r1_own.do, r2_adv.do, r2b_matatest.do (harness fact), r3_honest.do,
r4_boot.do (bootstrap/basic/subgroups/ame-accuracy), r5_ameoracle.do (the oracle
that caught finding A), r6_customvce.do (the six-path comparison that localised it).
