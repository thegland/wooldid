# hv-fastnose-doc-refute — running notes (skeptic pass)

Role: refute HV-FN1 / HV-FN2 / HV-FN3. All evidence below is from MY OWN data and
MY OWN do-files; I read the prober's files but ran none of them.

## My design (all probes)
40 states x 16 years (2005-2020) x 20 individuals = 12,800 rows. i = state,
t = year, ttre = adoption year (2010 / 2013 / 2016; 12 states never treated).
`dose` = state policy generosity, CONSTANT WITHIN STATE (max within-state SD
measured 1.9e-14), spread 0.5-3.0. DGP: y = 2 + a_s + .03*(yr-2005) + .40*dose*D + N(0,1).
With i = state this is exactly the help's own named d3 case (sthlp:902, "a dose
constant within each cohort"). Structural dose coefficient is 0.40, not zero.

## Probes
- `p1_fn1_binary.do` / `.log`  — FN1: binary att, analytic vs fastnose
- `p2_fn2_d3.do` / `.log`      — FN2/FN3: three arms, analytic / fastnose / safety(off)
- `p3_permute_and_ic.do`/`.log`— interactivecontrols screen twin; permute syntax miss (r198)
- `p4_permute_d3.do` / `.log`  — permute(shuffle) permreps(20) fast-mode default on d3
- `mkdata.do`                  — standalone generator (panel_states.dta)

## Results so far

### FN1 — CONFIRMED (the prober's refutation of the recon lead holds)
p1: fastnose arm gives `e(b_es_att)` rc 0, `e(numpre_att)` = 3, `e(V_es_att)` rc 111,
`e(V)` rc 111, `e(wdid_vcvstate)` = absent, colnames `pre4 pre3 pre2 post0..post3`.
Analytic twin: b rc 0, V rc 0, numpre 3, vcvstate analytic-blocks. Max abs difference
between the two coefficient vectors 4.441e-16. The caution "no horizon of the att event
study has an identified variance ... The estimates stand" prints in the fastnose arm.
So sthlp:1325's sentence is TRUE of fastnose on the pinned build, and known_bydesign §7's
last bullet (item 11 "looks unfinished") is refuted by measurement.
Paper trail: held_for_ruling_2026-09-13.md's last column is headed "recommendation",
not "ruling" — so (a) was the RECOMMENDATION; 6fe814d's own message says the code was
changed instead. The file header (added 40d46ac) states item 11 rests on that message
plus partial code evidence; my run supplies the missing code evidence.

### FN2 — behaviour reproduces exactly, framing does not
p2 A2 (robust fastnose) on my data: the R6 caution is ABSENT; main:ame = 0,
pretreatment:ame = 0, all seven ES_ame rows = 0 in the printed table and in
e(wdidmainresults); e(b_es_ame) posts with seven exact zeros (%21x +0.0...X-3ff);
e(numpre_ame) = 3; the printed ame message is the "The estimates stand" caution.
A1 (analytic robust): R6 caution prints in full, every ame cell missing,
e(b_es_ame) rc 111. att side identical across arms (.70128565 in both).

BUT:
* **sthlp:902, the `contreat()` option entry itself:** "...under default safety wooldid
  names this situation, reports the ame results as missing rather than as the
  regression's structural zeros ... Only the raw e(b) coefficient vector retains the
  structural zeros ... **and the check is skipped under fastnose.**"  The skip is
  documented in the option paragraph for the very option that creates the design.
* **A3 (robust safety(off)) measured: identical output.** ame rows print 0, e(b_es_ame)
  posts, and NO caution names them. So the prober's stated expectation — the zeros may
  reach output "only behind an explicit opt-out AND a warning" — is false on the warning
  half: safety(off) is equally silent. sthlp:1209 says so in advance: "the ame rows print
  the regression's structural zeros instead of missing values, and without the caution
  that names them. This is the same output an explicit safety(off) produces."
* Both arms print the generic "dropped for multicollinearity ... triggered 420 times".
=> PARTIAL, severity down from high to low.

### FN3 — CONFIRMED at low
p2 A1: e(b_es_ame) rc 111 while e(numpre_ame) = 3 posts. ado:4149 is unguarded.

### Negative result (my own probe, p3 B1/B2)
The interactivecontrols BS-3 screen twin (ado:12005) carries the identical
`"`fastnose'" == ""` gate, and it does fire that way: the caution "every treated
cohort-period slope term for the interactive control urate was dropped ... this run is
identical to controls(i.t#c.urate)" prints in the analytic arm and is ABSENT under
fastnose. NOT a finding: sthlp:1068 documents it verbatim — "(that check does not run
under fastnose or safety(off), where the reduction is silent)". Same pattern as FN2,
same documentation. This is a deliberate, documented, twice-repeated convention.

### p4 (permute) result — FN2 does NOT escalate
permute(shuffle) permreps(20) on the d3 design: ame rows reach the reportable table as
estimate 0, p_perm 1.000, [0,0] — but the engine prints "Caution: the permutation
distribution is degenerate (every completed draw equal) for: main_ame pretreatment_ame
ES_ame_pre4 ... ES_ame_post3; those permutation p-values are uninformative."
e(b_es_ame) rc 111, e(numpre_ame) empty (sthlp:283). Prober's REACH claim confirmed.

### Dating
git log -S on wooldidinstall/wooldid.sthlp: both the "the check is skipped under
fastnose" sentence (sthlp:902) and the interactivecontrols "does not run under fastnose
or safety(off)" sentence (sthlp:1068) enter at 81c284b, 2026-09-10 — a week before the
pinned build. Documented in advance, not patched after.

## FINAL VERDICTS
FN1 CONFIRMED (low) — the prober's refutation of the recon lead holds.
FN2 PARTIAL  (high -> low) — behaviour exact, framing refuted by sthlp:902 + the
             safety(off) arm; residue is an intra-help inconsistency in the fastnose
             option paragraph.
FN3 CONFIRMED (low).
New findings: none. Two clean negative results (IC screen twin; subgroup V_es_*).
