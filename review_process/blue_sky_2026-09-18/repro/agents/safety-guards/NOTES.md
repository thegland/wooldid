# safety() guard-by-guard probe -- NOTES

Build asserted: /private/tmp/.../bluesky/build/wooldid.ado  (p01_tokens.log:56-62)

## Guard map (pinned build, line -> gate variable)
`safetyinput` = the user's own typed value.  `safety` = WORKING value, rewritten to
`off` at ado:4829 when contreatwithin+contreat() are both given.

Gated on USER value (`safetyinput`) -- contreatwithin does NOT relax:
  718   customweightmultiplier() >= 0                      (r411)
  772   F13 onset-constancy within i                       (r459)
  1355  si8 granular-i caution + maxvar preflight          (caution / r907)
  1397  SI-A1 joint-test cluster screen                    (r459 / caution)
  5246  weights >= 0                                       (r411)
  5290  lattice() domain                                   (r411)
  5315  cfxlattice() domain + missing caution              (r411 / caution)
  5340  t-ttre commensurability                            (r411 / caution)
  6107/6116 permute()+subgroups design-population screens  (r459 / caution)
  6439  bootstrap(bayes) + cell-weighted estimand          (error / caution)
  12248/12272 sgdom subgroup integrality + unit-splitting  (r411/r459 / caution)
  11204 sgdeg subgroup-contrast degeneracy                 (blank vs report)

Gated on WORKING value (`safety`) -- contreatwithin DOES silently relax:
  4617  unconditionalse + cell-weighted estimand           (error)
  4623  unconditionalse + aweight                          (caution)
  4832  makecfxplots for non-att estimands                 (documented)
  4959  i must be positive non-zero integers               (r411)   <-- NOT in the help's
  4973  float-stored i/t/ttre above 2^24 caution           (caution) <-- working-value list
  836   esf: periods without a control unit                (caution)
  908   esf: sufficient pre/post periods                   (r459)
  955   no never-treated control group (also & ctw=="")    (r459)
  1713  wooldid_ccscreen: whole post-estimation collinearity family, incl. R6 ame void
  5788  wooldid_fcdeg collapse-degeneracy re-evaluation
  12086/12108 reltime first-period reference-forcing caution
  12148 (NOT gated -- always runs)

## Results so far
p01: safety(ON|OFF|Off|oFF|0|1|yes|no) all r(198) "safety() must be on or off";
     safety(off ) and safety( off) accepted.  The 500-line validation gap is BENIGN:
     the uppercase token is treated as safety-ON by the early checks (conservative)
     and the run then dies at ado:5473 with the right message.  CLEAN.
p02b: i coded 0 (0-indexed cohort codes) -- the ado:4959 refusal is skipped under
     contreatwithin+contreat, BUT the numbers are IDENTICAL to the 1-indexed run
     (att -.95083572 / ame -.92666667 both ways), and identical again under an
     explicit safety(off) on the binary path (att 1.5800311 == the 1-indexed
     default-safety answer).  _i==0 does NOT collide with the never-treated
     sentinel because never-treated rows always carry _t==0.  NO number moves.
     Non-integer i under contreatwithin: r(452) (uninformative crash, not silent).
     FLOAT CAUTION: ado:4973 fires twice without contreatwithin (log:922,928) and
     ZERO times with it.  <-- live doc-vs-code divergence, see FINDINGS.

## Later results
p03: region#t FE holding one treated state does NOT trip r(322) -- reghdfe removes the
     state's rows as singletons first; reported att 1.0714286 == design truth on the
     realized sample.  CLEAN.
p04: safety(off) vs default, 7 working configs: mreldif(b)=mreldif(V)=mreldif(M)=0.
     CFG5 (contreat+lattice(all 3), i=unit) is the one live r(322): default r(322),
     safety(off) returns att -2.026 and ame exactly 0 (structural zeros -- documented
     at sthlp:1209, so not filed).
p05: no-never-treated bypass gives the CORRECT number (1.6666667 == truth on e(sample),
     240/800 obs dropped and disclosed).  esf pre/post bypass -> bare r(2001).
     controls(policy2) collinear co-timed policy does NOT trip r(322).
p06: unconditionalse+att_i bypass SE 0.0505284 vs 0.0506031 robust vs 0.0491569 cluster
     (no caution printed under the bypass, but no visible damage on this design).
     FRACTIONAL ONSET bypass == the ttre+1 answer, manufacturing a 14-sigma pre-trend.
     Negative weights on missing-y rows: r(411) default, rc 0 under safety(off).
p07/p08: si8 maxvar preflight FALSE POSITIVE, at maxvar 2048 (1900 units, demand 2178)
     and at Stata's DEFAULT maxvar 5000 (4800 units, demand 5078).  Both refused
     r(907) under default safety and both rc 0 under safety(off).
p09: diagnosis pinned -- BOTH engine paths run at maxvar 2048 on the 1900-unit panel,
     so si8ilev+si8tlev are not variable demands at all.
     Cost-audit check: F13 still fires on the observed run with no engine (r459) and
     under bootstrap(standard) reps(5) (r459).  permute() leg NOT tested (my syntax
     error: permute() takes shuffle|probabilities).
p10: DECISIVE engine split, 2300 units, maxvar 2048, quoted demand 2578:
     default safety r(907); cluster(unit) safety(off) rc 0 N=18400 (reghdfe, FALSE
     POSITIVE); customvce(cluster unit) safety(off) r(907) (reg path, demand real).
     => the si8 formula must be engine-conditional, not "drop the level counts".
