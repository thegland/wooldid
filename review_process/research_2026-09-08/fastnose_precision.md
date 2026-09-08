# fastnose_grid: the c_lattice precision divergence

**Topic.** fastnose_grid c_lattice precision loss: why full-mode and fast-mode diverge 1e-9 relative on the AME family under lattice(all 3)

## Finding

The divergence has nothing to do with `nose`, with the cell collapse, or with aweights. It is caused
by the two legs fitting the stage-1 regression through **two different numerical routes**. On this
configuration the full leg trips the AMESE-R4 all-zero-VCE remedy (wooldid.ado:1394-1414) and
**refits with `regress` on explicit dummies** (355 columns); the fastnose leg skips that whole block
— it is gated on `if "`fastnose'" == ""` — and keeps the original `reghdfe` partialled-out fit (328
columns). Every reported number in the full leg therefore comes from a different solve than the
fastnose leg's, and the AME family (which reads the dose slopes) inherits ~1e-10 absolute round-off
from it while the ATT family (which reads a much better-determined combination) inherits only
~1e-13. The important inversion for the author: **the fastnose number is the more accurate one.**
The `regress` explicit-dummies route is intrinsically noisy on this design — merely reordering the
rows of the input data moves its own answers by 7.9e-11 relative on `main_ame` and 1.6e-8 relative
on `pretreatment_ame`, i.e. *more* than the full-vs-fast gap the component is flagging — whereas the
reghdfe route reproduces itself to ~1e-13/1e-14 under the same reordering, and a second,
independently parameterised reghdfe route (`ccnoabsorb`) agrees with it to 5e-14 while sitting
4.15e-11 away from the refit. The component's implicit assumption that the full leg is the reference
and fastnose the approximation is backwards in this class. Separately, the component header's
"~1e-15 abs / ~1e-13 rel" figure is not three orders wrong: it describes genuine nose-only paths
(the subgroup pairwise-contrast route, where `sgpost` is cleared under fastnose), and on the plain
contreat paths `nose` moves literally nothing — full and fastnose agree to exact bitwise zero once
the refit is taken out of the picture. A second, undocumented mechanism is firing, not a mis-sized
bound on the first.

## Mechanism

Chain, each link measured. (1) TRIGGER. `contreat()` + a cluster-type VCE arms the cell-containment
predicate (wooldid.ado:1373-1383): with `cluster(i)` and lattice, 40 treated i-t cells lie wholly
inside one cluster, so `amescrub` = 1 and reghdfe returns an all-zero e(V). `wooldid_amescreen,
mode(vzero)` (wooldid.ado:1395) then detects the zero VCE and, because the narrowness conditions
hold (single cluster var, no fe(), no ccc_absorb, 10+10 levels well under the 450 cap at
wooldid.ado:1405), the run refits at wooldid.ado:1407 with `reg y i.i i.t ... , cluster(i)` to
recover the SEs. That entire block is inside `if "`fastnose'" == ""`, so fastnose never sees it. (2)
WHY THE ROUTES DISAGREE AT ALL. Each treated cohort-by-period cell carries an intercept delta_g and,
under contreat, a dose slope beta_g. `lattice(all k)` splits each cell into k dose bins, giving a
(cell x bin) intercept/slope pair identified off the dose spread inside that bin of that cell. In
dta_cont each treated (i,t) cell holds 10 observations; at k=3 six of the 87 treated (cell,bin)
groups hold a single observation and 23 hold two or fewer. In such a group the intercept and the
slope are near-perfectly collinear — a ridge. Along that ridge the two solvers land in different
places: the individual stage-1 coefficients differ by up to 8.1e-8 absolute on the cell intercepts
and 8.4e-9 on the dose slopes. (3) WHY ONLY THE AME FAMILY CROSSES THE GATE. Both coefficient blocks
are noisy; the discriminator is the functional read off them. The ATT reads dydx(targetflag) =
delta_g + beta_g*xbar_cell — the fitted value at the cell's mean dose, i.e. the *well-determined*
direction across the ridge — so the intercept and slope errors cancel almost completely, down to
~1e-13 absolute. The AME reads beta_g alone, which is the ridge direction itself, so it cancels only
to ~1e-10. That is the whole 1000x gap between the two families, and it is why the pre-period AME
entries (magnitudes 0.006 to 0.04) cross a 1e-10 relative gate while the order-10 ATT entries do
not. (4) WHY IT SCALES WITH THE LATTICE. More bins means narrower within-(cell,bin) dose ranges and
more singleton groups, so a sharper ridge. Singleton (cell,bin) group counts and the worst AME
absolute divergence move together across lattice(all 2/3/4/5): 0 -> 4.5e-12, 6 -> 8.9e-11, 20 ->
1.95e-10, 35 -> 4.95e-10. With no lattice the refit does not even fire and the two legs agree to
exact zero. (5) WHAT IT IS NOT. Not `margins ... nose`: nose leaves dydx point estimates bit-
identical in five synthetic shapes and, inside wooldid, gives exact full-vs-fast zeros on every
contreat config where the refit does not fire. Not the cell collapse: `lattice()` is in `fcbad`
(wooldid.ado:4561) so `fcon` = 0 and both legs run margins on the full data. Not aweights: dropping
`[aw=w]` leaves the same three violations. Not a reghdfe convergence tolerance: the fastnose leg is
already converged at regtol(10) and does not move through regtol(14), and the full leg is regtol-
invariant because its numbers come entirely from the `regress` refit.

## Evidence

- REPRODUCTION (e1.do). Full vs fastnose on c_lattice: pretreatment_ame -0.016901229872367567 vs
  -0.016901229940361192 (abs 6.80e-11, rel 4.02e-09); ES_ame_pre2 abs 4.72e-11 rel 7.74e-09;
  ES_ame_pre1 abs 8.86e-11 rel 2.23e-09. Note main_ame also diverges 4.15e-11 on an order-1 value —
  it passes only because 4.15e-11 < the 1e-10 relative arm. The ATT family diverges only 1.5e-13
  absolute.
- `nose` IS NOT THE CAUSE (e2.do). `margins, dydx(x) df(100)` vs `... nose` after regress: bitwise
  identical (%21x hex equal, absdiff 0.000e+00) in all five shapes tested — plain `y x z`;
  `i.g##c.x`; `i.g##c.x [aw] subpop()`; the same plus a duplicate 'ringer' copy of x in the model;
  and `over(g)`.
- THE REFIT IS THE CAUSE (e3.do, verbose). Full leg prints both 'Caution: with standard errors
  clustered on i, 40 treated i-t cell(s) lie completely within a single cluster' and 'Caution: the
  estimation command returned an all-zero variance-covariance matrix; wooldid is re-estimating
  through the explicit-dummies fallback' (wooldid.ado:1406). Fastnose leg prints neither. Stored
  stage-1 model: 355 columns (full) vs 328 (fastnose).
- NECESSARY AND SUFFICIENT (e4.do, 6 configs with refit detection). c_lattice_robust (identical but
  `robust` instead of `cluster(i)`, so no containment screen and no refit): worst abs 3.55e-15,
  worst rel 5.20e-12, 0 violations — 9 of 12 e(b) entries agree to exact 0.000e+00. c_lattice_noaw
  (no aweight, refit still fires): the same 3 violations. c_basic and c_basic_aw (containment fires,
  refit does NOT — reghdfe returned a nonzero VCE): worst abs 3.55e-15, 0 violations, 11 of 12
  entries exactly equal. c_uz (refit DOES fire but no lattice): worst rel 6.65e-13, 0 violations.
- NOT A TOLERANCE ARTIFACT (e5.do). Full leg identical at regtol(9) and regtol(13) to all 17 digits
  — its numbers are the `regress` refit's, not reghdfe's. Fastnose leg: regtol(9) 1.0000063689476888
  -> regtol(10) 1.0000063689477148, then regtol(11), (12), (13), (14) all identical to regtol(10).
  reghdfe has converged; the residual 4.14e-11 gap to the refit does not shrink.
- ACCURACY ARBITER — ROW-ORDER NOISE FLOOR (e7.do). Re-running each route on the same data in three
  different row orders (gsort -i -t, sort intensity, random permutation). REGRESS REFIT route
  main_ame: 1.0000063689891556 / ...9543504 / ...9201171 / ...9106269 — internal spread 7.85e-11;
  pretreatment_ame internal spread 2.66e-10 (1.57e-08 relative); ES_ame_pre1 spread 5.06e-10
  (1.28e-08 relative). REGHDFE route under the same reorderings: main_ame spread 4.6e-14,
  ES_ame_pre2 spread 2.3e-14, pretreatment_ame spread 2.9e-13. The refit route's own
  irreproducibility exceeds the full-vs-fast gap the component flags; the reghdfe route's does not.
- ACCURACY ARBITER — INDEPENDENT ROUTE AGREEMENT (e7.do). `ccnoabsorb` moves the dose-control terms
  out of reghdfe's absorb() into its main equation (wooldid.ado:1183-1186), a differently-
  parameterised reghdfe solve. BNA_H main_ame = 1.0000063689476359 vs BHDFE 1.0000063689476888 —
  agree to 5.3e-14 — while both sit 4.15e-11 from the refit's 1.0000063689891556. Two controls
  confirm the refit is the full leg's whole story: `customvce(cluster i)` (which forces the same
  explicit-dummies `regress` path in BOTH legs, wooldid.ado:1327-1331) reproduces the full-leg
  numbers to 0 or 3.5e-16, and full+ccnoabsorb reproduces them exactly.
- MECHANISM AT THE COEFFICIENT LEVEL (e6b.do, _cons excluded as non-comparable across
  parameterisations). Of 326 shared nonzero stage-1 coefficients: 106 dose-slope (c.intensity)
  terms, worst divergence 8.383e-09 absolute / 2.849e-07 relative; 116 cell-intercept terms, worst
  8.128e-08 absolute / 2.607e-07 relative. Dose-slope coefficient magnitudes run 1.06e-03 to 1.96.
  Both blocks are noisy at 1e-8; the reported ATT still cancels down to 2.5e-13 while the reported
  AME only reaches 8.9e-11.
- CONDITIONING SCALES WITH THE LATTICE (e8.do + e10.do). No lattice: refit does not fire, AME
  divergence exactly 0.000e+00. lattice(all 2): refit fires, AME worst abs 4.50e-12, 2 violations.
  all 3: 8.86e-11, 3 violations. all 4: 1.95e-10, 4 violations (worst rel 1.51e-08). all 5:
  4.95e-10, 5 violations — main_ame itself now violates at 5.78e-09 relative. Singleton treated
  (cohort,period,bin) group counts on dta_cont over the same sweep: 0, 6, 20, 35 (groups of size
  <=2: 2, 23, 56, 79).
- DATA STRUCTURE (e9.do). dta_cont: 1000 obs, 10 cohorts x 10 periods x 10 obs per cell, one onset
  (ttre=5), intensity mean 9.99 sd 0.495. Every treated (i,t) cell has 10 rows, so at lattice(all 3)
  a (cell,bin) group can hold as few as 1.
- NO e() MARKER. wooldid.ado:4009-4041 posts wdid_fcol, wdid_fcskip, wdid_fastnose, wdid_vcespec and
  the rest, but nothing records that the R4 refit fired. There is no way for the parity component —
  or a user — to detect that the reported estimates came from the explicit-dummies fallback rather
  than from reghdfe.
- NOT A REGRESSION AND NOT COLLAPSE-RELATED. `lattice()` is on the fcbad list at wooldid.ado:4561,
  so fcon = 0 and the verbose line reads 'Wooldid is evaluating the margins calls on the full data:
  the default cell collapse does not apply with lattice()' in BOTH legs. The task's isolation
  already established the failure predates the current campaign; the refit block is long-standing.
- BOOTSTRAP/PERMUTE ARE ALREADY PROTECTED. Replications run under fastnose (wooldid.ado:126, 302),
  i.e. on the reghdfe route, while an observed run in this class reports the refit route. The rep-0
  fastnose comparator design (wooldid.ado:6701-6706 for permute, 7767-7769 for bootstrap) re-
  estimates the observed point under fastnose precisely so the comparison is bitwise consistent. The
  R4 refit is a concrete reason that design is load-bearing, not merely defensive.

## Recommendation

Not a code fix to the estimation path, and not a tolerance keyed to `lattice()`. Three steps, in
order. 1. ADD A NUMBER-NEUTRAL MARKER (the one code change I do recommend). One line beside the
other hidden posts at wooldid.ado:~4020, e.g. `ereturn hidden local wdid_amerefit "1"` set when the
wooldid.ado:1407 refit executes. It moves no number, and it is the only way anything downstream can
tell that the reported estimates came from the explicit-dummies fallback. Users deserve it too:
right now the only trace is a Caution line that a `qui` run swallows. 2. KEY THE COMPONENT'S
TOLERANCE ON THAT MARKER, NOT ON THE CONFIG NAME. The class is "R4 refit fired" = contreat() +
cluster-type VCE + containment + all-zero reghdfe VCE + single cluster var + no fe()/ccc_absorb +
<=450 i/t levels. `lattice()` and `[aw=w]` are incidental members: c_uz is in the class and passes
at 6.65e-13, c_basic_aw has lattice's aweight but is not in the class and agrees exactly. For
configs where the marker is set, relax the relative arm to 1e-6 while keeping the 1e-12 absolute
arm; worst measured full-vs-fast AME relative divergence over lattice(all 2..5) is 1.51e-08, so 1e-6
is ~66x headroom — and that headroom is warranted, because the refit route's own row-order
irreproducibility is already 1.6e-08 relative, i.e. the quantity being compared has no more digits
than that. Where the marker is not set, keep the existing 1e-10 arm: it holds at exact zero on every
non-refit contreat config measured. 3. KEEP TIGHT COVERAGE OF THE LATTICE CODE PATH by adding a
`c_lattice_robust` config (the same call with `robust` in place of `cluster(i)`) at the unrelaxed
tolerance. That exercises lattice + contreat + aweight + fastnose end to end and, measured, passes
at worst abs 3.55e-15 / worst rel 5.20e-12 with 9 of 12 entries bitwise equal — so the relaxation in
step 2 costs no real regression sensitivity on the lattice path. Also worth one line in
wooldid.sthlp's fastnose entry: the documented divergence there describes the nose-only paths, and
should say that a separate and larger divergence (up to ~1e-8 relative on ame-family estimates,
growing with lattice() bin count) arises when the all-zero-VCE refit fires, and that in that case it
is the FULL-mode estimate carrying the round-off, not the fast-mode one. HELD FOR THE AUTHOR, NOT
RECOMMENDED BY ME (number-moving, per the standing workflow rule): the refit exists to recover
standard errors, but it currently replaces the point estimates too. An alternative is to keep
reghdfe's b and take only V from the refit. That would move published full-mode ame numbers by
~1e-10 relative and would break the b/V-from-one-fit invariant, so it is a ruling, not a bug fix. My
own view is that ~1e-10 on a point estimate is statistically weightless and the coherence of b and V
is worth more, so I would leave the estimation code alone — but the author may weigh the
reproducibility argument differently now that the refit route is documented as the noisier one.

## Open questions

- Does the refit actually earn its cost here? The full leg reports 'standard errors were not identified for 6 of 12 reported estimates', so the refit recovered 6 SEs — presumably the ATT family, since amescrub nulls the ame se/t/p/ci cells regardless. I did not verify which 6, nor whether those recovered cluster-robust SEs are trustworthy under the very containment that zeroed the original VCE. If they are not, the refit's scope is worth revisiting on its own merits, independently of this parity finding.
- I did not establish an absolute reference value. The accuracy ranking (reghdfe better than the refit) rests on two strong but indirect arguments — the refit's own row-order irreproducibility, and the agreement of two differently-parameterised reghdfe solves to 5e-14 — not on a high-precision solve of the same normal equations. Both reghdfe variants share an algorithm family, so a common-mode error in reghdfe cannot be fully excluded from what I ran.
- The 1.5e-08 worst-measured relative divergence is measured, not bounded. It comes from one dataset (dta_cont, n=1000, 10x10 cells) at one lattice size. Nothing I ran establishes how bad the refit's round-off can get on a design with more singleton (cell,bin) groups, more cohorts, or a dose with a heavier tail. A recommended tolerance of 1e-6 is a judgment about headroom, not a proof.
- Whether other refit-class configs exist elsewhere in the battery that are currently passing by luck. amescrub requires contreat(), so the class is contreat-only, but I only surveyed the six configs in e4.do plus the lattice sweep — not the full suite.
- Whether the same b/V-route asymmetry shows up anywhere the full leg's e(V) is consumed downstream (honestdid, pretrendpower, adversarial): those read the refit's covariance against the refit's coefficients, which is internally consistent, but any cross-run comparison against a fastnose-derived quantity inherits the 1e-10 offset.

## Proposed simulation

No Monte Carlo is needed to settle the cause — it is settled by the direct experiments above. What a
follow-up SHOULD run, if the author wants the accuracy ranking nailed to an absolute reference
rather than to reproducibility arguments, is a deterministic numerical study, not a sampling study:
DESIGN A — quad-precision reference (settles open question 2). Take an instrumented COPY of
wooldid.ado in scratch (never the repo build) and, immediately after the wooldid.ado:1322 reghdfe
call, dump the realised design to disk: y, the weight, e(sample), and the tempvars appearing in the
verbose model line (ifvar __000000, cell _i __000004, relative period _t __000005, targetflag
__000006, lattice var __000007, ringer __000008, and the contreat-control group __000009). Outside
wooldid, materialise the full explicit-dummies design with `fvrevar`, drop collinear columns with
`_rmcoll`, and solve the weighted normal equations in Mata with `quadcross` (quad-precision
accumulation) plus `qrsolve`. Form the same AME functional — the estimand weights times the cell-by-
bin dose slopes — from (a) the quad solution, (b) reghdfe's b, (c) the refit's b. Report each
candidate's distance from the quad reference and, separately, the relative normal-equation residual
||X'W(y - Xb)|| / ||X'Wy||. Expected outcome given e7: reghdfe within ~1e-13 of the reference, the
refit ~4e-11 away. This is a single deterministic run, not a simulation. DESIGN B — bound the round-
off, which is where randomisation genuinely helps (settles open question 3). Fix the estimand and
vary only the numerics and the conditioning. Grid: n_cohorts in {10, 20, 40} x periods in {10, 20} x
obs-per-cell in {4, 10, 25} x lattice(all k) for k in {2,3,4,5,8} x dose sd in {0.1, 0.5, 2.0}
(relative to its mean, to sweep the within-(cell,bin) spread that sets the ridge sharpness). For
each design point draw 25 datasets; for each dataset run wooldid full and fastnose, and additionally
re-run each leg under 5 random row permutations. Record per replication: the full-vs-fast relative
divergence on every ame- and att-family e(b) entry; each route's own across-permutation spread; the
singleton and size<=2 (cell,bin) group counts; and whether the R4 refit fired. Report the 99th
percentile of the full-vs-fast ame relative divergence as a function of the singleton-group share —
that is the number a tolerance should be set from, and it also tests directly whether singleton
share is the right conditioning summary. Roughly 2250 wooldid pairs plus permutation re-runs; on
n<=5000 designs that is hours, not days, and it should be run in the background with the
4-concurrent-Stata cap respected.
