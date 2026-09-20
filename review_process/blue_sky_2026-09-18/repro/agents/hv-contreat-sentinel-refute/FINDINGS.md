# SKEPTIC report - hv-contreat-sentinel

All evidence is from MY OWN data and MY OWN do-files. I read the prober's scripts
but ran none of them as evidence.

My panel: hospital-year administrative shape. i = STATE (12), 30-60 hospitals per
state, T = 2008..2019, staggered expansion 2012/2014/2016, three never-expanding
states, dose = expo_h * 1{year >= expyr_s} (canonical exposure x post dose).
Self-contained builders in every repro file.

---

## VERDICTS ON THE SEEDED CLAIMS

### hv-F5-verify -> PARTIAL (severity low -> none-as-a-defect)
Repro: repro_seeded_claims.do (tags F5_*).
Reproduces exactly: 4 of 16 e(V) diagonals = -1.0000e+100 with b = exactly 0 on
pretreatment_ame and ES_ame_pre4/pre3/pre2; e(V_es_ame) diag = 0,0,0,x,x,x;
e(es_ndegen_ame) = 3. Also under plain iid (r1_own.log) and customvce(robust).

CORRECTION: the prober frames this as firing "outside the ruled AMESE R1
cell-containment geometry" - true, but the wrong governing document. The two
representations disagree BY DOCUMENTED DESIGN, in three places:
 - sthlp:903  "expected zero doses in untreated periods, as under
   contreatcontrols(untreatedZero), do not trigger this [dose-slope screen]"
 - sthlp:1325 degenerate horizons of e(V_es_type) "carry zero (structural-zero)
   or missing diagonal rows, flagged by a Caution and by e(es_ndegen_type)"
 - sthlp:1616-1617 e(V)'s placeholder "is stored as a large negative number ...
   code that reads e(V) cells directly ... should treat negative diagonal entries
   as not identified"
One export is documented to carry 0; the other is documented to carry a negative
placeholder. That is two documented conventions meeting, not a defect.

I also found a THIRD announcement the prober missed. Every affected run prints
FIRST: "Caution: a treatment period coefficient interacted with something (likely,
a term implementing a continuous treatment variable ...) was dropped for
multicollinearity ... Note: the multicollinearity caution above was triggered 41
times total." So three cautions fire, not two.

### hv-adv-ame-rank -> CONFIRMED (severity medium; disclosure-only, as claimed)
Repro: repro_seeded_claims.do (tags ADV_*).
I transcribed the alt==0/tt==0 branch of wooldid's OWN kernel wooldid_advW
(ado:14105-14155) - wooldid's Mata functions are not addressable from a caller's
session, see the harness note - and validated the transcription by reproducing
wooldid's PRINTED adversarial p to <1e-12 for BOTH families on BOTH VCEs first.
 - ame, cluster(hid): W(full 6x6) = +1.90489eed2e491X+001, q=2;
   W(post-only 3x3) = +1.90489eed2e491X+001, q=2. BITWISE EQUAL.
   my p = printed p = 0.210801747162294.
 - ame, iid: W = +1.a35bd8c4cc796X+000 both, q=2 both, p = 0.440917229438995.
 - att comparator, same runs: W full != W post, q=5 vs 2, p 0.534 vs 0.326
   (cluster) and 0.917 vs 0.696 (iid).
The premise is exact: the WHOLE first three ROWS of e(V_es_ame) - diagonal AND
the cross-covariances with the post block - are exactly 0, so the eigen-truncation
provably cannot retain any pre-period direction. The rank rule is documented
(sthlp:750-754, 389); what is missing is any signal that it FIRED, while the two
sibling pre-trend diagnostics on the same call do signal.

### hv-honestdid-Minvariant -> PARTIAL (severity high -> medium)
Repro: repro_seeded_claims.do (tags HONEST_*, CTRL_*).
The behaviour reproduces: delta(rm), cluster(hid), the ame CI is identical at
M = 0, 0.5, 1.0 (iid [-0.3430373248, 4.828910034] at all three; cluster
[-0.2242453055, 4.659319125]); att on the same run widens normally.

BUT THE ATTRIBUTION IS WRONG, and my control refutes it. I built
b = (0,0,0,1.2,1.0,0.8) with an exactly-zero 3x3 pre block and an ordinary post
block BY HAND, with no wooldid anywhere in the session, and fed it to honestdid:
CTRL_M_INVARIANT_no_wooldid=1 - [0.6186186186, 1.783783784] at M = 0, 0.5 and 1.0.
The same b with a non-degenerate pre block widens
(CTRL2_M_INVARIANT_nondegenerate_preblock=0). honestdid delta(sd) dies r(3498) on
my hand-built matrix exactly as on wooldid's. The M-invariance is honestdid's
behaviour given ANY degenerate pre block; wooldid computes none of it.

What wooldid actually does (the prober lists these, then discounts them):
 - sthlp:1321 makes "a covariance with no structurally degenerate rows (see
   e(es_ndegen_type))" a REQUIRED condition of the positional handoff;
 - sthlp:1322 states that geometry violations hand honestdid the wrong geometry
   SILENTLY - silence at the point of use is the documented character of a
   positional handoff;
 - ado:4160 prints, at estimation, a caution that NAMES THE CONSUMER: "...inspect
   the diagonal before joint inference - honestdid/lincom on those rows is not
   meaningful.";
 - e(es_ndegen_ame)=3 is posted and survives estimates save/use.
Documented precondition + named-consumer caution + machine-readable flag. "high"
overstates it. Medium: the residual defect is a UX one - the only protection is a
caution that has scrolled by the time the handoff is run.

NEW, and it sharpens the case: the OTHER handoff the same help advertises
(sthlp:802, Roth's pretrends) FAILS LOUDLY on identical inputs - PRET_ame_rc=198,
"PreTrendsPower(): Power level .5% too high ... PreTrendsBisect(): Invalid
starting points" - while working on att (PRET_att_rc=0, slope for 50% power
= .10215).

### hv-noVsentinel-marker -> CONFIRMED (severity low, as claimed)
Repro: repro_seeded_claims.do (tags SENT_*). Code: ado:4054-4085 displays the
count and posts nothing. After estimates save + clear all + estimates use: 4
sentinel diagonals present, e(es_ndegen_ame)=3 preserved, and the full ereturn
list (printed in the log) carries no scalar naming the 4.

---

## NEW FINDINGS (mine)

### A. [CRITICAL] contreatcontrols(untreatedZero) is silently ignored on every non-reghdfe estimation path
Repro: finding_A_untreatedzero_regresspath.do (self-contained; read the VERDICT_*
block at the end of the log).

With contreat() + contreatcontrols(untreatedZero), the regress-based path
(customvce(anything), unconditionalse) and the poisson path estimate the
contreatcontrols(BASIC) model instead. One panel, one seed:

  independent hand-built ORACLE (untreatedZero's own equation, sthlp:997)
                                          ame 0.9127941639   att -0.1141627684
  DEFAULT path + untreatedZero            ame 0.9127941639   att  0.8064147779
  robust (still reghdfe) + untreatedZero  ame 0.9127941639   att  0.8064147779
  customvce(ols) + untreatedZero          ame 0.4766147418   att  0.3661964684
  customvce(robust) + untreatedZero       ame 0.4766147418   att  0.3661964684
  customvce(cluster hid) + untreatedZero  ame 0.4766147418   att  0.3661964684
  unconditionalse cluster(hid) + uZ       ame 0.4766147418   att  0.3661964684
  DEFAULT path + BASIC                    ame 0.4766147418   att  0.3661964684

VERDICT_cv_equals_BASIC_ame = 1, VERDICT_cv_equals_BASIC_att = 1,
VERDICT_un_equals_BASIC_ame = 1 (reldif < 1e-9). Every post-treatment ES horizon
moves with them; the PRE-treatment rows are identical to 10 decimals, e(N) is
identical (8640), e(r2) is identical (0.236832847), e(df_r) is identical, and no
caution distinguishes the runs. The ratio here is 0.522 on ame and 0.454 on att;
it is NOT a fixed shrinkage - on a differently-seeded draw of the same design the
regress path returned main_ame = -0.3095 (r4_boot.log, section E), the wrong sign
relative to the design's known truth of +0.8. Caveat: the DEFAULT comparator was
not run on that particular draw, so the sign flip is against the known truth, not
against a measured default on the same data.

Controls that hold: a BINARY run agrees across the two paths exactly
(VERDICT_binary_paths_agree = 1); a dose nonzero in every period with basic agrees
across the two paths exactly (r6_customvce.log, section 4). So the defect is
specific to contreat() + untreatedZero + a non-reghdfe path.

POISSON: poisson + untreatedZero and poisson + basic return identical numbers
(att .02056812, ame .02649715). SCOPE HONESTLY: on this outcome BOTH poisson runs
printed "Caution: initial estimation of ppmlhdfe failed" and took the rc-3200
explicit-regressor fallback, so the measured equality is on that fallback; the
ppmlhdfe-absorb path carries the same term through contreatcontrolinserterABS by
the ado:1234 condition but I did not measure it, and I built no poisson oracle.
Claim for poisson: "no correct path" is code-implied plus the measured uZ==basic
equality on the fallback, not oracle-confirmed.

CAUSE, read in the code. wooldid.ado:1234
  if strpos("`contreatcontrols'","basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
      ... local ccringerhead  i.`null'##`ringerpoly' ...
inserts the pooled baseline dose control; its own comment says it is there to
"correctly manipulate the behavior of reg and poisson's colinearity drop rules"
for the basic-type controls. The untreatedZero branch at ado:1313
  if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
  }
rewrites contreatinj only - its comment says untreatedZero works by "just not
having contreatcontrols of any kind" - but never clears
contreatcontrolinserterFE/ABS, so on those paths the baseline dose slope that
untreatedZero exists to omit is put back. Under an exposure x post dose that term
is exactly collinear with the treated-cell dose terms, so the FIT is unchanged
(identical r2/rss/df) while the ame/att decomposition is not.

Realism: contreat() + contreatcontrols(untreatedZero) is the package's own worked
example (sthlp:1551); customvce() is the documented route to any vce wooldid does
not offer natively; unconditionalse is a documented feature; poisson with a
dose x post exposure is a standard count-outcome design. Nothing in the help says
untreatedZero behaves differently on those paths (sthlp:264, 270-273 describe
customvce() as a variance switch only), and no caution is printed.

### B. [low] the two advertised third-party handoffs disagree about refusing
On identical degenerate-ame inputs, the advertised pretrends handoff (sthlp:802)
refuses with r198 and the advertised honestdid handoff (sthlp:1559-1567) returns
an M-invariant robust set. Evidence in repro_seeded_claims.log (PRET_*, HONEST_*).

---

## PROBED AND CLEAN (negative results)

- bootstrap() does NOT remove the disclosure. bootstrap(standard) reps(100) +
  untreatedZero + cluster(hid): rebuilt e(V_es_ame) still carries 3 of 6 exact-zero
  diagonals, e(es_ndegen_ame)=3 still posted, the caution still fires,
  e(boot_ves_failed) empty, e(wdid_vcvstate)=bootstrap, 4 e(V) sentinels still
  present (r4_boot.log, section A). honestdid on the BOOTSTRAP matrix is
  M-invariant too (HON_BOOT_INVARIANT=1) - same honestdid property, same wooldid
  disclosure. Nothing is lost on the bootstrap path.
- PS-62 splice integrity after the partial-identification rescue: e(V)'s ES_ame
  post block vs e(V_es_ame)[4..6,4..6] gives mreldif 2.6e-11 (iid) / 2.9e-11
  (cluster) - not bitwise, but inside the project's 1e-10 indifference threshold.
  Not filed. (r2_adv.log, SPLICE_*.)
- subgroups() + untreatedZero exits r(459) on this design before any per-subgroup
  V_es_ame_s# is built, so that caution path is UNTESTED, not clean
  (r4_boot.log, section D).

## HARNESS NOTES (not wooldid defects; they will cost the next agent time)
- wooldid's own Mata functions are NOT in Mata memory after a successful run:
  mata: mata describe after qui wooldid ... lists only the caller's own functions
  plus a few tempname scalars (r2b_matatest.log). Transcribe wooldid_advW; you
  cannot call it.
- cap dis "X=" X prints NOTHING AT ALL when X is undefined - the literal is lost
  with the line - so sentinel prints built that way are invisible rather than
  wrong. Print from inside Mata with printf.

## ONE HYPOTHESIS I HELD AND DISPROVED (closes a hunting ground)
I expected the degenerate-row disclosure to live only in the analytic post block
(ado:4130-4165), so that under bootstrap() - where e(V_es_*) is REBUILT from the
replications and the pre-period ame is 0 in every replication - the caution and
e(es_ndegen_ame) would disappear while the exact-zero pre block survived. That is
FALSE. Measured (r4_boot.log section A): bootstrap(standard) reps(100) +
untreatedZero + cluster(hid) still prints the caution, still posts
e(es_ndegen_ame)=3, still carries the 4 e(V) sentinels, and leaves
e(boot_ves_failed) empty. The bootstrap path is covered.
