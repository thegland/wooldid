# SKEPTIC pass on the safety() findings -- NOTES

Build asserted in every log (`which wooldid` -> .../bluesky/build/wooldid.ado).
Host Stata reports `c(edition)=BE`, default `c(maxvar)=5000`.
TRAP found: `set maxvar` refuses (r4) while Stata matrices or Mata exist -- `clear`
alone is not enough; needs `clear matrix` + `clear mata` (or `clear all`).

## Verdict summary
| id | verdict | severity (theirs -> mine) |
|----|---------|---------------------------|
| SAFETY-01 | CONFIRMED (stronger than claimed) | high -> high |
| SAFETY-02 | CONFIRMED (doc divergence real; no number moves) | medium -> medium |
| SAFETY-03 | PARTIAL -- mechanism and "silently" both refuted | medium -> low |
| SAFETY-04 | CONFIRMED | low -> low |
| SAFETY-05 | PARTIAL -- real, but the rc is not r(2001) | low -> low |

## My own runs
- r1a_si8_defaultmaxvar.do  4800 units / 8 yrs / 20 treated, DEFAULT maxvar 5000:
  default safety -> r(907) quoting demand 5107; safety(off) -> rc 0, N=38400, att 1.5271631.
- r1b_si8_control.do  same data at maxvar 20000 under DEFAULT safety:
  att hex +1.86f429091af9eX+000 and se hex +1.3c63b9e9bc101X-007, BIT-IDENTICAL to the
  safety(off)-at-5000 run.  The refusal was purely false.
- REPRO_SAFETY01_si8_falsepositive.do  the self-contained version, three parts:
  P1 4000 units pass / 4900 units refuse with the SAME 20 treated units (the trigger is the
  NEVER-TREATED CONTROL count); P2 ATT_BITIDENTICAL=YES SE_BITIDENTICAL=YES; P3 at the
  default ceiling with 5200 units: default r(907), safety(off)+reghdfe rc 0 (false positive),
  safety(off)+customvce(cluster unit) r(907) (the demand is REAL on the reg path).
- r2/REPRO_SAFETY02: float caution printed 2x binary, 2x contreat-without-ctw, 0x with
  contreatwithin.  i coded 0: r(411) default, accepted under ctw, att/att_i/ame/ame_i ALL
  bit-identical to the 1-indexed control -> NO number moves.  Non-integer i -> r(452).
- r3/REPRO_SAFETY03: ttre+0.5 safety(off) == ttre+1 recode bit-for-bit, BUT the ES rows are
  labelled pre1.5 / pre.5 / post.5, and ttre+0.7 gives pre1.7 / pre.7 / post.3 -- there is no
  rounding at all.  ttre+0.3 dies r(452); esprelength(0) dies r(452).
- r4_own.do BLOCK E: no-control-period drop under safety(off) -- the caution is suppressed but
  "N = 360 (40 obs dropped from initial sample)" and e(N_dropped)=40 print either way, and the
  att is bit-identical.  sthlp's "never suppressed a disclosure" HOLDS.  Negative result.
- r4_own.do BLOCK F: after the mid-run r(907) under safety(off)+customvce, the dataset varlist
  is unchanged (F_LEAKED=NO) and `est dir` is clean.  Negative result.
- r5_f13_and_esf.do BLOCK G: F13 onset-constancy SURVIVES contreatwithin (r459 with no safety()
  typed).  Its safety(off) bypass moves no number on the DEFAULT pooled-reference spec, because
  the onset-onward-only rows land in the reference pool either way.
- r5 BLOCK H: the esfixedbaseperiod safety(off) degradation is r(198) "invalid syntax" with an
  ES window and r(461) "no observations in subpop()" without -- NOT r(2001).

## Added after the review pass
- r8_calib.do block D: calibration arm for SKEP-01.  20 treated of 60 units, i=unit,
  cluster(unit): SE/SD 0.508, coverage 0.733.  Against arm A (20 of 300): 0.212 / 0.367.
  The gradient confirms the control-to-treated-ratio mechanism.
- r8_calib.do block E: `bootstrap(standard) reps(2) cluster(unit)` on the 4800-unit panel
  ALSO exits r(907) quoting demand 5107 -- the si8 false positive refuses the very remedy
  sthlp:251 prescribes for i=unit.
- FINDINGS.md corrected: the "Use set maxvar" line in the r(907) output is Stata's own
  boilerplate, not part of wooldid's message.
