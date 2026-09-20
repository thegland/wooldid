# hv-permsave-comma-refute — SKEPTIC report

Pinned build asserted in every run (`which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado).
Own data, own do-files, own subdirectories w1/..w9/. I read the prober's do-files but never ran them.

---

## VERDICTS ON THE SEEDED FINDINGS

### hv-permsave-comma-01 — CONFIRMED (severity medium, unchanged)
Repro: r1_permsave_comma.do (workdir w1/).

| test | command | rc | effect |
|---|---|---|---|
| T-A | `permute(shuffle) permreps(20) permsave(pd, replace)` | **0** | writes `pd, replace.dta` (80,453 bytes, 21 obs, `_rep main_att onset_1..onset_90`) |
| T-B | `bootstrap(standard) reps(20) bootsave(bd, replace)` | **198** | refused: "bootsave() takes a plain file name and no comma; replace is a separate option of wooldid, so write bootsave(bd  replace) replace." Nothing written. |
| T-C | T-A again, verbatim | **602** | "permsave(pd, replace) - pd, replace.dta already exists in the working directory; add the replace option to overwrite it." (the user typed `replace`) |
| T-D | `permsave(pd) replace` | 0 | writes `pd.dta` |
| T-E | `use "pd, replace.dta"` | 0 | 21 obs, contents correct |

Mechanism as claimed. wooldid_permctl (ado:5868-5880) tests `.dta`, `/`, `\` and emptiness;
wooldid_bootctl (ado:6279-6298) tests those plus strpos(...,",") (dated RF-6, 2026-09-13).
The ado's own comment at 5866 — "permsave(): the name rules are bootsave()'s" — is falsified.
I checked the rest of the name-rule surface for asymmetry (r6, S4) and found none: a directory
part and an empty name are both r(198) on BOTH savers. The comma is the only divergence.

### hv-permsave-comma-02 — PARTIAL (absences real; stated mechanism wrong)
Repro: r7_W_lean.do (workdir w7/). r3_bootW.do was the heavy version (esfixedbaseperiod +
adversarial + subgroups); it hit its 15-minute deadline and was killed. r7 is the lean
replacement and answers the same question.

| run | boot_jointW | boot_jointWobs | boot_sgW | boot_sgWobs | jt_family | printed p |
|---|---|---|---|---|---|---|
| A boottstat, jointtests | PRESENT | absent | absent | absent | studentized | joint_horiz = .6923 |
| B same, no boottstat | PRESENT | PRESENT | absent | absent | fixedmetric | joint_horiz = .8462 |
| C boottstat + subgroups | PRESENT | absent | absent | absent | studentized | .7692, allsgssame = .0556 |
| D same, no boottstat | PRESENT | PRESENT | PRESENT | PRESENT | fixedmetric | .8462, .0385 |
| E permute() + subgroups | — | — | — | — | — | perm_jointW/Wobs/sgW/sgWobs ALL PRESENT |

r3_bootW.do's first run did flush before the kill, and supplies the adversarial row that r7
lacks: with `esfixedbaseperiod adversarial subgroups() jointtests bootstrap(standard) reps(30)
boottstat`, `R1_BOOTTSTAT|boot_advW|PRESENT` and `R1_BOOTTSTAT|boot_advWobs|absent`, with
`e(adversarial_pval_att)=.09677419` printed and `e(joint_att_allsgssame)=.` (r3_bootW.log:365-380).

The divergence from sthlp:1340 / sthlp:1391 is real and reproduces. Severity stays low.

CORRECTION — the prober's mechanism is wrong, and it matters for the fix. They wrote that
boot_jointWobs/boot_advWobs "exist at run time as the tempnames wdbjobs/wdbaobs (ado:365,
handed to wooldid_bootsave at ado:437)" and that they "measured the observed statistic in
row 0 of jointW_att = .2704913 in a bootsave file", concluding it is "a two-line posting
omission". wdbjobs/wdbaobs are NOT observed W statistics: they are the observed EVENT-STUDY
COEFFICIENT vectors that wooldid_bootjwobs (ado:15310) and wooldid_bootadvobs consume to
COMPUTE W_obs, which they return in a Mata-set local and then discard. The value the prober
measured in the bootsave file is the first pre-treatment observed coefficient, not a Wald
statistic — see N1 below. (I could not measure their own .2704913: none of the .dta files left
in their directory carries a jointW_ column at all — bdraws/b1/b2 are fixed-metric files whose
columns are `_rep main_att w_1..w_G` (r10_theirfile.do). But `wooldid_bootsavemk` has exactly one
path, so every jointW_* row-0 value in every bootsave file is JO[1,k]; a POSITIVE value like
theirs simply fails to trip the sign tell.) All three studentized gaps
are therefore the same kind of design gap they correctly identified for boot_sgW, not a
posting oversight; a two-line `ereturn matrix` of wdbjobs would post a coefficient vector
under a name that promises a Wald statistic.

### hv-permsave-comma-03 — CONFIRMED (severity low, unchanged)
Repro: r2_dots.do (workdir w2/).
- `bootstrap(standard) reps(20) dots(-5)` -> r(198) "dots() must be a nonnegative integer
  (0 suppresses the progress display)."
- `permute(shuffle) permreps(20) dots(-5)` -> rc 0, no progress line — identical output to
  dots(0), whereas dots(5) prints "Permutation replications (20): . (~14 s projected for 20
  replications)... done".
- Control that the mechanism is the early exit: rseed()'s range check WAS hoisted above it,
  and `permute(shuffle) rseed(-3)` -> r(198).

---

## NEW FINDINGS (mine, not seeded)

### N1 — bootsave() under boottstat writes the WRONG number in the `_rep = 0` row of `jointW_*` and `advW` [HIGH / certain]
Repro: r9_bootsave_jointW_repro.do (workdir w9/); first seen in r8_bootsave_jointW.do.

sthlp:557 (the bootsave(name) entry): "Under boottstat the replication joint-test statistics
are added as jointW_name and the adversarial statistic as advW, WITH THE OBSERVED STATISTICS
ON THE _rep = 0 ROW." ado:14905 says the same in the writer itself.

Measured, `wooldid y unit year ttre, att cluster(unit) jointtests esprelength(2)
espostlength(2) bootstrap(standard) reps(25) boottstat bootsave(bj) replace`:

    REPORTED_joint_horiz = .03846154        (= 1/26: no replication at or above W_obs)
    e(boot_jointW): 25 x 1  colnames = att
    jointW_att[_rep==0]  = -.09398543   (-1.80f6dee4819bfX-004)   <- NEGATIVE
    jointW_att reps      : n = 25, min = .06599616, max = 6.7419682   (all >= 0)
    e(b_es_att)[pre2]    = -.09398543   (-1.80f6dee4819bfX-004)   <- IDENTICAL
    RECONSTRUCTED p from the file = 1        (25 of 25 reps at or above row 0)
    REPORTED p                    = .03846154

and the adversarial twin (same shape with esfixedbaseperiod adversarial):

    ADV_pval_att  = .42307692            (= 11/26)
    advW[_rep==0] = -.17179509           <- NEGATIVE
    advW reps     : n = 25, min = .17107281, max = 37.347159

Mechanism, wooldid_bootsavemk (ado:14919-14931):

    JW = st_matrix(jwacc); jcs = st_matrixcolstripe(jwacc);
    JO = (jobs != "" ? st_matrix(jobs) : J(1, cols(JW), .))
    for (k = 1; k <= cols(JW); k++) {
      nm = substr("jointW_" + strtoname(jcs[k, 2]), 1, 32)
      idx = st_addvar("double", nm); st_store(., idx, (cols(JO) >= k ? JO[1, k] : .) \ JW[., k])
    }

jwacc is e(boot_jointW), ONE COLUMN PER TEST FAMILY (att, att_it, ame, ame_it). jobs is
wdbjobs, the OBSERVED COEFFICIENT VECTOR, one column per pre-treatment HORIZON per family
with the family key in the stripe's equation slot — the input wooldid_bootjwobs(fk, jobs,
wobs) reads to COMPUTE W_obs. Pairing them by the running index k writes coefficient element
k where family k's observed statistic belongs. With one family requested k = 1, so row 0 gets
the earliest pre-treatment coefficient. advW does the identical thing with aobs = wdbaobs.

Why it matters: the help sells the file as self-sufficient ("so the draws can be recentred
from the file alone"), and the observed statistic is the only thing in it that turns the 25
replication statistics into the printed p-value. A user doing what the help says gets
p = 1.000 instead of p = 0.038 for the joint pre-trend test, and 1.000 instead of 0.423 for
the adversarial test. Nothing in the file marks the gap; the only tell is that a Wald
statistic came out negative. The displayed table is unaffected — the run's own p-values come
from the correctly-computed-and-discarded wobs.

Not reachable without boottstat: the jwacc/awacc arguments at ado:437 are
cond("`boottstat'" != "", ...), so a fixed-metric bootsave file carries no jointW_* columns.

### N2 — dots(-1) collides with the internal "not given" sentinel [low / certain]
Repro: r2_dots.do, tests E/F/G/H.
- `bootstrap(standard) reps(20) dots(-1)` -> rc 0 AND the progress display prints, i.e. the
  explicit request is silently reinterpreted as the default.
- `dots(-1)` with NO engine -> rc 0, while dots(0) and dots(7) with no engine -> r(198)
  "dots() require(s) the bootstrap() option". Exactly one value of dots() escapes the
  stray-companion screen.
Cosmetic; no number moves.

---

## NEGATIVE RESULTS (attacked hard, clean)

1. permsave() FILE CONTENTS are correct (r4_permsave_contents.do, w4/). 80 clusters,
   40 replications, esprelength(2) espostlength(2):
   - 80 onset_* columns for 80 clusters; all 80 row-0 values equal the unit's own observed
     onset in the data (CHK1_bad=0);
   - the onset multiset (16/16/16 treated + 32 never-treated) is preserved in ALL 41 rows
     (CHK2_bad=0) — a shuffle never invents or loses an onset;
   - for every coefficient column, the file's SD over reps 1..R equals e(perm_sd)
     (bit-identical hex on 5 of 6 columns, 8.5e-16 relative on the sixth), the non-missing
     count equals e(perm_n), and the _rep = 0 value equals e(perm_comparator_b) bit-for-bit
     (CHK3_bad=0 CHK4_bad=0 CHK5_bad=0).
   This is the permutation-side analogue of N1 and it is CLEAN.
2. permsave() under permstrata() keeps the cluster keying (r5_strata_both.do, w5/). Regions
   were assigned so stratum order crosses ascending cluster order — the exact geometry
   ado:7500 documents as having once produced silently wrong stratified inference. Result:
   S1_ROW0_BAD=0, S1_MULTISET_BAD=0, S1_STRATA_BAD=0 — the onset multiset is preserved within
   every region in every replication, so no draw leaked across strata and no column is
   mislabelled.
3. The rest of the permsave/bootsave name-rule surface is symmetric (r6, S4):
   permsave(sub/dir), bootsave(sub/dir), permsave(.dta), bootsave(.dta) are all r(198).
4. Guards near this surface fire informatively: permonsetcohorts with a cohort spanning
   several i levels -> named r(459) (r5/r6); permute() with cluster() not nesting i -> named
   r(459) naming the remedy (r6); permute()+bootstrap() together -> r(198) "bootstrap() may
   not be combined with permute(); each orchestrates its own replication loop. Run them in
   separate wooldid calls." (r7, S3).
