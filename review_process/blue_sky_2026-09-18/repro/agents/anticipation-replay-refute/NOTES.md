# SKEPTIC pass on anticipation/replay/level findings (AR-1..AR-4)

## Done so far (no wooldid run yet)
- s1_render.do/.log — how Stata renders levels into macros.
  * `syntax [, Level(real 0)]` with 99.9 -> local = "99.90000000000001"; 99.99 -> "99.98999999999999";
    98.9 -> "98.90000000000001"; 79.9 -> "79.90000000000001"; but 99.91 -> "99.91" (clean).
  * `global g = c(level)` after `set level 99.9` -> "99.90000000000001"; after `set level 99.91`
    -> "99.91000000000001" (NOT clean!).
  * `set level 99.90000000000001` -> r(198), message "level() can have at most two digits after
    the decimal point".  **That message is STATA's, not wooldid's** (grep of wooldid.ado for the
    string returns nothing).
- s2_setlevel.do/.log, s2_map.csv — confirms the message source and a sampled map.
- s3_maps.do/.log, s3_full.csv — FULL enumeration of all 9000 legal 2-dp levels:
  * MAP A (the `level()` option path: parse literal -> macro -> `set level`): **1152/9000 fail**,
    first failure L=64.01. This matches AR-4's claim exactly (1152, all >= 64.01).
  * MAP B (the SESSION path: `set level L` -> `global = c(level)` -> restore): a STRICT SUPERSET.
    Over the 7848 rows testable in s3 it fails 1768; A-only = 0, B-only = 616. First map-B
    failure L=**32.22**, far below 64.01. e.g. L=99.91 and L=99.94 PASS as an option but POISON
    the restore when they are the session level.
  * HEX: `set level 99.91` stores +1.8fa3d70a3d70bX+006 while the parsed literal 99.91 is
    +1.8fa3d70a3d70aX+006 — one ulp apart. That is why the two maps differ.

## Next
- AR-4 wooldid runs (a)-(e) incl. the sticky-global follow-up.
- AR-3 + AR-2.
- AR-1.
- Novel: anticipation(k) vs shifted-ttre differential oracle.

## FINAL STATE (all runs complete)
Verdicts: AR-4 CONFIRMED/high (scope extended); AR-3 CONFIRMED/low; AR-2 CONFIRMED/low
(divergence is two lines, not one); AR-1 PARTIAL/low (behaviour confirmed, "omission not policy"
mechanism refuted - e() carries no flag to reconstruct the ame/jointtests cautions).
New: stale $WOOLDID_LEVEL0 masks the true rc and skips the shell's rc-gated cleanup epilogue on
every later parse-stage failure (medium).
Clean negatives: anticipation() == manual ttre shift bit-for-bit on 4 estimands / ES /
esfixedbaseperiod / coarse i; anticipation(-2) refused informatively; level() propagation exact
to 1.5e-15 on 260 gteffects rows; CS-2 holds under bootstrap(); replay with no data in memory.
See FINDINGS.md.

## CORRECTION (s11_mapb.do, typed-literal enumeration)
The earlier "2920 / 32%" estimate for map B was an inference, and it was WRONG. Measured exactly:
map A (option path) 1152/9000 (12.8%, first 64.01); map B (session path) 1768/9000 (19.6%,
first 32.22); map A is a strict SUBSET of map B (both=1152, A-only=0, B-only=616). `set level'
itself refuses none of the 9000. FINDINGS.md updated.
