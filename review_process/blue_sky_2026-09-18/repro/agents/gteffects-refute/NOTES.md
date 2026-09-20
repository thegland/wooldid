# SKEPTIC pass on GTE-1..GTE-5 — COMPLETE

Verdicts and evidence are in FINDINGS.md. Summary:

| id | verdict | severity |
|---|---|---|
| GTE-3 late-entry silent drop | CONFIRMED (two legs of the argument corrected) | medium |
| GTE-1 customweightmultiplier | PARTIAL — code matches its own paragraph; three OTHER paragraphs lack the qualifier | low |
| GTE-2 fastnose posts b_es but not V_es | CONFIRMED as a measurement; it is a negative result, not a defect | none |
| GTE-4 w_att_it populated when unrequested | CONFIRMED | low |
| GTE-5 timetrends se column not missing | CONFIRMED; the prober stopped one step short (see NEW-2) | low |

NEW (mine):
- NEW-1 [high] esfixedbaseperiod silently deletes any treated unit lacking an observation at
  exactly esrelativeto(). One missing outcome year moves the ATT from 2.00 to 1.02 on a DGP
  whose truth is 2.0. Repro: repro_esf_refperiod_drop.do
- NEW-2 [high] gteffects reports numerically-zero cell variances as real SEs (se ~1e-15,
  p down to 8e-258, zero-width CIs) although the help promises they are reported missing.
  Guard at ado:13496 is exact-zero only. Repro: repro_gtse_degenerate.do

Files: r01_low4 (GTE-1/2/4/5), r02_drops (GTE-3 + NEW-1), r03_degen (fastnose + NEW-2 discovery),
r04_statepanel (NEW-2 on a state panel), r05_map (NEW-2 reach), repro_esf_refperiod_drop,
repro_gtse_degenerate. Each has its .log next to it.

METHOD TRAP recorded for the next agent: `cap matrix Z = e(name)` returns _rc==0 even when
e(name) is not a matrix — Stata evaluates e(name) as a missing SCALAR and builds a 1x1 matrix.
Always check rowsof() or read `ereturn list`. This trap gave me the wrong answer on GTE-2 first
time round.
