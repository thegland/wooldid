# fp-estimands-refute — SKEPTIC notes

Assignment: refute (or confirm) four findings from the FIRST-PRINCIPLES ORACLE prober.
Rule: build my own repro from the claim's terms, on my own data.

## Plan
- A. FP-SEMIELAS-FE. Own DGP (unit-level heterogeneous intercepts) x {default, customvce(robust)}
     x {level, semielasticity}. Plus a wooldid-free mechanism demo
     (reghdfe+margins,eydx vs reg i.i i.t + margins,eydx). Homogeneous-level control arm.
     PRIOR SEARCH DONE: grep of review_process/ + wooldid_ado_review.md for eydx/semielas/xbd
     finds ONLY the poisson xb-scale remap (BS-6). The reghdfe-absorbed-FE issue is NOT
     on record anywhere -> if it reproduces it is novel.
- B. FP-CTZERO-ENGINE. Their DGP has `ui` generated AFTER expand => per-ROW noise (SD~5.8),
     so "low noise" is false. Rebuild clean (unit-level ui, sigma 0.3). Compare BOTH arms to
     the DGP truth 3.0*E[dose|post]. Prober's DEFAULT 4.358 is ~= truth; EXPLICIT -5.07 is the
     deviant => the framing may be inverted.
- C. FP-NOPRE-LATEENTRY. One tiny panel; confirm per-unit mint drop + messaging; grep help.
- D. FP-CWM-DOC. sthlp:220 (att phang) DOES name cwm; sthlp:1107 covers it generally; the
     three cell-weighted phangs do not. Doc-only.
- E. NEW: does wooldid_ccscreen (ado:11818) ever fire on the DEFAULT reghdfe engine? It scans
     e(b) colnames for "o."; reghdfe may strip omitted terms.

## Status
- (in progress)

## Status: COMPLETE (2026-09-17)
All four seeded findings adjudicated; FINDINGS.md has the full write-up.
- FP-SEMIELAS-FE  -> CONFIRMED high (+ 2 extensions: ey* contreat operators measured;
                     r(3598) Mata traceback when reghdfe's _cons+Xb goes negative)
- FP-CTZERO-ENGINE-> PARTIAL, culprit INVERTED, severity raised to high; mechanism is the
                     DF-1 ringer at ado:1233 surviving the untreatedZero branch at ado:1313;
                     new: poisson hits the DEFAULT engine the same way
- FP-NOPRE-LATEENTRY -> CONFIRMED low
- FP-CWM-DOC      -> CONFIRMED low (nuance: sthlp:220's att phang DOES name cwm)
Negatives: ccscreen is NOT blind on reghdfe (reghdfe keeps o. names); basic-controls engine
parity is bit-exact; dydx/level estimands agree to 1-2 ulp everywhere.
