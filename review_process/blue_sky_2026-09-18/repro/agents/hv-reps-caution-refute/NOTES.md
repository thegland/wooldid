# hv-reps-caution REFUTE — notes

Role: skeptic. Refute or confirm hv-reps-01 (reps() caution unreachable without cluster())
and hv-reps-02 (permute caution prints before an aborting refusal).

## Code reading done BEFORE any run (own reading, not the prober's)
- build/wooldid.ado:6459 `if "`clus'" != "" {` opens; :6466 two-cluster r459; :6473-6478 the
  reps() caution block; the branch closes at :6501. So the caution IS nested. Indentation of
  the `if `reps' > 0 {` is 2 spaces, same as the enclosing if -> reads top-level.
- bootctl signature :6194 has `clus(varlist)`; call site :4941 passes `clus(`cluster')`
  verbatim. No default substitution, so "no cluster()" == empty clus.
- Wording sweep for ANY differently-worded floor caution: only 5911 (permute), 6475/6476/6477
  (the nested block) and the FOOTER at 11784/11785. 11785 is guarded
  `else if "`wdbdsp'" != "wald"`, so wald-without-boottstat gets no footer line either.
- git show 813c414 -- wooldidinstall/wooldid.ado: the block is deleted at parent line 6011
  (top level) and re-inserted at 6233+ inside the cluster branch. Regression claim is
  structurally real.
- permctl: caution at :5911, `wordcount("`clus'") != 1` refusal at :5946-5949. Only ONE such
  refusal in the file (grep "permute() requires cluster"), so nothing earlier catches it.

## Own empirical run
refute_reps_caution.do — my own 70-unit x 12-year panel, 3 staggered cohorts + never-treated.
Legs A/B (percentile reps(10) -/+ cluster), C/D (wald+jointtests reps(18) -/+ cluster),
E (boottstat reps(10) no cluster), F/F2 (reps(39) -/+ cluster), G (permute no cluster),
H (permreps(19)), I (permute + cluster, healthy path).
STATUS: launching.

## RESULTS (all measured)
- hv-reps-01 CONFIRMED medium/certain. A/C/E/F silent (no cluster), B/D/F2 print.
- hv-reps-02 CONFIRMED low/certain (LEG_G: caution then r198).
- N1 (mine): caution still precedes the CALLER's r(459) in both engines (probe2 LEG_K/LEG_L).
- N2 (mine): sthlp:636 "neither engine has a hard minimum" contradicts sthlp:632-635 and
  ado:5890 (LEG_H: permreps(19) -> r198).
- NEGATIVE: no r(max) clobber from the moved block (probe2 LEG_J1/J2).
STATUS: DONE. FINDINGS.md written.
