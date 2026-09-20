# hv-reps-caution — NOTES (live)

## Claim under verification
The reps() p-value-grid caution (ado:6473-6478) is unreachable without cluster().

## Status of the static half: CONFIRMED, and it is a REGRESSION
- ado:6459 `if "`clus'" != "" {` ... ado:6501 `}` encloses ado:6473-6478 verbatim.
  Printed 6455-6502 with line numbers; the only `}` between 6468 and 6501 are 6478 (the
  caution block's own), 6499 and 6500 (the ttrevar/treated-cluster block). So the caution
  block IS inside the cluster branch. Its indentation (2 spaces, same as the enclosing
  `if`) is why it reads as top-level.
- Call site ado:4941 passes `clus(`cluster')` straight from the user option. NO default
  substitution, so with no cluster() the local is empty and the branch is skipped.
- Provenance: commit 813c414 "Fix three disclosure defects..." (2026-09-13, RF-2).
  `git show 813c414 -- wooldidinstall/wooldid.ado` removes the block from top level
  (parent file line 6014, between the CS-20 comment block and `* companion-value
  validation`) and re-adds it inside the cluster branch. The parent is saved at
  ./prev_wooldid.ado.
- The block it defeats is CS-20 (author ruling 2026-09-07, "caution, no floor"), whose own
  comment at parent:5996-6013 says reps() "carried neither, so reps(20) printed percentile
  p-values on a .10 grid in silence" and that wooldid_pfloor "prints only when a displayed
  p is exactly 0 ... so a low-rep run whose p-values all happen to be non-zero was told
  nothing at all". That is exactly the state restored for unclustered runs.
- Commit message claims "a healthy run still prints it exactly once" — a clustered run.

## Doc-vs-code
- sthlp:535: "A rank-based p-value cannot be finer than its own grid, and wooldid prints a
  Caution naming that grid whenever it is coarser than .05". Unconditional.
- sthlp:634-636: "reps() and permreps() both print a caution when the floor they imply sits
  above .05, naming the number". Unconditional.

## Permute side: CLEAN
- wooldid_permctl ado:5899-5912 sits at program top level (preceded by the permreps<20
  refusal 5885-5892, followed by `if "`permmethod'" == "probabilities"` at 5913). No
  cluster conditional wraps it. And permute() REQUIRES cluster() anyway (ado:5946-5949),
  so there is no unclustered permute leg to strand.
- permreps() also has a hard floor of 20 (ado:5889), bounding the hazard.

## TODO
- [x] empirical: DONE, confirmed 3/3 clustered vs 0/5 unclustered
- [ ] all three else-if legs (percentile / boottstat / bootdisplay(wald))
- [ ] show the 2/reps grid is materially present unclustered
- [ ] mechanical brace/condition sweep of every `Caution` in bootctl + permctl

## RESULT (2026-09-17, complete)
CLAIM CONFIRMED empirically and statically; see FINDINGS.md.
- repro_reps_caution.do/.log : 9 legs, reps(10). caution 0/5 unclustered, 3/3 clustered.
- probe2_boundary.do/.log    : reps(39)/reps(40) boundary correct when reachable; reps(39)
                               unclustered silent; permute-without-cluster prints its caution
                               before its own r(198) refusal (HV-REPS-02).
- sweep_cautions.py/sweep_all.txt : mechanical brace sweep of all 143 Caution lines. Only
                               ado:6475-6477 are subject-mismatched with their condition.
- CONSEQUENCE NUANCE: the routine footer at ado:11784-11785 DOES name the grid unconditionally
  on the percentile and boottstat families, so the unclustered loss there is the judgment
  ("above .05, cannot reject, use 40 reps"), not the number. Under bootdisplay(wald) without
  boottstat there is NO footer branch at all and wooldid_pfloor exits, so that leg is silent.
