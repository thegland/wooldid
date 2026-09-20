# hv-reps-caution — SKEPTIC's verdicts and own findings

Build asserted in every log: `which wooldid` ->
`.../scratchpad/bluesky/build/wooldid.ado`, `*! version 1.4 2026-08-19`.
All evidence below is from MY OWN data and MY OWN do-files:
- `refute_reps_caution.do` / `.log`  (70 units x 12 years, 3 staggered cohorts + 20 never)
- `probe2_ordering.do` / `.log`      (6-state panel; all-treated panel)

---

## hv-reps-01 — CONFIRMED (medium, certain). Mechanism and severity both hold.

Own measurement (refute_reps_caution.log, counts of `Caution: with reps` per leg):

| leg | command | caution | footer floor line |
|-----|---------|---------|-------------------|
| A | `bootstrap(standard) reps(10)` | **0** | yes (`resolution 2/10`) |
| B | A + `cluster(unit)`            | 1 | yes |
| C | `bootstrap(standard) reps(18) bootdisplay(wald) jointtests esprelength(2) espostlength(2)` | **0** | **none** |
| D | C + `cluster(unit)`            | 1 | none |
| E | `bootstrap(standard) reps(10) boottstat` | **0** | yes (`both 2/(n + 1)`) |
| F | `bootstrap(standard) reps(39)` | **0** | yes |
| F2| F + `cluster(unit)`            | 1 | yes |

Every leg rc 0. Four for four silent without `cluster()`, three for three with it, on
percentile / wald / boottstat / the sthlp:538 boundary alike. The `cond()` ladder picks the
right grid at the same time: leg D's text is the WALD branch, "cannot fall below 0.053
(1/(reps + 1))" = 1/19, while leg B/F2 print "a grid of 0.200 (2/reps)".

Structural cause verified independently of the prober: `ado:6459` opens
`if "`clus'" != "" {`, `:6473-6478` is the caution, the branch closes at `:6501`;
`wooldid_bootctl`'s signature (`:6194`) takes `clus(varlist)` and the call site (`:4941`)
passes `clus(`cluster')` with no default, so "no cluster()" is an empty macro.
`git show 813c414 -- wooldidinstall/wooldid.ado` deletes the block from parent-file line
6011 (top level) and re-inserts it at 6233+ inside the branch: the regression is real.

A wording sweep (not just the identical string) finds NO other floor disclosure:
`grep -n "grid of\|cannot fall below\|smallest studentized\|smallest attainable"` returns
only 5911 (permute), 6475-6477 (the nested block) and the footer 11784/11785.

SEVERITY STANDS AT MEDIUM, with the prober's own narrowing confirmed: under percentile and
under boottstat the footer still names the *number* (leg A, leg E), so what is lost there is
the judgment ("above .05", "nothing can reject", "use at least 40"). Under
`bootdisplay(wald)` + `jointtests` without `cluster()` (leg C) nothing at all is said: the
run printed `p-value for joint test that all pre-treatment att horizon estimates are 0:
0.63..` and the complete leg contains no statement of the 1/19 = .0526 floor — the footer
line at `ado:11785` is guarded `else if "`wdbdsp'" != "wald"` and `wooldid_pfloor` exits at
`ado:10783` on that same condition. sthlp:535 and sthlp:634-636 are unconditional.

## hv-reps-02 — CONFIRMED (low, certain).

`refute_reps_caution.log`, LEG_G: `permute(shuffle) permreps(20) rseed(105) dots(0)` with no
`cluster()` prints `Caution: with permreps(20) the smallest attainable per-coefficient
permutation p-value is 0.095 (2/(permreps + 1)), above .05, ...` and then
`Error: permute() requires cluster() with a single variable ...`, rc 198.
`grep -n "permute() requires cluster"` returns exactly one site (`ado:5947`), after the
caution at `ado:5911`, so nothing earlier intercepts it. Cosmetic; correctly filed low.

---

## NEW (mine, not seeded)

### N1 — RF-2's own rule is still broken in BOTH engines: the caution precedes the caller's r(459)
`probe2_ordering.log`, panel 2 (every unit eventually treated):
- LEG_K `bootstrap(standard) reps(10) cluster(unit)` -> prints
  `Caution: with reps(10) the coefficient p-values fall on a grid of 0.200 (2/reps...` and
  THEN `Error: you have no never-treated control group and no not-yet-treated control g...`,
  rc 459.
- LEG_L `permute(shuffle) permreps(20) cluster(unit)` -> same shape, rc 459.
- LEG_K2 (same panel, no engine) -> the r(459) alone, rc 459.
813c414 moved the block "below every aborting validation", but `wooldid_bootctl` is called at
`ado:4941` and that r(459) is raised at `ado:958`, inside `wooldid_worker` (defined at
`ado:531`), which runs downstream of the bootctl call - so the move only cleared the aborts
local to the subprogram. The commit therefore surrendered the
unconditionality of the disclosure (hv-reps-01) without achieving its stated goal on the
most common abort a real user hits.

### N2 — sthlp:636 contradicts sthlp:632-635 and ado:5890 on permute()'s hard minimum
Same paragraph: bolded, "permute() requires permreps() to be at least 20"; three lines later,
"neither engine has a hard minimum, since a handful of draws is a legitimate smoke test."
`ado:5889-5893` errors r(198) "permreps() must be at least 20". Measured, LEG_H of
`refute_reps_caution.log`: `permreps(19) cluster(unit)` -> that error, rc 198. sthlp:535 gets
it right ("without permreps()'s hard floor of 20"), so the trailing sentence is the wrong one.

---

## NEGATIVE RESULT (attacked, clean)
The moved block now sits between `gstats sum` and `if r(max) < 10` inside the cluster branch,
so I tested whether executing it clobbers `r(max)` and silences the cluster-count cautions.
It does not: `probe2_ordering.log` LEG_J1 (`reps(10)`, caution prints) and LEG_J2 (`reps(60)`,
it does not) BOTH print `Caution: only 6 clusters are available for the bootstrap` and
`Caution: only 3 treated cluster(s) are available for resampling`. No r()-clobber defect.
