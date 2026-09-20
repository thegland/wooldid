# SKEPTIC notes -- HV-EPOST-01 / HV-EPOST-02

Workdir: .../probes/hv-failed-ereturn-refute
Build asserted in probeA log (adopath ++ .../bluesky/build).

## Probe A (probeA_eresidue.do / .log) -- my own 40-county x 8-year panel,
## policy adopted AT THE REGION LEVEL (regions 1-2 in 2015, 3-4 in 2017, 5-8 never).

Good run: e(cmd)=wooldid, e(cmdline)="wooldid y unit year onset, att",
colsof(e(b))=1, e(N)=320, main_att=+1.8ddb543955ea6X+000 (=1.5546...).

Trigger battery (each from the same restored good e()):
| tag | spec                                        | rc  | e(cmd) after | kb  | e(N) |
|-----|---------------------------------------------|-----|--------------|-----|------|
| T1  | fe(regyear)  <- MY OWN, not on their list   | 322 | reghdfe      | 82  | 320  |
| T2  | poisson on an outcome with negatives        | 411 | wooldid KEPT | 1   | 320  |
| T3  | subgroups(region), 4 wholly-untreated levels| 0   | (succeeded)  | 5   | 320  |
| T4  | contreat(dose) contreatcontrols(untreatedZero)| 322| reghdfe      | 163 | 320  |
| T5  | parse typo `attt`                           | 198 | wooldid KEPT | 1   | 320  |

=> HV-EPOST-01's core mechanism REPRODUCES on my data with a trigger they never
   used (T1). T1 is MORE realistic than any of their four: a state/region-by-year
   FE on a panel whose policy varies at the region level.

WORSE THAN THEY SAID: after T1, e(N)=320 -- IDENTICAL to the good run's e(N).
Their case (c) had e(N)=456 vs 600, which an N sanity-check would catch. Here
nothing about the residue's scalars distinguishes it from the real thing.
e(b)[1,1] = exactly 0 (an omitted base cell). Harvest reads b=0, N=320.

Consequences measured: lincom on an internal colname -> rc 0 (silent, wrong
numbers); lincom _b[main_att] -> r(111); predict -> r(111) (their refutation of
the recon "silently predicts" sub-claim is correct).

## HV-EPOST-02 (SEC4 of probe A)
- post-regression failure, bare `wooldid` -> "Error: command line input should
  feature all of y i t ttre; too few variables specified." r(198)
- parse failure, bare `wooldid` -> designed message, r(198)
- after `est restore good`, bare `wooldid` -> designed message, r(198)
- SEC4d: VIRGIN state (marker dropped, e(cmd)=regress), bare `wooldid` ->
  THE SAME "too few variables" r(198). So the post-regression message is not a
  misdiagnosis; it is the generic invalid-call error, identical to what a user
  gets typing `wooldid` alone in a fresh session.
- SEC4e: with e(cmd)=regress and the marker set, the prober's proposed remedy
  (gate on WOOLDID_LASTRC alone) would print "the stored results predate that
  failed call" about a regress model -> the remedy is wrong.

## Baseline (SEC5, my own)
regress price mpg; then
 (a) `regress price nosuchvar` r111 -> e(cmdline) still "regress price mpg"
 (b) eclass ado: regress weight foreign then error 322 -> e(cmdline) becomes
     "regress price weight foreign"  (wooldid's shape)
 (c) same ado with `_estimates hold `h', restore nullok` -> "regress price mpg"
=> the idiom works and wooldid does not use it (0 hits for _estimates hold).

## Still to run
- Probe B: plot-save failure AFTER a completed estimation (saveplots without
  replace on a rerun) -> is a complete, valid e() then refused for replay with
  "the stored results predate that failed call"?
- Probe B: the realistic loop-over-specs harvest with est store.

## Probe B (probeB_plotreplay.do/.log) -- DONE
B1: `makeplots saveplots(pbstem)` run twice -> run 2 prints the FULL result table,
then r(602) "file pbstem_ES_att.gph already exists". e(cmd)=wooldid, kb=6, N=320,
e(wdidmainresults) present -> COMPLETE valid result. Bare `wooldid` ->
"the last wooldid call failed with r(602); the stored results predate that failed
call."  FALSE: they ARE that call's results. With `replace`: rc 0, replay rc 0.
B2: spec loop with `cap noi wooldid; est store m`k'' -> m2 is the reghdfe model
(kb=82, b=0); `estimates table m1 m2 m3, keep(main_att)` rc 0, blank middle column.
Sub-threshold note: B1 run1 vs run2 main_att differ by 1 ulp (...5f vs ...5e).

## Probe C (probeC_saveplots.do/.log) -- DONE
C1: bootstrap(standard) reps(12) + makeplots + saveplots, run twice -> run 2
completes all 12 reps AND the rewiring (e(boot_method)=standard, e(boot_reps)=12,
kb=6) and then dies r(602) at the plot save. Whole loop lost to a file name.
C2 mis-specified (makecfxplots needs lattice/cfxlattice) -> superseded by probe D.

## Probe D (probeD_histsave.do/.log) -- DONE, the headline
`histogramcohorteffects saveplots(pestem)` run twice -> run 2 r(602) at ado:12390
(reached from ado:1944, BEFORE the ereturn post at 4051): e(cmd)=reghdfe, kb=82,
e(N)=320, e(b)[1,1]=0, previous good wooldid results destroyed. So HV-EPOST-01's
realistic trigger is "run the do-file twice", not an exotic option combination.

## Consolidated self-contained repro: repro_saveplots_rerun.do/.log
PARTS A-D. Everything above in one file.

## VERDICTS
HV-EPOST-01 CONFIRMED, severity medium (unchanged), realism badly understated.
HV-EPOST-02 PARTIAL, severity medium -> LOW, not independent of 01, and the
            proposed remedy (LASTRC-only gate) is wrong.
NF-SKEP-01  NEW, medium: saveplots() has no pre-flight existence check; a do-file
            rerun loses a completed run (and a whole bootstrap loop); sthlp:696-697
            says "the run is not lost".
NF-SKEP-02  NEW, low: the replay guard calls complete, valid results of the failing
            call "the stored results predate that failed call".
