# SKEPTIC report -- HV-EPOST-01 / HV-EPOST-02, plus two findings of my own

All evidence below is from MY OWN do-files and MY OWN generated data; I read the
prober's script but never ran it. Build asserted in every log
(`adopath ++ .../bluesky/build`).

Files here:
- `probeA_eresidue.do` / `.log`  -- my repro of both seeded claims + my own baseline
- `probeB_plotreplay.do` / `.log` -- plot-save failure after the post; spec-loop harvest
- `probeC_saveplots.do` / `.log`  -- the same with `bootstrap(standard) reps(12)`
- `probeD_histsave.do` / `.log`   -- the same landing BEFORE the post
- `repro_saveplots_rerun.do` / `.log` -- consolidated self-contained repro for the
  two new findings

Panel used throughout: 40 counties in 8 regions x 8 years (2012-2019), policy
adopted **at the region level** (regions 1-2 in 2015, 3-4 in 2017, 5-8 never) --
an ordinary administrative shape, not a stress case. 320 rows.

---------------------------------------------------------------------------
## HV-EPOST-01 -- CONFIRMED (severity medium; realism badly UNDERSTATED)

Reproduced on my own data with a trigger the prober never used.

Good run: `e(cmd)=wooldid`, `colsof(e(b))=1`, `e(N)=320`,
main_att = `+1.8ddb543955ea6X+000`.

| trigger (each from the same restored good e())        | rc  | e(cmd) after | kb  | e(N) |
|-------------------------------------------------------|-----|--------------|-----|------|
| T1 `fe(regyear)` -- region-by-year FE on a region-level policy (MINE) | 322 | **reghdfe** | 82 | 320 |
| T2 `poisson` on an outcome containing negatives        | 411 | wooldid KEPT | 1   | 320 |
| T3 `subgroups(region)`, 4 wholly-untreated levels      | 0   | (succeeded)  | 5   | 320 |
| T4 `contreat(dose) contreatcontrols(untreatedZero)`    | 322 | **reghdfe**  | 163 | 320 |
| T5 parse typo `attt`                                   | 198 | wooldid KEPT | 1   | 320 |

Two corrections, both in the direction of *worse*:

1. **The realistic trigger is not on their list at all: rerunning a do-file.**
   `histogramcohorteffects` + `saveplots(stem)` writes its .gph at ado:12390,
   reached from ado:1944 -- long before wooldid's own `ereturn post` at ado:4051.
   Run the identical line twice without `replace` and the second run dies
   r(602) "file rrstem_hist_att_i.gph already exists" **with e(cmd)=reghdfe,
   kb=82, e(b)[1,1]=0**, and the user's previous good run (+1.b7fbc5505179aX+000 = 1.71868...) gone.
   See `repro_saveplots_rerun.do` PART C and `probeD_histsave.do`.
   Their four triggers are all "you asked for something unestimable"; this one
   is "you pressed run twice".

2. **The residue is harder to detect than they measured.** Their case (c) had
   e(N)=456 against the user's 600; on my panel every wrecked e() reports
   **e(N)=320, exactly the good run's N**. Nothing in the residue's scalars
   distinguishes it from a real wooldid result.

Consequences measured on my data (probeA SEC3, probeB B2):
- scripted harvest after `cap noi`: `e(b)[1,1]` = 0 exactly, `e(N)` = 320.
- `lincom _b[<internal colname>]` -> rc 0, full b/se/t/CI from the internal model.
- `lincom _b[main_att]` -> r(111); `predict` -> r(111). (The prober's refutation
  of the recon "silently predicts" sub-claim is correct.)
- The realistic vector, measured: a loop `foreach sp in "att" "att fe(regyear)"
  "att controls(x1)" { cap noi wooldid ...; est store m`k' }` stores the reghdfe
  model as `m2`; `estimates table m1 m2 m3, keep(main_att)` returns rc 0 and
  prints a blank middle column with no warning.
- No other residue: frame `default`, no `___wooldidfullmodel` left in
  `estimates dir`, no `_ppmlhdfe_d` (r111 on confirm).

My own baseline (probeA SEC5, not inherited from them):
- `regress price mpg`; `regress price nosuchvar` (r111) -> e(cmdline) unchanged.
- eclass ado `qui regress price weight foreign` then `error 322` -> e(cmdline)
  becomes "regress price weight foreign". That is wooldid's shape.
- the SAME ado with `_estimates hold `h', restore nullok` -> e(cmdline) back to
  "regress price mpg". The idiom works; wooldid has 0 hits for it.

Severity kept at **medium**: the failure itself is loud every time (wooldid
prints its own diagnostic before the rc), so nothing is wrong on screen; the
damage is to state the user did not ask about. `realistic` moves from "an
unusual option combination" to "any user who runs a do-file twice".

---------------------------------------------------------------------------
## HV-EPOST-02 -- PARTIAL (severity revised medium -> LOW; not independent)

What reproduces (probeA SEC4):
- post-regression failure, bare `wooldid` -> "Error: command line input should
  feature all of y i t ttre; too few variables specified." r(198)  [TRUE]
- parse failure, bare `wooldid` -> the designed staleness message  [TRUE]
- so the guard does fire only for failures that left e(cmd)=wooldid  [TRUE]

What is wrong with the claim:

1. **"Contradicts the design that comment states" -- no.** The branch is a
   *replay* branch; its job is to refuse a **misleading redisplay**. With
   e(cmd)=reghdfe there is nothing to redisplay: the branch would print a
   reghdfe model as a wooldid table, which is precisely what the gate prevents.
   The guard is not unreachable "for the class it was written for"; it is moot
   in a state where its function cannot arise. The real defect there is
   HV-EPOST-01, not the gate.

2. **"An error about a command line they never typed" -- not a misdiagnosis.**
   probeA SEC4d: with the marker dropped and e(cmd)=regress -- a virgin state in
   which no wooldid ever ran -- bare `wooldid` gives *the same* r(198)
   "too few variables specified". It is the generic invalid-call error for
   `wooldid` typed with no variables, identical in a fresh session.

3. **Their proposed remedy is wrong.** probeA SEC4e: gate on WOOLDID_LASTRC
   alone and a session where any earlier wooldid failed would answer a bare
   `wooldid` -- with `regress` results in e() the user deliberately put there --
   with "the stored results predate that failed call", and would hijack the
   generic invalid-syntax error for every later bare `wooldid`.

4. **It is a symptom, demonstrated directly.** probeA SEC4c: after the same
   post-regression r(322), `est restore good` (valid wooldid results back in
   e()) and bare `wooldid` prints "the last wooldid call failed with r(322);
   the stored results predate that failed call." The gate works the moment e()
   is not wrecked. Fix HV-EPOST-01 and HV-EPOST-02 disappears.

Severity **low**, `realistic` true but consequence is one less-informative
error message at a moment when an error is already correct.

---------------------------------------------------------------------------
## NF-SKEP-01 (NEW, mine) -- `saveplots()` never got the pre-flight file-name
## check the author added to bootsave()/permsave(), so rerunning a do-file loses a
## COMPLETED run -- and under bootstrap() the whole replication loop

`permsave()` and `bootsave()` check their file name for an existing file
**before** any work runs (ado:5866-5877, ado:6276-6298); the comment at
ado:6276-6278 says why: "Until 2026-09-05 only wooldid_bootsave checked them,
after the loop, so a name clash cost the whole run".

`saveplots()` has no such check anywhere (`grep -n "confirm file"` in the ado
returns only the four permsave/bootsave sites). Every `graph save` /`save` it
performs is bare, e.g. ado:13325 `qui graph save "`saveplots'_`graph'.gph",
`replace'`, ado:13339, ado:12390, ado:2494.

Measured (`repro_saveplots_rerun.do`):
- PART B, plain analytic run with `makeplots saveplots(rrstem)`: run 1 rc 0;
  run 2, identical line, **r(602) after the full result table has printed**.
- `probeC_saveplots.do` C1: with `bootstrap(standard) reps(12)`, run 2 completes
  all 12 replications and the rewiring (`e(boot_method)=standard`,
  `e(boot_reps)=12`, kb=6) and *then* dies r(602) at the plot save. The entire
  replication loop is lost to a file name -- the exact cost the bootsave()
  pre-flight was added to prevent.
- PART C, `histogramcohorteffects saveplots(rrstem)`: run 2 dies r(602) at
  ado:12390, *before* the post, so e() is left holding reghdfe (see above).
- PART D: both reruns succeed with `replace`. Trigger is only the missing option.

**This is NOT a claim that plot failures should not error.** ado:277-280 states
the design plainly: "A plot failure is an error of this call (its rc propagates,
the caller's frame is restored and the replay marker is set), as it is on a plain
run", and known_bydesign 6 rules the WOOLDID_LASTRC machinery as NFRPL-01/02
design. The claim is narrower: the file-exists case is the one trigger that is
trivially pre-flightable, the author already wrote that pre-flight twice for
exactly this reason, and saveplots() -- the option most likely to hit it, because
it writes on every ordinary rerun -- was left out. Under bootstrap()/permute()
the cost is the whole replication loop, which is verbatim what ado:6276-6278
says the 2026-09-05 change was for.

On the help (stated carefully, because the scope matters): sthlp:696-697 says
"If that write fails, typically because a file of the name already exists and
{opt replace} was not given, **the run is not lost**: a caution says the .dta is
missing, and the figures, the .gph files, the estimates and every stored result
are unaffected." Read strictly, "that write" is the `_cfx.dta` write named in the
preceding clause, and the sentence is internally coherent (it lists the .gph
files among the *unaffected* things). I did NOT measure that soft-fail -- by code
reading it is the `cap noi wooldid_cfxplots , ... dataonly` at shell ado:517-524.
So the accurate statement is: the help documents a soft failure for one
saveplots() write and is **silent** that the event-study .gph/.dta and the
histogram .gph writes hard-fail and lose the run. My inference, not the help's
words, is that sthlp:700 ("for those two the name is checked before the
replications run and an existing file is a hard error (r(602)) rather than a
silent skip", "those two" = permsave and bootsave) invites a reader to place
saveplots() in the soft class.

Severity **medium**, anchored on C1: the bootstrap case loses a completed
12-replication run. The plain analytic case alone would be low.
Realistic **true** (rerunning a do-file without `replace`).

---------------------------------------------------------------------------
## NF-SKEP-02 (NEW, mine) -- the replay guard misdescribes COMPLETE, VALID
## results as predating the failure

The mirror image of HV-EPOST-02, and the direction that actually misleads.

When the r(602) above lands after `ereturn post`, e() holds the complete result
of that very call -- `e(cmd)=wooldid`, all six columns, `e(wdidmainresults)`
present, `e(N)=320`, and its table already printed on screen. Bare `wooldid`
then answers (`repro_saveplots_rerun.do` PART B, `probeB_plotreplay.do`):

    Error: the last wooldid call failed with r(602); the stored results
    predate that failed call. Refit wooldid before replaying.

They do not predate it: they *are* it. Under `bootstrap(standard) reps(12)`
(probeC C1) the same message hides a fully rewired bootstrap result the user
just paid 12 replications for. Control: with `replace` the run succeeds and
bare `wooldid` replays at rc 0 (probeB B1_THIRD_REPLAY_RC=0).

**The refusal itself is design, and I am not disputing it**: ado:277-280 says a
plot failure "is an error of this call (its rc propagates, the caller's frame is
restored and the replay marker is set)", and ado:10731 shows the author knows a
mid-display abort then refuses the replay. The finding is the message TEXT: in
this class the assertion "the stored results predate that failed call" is
measurably false.

Measured cost, not just wording (probeC C1): in bootstrap run 2 the r(602) cut
in after `wooldid_display, part(1)` but before `part(boot)` -- the whole
bootstrap interpretation block that run 1 printed (the p-resolution note, the
`e(V)` / `e(df_r)` / esttab-coefplot-honestdid caveats) never appeared. So the
user paid 12 replications, has the fully rewired inference sitting in e(), was
never shown half of its documentation, and the one command that would print it
is refused with a false statement about what is stored.

Severity **low**, realistic **true**.

---------------------------------------------------------------------------
## Negative results / things I attacked that were clean
- `poisson` on an outcome containing negative values is caught at validation
  (r411) with e() intact -- not an e()-residue path.
- `subgroups(region)` with four wholly-untreated region levels ran to rc 0.
- A parse-stage failure never damages e() (T5), matching their control.
- `predict` after a wrecked e() errors r(111); it does not silently predict.
- After the wreck, no stray `___wooldidfullmodel`, no `_ppmlhdfe_d`, frame back
  to `default` -- the shell epilogue does its registered job; e() is the single
  unhandled piece of state, as claimed.
- `permsave()`/`bootsave()` name clashes ARE pre-flighted, so the analogous
  "lose the whole run to a file name" bug does not exist there. Read-only check:
  ado:168-176 (permutation) and 367-398 (bootstrap) do bracket their rep loops
  with `est sto`/`est restore`, so the rep engines do not leak e() either.
- Sub-threshold, NOT filed: two identical `wooldid` calls in one session on
  unchanged data differed in the last bit of main_att
  (+1.b74ac77ebbd5fX+000 vs +1.b74ac77ebbd5eX+000, probeB B1). 1 ulp, far below
  the project's rtol 1e-10 indifference threshold; most likely a sort-order
  dependent summation in the collapse. Noted, not claimed.
