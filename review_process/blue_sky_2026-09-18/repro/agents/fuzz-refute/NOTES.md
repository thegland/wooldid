# fuzz-refute NOTES (skeptic pass over FUZZ-01..08)

Build asserted: `which wooldid` -> .../bluesky/build/wooldid.ado, version 1.4 2026-08-19.
My own generator: `gen_mypanel.do` (no repo data, no repo generator).

## Probes run
- p0_orient.do / .log    : level parsing mechanism
- pA_label_level_replay_dots.do / .log : FUZZ-01, -03, -04, -06 (first pass; -06 leg invalid)
- pB_level_deep.do / .log: level attribution (wooldid vs reghdfe/margins), FUZZ-04 isolation,
                           FUZZ-06 retry (permute(reassign) invalid -> leg invalid again)

## Established so far

### FUZZ-01 CONFIRMED (with two corrections to the claim)
- binary mode, contreatelasticitytype(dydx|dyex|eydx|eyex): all rc 0, att and se bit-identical
  to baseline (%21x +1.fa285156713e0X+000 / +1.144ad4828c292X-004, SAME=11 on both).
- Header really does read "Main Results (Full Estimation Sample): Marginal Treatment Effects
  Presented as Elasticities" over a level ATT of 1.9772 (design truth tau=2).
- The bare-`wooldid` REPLAY reproduces the same header, so the lie is stored, not transient.
- CORRECTION 1: the prober says it "posts nothing to e(contreatelasticitytype)". It DOES post
  e(wdid_cetype) (hidden), which is exactly what the replay header reads. What it does not post
  is a visibly-named return.
- CORRECTION 2: "Every other contreat() companion is refused r(198) in that same position" is
  FALSE. Refused: lattice(), contreatcontrols(), contreatwithin, makecfxplots (4).
  Accepted rc 0 in binary mode: cfxlattice() (with an explicit "is being ignored" Caution),
  latticeignoreweights, ccnoabsorb, contreatpoly(2), cfxplottypes(), makecfxplotsbysg,
  cfxplotnose, semielasticity (8). So inertness is common; what is unique to
  contreatelasticitytype() is that it is inert AND rewrites the header.
- contreatelasticitytype(garbage) is refused r(198) even with no contreat(): the domain screen
  fires without a parent screen.

### FUZZ-03 CONFIRMED in substance, mechanism PARTIALLY WRONG, and much wider than reported
- level(99.9) -> r(198) "level() can have at most two digits after the decimal point";
  level(99.99) -> same. ALSO level(95.4) -> same (not in the prober's list).
- Accepted: 95 90 99 99.5 97.5 90.25 99.95 68.27 99.8 99.7 99.75 84.13 92.5.
- MECHANISM: NOT "binary-representable". Stata renders macro `=` results at 16 significant
  digits; 99.9 -> "99.90000000000001", 99.99 -> "99.99000000000001", 95.4 -> "95.40000000000001",
  while 99.95 -> "99.95" and 68.27 -> "68.27". Whatever round-trips short is fine.
- **NEW AND WORSE**: the mangling ALSO happens at `local wdlvl = c(level)` (ado:578). With
  `set level 99.9` in force and NO level() option at all, EVERY wooldid call fails r(198)
  with the same message: att robust, att robust safety(off), att unconditionalse,
  att robust fastnose -> all 198. `set level 95` -> rc 0 again.
  Plain `regress`, `reghdfe y x1, absorb(i t)` and `margins, dydx()` all rc 0 at c(level)=99.9,
  so this is wooldid's own re-render, not a Stata-wide restriction.
  => a legal session setting (profile.do `set level 99.9`) makes wooldid entirely unusable, and
  the message blames a level() the user never typed. This is a NEW finding beyond FUZZ-03.

### FUZZ-04 CONFIRMED, and the blast radius is LARGER than reported; one sub-claim REFUTED
- estimates store m1 -> rc 0; estimates replay m1 -> r(198) "Error: command line input should
  feature all of y i t ttre; too few variables specified."
- NEW: the failed replay sets the shell's WOOLDID_LASTRC=198, so the WORKING route (bare
  `wooldid`) is then refused too: bare wooldid rc 0 before the store, rc 0 after the store,
  r(198) after the failed `estimates replay`, and STILL r(198) after `estimates restore m1`
  (GLOBAL_LASTRC=[198]). The prober reported BARE2_RC=0 after restore; that does not reproduce.
- REFUTED sub-claim: "even a trailing blank defeats it". `wooldid  ` (two trailing blanks)
  reached the EMPTY-line branch - it printed the LASTRC refusal, which lives strictly inside
  that branch (ado:40-44). So trailing blanks do NOT defeat the gate; the prober's r(198) for
  the trailing-blank case was the poisoned-global refusal, misattributed.

### FUZZ-06 NOT YET TESTED VALIDLY
- pA used permute(shuffle) without cluster() -> refused for a different reason.
- pB used permute(reassign) -> not a legal method. Legal methods: shuffle | probabilities,
  and cluster() with one variable is required. Retry in pC.

### FUZZ-07 not to be executed: ado:558-562 is `net install wooldid, from(github...) replace`
  then `exit`, before wooldid_validate. Code is the proof; running it writes to the user's PLUS
  directory over the network.

## Still to do
- pC: FUZZ-06 valid legs; level grid mapped with STRING-built levels (pB's grid was circular -
  `local L = `k'/100` is itself 16-digit-rendered, so it found 0 mismatches); trailing blank
  from a clean state.
- pD: FUZZ-05 zero-noise invariance (the discriminating test) + |dATT|/se on noisy data.
- pE: FUZZ-08 pristine-twice, then reghdfe-direct sort test.
- pF: FUZZ-02 pass-through separation (plain margins, eydx on the same fitted model) and the
  r(3598) tempvar leg in isolation.

## FINAL STATE (all probes complete; FINDINGS.md is the deliverable)

Probe inventory (all in this directory, each with its .log):
  p0_orient.do              level parsing mechanism
  pA_label_level_replay_dots.do   FUZZ-01/-03/-04/-06 first pass
  pB_level_deep.do          level attribution (wooldid vs reghdfe/margins); FUZZ-04 isolation
  pC_grid_dots_blank.do     FUZZ-06 valid legs; trailing blank from clean state; level override
  pD_level_state.do         level grid mapped (1152/9000, 144/900); session-state corruption
  pE_eswindow.do            FUZZ-05 zero-noise battery + level poison persistence
  pF_level_trace.do         set trace proof of WHERE the level error fires (ado:576 and ado:70)
  pG_sort.do                FUZZ-08 pristine-twice, reverse sort, reghdfe direct
  pH_semiel.do / pH2_semiel_controls.do / pI_final.do   FUZZ-02 + closing probes
  repro/REPRO_*.do          nine self-contained repro files, each RUN and logged

Verdicts: F01 CONFIRMED (med), F02 CONFIRMED (med), F03 CONFIRMED/mechanism PARTIAL (med, unchanged),
F04 CONFIRMED + one sub-claim REFUTED (med, raised), F05 PARTIAL (low, lowered),
F06 CONFIRMED (low), F07 CONFIRMED from source only (low), F08 PARTIAL (low).
New: NEW-1 level/c(level) brick + state corruption (med); NEW-2 silent inert companions (low).

Key traps recorded for whoever comes next:
  * `permute()` needs cluster() with EXACTLY ONE variable and permreps() >= 20; methods are
    shuffle | probabilities only. Two separate probes wasted a run on this.
  * weights go BEFORE the comma: `wooldid y i t ttre [aw=w], att ...`. `, att ... [aw=w]` is
    "option [ not allowed".
  * `foreach x in a b { cmd }` on ONE line is r(198) in Stata ("code follows on the same line
    as open brace").
  * a failed `estimates replay` sets WOOLDID_LASTRC and poisons every later bare-`wooldid`
    replay in the session; any probe that tests the bare replay must `macro drop WOOLDID_LASTRC`
    or run it before anything that can fail.
  * `local x = <expr>` in Stata renders at 16 SIGNIFICANT DIGITS, so a grid built as
    `local L = `k'/100` is already mangled and comparing it to a program's echo finds nothing.
    Build decimal test values as STRINGS.
