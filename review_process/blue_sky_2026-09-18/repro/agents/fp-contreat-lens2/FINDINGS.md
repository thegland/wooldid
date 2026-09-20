# fp-contreat-B1 — LENS 2 (realism & impact) adjudication

**Verdict: PARTIAL. Severity downgraded high -> medium.** The engine split is real and I
do not dispute it; I add an independent confirmation on a spelling nobody tried
(cluster(cty), FINER than i, where the documented containment predicate cannot fire).
But the *impact framing* is wrong on the half of the table that carries most of the
claim's weight, and the remedy the claim proposes would, on that half, replace a loud
blank with a silent wrong number.

Why medium and not high, under this lens (a silent wrong number far outranks a loud
failure): the failure here is LOUD and clean end to end — three cautions print, e(V) is
ABSENT rather than zero (r(111)), test/lincom refuse with r(321), boottstat blanks every
column rather than dividing by zero. No wrong number reaches the user. What is lost is
valid att inference on a documented headline feature across 7 of 8 lattice spellings on a
realistic administrative panel, and the message the user is handed misdiagnoses the cause
— brief tier 4 (crash/loud failure on realistic input) plus tier 5 (documented-vs-actual),
not tier 1. The prober's "high" rested on "loses EVERY reported standard error"; half of
those standard errors were not valid to begin with.

## 1. What is confirmed (I did not re-litigate; both prior agents measured it)
On contreat()+lattice() panels the wooldid default (reghdfe) path returns max|e(V)| = 0
while plain regress on the identical fit returns a finite robust covariance. Prober: 4+
spellings, 3 fixtures. Skeptic: 6/6 on an independent fixture with identical rss and
6-decimal-identical margins-ame. My own default-path legs: 12/12 with lattice(all 3).
Accepted as fact.

## 2. My own measurement (L1_impact.do, 12 draws, own fixture, own DGP truth)
Three arms on the SAME draw, all through customvce(robust) — i.e. exactly the covariance
B1 proposes wooldid should adopt.

| row | arm              | mean    | MC SD   | mean SE | SE/SD | cover95 |
|-----|------------------|---------|---------|---------|-------|---------|
| att | A lattice(all 3) | 0.85649 | 0.12465 | 0.16695 | 1.339 | 12/12   |
| att | B lattice(reg3)  | 0.85063 | 0.12468 | 0.16565 | 1.329 | 12/12   |
| att | N no lattice     | 0.85362 | 0.12396 | 0.15974 | 1.289 | 12/12   |
| ame | A lattice(all 3) | 0.75096 | 1.18272 | 0.32450 | 0.274 |  5/12   |
| ame | B lattice(reg3)  | 0.63188 | 0.55366 | 0.15573 | 0.281 |  5/12   |
| ame | N no lattice     | 0.40053 | 0.13291 | 0.11370 | 0.855 | 11/12   |

true att (per draw) mean 0.89852 sd 0.03454; true ame = 0.400 exactly.
Arm B's lattice is a substantive grouping INDEPENDENT of dose with wide within-bin dose
variation — the friendliest, most natural applied use of lattice(). It behaves like the
dose-quantile lattice, so this is not a "misuse only" surface.

### The two rows are not the same casualty
- **att: a genuine, valid loss.** SE/SD ~1.33 — with 12 draws the SD estimate carries
  ~20% relative noise, so 1.33 is consistent with nominal, not demonstrably conservative;
  what is solid is 12/12 coverage in all three arms. The att point estimate shows no
  established bias either (0.856 vs true 0.899 = 1.2 MCSE). And the lattice does not
  degrade the att's sampling SD at all (0.1247 / 0.1247 / 0.1240 across the three arms).
  The default path throws away an SE that works. THIS is the user cost.
- **ame: not a loss — the blank is protective, though only by accident.** Under EITHER
  lattice the recovered SE is 3.6x too small — and 0.27 sits ~3.5 SDs of the 12-draw noise
  away from 1.0, so unlike the att ratio this one is real and the nominal 95% interval covers 5/12 = 42%. Without a lattice the
  same engine is fine (0.855, 11/12). Adopting the regress covariance wholesale would
  print a 3.6x-understated ame SE on every such run: a silent wrong number replacing a
  loud blank. Under this lens that is a worse outcome than the status quo.
  (Independently consistent with the skeptic's own B2 MC on the same design: reported se
  0.341 against realized SD 2.519, i.e. 7.4x too small. Neither prior agent connected B2
  to B1's remedy.)
- No bias established on the ame (arm B: 0.632 vs 0.400, MCSE 0.160 -> 1.45 sigma). The
  defect is variance, not location — matching the skeptic's B2 amendment.

## 3. Incidence and intermittency (my own, novel)
zero-V incidence on the default path across 12 draws of the SAME specification
(verified on L1mc.dta: defmissB is 8 ones, 4 zeros, 0 missings — the other 4 draws
returned rc 0 WITH standard errors, they did not error out):
    lattice(all 3)  12/12   deterministic
    lattice(reg3)    8/12   **INTERMITTENT**
The same command on the same DGP returns standard errors on 4 samples and none on 8.
This directly contradicts what wooldid PRINTS: "This is a property of the specification,
not a numerical accident" (ado:1567). It is not a property of the specification — it is
a property of the draw and of the engine.

## 3b. Incidence across the specs a user would actually write (L2_side.log S5, one draw)
N=1100, i=st (12 states), dose at county, 110 counties.

| spec | att se | ame se | df_r |
|------|--------|--------|------|
| lattice(all 3) robust        | MISSING | MISSING | 846 |
| lattice(all 3) cluster(cty)  | MISSING | MISSING | 109 |
| lattice(all 2) robust        | MISSING | MISSING | 913 |
| lattice(all 2) cluster(cty)  | .1391112 | .2217282 | 109 |
| lattice(reg3) robust         | MISSING | MISSING | 828 |
| lattice(reg3) cluster(cty)   | MISSING | MISSING | 109 |
| lattice(all 5) robust        | MISSING | MISSING | 779 |
| robust        (no lattice)   | .1622831 | .1167613 | 994 |
| cluster(cty)  (no lattice)   | .1301635 | .1164210 | 109 |

Seven of eight lattice spellings lose the whole table's inference on this single draw;
neither no-lattice spelling does.
**`cluster(cty)` is FINER than `i` = st.** The documented cell-containment geometry
(sthlp:1076, known_bydesign §1/§5) requires cluster() AT OR ABOVE i, so its predicate
cannot fire here — and the zero fires anyway, on two of the three clustered lattice
spellings. This is my own independent confirmation of B1's central negative claim, on a
spelling neither prior agent tried, and it is exactly the case known_bydesign §5 names as
"not documented and IS a finding".

## 3c. No silent-number leak anywhere I could find (L2_side.log S2/S3/S4)
- `mat V = e(V)` -> r(111): e(V) is genuinely ABSENT, not a zero matrix. So esttab,
  coefplot and honestdid cannot silently read se=0.
- `test` / `lincom` -> r(321). Loud.
- `bootstrap(standard) reps(30) boottstat` on the zero-V spec: RC=0, and the printed
  table carries `att 1.0594861 . . . .` and `ame 1.4315033 . . . .` — every inference
  column blank. boottstat studentizes by the analytic SE, which does not exist here, and
  it correctly refuses to invent one rather than dividing by zero. No wrong number.
  Cost note: the user who follows wooldid's own advice ("switch to bootstrap SEs") and
  types boottstat pays for 30 replications and gets a table of blanks. Plain
  bootstrap(standard) without boottstat DOES deliver — see §3e.

## 3d. What the DEFAULT path PRINTS when it does NOT zero out (L5_printed.do, 14 draws)
Same spec, lattice(reg3) robust, default (reghdfe) path only. Fired 9/14; on the other 5
it returned a full table.

    ame across all 14 draws : mean -0.19954,  realized MC SD 1.74093   (truth 0.400)
    att across all 14 draws : mean  0.89392,  realized MC SD 0.12127
    ame SE the DEFAULT PATH PRINTED on the 5 non-zeroing draws: mean 0.1132
    att SE the DEFAULT PATH PRINTED on the 5 non-zeroing draws: mean 0.1645
    cover95 of the PRINTED default-path ame interval: 2/5

CALIBRATION (important): 1.741 is the SD over ALL 14 draws and is driven by draw 9's
-6.109, which is a FIRING draw the engine never printed an SE for. The honest comparison
is within the 5 non-firing draws only — ame values 0.639, 0.254, 0.124, 0.809, 0.346,
SD ~0.28 against a printed SE of ~0.113, i.e. about **2.5x too small**, with the nominal
95% interval covering **2 of 5**. Still damning, and 2/5 coverage is the direct statement.
Observed but NOT established: the draws the engine zeroed out include the wildest ame
values (-6.109, -0.465, -0.193, -0.128); whether the engine failure is correlated with the
fragility that makes the SE wrong is not something 14 draws can settle, and I build
nothing on it.
The att SE it prints on the same draws (0.1645 against a realized SD of 0.1213) is fine.
This is the decisive impact fact for adjudicating B1: the all-zero covariance is the ONLY
thing standing between the user and a badly wrong ame standard error on this surface. It
is not protecting them on purpose — it is a numerical failure that happens to land on a
number that was wrong anyway. Adopting the regress covariance, as B1 proposes, would print
that wrong number on 14 samples out of 14 instead of 5.
(This is fp-contreat-B2's mechanism. B2 was measured on the customvce path; measuring it
on the shipping DEFAULT path is new here.)

## 3e. The remedy wooldid's own caution names DOES work (L5_printed.log part b)
`bootstrap(standard) reps(25)` on the zero-V spec: rc 0, honest numbers, and a caution
that the bootstrap covariance went to e(V_boot) rather than e(V) (rc 503).
    att  .92697543  se .13907921  t 6.67  p 0.000
    ame  4.6055154  se 5.5481107  t 0.83  p 0.240   CI includes zero
The bootstrap SE on the ame is large and the p-value correctly says nothing is known —
exactly the right answer on a fragile lattice ame. So the user is not stranded: the
remedy the caution names delivers correct inference, at the cost of runtime. That
materially reduces what they lose and is part of why I downgrade severity.

## 4. Would the user notice? YES — and that caps severity, but the message misleads
Verified in my own run (L2_side.log): the default zero-V path prints THREE cautions and
posts no e(V) at all:
  - "estimation returned an all-zero variance-covariance matrix. The issue is likely
     related to continuous treatment effect estimation coupled with use of analytic
     clustered SEs, with clusters occuring at a high level of aggregation..."
    -> FALSE here. There is no cluster() token anywhere; the run is bare `robust`.
  - "...This is a property of the specification, not a numerical accident: see the
     cell-containment caution above..." -> the referenced cell-containment caution never
     printed, and 8/12 intermittency falsifies "not a numerical accident".
  - "standard errors were not identified for 2 of 2 reported estimates, so no
     variance-covariance matrix e(V) was posted."
e(V) is genuinely ABSENT (`mat V = e(V)` -> r(111); test/lincom -> r(321)). So there is
NO silent-wrong-number leak downstream: no t=inf, no p=0. The failure is loud and clean.
The three remedies offered (impvalues dispersion, finer-than-cell clustering, conditional
robust SEs) never name `customvce(robust)`, which restores the valid att SE in one second.
(The misdirected caution text is fp-contreat-C, confirmed separately; cited here only as
the "what are they told" component of impact.)

## 5. What a real applied economist loses
A state-level `i` with county-level dose is not an exotic shape — it is close to FORCED,
because a unit-level `i` gives singleton treated cells and no within-cell dose variation.
lattice() is a headline documented feature with 11 latticetype spellings and its own
worked example (sthlp example 8). So: realistic, yes.
What they lose is the att row's standard errors, t, p and CI — the headline number of the
table, whose inference is legitimate — plus the event-study att horizons, plus e(V).
What they do NOT lose is anything valid on the ame side.
They then read a caution that blames a cluster() they did not specify, calls the outcome
a property of their specification when the same command works on other samples, and sends
them to a bootstrap rather than to the one-word option that recovers the att SE.

## 6. Remedy — narrower and safer than the one B1 proposes
B1 says "a wooldid-side remedy is cheap because the vzero screen already detects the
condition and wooldid already knows how to drive regress". Do NOT drive regress wholesale.
1. (cheap, no numbers move) Fix the vzero caution: it must not blame clustering when no
   cluster() is present, must not assert "not a numerical accident" (8/12 intermittency),
   and should name customvce(robust) as the att-side route with an explicit warning that
   the ame SE it returns is not trustworthy under a lattice.
2. (the real fix) Reinstate the att-FAMILY rescue only — the AMESE R4 refit removed
   2026-09-08 (ado:1546-1565) — and keep the ame side blank with its caution. My numbers
   say the premise recorded for that removal ("cell-level containment has already made the
   ame inference invalid") was written for the containment geometry and does not describe
   this surface at all — there is no cluster() here, so there is no containment — while
   its conclusion about the ame happens to hold here for a different reason (lattice
   fragility, 42% coverage). The rescue would have bought a valid att SE on this surface,
   which the ruling's cost-benefit never credited it with. The removal's cost objection (a solve over hundreds of collinear dummy
   columns) still stands and argues for gating the refit on the att family alone.

## 7. Honest limits of my evidence
- 12 draws. The 3.6x SE/SD gap and 5/12 coverage are far outside sampling noise, but the
  MC SDs themselves are coarse.
- All arms are bare `robust`. Part of arm N's 0.855 ame shortfall is the unabsorbed
  county effect (i=state) under robust, not the lattice; the lattice contributes the
  further factor of ~3 (0.855 -> 0.27). A cluster(cty) MC would separate the two cleanly
  and I did not run one.
- I did not find a design where the ame's analytic SE IS valid under a lattice. If one
  exists at much larger cell sizes, the ame SE loss would be a real loss there too.
- boottstat was checked and is clean (§3c). Whether PLAIN bootstrap(standard) — the
  remedy wooldid's own caution names — recovers usable numbers on a zero-V spec is in
  L4_boottstat.do; see NOTES.md for the result.
- S6 (i at the unit level with contreat()) returned rc 0 but posted nothing I could read,
  so I do NOT claim that an i above the unit is *forced* by contreat(); I claim only that
  the state-i / county-dose shape is the ordinary administrative one.

## Files
- L1_impact.do / .log   — the 12-draw MC, three arms, att and ame tracked separately
- L2_side.do / .log     — caution text, e(V) absence, downstream test/lincom
- L3_incidence.do / .log— incidence across realistic specs incl. cluster(cty)
- L4_boottstat.do       — the isolated boottstat leg
