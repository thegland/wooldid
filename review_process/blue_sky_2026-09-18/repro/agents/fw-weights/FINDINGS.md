# fw-weights — FINDINGS
Agent: adversarial prober, target "fweights across the whole surface", method METAMORPHIC.
Build: /private/tmp/.../scratchpad/bluesky/build/wooldid.ado (sha c1d8979). 2026-09-17.

---

## PREMISE INVERSION (report this first)

**`fw` is not a user-facing weight type in this build.** `wooldid.ado:628-632` refuses it:

    if "`weight'" == "fweight" & "`wdrepmode'" == "" {
      dis as error "Error: wooldid does not accept fweights. Specify aweights or pweights."
      error 198
    }
    if "`weight'" == "fweight" local weight aweight

`fw` appears in the `syntax` line only so the internal bootstrap replication call can pass an
fweight while `${WOOLDID_REPMODE}` is set (ado:9048), after which the demotion on the next line
keeps margins/gstats/gquantiles on aweight semantics. The help (`sthlp:57`) advertises
`[aw/pw=weights]` and never claims fweights, so there is no doc-vs-code divergence at the
surface. Verified in probe1: `[fw=k]` gives r(198) with that message on the plain parse, under
`oldsyntax`, and with a non-integer k.

The metamorphic identity ("an fweight of k must equal k physically duplicated rows") therefore
had to be tested through the vehicles the refusal itself prescribes — `[aw=k]` and `[pw=k]` —
plus weight-scale invariance and the two internal bootstrap legs. **The identity holds
essentially everywhere I could test it; the two defects below are things I found on the way.**

---

## WHAT CAME BACK CLEAN (the metamorphic core)

Fixture: 24 units x 8 periods x 3 rows per (i,t) = 576 rows, 3 staggered onsets + 6
never-treated, integer weight `k` in 1..4 varying WITHIN cell (so the weighted and unweighted
cell means differ). Frame B = `expand k` = 1440 rows. Files `p2_A.dta` / `p2_B.dta`,
generator inside `probe2_core.do`. Max relative difference on `e(b)`:

| surface | [aw=k] vs expand k |
|---|---|
| att | 9.8e-16 |
| att + att_i + att_it + att_itime | 9.8e-16 |
| ...same under `nofastcollapse` | 9.8e-16 |
| event study (pre 2 / post 3) | 5.4e-13 (maxabs 9.3e-15) |
| ...same under `nofastcollapse` | 5.4e-13 |
| `gteffects` | 9.8e-16 |
| `subgroups()` | 2.3e-15 |
| `cluster(i)` | 9.8e-16 |
| `cluster()` finer than i | 9.8e-16 |
| `controls()` | 0 |
| `fe()` | 3.1e-15 |
| `summarystats` | 9.8e-16 |
| `customweightmultiplier()` | 9.8e-16 |
| **`[pw=k]`** vs expand k | 9.8e-16 |

**The estimand weights themselves reflect the weight, not the row count.** `e(escomp)`:
`sumw_att`, `wshare_att`, `wshare_att_it`, `ncohorts`, `nunits`, `relyr` all match EXACTLY
between the two frames, and the scalars `e(escomp_refW)`=276 / `e(escomp_trtW)`=1080 are
identical while `refN`=108 vs 276 and `trtN`=432 vs 1080 correctly differ. `e(gteffects)`:
`estimate`, `w_att`, `w_att_it` match; only `se`, `p`, `lb/ub_95ci` and `nobs` move.

**Weight-scale invariance (ado:1654's claim) holds at ordinary scale.**
`[aw=k]` vs `[aw=1000k]`: b 3.3e-16, `e(V)` 0, whole table 4.0e-16. `[aw=k]` vs `[aw=k/7]`:
3.3e-16 / 0. `[pw=k]` vs `[pw=1000k]`: 3.3e-16 / 9.9e-16. With an event study + `jointtests`:
3.3e-13. With `cluster(i)`: 3.3e-16.

**`[aw=one]` with `one` identically 1 is BITWISE identical to no weight at all** (b, `e(V)`
and the whole main table all 0.000e+00) on att/att_i/att_it/att_itime and under `cluster(i)`;
1.9e-13 with an event study. The `iw`-when-unweighted margins plumbing at ado:633-639 is
consistent with the `aweight` path.

**The bootstrap engines carry the identity.** `[aw=k]` vs `expand k`, `cluster(i)`, common
`rseed(555)`, comparing the ENTIRE displayed table:
  - `bootstrap(standard)`, `att` (frequency-count leg): b 9.8e-16, table **5.3e-15**
  - `bootstrap(standard)`, `att_i` (physical re-identification leg): b 9.7e-16, table **5.2e-15**
  - `bootstrap(bayes)`, `att`: b 9.8e-16, table **2.9e-15** — clean despite the held F001
    `/W_c` split, because the same cluster set draws the same Dirichlet mass in both frames.
  - `boottstat`, `att` (physical leg): b 9.8e-16; table 4.7e-2 but **fully attributable**:
    `se_analytic` .0527255 vs .0503719 (the aweight N/df factor), t 25.76 vs 26.96,
    **p identical at .0327869**, and the studentized interval 1.21229/1.51022 vs
    1.21257/1.51061 — a 2e-4 difference. Percentile-t is invariant exactly as theory says.

**`permute()` carries it too.** `permute(shuffle) cluster(i) permreps(60) rseed(909)`,
`[aw=k]` vs `expand k`: b 1.3e-15, whole table 7.5e-15.

**ado:9040's measured claim reproduces.** Unclustered `boottstat` with no weight takes the
INTERNAL FWEIGHT count leg (ado:9048); with `[aw=one]` it takes the PHYSICAL-EXPANSION leg
(ado:379/9236). Same data, different legs, common `rseed(4242)`, `att`: b 0, whole table
**7.2e-16**, `e(V)` 1.0e-15. The "arithmetically IDENTICAL ... on all 400 pairs" comment is
accurate on this fixture.

**A leg pair that does NOT agree, and why it is not a defect.** The same comparison on `att_i`
gives 1.4e-1 at 50 reps, 3.6e-2 / 2.2e-2 at 150 reps on two seeds. The two runs take
*different* legs: no-weight takes the re-identification leg, whose draw is a plain `bsample`
(ado:9345), while `[aw=one]` takes the physical leg, whose draw is
`bsample, cluster(`clusv') idcluster(__wdbns)` (ado:9248). ado:9296 records in its own words
that "plain bsample diverges" from the cluster+idcluster form, so the two legs draw different
samples at a common seed. The gap shrinks with reps — Monte Carlo, not arithmetic. Worth
knowing only as a seed-reproducibility wrinkle: adding a weight of 1 changes which draws a
seed produces for a cell-weighted estimand.

**`${WOOLDID_REPMODE}` leak: no path found.** With the global set by hand, `[fw=k]` IS accepted
and `level(90)` is silently ignored (`e(level)`=95 against 90 on the clean run) — but all seven
`global WOOLDID_REPMODE 1` sites (7672, 7786, 9065, 9250, 9380, 9450, 9500) sit inside
`nobreak { ... cap ... global WOOLDID_REPMODE }`, and `wooldid_worker` consumes and clears the
global at entry (ado:538-539), so even a leak is eaten by the next call. Not filed.

**Parse edges, all as documented.** `[aw=0]` on 20 rows -> N 576->552; missing weight -> same;
negative -> r(411) (VAL-G06). Non-integer `fw` -> the fweight refusal (wooldid's message fires
before Stata's integrality check).

**Population-scale weights: the cell-weight totals are exact.** ado:1619/1620 build the
denominators `w_per_it` / `w_per_i` with `gegen ... = total(...)` and NO `double`, while the
escomp twin at ado:3833 writes `gegen double`, and commit `3cb0c41` ("B21: gen double for
ringer, interactive-control, and estimand-weight copies") converted the four `w_att*`
NUMERATORS to `gen double` and left these two denominators untouched — its own message scopes
itself to "the four w_att* estimand weights". That looked like an incomplete B21 fix. It is
not: **gtools' `gegen` defaults to `double`** where Stata's `egen` defaults to `float`
(measured directly in probe10: `GEGEN_TYPE=double`, `EGEN_TYPE=float`). So with county-scale
weights (2.5e5-9.3e6, cell totals ~1e8) and with survey weights summing to ~3e8, `[aw=pop]`
vs `[aw=pop/1e6]` vs `[aw=pop/1e3]` agree to 3.3e-16 - 9.3e-15 on all four estimands and on
the event study. A float denominator would have shown ~6e-8 here, six hundred times the
project's 1e-10 threshold. Recorded as a latent fragility only: those two lines are correct
by gtools' default, not by their own spelling, and a fallback to `egen` would silently
reintroduce B21 on the denominator side.

**`permute()` + `permtstat` settled (probe11).** The 4.7e-2 table gap there is the identical
DOF factor. `[aw=k]` vs `expand k` at `rseed(909)`, `permreps(60)`:
`se_analytic` .05272545083 vs .05037188006 and `t_perm` 25.757 vs 26.961 — but `p_perm`
.03278688525 on both, the `level()` null band 1.207006082/1.555633593 identical to every
printed digit, and `e(perm_p_att)`/`e(perm_p_att_i)` **bitwise identical**
(`+1.0c9714fbcda3bX-005`). Under `permute()` N and K are constant across replications, so the
aweight rescales every t by the same factor and the ranks are invariant. Exactly right.

**Zero weight == deleted row (probe7).** `[aw=w]` with `w` in {0,1} against the dataset with
the zero-weight rows physically dropped: identical `e(N)` (494 = 494), identical `e(df_r)`
(391), and b / whole table / `e(V)` agreeing to 3.0e-16; with the weight `k*keeprow` against
`[aw=k]` on the dropped frame, **bitwise 0**; and the `if` route (`if keeprow == 1`) also
bitwise 0. The markout at ado:683 and BS-37's exact-one `if` rule agree perfectly.

**A constant `customweightmultiplier()` is a no-op (probe7).** `cwm(1)` vs `cwm(2)` is bitwise
identical on all four estimands; `cwm(1)` vs no cwm is 3.3e-16.

**More surface, all clean on `e(b)` (probe7).** `impvalues cluster(i)` 8.1e-16 (and the WHOLE
table 1.4e-14, since the IM p-values are scale-free); `anticipation(1)` 3.9e-16;
`oostreatedcontrols(1)` 9.8e-16; `esfixedbaseperiod` 3.7e-13; `esrelativeto(-2)` 5.4e-13;
`timetrends` 1.1e-14; `unconditionalse robust` 1.1e-15; `interactivecontrols()` 8.4e-16.
Their SE columns move by the aweight N/df factor, as they must (df_r 455 vs 1319 etc.).

**poisson `[pw=k]` vs duplication (probe7/probe13).** `e(b)` 1.7e-15 / 2.2e-15; `aweight`
under `poisson` is refused r(198) as sthlp:710 says. (probe7's SE comparison there was void -
my first fixture's outcome was a deterministic `exp()` of the model, so ppmlhdfe fit it
exactly and every SE was machine zero, se 1.2e-11 and z 1.7e10, in BOTH frames. probe13 redoes
it on `rpoisson()` draws.)

**Weighted quantiles agree bitwise (probe15).** The lattice cuts the dose at WEIGHTED
quantiles (sthlp:922), the classic aw/fw divergence site. `_pctile dose [aw=k]`,
`_pctile dose` on the duplicated frame, and `_pctile dose [fw=k]` all return the same
median to the last bit: `+1.0cf1acc961b23X+001`.

**The finer-than-i machinery shipped 2026-09-16 carries the identity (probe9).**
`att_i` + `cluster()` FINER than `i` + `bootstrap(standard) reps(60) rseed(31415)`,
`[aw=k]` vs `expand k`: b 9.7e-16 and the WHOLE table 3.4e-15. The `boottstat` twin there
is 6.3e-2, the same analytic-SE DOF factor as every other boottstat pair.

**Continuous treatment (probe9/probe15).** `contreat() contreatcontrols(basic)`: b 3.4e-14;
with an event study 1.5e-13; under `nofastcollapse` 8.6e-12 (maxabs 7.9e-14).

---

# FINDINGS

## FW-01 — a weight written as a `+`/`-` EXPRESSION silently corrupts the cell-weighted estimands on the `customvce()` / `unconditionalse` paths
**severity medium · confidence certain · realistic true (three-way conjunction, stated below) · FOUND BY ME**

TWO DISTINCT BEHAVIOURS, do not conflate them: on the DEFAULT estimation path a weight
expression ERRORS (reghdfe's own guard); on the `customvce()` / `unconditionalse` paths it is
ACCEPTED and mis-binds. A verifier who runs only `wooldid ... [aw=a+b], att_it` will see
r(198) - that is the masked case, not a refutation.
repro: `/private/tmp/.../probes/fw-weights/repro_wexpr.do`

`wooldid_worker` keeps the raw weight-expression TEXT in `` `w' `` (ado:607
`local w \`exp'`) and pastes it into arithmetic:

    ado:1619  qui gegen `w_per_it' = total(`esamp' * `w') , by(`i' `t')
    ado:1630  qui gen double `w_att_i'     = `w' `customweighttext'/`w_per_i'
    ado:1637  qui gen double `w_att_itime' = `w' `customweighttext' /(`w_per_it'*`t_per_i')
    ado:1644  qui gen double `w_att_it'    = `w' `customweighttext'/`w_per_it'
    ado:1650  qui gen double `w_att'       = `w' `customweighttext'

With `[aw = wf+wm]`, line 1644 evaluates `wf + wm/W_it`, not `(wf+wm)/W_it`. `att` has no
denominator, so it is exact; `att_it`, `att_i` and `att_itime` are wrong. `customweighttext`
is `* cwm`, so with `customweightmultiplier()` even `att` becomes `wf + wm*cwm`.

OBSERVED (probe17, `wf+wm` constructed to equal `wtot` exactly, `customvce(robust)`, no
`safety(off)` needed):

| | `[aw=wtot]` | `[aw=wf+wm]` |
|---|---|---|
| att | 1.35806635445 | 1.35806635445 |
| att_it | 1.35539110268 | **1.35837510195** |
| att_i | 1.36654660882 | **1.35953352747** |
| att_itime | 1.36252940698 | **1.35932150079** |

max relative 5.13e-3, max absolute 7.01e-3 — five hundred thousand times the project's 1e-10
threshold. `customvce(cluster i)` moves 2.20e-3; `unconditionalse robust safety(off)` moves
5.13e-3. No error, no caution, rc 0.

EXPECTED: identical, because `wf+wm == wtot` on every row (asserted in the repro).

WHY IT IS NOT SEEN ON THE DEFAULT PATH: reghdfe refuses the expression first — traced in
probe16, `reghdfe.Estimate` runs `unab exp : wf+wm, min(1) max(1)` and dies "+ invalid name"
r(198). So on the default path this is instead **an uninformative error**: `[aw=pop/1000]`
gives r(198), `[aw=2*k]` gives "variable 2*k not found" r(111), `[aw=int(k)]` and
`[aw=cond(k>2,k,1)]` give r(101), `[aw=(ka+kb)]` r(198) — none of which mentions weights, and
`safety(off)` does not change any of them. Plain `regress y post i.i i.t [aw=wf+wm]` accepts
the same expression (rc 0), which is why a user would write it. `customvce()` and
`unconditionalse` estimate with plain `regress`, which is exactly where the guard is missing
and the bad arithmetic runs.

FIX SHAPE: materialise the expression once — `tempvar wv ; qui gen double \`wv' = (\`exp')`
— and use `` `wv' `` everywhere. That fixes the arithmetic AND turns every message above into
a working run.

## FW-02 — pweights give a robust covariance but the run is labelled OLS in `e(vce)`
**severity medium · confidence certain · realistic true · FOUND BY ME**
repro: `/private/tmp/.../probes/fw-weights/repro_pw_vce_label.do`

`sthlp:268` says pweights "imply robust standard errors" and that is true of the NUMBERS:
with `[pw=svywt]` and no vce option, `e(V)` is bitwise the covariance the same run produces
under `robust`.

    A  [pw=svywt]              se = 0.0919359312   e(vce)=[ols]     e(vcetype)=[.]
    B  [aw=svywt], robust      se = 0.0919359312   e(vce)=[robust]  e(vcetype)=[Robust]
    C  [aw=svywt]  (iid)       se = 0.0817970344   e(vce)=[ols]     e(vcetype)=[.]
    se(A) == se(B) exactly (|diff| 0.0000e+00);  se(A) != se(C) (|diff| 1.0139e-02)

So the covariance is B's and the label is C's. The header prints the blank
`Standard Errors: ` as well. The contract this breaks is stated in wooldid's own comment at
ado:4283 — *"e(vce) and e(vcetype) name what e(V) IS"* — and the fallthrough at ado:4301,

    else ereturn local vce = cond("`poisson'" == "poisson", "oim", "ols")

never looks at the weight type. ATTRIBUTION: this is wooldid's, not inherited — reghdfe on
the identical pweighted fit posts `e(vce)`="robust", `e(vcetype)`="Robust", and wooldid
replaces that with "ols".

WHY IT MATTERS: `e(vcetype)` is exactly what estout/esttab and outreg2 read to decide whether
to print "Robust standard errors" under a table, so a pweighted wooldid result is published
labelled as OLS. Survey data is the canonical pweight use case.

NOT the same as PS-33 (the ruled cosmetic: the iid default's header carries no label). Here
the vce is NOT the iid default and a machine-readable `e()` macro makes a false positive
claim. The one-line fix is to add the weight-type branch beside the `robust`/`cluster` ones
already at ado:4288-4293.

## FW-03 — the fweight refusal prescribes a remedy that does not reproduce fweight inference, and nothing says so
**severity low · confidence certain · realistic true · SEEDED (the refusal itself is the seed's premise)**
repro: `/private/tmp/.../probes/fw-weights/repro_collapsed_admin.do`

`[fw=n]` is refused (ado:628-632) with *"Specify aweights or pweights."* Neither substitute
reproduces what the fweight would have given, because neither carries the row count.

Fixture: individual administrative microdata, 30 states x 12 years x ~25 people = 8941 rows,
collapsed the ordinary way to 360 (state, year) cells with a mean outcome and a headcount `n`.

    microdata,      cluster(state):  b .38464638  se 0.02936881  N 8941
    [fw=n] on cells:                 r(198) refused
    [aw=n] on cells, cluster(state): b .38464638  se 0.03673877  N 360   se ratio 1.2509
    [pw=n] on cells, cluster(state): b .38464638  se 0.03673877  N 360   se ratio 1.2509
    microdata,      iid:                        se 0.04260507
    [aw=n] on cells, iid:                       se 0.03885874        se ratio 0.9121
    [pw=n] on cells, iid:                       se 0.04068764        se ratio 0.9550

The point estimate is exactly right (reldiff 0.0000e+00 — the identity this whole probe
verified). The standard error is 25% too LARGE clustered and 9% too SMALL iid, and `e(N)`
reports 360 where the analysis sample is 8941. There is no way to get the microdata answer out of the
collapsed panel WITHOUT physically expanding it (`expand n` reproduces it exactly - that is
the metamorphic identity this whole probe verified).

This is aweight semantics, which are Stata's and not wooldid's, so it is not a wrong number
in wooldid's own terms. What is wooldid's is (a) the remedy the message names is not
equivalent and the message does not say so, and (b) `sthlp:57` documents `[aw/pw=weights]`
and the help nowhere mentions frequency weights, count data, or collapsed panels — so a user
who has not read ado:620 has no way to learn that the substitution changes their inference.
A sentence in the error message ("aweights and pweights do not carry the row count, so the
standard errors and e(N) will be those of the collapsed panel; expand the data if you need
frequency-weight inference") would close it.
