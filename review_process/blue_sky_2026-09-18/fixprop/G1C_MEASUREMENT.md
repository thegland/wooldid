# R18 / G1c — the allocation bisection the proposal asked for, run

The proposal (`proposals_GUARD.md:187`) offered G1c as a question, not a patch, and named the
measurement that would settle it: *"for each of {reghdfe, ppmlhdfe, reg, poisson} x {binary,
contreat, contreat+lattice}, bisect `si8cells` against a fixed `maxvar` and record the break
point."* The ruling that came back was **"Tighten to the measured allocation"**, on an option
whose text reads *"Ships with the primary G1a fix"* — so R18 is G1a **and** G1c, and the frame
for reading every break point is G1a's, not the pinned formula's.

That distinction decides the answer. Read through the pinned bound the overhead is
`si8ilev + si8tlev + 200`; read through G1a's, on the absorbing path it is just `200`, because
the `i` and `t` levels go to `absorb()` and are never design-matrix columns. The same break
points mean different multipliers in the two frames.

`maxvar 2048` throughout, `safety(off)` on every measuring leg so wooldid's own guard can never
preempt Stata's allocator. Artefacts: `scratchpad/g1c/{split,skew,mult,units,ctreat,recheck,regfp}.do`
and their logs.

## 1. The axis is `si8cells`, not the event window

A design can have a wide window and few treated cells, or the reverse:

| arm | si8cells | cohorts x ALL relative periods | rc |
|---|---|---|---|
| E | 300 | 2000 | **0** |
| G | 1700 | 2000 | 907 |

E allocates fine with 2000 cohort-by-period combinations in its window; G dies with the same
2000. The allocation tracks the TREATED cells the bound already counts. G1c's axis is right.

## 2. Units are not allocated on the absorbing path

`si8cells` held at 800, unit count raised. Under the PINNED formula a `2*` bound would demand
2420 and 3020 here and refuse both:

| engine | units | rc |
|---|---|---|
| reghdfe | 140 | 0 |
| reghdfe | 600 | **0** |
| reghdfe | 1200 | **0** |
| ppmlhdfe | 140 / 600 | 0 / 0 |

This is G1a's premise, measured directly, and it is why the multiplier has to be read in G1a's
frame. It also reproduces G1a's own headline: 5,200 units x 50 cells runs clean at `maxvar 5000`
where the pinned bound demands 5458.

## 3. The multiplier, bisected with overhead constant

140 units and 20 periods in every arm, onset slid to move `si8cells` alone:

| si8cells | 300 | 500 | 700 | 800 | **900** | **1000** |
|---|---|---|---|---|---|---|
| reghdfe | 0 | 0 | 0 | 0 | **0** | **907** |
| ppmlhdfe | 0 | 0 | 0 | 0 | **0** | **3598** |

In G1a's frame the break gives `m*900 + 200 <= 2048` and `m*1000 + 200 > 2048`, i.e.

> **1.848 < m <= 2.053** — 2 sits inside it, and is the only integer that does.

The proposal's own fixture (`run_bound2.log`, break between 900 and 950) gives
`1.945 < m <= 2.053` in the same frame. The two agree, and the proposal's "the true demand is
~= 2*si8cells + ~150" was right.

## 4. But NOT 2x on top of the contreat doubling

The pinned bound already charges `si8cells*(1 + contreat)`. Multiplying that by 2 would charge
**4** per cell under `contreat()`. Measured, absorbing path, same constant overhead:

| si8cells | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900 | **1000** |
|---|---|---|---|---|---|---|---|---|---|
| reghdfe | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **907** |
| ppmlhdfe | 0 | 0 | 0 | 0 | 0 | 0 | - | - | - |

`contreat()` breaks in the **same place as binary**, between 900 and 1000, so its multiplier is
the same 2 rather than 4. A 4x charge would demand 2200 at `si8cells = 500` and refuse a run
that finishes - the exact false-positive class G1 exists to remove. **So the multiplier is a
flat 2 on the absorbing path, REPLACING the `(1 + contreat)` factor rather than multiplying
it.** `contreat()` at 600 units and `si8cells = 500` also runs, so units buy nothing there
either.

## 4b. Where the constant sits

Binary reghdfe passes at 900 and fails at 927 and 945, so `2*900 + c <= 2048` and
`2*927 + c > 2048`:

> **194 < c <= 248**

`+200` is inside that bracket, near its floor. It is therefore a fair charge rather than a round
number, and it must NOT be reduced: anything at or below 194 would under-charge and re-open the
late r(907). The residual false-positive band it leaves is `si8cells` in {925, 926} - two cells
wide at the very edge of `maxvar`, and 927 itself was measured to fail.

## 5. The non-absorbing path is left alone, and that is a measurement result

On `unconditionalse` (the reg / poisson path, where G1a keeps both level counts):

| arm | si8cells | units | rc |
|---|---|---|---|
| reg | 900 | 140 | 0 |
| reg | 1000 | 140 | **907** |
| reg | 800 | 1200 | 0 |
| reg | 800 | 1240 | **907** |

A single linear rule `a*cells + b*units + c` *can* be fitted to these four - `a = 2, b = 0.3,
c ~= 80` does it - but `b = 0.3` has no mechanical meaning, and the honest reading is that this
path is governed by **two** constraints rather than one: it breaks at 1000 cells exactly where
the absorbing path does, with units irrelevant to that break, and separately at
`cells + units + periods` near 2040, which is a charge of one column per unit. Demand there is
the larger of two terms, not a single product. With no one rule to fit, there is no "measured
allocation" on this path to tighten to, and the bound stays exactly as it is. `poisson` tracks
`reg` on every arm run.

**A residual false positive on this path, measured and NOT fixed.** At `si8cells = 800` with
1200 units the pinned bound demands 2260 against `maxvar 2048` and refuses, while the same
specification under `safety(off)` **completes** (measured both ways: `r(907)` with the guard on,
`rc 0` with it off). G1a's reg branch encodes "the levels are written out, one column each",
and that is not what the allocator does. R18 does not change it: nothing measured here says
what the right charge on that path would be, and a wrong tightening is worse than a known gap.

## 6. What ships, and every measurement it is checked against

Absorbing path: **`si8demand = 2*si8cells + 200`** (G1a drops the two level counts, G1c sets the
multiplier). Non-absorbing path: unchanged. Checked against every absorbing arm measured:

| arm | si8cells | bound | maxvar | verdict | actual | |
|---|---|---|---|---|---|---|
| mult 300 | 300 | 800 | 2048 | allow | rc 0 | ok |
| mult 500 | 500 | 1200 | 2048 | allow | rc 0 | ok |
| mult 700 | 700 | 1600 | 2048 | allow | rc 0 | ok |
| mult 800 | 800 | 1800 | 2048 | allow | rc 0 | ok |
| mult 900 | 900 | 2000 | 2048 | allow | rc 0 | ok |
| mult 1000 | 1000 | 2200 | 2048 | **refuse** | r(907) | **caught early** |
| split A | 1000 | 2200 | 2048 | **refuse** | r(907) | **caught early** |
| split B | 500 | 1200 | 2048 | allow | rc 0 | ok |
| split C | 950 | 2100 | 2048 | **refuse** | r(907) | **caught early** |
| split D | 1050 | 2300 | 2048 | **refuse** | r(907) | **caught early** |
| skew E | 300 | 800 | 2048 | allow | rc 0 | ok |
| skew F | 180 | 560 | 2048 | allow | rc 0 | ok |
| skew G | 1700 | 3600 | 2048 | **refuse** | r(907) | already refused by the pinned guard too |
| units 600 / 1200 | 800 | 1800 | 2048 | allow | rc 0 | ok |
| contreat 200-700 | 200-700 | 600-1600 | 2048 | allow | rc 0 | ok |
| G1 headline | 50 | 300 | 5000 | allow | rc 0 | **false refusal removed** |

**Zero false positives on twenty arms, four late r(907)s converted into an early named refusal
(mult 1000, split A, split C, split D), and two false refusals removed (units 1200, and the G1
headline).** The pinned bound classifies six of these wrongly; the new one classifies none wrongly.

## 7. Verified on the landed build

The G1 headline rebuilt at `maxvar 5000`, 5,200 units, 20 ever-treated: default safety now returns
**rc 0 in 1.18 s** where the pinned build refused with "requires at least 5458", and
`e(wdidmainresults)` is **bit-identical in `%21x`** across all seven cells to the same run under
`safety(off)`. The `customvce(cluster unit)` leg still refuses, as it must: that path keeps the
level counts, so its demand is unchanged in form.
