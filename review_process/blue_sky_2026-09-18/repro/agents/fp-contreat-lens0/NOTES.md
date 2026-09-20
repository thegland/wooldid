# fp-contreat-E — LENS 0, INDEPENDENT REPRODUCTION

Repro: `L0_untreatedZero.do` / `L0_untreatedZero.log` (single run, all 16 wooldid calls rc=0).
Build asserted: `.../bluesky/build/wooldid.ado`, version 1.4 2026-08-19. Suite not used.

## VERDICT: CONFIRMED, severity critical. Reproduces exactly on my own DGP, in closed form.

My DGP is my own: `i` = 8 groups, 15 sub-units each, t = 1..9, onsets 4/5/6/7/5 for groups 1-5,
groups 6-8 never treated, dose exactly 0 off treated cohort-periods (the untreatedZero premise),
dose ~ U(0.5,3) inside them. 1080 obs, 23 treated (i,relyr) cells of 15 obs each.

## 1. Homogeneous exact-fit design (truth: ame = 0.80 under ANY cell weighting)

| run | att | ame |
|---|---|---|
| `robust` (default reghdfe) | 3.3725563846 | **0.8000000000** |
| `ccnoabsorb robust` | 3.3725563846 | **0.8000000000** |
| `customvce(robust)` | 1.9869565217 | **-0.0000000000** |
| `customvce(cluster grp)` | 1.9869565217 | **-0.0000000000** |
| `unconditionalse robust` | 1.9869565217 | **-0.0000000000** |

The three gated paths return exactly ZERO for a true effect of +0.80. Mechanism check: with a
homogeneous beta the inserted pooled slope absorbs it and every cell slope becomes 0.

ATT check in closed form: mean cell intercept over treated obs = **1.9869565217** (python,
independent of Stata). That is exactly what the gated paths report as the ATT — i.e. the ENTIRE
dose contribution (0.8 * mean dose = 1.3856) has been stripped out of the headline ATT.
**So the ATT row is corrupted too, not just the ame row.**

## 2. Same design with N(0,0.3) noise — the realistic case

Same data, same specification, only the variance option differs:

```
default robust :  att 3.4056883 (se .0349)   ame  .80437912 (se .0244, p 0.000)   [truth .80]
customvce(rob) :  att 1.8794300 (se .2792)   ame -.07758669 (se .1568, p 0.621)
```

The wrong answer carries an entirely ordinary standard error and reads as a clean, publishable
NULL ("dose has no effect, p = 0.62") when the truth is +0.80 — 5.6 of its own SEs away. The ATT
moves 3.41 -> 1.88, 5.5 of its own SEs. Nothing in the printed table looks broken.

## 3. Heterogeneous exact-fit design — the identity, in closed form

True cell slopes bcell = 0.50 + 0.13*grp + 0.07*relyr, all 23 distinct.

- closed-form truth (python, cell-obs mean) = **1.0013043478**; wooldid `robust` = 1.0013043478 (exact)
- `customvce(robust)` and `unconditionalse robust` = **-0.4286956522**
- truth − bcell(grp=5, relyr=4) = 1.0013043478 − 1.4300000000 = **-0.4286956522** (exact)

So the reported number is the truth minus the true slope of whichever cell the collinearity screen
happened to drop. Independently reproduced; I never looked at the other agents' numbers to get it.

## 4. THE CONTROL NOBODY RAN — kills the "the engines just differ" defence

Separate design, dose defined and nonzero EVERYWHERE, `contreatcontrols(basic)` (where the pooled
baseline control is legitimate and asked for):

| run | att | ame |
|---|---|---|
| `robust` | 3.4302233064 | 0.8000000000 |
| `customvce(robust)` | 3.4302233064 | 0.8000000000 |
| `unconditionalse robust` | 3.4302233064 | 0.8000000000 |

Bit-identical to 10 decimals. `customvce()` changes only the variance estimator on an ordinary
specification, exactly as documented. The divergence is specific to `untreatedZero`.

## 5. Poisson — the prober's mechanism wording is WRONG, the skeptic's is right

- Homogeneous (true semi-elasticity 0.25 everywhere): reported ame = **3.968e-17, p = 0.999**;
  att = 0.2 exactly (my post intercept; the 0.25*dose part stripped).
- Heterogeneous (bcell = 0.10+0.04g+0.02r): closed-form truth **0.2513043478**, reported
  **-0.1286956522** = truth − bcell(5,4)=0.38, exact. **Not zero.**

So "driven to zero regardless of the truth" is the homogeneous special case only. The visible
mechanism in the log is the fallback the skeptic named:
`Caution: initial estimation of ppmlhdfe failed, possibly due to inclusion of an absorbed
categorical fixed effect with only 1 level.` — that one-level term is wooldid's OWN
`i.<null>##c.(ringer)`. The fallback message then says *"In principle, this fallback solution should
usually yield the same results"*, which is false here.
There is NO uncorrupted poisson path: the gate at ado:1234 fires on `poisson` alone.

## 6. What the user sees

A generic multicollinearity caution fires on exactly the 9 corrupted runs and on none of the 7
correct ones — so wooldid does know something was dropped. But it is the same generic caution used
elsewhere, it says only *"This could imply your results are not identified, depending on exactly
what was dropped"*, and the table below it prints the corrupted att and ame as ordinary numbers
with ordinary SEs and p-values. Nothing tells the user the point estimate is wrong.

## Scope / realism
`contreat() contreatcontrols(untreatedZero) poisson robust` on a count outcome that is zero before
treatment is sthlp:1551's own worked example — no exotic option, and the gate fires on `poisson`
alone. sthlp:1033-1034 states untreatedZero means g() = 0 and sthlp:986 says it "cannot be paired
with any other type", so the user cannot even ask for the control that is inserted.
