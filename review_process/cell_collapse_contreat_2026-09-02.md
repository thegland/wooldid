# Cell collapse and contreat(): does the cell mean of the dose make the collapse exact? (2026-09-02)

Assessment requested by the author on 2026-09-01: "assess the feasibility of extending our
approach to continuous treatments. Assume no lattice() or interactive() controls or anything
that also breaks the binary case. Would it work if after preserving, we assigned the first row
of each i-t cell's continuous treatment variable to the cell mean of the continuous treatment
variable?" The prototype lives in the session scratch directory and is reproduced in full in
section 9.

**Status: RULED AND BUILT, 2026-09-02.** The author ruled after reading sections 1-8: build the
collapse for `contreat()` with the gates below; gate `poisson` + `contreatelasticitytype(dyex)`;
take the one-token `wsum != 0` guard; and treat `customweightmultiplier()` as caveat emptor, so
the Σw = 0 detector of section 5.1 stays built but unshipped. That is what is now in
`wooldidinstall/wooldid.ado` - the shipped code differs from the section-9 prototype only by
those two lines, plus rewritten comments - together with a rewritten refusal message and a
rewritten help section ("What turns the collapse off: the complete list") that enumerates every
route by which the collapse does not run. Section 6a answers the four questions the author raised
on reading this document. Companion to
`cell_collapse_2026-08-31.md`, `cell_collapse_followup_2026-09-01.md` and
`cell_collapse_subgroups_2026-09-01.md` (route P, whose working-tree ado is the base of the
prototype).

**Headline.** Yes - the cell mean is the right answer to the author's question, with four
refinements to it plus one implementation proviso (the ringer copy of the dose is averaged the
same way; see (vi)), and with those in place the collapse is exact for `contreat()` over the whole
domain the help file already documents for the binary collapse. The refinements: the mean must be
the **estimand-weighted** mean, not Σx/n (i); it must be formed as Σ(w·x)/Σw and written onto
**every** row before the row pick, not by `gcollapse (mean) [aw=w] merge`, which leaves a
zero-weight first row holding its raw dose (ii); the user's cohort variable `i` must join the
collapse key (iii); and the dose must be **recast to double** before the mean is stored (iv).
Each was tested by deleting it from the prototype: each deletion produces divergence, and three
of the four produce divergence in numbers wooldid prints - a subgroup pre-trend p moving
0.0255 -> 0.0527 with its CI bound changing sign (i), a pre-trend horizon moving 0.167 -> 0.0022
(ii), a joint pre-trend p moving 0.0956 -> 0.1112 (iii). The form the collapse would have had
without the author's proposal - the surviving row keeping its own dose, which is what the shipped
gate was written against - fails on 9 of 9 contreat configs, and not marginally: on the smallest
possible request the ATT reads 1.3775620504 against 1.3836141895 with its variance overstated
11.5%, p-values cross 0.05 in both directions, and six coefficients flip sign. So the mean is not
an optimisation of the collapse; it is what makes it correct.

With the four refinements in place, 328 paired collapsed-vs-full runs over seven adversarial
angles (estimands x weights x dose storage; model specification; every event-study, pre-trend and
diagnostic read; subgroups; gates and nonlinear forms; VCE, bootstrap and permute; data shapes)
give **no divergence on the documented domain**. Comparing every `e()` object in `%21x` over the
305 of them that both produced dumps and are not accounted for individually below - 846,157
matrix cells, 8,881 scalars and every macro - 8.7% of matrix cells differ in their last bits and
**not one differs by more than 1e-9 relative** at any magnitude down to 1e-12; 40 configs are bit-identical throughout; the
only scalars moving more than 1e-9 are four joint-test p-values of order 1e-38 and 1e-21. The
largest relative move on a cell of ordinary size is 1.9e-10 (a p-value of 1.2e-8), the largest
absolute move anywhere is 4.8e-10 (on a t statistic of 39.2498); `e(b)` tops out at 1.1e-10
relative and `e(V)` at 2.6e-11.

Three things sit outside that. A *signed* `customweightmultiplier()` under `safety(off)` whose
cell weight sum is negative (repaired by a one-token guard change, measured as a no-op elsewhere)
or exactly zero with non-zero rows (a limit of the summed-weight representation; a 15-line
detector is built and priced, and the author has ruled that cwm is caveat emptor, so it stays on
the shelf). `poisson` with `contreatelasticitytype(dyex)`, where margins differentiates
`expression(predict(xb))` numerically and the collapse perturbs the ame variance by the same order
as the input row order does - gate it or document it. And one regime found by trying to break the
prototype rather than to confirm it: when a cell's doses cancel to four or more digits (within-cell
spread >= ~1e4 times the cell mean) the two arms part company - and every anchor says the
**uncollapsed** arm is the one that is wrong, to the point that at spread/mean 1e10 the full-data
run moves a reported p-value from 0.062 to 0.573 on row order alone while the collapsed run holds
its analytic value to ten digits (section 5.3). That regime is pathological fixture territory and
argues for no code, but it does mean a battery MISS there is not evidence against the collapse.

The prize is the binary prize: 1.27-1.44x on a 4k-row admin sample, 2.2-4.2x at 30k rows,
3.7-4.1x at 145k rows, with `mreldif(b)` ~1.7e-15.

---

## 1. Why contreat() was gated, and what changes

The shipped gate (`wooldid.ado` 4694: `if "\`contreat'" != "" local fcbad \`fcbad' contreat()`)
was written when the first collapse measurement showed contreat runs "coefficients shifted,
one sign flip". That measurement was of the *first-row* collapse: the surviving row of each
cell kept its own dose, so every margins read that involves the dose saw one draw of a
within-cell-varying regressor instead of the cell. The author's question is whether replacing
that draw with the cell mean repairs it. It does, for a reason worth writing down because it
decides which options must stay gated.

**Every collapsible read is linear in the dose within a cell.** With `contreat(x)` the model
carries, for each cohort-by-period treatment cell g, an intercept term δ_g·D and a slope term
β_g·D·x (plus the `contreatcontrols()` terms, which are absorbed or explicit regressors that
margins does not differentiate). Inside the margins window the reads are:

| read | quantity at a row of cell g with dose x | within-cell form |
|---|---|---|
| `dydx(targetflag)` (att, att_it, att_i, att_itime, ES horizons, pre-trend rows, joint tests, adversarial line) | δ_g + β_g·x | linear in x |
| `dydx(contreat)` (the `ame` rows, `ame_it`, `V_es_ame*`) | β_g | constant |
| `dyex(contreat)` | β_g·x | linear in x |
| `eyex(contreat)` | β_g·x / ŷ | nonlinear (ŷ depends on x) |
| `eydx(contreat)` | β_g / ŷ | nonlinear |
| `contreatpoly(k>1)` | β1·x + β2·x² + ... | nonlinear |
| `lattice()` / `cfxlattice()` | any of the above, `over()` a binning of x | splits the cell on a non-key variable |

A margins read over a subpop is Σ_rows w_r·(a_g + b_g·x_r) / Σ_rows w_r. Grouping rows by
cell, Σ_g [a_g·W_g + b_g·Σ_r w_r x_r] = Σ_g W_g·(a_g + b_g·x̄_g) with x̄_g = Σ w x / Σ w. The
right-hand side is precisely the read on a collapsed dataset holding one row per cell with
weight W_g and dose x̄_g. The delta-method Jacobian is the derivative of the same expression
in (a_g, b_g), so `e(V)` matches too. Three conditions make this exact rather than
approximate: (1) a_g and b_g are constant inside a collapse cell, i.e. the key contains every
variable that indexes the treatment terms and every variable a `subpop()`/`over()` splits on;
(2) the weight used to form x̄_g is proportional, inside the cell, to the weight the read
uses; (3) the read is one of the linear rows of the table. Each of the prototype's provisos is
one of these conditions, and each gate keeps out a read that violates (1) or (3).

**What the prototype does** (full diff in section 9; ~70 lines against the route-P working
tree):

- Collapse site (shipped 1835): when `contreat` is set, the user's cohort variable `i` joins
  the key (`local fckeyx \`i'`) and `contreat(\`contreat' \`ringer')` is passed to the worker.
  The ringer (wooldid's double copy of the dose) is averaged too, because wherever it reaches the
  regression its level shifts the round-off of every prediction. It reaches it through
  `contreatcontrolinserterFE`, built at 1123-1141 under `contreatcontrols(basic)` without a
  lattice, under `byTreated`/`byTtre`/`byCohort`/`byTime`, and unconditionally under
  `unconditionalse`, `customvce()` or `poisson`; `ccnoabsorb` (1157-1159) moves those terms out of
  the absorb into the explicit regressor list. `contreatcontrols(untreatedZero)` **alone** puts the
  ringer nowhere - it only rewrites `contreatinj` (1152-1154) - so the prototype's own comment,
  which names `untreatedZero` as a ringer path, is wrong on that point (section 7).
- Worker `wooldid_fcol` (shipped 5133): recasts each averaged variable to `double`; forms
  Σ(w·x) and Σw per key cell with `gcollapse (sum) ... merge` where w is the *first existing*
  estimand weight (`w_att`, else `w_att_it`, ...); writes Σ(w·x)/Σw onto every row of the cell
  `if wsum > 0 & wsum < .`; then the existing row pick keeps the first row and sums the
  estimand weights as before.
- Gate (shipped 4694 replaced): `contreat()` itself no longer refuses. Refused with
  contreat: `contreatpoly(>1)`, `lattice()`, `makecfxplots`/`cfxlattice()`,
  `contreatelasticitytype(eyex|eydx)` (dyex allowed), and the contreat variable appearing
  inside `controls()`, `coarsecohortcontrols()`, `fe()` or `ccc_absorb()` (token match after
  stripping factor-variable prefixes; gate name `contreat-in-controls`). The pre-existing
  gates (`interactivecontrols()`, `unconditionalse`, `semielasticity`, `impvalues`,
  `poisexpresults`, `customvce()`) are unchanged. Poisson under `dydx`/`dyex` is allowed
  because wooldid already forces the linear index there (shipped 1804-1807: `predict(xb)` under
  `dydx`, `expression(predict(xb))` under `dyex`).

## 2. Provisos and why each is necessary

The provisos were tested two ways: the prototype was probed with every one in place (section
4), and four of the six were removed in turn in a variant ado to show the failure each prevents.
Proviso (v) is the base idea rather than a refinement and was tested by removing the mean
altogether; proviso (vi) rests on the code path plus sufficiency evidence, not on a deletion
test.

| proviso | what breaks without it | evidence |
|---|---|---|
| (i) estimand-weighted mean, not Σx/n | any aw/pw/cwm that varies *inside* a cell, on every read of the att family (not just the elasticity forms): estimates move 0.15% with the fixture's own `analytic_weight` and 15% with a strong weight, and one measured subgroup pre-trend crosses 0.05 with its CI bound changing sign | V1a(b), section 2.1 |
| (ii) Σ(w·x)/Σw written on every row before the pick | a zero-weight first row keeps its raw dose (`gcollapse (mean) [aw=w] merge` writes missing on zero-weight rows, and the pick row is the cell's first row whether or not its weight is zero) | V1b(f), section 2.2 - reproduced on the first prototype's form: a pre-trend horizon moves p = 0.167 -> 0.0022; re-verified fixed in round 2 from six angles (A1: cwm 0 on the first row with pw, with within-cell aw and no att, with a 1e2 dose offset, combined with an if() inside cells; cwm 0 on the last row; cwm 0 on a whole cohort - all PASS) |
| (iii) user cohort `i` in the key | the pooled reference-period cells (`_i==0`) are read under an estimand weight that is not proportional, across the cohorts sharing an onset, to the weight the mean was built with; the whole-sample pre-treatment joint-test p-values then move (one measured case crosses 0.10) | V1a(a), section 2.3 |
| (iv) recast to double inside the preserve | a float dose stores the cell mean at float precision, injecting ~2.4e-8 relative into every dose-scaled read (no measured p-crossing or CI sign flip, but it breaks bit-reproducibility and turns 12 per-subgroup joint p-value *macros* into structural differences) | V1a(c), section 2.4 |
| (v) the mean at all (vs the shipped pre-mean form, where the surviving row keeps its own dose) | the shipped "coefficients shifted, one sign flip" - and worse: 9 of 9 contreat configs diverge, the plain `att contreat()` ATT moves from 1.3836141895 to 1.3775620504 with its variance overstated 11.5%, p crosses 0.05 in both directions, six coefficients flip sign | V1a(d), section 2.5 |
| (vi) the ringer copy averaged with the dose | its level enters `xb` wherever `contreatcontrolinserterFE` is non-empty (1123-1141: `basic` without a lattice, `byTreated`/`byTtre`/`byCohort`/`byTime`, or any of `unconditionalse`/`customvce()`/`poisson`), and `ccnoabsorb` moves those terms into the explicit regressor list | code path (1123-1159) plus sufficiency: A2's `ccnoabsorb`, `byCohort`, `byTime` and every poisson config PASS. Not tested by deletion, and every poisson config on this fixture took the ppmlhdfe fallback, so the primary absorbing ppmlhdfe path with contreat is unmeasured |

**What "estimand-weighted" means, exactly** (the weights are built at `wooldid.ado` 1426-1457):
`w_att = w · cwm` - the user's sampling weight times `customweightmultiplier()` - while
`w_att_it = w · cwm / w_per_it`, `w_att_i = w · cwm / w_per_i` and
`w_att_itime = w · cwm / (w_per_it · t_per_i)`. Every denominator there is a total taken **by
(i,t)** or **by (i, posttreatment)**, so it is constant inside an (i,t) cell. The worker weights
the dose mean by the *constructed* estimand weight (the first one that exists), but since all
four are `w · cwm` divided by a within-cell constant, and since proviso (iii) puts `i` in the key
so that every collapse cell sits inside one (i,t) cell, **the estimand-weighted mean and the
plain `w · cwm`-weighted mean are the same number**. The distinction only bites where a collapse
cell spans more than one (i,t) cell - which is precisely the pooled reference rows that proviso
(iii) removes (section 2.3).

Why the *first existing* estimand weight suffices for (i): inside a genuine (i,t) cell every
estimand weight is w × cwm × (a per-cell factor), so all of them give the same weighted mean.
A1 exercised every branch of the selection (att only; att_it only; att_i only; att_itime only;
att_i att_itime without att; all four) with a weight that varies strongly within cells
(`0.05 + 5·runiform()` as aw and as pw) and with cwm: all PASS. The proportionality is exactly
what fails when cohorts pool, which is what makes proviso (iii) load-bearing - see section 2.3,
where dropping `i` from the key makes a read weighted by `w_att_it` consume a mean built under
`w_att` across three pooled cohorts.

### 2.1 Proviso (i) measured: the weighted mean, not Σx/n

V1a(b) ran a variant whose only change is that the worker's weight source is left empty, so the
cell mean is Σx/n, over thirteen configs on `continuous_control_dense`: **10 MISS, 3 PASS**. The
condition is a clean product of two factors: the proviso is load-bearing exactly when the
estimand weight w × cwm varies *within* a collapse cell **and** the dose varies within that same
cell. Remove either and the two means coincide by construction.

- Smallest breaking configuration: `Vb_dyex_min_att_wbig` - `att` only, `robust`,
  `esprelength(1) espostlength(1)`, no cluster, no joint tests, no subgroups, one weight
  (`aw=wbig`): `pretreatment:att` 0.0743442818 vs 0.0799497572 (7.0% relative).
- With the fixture's own stock weight rather than a constructed one (`Vb_dyex_aw_stock`,
  `aw=analytic_weight`, within-cell sd 0.069-0.080): 133 matrix cells move beyond 1e-9, worst
  estimate 0.15% (`pretreatment:att` -0.0317309228 vs -0.0316836921). This is the realistic
  floor, not the constructed one.
- Decision-changing case (`Vb_dyex_sg3_aw_wbig`): the subgroup-1 `ES_att:pre1` row reads
  estimate 0.0673367818 / p **0.0527** / CI [-0.00080, 0.1355] under the unweighted mean against
  0.0780380401 / p **0.0255** / CI [+0.0097, 0.1464] on the full data - the CI lower bound
  changes sign and the pre-trend test crosses 0.05 in the *masking* direction. The `att_it` row
  of the same table crosses independently (p 0.0693 vs 0.0354).
- Closest non-breaker: `Vb_ctl_dyex_aw_wcell` - the identical request with `aw=wcell`, a weight
  that varies *across* cells but is constant within them (within-cell sd 2.8e-16): worst residue
  2.4e-13, zero cells beyond 1e-9. That single substitution turns a 15% MISS into a PASS.

The angle also corrected its own pre-registration: the two configs predicted to PASS under
`dydx` (on the argument that `dydx` is level-blind) both MISSed, and their att-row moves
reproduce their `dyex` twins' to the printed precision while their ame rows sit at <= 3.7e-15.
The reason is the regression specification - the treatment enters as
`i._i#i._t#i.targetflag##c.contreat`, so `dydx(targetflag)` = β_level + β_slope·x is
level-sensitive under **every** `contreatelasticitytype`, the default included; only
`dydx(contreat)` is level-blind. So proviso (i) is load-bearing on ordinary runs, not just on
the elasticity forms.

### 2.2 Proviso (ii) measured: the first prototype's `gcollapse (mean) [aw=w] merge`

The defect round 1 found was reproduced deliberately, on a variant ado that restores the
original mean form (`ct/V1b/wooldid_v1_oldmean.ado`; everything else identical to the
prototype), over four configs on `continuous_control_dense` with all four estimands,
`cluster(person_id) esprelength(2) espostlength(2) jointtests` (`adv2/ctV1b_f/RESULT.txt`):

| config | verdict | what moved |
|---|---|---|
| `V1bf_w0f_cwm` (cwm = 0 on the first row of every cell) | MISS, 160 cells | `ES_att:pre1` estimate 0.0334412 -> 0.0780537, se 0.0240993 -> 0.0250416, **p 0.167405 -> 0.0022092**, 95% CI [-0.0141957, 0.0810781] -> [+0.0285543, 0.127553]; `pretreatment:att` -0.03559444 -> -0.00287964 (92% of its own magnitude); `main:att` 1.38151 -> 1.37451 |
| `V1bf_w0f_pw_cwm` (same, with pw) | MISS, 160 cells | `ES_att:pre1` p 0.161364 -> 0.00168835, CI lower -0.0135941 -> +0.0303515 |
| `V1bf_w0last_cwm` (cwm = 0 on the **last** row of every cell) | PASS | worst `e(b)` 4.7e-14 - the pick row has weight, so the old form is exact |
| `V1bf_nocwm_aw` (no zero weights) | PASS | worst `e(b)` 2.8e-14 |

That is the sharpest statement of why (ii) is a proviso and not a detail: the failure mode is
silent and it lands on inference, not just on the last decimal - a pre-trend horizon that the
full data calls insignificant (p = 0.17) is reported as significant at the 1% level, with a
confidence interval that no longer contains zero.

**How the pick row gets a zero weight at all** (measured 2026-09-02, `ct/Q/qC_zeroweight.log`,
`qD_nsample.log`). The worker's first act is `keep if \`varlist' == 1` on the e(sample) marker
(shipped 5136), and the surviving row is the minimum original row number *within the sample*
(`(min) pick = ord` at 5147), chosen without reference to any weight. So the hazard needs a row
that is **in** e(sample) and carries zero estimand weight - and there is exactly one route to
that. A zero *sampling* weight does not do it: Stata drops zero-`aw` observations from the
estimation sample, measured as `e(N)` falling from 1728 to 1584 when 144 first-rows were given
`aw = 0`, and on that fixture the **old** mean form reproduces the fixed form to the bit
(`O1` = `P1`, `mreldif(b)` 8.65e-16 both). `customweightmultiplier()` is different: it multiplies
the estimand weight without touching e(sample) (`e(N)` stays 1728), so a cwm of zero leaves the
row in the collapse with weight zero, first in line to be picked. Since every other in-sample row
carries `w · cwm > 0` for the first existing estimand weight (1426-1457), **proviso (ii) is
reachable only through `customweightmultiplier()`**. It is still the right way to write the code -
the fixed form costs nothing and is correct on both routes - but its severity should be read
against the author's standing position that cwm is caveat emptor (section 7). The current
prototype's form (Σ(w·x)/Σw onto every row, then the pick) passes the two of these configs that
round 2 re-ran on it verbatim - `A1r_w0f_pw_cwm` and `A1r_w0last_cwm` carry exactly the prepcode
and options of `V1bf_w0f_pw_cwm` and `V1bf_w0last_cwm` - and the six A1 variants of the same
shape.

### 2.3 Proviso (iii) measured: why `i` has to be in the key, and where it bites

`_i`, wooldid's internal cohort index, is the user's `i` on every row inside the estimation
window and **0 elsewhere** (`wooldid.ado` 628-629, with an identical twin at 752-753 on the
`esfixedbaseperiod` branch: `gen double _i = i if _t > 0 & _t < . & <the user's if>`, then
`replace _i = 0 if _i == . & <the user's if>`). So `i` in the key changes exactly one thing: the reference-period
rows, which otherwise pool across every cohort that shares an onset (`ttre` is in the key, so
cohorts 1-3, 4-6, 7-9, 10-12 of `continuous_control_dense` pool in pairs of three).

V1a(a) ran the fourteen pre-registered configs on a variant with that one line removed. All
fourteen come back PASS from the harness - and that is the harness's blind spot, not a clean
bill of health: the differ never lets an `e()` scalar's value reach a verdict. Reading the `%21x` dumps
directly, three of the four pre-registered predicted-MISS configs move a number wooldid prints:

- `V1a_noi_all_uzv_dyex_if` (no user weights at all, just an `if` that thins cohorts 1 and 4):
  the printed "p-value for joint test that all pre-treatment ame_it horizon estimates are 0"
  reads **0.1112** collapsed against **0.0956** on the full data - across the conventional 0.10
  line.
- The smallest breaking configuration is `V1a_noi_attit_uzv_dyex_wbig`
  (`att att_it contreat(dose_vary) contreatcontrols(untreatedZero) contreatelasticitytype(dyex)
  cluster(person_id) esprelength(2) espostlength(2) jointtests`, `aw=wbig`):
  `e(joint_Choriz_i)` = 0.0031356104 vs 0.0031085516, 8.6e-3 relative.
- The closest **non**-breaker is the same call with `att_i att_itime` in place of `att att_it`:
  every `e()` scalar bit-identical.

That minimal pair names the mechanism. The worker builds the cell mean under the *first
existing* estimand weight - `w_att` when `att` is requested - while the pre-trend joint-test
read (prototype 2675) weights by `w_att_it` over a subpop that includes the `_t==0` reference
bucket. Inside a genuine (i,t) cell the two weights are proportional and the mean is exact;
across pooled cohorts they are not. Drop `att` and `w_att_it` becomes the weight the mean is
built with, and the divergence vanishes. Attribution is closed at the other end too: re-running
the three breaking configs on the prototype (`i` in the key) reproduces the full-data
`e(joint_Choriz_i)` **bit for bit** in all three, with worst scalar residue 4.5e-12.

Two riders. The damage is confined to the whole-sample pre-treatment joint-test p-values
(`e(joint_Choriz_i)`, `e(joint_horiz_i)`): across all fourteen variant configs the main results
table moves by at most 5.1e-13 relative and `e(b)`/`e(V)` by 5.4e-14/6.9e-15, so this is not
"the collapse changes the estimates". And the ingredient list is narrower than "contreat": the
breakage needs `contreatcontrols(untreatedZero)` (so the reference-period dose slope is
estimated rather than absorbed) *plus* a read whose weight differs from the mean-building
weight; the `basic`-controls twin passes at 2.9e-12 and the `ccnoabsorb` twin at 1.5e-15. It is
not purely a `dyex` phenomenon either - under `untreatedZero` the dose enters as
`c.dose#i._i#i._t#i.targetflag`, so the att-side derivative is level-sensitive too and
`e(joint_horiz_i)` moved 5.4e-4 on a plain-`dydx` config.

Thirteen of the fourteen pre-registered predictions were right. The fourteenth,
`V1a_noi_sg3_uzv_dyex_wbig`, is a mechanism-confirming null: `sg_three` is a cohort surrogate on
this fixture (cohorts 1,4,7,10 -> SG A, 2,5,8,11 -> SG B, 3,6,9,12 -> SG C) and `subgroups` is
already in the collapse key, so it silently reinstates the proviso under test.

### 2.4 Proviso (iv) measured: the recast to double

V1a(c) deleted the `recast double` line and ran seventeen configs over two fixtures - twelve
with a float dose (`dta_cont`'s native float `intensity`, and `continuous_control_dense` doses
forced to float) and storage-type twins holding the identical numeric values in double. Nine of
the twelve float configs MISS; every double twin PASSES at 5.7e-14..1.7e-11. The smallest
breaker is `B_dydx_f_basic`: `dta_cont`, `att att_it contreat(intensity)
contreatcontrols(basic)`, no weights, no cluster, no subgroups, default `dydx` - 70 cells beyond
1e-9, worst 5.6e-7 relative on an `ES_att:pre1` t statistic. Within that same run every
AME-family row holds at last-place noise (<= 2.5e-16) while every dose-scaled ATT row moves, up
to 5.6e-7 - the same level-sensitivity mechanism as section 2.1. The mechanism control is decisive: a float dose that
is *constant within every cell* (`A_dyex_f_cellconst`) is bit-identical over 840 cells, because
the double mean then rounds back to exactly the same float - float storage alone is not enough,
the dose has to vary inside the cell.

The magnitudes match float theory: the measured primitive damage (a double cell mean of a float
dose rounded back to float) is 2.4e-8 median, 5.6e-8 max relative, and the observed `e(b)`/`e(V)`
moves land at 1e-8..4e-8, reaching 5.6e-7 only on near-zero pre-period rows.

**What the recast is and is not doing** (measured 2026-09-02, `ct/Q/qF_recast.do|log`). It does
not change how the mean is *computed* - that is double either way (`gen double` for `w*x`, a
`gcollapse (sum)`, then a double division). It changes how the result is *stored back onto the
row*. Stata does not refuse or truncate an out-of-type `replace`; it silently **promotes the
storage type**, and the promotion stops at float for a `byte`, `int` or `float` variable (measured:
`replace` of 3.4 into a `byte` gives 3.4000001, the float value; into a `long` it gives 3.4
exactly, the double value). So without the recast a `byte`- or `float`-stored dose carries its cell
mean at float precision. Re-measured on one fixture with the same request and three storage types
for the same dose values:

| dose storage | with the recast (shipped) | without it |
|---|---|---|
| `byte` | `mreldif(e(b))` 8.6e-16, table 1.8e-14 | `mreldif(e(b))` 3.5e-09, table 6.1e-08 |
| `float` | 8.6e-16, table 1.5e-14 | 1.4e-09, table 2.7e-08 |
| `double` | 8.6e-16 | 8.6e-16 (the recast is a no-op) |

So it is a rounding question, but at *float* ulp rather than double ulp - seven orders of
magnitude, which is the difference between "no digit wooldid displays is affected" and a reported
`main:att` reading 1.381177483710 against the full data's 1.381177483150. Two things follow beyond
the magnitude. The help file's documented residue bound (about one ulp of the linear prediction,
variances at ~1e-16 relative) would simply be false for contreat runs without it. And on the
subgroup route the move changes a comparison qualitatively rather than numerically, because twelve
per-subgroup joint p-values are stored as e() *macros*: a residue-sized numeric change becomes a
string difference. What the recast cannot do is touch the regression - the collapse runs after
estimation, inside the preserve - so coefficients, the variance matrix and the clustering are
unaffected either way. **The honest
qualification** - which the angle insisted on - is that across all seventeen configs and all
24 result matrices carrying a p column there were **zero** p-value crossings of 0.01/0.05/0.10
and **zero** CI sign flips; the worst absolute p move anywhere was 1.66e-7. So the case for the
recast is a residue-budget and bit-reproducibility case, plus one structural consequence: on the
subgroup route (`A_dyex_f_sg3_wbig2`) twelve per-subgroup joint p-values are stored as e()
*macros* (strings), so a residue-sized move there registers as a structural difference
(`joint_Choriz_i_s1` 3.76145473265e-08 vs 3.76145875442e-08) - and the double twin has none.
(The scalar `joint_pre_att_it_allsgssame` moves too, 0.0001533690374492 vs 0.0001533685243506,
but a scalar sits below the verdict rather than counting as structural.) Round 2 closes the sufficiency side that this angle could only reason about: the
prototype *with* the recast ran sixteen distinct `dta_cont` configs on the natively-float
`intensity` (A3, A4, A7) and agreed with the full data to <= 1.1e-10 (worst `A7_dc_if`,
`ES_att:pre2` = 3.23435899e-05; the separately diagnosed `ctA4_unit_sg2_jt` residue MISS sits at
3.6e-9).

### 2.5 Proviso (v) measured: the pre-mean form the shipped gate was written against

V1a(d) ran the collapse with the dose left alone - the surviving row keeps its own draw. This is
not the author's proposal (which is the cell mean); it is the form the collapse had before the
proposal, and the one the shipped gate's "coefficients shifted, one sign flip" was measured on.
Ten paired configs on three datasets,
one Stata at a time. The binary control PASSES (5.2e-13), and **all nine contreat configs
MISS**: with and without weights, with and without cluster / event study / joint tests, OLS and
poisson, `dydx` and `dyex`. The smallest breaking configuration is
`att contreat(dose_vary) contreatcontrols(basic)` and nothing else: the reported ATT is
1.3775620504 against 1.3836141895, its standard error is overstated 5.6% (0.0145868807 vs
0.0138113285) and `e(V)[1,1]` by 11.5%. The error is not conservative - elsewhere the first-row
form *understates* the standard error. In reported terms it is decision-changing: a p-value
crosses 0.05 in five configs at eight row-level sites, in both directions (a fabricated
pre-trend, `ES_att:pre1` p 0.0022 vs 0.1243; and a hidden one, poisson `pretreatment:att`
p 0.1569 vs 0.0093). The shipped text's "one sign flip" is an understatement: on `dta_cont`
six `e(b)` estimates flip sign, with one t statistic going from -0.4511 to +0.3920.

The only reads that survive the first-row form are `dydx(contreat)`/`ame` under `dydx`
(<= 9.9e-15 across five configs), which are cell-constant and so do not see the dose level at
all - and they break too the moment `contreatelasticitytype(dyex)` makes them level-linear
(0.68-0.78 relative). Since wooldid always emits the `dydx(targetflag)` family alongside them,
and that family is δ_g + β_g·x, every contreat run is corrupted somewhere. The closing argument
is arithmetic rather than empirical: the first-row draw equals the cell mean only if the dose
is constant within the cell - in which case β_g·D·x is collinear with δ_g·D and β_g is not
identified - or if every cell is a singleton, in which case the cell-saturated model is exactly
fitted. There is no non-degenerate contreat configuration in which the pre-mean form is exact.

## 3. What the gates keep out (measured, gate-off variant)

V1b(e) removed the eyex/eydx, `contreatpoly(>1)` and contreat-in-controls gates from the
prototype and ran twelve configs so the collapse engaged on them (`adv2/ctV1b_e/RESULT.txt`;
`continuous_control_dense`, robust or `cluster(person_id)`, `esprelength(2) espostlength(2)`,
worst relative move in the arm-C-vs-U comparison):

| config | form | worst e(b) | worst e(V) | verdict |
|---|---|---|---|---|
| `V1be_ctl_dx1` | `controls(c.dose_vary#c.x1)` | 0.0063 | 0.0086 | MISS, 110 cells - gate needed |
| `V1be_ctl_dx1_ccnoabs` | same with `ccnoabsorb` | 0.0063 | 0.0086 | MISS, 110 cells |
| `V1be_ctl_isgp` | `controls(i.sgp#c.dose_vary)`, `sgp = mod(person_id,3)+1` (varies *within* a cell) | 0.0017 | `V_es_ame` 0.083 | MISS, 116 cells |
| `V1be_eyex` | `contreatelasticitytype(eyex)` | 0.00029 | 0.0011 | MISS, 130 cells |
| `V1be_eydx` | `contreatelasticitytype(eydx)` | 1.1e-5 | 6.8e-5 | MISS, 130 cells |
| `V1be_pois_eyex` | poisson + `eyex` | 0.0029 | 0.0042 | MISS, 124 cells |
| `V1be_poly2` | `contreatpoly(2)` | 0.35 | 0.27 | MISS, 130 cells |
| `V1be_poly3` | `contreatpoly(3)` | 0.84 | - | MISS, 48 cells |
| `V1be_ctl_isg3` | `controls(i.sg_three#c.dose_vary)`, `sg_three` constant within a cohort | 2.5e-16 | 9.4e-16 | PASS - the gate is conservative here |
| `V1be_ctl_plain` / `_cplain` | `controls(dose_vary)` / `controls(c.dose_vary)` (plain extra slope) | 5.4e-14 | 1.6e-15 | PASS - conservative |
| `V1be_dyex_ctl` | `contreatelasticitytype(dyex)` (the allowed elasticity form) | 2.3e-16 | 1.4e-15 | PASS (control) |

So every refused form is refused for cause, at magnitudes from 1e-5 to 0.8 relative, and the
rows separate the mechanism cleanly.

*Which reads break.* Under contreat-in-controls it is only the **ame** family: in
`V1be_ctl_dx1` the worst att-family `e(b)` move is 9.8e-15 while `ES_ame_pre1` moves
1.479653786 -> 1.489015445 (0.63%); in `V1be_ctl_isgp` the att family holds at 2.4e-14 while
`ES_ame_it_pre1` moves 0.3282987778 -> 0.3288642725. That is exactly what the section 1 table
predicts: `dydx(contreat)` picks up the control term's slope, which varies inside the cell,
whereas `dydx(targetflag)` does not. Since wooldid always reports the ame rows when `contreat()`
is given, the gate has to refuse the whole run, but it is worth knowing that the att side is
never the casualty here.

*Which interactions break.* `i.sgp#c.dose_vary` (MISS) against `i.sg_three#c.dose_vary` (PASS)
isolates it: the gate is needed when the variable interacted with the dose varies *inside* a
cell, and `sg_three` does not (each cohort of `continuous_control_dense` sits in exactly one
`sg_three` level, verified by crosstab, so the interaction is cell-constant). The token gate
cannot tell those apart, so it refuses both, and it likewise refuses a dose entered as a
stand-alone linear control, which is also exact. In the other direction the gate matches the
*name* of the contreat variable, so a copy of the dose inside a control term
(`A2n_ctl_dosecopy_x1` with `controls(c.x1#c.dose_copy)`, `A2n_ctl_dosecopy_t` with
`controls(c.dose_copy#i.t)`) is not gated - and both PASS, `e(b)` agreeing to 2.9e-14, because
margins does not differentiate through the copy. Both asymmetries are acceptable: the first
costs the speed-up on unusual specifications, the second is exact anyway.

Gate audit (A5, 42 configs): all 24 gated configs report `fcol C=0` with C bit-identical to U;
every token in the refusal list fired at least once. Stand-alone explicit requests:
`contreatpoly(2) fastcollapse` -> r(198) (also with `safety(off)` - no escape);
`contreat(dose_vary) fastcollapse` (poly 1) -> rc 0, `e(wdid_fcol)=1`; `eyex fastcollapse` ->
r(198); `dyex fastcollapse` -> rc 0, fcol 1; `lattice(all 2) fastcollapse` -> r(198);
`poisson poisexpresults fastcollapse` -> r(198); `poisson dyex fastcollapse` -> rc 0, fcol 1.
A3 (61 configs) confirms that none of the event-study/diagnostic options (`jointtests`,
`adversarial`, `altmethod_adversarial`, `esfixedbaseperiod`, `esrelativeto`, `pretrendpower`,
`histogramcohorteffects`, `summarystats`, `makeplots`/`saveplots`, `subgroups`,
`suppresssgprimaryeffects`, `bootstrap`, `permute`, `fastnose`, `contreatwithin`) needs a gate:
the collapse engaged on all of them and all PASS. `gteffects` is refused by wooldid itself with
contreat (r(198), both arms).

## 4. The probe: 328 paired runs on the fixed prototype

Harness as in the subgroups assessment: arm C runs the prototype (collapse engages where the
gate allows), arm U the same config with `fullprecisionslow`; every `e()` scalar, macro and
matrix cell is dumped in `%21x` and compared (structure: rc, dims, names, present-vs-missing,
Cautions; cells: a move that is both > 1e-9 relative and > 1e-12 absolute is a MISS -
`differ.py:326-327`). A first round over the same seven
angles (290 distinct configs, 276 of which produced dumps) found the zero-weight-first-row defect
of section 2.2 and the
contreat-in-controls gate gap; both were fixed and every angle was then re-run from scratch.
Those re-runs are the tally below (config counts are distinct config names; where a config was
re-run after a config-file mistake, only its final verdict is counted):

| angle | configs | PASS | MISS | NOTENGAGED | BOTH_ERR | worst PASS residue (relative) |
|---|---|---|---|---|---|---|
| A1 estimands x weights x dose storage | 37 | 33 | 2 (off-domain, section 5) | 0 | 2 | 1.87e-10 on the `ES_att_it:pre2` p-value 1.216355047e-08 (abs 2.3e-18), in the config with a 1e2 dose offset and cwm zeros; `e(V)` <= 2.4e-11 |
| A2 contreatcontrols / controls / FE / timing options | 48 | 43 | 0 | 3 (contreat-in-controls, expected) | 2 (wooldid's own refusals) | 2.16e-11 on the `main:ame` p-value 7.409846658e-07 (abs 1.6e-17) |
| A3 event-study, pre-trend, diagnostic reads | 61 | 56 | 0 | 0 | 5 (config errors / documented refusals) | 1.65e-11 on `ES_att_pre1` = 2.15100705e-04 on the float unit-level fixture (abs 3.6e-15) |
| A4 subgroups | 44 | 39 + 3 FLAG_SEFAIL | 2 (both residue, diagnosed) | 0 | 0 | 3.95e-11 on `ES_att_SG2_pre4` = 8.994515445e-05 (abs 3.6e-15); a 2.6e-82 pairwise p-value moves 6.3e-11 relative (1.6e-92 absolute) |
| A5 gates and nonlinear forms | 42 | 13 | 5 (all poisson+dyex, section 5) | 24 (expected) | 0 | 4.5e-12 on the `pretreatment:att` CI bound -9.89540101e-05 (abs 4.5e-16) |
| A6 VCE and resampling | 48 | 39 + 2 FLAG_SEFAIL | 0 | 0 | 7 (config errors / documented refusals) | 3.66e-11 on `bias_boot` ES_att_pre2 = -4.84683241e-05 (abs 1.8e-15) |
| A7 data shapes and edge cells | 48 | 48 | 0 | 0 | 0 (one config-file mistake errored in both arms, PASS on re-run) | 4.15e-10 on the `ES_att:pre2` p-value 1.054289332e-11 under a deliberate +100 dose offset (abs 4.4e-21); on cells of ordinary size, 1.1e-10 on `ES_att_pre2` = 3.234358994e-05 |

Three scope facts belong with that table. The A1 row counts the 37 configs run on the prototype;
a 38th, `A1r_neg_wsumneg2_sg`, ran on the signed-guard variant and PASSes (section 5.1). Every
paired `%21x` comparison ran on a small fixture - the largest `e(N)` anywhere is 3,456, the median
1,728, and 287 of the configs are `continuous_control_dense` - so the evidence that the collapse
scales is the timing runs of section 6 (`mreldif(b)` <= 4.3e-15 at 145,335 rows), not this table.
And every poisson config in the probe took the ppmlhdfe fallback on this fixture (100 "initial
estimation of ppmlhdfe failed" hits across the round-2 C-arm logs), so the primary absorbing
ppmlhdfe path with contreat is unmeasured; the ringer average is a no-op there by construction,
but that is reasoning, not measurement.

FLAG_SEFAIL rows are the route-P subgroup degeneracy screen firing on cluster-at-cohort /
cluster-at-unit runs that wooldid itself reports as SE-degenerate. Two of the three take the `fallback` branch: `ctA4_relyr_sgr` is bit-identical to the full run
and `ctA4_sg3_cohcl_nojt` agrees to residue (332 of its 10,448 cells move, worst 4.6e-12
relative, none beyond 1e-9). The third, `ctA4_sg3_cohcl_safoff_jt`, is `safety(off)`, so by design the screen
does **not** fall back: it proceeds with a caution, and there the collapsed run does differ -
`e(pwsgs)` rows 2 and 8 (the `sg_three` pairwise contrasts) are missing under the collapse
where the full run reports them, and 12 cells move beyond 1e-9, the worst being
`e(pwsgs)[10,2]` = 7.958926556e-09 vs 7.965734878e-09 (8.6e-4 relative on a standard error of
8e-9). That is the documented `safety(off)` contract from the subgroups assessment, not a new
finding: the contrast is numerically singular (the full run reports those rows with an SE of
1.6e-10 and a t of 5.2e+08) and wooldid prints its "standard errors were not identified" caution
in both arms.
Apart from that one run, rc, cautions and notes were textually identical between arms in every
config of every angle (scripted comparison), and the only other structural difference anywhere
in round 2 is the `A1r_w0h_cwm_pw` hopeless config, where both arms error with different
return codes.

Aggregating the dumps directly rather than trusting the per-run summaries - which matters,
because the harness's differ never lets an `e()` scalar's *value* reach a verdict (only its
presence or absence), screening scalar values at 1e-6 against a denominator floored at 1e-12
(`differ.py:300-302`), so a PASS verdict alone does not mean "no reported number moved" (the
V1a(a) probe of section 2.3 exploits exactly that gap):
all 328 configs dumped in both arms; ten are accounted for individually elsewhere on this page
(the two signed-cwm and five poisson+dyex configs of section 5, the two diagnosed A4 MISSes above
and the `safety(off)` caution run) and thirteen more returned an empty `e()` because both arms
refused the request, leaving **305** configs compared object by object in `%21x`:
**846,157 matrix cells, 8,881 `e()` scalars and every `e()` macro**. Of the matrix cells 73,503 (8.7%) differ in
their last bits and **none** differs by more than 1e-9 relative once the cell's own magnitude
is 1e-12 or larger; 40 configs are bit-identical cell for cell. Of the scalars 460 differ and
only four move more than 1e-9 relative - all four are joint-test p-values of magnitude 4.9e-38
or 2.1e-21 (worst 1.4e-8 relative on 4.87317690942e-38), which print as 0.0000 in either arm.
The string-valued per-subgroup joint p-values (`e(joint_horiz_s*)` and friends) differ in 47
places, all at <= 7.2e-12 relative; the largest absolute move is 4.7e-13, on a p-value of 0.374,
so no printed digit changes. No other macro differs. The extreme values quoted in the table
above are all near-zero cells: the largest relative move on a cell of magnitude >= 1e-8 is
1.87e-10 (`A1r_d_big2_cwm0f`, an `ES_att_it:pre2` p-value of 1.216355047e-08, which prints the
same ten digits in both arms), and the largest absolute move anywhere is 4.77e-10, on the
t statistic 39.249807005 of the same config's `ES_att:post0` row - the tenth significant digit
of a t statistic, under a deliberate 1e2 dose offset chosen to maximise cancellation.

The two A4 MISSes were re-run by independent diagnosers and are harness-threshold residue:

- `ctA4_unit_sg2_jt`: a binary ES pre1 coefficient of -1.97e-6 for subgroup 2 on the
  float-precision `dta_cont` fixture moves 7.1e-15 absolute (two 2^-48 quanta; 3.6e-9
  relative only because the cell is 2e-6). The read is the per-subgroup horizon margins
  (shipped ~2760); the full-data arm itself varies by up to 7.1e-15 at this cell under row
  permutations (two shuffles give C == U bit-identical, reverse sort gives one quantum).
- `ctA4_pois_sg3_jt`: the `pre_ame_it:2vs1.sg_three` pairwise contrast (shipped 3361) is an
  exact-zero contrast (both subgroup pre-period AMEs are -.0305073 under the ppmlhdfe
  fallback), se 2e-9 is cancellation residue, and the C-U move of 3.5e-18 absolute is 1.7e-9
  of that se; the U arm moves 2e-15 at the cell under row order alone. Not a branch flip:
  r(p) present in both arms, pairwise rc identical, se bit-identical.

Across the 44 subgroup configs the contreat subgroup section - which the route-P degeneracy
screen does not cover - never flipped a branch: `r(p)` present-vs-missing, `pwcompare` rc, the
rank check and the `o.` markers agreed between arms in every config except the `safety(off)`
caution run above, where the screen had already fired and wooldid reports the SEs as
unidentified.

Resampling (A6): bootstrap `boot_completed`, `boot_z0`, `boot_navail` and the draw sequence
are bit-identical - the collapse changes no draw - while the quantities built from the
replications carry the collapse's residue through 30-200 reps: worst over the 48 configs,
`bias_boot` 3.7e-11 (on a 4.8e-5 cell), `V_boot` 2.6e-11, `wdidsgresults_*` 5.4e-12,
`bootdist` 1.9e-12, `boot_ci_normal` 1.3e-12, `boot_ci_bc`/`_percentile` 6.7e-13,
`V_analytic`/`V_es_analytic_*` <= 8.2e-15. On the permute side the counts and accept/reject
tallies are identical and `permdist` <= 1.2e-14, `perm_comparator` <= 1.2e-15, across
`permute(shuffle)`, `permute(probabilities)`, `permstrata()`, `permtstat`, `permtci` and the
subgroup form. (`permexact` appears in no contreat config and is therefore *not* covered by this
evidence; `permtci` is - `A6_pm_tci` PASSes.) Under the subgroup fallback the reps follow the
observed run's fallback.

## 5. Where the arms disagree

### 5.1 Signed `customweightmultiplier()` under `safety(off)` (A1, two configs)

With cwm = -1 on some rows of one cell, the cell's estimand-weight sum can be negative or
exactly zero:

- Σw < 0 (`A1r_neg_wsumneg2`, 9 rows at -1 and 3 at +1 in one cell): the weighted mean is well
  defined but the prototype writes it only `if wsum > 0`, so the pick row keeps its raw dose.
  The main att family moves 5.5e-4 relative (1.38099003327 against 1.38174568463), `e(V)` up to
  1.2e-2, and the worst reported cell is the `ES_att:post1` standard error, 0.0312441770 against
  0.0314276695. Changing the guard to
  `wsum != 0` makes the config PASS (b residue 5.4e-14, V 1.3e-15; variant
  `A1/wooldid_ct_signedguard.ado`). The guard change is a **no-op everywhere else**: the 37
  A1 configs re-run on the signed-guard variant (`adv2/ctV4_a..d`) reproduce 36 of 37 census
  lines *and* verdicts byte-identically - same cell counts, same worst residues, same
  matrices, including both BOTH_ERR configs and the Σw = 0 MISS below - and the 37th is
  `A1r_neg_wsumneg2` going MISS -> PASS. This is what a free fix looks like: `wsum > 0` and
  `wsum != 0` can differ only where the sum is negative, which requires `safety(off)`.
- Σw = 0 with non-zero rows (`A1r_neg_wsum0`, 6 rows at +1 and 6 at -1 in cell cohort 5 / t 8):
  the full-data read contributes b_g·Σ w_r x_r ≠ 0 from that cell, which no (summed weight,
  weighted mean) pair can represent - the cell instead enters every read with weight zero. The
  main att family moves 2.7e-4 relative (1.38252203074 against 1.38289179654), `e(V)` 1.7e-3,
  `ES_att_post1` 1.24123837402 against 1.24368807244, and the worst reported cell is that
  horizon's CI lower bound, 1.18107330675 against 1.18347063323. The ame (slope) reads stay at
  residue throughout (worst 1.7e-16 on `e(b_es_ame*)`/`e(V_es_ame*)`, 2.1e-15 on the ame rows of
  the reported table, in both signed-weight configs), as the mechanism predicts: only the reads that carry the dose *level* can see a mis-set pick-row dose. This is an inherent limit:
  the only exact treatment is to refuse to collapse (fcol=0 or the route-P fallback) when any
  cell has Σw = 0 with a non-zero row. The predicate must be Σw = 0 **and some row of the cell
  carries a non-zero weight**: a cell whose rows are all zero-weight is legitimate and common
  (a cohort or a horizon zeroed by `customweightmultiplier()`), contributes nothing to either
  arm, and must still collapse - `A1r_w0cohort_cwm`, `A1r_w0allbut1_cwm`, `A1r_w0f_pw_cwm` and
  `A1r_w0last_cwm` all engage and PASS today and a bare `wsum == 0` test would refuse them
  all. Both directions were built and measured (V4(z), variant `ct/V4z/wooldid_ct_zerodetect.ado`,
  run `adv2/ctV4_z`): the predicate `wsum == 0 & wsrc != 0 & !missing(wsrc)` is **15 executable
  lines**, and with it `A1r_neg_wsum0` comes back `fcol C=0` with all 1392 dumped matrix cells
  bit-identical between arms and zero structural differences, while `A1r_neg_wsum0_pw` (the same
  signed pattern with `pw` breaking the cancellation) stays engaged and PASSes at 2.4e-13 and the
  four legitimate all-zero-weight-cell configs reproduce their pre-detector census lines byte for
  byte. Cost, measured at 2,000,000 rows and 20,000 cells: one extra `count if` pass at 0.015 s
  against 0.124 s for the `gcollapse` beside it - about 8% of the contreat-mean block, invisible
  next to the margins calls the collapse exists to accelerate - and it reuses the `wsum` the mean
  already materialises, adding no sort and no second collapse. Two limits worth writing into the
  code comment: the test is exact-zero only (weights that cancel to +-1e-17 rather than to 0
  still slip through, and the signed guard then divides by that residue), and the refusal reason
  reaches only the verbose line, not `e()`, unless a production version surfaces it in
  `e(wdid_fcskip)`.

Neither is reachable on the documented domain: wooldid refuses a negative cwm without
`safety(off)`, and with non-negative weights Σw = 0 means every row of the cell has weight
zero, so the cell contributes nothing in either arm (`A1r_w0cohort_cwm`, `A7` zero-weight
cells: PASS). A7's `cwm = -1 on one row per cell` under `safety(off)` keeps every cell sum
positive and PASSES. Per the pathological-data doctrine this is a document-or-guard decision,
not a build item: the one-token guard is free, and the Σw = 0 refusal is a seven-line count
inside the worker plus an eight-line restore path in the caller - the fifteen lines measured
above.

### 5.2 `poisson` + `contreatelasticitytype(dyex)`

(A5, 5 of 5 variants: plain, `dose_common`,
if()+cluster, jointtests+att_it, pw). Only the ame (dose) reads move, and only through their
variance: `e(V)` up to 3.5e-6 relative (`pretreatment_ame` 5.484551479e-06 vs 5.484570661e-06),
`e(V_es_ame*)` up to 3.2e-6, and in the reported table a `pretreatment:ame` CI bound of
-2.522112643e-05 vs -2.522825531e-05 - a 2.8e-4 relative move on a reported number. The point estimates barely move (`e(b)` <= 6.4e-10 relative, `e(b_es_ame)`
<= 4.2e-10), and the binary att rows computed in the same runs agree to 1.5e-13
(`e(V_es_att)` to 2.3e-15).
Mechanism: under poisson without `poisexpresults` the dyex read is `margins, dyex(x)
expression(predict(xb))`, and margins differentiates an `expression()` numerically (in x for
the estimate, in the coefficients for the Jacobian), so the round-off pattern of the numerical
derivative changes when a within-cell-varying x is replaced by its mean. Calibration
(`ct/A5/A5_calib2`, same config, same Stata): merely **permuting the rows of the full data**
moves the same cells by `e(V)` 1.28e-6, `e(b)` 4.2e-10 and a reported table cell 1.08e-5,
against the collapse's 3.05e-6 / 4.2e-10 / 1.25e-5 - so the collapse perturbs this read by the
same order of magnitude as the input row order does (at most 2.4x it), and with x constant
within cells the collapsed numeric read is bit-identical to the full one. The regress dyex
path (analytic chain rule) and the poisson dydx/att reads (`predict(xb)`, analytic) PASS at
<= 4.5e-12 (`e(b)` <= 1.5e-13). Options: gate `poisson` + `dyex` (one more token; costs the speed-up on a rare
combination) or document the residue - up to the fourth significant digit of a near-zero CI
bound, which is more than the help file's own contract allows: `wooldid.sthlp` 291-300 says
results "agree to machine precision but are not bit-identical" and quantifies the variance side as
"about 1e-16 relative regardless of the scale of the predictions". The per-rep cost under bootstrap/permute is unaffected either way.

### 5.3 A third disagreement, in the opposite direction: extreme within-cell cancellation

The refutation angle (V1b(r)) was asked to break the weighted-mean claim on the documented
domain with attacks the seven probe angles could not run. Its strongest attack was a
catastrophic-cancellation ladder: a dose `1 + S·(z - z̄_cell)` whose within-cell spread is S
times its cell mean, unweighted, with the estimand made analytically S-invariant by
construction, run at S = 1, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e12 (runs `adv2/ctV1b_r` and
`adv2/ctV1b_r2`, 19 configs, plus row-order twins at 1e6, 1e8 and 1e10 and two controls). The
harness reports eleven MISSes from S = 1e4 up (the eight rungs at or above 1e4, plus the three
row-order twins). **The direction of fault is the finding**: it is the
*uncollapsed* arm that moves.

| S (spread/mean) | `main:att` estimate C vs U | collapsed `main:att` se | full-data `main:att` se | collapsed `main:ame` se | full-data `main:ame` se |
|---|---|---|---|---|---|
| 1 | 1.6e-16 | 0.1302380935 | 0.1302380935 | 0.183013 | 0.183013 |
| 1e4 | 6.5e-14 | 0.1302380935 | 0.1302380872 | 1.83013e-05 | 1.83013e-05 |
| 1e5 | 6.1e-13 | 0.1302380935 | 0.1302379502 | 1.83013e-06 | 1.83013e-06 |
| 1e6 | 1.3e-11 | 0.1302380935 | 0.1302372095 | 1.83013e-07 | 0.00183469 |
| 1e7 | 5.7e-11 | 0.1302380935 | 0.1302673181 | 1.83013e-08 | 0.000393993 |
| 1e8 | 5.1e-10 | 0.1302380936 | 0.1306997719 | 1.83013e-09 | 0.30165 |
| 1e9 | 8.7e-09 | 0.1302380927 | 0.1187200712 | 1.83013e-10 | 5.4609 |
| 1e10 | 4.0e-08 | 0.1302380882 | 0.0824266096 | 1.83013e-11 | 1.40567 |
| 1e12 | 9.2e-06 | 0.1302396883 | 0.0674313401 | 1.83007e-13 | 4178.73 |

The estimand is S-invariant by construction (`main:att` = 0.3429364417 analytically), and the
dose AME is exactly 0.5 with an analytic standard error of order 1e-9 x 1/S. Read the table
across: the **collapsed** columns are the analytic answer - the att se holds 0.1302380935 to ten
digits from S = 1 through S = 1e10, and the ame se follows the exact 1/S law over twelve decades -
while the **full-data** columns drift away from it from S = 1e6 on.

The DGP fixes the dose AME at exactly 0.5 with an analytic standard error of order 1e-9. At
S = 1e8 the collapsed arm reports 0.499999997349 with CI [0.499999993731, 0.500000000966]; the
full-data arm reports 0.500013987294 with se 0.301649731716, p 0.0995908456619 and CI
[-0.096254683896, 1.09628265848] - insignificant at 5%. Four independent anchors put the fault
on the uncollapsed side. (1) The collapsed ame se follows the exact analytic 1/S law to six
significant figures through S = 1e10 (five at S = 1e12) while the uncollapsed one stops scaling
at S = 1e6. (2)
The collapsed att se holds its S-invariant value 0.1302380935 to ten digits through S = 1e10
while the uncollapsed one drifts to 0.0824266096. (3) **Row-order twins** - the same data in a
different order, nothing else changed, compared within one arm - move the *uncollapsed* arm at
S = 1e8 by 1.09e-2 on the att se, 0.213 on the ame se and 0.488 on the ame p (0.0995908457 ->
0.194434684), while the collapsed arm moves 5.7e-10, 1.8e-10 and 0 on the same three; at
S = 1e10 the uncollapsed arm moves the reported `ES_att:post0` p-value from 0.0615285935 to
0.572652779 - across 0.10 on row order alone - where the collapsed arm moves it 5.6e-06. (4)
Both arms agree to 1.7e-15 relative at S = 1 - last-place residue, with the att se
bit-identical - and the collapsed values are the exact analytic continuation of that agreed
point. The two structural differences in this
family are the collapsed arm alone printing wooldid's "reported standard errors appear to be
identical" caution - which fires precisely because the collapsed SEs are correct and nearly
deterministic while the uncollapsed ones are noise.

Controls close it: `V1br_magctl` gives the dose the same 1e10 *magnitude* with the cancellation
removed (cell mean 1e10, ratio 0.3) and agrees to 1.7e-15, so magnitude is innocent; the scale
probe (an ordinary dose in [-0.5, 0.5] with a cell mean of 1e-9, ratio 1e9) reproduces the
S = 1e9 rung, its `main:att` estimate moving 7.65e-09 against the rung's 8.70e-09, so the ratio
is the whole story; and in every rung the
pre-treatment rows, whose cells carry no dose channel in this DGP, stay exact at <= 1.7e-15.

Three consequences. First, nothing here argues for code: a cell whose signed doses cancel to four
or more digits is a fixture regime, not a data regime, and it sits under the pathological-data
effort cap. Second, and worth writing into the doctrine: **in this regime a battery config will
MISS without the collapse being at fault**, because the harness's "U is truth" assumption fails -
the full-data arm moves a reported p from 0.062 to 0.573 on row order alone. Third, the honest version of
the headline is not "the collapse is exact and the full run is broken": at ratio 1e4 the
full-data se error is 4.9e-8 relative, immaterial; it becomes material (>0.1% of a reported se)
only from ratio ~1e8. The angle's own two nulls belong here too: the first-existing-weight attack
found the four estimands moving 1.11e-09 / 8.69e-10 / 1.22e-09 / 1.25e-09 - a 1.44x spread
against a 10x refutation criterion, so proviso (i)'s weight selection survives its sharpest test;
and the joint-test amplification attack is *inconclusive by construction* rather than passed,
because any fixture that amplifies the dose channel also makes the dose slope near-deterministic,
so the pre-period cells the joint tests read carry no dose channel (the four joint scalars are
bit-identical at every cancellation rung, S = 1e8 and S = 1e10 included, where reported SEs move
by tens of percent; at the S = 1 control three of the four differ in their last bits).

## 6. The prize

Stand-alone timing (`ct/V3/v3_timing.log`, `v3_scale.log`): one Stata, `timer` around each
call, collapsed run first then `fullprecisionslow`, medians of three; `y2 cohid year ttre` on
`dta_s_admin` with a synthetic dose (`set seed 4242; dose = 0.5 + runiform();
y2 = y + 0.3*dose*post`), `cluster(state) esprelength(2) espostlength(2)`, contreat runs with
`contreat(dose) contreatcontrols(basic)`. The admin sub-sample keeps one fips in seven; the
x5 sample is `admin_full` expanded five times with jitter. The one row that is not on
`dta_s_admin` is the `continuous_control_dense` run, which uses `cluster(person_id)` and a
3-pre/3-post window:

| data (`e(N)`) | run | collapsed | full | ratio | mreldif(b) |
|---|---|---|---|---|---|
| admin sub-sample (4,162) | att, contreat | 0.61 s | 0.78 s | 1.27 | 1.7e-15 |
| | att att_it, contreat | 0.94 | 1.28 | 1.35 | 1.7e-15 |
| | att att_it, contreat, subgroups(region), cohreg | 13.45 | 17.94 | 1.33 | 1.7e-15 |
| | att att_it, contreat, subgroups(region), safety(off) | 3.79 | 5.35 | 1.41 | 1.8e-15 |
| | att, contreat, adversarial jointtests safety(off) | 0.78 | 1.12 | 1.44 | 1.7e-15 |
| | att att_it, contreat, histogramcohorteffects | 1.75 | 2.30 | 1.32 | 1.7e-15 |
| | att att_it, **binary** (comparison) | 0.40 | 0.57 | 1.43 | 1.7e-15 |
| continuous_control_dense (1,728) | all four estimands, contreat, jointtests adversarial, cluster(person_id) | 4.98 | 6.56 | 1.32 | 8.6e-16 |
| admin full (29,067) | att att_it, contreat | 1.76 | 4.83 | 2.75 | 1.8e-15 |
| | att att_it, binary | 0.70 | 1.56 | 2.23 | 1.8e-15 |
| | att att_it, contreat, subgroups(region), safety(off) | 4.90 | 20.76 | 4.24 | 1.8e-15 |
| | att, contreat, adversarial jointtests safety(off) | 1.74 | 4.05 | 2.32 | 9.7e-16 |
| admin x5 (145,335) | att att_it, contreat | 6.51 | 26.93 | 4.14 | 4.3e-15 |
| | att att_it, binary | 2.23 | 8.13 | 3.65 | 1.8e-15 |

`mreldif(V)` <= 3.6e-16 throughout; `e(wdid_fcol)` 1/0, guard empty in every pair. The
contreat prize tracks the binary prize at every size (the regression dominates small samples;
the margins window dominates large ones), and the extra worker cost (one `gcollapse (sum)` of
two variables per collapse) is invisible at this scale.

## 6a. Two questions the author raised on reading this (2026-09-02)

**Does the binary collapse survive unequal cell sizes and per-observation weights?** Yes, and it
was re-measured rather than assumed (`ct/Q/qA_binary.do|log`, the shipped ado, ten paired
configs on a deliberately unbalanced `continuous_control_dense`: 144 cells holding 3 to 12
observations each, `aw` cell totals ranging 3.50 to 40.88, and a lognormal weight whose cell
totals range 0.27 to 12,065). Every case agrees to residue - `mreldif(e(b))` <= 1.7e-15,
`mreldif(e(V))` <= 8.7e-19, `mreldif(e(wdidmainresults))` <= 3.9e-14:

| case | mreldif(b) |
|---|---|
| unequal cell sizes, no weights | 8.2e-16 |
| balanced cells, per-observation `aw` | 8.7e-16 |
| unequal cells + per-observation `aw` | 8.7e-16 |
| unequal cells + lognormal `aw` (44,000x spread in cell totals) | 8.9e-16 |
| unequal cells + `pw` | 8.7e-16 |
| unequal cells + weight constant *within* cells (control) | 8.5e-16 |
| unequal cells + `aw` + `subgroups()` + `adversarial` | 8.9e-16 |
| unequal cells + `aw` + cwm zero on every cell's **first** row | 8.7e-16 |
| unequal cells + `aw` + cwm zero on every cell's **last** row | 1.7e-15 |
| an `if()` that drops the first rows of some cells (so the pick row differs from the full-data first row, and out-of-sample rows exist) | 8.8e-16 |

The mechanism is why: the read is Σ_r w_r·d_g / Σ_r w_r with d_g cell-constant, so grouping gives
Σ_g W_g·d_g / Σ_g W_g with W_g = Σ_r w_r. Unequal counts and unequal weights are both absorbed
into the summed weight, which is exactly what `gcollapse (sum)` builds. The zero-weight first row
is harmless in the binary case for the same reason: the surviving row supplies only key values
(identical across the cell by construction) and the summed weights, never a value that varies
inside the cell. That last part is what changes under `contreat()`, where the row also supplies
the dose - hence proviso (ii).

**Is there a speed-up with no event-study horizons?** Yes - `esprelength()`/`espostlength()`
default to 0, so a vanilla `att contreat()` run produces only the main and pretreatment att and
ame rows, and the collapse still pays for itself (`ct/Q/qB_vanilla.do|log`, medians of three,
same admin fixtures as section 6; both admin fixtures collapse to 60 cells, so 4,258 rows -> 60
and 29,667 -> 60):

| fixture | request (no ES horizons) | collapsed | full | ratio |
|---|---|---|---|---|
| admin sub-sample (`e(N)` 4,162) | `att`, contreat | 0.295 s | 0.366 s | **1.24** |
| | `att`, binary | 0.139 | 0.152 | 1.09 |
| | `att att_it`, contreat | 0.363 | 0.449 | 1.24 |
| | all four estimands, contreat | 0.501 | 0.717 | 1.43 |
| admin full (`e(N)` 29,067) | `att`, contreat | 0.801 | 1.003 | **1.25** |
| | `att`, binary | 0.306 | 0.378 | 1.24 |
| | all four estimands, contreat | 1.033 | 2.087 | 2.02 |
| admin x5 (`e(N)` 145,335) | `att`, contreat | 3.369 | 4.490 | **1.33** |
| | `att`, binary | 1.106 | 1.349 | 1.22 |

So the cell-mean construction does not cancel the saving: the vanilla contreat gain is 1.24x at
4k rows, 1.25x at 29k and 1.33x at 145k, and it is *larger* than the binary vanilla gain at the
small size (1.09x) because `contreat()` adds the ame margins reads - more margins work to
accelerate against the same one-off collapse cost. The ratio scales with the number of margins
calls rather than with the row count: adding the other three estimands takes the 29k run from
1.25x to 2.02x, where adding a 2/2 event-study window takes it to 1.77x. The one regime where the
collapse is close to a wash is the smallest one with the fewest reads - binary, `att` only, 4k
rows, where the minimum-of-three ratio is 0.97. On these runs the collapsed and full arms are
bit-identical in `e(b)` (`mreldif` 0.000e+00 for the three vanilla `att` rows).

## 7. Loose ends to fix if this is built

1. **Stale refusal text** (A5). The `fastcollapse is incompatible with:` block still says
   "contreat() and interactivecontrols() interact the treatment dummy with a variable that
   varies inside a cell, so the effect is no longer cell-constant (measured: coefficients
   shifted, one sign flip)" - which the prototype contradicts - and has no bullet for the tokens
   that actually trigger it under contreat (`contreatpoly(>1)`, `lattice()`,
   `makecfxplots/cfxlattice()`, `contreat-in-controls`); the elasticity bullet does not say that
   dyex is allowed while eyex/eydx are refused. The help file's fastcollapse paragraph needs the
   same edit, plus the residue note for poisson+dyex if that is documented rather than gated.
2. **The prototype's own comments** need two corrections before they ship: the collapse-site
   comment names `untreatedZero` as one of the paths that put the ringer in the model, which it is
   not (section 1); and the worker comment's "measured: att-type reads moved by up to 0.045
   absolute" should carry the sharper V1b(f) measurement of section 2.2 (a pre-trend horizon moving
   p 0.167 -> 0.0022).
3. **Cosmetic** (A2, A5): the verbose skip line reads "does not apply withcontreatpoly(>1)"
   (no space: Stata trims the leading blank of `\`fcskip'`); the same missing space is in
   "incompatible with:". Pre-existing message shape, shared by every gate name.
4. **Guard and detector** (section 5) - **author-ruled 2026-09-02: not a build item.** The
   author's position is that `customweightmultiplier()` is caveat emptor, and every case in
   section 5.1 needs a signed cwm under `safety(off)`. The one-token guard change
   (`wsum > 0` -> `wsum != 0`) is still worth taking because it is measured as a no-op on 36 of
   37 configs and repairs the 37th at zero cost; the 15-line Σw = 0 detector
   (`ct/V4z/wooldid_ct_zerodetect.ado`, predicate `wsum == 0 & wsrc != 0 & !missing(wsrc)`) is
   built and priced and can sit on the shelf. Record instead a one-line doctrine note that
   signed cwm under `safety(off)` is outside the collapse's domain.
5. **Poisson + dyex**: gate or document (section 5). The measurement moved from "1e-6 on a
   variance" to "the fourth significant digit of a reported CI bound", which argues for the
   gate.
6. **The contreat subgroup section is not covered by the route-P degeneracy screen**
   (`wooldid_fcdeg` screens the binary subgroup reads). A4 found no branch flip in 44 configs
   and the only exact-zero contrasts (duplicate-cohort probes, the ppmlhdfe-fallback AME
   coincidence) differ between arms only by residue on cells whose se is itself residue. The
   one case where the arms genuinely part company is the `safety(off)` caution branch of the
   *binary* screen (`ctA4_sg3_cohcl_safoff_jt`, section 4) - which is the contract the author
   already ruled on for subgroups, not a contreat question. If the screen is extended to the
   contreat reads, the cheap hook is the same 1e-9 contrast-variance band the binary section
   uses, applied to the `pwcompare(effects)` reads at shipped 3322/3361.
7. **Never exercised under the prototype**: `cleanmatrices`, `oldsyntax` and `dots()` appear in no
   contreat config. `cleanmatrices` is the one worth a pin - the subgroups assessment hit an
   `r(3499)` do-file artefact around it - and the primary absorbing ppmlhdfe path (section 4) is
   the other measurement gap.
8. **The help file's residue paragraph** (`sthlp` 291-300) quantifies the collapse's residue as
   ~1e-16 relative on variances "regardless of the scale of the predictions" and ~1e-13 on
   near-zero coefficients, with no carve-out. The contreat measurements exceed both bounds (`e(V)`
   to 2.6e-11, a near-zero p-value to 4.2e-10, and poisson+dyex to 2.8e-4 if it is not gated), so
   that paragraph needs re-quantifying whether or not poisson+dyex is gated.
9. **Tests to pin**: a contreat F-series in `test/test_fastcollapse.do` (engaged PASS on
   dydx/dyex/poisson-dydx; NOTENGAGED on each gate token; the zero-weight-first-row config;
   the float-dose config; a jointtests+adversarial config that exercises the pooled reference
   cells), and the suite's `fcol_*` pins extended with a contreat member. The 328-config
   harness inventory (section 10) is the source for the pin configs.

## 8. Recommendation

Build it, as the prototype stands, with the section-7 items folded in. Three things carry that
recommendation. The arithmetic is exact by construction on the documented domain, and the reason
is written down in section 1 rather than asserted: every collapsible read is linear in the dose
within a cell, so a summed weight and a weighted mean reproduce it exactly. The probe found
nothing on that domain outside the section-5 disagreements - 305 dump-complete configs compared
object by object, 846,157 matrix cells, 8,881 scalars, no cell beyond last-place residue, with
poisson+dyex the one on-domain exception and a ruling of its own. And each proviso and each gate
is shown to be load-bearing by deleting it and measuring the damage, so the design has no
ornamental parts: drop the weighting and a subgroup pre-trend p crosses 0.05; drop `i` from the
key and a joint pre-trend p crosses 0.10; drop the recast and twelve joint p-values become
structural differences; drop the mean itself and the ATT is wrong in its third significant digit.
Section 6a adds the last practical question: the collapse pays for itself even on a vanilla run
with no event-study horizons (1.24x at 4k rows, 1.25x at 29k, 1.33x at 145k), so there is no size
or request shape at which the cell-mean construction eats the saving.

One question needed the author's ruling rather than more work; it has been answered, and one
remains:

1. ~~**`safety(off)` + signed `customweightmultiplier()`.**~~ **Ruled by the author 2026-09-02:**
   `customweightmultiplier()` is caveat emptor, so this is a documentation item, not a build one.
   Take the free `wsum != 0` guard anyway (a measured no-op on 36 of 37 configs that repairs the
   37th); leave the 15-line Σw = 0 detector built but unshipped.
2. **`poisson` + `contreatelasticitytype(dyex)`.** The residue here is not last-place: a
   reported CI bound moves 2.8e-4 relative, because margins differentiates
   `expression(predict(xb))` numerically. It is of the same order as the perturbation the input
   row order already produces on that read, so it is arguably the model's noise rather than the
   collapse's - but it is larger than anything else in this assessment. I recommend the gate (one
   token), so the help file's residue paragraph does not need an exception clause.

With the cwm question ruled and poisson+dyex gated, the build is the section-9 diff plus one gate
token, the guard token, and the section-7 documentation edits - and the help file's
existing promise ("results agree to machine precision but are not bit-identical", `sthlp` 291-300)
holds for contreat as written, with no exception clause.

One thing not to build: the extreme-cancellation regime of section 5.3. It is worth a line in
this repo's testing doctrine rather than in the ado - when a cell's doses cancel to four or more
digits, the two arms disagree because the *full-data* margins accumulation is unstable there, so
a future battery MISS in that regime is not evidence against the collapse and should be
adjudicated against the analytic value, not against `fullprecisionslow`.

## 9. The prototype diff

Against the route-P working tree of `wooldidinstall/wooldid.ado` (line numbers are the shipped
file's):

```
1835c1835,1847
<     wooldid_fcol `esamp' , key(`_i' `_t' `targetflag' `ifvar' `t' `ttre' `subgroups') wvars(`w_att' `w_att_it' `w_att_i' `w_att_itime')
---
>     * PROTOTYPE (contreat feasibility): with a continuous treatment the pooled reference-period
>     * cells (`_i'==0) enter the pre-trend joint tests under estimand weights that are not
>     * proportional across cohorts, so the user's cohort variable joins the key; the cell's
>     * weighted mean of the treatment variable is then written onto the surviving row
>     local fckeyx
>     local fccont
>     if "`contreat'" != "" {
>       local fckeyx `i'
>       * the ringer copy is averaged too: it enters as an explicit regressor under ccnoabsorb,
>       * untreatedZero and the ppmlhdfe fallback, where its level shifts xb round-off
>       local fccont contreat(`contreat' `ringer')
>     }
>     wooldid_fcol `esamp' , key(`_i' `_t' `targetflag' `ifvar' `t' `ttre' `subgroups' `fckeyx') wvars(`w_att' `w_att_it' `w_att_i' `w_att_itime') `fccont'
4694c4706,4726
<   if "`contreat'" != ""                     local fcbad `fcbad' contreat()
---
>   * PROTOTYPE (contreat feasibility): contreat() itself is allowed; only the reads that are
>   * nonlinear in the treatment variable, or that split cells on a variable outside the key, are refused
>   if "`contreat'" != "" {
>     if `contreatpoly' > 1                                                       local fcbad `fcbad' contreatpoly(>1)
>     if `"`lattice'"' != ""                                                      local fcbad `fcbad' lattice()
>     if "`cfxlattice'" != "" | "`makecfxplots'" != ""                            local fcbad `fcbad' makecfxplots/cfxlattice()
>     if "`contreatelasticitytype'" == "eyex" | "`contreatelasticitytype'" == "eydx" local fcbad `fcbad' contreatelasticitytype(`contreatelasticitytype')
>     * the treatment variable inside a user control term (controls(c.dose#c.x) and the like)
>     * makes the treatment derivative vary within a cell on a variable the collapse does not
>     * average, so any such term refuses the collapse
>     local fcctrl `"`controls' `coarsecohortcontrols' `fe' `ccc_absorb'"'
>     foreach tok in "#" "(" ")" "," {
>       local fcctrl = subinstr(`"`fcctrl'"', "`tok'", " ", .)
>     }
>     local fcctrlhit 0
>     foreach w of local fcctrl {
>       local w2 = regexr("`w'", "^[^.]*\.", "")
>       if "`w2'" == "`contreat'" local fcctrlhit 1
>     }
>     if `fcctrlhit' local fcbad `fcbad' contreat-in-controls
>   }
5133c5165
<   syntax varname(numeric) , KEY(varlist) [WVARS(string)]
---
>   syntax varname(numeric) , KEY(varlist) [WVARS(string) CONTreat(varlist numeric)]
5144a5177,5207
>   * PROTOTYPE (contreat feasibility): the surviving row carries the cell's estimand-weighted
>   * mean of the continuous treatment. Every estimand weight is proportional to w x cwm inside
>   * a cell, so the first existing weight variable serves as the analytic weight. The variable
>   * is recast to double first: a float store of the mean injects ~1e-8 into every read.
>   * The mean is formed as sum(w*x)/sum(w) and written onto EVERY row of the cell, not via
>   * gcollapse (mean) [aw=w] merge: that form writes MISSING onto zero-weight rows, and the
>   * surviving row is the cell's first row whether or not its weight is zero, so a
>   * customweightmultiplier() of exactly 0 on a cell's first row left the raw dose there
>   * (measured: att-type reads moved by up to 0.045 absolute). A cell whose weights are all
>   * zero keeps its raw values; it carries no weight in any read of either arm.
>   local cmv
>   if "`contreat'" != "" {
>     tempvar wsum
>     local wsrc
>     if "`wv'" != "" local wsrc `: word 1 of `wv''
>     foreach cv of local contreat {
>       local ctype : type `cv'
>       if "`ctype'" != "double" qui recast double `cv'
>       tempvar cm_`cv'
>       if "`wsrc'" != "" qui gen double `cm_`cv'' = `cv' * `wsrc'
>       else               qui gen double `cm_`cv'' = `cv'
>       local cmv `cmv' `cm_`cv''
>     }
>     if "`wsrc'" != "" qui gen double `wsum' = `wsrc'
>     else               qui gen double `wsum' = 1
>     qui gcollapse (sum) `cmv' `wsum', by(`key') merge replace
>     foreach cv of local contreat {
>       qui replace `cv' = `cm_`cv'' / `wsum' if `wsum' > 0 & `wsum' < .
>     }
>     qui drop `cmv' `wsum'
>   }
```

(Recommended edits per section 7: `wsum != 0 & wsum < .`; a `count if wsum == 0 & anynz`
refusal; `poisson` + `dyex` in the gate list.)

## 10. Artifacts

Session scratch directory (`/private/tmp/claude-501/.../scratchpad`, not in the repo):

- `ct/wooldid_ct.ado` - the prototype. Variants, each a surgical edit against it, every one
  diffed before use (five are single lines; `gateoff` removes three gate lines, `oldmean` restores
  the earlier mean block, `zerodetect` adds the detector): `ct/A1/wooldid_ct_signedguard.ado` (`wsum != 0`),
  `ct/V4z/wooldid_ct_zerodetect.ado` (signed guard + the Σw = 0 detector, with `variant.diff`),
  `ct/V1a/wooldid_v1a_{noi,unw,norecast,firstrow}.ado` (provisos iii, i, iv, v),
  `ct/V1b/wooldid_v1_{gateoff,oldmean}.ado` (gates off; the first prototype's mean form).
- Per-angle working directories with the config CSVs, findings notes and Stata logs:
  `ct/V1a/` (+ `FINDINGS_a.txt`), `ct/V1a_b/`, `ct/V1a_c/`, `ct/V1a_d/`, `ct/V1b_r/`,
  `ct/V4z/` (+ `cost.do`).
- `ct/A1..A7/*.csv` - the 328 round-2 configs; `adv2/<run>/RESULT.txt` with `run_C/`, `run_U/`
  `%21x` dumps and `wd_C/`, `wd_U/` logs for every run (`A1r_a..e`, `ctA2_a..c`, `ctA3_e..h2`,
  `ctA4_a..c`, `A5_e..h`, `ctA6_a2,b2,c2,c3`, `A7_r1..r6`, `ctD_*`, `ctV1a_*`, `ctV1b_*`,
  `ctV4_*`). The round-1 runs (`A1_a..d`, `ct_smoke*`, `A2_a..d`, `ctA3_a..d`, `A4_a..c2`,
  `A5_a..c`, `ctA6_a,b`, `A7_a`; 290 distinct configs, 276 with dumps) are in the same directory.
- `ct/A5/A5_explicit.log` (stand-alone gate requests and the refusal text), `A5_calib*.log`
  (poisson numeric-derivative calibration); `ct/A3/saveplots_cmp.log` (saved ES plot data
  comparison); `ct/D_*/` (diagnoser reproductions with print-hook variants).
- `ct/Q/qA_binary.do|log` (the binary unequal-cells / heterogeneous-weight probe of section 6a),
  `qB_vanilla.do|log` (the no-event-study timing), `qC_zeroweight.do|log` and `qD_nsample.do|log`
  (the zero-sampling-weight-versus-cwm question of section 2.2).
- `ct/V3/v3_timing.do|log`, `v3_scale.do|log` (timing; the trailing `r(4012)` in the scale log
  is an optional `profiler` block that this Stata build lacks - all SUMMARY lines precede it).
- `ct/r2_results.json` - the seven angle reports and two diagnoser reports as returned.
