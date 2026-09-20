# SKEPTIC REPORT — hv-lat-01..04 (maxvar preflight / lattice)

Workdir: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-lattice-maxvar-refute

## NEW FINDING (mine) — sk-maxvar-fp: the preflight is NOT a lower bound and DOES fire falsely

**Repro:** `s1_falsepos.do` (maxvar 2048 arm) and `s2_falsepos_default.do` (untouched
default maxvar 5000 arm), controls in `s3_controls.do`.

**Claim.** `si8demand = si8cells*(1+contreat) + si8ilev + si8tlev + 200` (ado:1371) is
documented as a LOWER bound on the regression's variable demand, and the documented
justification for the whole design is ZERO FALSE POSITIVES:

> ado:1342-1348: "a lower bound on the regression's variable demand, so an early refusal
> can only fire when the late r(907) was certain anyway."
> open_defects.md:74-79: "The design accepts false negatives ... in order to guarantee zero
> false positives, because a false positive would refuse a specification that would actually
> have estimated fine."
> known_bydesign.md section 5: same.

That guarantee is FALSE on the default estimation path. `ado:1476` is
`reghdfe y ... i._i#i._t#i.targetflag ... , absorb(i.`i' i.`t' `fe')`. The `i` and `t`
levels are ABSORBED — they are never design-matrix columns. So `si8ilev + si8tlev` is not
part of any bound; it is additive slack that grows one-for-one with the number of
NEVER-TREATED units.

**The invariant is asserted in FOUR project documents, and none of them records this
direction.** ado:1342-1348; `open_defects.md:74-79`; `known_bydesign.md` section 5;
`research_2026-09-16/guard_failure_modes.json:59` ("it cannot convert a would-be-successful run
into a refusal"). `grep -n "ilev|absorb|over-pred|false positive" guard_cost_2026-09-16.md`
returns only the UNDER-prediction direction (line 342) and the restatement of the invariant
(line 406). The over-prediction direction is nowhere in the record.

**Measured, at Stata's UNTOUCHED DEFAULT `c(maxvar)` = 5000** (s2_falsepos_default.log):
5200 units x 12 years, 20 treated (onsets 2012-2015), `wooldid y unit year ttre, att`
- DEFAULT safety: `Error: this specification requires at least 5542 variables (about 130
  treated cohort-by-period interaction terms) but maxvar is currently 5000, so estimation
  would fail with r(907).` -> rc 907, nothing estimated.
- IDENTICAL data, IDENTICAL maxvar, only `safety(off)` added: **rc 0 in 3.75 seconds**,
  att = 1.3745078 (truth 1.5), N = 62400.
- Same spec, DEFAULT safety, `set maxvar 8000`: rc 0, att = 1.3745078,
  hex +1.5fdfbe7391aecX+000 — identical to the safety(off) run, so `safety(off)` changed
  nothing about the estimate; the refusal was purely the variable count.
- The wooldid-shaped design, hand-built and fitted on the same data
  (`reghdfe y i._i#i._t#i.tf, absorb(i.unit i.year)`, s3_controls.log): **F(220, 56969)**,
  i.e. 220 estimated columns out of 240 (20 omitted for collinearity). Against a claimed
  "at least 5542 variables". Overshoot **>= 23x**, and the run needed under 5% of maxvar.
  (My hand-build is a generous UPPER proxy: it sets `_i` = unit on every row of a treated unit,
  whereas wooldid sets `_i` = 0 and `_t` = 0 outside post-treatment cells (ado:928-929, 1009), so
  wooldid's own design carries only the 130 nonzero cells plus the reference. The true overshoot
  is therefore larger than 23x, and I report the conservative figure.)

**The same at maxvar 2048** (s1_falsepos.log), which is Stata/BE's FIXED ceiling that a BE
user cannot raise at all: 2200 units, same 20 treated. Guard: "requires at least 2542
variables (about 130 ...)" -> rc 907. `safety(off)` at the same 2048: rc 0 in **1.91 s**,
att = 1.4700209.

**The controlled contrast.** The 2200-unit and 5200-unit panels have the SAME 20 treated
units, so `si8cells` = 130 on both and the design matrix is the same ~240 columns on both.
si8demand moves 2542 -> 5542, a change of exactly 3000 = the change in `si8ilev`. The term
that decides the refusal is the one that costs nothing.

**Two more realism anchors, both at the untouched default maxvar 5000.**
- `s8_zip_contreat.do` -- **continuous treatment**, 5,000 ZIP codes x 12 years, 60 exposed,
  `att contreat(exposure) contreatcontrols(basic)`, no lattice, no polynomial, no interactive
  controls. DEFAULT safety: *"requires at least 5992 variables (about 390 treated
  cohort-by-period interaction terms)"* -> **rc 907**. Identical command + `safety(off)`, same
  maxvar: **rc 0 in 68.18 s**, estimate .63452766 (truth .60), N = 60000. The design needs about
  390 x 2 + 1 = 781 columns; the guard claimed 5992. Overshoot ~7.7x. (Formula check:
  390*2 + 5000 + 12 + 200 = 5992, exact.)
- `s7_county_scale.do` -- 3,143 US counties x 18 years, 100 exposed, same options. DEFAULT
  safety: *"requires at least 5761 variables (about 1200 treated cohort-by-period interaction
  terms)"* -> **rc 907** (1200*2 + 3143 + 18 + 200 = 5761, exact). **Honest limit: the
  `safety(off)` arm did NOT finish inside my 500 s deadline and was killed**, so at county scale
  with a continuous dose I did not prove the run is feasible -- that design really is expensive
  (~2,401 columns on 56,574 rows). I report ARM 1 only and rest the finding on s1/s2/s3/s8.

**Realistic?** Yes. `i` = unit with a few thousand units and a handful of treated ones is an
ordinary administrative DiD: hospitals (~6k), school districts (~13k), nursing homes (~15k),
ZIP codes (~41k), establishments. At default maxvar the threshold is about 4,650 units; on
Stata/BE it is about 1,700. wooldid's own help (sthlp:1163) and its own caution steer the
user toward `fe()` + a coarser `i`, which is the right advice for RUNTIME and is presented
here as if it were required for the variable limit, which it is not.

**Path note (the precise statement).** On the `unconditionalse` / `customvce()` path the
command is `reg y i.`i' i.`t' `fe' ...` (ado:1490), where the `i` and `t` levels ARE real
columns. The formula is right for that path and wrong for the default reghdfe/ppmlhdfe path
— which is the path essentially every run takes.

**The block contradicts itself on the same data.** One line above the refusal, the same
`if "`safetyinput'" != "off"` block prints the granular-`i` caution: *"runtime grows faster than
linearly in that cell count ... so a granular unit can dominate runtime long before it approaches
Stata's variable limits."* It then refuses, on the variable limit, a run that takes 3.75 seconds.
Verbatim, in that order, in `s2_falsepos_default.log`.

**A harm beyond the message.** The refusal's named remedy is "swap `i' to a coarser level of
aggregation." Following it changes what is estimated: for the cell-weighted estimands
(`att_i`, `att_it`, `att_itime`) the cells themselves are redefined, so the estimand changes
outright; for the observation-weighted `att` the within-cell weighting changes. A spurious
refusal that steers the user to coarsen `i` pushes them off the quantity they asked for, for no
reason.

**Severity: medium.** It refuses a correct, fast, realistic specification, with a message
whose number is wrong by more than an order of magnitude and whose named remedy is
unnecessary. The escape hatches are `set maxvar` (impossible on BE) or `safety(off)`, which
disables ~60 other checks. Not critical: nothing silently wrong is computed, and the failure
is loud.

**Direction relative to the seeded findings:** hv-lat-01/02 say the bound is too SMALL. This
says the bound is ALSO too LARGE, on a different and more common axis, and that the property
the retraction relied on to close hv-lat-01's mechanism does not hold.


---

# VERDICTS ON THE SEEDED FINDINGS

## hv-lat-01 -> PARTIAL (severity revised medium -> low)

**Mechanism half: REDISCOVERY of a closed item.** "si8demand under-predicts, the run passes the
preflight silently and dies on Stata's own r(907)" is exactly:
- `known_bydesign.md` section 5 ("a deliberate LOWER bound ... Opened as a defect and retracted
  the same day, 2026-09-16");
- `review_process/open_defects.md:63-85`, the retraction itself;
- `review_process/guard_cost_2026-09-16.md:172-175` ("did not fire on a case that then failed
  with Stata's own ... cells=2004 ilev=200 tlev=20 demand=2424 against c(maxvar)=5000"), and
  line 342 ("the pre-emption does not reliably work ... under-predicts by more than 2x").
The prober's contribution on this half is a larger multiplier (4.6x on the lattice axis vs the
project's 2x on the plain axis). That is a sharper measurement of a ruled item, not a new defect.
The prober acknowledged the ruling and pivoted; the pivot is what survives.

**Doc half: SURVIVES, and is genuinely new.** The retraction is SILENT on the help text, and:
- `sthlp:1167`: "...if the count of variables the specification requires exceeds the current
  maxvar setting, wooldid halts with an informative error before estimation rather than failing
  mid-run with r(907)." No hedge.
- `sthlp:1207` lists it under "refuses, and safety(off) proceeds" alongside hard domain screens.
- `grep -c "lower bound" wooldid.sthlp` = **0**. The qualification exists only in an ado comment.
That is a documented-vs-actual divergence (brief category 5) and is the residue I confirm.

**My own reproduction, my own data** (`s4_lattice_own.do`, 40 states x 10 counties x 18 years,
i = state, continuous exposure, DEFAULT maxvar):
- si8demand PREDICTED 1030 (cells=386, ilev=40, tlev=18) -- preflight silent, correctly so here.
- realized nonempty design columns: no lattice ~775; `lattice(all 4)` ~3685; `lattice(all 5)` ~4191.
- under-prediction factor 3.6x and 4.1x on my shape (prober: 4.5x-5.2x on theirs). Confirmed.
- RUN 1 (no lattice) rc 0 in 13.15 s.
**End-to-end on my own data, my own shape** (`s6_lattice_big.do`, 50 states x 8 counties x 22
years = 8,800 rows, 46 treated, i = state, continuous exposure, UNTOUCHED default maxvar 5000):
- si8demand PREDICTED **1568** (cells=648, ilev=50, tlev=22) -> preflight SILENT.
- realized nonempty columns: no lattice ~1299; `lattice(all 5)` ~**6827**. Under-prediction 4.35x.
- RUN 1 (no lattice): rc 0 in **44.00 s**.
- RUN 2 (same + `lattice(all 5)`): **rc 907**, and the user sees only wooldid's fastcollapse
  caution followed by Stata's canned "maxvar too small ... it is currently 5000. Use set maxvar".
- 1568 < 5000 < 6827. hv-lat-01's arithmetic reproduces on a panel the prober never used.

**CORRECTION to hv-lat-01/hv-lat-03's "mid-run" framing.** RUN 2 took **0.194 seconds** end to
end (`timer list 2`). The r(907) is raised at factor-variable expansion, before any fitting, so
the cost of the missing preflight is a worse MESSAGE, not a wasted expensive estimation. The help
sentence "rather than failing mid-run with r(907)" (sthlp:1167) and the prober's "dies mid-run"
both overstate the practical harm -- no fitting is wasted, and what is lost is wooldid's
pre-estimation prep (validation, `_i`/`_t` construction, the three `gegen tag` passes, collapse
setup), 0.19 s on this 8,800-row panel and proportionally more on a large one. The dataset
survives intact (_N=8800, 11 vars),
`WOOLDID_LASTRC` = 907, and a bare `wooldid` replay is refused with
"Error: the last wooldid call failed with r(907); the stored results predate that failed call."
(r(198)) -- confirming hv-lat-03's side claims.

**Severity: LOW, not medium.** The mechanism is ruled and the ruling is explicit; what is left is
one unhedged help sentence plus a table row. The failure is loud, and Stata's own r(907) text
names `set maxvar`.

## hv-lat-02 -> PARTIAL (severity revised medium -> low)

True, and true by construction: `interinject` (ado:1171), `coarsecontrolinserter` (ado:1138) and
`ccc_absorbinserter` (ado:1145) are all appended to the same `reghdfe` call at ado:1476, and
si8demand reads none of them. But this is the SAME root as hv-lat-01 -- the bound models one term
family -- and the same ruling covers it. It is fix-scoping information, not a second defect.

**My own arms** (`s5_family.do`: 30 states x 8 counties x 14 years = 3,360 rows, 26 treated,
`set maxvar 2048`, NO lattice anywhere; si8demand = **642** on every arm, cells=199, ilev=30,
tlev=14, so the preflight is silent throughout):
- ARM 1 `att contreat(exposure) contreatcontrols(basic)`                       -> **rc 0**
- ARM 2 same + `interactivecontrols(x1 x2 x3 x4 x5 x6 x7 x8)`                   -> **rc 907**
- ARM 3 same + `coarsecohortcontrols(i.reg)` (12-level region)                  -> **rc 907**
- ARM 4 plain binary `att` + `interactivecontrols(x1..x8)` (si8demand 443)      -> rc 0
Eight demographic controls and a 12-level region factor are ordinary specifications, and neither
moves si8demand by one variable. Confirmed on my own data. Severity low, as a scoping note.

## hv-lat-03 -> CONFIRMED (severity low, as filed)

Independently confirmed as a by-product of my own runs: when the preflight DOES fire, Stata's
canned "maxvar too small ... You need to increase maxvar; it is currently N. Use set maxvar"
text prints immediately after wooldid's two lines (verbatim in `s1_falsepos.log` and
`s2_falsepos_default.log`). So the remedies are indeed SPLIT between Stata and wooldid, and under
`lattice()` only Stata's half survives. The prober's own downgrade of the seeded claim's
"uninformative" is correct. No correction needed.

## hv-lat-04 -> REFUTED (severity none)

The claim is that at the cell counts where the variable ceiling is first breached under
`lattice()`, runtime is still "tens to low hundreds of seconds", so the help's runtime caution
steers the user away from the real constraint. My own measurement points the other way:
- 40 states x 10 counties x 18 years, `lattice(all 4)`, ~3,685 realized columns -- BELOW the
  default ceiling -- was **still running after about 4 minutes** and I killed it myself to free a
  Stata slot (it was NOT deadline-killed, and I did not let it finish), against **13.15 s** for
  the identical spec without the lattice in the same session. So **>= 18x and climbing**, on a
  7,200-row panel, under 8-process contention.
Column ratio 3685/775 = 4.75x, and 4.75^1.9 ~ 19 -- my >=18x-and-unfinished ratio is consistent
with sthlp:1272's "lattice() is the steepest axis, at roughly the 1.9th power of K" rather than
merely not contradicting it. At column counts approaching the ceiling, runtime is emphatically
not unremarkable. The prober measured two points (i = 25 and i = 50) at which
their own report concedes the caution never printed at all, and concedes the sentence is an
asymptotic claim they did not test. A caution that does not print cannot mislead. REFUTED.

---

# WHY THE TWO HALVES ARE ONE PROBLEM (for whoever scopes the fix)

The 2026-09-16 retraction closed hv-lat-01's mechanism on this reasoning
(`open_defects.md:74-79`): *"The design accepts false negatives ... in order to guarantee zero
false positives, because a false positive would refuse a specification that would actually have
estimated fine ... trading recall for precision is the right way round."*

**There is no precision to trade for.** On the default estimation path the formula is wrong in
both directions at once:
- **too small** by every regressor family the bound does not model -- `lattice()` (ado:1087/1117),
  `contreatpoly()` (ado:1107-1116), `coarsecohortcontrols()` (ado:1138), `ccc_absorb()` (ado:1145),
  `interactivecontrols()` (ado:1171), plain `controls()`. Measured 3.6x-4.35x under lattice on my
  two panels; the project's own record has >2x on a plain panel.
- **too large** by `si8ilev + si8tlev + 200`, which `reghdfe ... absorb(i.`i' i.`t' `fe')` never
  turns into columns. Measured >=23x over on a 5,200-unit panel.
So the guard has neither recall nor precision, and the premise of the retraction does not hold.
The retraction should be reopened, and any fix should address both directions in one change.

Fix shapes, without prescribing one (author's call):
(a) split by path -- on reghdfe/ppmlhdfe drop `si8ilev + si8tlev`; keep the current formula only
    for the `reg`/`poisson` path at ado:1490, where those levels really are columns;
(b) count from the assembled regressor macros (`contreatinj`, `interinject`,
    `coarsecontrolinserter`, `contreatcontrolinserter`), as hv-lat-02 proposes;
(c) wrap the estimation command in `capture` and translate rc 907 into wooldid's own message
    naming the option that blew the budget (cheapest, and it is exact by construction);
(d) delete the arm -- `guard_cost_2026-09-16.md:339-341` already observes that `si8cells` and
    `si8tlev` exist *only* to form `si8demand`, so the >100-levels caution survives without it.
And separately, hedge `sthlp:1167` / `sthlp:1207`, which today promise a guarantee no version of
this formula can deliver.
