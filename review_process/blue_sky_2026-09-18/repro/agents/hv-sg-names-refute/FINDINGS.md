# hv-sg-names-refute -- SKEPTIC findings

Pinned build asserted in every log (`which wooldid` ->
`.../scratchpad/bluesky/build/wooldid.ado`, `version 1.4 2026-08-19`).
All evidence is from MY OWN do-files on MY OWN generated data. I never ran the
original prober's `p1_names.do`.

Probe files (all self-contained, all with `.log` beside them):
- `r1_sgnames.do`    -- subgroup e(b)/e(V) names, help-name vs true-name addressing, V splice, adversarial scalars
- `r1b_consumers.do` -- coefplot counts, e(wdidsgresults_1) rows, bootstrap sg table keys
- `r2_esnames.do`    -- att_i / att_itime event-study names + the numeric consequence
- `r3_sharpen.do`    -- bit-identical shared ES block across att_it/att_i/att_itime
- `r4_ame_and_sg.do` -- contreat()/ame naming; subgroups + att_i; bootstrap column names
- `r5_prefix.do`     -- NEW: prefix-nested names + Stata `_b[]` abbreviation
- `r6_varabbrev.do`  -- NEW: `set varabbrev off` is the switch that blocks it

---

## VERDICT hv-sgname-01 -- CONFIRMED (severity medium, unchanged)

Reproduced exactly on my own 90-unit x 10-year panel (N=900, onsets
2004/2006/2008, 36 never-treated, sg in {1,2}).  e(b) has 48 columns, every
one FLAT (coleq = `_` for all 48 -- the one refutation angle that could have
bitten: there is no equation structure that would let [main_s1]att or
[main]att_s1 work).  The names are

    main_att main_att_it main_att_i main_att_itime
    pretreatment_att ... pretreatment_att_itime
    ES_att_pre2/pre1/post0/post1   ES_att_it_pre2/.../post1
    main_SG1_att main_SG1_att_it main_SG1_att_i main_SG1_att_itime
    pretreatment_SG1_att ... pretreatment_SG1_att_itime
    ES_att_SG1_pre2/.../post1      ES_att_it_SG1_pre2/.../post1
      (SG2 block identically)

`grep -c '_SG' wooldid.sthlp` == 0.  (Bare `SG` occurs 11 times, always as the
prose abbreviation for "subgroup".)

Addressing (r1 log 525-648): every help-derived name errors r(111) and every
true name works --

| reference | rc |
|---|---|
| `_b[main_att_s1]` | 111 |
| `_b[ES_att_post0_s1]` | 111 |
| `_b[ES_att_s1_post0]` | 111 |
| `test main_att_s1 = main_att_s2` | 111 |
| `lincom main_att_s1` | 111 |
| `esttab ., keep(*_s1*)` | 111 |
| `esttab ., drop(*_s1*)` | 111 |
| `_b[main_SG1_att]` | 0 |
| `_b[ES_att_SG1_post0]` | 0 |
| `test main_SG1_att = main_SG2_att` | 0 |
| `lincom main_SG1_att` | 0 |
| `esttab ., keep(*SG1*)` | 0 |

The silent path reproduces and is exactly one path, measured properly
(r1b log 210-267; coefplot's generate() variables are `cpAb`, not `cpA_b`,
which is why r1's first attempt returned missing):

| call | plotted coefficients |
|---|---|
| `coefplot` | 36 |
| `coefplot, drop(*_s1*)` | 36 -- SILENT NO-OP, rc 0, no message |
| `coefplot, drop(*SG1*)` | 24 |
| `coefplot, keep(*SG1*)` | 12 |
| `coefplot, keep(*_s1*)` | rc 0 but prints `(nothing to plot)` |

### Corrections / additions to the prober's account (none are refutations)
1. The `_SG` infix really is used consistently by the code, so their "the fix
   belongs in the help, not the code" holds against my independent check:
   e(bootdist) and e(V_boot) carry main_SG1_att, ES_att_SG1_post0, ...
   (r4 log 702-739), and the event-study covariance SPLICE lands --
   e(V)["ES_att_SG1_pre2","ES_att_SG1_post0"] equals
   e(V_es_att_s1)["pre2","post0"] to the bit on 10 of 16 cells and to 1 ulp on
   the other 6 (r1 log 803-818; e(V_es_att_s1) is itself asymmetric at the ulp
   level and ereturn post symmetrizes -- inside rtol 1e-10, not a finding).
2. Their count of THREE key TYPES stands, and the types reach more objects
   than they listed.  `_SG<v>_` infix on e(b)/e(V)/e(bootdist)/e(V_boot)
   (UNDOCUMENTED); `_s<v>` suffix on b_es_*/V_es_*/escomp_* (sthlp:1315) and
   also on e(adversarial_pval_att_s1) (sthlp:1398 -- documented, not a fourth
   regime); bare `_<v>` on wdidsgresults_* (sthlp:1318) and also on
   e(boot_wald_sgresults_1) / e(boot_pctl_sgresults_1).  Only the infix is
   undocumented.  Verified: e(b_es_att_s1) rc 0 / e(b_es_att_SG1) rc 111;
   e(wdidsgresults_1) rc 0 / _s1 and _SG1 both rc 111;
   e(boot_wald_sgresults_1) rc 0 / _s1 rc 111.
3. The failure is LOUD on every consumer the help names except coefplot's
   drop().  esttab, drop(*_s1*) errors 111 (a non-matching drop pattern is not
   a no-op in estout), and coefplot, keep() prints a message.  The one silent
   path is coefplot's own semantics for a non-matching drop() pattern, not
   something wooldid does.  That is why medium, not high, is the right
   severity -- but see NEW FINDING sk-prefix-01, the genuinely silent sibling.
4. Their "aggravating" note is accurate: the printed table shows
   "Results for Subgroup sg == 2:" with equation-style rows `main | att`,
   `ES_att | post0` -- no e(b) name is ever displayed, so the help is the
   user's only stated route to the name (r1 log 274-320).

---

## VERDICT hv-sgname-02 -- CONFIRMED (severity medium, unchanged), scope WIDER

Own panel (N=1920): 5 treated cohort-groups, 3 of 40 units with true effect 0.4
and 2 of 6 units with true effect 2.5, plus 60 never-treated units in 2 groups;
i = the cohort-group.

`wooldid y grp year ttre, att_i esprelength(2) espostlength(2) robust` posts
exactly 6 columns: main_att_i, pretreatment_att_i, ES_att_it_pre2,
ES_att_it_pre1, ES_att_it_post0, ES_att_it_post1.  att_itime alone gives the
identical set with main_att_itime.  _b[ES_att_i_post0] r(111),
_b[ES_att_itime_post0] r(111), _b[ES_att_it_post0] rc 0.  e(b_es_att_i) r(111);
e(b_es_att_it) rc 0.  The shared block is BIT-IDENTICAL across the three
estimands: ES_att_it_post0 = +1.39e70e8795764X+000 under att_it, under att_i and
under att_itime (r3 log 141-149), so the code is right and sthlp:1315 describes
it correctly.

The numeric consequence reproduces on my own design (r2 log 313-378):
`wooldid ..., att att_i` gives ES_att_post0 = 0.5787 (design truth,
cell-weighted, (120*0.4+12*2.5)/132 = 0.5909) and ES_att_it_post0 = 1.2262
(design truth, equal-weighted, (3*0.4+2*2.5)/5 = 1.24) -- a factor of 2.12.
A user who asked for att_i, found no ES_att_i_post0 and read the ES_att_* rows
as "the event study" would be 53% low.  As the prober said, that last inference
step is the user's; the posted numbers are right.

### Additions (the scope is wider than the claim states)
- The prober's "receives event-study rows labelled for an estimand they never
  asked for" is confirmed AND is visible in the PRINTED table, not only in
  e(b): on att_i alone the event-study block is headed ES_att_it while main and
  pretreatment correctly read att_i (r2 log 176-194).  A user who never touches
  e(b) still sees an estimand they did not request.
- The same false promise fires for the ame family, at a THIRD help location.
  See sk-amename-02.

---

## NEW FINDING sk-prefix-01 [medium / probable / realistic] -- SILENT WRONG
## NUMBER: estimand names are prefix-nested, so Stata's unique-prefix
## resolution silently redirects `_b[main_att]` to `main_att_i`

main_att < main_att_i < main_att_it < main_att_itime;
pretreatment_att < pretreatment_att_i < ...; main_ame < main_ame_i < ...;
main_SG1_att < main_SG1_att_i < ....  Stata resolves a coefficient reference by
unique prefix.  Demonstrated DIRECTLY on wooldid's e(b), where no dataset
variable shares any of these names (r6 log: `confirm variable main_att` and
`confirm variable main_att_i` both r(111)) and where `lincom main_att` echoes
the expanded `main_att_i` in its "( 1)" line.  So on a fit that posts exactly
one member of a family, a reference to a DIFFERENT member of that family
resolves silently to the one that is there.  (r5 leg 0 shows the same
resolution for ordinary regressors, but that leg is dataset-variable
abbreviation and is NOT the evidence for the wooldid case.)

Measured (r5_prefix.do, 192-unit design, log 177-322):

| fit | reference | value returned | rc |
|---|---|---|---|
| att alone | _b[main_att] | 0.47640617 (correct) | 0 |
| att_i alone | _b[main_att] | 1.21489080 (= main_att_i) | 0 |
| att_i alone | _b[main_att_i] | 1.21489080 | 0 |
| att_i alone | _b[pretreatment_att] | -0.00952806 (= pretreatment_att_i) | 0 |
| att + contreat | _b[main_ame] | +0.0772947 (correct) | 0 |
| att_i + contreat | _b[main_ame] | -0.02992339 (= main_ame_i) | 0 |
| att_i + subgroups() | _b[main_SG1_att] | 1.4232655 (= main_SG1_att_i) | 0 |

Ratio wrong/right on the att pair: 2.55.  On the ame pair the two numbers have
OPPOSITE SIGNS (+0.077 vs -0.030).

- `di _b[main_att]` / `scalar b = _b[main_att]` is COMPLETELY SILENT.
- `lincom main_att` and `test main_att = 0` return rc 0 but DO echo the expanded
  name in the "( 1)  main_att_i = 0" line (r5 log 210-231) -- a visible trace,
  easily missed and absent under `quietly`.
- `esttab ., keep(main_att)` errors 111 ("coefficient main_att not found").
- The misdirection fires ONLY when exactly one member of the family is posted:
  with att_i att_it together _b[main_att] is ambiguous and errors r(111)
  (r5 leg 3); with att att_i together the exact match wins and returns the
  correct 0.47640617 (r5 leg 4).

Why this is wooldid's problem and not merely Stata's: abbreviation is normally
harmless because a prefix names the SAME object (pri -> price).  wooldid is the
unusual case where the prefix of one estimand's name is another ESTIMAND -- a
different quantity with a different weighting.  sthlp:1305 presents main_att
first, as "the headline average treatment effect on the treated", and then says
the other estimands use "the same names with att replaced by ...", which invites
exactly the reference that misfires.  A user who switches a spec from att to
att_i and leaves `di _b[main_att]` in the post-estimation block gets a
2.5x-different number with no error, no caution and no trace.

Not in known_bydesign.md; its section 6 says explicitly "a missing hidden twin
is documented; a WRONG or COLLIDING name is not."

MECHANISM AND REMEDY (r6_varabbrev.do, log 114-161).  The redirection is
governed by Stata's `set varabbrev`, which applies to coefficient names as well
as dataset variables, and whose DEFAULT IS `on`:

    set varabbrev on  (default)  ->  _b[main_att] = 1.2422169  rc 0   (= main_att_i)
    set varabbrev off            ->  _b[main_att]  r(111);
                                     lincom main_att  r(111);
                                     _b[main_att_i] = 1.2422169  rc 0

No dataset variable is involved: `confirm variable main_att` and
`confirm variable main_att_i` both r(111) on the same data.  So the user-side
remedy exists (`set varabbrev off`) but is off the default path, and a wooldid
caution or a help note naming it would be cheap.  A code-side remedy would be
to stop making one estimand's name a strict prefix of another's (e.g.
`main_att_c` for the cell-weighted att), which is a number-moving rename and a
ruling for the author.

Repro: r5_prefix.do (the measurement), r6_varabbrev.do (the mechanism)

---

## NEW FINDING sk-amename-02 [low / certain] -- the att->ame naming
## promise at sthlp:239 fails for the ame event study exactly as sthlp:1305-1306
## does; this help location was not identified

sthlp:239: "Under contreat() every requested estimand is therefore reported
twice ... once as the effect of the binary onset (att, att_it, att_i,
att_itime), and once as the average marginal effect of the dose (ame, ame_it,
ame_i, ame_itime).  Both families are displayed, and both are stored under the
same names with att replaced by ame."

Measured (r4_ame_and_sg.do, log 203-243): `wooldid y grp year ttre, att_i
contreat(dose) contreatcontrols(basic) esprelength(2) espostlength(2) robust`
posts 12 columns:

    main_att_i  pretreatment_att_i  main_ame_i  pretreatment_ame_i
    ES_att_it_pre2 pre1 post0 post1   ES_ame_it_pre2 pre1 post0 post1

The LEVEL rows honour sthlp:239 (main_ame_i exists, rc 0) but the event-study
rows do not: _b[ES_ame_i_post0] r(111), _b[ES_ame_itime_post0] r(111),
_b[ES_ame_it_post0] rc 0; e(b_es_ame_i) r(111), e(b_es_ame_it) rc 0.  With att
requested instead, the ame ES block is ES_ame_* (r4 log 325-337), so the ame
family has exactly the same two-regime split the att family has.

`ame_i` is not a user-facing option ("option ame_i not allowed", r3 log 165-172)
-- the ame family comes along with the att family under contreat() -- so the
user's only route to the ame_i event-study name is the help, and the help gives
a name that does not exist.  sthlp:1315's correct account names only
att_it / ame_it, so sthlp:239 and sthlp:1306 are both wrong against it.

Repro: r4_ame_and_sg.do

---

## Negative results (things I attacked that are CLEAN)
- The _SG infix is internally consistent.  e(b), e(V), the spliced ES blocks in
  e(V), e(bootdist) and e(V_boot) all use it, and the splice is numerically
  right.  No mismatched or colliding _SG name found.
- No equation structure hides a `_s#` addressing route.  All 48 e(b) columns
  have coleq = `_`.
- The shared ES block is not a numeric bug.  ES_att_it_post0 is bit-identical
  under att_it, att_i and att_itime, and matches the equal-weighted design truth
  (1.2262 vs 1.24 with sd-0.20 noise on 5 groups), while ES_att_post0 matches
  the cell-weighted truth (0.5787 vs 0.5909).  The estimands genuinely coincide
  horizon by horizon, as ado:6828-6835 says.
- e(adversarial_pval_att_s1) matches sthlp:1398.  Not an undocumented regime.
- Negative and long subgroup values: not probed -- sthlp:858 refuses
  negative/non-integer sg values on treated rows and ado:12283 already cautions
  at >=8 digits (already ruled).
