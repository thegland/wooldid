# FINDINGS — e() return contract (agent: ereturn-contract)

Method: contract testing against `recon/ereturn_map.md` plus automated assertions over a broad
configuration sweep (`ckc.doi` + `ckc2.doi`, ~25 invariants, ~70 configurations on 8 generated
panels). Every finding below has a self-contained repro do-file in this directory.

STATUS: COMPLETE.  ~70 wooldid configurations on 10 generated panels, 15 probe do-files and
logs in this directory.  Findings EC-1 … EC-9 below; self-contained repros are
`repro_EC1_bootrank.do`, `repro_EC2_level.do`, `repro_EC4_singular.do`, `repro_EC5_stale.do`,
`repro_EC7_sgnames.do`, `repro_EC8_restore.do` (each runs standalone: generates its own data,
`adopath ++` the pinned build, no dependency on anything else here).

SUMMARY
| id | title | severity | seeded? |
|---|---|---|---|
| EC-1 | bootstrap `e(V)` is rank-deficient at EVERY `reps()`, so `reps >= K` cannot fix it; at `reps(5)` `test` silently drops 3 of 7 constraints and returns p = 1e-16 | high | seeded (extends F7) |
| EC-2 | `level(99.9)` and `set level 99.9` make wooldid exit r(198) saying "at most two digits after the decimal point" | medium | no |
| EC-3 | a run that reports r(198) can leave a COMPLETE, CORRECT wooldid `e()`, and the replay refusal then says something false about it | medium | partly (F1 push-on) |
| EC-4 | `e(V_es_<fam>)` singular with a strictly positive diagonal; `e(es_ndegen)` is a diagonal-only scan, so the help's "its absence marks a clean covariance" is false | medium | no |
| EC-5 | the assigned failure-site enumeration: only 3 of 20 sites destroy the prior `e()`, one of them (`ccc_absorb`) is new, and the replay-staleness message never fires for exactly those 3 | medium | seeded (assigned) |
| EC-6 | `[aw=]`/`[pw=]` leave no trace in `e()` and `e(vce)` reads `ols` under pweights although the SEs ARE robust | medium | seeded (F6 push-on) |
| EC-7 | the help's `_s#` subgroup coefficient suffix does not exist in `e(b)` (`lincom main_att_s1` -> r(111)) | medium | seeded (F2) |
| EC-8 | after any failed wooldid, `estimates restore` gives a valid result set that wooldid refuses to replay | medium | no |
| EC-9 | `estimates save`/`use` is a lossless round trip, hidden entries included | CLEAN NEGATIVE | no |


## Corrections to the recon map (contradictions, per the coordinator note)
- **C1. recon F11's `e(roth_slope_ul_<type>)` "documented nowhere" is FALSE.** The help documents
  the pair at `sthlp:1381` as `{it:e(roth_slope_ll/ul_{it:type})}`; recon's
  `grep -c roth_slope_ul` misses the slash form. Both `roth_slope_ll_` and `roth_slope_ul_` are
  documented. Do not file.
- **C2. recon F10's "9991-9994 have no strlen guard" is FALSE.** `ado:9990` opens
  `if strlen("boot_pctl_sgresults_`sg'") <= 32 {` and that brace encloses the wald/pctl/bc/bca
  postings through 9995. They are guarded (conservatively, by the longest sibling's length), so
  the help's "skipped silently" promise at `sthlp:1293-1299` holds for all four.

## MEASURED RESULTS

### EC-1 (EXTENDS the seeded F7) — the bootstrap `e(V)` is rank-deficient at EVERY reps()
`p07_atomic.do`, `p07_atomic.log`, `p07b_rank.do`. Spec: 120 units x 12 periods, four
estimands + `esfixedbaseperiod esprelength(3) espostlength(3)`, `cluster(clu)`,
`bootstrap(standard)`, K = 20 coefficients, zero failed replications:

| reps | boot_completed | rank(e(V)) | K  | `test` of 7 coefficients | `test` of 3 ES coefficients |
|------|----------------|-----------|----|--------------------------|------------------------------|
| 5    | 5              | **4**     | 20 | p=1.155e-16 **on df=4** (3 constraints silently dropped) | p=6.624e-19, F=587.6 |
| 8    | 8              | **7**     | 20 | p=2.248e-37 on df=7      | p=3.437e-15, F=234.8 |
| 12   | 12             | **10**    | 20 | p=5.140e-16 on df=7      | p=1.571e-13, F=154.8 |
| 20   | 20             | **10**    | 20 | p=5.152e-12 on df=7      | p=8.531e-13, F=128.5 |
| 40   | 40             | **10**    | 20 | p=3.212e-12 on df=7      | p=1.768e-12, F=118.5 |

**The decisive control (`p07b_rank.log`):** the deficiency is caused by the ESTIMAND SET, not
by `reps`. On the same panel, `reps(60)`:
```
RANK|ANALYTIC_att_only |K=8 |rank(e(V))=8 |   RANK|BOOT60_att_only|K=8|rank(e(V))=8|completed=60
RANK|ANALYTIC_att_atti |K=16|rank(e(V))=16|
RANK|ANALYTIC_all4     |K=20|rank(e(V))=20|   (from p07: BOOT all4 K=20 -> rank 10 at reps 12/20/40)
```
`att` alone (K = 8) gives a FULL-rank bootstrap `e(V)` at 60 replications; adding `att_it`,
`att_i` and `att_itime` (K = 20) drops it to 10 and it never recovers.

Two things recon's F7 did not have:
1. **The rank SATURATES at 10 of 20 and never reaches K.** It is not an artefact of
   `reps < K`: at `reps(40)` — the number wooldid's own caution recommends — `e(V)` is still
   rank 10. The coefficient vector is structurally ~10-dimensional (main_att == main_att_it
   numerically, the pooled rows are weighted averages of the ES rows), so the replication
   covariance is CORRECTLY singular. Recon's proposed remedy "does the caution want
   `reps >= colsof(e(b))`" therefore would not fix it. Nothing marks the rank: `e(properties)`
   is `b V`, `e(boot_v_posted)` is 1, and the one caution printed is about the p-value GRID.
   (The analytic `e(V)` is rank 20 of 20 on the same spec, but that is just the ruled PS-62
   `diag(se^2)` construction, which cannot be singular — not news. The news is that
   `reps >= K` is not the boundary, so the remedy recon proposed cannot work.)
2. **At reps(5) Stata silently answered a different hypothesis.** `test` of seven named
   coefficients returned **df = 4**; the scripted caller reads `r(p) = 1.155e-16` and has no
   way to know three of its seven constraints were dropped. The F for the 3-coefficient test
   falls monotonically 587.6 -> 234.8 -> 154.8 -> 128.5 -> 118.5 as reps goes 5 -> 40, i.e. a
   five-replication run overstates the F by ~5x. A user is explicitly invited into this:
   `ado:6235` records "reps() is mandatory (author ruling: no floor - low-rep runs are
   legitimate)" and the caution text calls small reps "a legitimate way to time a specification".

### EC-2 — `level(99.9)` (and `set level 99.9`) make wooldid exit r(198) with a message that contradicts what the user typed
`p03_sweepA.log` A13 (discovery), **`p09_level.do` / `p09_level.log` (the grid)**.

Twelve levels, each put to `regress ..., level(L)`, to `set level L`, and to wooldid:
```
LEVELGRID|  90  |regress=0|setlevel=0|wooldid=0  |elevel=90    |lb_90ci
LEVELGRID|  95  |regress=0|setlevel=0|wooldid=0  |elevel=95    |lb_95ci
LEVELGRID|  99  |regress=0|setlevel=0|wooldid=0  |elevel=99    |lb_99ci
LEVELGRID| 97.5 |regress=0|setlevel=0|wooldid=0  |elevel=97.5  |lb_97p5ci
LEVELGRID| 99.5 |regress=0|setlevel=0|wooldid=0  |elevel=99.5  |lb_99p5ci
LEVELGRID| 68.5 |regress=0|setlevel=0|wooldid=0  |elevel=68.5  |lb_68p5ci
LEVELGRID| 99.9 |regress=0|setlevel=0|wooldid=198|             |         <-- r(198)
LEVELGRID| 90.1 |regress=0|setlevel=0|wooldid=198|             |         <-- r(198)
LEVELGRID| 99.99|regress=0|setlevel=0|wooldid=198|             |         <-- r(198)
LEVELGRID| 99.95|regress=0|setlevel=0|wooldid=0  |elevel=99.95 |lb_99p95ci
LEVELGRID| 68.3 |regress=0|setlevel=0|wooldid=0  |elevel=68.3  |lb_68p3ci
LEVELGRID| 92.5 |regress=0|setlevel=0|wooldid=0  |elevel=92.5  |lb_92p5ci
```
The message is `level() can have at most two digits after the decimal point`; 99.9 has ONE.
`regress y t, level(99.9)` prints a "[99.9% conf. interval]" header on the same data, and
`set level 99.9` is accepted by Stata itself. So this is wooldid-specific.

**The second route is worse: it fires with no `level()` option at all.**
```
. set level 99.9                                         // Stata accepts it
. wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3)
level() can have at most two digits after the decimal point
level() can have at most two digits after the decimal point
r(198)                                     // and e() still holds the PREVIOUS run: elevel=92.5
```
A user whose session sets `set level 99.9` cannot run wooldid at all, and is told about an
option they never typed.

**Mechanism, MEASURED (`p14_controls.log`).** `ado:551` parses `Level(real 0)`; `ado:574-576`
does `if `level' != 0 { global WOOLDID_LEVEL0 = c(level) ; set level `level' }`. Stata expands
a macro holding a real with `%18.0g`, and that is exactly the discriminator:
```
  level   macro expansion of `local x = <level>'      set level `x'
   90     90                                          rc 0
   97.5   97.5                                        rc 0
   99.5   99.5                                        rc 0
   68.3   68.3                                        rc 0
   99.95  99.95                                       rc 0
   92.5   92.5                                        rc 0
   99.9   99.90000000000001                           rc 198   <-- three decimals appear
   90.1   90.09999999999999                           rc 198
   99.99  99.98999999999999                           rc 198
```
The failing set is precisely the set whose `%18.0g` rendering needs more than two decimals.
The no-option route reaches the same line because `wooldid_validate` at `ado:5507` does
`if `level' == 0 local level = c(level)`, converting the session level into the same real.
Passing the user's literal text through, or `set level `=round(`level',0.01)'`, closes both.
Severity: medium — a crash with a misleading message on ordinary input, and it also blocks the
`lb_99p9ci`/`ub_99p9ci` interval-column naming rule the help documents at `sthlp:1310`.

### EC-6 (NEW; answers recon F6's open push-on) — `[aw=]` and `[pw=]` leave NO trace in `e()`, and `e(vce)` says `ols` under pweights
`p09_level.do` / `p09_level.log`:
```
. wooldid y id t ttre [aw=awt], att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3)
  AW_RC=0  wtype=  wexp=  vce=ols
. wooldid y id t ttre [pw=pwt], att att_i att_it att_itime ...
  PW_RC=0  wtype=  wexp=  vce=ols  vcetype=
. [fw=...] -> r(198) (documented refusal)   [iw=...] -> r(101)
```
`e(wtype)` and `e(wexp)` are the two entries every Stata estimation command posts to record its
weights, and `esttab`/`estimates table`/any replication script read them; wooldid posts
neither, so a stored `.ster` cannot tell a weighted run from an unweighted one. Worse, under
`[pw=]` — which implies a robust VCE in every Stata command — `e(vce)` reads `ols` and
`e(vcetype)` is empty, so nothing in `e()` says the standard errors are (or are not) robust.
This is the exact open push-on recon F6 listed ("pweights (imply robust - does `e(vce)` say
so?)"): the answer is **no**. Compare `cluster()`, which does post `e(vce)=cluster`,
`e(vcetype)=Robust`, `e(clustvar)` and `e(N_clust)`.

**The NUMBERS are right; only the labels are wrong** (`p14_controls.log`, same data):
```
SE|pweight       |att se = 0.210058487 | vce=ols    | vcetype=       | wtype=
SE|aweight+robust|att se = 0.210058487 | vce=robust | vcetype=Robust | wtype=
SE|aweight       |att se = 0.176744105 | vce=ols    | vcetype=       | wtype=
SE|noweight      |att se = 0.167907539 | vce=ols    | vcetype=       | wtype=
```
`[pw=]` reproduces `[aw=] robust` to nine decimals, so the estimator is doing the right thing —
`e(vce)`/`e(vcetype)` simply fail to say so, and `e(wtype)`/`e(wexp)` are missing on every
weighted run. Cause: the single `ereturn post` at `ado:4051` is
`ereturn post `eb' `ev', esample(`esamp') depname(`y') obs(`N') dof(`df_r')` — it has no
weight clause, and nothing later posts `e(wtype)`/`e(wexp)`.

### EC-7 (CONFIRMS the seeded recon F2, with a one-line repro) — the help names `e(b)` subgroup coefficients with a suffix the code never produces
`p09_level.do` / `p09_level.log`. `sthlp:1301-1306` is the paragraph that tells users what
`esttab`, `coefplot`, `lincom` and `test` address a coefficient by, and it ends:
*"Under `subgroups()` each subgroup's rows take the suffix `_s#`."*  `_SG` appears ZERO times
anywhere in `wooldid.sthlp`.
```
. wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) ///
        espostlength(3) subgroups(sg) cluster(clu)
. lincom main_att_s1        -> r(111)      <- what the help tells you to type
. lincom main_SG1_att       -> rc 0, est = 3.0590    <- what the code produces (ado:4015)
. test main_SG1_att = main_SG0_att -> rc 0
```
Every OTHER subgroup-keyed name does use the documented suffix, which is what makes the one
exception so easy to trip over: `e(b_es_att_s1)` (matrix), `e(numpre_att_s1)` (scalar),
`e(wdidsgresults_1)`, `e(escomp_s1)` all exist under `_s#`.

### EC-3 (NEW) — a run that reports r(198) can leave a COMPLETE and CORRECT wooldid `e()`, and the replay refusal then tells the user something false about it
`p10_partial.do`, `p10_partial.log`, `partial_names.txt`.

recon F1 enumerated two post-failure `e()` states (the prior e() untouched; reghdfe's internal
saturated model). There is a third, and I pinned its exact site:

```
. wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) ///
        subgroups(sg8) jointtests adversarial cluster(clu)          // sg8 = 10000000/10000001
  ... the full main table, the ES table and the joint tests all PRINT ...
  adversarial_pval_att_it_s10000000 invalid name
  r(198)
. di "`e(cmd)'"          ->  wooldid
. di colsof(e(b))        ->  24
  e() name set is BYTE-FOR-BYTE the same 77 entries as the identical run with 7-digit
  subgroup values, which exits rc 0.        (GOOD_NSCALARS=34 == FAIL_NSCALARS=34)
. lincom main_att        ->  rc 0, est = 3.0946, se = .30871     (correct numbers)
. estimates store X      ->  rc 0
. wooldid                ->  "Error: the last wooldid call failed with r(198); the stored
                              results predate that failed call. Refit wooldid before
                              replaying."   r(198)
. macro drop WOOLDID_LASTRC
. wooldid                ->  r(198) again
```
**The crash site is a READ, not a post.** `ado:11577` in `wooldid_display` is
`if e(adversarial_pval_`typetext'_s`sgcase') < .`, and `adversarial_pval_att_it_s10000000`
is 33 characters, so the *reference* `e(<name>)` is a hard r(198). The matching posting site
(`ado:4349`) never ran for `att_it` — it is gated on `cap confirm matrix` — so nothing is
missing from `e()`; the display loop simply builds a name it can never read. This is why the
`sthlp:1293-1299` "skipped silently when the name would run past 32 characters" contract does
not save this case: there is no twin to skip.

Two consequences, of which the second is the general one:
1. The estimation is COMPLETE and every stored result is correct, yet the call reports an
   error and (because `WOOLDID_LASTRC` is set) the results become unreachable through
   wooldid's own replay. `macro drop WOOLDID_LASTRC` does not rescue them either — replay
   re-enters the same display and dies again. (recon F4's `macro drop` defeat IS real for
   pre-post failures — measured in `p13_replay.log`: `R4_REPLAY_AFTER_MACRODROP_RC=0`, the
   stale results redisplay — it just does not apply when the failure was in the display.)
2. **The replay refusal's message is false for every failure that happens after
   `ereturn post` (`ado:4051`).** It says "the stored results predate that failed call"; here
   they ARE that call's results, complete and correct. Any failure in the display or plot
   stage produces the same false statement.

Scope note, stated honestly: the r(198) ITSELF is inside the ruled 8-digit-subgroup caution
family (`wooldid_solved_issues.md:180`; `sthlp` names `adversarial` as one of the exposed
options), so do not re-file that. What is new is (a) the site is a READ in `wooldid_display`,
which a `strlen` guard there would turn into one skipped display line instead of a crash, and
(b) the false replay message, which is not specific to long names.

### EC-4 (NEW) — `e(V_es_<fam>)` is SINGULAR with a strictly positive diagonal on an ordinary clustered design, and `e(es_ndegen_<fam>)` cannot see it, contradicting the help's "its absence marks a clean covariance"
Discovery `p06_ragged.do`; isolation `p11_singular.do`; **self-contained repro
`p12_ndegen.do` / `p12_ndegen.log`**.

The help, in the `e(b_es_*)/e(V_es_*)` entry, promises:
> "...horizons whose standard errors are not identified then carry zero (structural-zero) or
> missing diagonal rows, flagged by a Caution and by `e(es_ndegen_<type>)` (the count of
> degenerate rows, posted only when such rows exist — **its absence marks a clean covariance**)."

`ado:4130-4136` computes that count by scanning **only the diagonal**:
```stata
forval f2d = 1/`f2k' {
  if `f2vm'[`f2d',`f2d'] > 0 & `f2vm'[`f2d',`f2d'] < . local f2ok 1
  else local ++f2nbad
}
```
A block that is rank-deficient *off* the diagonal is therefore never counted, never cautioned,
and `e(es_ndegen_*)` is never posted — which the help then reads as proof the matrix is clean.

**MEASURED (`p12_ndegen.log`).** One panel: 120 units x 12 periods, three treated cohorts at
t = 4/7/10 plus never-treated, `clu = ceil(id/6)` (20 clusters). One window:
`esfixedbaseperiod esprelength(3) espostlength(5)` -> H = 8. Only the VCE changes:

| run | `e(V_es_att)` | rank | eigenvalues (sorted) | min diagonal | `e(es_ndegen_att)` | `cholesky()` |
|---|---|---|---|---|---|---|
| `cluster(clu)` | 8x8 | **6** | -8.35e-17, -8.37e-18, 4.09e-04, 2.60e-03, 8.44e-03, 3.09e-02, 7.89e-02, 1.77e-01 | 1.57e-02 | **absent** | **36 missing cells**, rc 0 |
| (default, iid) | 8x8 | 8 | 1.43e-02 … 1.85e-01 | 3.06e-02 | absent | 0 missing |
| `robust` | 8x8 | 8 | 3.58e-03 … | 1.93e-02 | absent | 0 missing |
| `cluster(clu)`, `att` alone | 8x8 | **6** | identical to row 1 | 1.57e-02 | **absent** | 36 missing |
| `cluster(clu)`, `esprelength(2) espostlength(2)` (H=4) | 4x4 | 4 | 3.60e-04 … | 2.09e-02 | absent | 0 missing |

So it is the **cluster-robust VCE plus a window wider than the design independently
identifies** that makes the block singular — not the estimand set (row 4: `att` alone gives the
identical rank-6 matrix) and not a degenerate diagonal (every diagonal entry is >= 0.0157) and
not asymmetry (`VASYM = 0` on every run in the whole sweep). Two eigenvalues sit at ~1e-17.
The same geometry at subgroup level is in `p11_singular.log`
(`BLK|A_shortpre6|att|H=12|rank=9`, `BLK|B_sgonecohort|att|H=8|rank=6`), and it propagates to
`e(V)` itself: `VTOT|C_nosg|K=24|rank(e(V))=20`, `VTOT|A_shortpre6|K=76|rank(e(V))=68`.

**The control the finding needs, MEASURED (`p14_controls.log`).** The help says the design that
SATISFIES its conditions is "a coarser `i` (cohorts, or treatment-onset groups) with unit fixed
effects supplied through `fe()`". Same data, same window:
```
SH|i_is_unit_cluster      |H=8|rank=6|mindiag=1.57e-02|cholmiss=36|es_ndegen=.|rank(e(V))=20/24
SH|i_is_cohort_fe_cluster |H=8|rank=8|mindiag=2.35e-02|cholmiss= 0|es_ndegen=.|full rank
SH|i_is_cohort_fe_robust  |H=8|rank=8|mindiag=3.19e-02|cholmiss= 0|es_ndegen=.|full rank
SH|i_is_cohort_fe_iid     |H=8|rank=8|mindiag=3.16e-02|cholmiss= 0|es_ndegen=.|full rank
```
So the singularity is confined to **`cluster()` + a unit-level `i`** — which is exactly the
geometry the help's "Design conditions for every consumer of `e(V_es_*)`" block already flags
in prose ("a treated (i, relative-period) cell lies wholly inside one cluster"; "the covariance
carries control-cluster noise only"; "neither condition below is visible in the matrix
itself"). **That downgrades this from a statistical surprise to a marker defect**, and I grade
it MEDIUM accordingly.

**Honest scope.** The *statistical* hazard is documented. What is NOT is (a) the
machine-readable claim that
the absence of `e(es_ndegen_<type>)` marks a clean covariance — here it marks a rank-6-of-8
matrix — and (b) that the degeneracy screen is a diagonal-only scan, so no rank or eigenvalue
check stands between a singular block and `pretrends`, `pretrendpower()`
(`ado:10444-10445` feeds it `e(b_es_*)`/`e(V_es_*)` directly), `lincom`, `test`, or any user
`cholesky`/`invsym`. `pretrends power` ran rc 0 on the singular block (`p12_ndegen.log`);
**`honestdid` was NOT successfully tested** — both of my attempts were my own syntax error
(`vcv()` resolved to `e(V)` and errored on conformability), so I make no claim about it.
Adding `rank()` or a smallest-eigenvalue test beside the existing diagonal loop at
`ado:4130-4136` would close the marker gap. A unit-level `i` with `cluster()` is, on the
project's own evidence, the most common way this command is run.


### EC-8 (NEW) — after any failed wooldid, `estimates restore` gives you back a valid result set that wooldid then refuses to replay, with a false explanation
`p13_replay.do` / `p13_replay.log`. `sthlp:1299` states "The commands `est restore` and
`est save` are compatible with `wooldid` for saving its output".
```
. wooldid ... (good)          -> rc 0
. estimates store GOODM
. wooldid y id t ttre, att subgroups(nosuchvar)    -> r(111);  $WOOLDID_LASTRC = 111
. estimates restore GOODM     -> rc 0,  e(cmd)=wooldid,  K=16      (a VALID result set)
. di "$WOOLDID_LASTRC"        -> 111                                (still set)
. wooldid                     -> r(198) "Error: the last wooldid call failed with r(111);
                                 the stored results predate that failed call."
```
The same after a post-regression failure (`P1_FAIL_RC=322`, `e(cmd)=reghdfe`; restore brings
back K=16; `P3_REPLAY_AFTER_RESTORE_RC=198`). The refusal is keyed on a session global
(`ado:40-46`) that no `estimates restore` can clear, so the guard cannot tell "stale" from
"deliberately restored", and its message asserts something untrue about what is in `e()`.
The natural fix is to clear `WOOLDID_LASTRC` when `e(cmd)=="wooldid"` results arrive from
somewhere other than the failed call — or to key the guard on something `estimates restore`
carries, e.g. comparing `e(cmdline)` with the failed call's.

Related measurements in the same log, for the record:
* `$WOOLDID_LASTRC` survives `clear` (M2) and `discard` (M3), and `macro drop _all` clears it
  (M4) — so recon F4's `macro drop` route is confirmed: after a parse-stage failure,
  `macro drop WOOLDID_LASTRC` then bare `wooldid` **replays the stale results, rc 0**
  (`R4_REPLAY_AFTER_MACRODROP_RC=0`, and `e(cmdline)` still shows the OLD command line).

### EC-9 (CLEAN NEGATIVE, recon §8 listed it as untested) — `estimates save` / `estimates use` is a lossless round trip for wooldid, hidden entries included
`p13_replay.log`, block 5. Saving a 4-estimand + fixed-base event-study run to a `.ster` and
reading it back into an empty session loses **nothing**:
```
T4_LOST_SCALARS:      (empty)
T5_LOST_MACROS:       (empty)
T6_LOST_MATRICES:     (empty)
T7_RESTORED_CMD=wooldid  K=16
T9_REPLAY_AFTER_ESTUSE_RC=0     <- bare `wooldid' replays the restored run
```
Every hidden entry (`wdid_vcvstate`, `wdid_vblocks`, `V_analytic`, `df_r_analytic`, …) survives,
and the replay path works off the restored e(). The only thing that does not survive is
`e(sample)` (`T8_ESAMPLE_AFTER_USE=0`), which is Stata's own behaviour when the data is gone.
