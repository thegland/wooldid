# Mini Wooldid Assessment — 2026-08-23

Static-analysis pre-release sweep ahead of the full testing-and-assessment run, per the
work file's "Mini Wooldid Assessment" chunk. Method: a 10-lens static-analysis fan-out
(five close-read sections over shell/worker/validate, plus subprogram, permute+mata,
display/message, docs-drift, and cross-cutting-hygiene lenses; ~48 raw findings, deduped
and centrally adjudicated against the campaign ledger). No testing-suite runs were used
for the assessment itself; fixes were gated by a parent-vs-edited two-build smoke
harness (ad hoc, scratchpad-only) plus targeted static verification of every claim.
Out of scope by direction: technical-variable overwrite issues, deep statistical
correctness (adversarial calibration, RI validity — next phase), performance.

## 1. Implemented (commits, this branch)

- **`a5abbed` gitignore coverage** — the nine tranche-16 test suites write log-scrape
  artifacts under ten new name families (`advstat_*`, `cfxlat_*`, `cfxtok_*`,
  `esfsup_*`, `perm2_*`, `ptpd_*`, `ptpr_*`, `sgdom_*`, `vald_*`, `w3f_*.txt`) that
  `test/.gitignore` predated; all 26 previously-untracked `.txt` files verified
  regenerable (every reference is a `log using ... replace` or a same-test `fileread`
  scrape). Root `.gitignore` added carrying only `.DS_Store`.
- **Comment-accuracy pass** — zero code change: the shell replay comment now records
  the deliberate NFRPL-01 exception; the sample-membership markout note contradicted
  the loop below it since the original upload; the F6 gate comment still described the
  commit-A transitional state ("engine lands in a follow-up commit"); the F7 subprogram
  banner was stranded above `wooldid_permctl` when F6 inserted the permute programs
  (relocated to `wooldid_ptpcheck`); `wooldid_permx`'s header omitted `wdp_caut0`;
  typos (usge ×2, misspecifeid, doens't, udpated, eachother ×4, "initialize some
  times").
- **Dead-code removals** — eleven provably-unread vestiges, each grep-verified:
  worker `relyrlst` (typo'd sibling of `relyrlist`), `converged` (captured, never
  read or posted; the ship-e(converged) question stays with the author), the vestigial
  `sgnum` init/increment pair in the e(b) stacking block; permrun's pre-rep-0
  comparator-seeding loop (dead store since NFPERM-10 made rep-0 mandatory; its
  ULP/SCALARS rationale moved into the rep-0 comment the shell already
  cross-references), the reader-less cluster count (`wdpnclus`), unused tempname
  `distk`, and the empty placeholder `forval` loop; esplots' two dead
  `initialflagt` stores; mata unused declarations (permx: `D`, `blockperms`,
  `counts`, `S`, `j`; permx_enum: `sub`, `ustr`).
- **Console-message and format fixes** — user-visible text only, no rc or numeric
  changes: multicolinearity→multicollinearity (3 messages + 2 comments) and the
  trailing "warning"→"caution" on its count line; "Error:" prefix added to the
  non-integer-subgroup r(411) message; subject–verb repair of the empty-sglist r(459)
  message; "pre-treatment units"→"periods" in the esf insufficient-sample note
  (timetrends consumes periods, not units); "counter intuitive"→"counterintuitive"
  ×2; "safety = off"→"safety(off)" (sole non-canonical spelling); "suppress this
  warning message"→"suppress this caution"; verbose-stream fixes ("for the full
  sample", "in preparation for", the contreat completion line now gated on contreat
  so binary-only runs no longer announce continuous-treatment work, the per-estimand
  AME progress line now names its estimand instead of printing identically 4×); cfx
  plot-title capitalization aligned to its six siblings ("Post-treatment"/
  "Pre-treatment"); `%-8.4f`→`%-6.4f` on the lone divergent SG joint-test p line; R2
  header separator spacing normalized. Verified: none of these strings is a
  fileread anchor in any `test/*.do`.
- **esplots subgroup/timetrends sibling line** — the subgroup fixed-base branch
  omitted the timetrends min-adjustment its main-block sibling has
  (`if "timetrends" ... local min = esrelativeto - 1`), so under esfixedbaseperiod +
  timetrends + subgroups + makeplots with esrelativeto below the reported window, the
  second planted reference point (at esrelativeto−1) fell outside the subgroup plots'
  xlabel range and lost its tick/asterisk. Inserted the sibling line verbatim; the
  planted data row and asterisk logic already existed in the branch. Plot-label-only
  change; probed on the firing combo and non-firing siblings.
- **Docs** — sthlp reserved-names sentence singularized: `__marginsmodel` machinery
  was deleted by 4fcf127 (NFADV), so the claim that wooldid "overwrites and removes"
  that name on every run was stale; `___wooldidfullmodel` half kept verbatim.
  Citation-year alignment completed in the direction the author set in acd9f7e
  (pkg/README → Wooldridge 2025, Empirical Economics): the sthlp's 7 in-text
  "Wooldridge (2021)" mentions, the sthlp Citations entry (was: Working Paper 2021 +
  SSRN link), and the CITATION.cff abstract now match. CITATION.cff's
  `date-released`/`preferred-citation` (2023 SSC release) deliberately untouched —
  refresh-at-release item.

## 2. Tickets for the author (found, verified, NOT implemented — behavior-changing)

- **T-M1: `evertreated` coding inverted between siblings** (worker L992 vs L850, both
  since the original upload). The contreatcontrols-block fallback defines
  1 = never-treated while the lattice-block sibling defines 1 = ever-treated; the same
  `contreatcontrols(byTreated)` request gets opposite codings depending on whether
  latticetype also contains byTreated (the L990 guard decides which definition the
  consumers at L1023-24 receive). Impact is confined: inside absorb() the variable is
  grouping-only (bit-identical either way); on native-factor paths the pooled
  `i.null##ringerpoly` term from L1014 equalizes spans; the movable corner is
  ccnoabsorb-without-pooled-term, where base-level/collinearity-drop choices — and
  hence numbers — can shift. Related: **L1014's pooled-term condition omits
  ccnoabsorb** even though its comment says it should cover any FE object kicked
  outside an HDFE absorb (L1049-52 does exactly that) — the two should be fixed
  together, with pins. Numeric change ⇒ graduation + suite pinning required.
- **T-M2: continuous-SG joint pre-trend p stored as missing** (L2685/L2696 read
  `r(p)[1,2]` where every continuous sibling reads `[1,1]`; binary sites correctly
  read `[1,2]` for the two-column dydx layout). Empirical evidence already in the
  blessed baseline: `cont_sg.csv` has `joint_Choriz_s*` = "." while the same run's
  overall continuous `joint_Choriz` and binary-SG `joint_horiz_s*` carry real p's;
  the conditional-ereturn gate proves the assignment executed and stored literal
  missing. Old wooldid shares the bug, so old-vs-new comparison can never flag it. A
  prior record of this question ("logged below-the-line" during B6 verification)
  dangles into a since-rewritten review file. Fix is likely `[1,2]`→`[1,1]` at both
  sites but must be probed for the actual r(p) layout, and `cont_sg`/`cont_sg_gap`
  baselines pin the defective "." ⇒ graduation + re-bless required.
- **T-M3: adversarial horizon-gate disjunct is provably redundant** (L682):
  `(espre + espl + 1 < 3 & timetrends != "")` is implied by the standalone
  `espl == 0 | espre == 0` clauses for the non-negative integers that reach it.
  Either delete it (behavior-neutral) or — if the original intent was a *stricter*
  timetrends threshold (cf. the L633 message that timetrends raises required
  pre-periods) — strengthen it, which would be an rc-change. Author's call; code left
  untouched.
- **T-M4: latent r(b)-staleness window in the SG-binary adversarial block**
  (L2436-38 reads `r(b)` after the conditional r()-clobbering `glevelsof` fallback at
  L2399). Currently provably inert (the fallback fires only when colsof ≤ 2, and both
  clobbered and unclobbered paths fail the `> 2` gate identically, printing the same
  skip Caution); the sg-continuous and both main siblings have no clobber in their
  windows. Hazard-hygiene candidate for the deep phase: capture colsof while r(b) is
  authoritative (the f2bmat/f2vmat pattern three lines up).

## 2b. Author rulings on the tickets (2026-08-23) and their outcomes

- **T-M1: ignore** (author ruling; item stands recorded above for any future pass).
- **T-M2: author challenged the analysis ("checking the actual matrices would show
  otherwise"); the matrix check CONFIRMS the defect.** Through an instrumented
  scratchpad build on the exact battery `cont_sg` configuration (frozen
  `dta_cont.dta` + identical option string), the r(p) returned by the sg-continuous
  overjoint contrast is **1×1** with the real joint pre-trend p at [1,1] — sg1:
  6.209e-06, sg2: 4.109e-07 — while the code reads [1,2]. Stata's expression-context
  matrix subscripting silently returns missing for out-of-range cells (probed:
  `A[1,2]` on a 1×1 matrix gives "." with rc 0), so `e(joint_Choriz_s*)` stores "."
  and the console prints "." where those p-values exist. In the same run the main
  continuous site (structurally identical margins, reads [1,1]) returned .00308048
  into e(joint_Choriz), and the main binary site's r(p) is 1×2 with its p at [1,2]
  (base column missing) — confirming the per-shape index convention the siblings
  follow. The fix is [1,2]→[1,1] at the two sg-continuous sites; it changes e()
  values ("." → real p) and `test/baseline/cont_sg.csv` pins the defective "."
  (sole pinning config). Note the full old-vs-new assessment CANNOT catch this
  class on its own: the old build shares the bug, so both sides show matching "."
  cells.

  **AUTHOR RULING (2026-08-23): DEFER — address during the full assessment run,
  under a mandated evidence protocol.** Second-round verification (author requested
  full-run proof, skeptical of the toy-vs-full distinction): two cases × two builds
  — the battery `cont_sg` configuration and a purpose-built well-structured case
  (8 coarse cohorts × 20 units × 12 periods, unit-level dose varying within
  cohort-period cells, clean parallel pre-trends by construction, subgroup-specific
  post effects), each run through the stock build (untouched standard console
  output) and an instrumented build (r(p)/r(b) dumped mid-run immediately before
  the exact read lines). Every continuous r(p) — main and per-sg, 2- and 3-horizon
  — is 1×1 with the real joint p at [1,1] (clean case per-sg p's .053/.469; battery
  case 6.2e-06/4.1e-07); the binary sites' r(p) is 1×2 with the vacuous
  0b.targetflag base column at [1,1] and the real p at [1,2], which is the whole
  sibling index convention. Stock output prints "." on the per-sg ame joint lines
  while att twins are real; stock values match the blessed baseline bit-for-bit.

  **REQUIRED PROTOCOL for the fix wave (per the author):** before shipping
  [1,2]→[1,1] at the two sites, produce dumps of the FULL wooldid console output
  AND the relevant matrices (r(p) and r(b) at the sg-continuous joint-test sites,
  captured mid-run via an instrumented scratch build) from AT LEAST 3 different
  example dataset runs, including datasets rich enough to rule out data-input
  degeneracy as the driver of anything observed; then implement, re-bless
  `cont_sg.csv` (sole pinning baseline), and document the resulting old-vs-new
  suite diff explicitly (the comparison alone cannot flag an old-shared defect).
- **T-M3: ignore for now** (author ruling; revisit in the efficiency pass).
- **T-M4: author graduated; IMPLEMENTED** — finaladvtest is now captured beside
  matrixstartgrab while r(b) is authoritative, with the late block reduced to the
  empty/dot fallback. Two-build verification: byte-identical %21x e() dumps and
  caution streams on a normal sg+esf+adversarial config AND a crafted
  fallback-firing config (single in-window-horizon subgroup, firing proven in both
  builds by the one-column e(b_es_att_s2) signature).

## 3. Finder claims refuted during adjudication

- **Permutation count-string "e+ 15" interior-space claim**: mata's `strofreal`
  strips leading blanks (probed: `strofreal(36, "%3.0f")` = "36", and the baseline's
  `perm_npossible_str` = `1.1214e+36`, no space). Stata *display* format semantics
  don't carry over to mata. No defect.

## 4. Deferred observations (next phase's radar, from the finder wave)

Statistical/deep: permute reps strip subgroups() (v1 design) — check the observed
statistic and permuted null come from the same specification; esf rep-parity latency
(the shell feeds the already-incremented esprelength back through validate on rep
lines — inert today, model-side consumers would break parity); Break/rc-1 inside the
rep *loop* is treated as a failed rep rather than propagated (rep-0 propagates; the
ledger's Break-mid-loop note reflects the old reading); unit id 0: a treated unit
with user id 0 gets `_i = 0` and folds into the pooled never-treated factor level
with rc 0 (realistic for 0-indexed ids — deserves a probe + screen decision);
`gegen mint = min(t)` target type under >2^24 time codes; E1b's jointmarginsdof
hardening covers only the two binary whole-sample joint tests (contreat + SG sites
re-create the latent trap if the impvalues/fastnose joint-test gates are ever
relaxed); SG 'cont'-round row-grab by position (assumes sglist order of the over()
post matrix); F2 partial-post asymmetry (b_es ships silently without V_es when only
the V post fails); w_per_it gegen'd float vs wooldid_ptp's double recompute (~1e-7
diagnostic-only drift in hbar_it); NIS bisection leaves a cell silently missing on
non-convergence (unlike the k>4 omission, no note); permexact all-reps-failed corner
would hit `J(0,k,.)` instead of a clean refusal (argued unreachable).

Cosmetic/docs candidates not acted on: orphan-modifier options silently inert
(contreatpoly without contreat; cfxplottypes/cfxplotnose/makecfxplotsbysg without
makecfxplots — the ticketed poisexpresults class; cfxlattice already gets a Caution);
verbose "including event study estimates" prints under espre(0)+espost(0); the
contreat joint-test block re-runs its margins identically up to 4× (inside the
estimand loop; binary twin is outside — performance only); f9tab padding rownames
label row 3 "t" on poisson paths (internal-only); sthlp labels e(joint_horiz_i)
"(ATT_i)" though it's computed on the shared att_it-family curve; the pooled-regime
planted-zero display row isn't covered by the sthlp's reference-period-zeros
sentence; escomp nclusters is silently all-missing under multiway clustering (docs
gap); shared-curve labeling on att_i/att_itime-only runs (NFPTP-03 added the
stand-in note only to the Roth block); part(roth) replay prints nothing where a live
run showed the cap-noi error (already-broken-diagnostics corner); SG joint tests post
as e() macros while full-sample counterparts are scalars (known ledgered
escalation); CITATION.cff is not on the acd9f7e release checklist (its 2023 dates
must be refreshed at actual release); permrun's wdpo1-4 comparator tempnames cap at
four estimands (implicit coupling).

Added post-ruling (from the T-M2 verification runs): the continuous SG-equality
joint tests print "." alongside T-M2's per-sg lines — "all SGs have the same ame"
and "...pre-treatment ame" are missing-valued in BOTH probe cases while their att
twins are real, and the blessed cont_sg baseline pins those scalars as MISSING
(%21x `+1.0…X+3ff` decodes to Stata missing, not 1.0 — longstanding state, not a
regression). The clean case's subgroup dose effects differ by construction, so
rejection, not ".", is the expected outcome there. Distinct site family from T-M2
(the E1 chain/contrast machinery, not the [1,2] read); undiagnosed; deep-phase item.

Standing-ticket correction: `test/suite/REPORT_2026-08-18.md` ticket 4 is now
two-thirds stale — `contreatelasticitytype(notAType)` and `safety(nonsense)` became
hard r198 screens in the VAL-G21 closures; only poisexpresults-without-poisson
remains a silent accept.

## 5. Post-commit adversarial verification

One adversarial diff reviewer attacked all five ado/docs commits with a refute-first
mandate: **5/5 UNREFUTED**. Highlights: the comment-only commit proven by
comment-stripped byte-identity of parent vs commit blobs (no block comments,
continuations, or #delimit exist to smuggle code into `*` lines); the permrun
seeding-loop deletion proven by a full-program read (zero wdpo references before the
rep-0 overwrite; every intervening exit is a hard error; overwrite-failure exits fire
before any read); the r()-state attack on the wdpnclus deletion found every r() read
downstream dominated by its own setter; the contreat message gate proven by
brace-depth tracking to be literally the enclosing block's own condition; the esplots
insert byte-identical to its production sibling after whitespace strip, with `min`'s
only downstream consumer the twoway xlabel() option. Residual observation (pre-existing,
harmless direction, left as-is): the sthlp reserved-name sentence says the scratch
estimate is removed "on every run, including runs that end in an error", but a run
failing before the store neither overwrites nor removes a pre-existing user estimate
(the shell drop is flag-gated on the store having happened) — the warning overstates
only in the direction that makes users more cautious.

## 6. Section verdicts (condensed)

Every lens returned "clean after 16 tranches" modulo the items above: shell/validate
syntax-line parity and the c_local closure verified complete; the four estimation
blocks' E1/F9/NFADV machinery matches comments and siblings; wooldid_advstat's rclass
contract verified spelled identically at all four call sites with correct
gapskip/degenerate handling; permute screens match engine contracts; display/replay
parity holds by construction; the post-7/28 docs surface was clean apart from the two
items fixed; all 14 programs brace-balanced; no duplicated adjacent statements; r()
discipline solid except T-M4.
