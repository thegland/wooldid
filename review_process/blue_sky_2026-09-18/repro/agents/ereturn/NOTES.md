# ereturn recon — NOTES (live, 2026-09-16)

## DELIVERABLE
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/ereturn_map.md` — written, being revised.

## FILES HERE
- `mkdata.doi` (80 units x 10 periods panel) -> `wd_base.dta`
- `dump.doi` -> `edump <tag> <rc>` writes the full e() contract of the current results to `edump.txt`
- `probeA.do` 28 analytic configs | `probeB.do` 18 event-study configs | `probeC.do` contreat+name
  stress+staleness+sentinel | `probeD.do` boot/perm | `probeE.do` sentinel row dump + consumers |
  `probeF.do` permute + BLANKING-GUARD AUDIT (G01/G02/G03) + bootstrap reps(8) rank (H01)
- `edump.txt` consolidated empirical surface; `probe?.log` full logs

## STATE
- probeA/B/C/D/E: DONE. probeF: RUNNING (G01 plain done; G02 boot, G03 perm, H01 pending).
- permute configs D06-D10 and F01-F04 ALL FAILED r(459): permute needs cluster() to nest i AND
  clusters to be subsets of treatment-onset groups. My `clu = ceil(id/4)` mixes onsets.
  FIX FOR A SUCCESSOR: build cohorts in contiguous id blocks (grp = floor((id-1)/20)) and set
  clu = grp*10 + ceil(mod(id-1,20)/4). NOT YET RUN -> §4b of the map is static-only.

## INCIDENT
`probeE.do` ran `wooldid , update`, which performed a LIVE net-install check of the user's installed
wooldid. It reported "all files already exist and are up to date" (no change made). Do not repeat.

## DE-DUP against repo records (done)
- **-1e100 sentinel: KNOWN AND RULED.** wooldid_solved_issues.md:262 + review_process/
  held_for_ruling_2026-09-09.md:11,92-108 + rulings_2026-09-07.md:128. Verified by the author on
  ereturn display / esttab / lincom / nlcom / test / sqrt. NOT verified on matrix-level consumers.
- **8-digit subgroup name crash: KNOWN AND RULED** as caution-only, wooldid_solved_issues.md:180.
- **Silent colname padding: RESOLVED INSTANCE ON RECORD**, wooldid_solved_issues.md:338 — Stata pads
  a SHORT namelist by repeating the last name, rc 0, no warning; `_b[ES_att_SG0_post1]` returned the
  post0 estimate. THE MECHANISM IS THE TEMPLATE; hunt siblings at unguarded colname sites.
- **Analytic joint/adversarial p-values surviving boot/permute: ON RECORD as partly not_covered**
  (bootstrap_coverage_2026-09-02.md:158 "analytic (cautioned)"; permute_coverage_2026-09-02.md:147,150
  "not_covered"). Those records PREDATE the 2026-09-05 bootjnt/permjnt work -> re-measure (probeF G).
- **No record found** for: reghdfe's e() left behind after a failed run; `macro drop` defeating the
  replay guard; `e(clustvar)` absent under `customvce(cluster ...)`; `e(roth_slope_ul_*)` undocumented;
  the help's `_s#` subgroup-suffix claim vs the actual `_SG#_` infix.
- Help DOES warn about cross-estimate `test`/`lincom` (sthlp 1408/1410/1468) -> that is a UX finding
  (warning is buried in the e(wdid_vcvstate) synopt), not doc-vs-actual.

## FINAL RESULTS (all probes complete except probeG P04)
- **F7 (new, CRITICAL)**: bootstrap `reps(8)` on K=20 -> `rank(e(V)) = 7`, posted as a full VCE with
  `e(properties)= b V`. `test` returns p = 1.46e-14 (df 3) / 9.9e-28 (df 7) / 5.4e-10 (df 2) on two
  independent panels. Only caution printed is about the p-value GRID. `reps()` has no floor by author
  ruling (6053-6056); `permreps()` does (min 20, ado:5890).
- **F8 CLOSED both engines**: bootstrap AND permute rewire every joint/adversarial p to a k/(R+1)
  share; adversarial is blanked to missing under bootstrap. The 2026-09-02 coverage records are stale.
  RESIDUE: the plain run blanks all three joint_*_allsgssame (sgdeg rank guard) while the permute run
  reports finite equality p's for the same contrasts -> does permute consume the degenerate analytic SE?
- **permute default surface measured** (P01): no e(V), no b_es_*/V_es_*/numpre_* (f2bonly carve-out at
  ado:1990 excludes permute from the fastnose ruling), e(vce)/e(vcetype) blanked but e(clustvar)/
  e(N_clust) kept, e(permdist) is 40x4 (estimands) while e(permdist_b) is 40x20 (all coefficients).
- Trap that bit me: after a failed resampled run, F1 leaves the previous run's e(); my first
  blanking comparison (probeF GPERM lines) read the bootstrap run's values. Always print rc + cmdline.
