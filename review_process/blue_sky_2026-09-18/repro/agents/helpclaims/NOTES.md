# helpclaims recon — NOTES

Task: build the HELP-FILE TESTABLE-CLAIM EXTRACT at
/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/help_claims.md

## State
- Read the WHOLE sthlp (1664 lines) in 9 chunks; de-SMCL'd copies kept as chunk2..chunk9.txt here.
  (chunk1 = lines 1-260, saved by the harness at
   /Users/hegland/.claude/projects/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/tool-results/bo32s93pz.txt)
- Cross-checked a handful of claims against the pinned .ado (read-only greps only). No Stata run yet.

## Confirmed-by-code findings so far (ado line refs on the pinned build)
1. **e(b) subgroup names are an INFIX `_SG<value>_`, not the suffix `_s#` the help claims.**
   ado:4015 `local rownamespool_sginj = subinstr("`rownamessg`sgcase''",":","_SG`sgcase'_",.)`
   plus ado:10910, 11061, 14638, 16941 all build `ES_<base>_SG<sg>_<label>`.
   sthlp:1307 says "Under subgroups() each subgroup's rows take the suffix _s#." FALSE.
2. **No ES_att_i_* / ES_att_itime_* rows exist**; att_i/att_itime event-study rows are named
   ES_att_it_*. ado:3556-3566 + ado comment 6828-6832. sthlp:1307 implies they exist;
   sthlp:848-851 and sthlp:1311 say they do not. Internal contradiction.
3. e(b) names are FLAT (colon -> underscore at ado:4006); e(wdidmainresults) keeps the colon
   equation form. Both stated correctly in the help (sthlp:1303-1307, 1311).
4. **fweights**: the worker's syntax accepts `[aw fw pw /]` (ado:555) but explicitly refuses a
   user fweight at ado:628-630 ("wooldid does not accept fweights"), demoting fw->aw only in
   replication mode. The help NEVER mentions fweights (Syntax line sthlp:57 says `[aw/pw=weights]`).
5. `set emptycells drop` at ado:602; shell snapshots c(emptycells) at ado:62 and restores on the
   error path at ado:84. Program-scoped otherwise.
6. `set level` IS mutated for the duration of a call and restored by the shell via
   $WOOLDID_LEVEL0 (ado:~66-70). Undocumented in the help.
7. simultaneousbands analytic engine: caller rngstate saved ado:10892, restored ado:10963;
   the "fixed internal seed" is actually per-family, 20260907 + hash(family name) mod 1e6
   (ado:10936-10937).
8. e(escomp) wshare_att / wshare_att_it are normalised over MATCHED horizons only (ado:3888-3893),
   so "each column sums to 1" is by construction — but a horizon with no matching raw row keeps a
   zero row.

## Next
- Write help_claims.md (the deliverable). Rank by where bodies are buried.

## DELIVERABLE WRITTEN (complete)
/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/help_claims.md
1806 lines / 200K. Structure:
  §0 how to use; §1 twenty ranked first checks (H-01..H-20); §2 full extract by help section (2.1-2.18);
  §3 ten internal contradictions C1-C10 with both line numbers; §4 twelve silences S1-S12;
  §5 name-length arithmetic table; App A measured e() inventory; App B what I did not cover.

## Empirical confirmations added after the code reads (recon1.do / recon1.log)
- Stata 19.0 MP. hdidregress / csdid / did_imputation ALL present (which rc 0).
- e(b) subgroup names measured: main_SG1_att, ES_att_SG1_pre3, pretreatment_SG1_att_i. No "_s1" suffix. H-01 PROVEN.
- att_i run colnames measured in full: main_att_i pretreatment_att_i ES_att_it_pre3 ES_att_it_pre2
  ES_att_it_post0 ES_att_it_post1. No ES_att_i_*. H-02 PROVEN.
- fw -> r(198) bespoke refusal; iw -> r(101). Help never mentions fweight at all (S1).
- Clean on a plain run: count if e(sample)==e(N)==360; c(level)/c(emptycells)/c(rngstate)/sortedby/
  datasignature all unchanged; rerun bit-identical e(b).
- e(escomp) 4x8 vs 8 ES rows in e(wdidmainresults) 12x7  => the "one-for-one" alignment claim is 1:2 (H-16).
- e(wdid_vblockn)=36 == 6 families x C(4,2). e(wdid_espre)="3" when esprelength(2) requested.
- Undocumented behaviour-changing globals found in code: WOOLDID_BOOTNOSUB (suppresses the boottstat
  cluster substitution => reintroduces the 2-1/G ~ sqrt(2) SE inflation, silently), WOOLDID_BOOTDEBUG.
  ado:9367-9378. Help names only WOOLDID_LASTRC. -> S3.
- ado:9364-9366 (dated 2026-09-16): the P7 consistency assert runs only "on every replication that posts
  e(N_clust)" -> unchecked paths exist. -> S4.

## Status: COMPLETE. Nothing left half-written.

## Post-advisor revision (final)
- grep -c '_SG' wooldid.sthlp == 0 (the string appears NOWHERE in the help); '_s#' appears 5x. H-01 airtight.
- Added C11: the shipped Examples (0, 7, 11, 12, 13) all use cluster(i), the exact configuration
  sthlp:67/255/1331 document as undercovering and as invalidating every e(V_es_*) consumer. Doc-side twin
  of the confirmed 0.79-0.88 coverage failure.
- Added S13: sthlp:2 version stamp "1.4 2026-08-19" is stale (body documents 09-05..09-08 and permonsetcohorts 09-16).
- Added literal wooldid commands to H-03, H-04, H-05, H-06, H-10, H-11, H-13.
- Recalibrated H-16 down to LOW (the 4-vs-8 count is certain but "the ES rows" is a reading).
- Promoted the permonsetcohorts rep-0-comparator check with an explicit command.
File is final at 1868 lines.
