
*! version 1.4 2026-08-19  wooldid
* Author: Thomas Hegland, Agency for Healthcare Research and Quality
* Contact: thomas.hegland@ahrq.hhs.gov, @thomas_hegland (twitter), thomashegland.com
* This code file is the work of the author. The code comes with no endorsement, guarantee, or warranty
* from the Agency for Healthcare Research and Quality, the Department of Health
* and Human Services, or the United States government. Any views expressed within
* the code file are also those of the author, and not necessarily that of any of
* the aforementioned governmental entities.

cap prog drop wooldid
* properties(svyr) is a DISPLAY-ROUTING declaration, not an svy-support claim: Stata's
* bootstrap/jackknife prefixes route their post-replication output through
* _prefix_display only when the command's properties intersect {svyb,svyj,svyr};
* otherwise they attempt a bare `wooldid, <diopts>' replay, which parses as a
* too-few-variables r(198) here and aborts the prefix AFTER a fully successful
* replication run. svyr is the narrowest property that fixes the routing. svy itself
* remains unsupported and fails loudly (measured: svy: wooldid -> r(101)).
prog def wooldid, eclass properties(svyr)
version 16.0
  * F11: replay - wooldid typed alone after a successful wooldid estimation redisplays
  * the stored results from e() (tables and diagnostics; plots are not redrawn, and no
  * options are accepted on replay). gated on a STRICTLY EMPTY command line, so every
  * other invocation - including "wooldid, update" and an option typo like "wooldid,
  * att" - falls through to the worker and behaves exactly as before
  if `"`0'"' == "" & "`e(cmd)'" == "wooldid" {
    * refuse to replay across a failed wooldid call: WOOLDID_LASTRC is a session
    * global (namespaced like WOOLDID_CLEANUP) set to the failing rc in the shell's
    * rc-gated epilogue below whenever a wooldid invocation fails - at the parse
    * stage or anywhere later - and cleared there on any successful completion.
    * (one deliberate exception, per NFRPL-01: a post-estimation permutation-phase
    * failure restores the still-valid observed results and clears the marker on its
    * way out, so those results remain replayable.)
    * while it is set, the stored e() results predate the user's latest (failed)
    * attempt, so redisplaying them silently would be misleading. exiting here
    * leaves the marker set: replay stays refused until a wooldid call completes
    * successfully. shell-level state only: never posted to e(), invisible to the
    * S5 c_local closure, and this check lives strictly inside the empty-line
    * branch, so "wooldid, update" and every other non-empty call is untouched
    if "${WOOLDID_LASTRC}" != "" {
      if "${WOOLDID_LASTRC}" == "postloop" dis as error "Error: the last wooldid call failed after its replication loop, while rewiring the stored results to the permutation or bootstrap inference; the stored results are partially rewired. Refit wooldid before replaying."
      else dis as error "Error: the last wooldid call failed with r(${WOOLDID_LASTRC}); the stored results predate that failed call. Refit wooldid before replaying."
      exit 198
    }
    wooldid_display , part(1)
    wooldid_display , part(roth)
    wooldid_display , part(2)
    if "`e(perm_method)'" != "" dis as result " "
    wooldid_display , part(perm)
    if "`e(boot_method)'" != "" dis as result " "
    wooldid_display , part(boot)
    dis as result  " "
    exit
  }
  * thin shell (B18/B22): run the worker under capture so that when a run fails midway,
  * the caller's frame is restored and this run's scratch estimates are removed before
  * the error is re-raised. the epilogue is rc-gated: successful runs are untouched by
  * construction. the worker registers what it actually created this run in the
  * WOOLDID_CLEANUP global (keep / fullmodel / ppml) so the epilogue never removes a
  * user-owned object the worker did not create. e() posted by the worker survives the
  * shell; rc and error messages pass through exactly (exit, not error - error reprints)
  local wdsh_ec "`c(emptycells)'"
  qui frame
  local wdsh_frame `r(currentframe)'
  cap macro drop WOOLDID_CLEANUP
  capture noisily wooldid_worker `0'
  local wdsh_rc = _rc
  * F054: the worker sets level() for the duration of the call; restore the caller's
  if "$WOOLDID_LEVEL0" != "" {
    set level $WOOLDID_LEVEL0
    macro drop WOOLDID_LEVEL0
  }
  nobreak {
    local wdsh_flags "$WOOLDID_CLEANUP"
    cap macro drop WOOLDID_CLEANUP
    * replay-staleness marker consumed by the replay branch above: any failed
    * invocation (parse-stage or later) records its rc; any successful completion
    * (estimation or update) clears it (assigning an empty value drops a global)
    if `wdsh_rc' != 0 global WOOLDID_LASTRC `wdsh_rc'
    else global WOOLDID_LASTRC
    if `wdsh_rc' != 0 {
      * set emptycells is program-scoped (stata restores it automatically at program
      * exit); this restore is belt-and-braces against version drift only
      if "`c(emptycells)'" != "`wdsh_ec'" set emptycells `wdsh_ec'
      * stata auto-drops the worker's tempname frames on error but then leaves the
      * user in the default frame, not the frame wooldid was called from
      qui frame
      if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
      if strpos("`wdsh_flags'", "keep") == 0 {
        if strpos("`wdsh_flags'", "fullmodel") != 0 cap est drop ___wooldidfullmodel
        if strpos("`wdsh_flags'", "repmodel") != 0 cap est drop ___wooldidrepmodel
        * frame-qualified: if the failure destroyed the snapshot frame, the restore above
        * cap-fails and an unqualified drop would act on whatever frame is current -
        * deleting an unrelated variable of this name there. qualified, it cap-fails
        * benignly instead
        if strpos("`wdsh_flags'", "ppml") != 0 cap frame `wdsh_frame': drop _ppmlhdfe_d
      }
    }
  }
  if `wdsh_rc' != 0 exit `wdsh_rc'

  * F6: permutation inference - runs only after a fully successful estimation. the
  * worker set wdp_on and (since 2026-09-05) deferred its table display and plots, as it
  * does under bootstrap(); re-running wooldid_validate materializes every normalized
  * option here (the S5 c_local blanket), from which the rep command line is built.
  * the observed results are stored, the rep loop re-runs wooldid quietly on permuted
  * onset assignments accumulating the FULL rep e(b) name-aligned (wooldid_permrun),
  * the observed results are restored, every reported inference surface is rewired to
  * the permutation distribution (wooldid_permwire, wooldid_permjnt, wooldid_permcellsx),
  * and the (e()-driven) display prints the rewired tables, the deferred plots and the
  * part(perm) footer. Author rulings 2026-09-05: only permutation p-values are ever
  * shown under permute(); joint, equality and adversarial tests are permutation tests
  * of a fixed-metric statistic; the intervals are shift intervals; the observed run is
  * a fast (no-SE) run unless permtstat or permanalyticse
  if "`wdp_on'" == "1" {
    qui wooldid_validate `0'
    * permtci implies permtstat (resolved in validate as well; kept for the rep string)
    if "`permtci'" != "" & "`permtstat'" == "" local permtstat permtstat
    * rep option string: keep everything that shapes the point estimates - INCLUDING
    * subgroups (their rows are permuted too, 2026-09-05) and the cfx plot requests
    * (their cells are harvested by key); strip display/diagnostic/plot work and the
    * permutation options themselves; reps run under fastnose (only coefficients are
    * consumed) UNLESS permtstat, in which case each replication also produces a VCE so
    * its own standard errors can form its t-ratios
    local wdp_ro `att' `att_it' `att_itime' `att_i' `esfixedbaseperiod' `timetrends' `poisson' `poisexpresults' `robust' `ccnoabsorb' `latticeignoreweights' `emptycellsoverride'
    if "`permtstat'" == "" local wdp_ro `wdp_ro' fastnose
    if "`cluster'" != "" local wdp_ro `wdp_ro' cluster(`cluster')
    if `"`customvce'"' != "" local wdp_ro `wdp_ro' customvce(`customvce')
    if "`fe'" != "" local wdp_ro `wdp_ro' fe(`fe')
    if `"`coarsecohortcontrols'"' != "" local wdp_ro `wdp_ro' coarsecohortcontrols(`coarsecohortcontrols')
    if "`ccc_absorb'" != "" local wdp_ro `wdp_ro' ccc_absorb(`ccc_absorb')
    if "`controls'" != "" local wdp_ro `wdp_ro' controls(`controls')
    if "`interactivecontrols'" != "" local wdp_ro `wdp_ro' interactivecontrols(`interactivecontrols')
    if "`customweightmultiplier'" != "" local wdp_ro `wdp_ro' customweightmultiplier(`customweightmultiplier')
    * a cell-collapse fallback in the observed run (e(wdid_fcguard), see wooldid_fcdeg)
    * extends to the reps, as in the bootstrap block. NB the screen reads a subgroup
    * margins r(V); under the fast-mode default only the subgroup-by-horizon call
    * carries one (BS-26), so the rule bites there and under permtstat/permanalyticse
    local wdp_fb = cond("`e(wdid_fcguard)'" == "fallback", "fallback", "")
    if "`fastcollapse'" != "" & "`wdp_fb'" == "" local wdp_ro `wdp_ro' fastcollapse
    if "`fullprecisionslow'`wdp_fb'" != "" local wdp_ro `wdp_ro' fullprecisionslow
    if "`contreat'" != "" local wdp_ro `wdp_ro' contreat(`contreat') contreatpoly(`contreatpoly')
    if `"`contreatcontrols'"' != "" local wdp_ro `wdp_ro' contreatcontrols(`contreatcontrols')
    if `"`lattice'"' != "" local wdp_ro `wdp_ro' lattice(`lattice')
    if "`makecfxplots'" != "" {
      local wdp_ro `wdp_ro' makecfxplots
      if "`cfxlattice'" != "" local wdp_ro `wdp_ro' cfxlattice(`cfxlattice')
      if "`cfxplottypes'" != "" local wdp_ro `wdp_ro' cfxplottypes(`cfxplottypes')
      if "`makecfxplotsbysg'" != "" local wdp_ro `wdp_ro' makecfxplotsbysg
      if "`cfxplotnose'" != "" local wdp_ro `wdp_ro' cfxplotnose
    }
    if "`semielasticity'" == "eydx" local wdp_ro `wdp_ro' semielasticity
    if "`contreat'" != "" & "`contreatelasticitytype'" != "" local wdp_ro `wdp_ro' contreatelasticitytype(`contreatelasticitytype')
    if "`subgroups'" != "" local wdp_ro `wdp_ro' subgroups(`subgroups')
    * the RAW user esprelength (validate's esf normalization re-fires on every nested
    * re-parse) and the user's OWN safety value, as the bootstrap block emits them
    if `esprelengthraw' != 0 local wdp_ro `wdp_ro' esprelength(`esprelengthraw')
    if "`esposttrip'" == "no" local wdp_ro `wdp_ro' espostlength(-999)
    else if `espostlength' != 0 local wdp_ro `wdp_ro' espostlength(`espostlength')
    local wdp_ro `wdp_ro' esrelativeto(`esrelativeto')
    if "`safetyinput'" != "" local wdp_ro `wdp_ro' safety(`safetyinput')
    if `regtol' != 9 local wdp_ro `wdp_ro' regtol(`regtol')
    if `oostreatedcontrols' != 0 local wdp_ro `wdp_ro' oostreatedcontrols(`oostreatedcontrols')
    if `anticipation' != 0 local wdp_ro `wdp_ro' anticipation(`anticipation')
    local wdp_est `att' `att_it' `att_i' `att_itime'
    if "`wdp_est'" == "" local wdp_est att
    local wdp_ts = ("`permtstat'" != "")
    tempname wdpobs wdpdist wdpcomp wdpobsb wdpacc wdpsacc wdpcmpb wdpcmpse wdpcell wdpcell0 wdpasg
    qui est sto `wdpobs'
    matrix `wdpobsb' = e(b)
    cap noi wooldid_permrun , y(`y') i(`i') t(`t') ttrevar(`ttre_asgiven') clus(`cluster') permute(`permute') permreps(`permreps') `permexact' permstrata(`permstrata') permprobs(`permprobs') `permtstat' repopts(`wdp_ro') estimands(`wdp_est') wtype(`wdp_wt') wexp(`wdp_wx') ifr(`wdp_if') distname(`wdpdist') compname(`wdpcomp') obsb(`wdpobsb') accname(`wdpacc') seaccname(`=cond(`wdp_ts', "`wdpsacc'", "")') cmpbname(`wdpcmpb') cmpsename(`=cond(`wdp_ts', "`wdpcmpse'", "")') cellacc(`=cond("`makecfxplots'" != "", "`wdpcell'", "")') cell0(`wdpcell0') dots(`dots') rseed(`rseed') asgname(`=cond("`permsave'" != "", "`wdpasg'", "")')
    local wdp_rc = _rc
    qui est restore `wdpobs'
    cap est drop `wdpobs'
    global WOOLDID_LASTRC
    if `wdp_rc' != 0 exit `wdp_rc'
    * from here until the display completes, e() is being rewired to the permutation
    * inference: an uncaught failure in between leaves a half-rewired e() that the replay
    * branch must refuse (the marker is consumed there and cleared at the end of this
    * block; a plot failure records its rc the same way). Added 2026-09-05
    global WOOLDID_LASTRC postloop
    * run-level permutation results on the restored observed results
    ereturn local perm_method `permute'
    if "`rseed'" != "" ereturn local perm_seed `rseed'
    ereturn local perm_rngstate `wdp_r_rng0'
    ereturn hidden local wdid_permest `wdp_est'
    ereturn scalar perm_reps = `wdp_r_reps'
    ereturn scalar perm_completed = `wdp_r_done'
    ereturn scalar perm_failed = `wdp_r_fail'
    ereturn scalar perm_npossible = `wdp_r_npd'
    ereturn local perm_npossible_str `wdp_r_nps'
    ereturn local perm_stat = cond(`wdp_ts', "t", "b")
    local wdp_M = cond("`permexact'" != "", "`wdp_r_reps'", ".")
    * BS-23 (author ruling 2026-09-07): the kernels switch the share convention on
    * MEMBERSHIP - is the observed statistic one of the reference values - rather than on
    * completeness. 0 when there is no observed row to point at.
    local wdp_OR = cond("`permexact'" != "" & "`wdp_r_obsrow'" != "", "`wdp_r_obsrow'", "0")
    * rewire the coefficient tables and post the stripe-aligned surfaces (eclass sub)
    wooldid_permwire , acc(`wdpacc') obsb(`wdpobsb') cmpb(`wdpcmpb') seacc(`wdpsacc') cmpse(`wdpcmpse') level(`level') m(`wdp_M') orow(`wdp_OR') done(`wdp_r_done') fail(`wdp_r_fail') blank(`wdp_r_blank') ts(`wdp_ts') method(`permute') reps(`wdp_r_reps') ase(`permanalyticse') clus(`cluster') gt(`gteffects') absp(`permabsp') side(`permside') rawp(`permrawp')
    * the legacy per-estimand scalars, read from the posted kernel by name so they are
    * the SAME numbers the tables carry (one convention: equal-tailed p, 5th-95th null
    * band, permutation SD; on the t-ratio scale under permtstat)
    tempname wdpk
    matrix `wdpk' = e(perm_kernel)
    foreach wdpe in `wdp_est' {
      local wdpc = colnumb(`wdpk', "main_`wdpe'")
      if `wdpc' < . {
        ereturn scalar perm_p_`wdpe' = `wdpk'[3, `wdpc']
        * CS-57 (author ruling 2026-09-08): rows 11-12 are the signed level() quantiles of
        * the NULL DISTRIBUTION - the region the statistic occupies under the null - not an
        * interval for the estimand. Naming them lb/ub made one word mean two things in a
        * single run, since e(perm_ci) and every e(boot_ci_*) use ll/ul for a genuine
        * interval. The band-specific names are canonical from here; the legacy pair stays
        * for compatibility and is documented as deprecated.
        ereturn scalar perm_bandlo_`wdpe' = `wdpk'[11, `wdpc']
        ereturn scalar perm_bandhi_`wdpe' = `wdpk'[12, `wdpc']
        ereturn scalar perm_lb_`wdpe' = `wdpk'[11, `wdpc']
        ereturn scalar perm_ub_`wdpe' = `wdpk'[12, `wdpc']
        ereturn scalar perm_sd_`wdpe' = `wdpk'[6, `wdpc']
        if "`permtci'" != "" {
          * inverted-t CI (permtci, retained): observed b +/- (level() quantile of the
          * permuted |t|) x observed se. The shell's level is 0-means-95: resolve it here
          local wdptb = _b[main_`wdpe']
          local wdpts = _se[main_`wdpe']
          local wdptclv = cond(`level' == 0, c(level), `level')
          tempname wdpq
          mata: st_numscalar("`wdpq'", wooldid_permq(sort(abs(select(st_matrix("`wdpacc'")[., `wdpc'] :/ st_matrix("`wdpsacc'")[., `wdpc'], (st_matrix("`wdpacc'")[., `wdpc'] :< .) :& (st_matrix("`wdpsacc'")[., `wdpc'] :> 0) :& (st_matrix("`wdpsacc'")[., `wdpc'] :< .))), 1), `wdptclv'/100))
          if `wdpts' < . & `wdpts' > 0 & `wdpq' < . {
            ereturn scalar perm_tlb_`wdpe' = `wdptb' - `wdpq'*`wdpts'
            ereturn scalar perm_tub_`wdpe' = `wdptb' + `wdpq'*`wdpts'
          }
          else dis as text "Caution: the inverted-t CI for `wdpe' is unavailable (no usable observed standard error)."
        }
      }
    }
    * CS-58 (author ruling 2026-09-08): e(permdist) is NOT the twin of e(bootdist). It
    * holds one column per requested estimand on the e(perm_stat) scale (the t-ratio under
    * permtstat), where e(bootdist) is the full coefficient-draw matrix stripe-aligned to
    * e(b). The true twin is e(permdist_b), posted by wooldid_permwire. e(permdist) is
    * retained and documented as deprecated in favour of it.
    if `wdp_r_done' <= 10000 cap ereturn matrix permdist = `wdpdist'
    * NFPERM-12 (lean kernel): the estimand comparators actually used for the legacy
    * scalars, on the e(perm_stat) scale (the full vector is e(perm_comparator_b))
    cap ereturn matrix perm_comparator = `wdpcomp'
    * joint tests, subgroup equality, adversarial and contrasts (eclass sub)
    wooldid_permjnt , acc(`wdpacc') obsb(`wdpobsb') cmpb(`wdpcmpb') level(`level') m(`wdp_M') orow(`wdp_OR') ts(`wdp_ts') jt(`jointtests') adv(`adversarial') advalt(`altmethod_adversarial') advesrel(`esrelativeto') advtt(`timetrends') esf(`esfixedbaseperiod') absp(`permabsp') side(`permside') rawp(`permrawp') sd(`permstepdown')
    * the continuous-treatment plot cells (eclass sub)
    if "`makecfxplots'" != "" & `wdp_ncfx' > 0 {
      local wdp_cfxnames
      forval wdk = 1/`wdp_ncfx' {
        local wdp_cfxnames `wdp_cfxnames' `wdb_cfxname`wdk''
      }
      wooldid_permcellsx , cell(`wdpcell') cell0(`wdpcell0') ncfx(`wdp_ncfx') names(`wdp_cfxnames') level(`level') m(`wdp_M') orow(`wdp_OR') ts(`wdp_ts') lattype(`latticetype') cfxlat(`cfxlattice') absp(`permabsp') side(`permside') rawp(`permrawp')
    }
    * permsave(): the replication draws as a .dta in the working directory
    if "`permsave'" != "" wooldid_permsave , file(`permsave') acc(`wdpacc') obsb(`wdpobsb') seacc(`=cond(`wdp_ts', "`wdpsacc'", "")') aseobs(`=cond(`wdp_ts', "e(perm_ase_obs)", "")') asg(`wdpasg') clusvals(`wdp_clusvals') obsons(`wdp_obsons') seed(`rseed') rng(`wdp_r_rng0') `replace'
    * deferred display: the worker skipped these under permute(); they read only e(),
    * which now carries the permutation inference (identical formatting, free replay)
    wooldid_display , part(1)
    if "`makecfxplots'" != "" {
      local wdcfopts
      if "`saveplots'" != "" local wdcfopts saveplots(`saveplots')
      cap noi wooldid_cfxplots , `wdcfopts' `replace'
      local wdsh_prc = _rc
      if `wdsh_prc' != 0 {
        qui frame
        if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
        global WOOLDID_LASTRC `wdsh_prc'
        exit `wdsh_prc'
      }
    }
    wooldid_display , part(2)
    * event-study plots, deferred from the worker so they read the REWIRED tables: their
    * bands are the permutation shift intervals the table prints. A plot failure is an
    * error of this call (its rc propagates, the caller's frame is restored and the
    * replay marker is set), as it is on a plain run
    * CS-87 (consistency sweep 2026-09-06): the worker draws the event-study plot only when
    * horizons exist (esprelength/espostlength > 0 after the panel clips them); the deferred
    * engine call must apply the same guard, else wooldid_esplots runs on an empty frame,
    * dies r(2000) and loses the whole run. The horizons are visible here as ES_ columns of
    * the (rewired) e(b)
    if "`makeplots'" == "makeplots" & strpos("`: colfullnames e(b)'", "ES_") > 0 {
      local wdesopts
      if "`saveplots'" != "" local wdesopts saveplots(`saveplots')
      cap noi wooldid_esplots , `wdesopts' `replace'
      local wdsh_prc = _rc
      if `wdsh_prc' != 0 {
        qui frame
        if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
        global WOOLDID_LASTRC `wdsh_prc'
        exit `wdsh_prc'
      }
    }
    wooldid_display , part(perm)
    * the one surface still drawn from the analytic covariance (mirrors the bootstrap block)
    if "`histogramcohorteffects'" != "" dis as text "Caution: histogram outputs are drawn from the stage-1 ANALYTIC covariance during estimation; the reported tables, the cfx matrices, the event-study and continuous-treatment plots and the e() inference above are permutation-based."
    dis as result " "
    global WOOLDID_LASTRC
  }

  * F-boot: bootstrap-SE orchestration - runs only after a fully successful estimation.
  * the worker set wdb_on and deferred its table display; the rep loop lives in
  * wooldid_bootrun and the e() rewiring in wooldid_bootwire (compiled-size reasons).
  * the observed results are stored, the rep loop runs nested fastnose wooldid calls
  * under reweighted [aw/pw], the observed results are restored, the reported
  * inference surfaces are rewired to the bootstrap covariance, and the (e()-driven)
  * display prints the rewired tables plus part(boot)
  if "`wdb_on'" == "1" {
    qui wooldid_validate `0'
    * rep option string: keep everything that shapes point estimates - INCLUDING
    * subgroups (their estimands get bootstrapped) and contreatwithin (R2 finding 5:
    * omitting it silently runs the non-within model in reps because validate's
    * internal rewrite emits safety(off)); strip display/plot/diagnostic work and the
    * bootstrap options themselves. reps run under fastnose (only e(b) is consumed) UNLESS
    * boottstat is requested, in
    * which case each replication also has to produce a VCE so its own standard errors can
    * studentize its coefficients
    local wdb_ro `att' `att_it' `att_itime' `att_i' `esfixedbaseperiod' `timetrends' `poisson' `poisexpresults' `robust' `ccnoabsorb' `latticeignoreweights' `emptycellsoverride' `contreatwithin'
    if "`boottstat'" == "" local wdb_ro `wdb_ro' fastnose
    if "`cluster'" != "" local wdb_ro `wdb_ro' cluster(`cluster')
    if "`fe'" != "" local wdb_ro `wdb_ro' fe(`fe')
    if `"`coarsecohortcontrols'"' != "" local wdb_ro `wdb_ro' coarsecohortcontrols(`coarsecohortcontrols')
    if "`ccc_absorb'" != "" local wdb_ro `wdb_ro' ccc_absorb(`ccc_absorb')
    if "`controls'" != "" local wdb_ro `wdb_ro' controls(`controls')
    if "`interactivecontrols'" != "" local wdb_ro `wdb_ro' interactivecontrols(`interactivecontrols')
    if "`customweightmultiplier'" != "" local wdb_ro `wdb_ro' customweightmultiplier(`customweightmultiplier')
    * a cell-collapse fallback in the observed run (e(wdid_fcguard), see wooldid_fcdeg)
    * extends to the reps: under fastnose the screen has no r(V) to read, so they must be
    * sent to the full data outright or the bootstrap covariance of a degenerate subgroup
    * contrast is built from collapsed replicates while the observed inference is not
    * (2026-09-01 probe, A4)
    local wdb_fb = cond("`e(wdid_fcguard)'" == "fallback", "fallback", "")
    if "`fastcollapse'" != "" & "`wdb_fb'" == "" local wdb_ro `wdb_ro' fastcollapse
    if "`fullprecisionslow'`wdb_fb'" != "" local wdb_ro `wdb_ro' fullprecisionslow
    if "`contreat'" != "" local wdb_ro `wdb_ro' contreat(`contreat') contreatpoly(`contreatpoly')
    if `"`contreatcontrols'"' != "" local wdb_ro `wdb_ro' contreatcontrols(`contreatcontrols')
    if `"`lattice'"' != "" local wdb_ro `wdb_ro' lattice(`lattice')
    * cell-level surfaces the author ruled must be bootstrapped, not blanked (2026-09-04):
    * the replications compute them too, and the loop harvests them by cell key
    if "`gteffects'" != "" local wdb_ro `wdb_ro' gteffects
    if "`makecfxplots'" != "" {
      local wdb_ro `wdb_ro' makecfxplots
      if "`cfxlattice'" != "" local wdb_ro `wdb_ro' cfxlattice(`cfxlattice')
      if "`cfxplottypes'" != "" local wdb_ro `wdb_ro' cfxplottypes(`cfxplottypes')
      if "`makecfxplotsbysg'" != "" local wdb_ro `wdb_ro' makecfxplotsbysg
      if "`cfxplotnose'" != "" local wdb_ro `wdb_ro' cfxplotnose
    }
    if "`semielasticity'" == "eydx" local wdb_ro `wdb_ro' semielasticity
    if "`contreat'" != "" & "`contreatelasticitytype'" != "" local wdb_ro `wdb_ro' contreatelasticitytype(`contreatelasticitytype')
    if "`subgroups'" != "" local wdb_ro `wdb_ro' subgroups(`subgroups')
    * A8: emit the RAW user esprelength - validate's esf normalization re-fires on
    * every nested re-parse, so embedding the mutated value widens the window per rep
    if `esprelengthraw' != 0 local wdb_ro `wdb_ro' esprelength(`esprelengthraw')
    if "`esposttrip'" == "no" local wdb_ro `wdb_ro' espostlength(-999)
    else if `espostlength' != 0 local wdb_ro `wdb_ro' espostlength(`espostlength')
    local wdb_ro `wdb_ro' esrelativeto(`esrelativeto')
    * emit the user's OWN safety value (the post-validate working copy may have been
    * rewritten to off by contreatwithin; the rep re-applies the same rewrite itself)
    if "`safetyinput'" != "" local wdb_ro `wdb_ro' safety(`safetyinput')
    if `regtol' != 9 local wdb_ro `wdb_ro' regtol(`regtol')
    if `oostreatedcontrols' != 0 local wdb_ro `wdb_ro' oostreatedcontrols(`oostreatedcontrols')
    if `anticipation' != 0 local wdb_ro `wdb_ro' anticipation(`anticipation')
    tempname wdbobs wdbobsb wdbcomp wdbacc wdbsacc wdbpwsacc wdbjwacc wdbjobs wdbaseobs wdbawacc wdbaobs
    tempname wdbcell
    qui est sto `wdbobs'
    matrix `wdbobsb' = e(b)
    * the physical re-id leg is taken for the cell-weighted estimands (their bootstrap is
    * only correct when a twice-drawn cluster contributes its cells twice) AND, since
    * 2026-09-04, for boottstat under cluster(): each replication's own standard error
    * has to cluster on the RESAMPLED cluster identity, and only a physical expansion can
    * carry one - frequency weights on the original rows cannot say "this cluster is here
    * twice as two clusters", so under the count leg a twice-drawn cluster's copies were
    * pooled into one score and the cluster-robust meat took m^2 where it needed m.
    * Measured inflation E[m^2]/E[m] = 2 - 1/G, ~sqrt(2) in the standard error, which
    * broke percentile-t pivotality. Author ruling 2026-09-04
    local wdb_reid = cond("`bootstrap'" == "standard" & ("`att_it'`att_i'`att_itime'" != "" | ("`boottstat'" != "" & "`cluster'" != "")), "reid", "")
    * UNCLUSTERED + user weight + boottstat: physical expansion WITHOUT re-identification.
    * A frequency weight cannot compose with a user aweight/pweight (their product w*m is
    * neither a count nor a sampling weight), but duplicating the drawn rows physically and
    * letting each copy carry the user's own weight represents both facts. No unit re-id:
    * the resampling unit is the observation, so re-identifying by draw would fragment every
    * unit into single-period units. Duplicate (i,t) rows are tolerated by the observation-
    * weighted estimands (att/ame); the cell-weighted ones are refused without cluster()
    * at parse and never reach this leg. Author ruling 2026-09-04
    local wdb_phys = cond("`bootstrap'" == "standard" & "`boottstat'" != "" & "`cluster'" == "" & "`wdp_wt'" != "" & "`poisson'" == "", "phys", "")
    tempname wdbW wdbjk
    cap noi wooldid_bootrun , wmat(`wdbW') jkacc(`=cond("`bootdisplay'" == "bca", "`wdbjk'", "")') y(`y') i(`i') t(`t') ttrevar(`ttre_asgiven') mode(`bootstrap') `wdb_reid' `wdb_phys' reps(`reps') dots(`dots') rseed(`rseed') bootalpha(`bootalpha') clus(`cluster') repopts(`wdb_ro') wtype(`wdp_wt') wexp(`wdp_wx') ifr(`wdp_if') pois(`poisson') obsb(`wdbobsb') accname(`wdbacc') compname(`wdbcomp') seaccname(`=cond("`boottstat'" != "", "`wdbsacc'", "")') pwseaccname(`=cond("`boottstat'" != "", "`wdbpwsacc'", "")') jwaccname(`=cond("`boottstat'" != "", "`wdbjwacc'", "")') jobsname(`=cond("`boottstat'" != "", "`wdbjobs'", "")') awaccname(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbawacc'", "")') aobsname(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbaobs'", "")') advesrel(`esrelativeto') advtt(`timetrends') advalt(`altmethod_adversarial') cellacc(`=cond("`gteffects'`makecfxplots'" != "", "`wdbcell'", "")')
    local wdb_rc = _rc
    qui est restore `wdbobs'
    cap est drop `wdbobs'
    global WOOLDID_LASTRC
    if `wdb_rc' != 0 exit `wdb_rc'
    * post-loop rewiring marker (see the permutation block): cleared at the end of this block
    global WOOLDID_LASTRC postloop
    * A10: the re-id leg duplicates clusters but user covariates keyed to unit
    * identity keep their original codes across copies (the official idcluster
    * convention shares this) - say so once
    if "`wdb_reid'" != "" {
      local wdb_ikeyed
      foreach wdbcv in `fe' `interactivecontrols' `controls' `coarsecohortcontrols' `ccc_absorb' {
        if "`wdbcv'" == "`i'" local wdb_ikeyed 1
      }
      if "`wdb_ikeyed'" != "" dis as text "Caution: a covariate specification names the i variable `i'; under resampling, duplicated clusters' copies share their original `i' codes inside fe()/controls() terms (the official bootstrap idcluster convention). Only the cell/aggregation identities are re-identified."
    }
    * bootdf() default, resolved here because it needs the resampling-cluster count that
    * wooldid_bootrun has just returned: standard normal when the run is unclustered, and
    * t with G-1 degrees of freedom when it is clustered, G being the number of resampling
    * clusters. the bootstrap standard error is estimated rather than known, and where a
    * cluster count exists it is the natural denominator for that correction; with no
    * cluster() the resampling unit is the observation, so a G-1 rule would be N-1 and
    * indistinguishable from the normal. an explicit bootdf() always wins (author ruling
    * 2026-09-02; the value is validated in wooldid_bootctl, and both branches here are
    * valid by construction)
    if "`bootdf'" == "" local bootdf = cond("`cluster'" != "", "clusters", "normal")
    * rewire the posted inference surfaces to the bootstrap covariance (eclass sub)
    wooldid_bootwire , mode(`bootstrap') reps(`reps') level(`level') bootalpha(`bootalpha') clus(`cluster') acc(`wdbacc') obsb(`wdbobsb') comp(`wdbcomp') `ties' done(`wdb_r_done') fail(`wdb_r_fail') blank(`wdb_r_blank') rseed(`rseed') bootdf(`bootdf') bdisp(`bootdisplay') nclus(`wdb_r_nclus') seacc(`=cond("`boottstat'" != "", "`wdbsacc'", "")') aseout(`=cond("`boottstat'" != "", "`wdbaseobs'", "")') jkacc(`=cond("`bootdisplay'" == "bca", "`wdbjk'", "")') sb(`simultaneousbands')
    * bootstrapped joint tests, subgroup contrasts and CI flavors (eclass sub)
    wooldid_bootjnt , level(`level') acc(`wdbacc') obsb(`wdbobsb') comp(`wdbcomp') jt(`jointtests') `ties' done(`wdb_r_done') bootdf(`bootdf') bts(`boottstat') bdisp(`bootdisplay') nclus(`wdb_r_nclus') clus(`cluster') pwseacc(`=cond("`boottstat'" != "", "`wdbpwsacc'", "")') jwacc(`=cond("`boottstat'" != "", "`wdbjwacc'", "")') jobs(`=cond("`boottstat'" != "", "`wdbjobs'", "")') aseobs(`=cond("`boottstat'" != "", "`wdbaseobs'", "")') seacc(`=cond("`boottstat'" != "", "`wdbsacc'", "")') awacc(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbawacc'", "")') aobs(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbaobs'", "")') adv(`adversarial') advalt(`altmethod_adversarial') advesrel(`esrelativeto') advtt(`timetrends') jkacc(`=cond("`bootdisplay'" == "bca", "`wdbjk'", "")')
    * bootsave(): the replication draws as a .dta in the working directory (author request
    * 2026-09-04). One row per replication (_rep = 1..reps; a failed replication keeps its
    * row, all missing) plus a _rep = 0 row carrying the observed estimates, so the file
    * can recentre its own draws; under boottstat also the replication and observed analytic
    * standard errors (se_*) and, with subgroups, the replication pairwise-contrast se
    * (pwse_*). Read from the accumulators, not e(bootdist), so the 10,000-row cap on the
    * posted matrix does not apply
    * the joint/adversarial accumulators were MOVED into e(boot_jointW)/e(boot_advW) by
    * wooldid_bootjnt (ereturn matrix moves its source), so bootsave reads them from e()
    if "`bootsave'" != "" wooldid_bootsave , file(`bootsave') acc(`wdbacc') obsb(`wdbobsb') seacc(`=cond("`boottstat'" != "", "`wdbsacc'", "")') aseobs(`=cond("`boottstat'" != "", "`wdbaseobs'", "")') pwseacc(`=cond("`boottstat'" != "", "`wdbpwsacc'", "")') wmat(`=cond("`cluster'" != "", "`wdbW'", "")') wobs(`=cond("`bootstrap'" == "bayes", cond("`bootalpha'" == "size", 1/max(`wdb_r_nrows', 1), 1/max(`wdb_r_nclus', 1)), 1)') jwacc(`=cond("`boottstat'" != "", "e(boot_jointW)", "")') jobs(`=cond("`boottstat'" != "", "`wdbjobs'", "")') awacc(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "e(boot_advW)", "")') aobs(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbaobs'", "")') seed(`rseed') rng(`wdb_r_rng0') `replace'
    ereturn local boot_rngstate `wdb_r_rng0'
    * deferred display: worker skipped these under bootstrap; they read only e(),
    * which now carries the bootstrap inference (identical formatting, free replay)
    wooldid_display , part(1)
    * Roth pre-trend power on the BOOTSTRAP covariance. Author ruling 2026-09-04 lifting
    * the plan-A13 refusal: e(V_es_*) has just been replaced with the bootstrap version,
    * so running the diagnostic here feeds Roth's own formula a bootstrap covariance
    * rather than an analytic one. The two hbar scalars were computed at worker level
    * (they need tempvars that no longer exist) and are handed back in. Placed between
    * part(1) and part(2) to match the ordering the non-bootstrap path uses
    if "`pretrendpower'" != "" {
      cap noi wooldid_ptp , levels(`pretrendpower') esrel(`esrelativeto') dfr(`=cond(e(df_r) < ., string(e(df_r)), ".")') pois(`poisson') poisexp(`poisexpresults') hba(`=cond(e(wdid_ptphbar_att) < ., string(e(wdid_ptphbar_att)), ".")') hbi(`=cond(e(wdid_ptphbar_it) < ., string(e(wdid_ptphbar_it)), ".")') level(`=cond(`level' == 0, c(level), `level')')
      local wdsh_prc = _rc
      if `wdsh_prc' != 0 {
        qui frame
        if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
        global WOOLDID_LASTRC `wdsh_prc'
        exit `wdsh_prc'
      }
    }
    * cell-level surfaces rewired to the bootstrap (ruling 2026-09-04): e(gteffects) columns
    * se/p/lb/ub and each e(cfx*) matrix's se/t/pvalue/ll/ul, per bootdisplay(); studentized
    * under boottstat from the replications' own cell standard errors
    if "`gteffects'`makecfxplots'" != "" {
      * the cfx accumulators are numbered by wooldid_bootrun's existence count and named
      * in its c_locals; pass the names explicitly (as the permute block does) rather than
      * re-deriving them from e(wdid_cfxmats) by position (F046, 2026-09-06)
      local wdb_cfxnames
      forval wdk = 1/`wdb_ncfx' {
        local wdb_cfxnames `wdb_cfxnames' `wdb_cfxname`wdk''
      }
      wooldid_bootcellsx , cell(`wdbcell') ncfx(`wdb_ncfx') names(`wdb_cfxnames') bdisp(`bootdisplay') bts(`boottstat') level(`level') bootdf(`bootdf') nclus(`wdb_r_nclus') `ties'
    }
    * continuous-treatment plots, deferred like the event-study plots (author ruling
    * 2026-09-05): drawn from the rewired e(cfx*) matrices, so their bands are the
    * displayed bootstrap family's; the worker drew nothing and saved nothing under bootstrap()
    if "`makecfxplots'" != "" {
      local wdcfopts
      if "`saveplots'" != "" local wdcfopts saveplots(`saveplots')
      cap noi wooldid_cfxplots , `wdcfopts' `replace'
      local wdsh_prc = _rc
      if `wdsh_prc' != 0 {
        qui frame
        if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
        global WOOLDID_LASTRC `wdsh_prc'
        exit `wdsh_prc'
      }
    }
    wooldid_display , part(2)
    * event-study plots, deferred from the worker so they read the REWIRED tables. The
    * bands are therefore whichever family bootdisplay() selected - percentile, BC or
    * wald - and match the printed table exactly. A plot failure is an error of this call
    * CS-87 (consistency sweep 2026-09-06): the worker draws the event-study plot only when
    * horizons exist (esprelength/espostlength > 0 after the panel clips them); the deferred
    * engine call must apply the same guard, else wooldid_esplots runs on an empty frame,
    * dies r(2000) and loses the whole run. The horizons are visible here as ES_ columns of
    * the (rewired) e(b)
    if "`makeplots'" == "makeplots" & strpos("`: colfullnames e(b)'", "ES_") > 0 {
      local wdesopts
      if "`saveplots'" != "" local wdesopts saveplots(`saveplots')
      cap noi wooldid_esplots , `wdesopts' `replace'
      local wdsh_prc = _rc
      if `wdsh_prc' != 0 {
        qui frame
        if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
        global WOOLDID_LASTRC `wdsh_prc'
        exit `wdsh_prc'
      }
    }
    wooldid_display , part(boot)
    * the list shrank on 2026-09-05: gteffects cells and the cfx matrices are now rewired to
    * the bootstrap (wooldid_bootcellsx), the adversarial p is bootstrapped under boottstat and
    * absent otherwise, and event-study plots are drawn after the rewiring
    if "`histogramcohorteffects'" != "" dis as text "Caution: histogram outputs are drawn from the ANALYTIC covariance during estimation; the reported tables, the gteffects cells, the cfx matrices, the event-study and continuous-treatment plots and the e() inference above are bootstrap-based."
    dis as result " "
    global WOOLDID_LASTRC
  }
  * P7 accident 2 (author ruling 2026-09-08): the analytic twin of the engine branches' own
  * wooldid_cfxplots call. Runs only when the worker recorded a stem, which it does only on an
  * analytic makecfxplots run with saveplots() - so neither engine reaches this, and a run that
  * asked for no plot data is untouched. cap, but NOT silent: a failure to write the data file
  * must not lose an estimation that has already succeeded and printed - which is not true of
  * the engine branches, where the same call also DRAWS the figures - and equally must not be
  * swallowed. P6's own published policy is that no accepted option is a silent no-op, and a
  * saveplots() that quietly produced the .gph files and not the .dta would be exactly that.
  if "`e(wdid_cfxsave)'" != "" {
    cap noi wooldid_cfxplots , saveplots(`e(wdid_cfxsave)') dataonly `e(wdid_cfxrepl)'
    local wdsh_dc = _rc
    if `wdsh_dc' != 0 dis as text "Caution: the continuous-treatment plot-data file `e(wdid_cfxsave)'_cfx.dta could not be written (r(`wdsh_dc')). The figures and their .gph files are unaffected, as is every estimate and every stored result; only that one machine-readable record of the plotted points is missing. The usual cause is an existing file of that name without the replace option."
  }
end

cap prog drop wooldid_worker
prog def wooldid_worker, eclass
version 16.0
  * internal replication marker, set by wooldid_permrun/wooldid_bootrun around each
  * nested call (inside nobreak, so it is always cleared) and consumed READ-ONCE here.
  * Replication results are read for their coefficients alone, so blocks that only fill
  * diagnostic e() returns can be skipped there. Read-once means even a leaked marker -
  * which the setters' nobreak discipline already prevents - affects at most one run.
  local wdrepmode "${WOOLDID_REPMODE}"
  global WOOLDID_REPMODE
  * the stored full model. Replications (nested wooldid calls under bootstrap() or
  * permute()) store theirs under a separate name and always drop it at exit, so they
  * never overwrite - or, at their exit, delete - the observed run's ___wooldidfullmodel
  * that verbose retains for the user (until 2026-09-05 a verbose bootstrap run left the
  * LAST replication's model under that name)
  local wdfm ___wooldidfullmodel
  if "`wdrepmode'" != "" local wdfm ___wooldidrepmodel
  * Every Stata replication prefix (bootstrap/jackknife/simulate/permute) runs its
  * replication loop under `set coeftabresults off', which silently removes margins
  * r(table) - and this worker harvests r(table) at 16+ sites, dying r(503)
  * "conformability error" on every replicate. The setting is program-scoped, so this
  * line is automatically undone at worker exit and is a no-op in a default session.
  set coeftabresults on

* S5: this line must remain BYTE-IDENTICAL to the syntax line in wooldid_validate
syntax [varlist(numeric default=none)] [if/]  [aw fw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist numeric) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose fastcollapse fullprecisionslow GTEffects PERMute(string) PERMREPS(integer 0) PERMEXact PERMTStat PERMTci PERMSTRata(varname numeric) PERMPRobs(string asis) PERMSAve(string) PERMANalyticse PERMABsp PERMSide(string) PERMRAWp PERMSTEPdown SIMULTaneousbands BOOTstrap(string) reps(integer 0) ties Level(real 0) dots(integer -1) rseed(string) BOOTALPha(string) BOOTDf(string) BOOTDISPlay(string) BOOTTStat BOOTSAve(string)]

* Execute update before running main program
  if "`update'" == "update" {
    net install wooldid, from(https://raw.githubusercontent.com/thegland/wooldid/main/wooldidinstall/) replace
    dis as text "Update complete; please restart stata or use 'discard' command to load the updated version of wooldid. Other dependencies (reghdfe, ftools, gtools, ppmlhdfe (optional)) must be updated manually if necessary."
    exit
  }

  * S5: all parse-time validation and option normalization lives in wooldid_validate,
  * which re-parses `0' with a BYTE-IDENTICAL copy of the syntax line above (guarded by
  * test/check_validate_sync.py) and returns its option/derived locals via c_local.
  * called before cleanmatrices, so a rejected call no longer clears user matrices/mata
  wooldid_validate `0'
  * F054 (quality review 2026-09-05, author ruling 2026-09-06): level() applies to the
  * analytic tables too. Every analytic interval comes from margins' r(table) rows at
  * c(level), so the level is set for the duration of this call and the shell restores
  * the caller's value on every exit path (WOOLDID_LEVEL0); replications inherit the
  * outer call's setting and never touch it
  if `level' != 0 & "`wdrepmode'" == "" {
    global WOOLDID_LEVEL0 = c(level)
    set level `level'
  }
  local wdlvl = c(level)
  local wdlvltxt = cond(mod(`wdlvl',1)==0, strofreal(`wdlvl'), subinstr(strofreal(`wdlvl',"%9.0g"), ".", "p", .))


  ***************************************************************************************************
  ******************************  Initialize and Process Command Input ******************************
  ***************************************************************************************************
  * This section converts user input into actionable information and checks for validity of the user input.



  * Clean the workspace if cleanmatrices is specified - this is helpful b/c sometimes reghdfe enters into a failure state
  * that leaves junk matrices floating around in stata and mata which prevent further estimation. cleanmatrices should fix that.
  if "`cleanmatrices'" != "" {
    clear mata
    matrix drop _all
  }

 * Unless disabled by the user, set emptycells to drop -- greatly improves performance since structure of the reghdfe matrix
 * involves lots of empty cells in the interaction term.
 if "`emptycellsoverride'" == "" {
    * note: set emptycells is program-scoped - stata automatically restores the caller's
    * value when wooldid exits, including on error paths (cf. StataCorp's anova.ado);
    * no manual capture/restore is needed here (B18 was verified illusory, 2026-07-26)
    set emptycells drop
 }


  * extract the regression weights variable and setup macros for handling weights; if no weight variable, create one. require weights to be non-negative everywhere
  local w `exp'

  * F6: hand the shell what it needs to orchestrate the permutation rep loop (the
  * loop lives in the shell for compiled-size reasons; the c_local blanket excludes
  * the if/weight pieces, so they are echoed here)
  if "`permute'`bootstrap'" != "" {
    if "`permute'" != "" c_local wdp_on 1
    if "`bootstrap'" != "" c_local wdb_on 1
    c_local wdp_if `"`if'"'
    c_local wdp_wt `weight'
    c_local wdp_wx `w'
  }
  local weightinserter [`weight' = `w']
  * fw is parseable only so the bootstrap replication path can pass an fweight; a user
  * typing it directly is refused. The demotion back to aweight on the next line is
  * LOAD-BEARING and measured: `weightinserter' was built one line above, so the fweight
  * reaches the REGRESSION and nothing else, while margins/gstats/gquantiles keep today's
  * aweight semantics. Without the demotion, customweightmultiplier() makes the margins
  * weight w*cwm non-integer even when w itself is integral and every replication dies
  * with r(459) - cwm is deliberately NOT folded into `w' (see the `local w `exp'' above),
  * which is exactly why splitting them works
  if "`weight'" == "fweight" & "`wdrepmode'" == "" {
    dis as error "Error: wooldid does not accept fweights. Specify aweights or pweights."
    error 198
  }
  if "`weight'" == "fweight" local weight aweight
  if "`w'" == "" {
    tempvar w
    qui gen `w' = 1
    local weightinserter
    local weight iw
    * when weights are absent, we need to tell what kind of weights margins should use
    * iw in margins is treated like fw w/o the N adjustment
    * that's desirable for our purposes
  }





  ***************************************************************************************************************
  ******************************  PRODUCE IN-SAMPLE INDICATOR / HANDLE MISSINGS  ********************************
  ***************************************************************************************************************


  * create a variable signifying sample membership and update it to reflect restrictions imposed by the user in "`if'"
  * as well as the following restrictions: non missing i & t identifiers, non missing outcomes, non-missing weights, weights are positive and exceed zero, non missing controls, non missing cluster variables,
  * note: the loop below extends the non-missingness restriction to the variables inside fe(), controls(), and the other passed varlists (factor notation stripped)
	tempvar ifvar
  if "`if'" != ""{
		qui gen `ifvar' = `if'
  }
  else{
		qui gen `ifvar' = 1
  }

	* BS-37 (bug sweep 2026-09-06, author ruling 2026-09-07 "keep equals-one, fix the
	* reporting"): the exact-one rule STAYS. It is what keeps a bare flag variable's
	* missing and out-of-range values out of the sample, and all six forwarded copies of
	* the expression share it (dropcount below, wooldid_validate's ifassert, permctl's
	* insmp, permrun's and bootrun's smp), so the estimation sample, the permutation
	* design and the bootstrap frame agree. What was wrong was the silence: rows the rule
	* alone removed drew no message and were counted nowhere - "if flag" with flag in
	* {1,2} lost the flag == 2 rows under a header reading "0 obs dropped". Count off
	* `ifvar' rather than off (`if') so the count is exactly what the rule did, and count
	* BEFORE the markout replace below, which would zero these rows for an unrelated
	* reason. wdvifexcl is initialised here so the fold at the dropped-observation site is
	* well-formed on every path. Body in a subprogram: the worker is at the size ceiling
	local wdvifexcl 0
	if "`if'" != "" wooldid_ifrule , ifv(`ifvar')

	* save this for counting to see if we drop observations
	local dropcount ((`if') == 1)
	if "`if'" == "" local dropcount (1 == 1)

  * i, t, y, and the weight variable must all be non missing. weights must exceed 0 to be included. a loop below imposes similar requirements for control variables and the continuous treatment variable
  qui replace `ifvar' =  0 if `i' >= . | `t' >= . | `y' >= . | `w' >= .  | `w' <= 0

  * sanitize controls passed that may contain factor variable notation
  local loopablecontrols
  foreach factorvarspossible in fe controls coarsecohortcontrols      {
    local processtemp ``factorvarspossible''
    local processtemp = subinstr("`processtemp'","#"," ",.)
    local processtemp = subinstr("`processtemp'","i."," ",.)
    local processtemp = subinstr("`processtemp'","c."," ",.)
    local loopablecontrols `loopablecontrols' `processtemp'
  }

  local checkclus
  if "`cluster'" != "`i'" local checkclus `cluster'

  local checklattice
  if wordcount("`lattice'") == 1   local checklattice `lattice'

  * normally, it is ok to have missing values in the sg variable -- we only produce estimates for flagged sg groups
  * not so if using sg to do a triple diff - here, missing values actually are an issue
  * that said, we aren't really implementing that at this time
  foreach mustbenonmissingvar in `checkclus'  `loopablecontrols' `ccc_absorb' `customweightmultiplier' `contreat' `checklattice'   `interactivecontrols' {
    qui replace `ifvar' = 0 if `mustbenonmissingvar' >= .
  }


  local ifstatement `ifvar' == 1

  * NFCOMP-02: the customweightmultiplier() domain is documented as >= 0 but was never
  * enforced - a negative value silently corrupts weight masses, estimates, and the
  * escomp shares under the default iweight plumbing (probed rc 0, no message). checked
  * on the estimation-aware domain (ifstatement folds the markout, so negative values
  * confined to deterministically-discarded rows stay tolerated - the G15A precedent);
  * missing values pass (already marked out above). safetyinput-gated per the F13/S5
  * gated-screen convention for data-integrity checks
  if "`customweightmultiplier'" != "" & "`safetyinput'" != "off" {
    cap assert `customweightmultiplier' >= 0 if `ifstatement'
    if _rc != 0 {
      dis as error "Error: the customweightmultiplier() variable contains negative values on in-sample observations; the documented domain is non-negative (missing values are dropped). To ignore this error, specify safety(off)."
      error 411
    }
  }


  *********************************************************************************************************************************
  ******************************  PROCESS EVENT STUDY INFO, ESTIMABILITY CHECKS, FORM TREATMENT DUMMIES ***************************
  *********************************************************************************************************************************
  * In this section, we form treated_cohort X treated_period dummies
  * and also figure out how large the panel is and how many observations we have across specific time periods - imposing
  * estimability constraints where appropriate.
  * We do this differently depending on if we are using a fixed or an pooled reference period.

  gstats sum `t' if `ifstatement', meanonly
  local maxyearobserved = r(max)
  local minyearobserved = r(min)
  tempvar ttre

  * create a modified version of the period-of-treatment-onset variable ttre
  * the modified version will treat cohorts with treatment onsets a set distance (specified in oostreatedcontrols)
  * from the end of the sample identically as to units that are never treated at all
  * F4 anticipation(): treat the effective onset as ttre - #, shifted ONCE here so every
  * downstream consumer (window screens, _t construction, summarystats, margins subpops)
  * sees the shifted onset. the shift applies to BOTH the value and the out-of-sample
  * window screen (a value-only shift would leave out-of-window anticipating cohorts in
  * the control pool). the single-expression form is bit-identical to the hand-edit
  * workaround (gen ttre2 = ttre - #); with anticipation absent/0 the emitted command is
  * unchanged from before the option existed
  local anticshift
  if `anticipation' > 0 local anticshift - `anticipation'
  qui gen double `ttre' = `ttre_asgiven' `anticshift' if `ttre_asgiven' `anticshift' <= `maxyearobserved' + `oostreatedcontrols' & `ifstatement'

  * F13: the _i/_t construction below assumes cohort-level timing - exactly one treatment onset
  * per unit. rows of the same i carrying different onsets, or a mix of missing and non-missing
  * onsets, are silently re-timed row by row: part of the unit joins the control pool or a
  * second cohort, shifting event-study horizons and pre-trends with no error. detect on the
  * internal `ttre' (the window screen above normalized extended missings to sysmiss) over the
  * screened estimation sample - exactly the rows the construction reads. all missing flavors
  * fold to one sentinel value so all-missing (never-treated) units pass while mixed
  * missing/nonmissing within a unit fails
  if "`safetyinput'" != "off" {
    tempvar f13ttremin f13ttremax
    qui gegen double `f13ttremin' = min(cond(`ttre' >= ., -8.9884656743115785e+307, `ttre')) if `ifstatement', by(`i')
    qui gegen double `f13ttremax' = max(cond(`ttre' >= ., -8.9884656743115785e+307, `ttre')) if `ifstatement', by(`i')
    cap assert `f13ttremin' == `f13ttremax' if `ifstatement'
    if _rc != 0 {
      tempvar f13badi
      qui gen double `f13badi' = `i' if `f13ttremin' != `f13ttremax' & `ifstatement'
      qui gstats sum `f13badi', meanonly
      local f13badid = strtrim(string(r(min), "%18.0g"))
      dis as error "Error: the treatment onset variable `ttre_asgiven' is not constant within `i' (smallest offending `i': `f13badid')."
      dis as error "wooldid requires exactly one treatment onset per group defined by `i': rows of the same `i' carrying different onsets, or a mix of missing and non-missing onsets, are silently split across cohorts, and part of the group can be used as its own control."
      dis as error "Common causes: an onset recorded only from the onset period onward (fill it in on every in-sample row of a treated group, including pre-treatment rows); underlying units with different onsets combined under one `i' value; or IDs stored as float, which rounds distinct large IDs together (regenerate `i' as double from its source; recasting cannot recover the lost values). Alternatively, restrict the estimation sample with if()."
      dis as error "To ignore this error, specify safety(off)."
      error 459
    }
    qui drop `f13ttremin' `f13ttremax'
  }

  * prepare the _i and _t variables, which are cohort and time variables where non-treated-unit non-treated-times are all pooled
  * noting that "non-treated-times" should be extended here to actually mean "any time period for treated units in its reference period"
  tempvar _i
  tempvar _t

  * we begin by handling the case where the user wants the comparison period to be the beginning of the sample, possibly multiple such periods
  if "`esfixedbaseperiod'" != "esfixedbaseperiod" {
    * relative-time construction lives in wooldid_reltime (worker compiled-size headroom,
    * 2026-09-07): identical behaviour; _t, _i, ifvar and ttre are mutated in place under
    * the caller's tempvar names and the four window locals come back by c_local
    wooldid_reltime , ivar(`i') tvar(`t') ttrevar(`ttre') rawttre(`ttre_asgiven') ifv(`ifvar') relt(`_t') relid(`_i') epre(`esprelength') epost(`espostlength') antic(`anticipation') ifs(`ifstatement') tt(`timetrends') sfty(`safety') ct(`contreat') cw(`contreatwithin')
  }
  else{
    * enter the es fixed base period case - esfixedbaseperiod
    * by default, estimates are relative to the first time period prior to treatment (-1);
    * values above -1 are rejected in wooldid_validate (formerly they were silently clamped
    * to -1 here, running with a different reference period than requested)
    * we also modify esprelength to stretch to the relative period
    if `esrelativeto' < -1 & `esprelength' > 0 & `esprelength' < abs(`esrelativeto') {
      local esprelength = abs(`esrelativeto') -1
    }
    if `esprelength' >= abs(`esrelativeto') & "`timetrends'" == "timetrends" local ++esprelength

    local refperiod = `esrelativeto'

    qui gen `_t' = `t' - `ttre' if `ifstatement'


    * VAL-G10B: the calendar-support screen runs FIRST, before the reference checks.
    * support-first is a fixed point: support is furnished only by never-treated rows,
    * which the reference screens never drop, so running support first cannot create new
    * unsupported periods. the old order could delete a group's only reference row AFTER
    * that group had passed the reference check, sending a structurally broken sample
    * downstream to die opaquely at the dropped-dummy screen (r322) - or run silently
    * wrong under safety(off) on some shapes
    * produce a check that each relative time effect is estimable in that there is an untreated cohort for it
    * controlcheck2 will be 1 if there is at least one in sample non treated unit in a given period
    * this logic works b/c we are in esfixedbaseperiod land and we really do care about each period - we aren't pooling anything
    * as above, don't do this if only using within-treated cohort continuous treatment variation
    if "`contreatwithin'" == "" {
      tempvar controlcheck
      tempvar controlcheck2
      qui gen `controlcheck' = `ttre' ==. & `ifstatement'
      qui gegen `controlcheck2' = max(`controlcheck') if `ifstatement' , by(`t')
      qui replace `ifvar'  = 0 if `controlcheck2' == 0
      if "`safety'" != "off"  {
        cap assert `controlcheck2' != 0
        if _rc != 0 {
          dis as text "Caution: estimability checks found time periods exist in sample without a control unit, causing them to be dropped."
          if "`contreat'" != "" dis as text "If intending to only use within treated-cohorts variation in your continuous treatment variable, specify option contreatwithin to avoid these checks."
        }
      }
      qui drop `controlcheck' `controlcheck2'
    }

    * update _t and ttre
    qui replace `_t' = . if `ifvar' == 0
    qui replace `ttre'  = . if `ifvar' == 0


    * here we implement a check to see whether or not the reference period exists for a given cohort
    * if not, we do not use the cohort (support-dropped periods no longer count as observed
    * reference rows: ttre was set missing above on support-dropped rows)
    tempvar refperiodobserved_temp
    tempvar refperiodobserved
    qui gen `refperiodobserved_temp' = `_t' == `refperiod' if `ifstatement' & `ttre' != .
    qui gegen `refperiodobserved' = max(`refperiodobserved_temp') if `ifstatement' & `ttre' != ., by(`i')
    qui drop `refperiodobserved_temp'
    qui replace `ifvar' = 0 if `refperiodobserved' == 0  & `ttre' != .
    * F4: attribute reference-period drops caused only by the anticipation shift
    if `anticipation' > 0 cap noi wooldid_anticesf , ivar(`i') tvar(`t') rawttre(`ttre_asgiven') refobs(`refperiodobserved') esrel(`esrelativeto') antic(`anticipation') ttrevar(`ttre') ifv(`ifvar')
    * update _t and ttre
    qui replace `_t' = . if `ifvar' == 0
    qui replace `ttre'  = . if `ifvar' == 0


    * if using time trends, you need to see a reference period block that is 2 periods in length; drop units where this is impossible
    if "`timetrends'" == "timetrends" {
      tempvar ttrefperiodobserved_temp
      tempvar ttrefperiodobserved
      qui gen `ttrefperiodobserved_temp' = `_t' == `refperiod' -1 if `ifstatement' & `ttre' != .
      qui gegen `ttrefperiodobserved' = max(`ttrefperiodobserved_temp') if `ifstatement' & `ttre' != ., by(`i')
      qui drop `ttrefperiodobserved_temp'
      qui replace `ifvar' = 0 if `ttrefperiodobserved' == 0  & `ttre' != .
    }
    * update _t and ttre
    qui replace `_t' = . if `ifvar' == 0
    qui replace `ttre'  = . if `ifvar' == 0

    * grab information about range of _t
    qui gstats sum `_t' if `ifstatement', meanonly
    * abs() here and `max(0, -r(min))' in wooldid_reltime's pooled twin are deliberately
    * different, and this form is the safe one only under esfixedbaseperiod: the
    * reference-period screen above zeroes ifvar for every unit that does not OBSERVE
    * _t == esrelativeto, and esrelativeto is rejected above -1, so every treated unit
    * reaching this line has a realized _t <= -1 and r(min) <= -1. The two expressions
    * therefore agree here, and the BS-33 failure mode the twin guards against (an
    * onset-period gap or an off-grid onset leaving r(min) POSITIVE, so abs() invents a
    * pre-window nothing realizes) is unreachable. Do not "fix" this one, and do not read
    * the difference as the BS-33 fix having been applied to only one of two sites
    local min_t = abs(r(min))
    local max_t = r(max) + 1


    * set _t into all positive terms to allow usage in regressions
    qui replace `_t' = `_t'+ `min_t' + 1

    *redefine the ref period macro for future use and then do a little bit of housekeeping
    * LOAD-BEARING: wooldid_advstat plants the adversarial reference row(s) at
    * t == esrelativeto + min_t + 1 via this same identity (relyr = t - min_t - 1);
    * if this definition moves, that subprogram must move with it
    local refperiod = `esrelativeto'+`min_t'+1
    qui drop `refperiodobserved'
    if "`timetrends'" == "timetrends" qui drop `ttrefperiodobserved'


    * check that there is a real sample in that someone, somewhere, has a pre-period and a post period
    if "`safety'" != "off"   {
      qui  gstats sum `_t'  if `ifstatement' & `ttre' != ., meanonly
      local failure 0
      cap assert r(N) != 0
      if _rc != 0 local failure 1
      cap assert r(min) <= `refperiod'
      if _rc != 0 local failure 1
      cap assert r(max) > `refperiod'
      if _rc != 0 local failure 1
      if `failure' == 1 dis as error  "Error: cannot find sufficient pre- and post-treatment periods for the treated units."
      if `failure' == 1 & "`timetrends'" == "timetrends"  dis as error  "Note that inclusion of cohort-specific time trends increases the required number of pre-treatment periods."
      if `failure' == 1 & `anticipation' > 0 dis as error  "Note: anticipation(`anticipation') shifted every effective onset earlier; a shift larger than the available pre-period window can leave no estimable treated periods."
      if `failure' == 1 error 459
    }
    * send _t to 0  when we want it to be  a comparison period
    qui replace `_t' = 0 if (`_t' == . | `_t' == `refperiod' ) & `ifstatement'
    if "`timetrends'" == "timetrends" qui replace `_t' = 0 if `_t' == `refperiod' - 1 & `ifstatement'

    * produce modified identifiers
    * _i is unit identifiers that are positive integers if i is ever treated in sample and 0 otherwise
    qui gen double `_i' = `i' if  `_t' > 0 & `_t' < . & `ifstatement'
    qui replace `_i' = 0 if `_i' == . & `ifstatement'


    **** gather some event study info
    local espre = min(`min_t', `esprelength')
    local espl = min(`max_t', `espostlength')
    * end data setup work in the case using  a fixed base period

    * when making horizons, non-fixed base will always present all horizons flagged by _t, since it only generates flags if they will be used
    * not so when picking a fixed base -- _t will take on many values the user may or may not be interested in reporting
    * in that case, horizon = 1 might not be of any interest
    * so what you want to start at is   1+`min_t' - `espre'
  }

  * F4: warn when an effective onset lands on a period not observed in the data (both
  * reference regimes; the checks live in a subprogram to conserve worker size headroom)
  if `anticipation' > 0 cap noi wooldid_anticgrid , tvar(`t') ttrevar(`ttre') ifv(`ifvar')

  * this is needed for correctly rendering tables / looping through the right number of horizons
  local horizonstarter = 1 + `min_t' - `espre'
  if "`esfixedbaseperiod'" == "" local horizonstarter 1


  ******* given the above processing and estimability checks, do we still have a control group?
  * note that the checks above should already trigger if you don't have a treatment group
  * do not check for contreatwithin b/c you don't expect a control group in that case
  if "`safety'" != "off"  & "`contreatwithin'" == "" {
    qui count if `ttre' == . & `ifstatement'
    if r(N) == 0 {
      dis as error  "Error: you have no never-treated control group and no not-yet-treated control group whose treatment falls outside the sample window."
      dis as error  "To ignore this error, specify safety(off) or modify the sample window manually."
      if `anticipation' > 0 dis as error  "Note: anticipation(`anticipation') moves effective onsets earlier, which can pull cohorts treated after the sample window into it; a larger oostreatedcontrols() or a smaller anticipation() may restore a control group."
      error 459
    }
  }

  * if number of requested horizons is too short, don't do adversarial pre-trend testing
  * the + 1 is to account for the 0 reference period.
  * CS-38 (author ruling 2026-09-08): say so. pretrendpower() refuses at parse time when its
  * prerequisites are missing, while this cleared its option silently and the run completed
  * with no sign the requested test had been dropped. A window the user asked for but the
  * DATA cannot support is a caution here; a window that could never have supported it is
  * refused in wooldid_validate before estimation starts.
  if (`espre' + `espl'+ 1 < 3 & "`timetrends'" != "") | `espl' == 0 | `espre' == 0 {
    if "`adversarial'" != "" dis as text "Caution: adversarial pre-trend testing was dropped: the estimation sample realizes `espre' pre-treatment and `espl' post-treatment horizon(s), and the test needs at least one of each (at least three coefficients in total under timetrends). Widen esprelength()/espostlength() to a window the data can support, or drop adversarial."
    local adversarial
  }

  * B16 guard: under esfixedbaseperiod + altmethod_adversarial the reference line is fit on
  * the pre-period coefficients plus the planted reference row; with esrelativeto(-1) and a
  * single pre horizon that fit sample is one row, so no line is estimable
  if "`adversarial'" != "" & "`altmethod_adversarial'" != "" & "`esfixedbaseperiod'" != "" & `esrelativeto' == -1 & `espre' <= 1 {
    dis as text "Caution: adversarial pre-trend testing was cancelled: with esfixedbaseperiod, esrelativeto(-1), and a single pre-treatment horizon, only the reference row is available to fit the adversarial reference line."
    local adversarial
    local altmethod_adversarial
  }

  * NFADV-01/02: the adversarial p-value is esf-scoped (like pretrendpower); under the
  * pooled reference the fitted line is still computed and returned, without inference
  if "`adversarial'" != "" & "`esfixedbaseperiod'" == "" dis as text "Caution: the adversarial pre-trend p-value requires esfixedbaseperiod (the pooled reference has no single omitted period); the fitted line is reported without a p-value."

  ****** after the above estimability checks and processing, does the SG variable have surviving members?
  * while we're at it, figure out how many levels the sg variable takes on during the treatment period - will be useful down the line

  * subgroup-domain validation lives in wooldid_sgdom (worker compiled-size headroom,
  * 2026-09-07): identical behaviour, gate kept here, sglist returned by c_local
  if "`subgroups'" != "" wooldid_sgdom , sg(`subgroups') ivar(`i') tvar(`t') ttrevar(`ttre') ifs(`ifstatement') si(`safetyinput') `makecfxplotsbysg'






  *********************************************************************************************************************************
  ***********************************  PRODUCE LATTICE FOR CONTREAT, CONTINUOUS TREATMENT VARIABLES  ******************************
  *********************************************************************************************************************************

  * first, create targetflag, a variable that is 1 in relative periods for which we are interested in producing treatment effects, and 0 otherwise
  * this ends up being needed in various places as a matter of convenience in building the regressions as well as elsewhere in making these variables
  tempvar targetflag
  qui gen `targetflag' = `_t' > 0 & `_t' < . if `ifstatement'



  * prepare the continuous treatment materials  - this must be done after ttre is fully updated and bad treatments are missinged out
  * lattice is a thing that allows piecewise estimation of the effect of your continuous treatment
  * you can pass through any lattice variable you like, or autogenerate one
  * if you autogenerate, the key concepts here are:
  * latticenum tells us how many chunks to split the sample into for the piecewise continuous treatment function
  * latticetype -- specifies how to compute quantiles. all implies over the whole sample, otherwise they are calculated within subsets of the sample specified by the by variable.
  *     e.g., if baseX is specified, quantiles are formed in a base period and the boundaries in the base period are used to define the lattice everywhere
  * the code block below handles autogenerating: it is largely a matter of generating percentiles within the correct group of variables
  * note that piecewise functions threaten the concept of derivatives, so don't expect clean interpretation of AMEs when using a lattice

  local latticevar
  local latticepart2hash
  local latticepart1hash
  local contreatinj

  * note that this thing will first work on making lattices
  * then work on making the continuous treatment variable text will add it to the regression
  if "`contreat'" != "" {
    * work on the lattice
    if "`lattice'" != ""  {
      * wordcount = 2 implies generate a lattice; wordcount of 1 means you just had the name of the lattice var directly passed through
      if wordcount("`lattice'") == 2 {
        * ignore weights in lattice percentile computation - useful in some situations where weights are not very uniform
        local latticeweightinserter `weightinserter'
        if "`latticeignoreweights'" == "latticeignoreweights" local latticeweightinserter

        * begin specifying commands to construct the different types of lattice variables / create ancillary supporting variables
        tempvar latticevar
        local byable
        local xtileif
        if "`latticetype'" == "baseTreated" local xtileif & `ttre' != .
        if "`latticetype'" == "baseTreatedPost" local xtileif & `t' >= `ttre' & `ttre' != .
        if "`latticetype'" == "baseTreatedPre" local xtileif & `t' < `ttre' & `ttre' != .
        if "`latticetype'" == "baseTreatedRef" local xtileif & `t' < `ttre' & `ttre' != . & `_t' == 0
        if "`latticetype'" == "byCohort" local byable `i'
        if "`latticetype'" == "byTime" local byable `t'
        if "`latticetype'" == "byCohortByTime" local byable `i' `t'
        if  "`latticetype'" == "byTreated" | "`latticetype'" == "byTreatedByTime" {
          tempvar evertreated
          qui gen `evertreated' = `ttre' != .
        }
        if   strpos("`latticetype'", "byTtre")!=0 {
          tempvar modified_ttre
          qui gen double `modified_ttre' = `ttre'
          qui replace `modified_ttre' = `minyearobserved' if `ttre' == .
          * minyearobserved is just any value not appearing in ttre; ttre is guaranteed to not include that value
        }
        if "`latticetype'" == "byTreated" local byable `evertreated'
        if "`latticetype'" == "byTreatedByTime" local byable  `evertreated' `t'
        if "`latticetype'" == "byTtre" local byable `modified_ttre'
        if "`latticetype'" == "byTtreByTime" local byable  `modified_ttre' `t'

        * make the 'base' style lattices - these compute percentiles in the specified base period (treated units only, treated units post treatment, treated units pre treatment)
        * base type needs to be handled special b/c we code this variable with a loop through the grabbed list of percentiles
        * for the others, we just use a percentile generator to create them under the specified conditions
        if "`latticetype'" == "baseTreatedPost" | "`latticetype'" == "baseTreated" | "`latticetype'" == "baseTreatedPre" | "`latticetype'" == "baseTreatedRef" {
          gquantiles `contreat' `latticeweightinserter' if `ifstatement' `xtileif', _pctile   nquantiles(`latticenum')
          gen `latticevar' = 1
          forval latticecase = 2(1)`latticenum'{
            local percentiletograb = `latticecase'-1
            qui replace `latticevar' = `latticevar' + 1 if `contreat' > r(r`percentiletograb')
          }
        }
        else {
          if "`latticetype'" != "all" qui gquantiles `latticevar' = `contreat' `latticeweightinserter' if `ifstatement' , xtile nquantiles(`latticenum') by(`byable')
          if "`latticetype'" == "all"  qui gquantiles `latticevar' = `contreat' `latticeweightinserter' if `ifstatement' ,   nquantiles(`latticenum') xtile
        }
        qui replace `latticevar' = 0 if `latticevar' == . & `ifstatement'
        * null lattice points are sent to 0
      }
      else{
        local latticevar `lattice'
      }
      * build some macros for inserting the lattice variable where needed
      local latticepart2hash ##i.`latticevar'
      local latticepart1hash #i.`latticevar'
    }
    else{
      local latticepart2hash
      local latticepart1hash
    }

    * lattice is done, work on making the continuous treatment variable work
    * finish up making contreatadder - the code chunk that injects into the regression statement the correct variables needed for estimating continuous treatment effects
    * an important part of this is partialling out the base effect of the continuous treatment variable in a fashion that is invisible to margins
    * that is to say, when we control for the continuous treatment variable outside of the treated/pre-treated observations, we don't want that to contribute to treatment effect estimation
    * we can do that by creating a ringer variable that is equal to contreat but with a different name, or by pushing stuff into the absorb() command
    * we handle requests for contreat polynomials using c.#c. interactions so that the polynomial terms are legible to margins; we don't care about doing so for the ringer variable though
    * the ringer copies must be double: gen defaults to float, and a float-truncated copy of a
    * double contreat is only approximately collinear with it, weakening the engineered
    * collinearity drop that keeps the base effect out of the margins-visible terms
    tempvar ringer
    qui gen double `ringer' = `contreat'
    local contreatpolypart c.`contreat'
    local ringerpoly c.(`ringer'
      if `contreatpoly' > 1 {
        forval polycount = 2(1)`contreatpoly' {
          local contreatpolypart `contreatpolypart'##c.`contreat'
          tempvar ringer`polycount'
          qui gen double `ringer`polycount'' = `ringer'^`polycount'
          local ringerpoly `ringerpoly' `ringer`polycount''
        }
      }
    local ringerpoly `ringerpoly')
    local contreatinj   `latticepart2hash'##`contreatpolypart'
    local ringerinj   `latticepart1hash'#`ringerpoly'
    * the macros above will be used for constructing the contreatcontrols below and the main regression statement
    * this parens below closes the if contreat == contreat if statement
  }




  *******************************************************************************************************************
  ******************************  HANDLE CONTROLS, CRITICALLY, THE CONTREATCONTROLS  ******************************
  *******************************************************************************************************************

  ************ HANDLE REGULAR CONTROLS ****************
  * direct controls and direct fe don't require modification

  * controls for coarse case - i.e., where we interacted them with the treatment dummies. ccc_absorb puts them into absorb - only viable for continuous controls of course
  * the astute eye will notice that coarsecontrolinserter w/ categorical controls should partial out the main treatment effects
  * this does not occur b/c this statement is placed toward the end of the regression statement, meaning multicollinearity forces drops among the controls
  local coarsecontrolinserter
  local ccc_absorbinserter
  if "`coarsecohortcontrols'" != "" local coarsecontrolinserter i.`_i'#i.`_t'#(`coarsecohortcontrols')

  if "`ccc_absorb'" != "" {
    local ccc_absorbinserter i.`_i'#i.`_t'#c.(`ccc_absorb')
  }

  * handle interactive controls - these are closer to the object that wooldridge covers in his paper for covariates
  if "`interactivecontrols'" != "" {
    * center the variables by cohort i
    local interlist
    local iter 1
    foreach intercon in `interactivecontrols' {
      tempvar ic`iter'
      qui gen double `ic`iter'' = `intercon'
      local interlist `interlist' `ic`iter''
      local ++iter
    }
    qui gstats transform (demean) `interlist'    `weightinserter'   if   `ifstatement', by(`i') replace

    * there are two sets of objects to introduce:
    * interactions between the time variable and the raw interactive controls
    * interactions between the de-meaned interactive controls and the dummies for treated units' treatment periods
    * SI-20 (author-ruled, 2026-09-05): the cell interaction is built with 1.targetflag, not
    * i.targetflag, so the reference cell (_i = _t = targetflag = 0, the pooled untreated
    * observations) gets NO slope column. That column is exactly collinear with i.t#c.x plus
    * the treated-cell slopes (modulo the absorbed cohort means), so with i.targetflag the
    * engine's collinearity drop decided which column survived and the reported horizons,
    * their standard errors and the joint tests were solver- and order-dependent (the
    * reference-cell margin was non-zero wherever the covariate varies within cohort over
    * time, under weights, or on an unbalanced panel). With 1.targetflag the reference-cell
    * slope is partialled out through i.t#c.x, the reference margin is a structural zero
    * (0bno.) and the horizons are estimable functions of the model. (The earlier form
    * mirrored a term in some of Wooldridge's example files; it is not harmless here.)
    local interinject i.`_i'#i.`_t'#1.`targetflag'#c.(`interlist')  i.`t'#c.(`interactivecontrols')
  }


  * handle the time trend option
  local trendinjector
  * the trend enters through a mean-centered double copy of t: with raw calendar codes
  * the demeaned i.i#c.t columns lose all precision and the trends go numerically inert
  * (silently - rc 0, no caution; measured as low as t ~ 5e5, panel-dependent), and even
  * at 4-digit years the raw form drifts under affine shifts of t where the centered
  * form is exactly shift-invariant. centering is a pure reparameterization absorbed by
  * the unit intercepts; measured movement on ordinary panels is machine precision.
  * t is used as a continuous regressor NOWHERE else - every other formula use is
  * factor/absorb/grouping and must stay raw
  if "`timetrends'" == "timetrends" {
    tempvar tcen
    qui gstats sum `t' if `ifstatement', meanonly
    qui gen double `tcen' = `t' - r(mean)
    local trendinjector i.`i'##c.`tcen'
  }
	*note that we dont absorb this - pooled continuous variabels is a tricky matter sometimes


  ************** CREATE CONTREAT CONTROLS ***********
  * this covers how we control for continuous treatment/treatment intensity outside of treated/pre-treatment observations
  * we probably want to control for the continuous treatment variable in a flexible way, possibly in a way that varies by period or by cohort or by period separately in the treatment and control groups
  * this section enables that by creating appropriate controls (using alternate varnames or absorb, so as to remain invisible to margins)
  * with interactions to allow the specified degree of flexibility - eg, have controls vary by period vs by cohort vs etc.

  if "`contreatcontrols'" != ""   {
    * setup missing helper variables (only if not created earlier)
    if strpos("`contreatcontrols'" , "byTreated" ) != 0 & strpos( "`latticetype'" , "byTreated") == 0  {
      tempvar evertreated
      qui gen `evertreated' = `ttre' == .
    }
    if strpos("`contreatcontrols'" , "byTtre" ) != 0 & strpos("`latticetype'", "byTtre")==0{
      tempvar modified_ttre
      qui gen double `modified_ttre' = `ttre'
      qui replace `modified_ttre' = `minyearobserved' if `ttre' == .
      * minyearobserved is just any value not appearing in ttre; ttre is guaranteed to not include that value
    }

    * initialize. ccringerhead / ccringerheadabs / ccringertail are the pieces the two
    * spellings are composed from at the end of this block (see DF-1 there)
		local contreatcontrolinserterFE
		local contreatcontrolinserter
		local contreatcontrolinserterABS
		local ccringerhead
		local ccringerheadabs
		local ccringertail

    * handle basic, byTreated, byTtre, and byCohort in order since these are successively more intensive and supersede each other

    * basic just partials out the baseline effect of your contreat control. implicitly, this effect is pooled in your control group and reference period
    * this text initializes contreatcontrolinserterfe whenever basic is requested - with or without a lattice (SI-20, author-ruled 2026-09-05:
    *     the block used to be skipped when a lattice was present, so plain basic + lattice() inserted NO baseline dose control on the default
    *     path, the base level of contreat stayed margins-visible, and every ame-family row - main, pretreatment, each ES_ame horizon - was off
    *     by the baseline dose slope; with a lattice the lattice##ringer twin below is added so the base slope is partialled out within each lattice level)
    * and also initializes it this way if there is any chance that this FE object might get kicked outside of an HDFE absorb, in which case this section being included
    * is necessary in order to correctly manipulate the behavior of reg and poisson's colinearity drop rules (nothing in an HDFE is dropped b/c of colinearity with something not in the FE)
		* re: the nature of the manipulation, specifying null##ringer and lattice##ringer (as opposed to ringer alone or lattice#ringer) -- and doing so explicitly, rather than
		*     implicitly like in the other control types -- ensures that the base level of contreat gets partialled out. if not specified, it does not always get partialled out, and
		*     instead one of the contreat control coefficients gets dropped for colinearity instead
    if strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
      tempvar null
      qui gen `null' = 1
      local ccringerhead  i.`null'##`ringerpoly'
      local ccringerheadabs  i.`null'##`ringerpoly'
			if "`latticepart2hash'" != "" {
			  local ccringerhead `ccringerhead'  i.`latticevar'##`ringerpoly'
			  local ccringerheadabs  i.`latticevar'##`ringerpoly'
			}
			*using latticepart2hash as a check here is arbitrary and a different lattice byproduct could be used
    }

    * built on top of the initialized contreatcontrolisnerterFE -- add lattice components only if useful
    if strpos("`contreatcontrols'",  "byTreated") != 0  local ccringertail `ccringertail' i.`evertreated'`ringerinj'
    if strpos("`contreatcontrols'",  "byTreated") != 0 & "`latticepart2hash'" != "" local ccringertail `ccringertail'  i.`evertreated'`latticepart2hash'
    if strpos("`contreatcontrols'",  "byTtre") != 0  local ccringertail `ccringertail' i.`modified_ttre'`ringerinj'
    if strpos("`contreatcontrols'",  "byTtre") != 0  & "`latticepart2hash'" != "" local ccringertail `ccringertail'  i.`modified_ttre'`latticepart2hash'
    if strpos("`contreatcontrols'",  "byCohort") != 0  local ccringertail `ccringertail'  i.`i'`ringerinj'
    if strpos("`contreatcontrols'",  "byCohort") != 0  & "`latticepart2hash'" != "" local ccringertail `ccringertail'   i.`i'`latticepart2hash'

    * handle by period - this is in addition to, not in lieu of the others
    if strpos("`contreatcontrols'",  "byTime") != 0  local ccringertail `ccringertail' i.`t'`ringerinj'
    if strpos("`contreatcontrols'",  "byTime") != 0  & "`latticepart2hash'" != "" local ccringertail `ccringertail'     i.`t'`latticepart2hash'

    * ---- DF-1 (author ruling 2026-09-13): TWO spellings of this control set, because the
    * engines consume it two different ways.
    * `ringer' is a bit-exact double copy of the dose, so the pooled i.`null'##c.(ringer)
    * set is a numerical twin of columns the margins-visible block already carries through
    * its ## chain (c.`contreat', and with a lattice i.lattice and i.lattice#c.`contreat').
    * It adds no rank. Inside absorb() it nonetheless adds a slope that reghdfe COUNTS in
    * df_a, so the residual df came out too small and every att-row se, t, p and interval,
    * plus rmse, F, r2_a and ll_0, carried a spurious sqrt(df_old/df_new) of about
    * 0.14-0.25%. Measured three ways: the rank of the design is unchanged; a reg with
    * explicit dummies agrees with every spelling except the shipped one; and across the
    * suite 36 configs move the residual count with a residual sum of squares differing by
    * exactly 0 or at most 3e-15. Point estimates, ll, rss, r2 and N do not move.
    * Measured again on dta_cont while building this (df_a read back as rss/rmse^2):
    * basic + lattice(all 3) goes 753 -> 754, which is the count the ccnoabsorb, customvce
    * and unconditionalse spellings of the SAME specification already reported, so the fix
    * makes the absorbing path agree with the explicit ones rather than inventing a third
    * answer; basic + byTime goes 890 -> 891; basic + byTime + lattice goes 703 -> 709.
    * The 48-config battery moves exactly its two lattice configs and, at rtol 1e-8, in
    * rmse alone - no b, ll, rss, r2 or N anywhere.
    * THE RULE: absorb the FINEST ringer partition requested and only that one. The pooled
    * i.`null'## set is redundant as soon as ANY other ringer set is absorbed (the lattice
    * twin, or any byX - byX carries the lattice inside `ringerinj'), and the lattice twin
    * is redundant as soon as a byX one is. It is correct and free only on plain basic with
    * neither a lattice nor a byX. `ccringertail' is non-empty exactly when a byX was asked
    * for, because every byX adds its ringer term, which is what makes it the test here.
    * The manipulation the long comment above describes is NOT given up: the finest set
    * that survives partials the base dose slope out exactly as the pooled one did, within
    * a finer cell, and it is still spelled explicitly.
    * A caveat worth keeping: the suite's 36-config sweep was over basic + lattice configs
    * (only one suite config carries a byX at all), so the byX arms are thinner evidence.
    * On basic + byTime the two spellings agree to 3e-9 in rss and 1e-10 in b rather than
    * the 2e-15 the lattice arm shows, and basic + byTime + lattice on dta_cont is a
    * degenerate design in BOTH builds - its own coefficients move by O(1) between
    * regtol(9) and regtol(14) before the change as well as after, and its clustered VCE
    * comes back all zero - so nothing there is evidence about this edit either way.
    * ENGINE-CONDITIONED, and this is the part to keep: the thinning applies ONLY to what
    * goes inside absorb(). On every path that passes these terms as explicit REGRESSORS -
    * ccnoabsorb, unconditionalse, customvce, and the ppmlhdfe rc-3200 fallback - the
    * i.`null'## spelling is load-bearing for the collinearity DROP ORDER, since nothing in
    * an HDFE absorb is dropped for collinearity with something outside it but an explicit
    * regressor certainly is; those paths already report the correct residual count, and
    * they keep reading `contreatcontrolinserterFE', which is the full set exactly as before
    local contreatcontrolinserterFE `ccringerhead' `ccringertail'
    if "`ccringertail'" != "" local ccringerheadabs
    local contreatcontrolinserterABS `ccringerheadabs' `ccringertail'


    * in terms of interpreting these control measures, note that allow each cohort to have its own contreat control function can have some interesting consequences
    * more aggressive controlling can reduce the amount of cross-group comparisons that occur
    * adding byTime implicitly says "compare the contreat effect in treated units in period 1 to controls in period 1"
    * but allowing contreat effects to vary by cohort and period simultaneously would not have this effect


    * if always zero when not treated, the only thing to do is not control for anything! this results in no contreat controls
    * it is implemented via modifying the contreatinj statement and just not having contreatcontrols of any kind
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
    }

  }

  if "`ccnoabsorb'" == "ccnoabsorb"   {
		local contreatcontrolinserter `contreatcontrolinserterFE'
		local contreatcontrolinserterFE
		* and with it the absorb spelling (DF-1): under ccnoabsorb the whole set is an
		* explicit regressor, so nothing is left for absorb() to take
		local contreatcontrolinserterABS
	}






  ***************************************************************************************************
  ********************************         ESTIMATE THE MODELS        *******************************
  ******************************* (And Extract Information From Them) *******************************
  ***************************************************************************************************

  if "`verbose'" == "verbose" dis as text "Wooldid is now estimating its regression."

  * handle SEs
  * SI-8 (author-ruled 2026-08-25): granular-i guidance + maxvar preflight, both
  * skipped under safety(off). Counting is tag-based (no level lists in macros - the
  * huge-cardinality lesson). The preflight uses the EXACT treated interaction-cell
  * count (distinct (_i,_t) pairs, x2 under contreat for the dose slopes), a lower
  * bound on the regression's variable demand, so an early refusal can only fire
  * when the late r(907) was certain anyway.
  * P6 (author ruling 2026-09-07): gate on the USER's safety() value, not the working copy.
  * contreatwithin with contreat() rewrites the working copy to off (wooldid_validate), which
  * silently disabled this preflight for a user who never asked for safety(off) - and the
  * preflight is the difference between an early named error and a late r(907) mid-estimation.
  * The relaxation contreatwithin needs is of the estimability screens, not of a variable count.
  if "`safetyinput'" != "off" {
    tempvar si8tag
    qui gegen byte `si8tag' = tag(`i') if `ifstatement'
    qui count if `si8tag' == 1
    local si8ilev = r(N)
    cap drop `si8tag'
    qui gegen byte `si8tag' = tag(`_i' `_t') if `ifstatement' & `_i' != 0
    qui count if `si8tag' == 1
    local si8cells = r(N)
    cap drop `si8tag'
    qui gegen byte `si8tag' = tag(`t') if `ifstatement'
    qui count if `si8tag' == 1
    local si8tlev = r(N)
    cap drop `si8tag'
    local si8demand = `si8cells'*(1 + ("`contreat'" != "")) + `si8ilev' + `si8tlev' + 200
    if `si8ilev' > 100 {
      dis as text "Caution: the cohort variable `i' takes `si8ilev' distinct values. wooldid builds one interaction term per cohort-by-period cell, and runtime grows faster than linearly in that cell count - roughly its 1.3rd to 1.7th power, and about its 1.9th under lattice() - so a granular `i' can dominate runtime long before it approaches Stata's variable limits. Sample size itself is comparatively cheap. Coarsening `i' is normally the single largest speedup available: specify unit-level fixed effects via fe() and swap `i' to a coarser level of aggregation, such as cohorts or treatment onset groups. (safety(off) suppresses this check.)"
    }
    if `si8demand' > c(maxvar) {
      dis as error "Error: this specification requires at least `si8demand' variables (about `si8cells' treated cohort-by-period interaction terms) but maxvar is currently " c(maxvar) ", so estimation would fail with r(907)."
      dis as error "Consider specifying unit-level fixed effects via fe() and swapping `i' to a coarser level of aggregation, such as cohorts or treatment onset groups. (safety(off) suppresses this check.)"
      error 907
    }
  }

  * GROUP B: byCohort point-identification screen (see wooldid_amescreen mode bycoh).
  * Pre-estimation, always-on when byCohort+contreat are active (the pre-treatment-
  * period sufficiency errors' sibling, one layer down: rows can exist while the
  * baseline dose slope's identifying VARIATION does not).
  if strpos("`contreatcontrols'", "byCohort") != 0 & "`contreat'" != "" {
    tempvar amebtouse
    qui gen byte `amebtouse' = (`ifstatement')
    wooldid_amescreen , mode(bycoh) tcell(`_t') esamp(`amebtouse') contreat(`contreat') iorig(`i') tvar(`ttre')
    if r(fail) == 1 error 459
    cap drop `amebtouse'
  }

  * SI-A1 joint-test screen (author-ruled 2026-08-25): with analytic clustered SEs
  * clustered at the level of the cohort variable i or above, the saturated
  * cohort-by-period cell dummies force each cell's residuals to sum to zero within
  * its own cluster, so the CRVE cannot capture the within-cluster cross-horizon
  * variation the joint tests depend on (measured rejection rates 11-66% at nominal
  * 5%, persisting at large cohort counts and large cells). Runs under BOTH safety
  * modes and branches: error under default safety, caution + proceed under
  * safety(off) - so it gates on safetyinput, not safety (contreatwithin overwrites
  * the latter internally). customvce cluster-type counts as clustered (the Group A
  * containment screen's convention); adversarial is deliberately not gated here.
  * F-boot carve-out: under bootstrap the analytic joint tests are never computed
  * (their p-values are Wald tests on the bootstrap covariance, which this screen's
  * own error text recommends), so the screen does not apply
  * permute() carve-out (2026-09-05): the same - under permute() the joint tests are
  * permutation tests built by the shell, the other remedy this screen's text recommends
  * PS-64 follow-through (audit 2026-09-10): pretrendpower() is the SECOND entry point to the
  * same statistic - wooldid_ptp computes its joint pre-trend test from the identical analytic
  * e(V_es_*) - and the ruling says one statistic, one rule, in BOTH. The screen fired only
  * under jointtests, so a pretrendpower() run at cohort-level clustering printed the oversized
  * test with nothing said. The gate now covers either entry point; the message names whichever
  * one asked for it.
  if ("`jointtests'" != "" | "`pretrendpower'" != "") & "`bootstrap'`permute'" == "" {
    local jtclusvars
    if "`cluster'" != "" local jtclusvars `cluster'
    local jtcvw = word(`"`customvce'"', 1)
    if "`jtclusvars'" == "" & strlen("`jtcvw'") >= 2 & strpos("cluster", "`jtcvw'") == 1 {
      local jtclusvars : subinstr local customvce "`jtcvw'" "", word
    }
    local jtcoarse
    foreach jtcv of local jtclusvars {
      cap confirm variable `jtcv'
      if _rc != 0 continue
      tempvar jtt
      qui gegen byte `jtt' = tag(`i') if `ifstatement'
      qui count if `jtt' == 1
      local jtni = r(N)
      cap drop `jtt'
      qui gegen byte `jtt' = tag(`i' `jtcv') if `ifstatement'
      qui count if `jtt' == 1
      local jtnic = r(N)
      cap drop `jtt'
      if `jtni' > 0 & `jtnic' == `jtni' & "`jtcoarse'" == "" local jtcoarse `jtcv'
    }
    if "`jtcoarse'" != "" {
      * PS-64 (author ruling 2026-09-08): compute the test with a named caution instead of
      * refusing. The comparator packages run their pre-trend tests at panel-level
      * clustering with no gate, and wooldid's own help twice recommends the clustering
      * that trips this condition, so a hard error sent users into a documented dead end.
      * The caution states what is wrong with the covariance; the number is the user's to
      * weigh. This supersedes the earlier ruling to suppress the same statistic under
      * pretrendpower(): one statistic, one rule, in both entry points.
      dis as text "Caution: `=cond("`jointtests'" != "", "jointtests", "pretrendpower()")' was combined with analytic clustered SEs clustered at the level of the cohort variable i or above (`i' is nested within the clustering variable `jtcoarse'). The saturated cell dummies leave each cell's residuals with no within-cluster variation to identify, so the joint test's covariance is degenerate and its p-value should not be read as a calibrated test. wooldid's randomization inference (permute()) or bootstrap-based clustered SEs give a joint test that is valid under this design, as does clustering strictly finer than the cohort-by-period cells."
    }
  }

  local vcespecifier
  if "`robust'" == "robust" local vcespecifier vce(robust)
  if "`cluster'" != "" local vcespecifier cluster(`cluster')

  * Regression!
  * unconditionalse and customvce punts us to using reg or poisson; otherwise, use reghdfe or ppmlhdfe
  * the transition is fairly easy in that you can get the same results by placing any absorbed FE before the treatment variable specification in the reg equation
  * since colinear terms get dropped in the order in which they appear (mostly)

  *Re regtol - it is set to 9 by default, is is the reghdfe default.
  if "`unconditionalse'" == "" & "`customvce'" == "" {
    * ppmlhdfe specification is different from reghdfe b/c it is more hostile to partialling out continuous variables
    * ppmlhdfe requires d for predict, which margins requires
    if "`poisson'" == "poisson" {
      if "`verbose'" != "" dis "ppmlhdfe `y' `contreatcontrolinserter' `trendinjector'    i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'   `contreatcontrolinserterABS'     `ccc_absorbinserter'  ) `vcespecifier' tol(1e-`regtol') d "
      qui cap ppmlhdfe `y' `contreatcontrolinserter' `trendinjector'    i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'    `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'    `contreatcontrolinserterABS'  `ccc_absorbinserter' ) `vcespecifier' tol(1e-`regtol') d
      if _rc == 3200 {
        dis as text "Caution: initial estimation of ppmlhdfe failed, possibly due to inclusion of an absorbed categorical fixed effect with only 1 level."
        dis as text "Falling back to estimation with ppmlhdfe, but with any included time-trend controls, ccc_absorb controls, and/or continuous treatment controls pulled out of the absorb() statement and placed in the main regression."
        dis as text "In principle, this fallback solution should usually yield the same results, but will be computationally less efficient."
        dis as text "If estimation still fails, the issue could be either with the fe() specified or a result of having few cohorts/time-periods."
        qui ppmlhdfe `y' `contreatcontrolinserter' `contreatcontrolinserterFE'  `ccc_absorbinserter'  `trendinjector'  i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'   ) `vcespecifier' tol(1e-`regtol') d
      }
      else {
        if _rc != 0 error _rc
      }
    }
    else{
      if "`verbose'" != "" dis "reghdfe `y' `contreatcontrolinserter' `trendinjector' i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'   `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'     `contreatcontrolinserterABS'   `ccc_absorbinserter'  ) `vcespecifier' tol(1e-`regtol')  "
      qui reghdfe `y' `contreatcontrolinserter' `trendinjector'  i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'     `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'     `contreatcontrolinserterABS'   `ccc_absorbinserter'  ) `vcespecifier' tol(1e-`regtol')
    }
  }
  else{
    * run the reg/poisson version, places all absorbed FE at the front of the reg equation
    if "`customvce'" != "" local vcespecifier vce(`customvce')
    local regtype reg
    if "`poisson'" == "poisson" local regtype poisson
    if "`verbose'" != "" dis "`regtype' `y' i.`i' i.`t' `fe'   `trendinjector' `contreatcontrolinserter' `contreatcontrolinserterFE'   `ccc_absorbinserter'   i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' ,  `vcespecifier'  "
    qui `regtype' `y' i.`i' i.`t' `fe'   `trendinjector' `contreatcontrolinserter' `contreatcontrolinserterFE'   `ccc_absorbinserter'   i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'    `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' ,  `vcespecifier'
  }

  * save the model; will be useful later
  * register the data-variable obligation with the shell's error-path cleanup epilogue
  * BEFORE storing: _ppmlhdfe_d already exists in the data at this point, and est sto
  * itself can fail (a real trigger: the estimates directory at Stata's 300-store
  * limit), which previously leaked the variable through the unregistered window
  if "`poisson'" == "poisson" & "`unconditionalse'" == "" & "`customvce'" == "" global WOOLDID_CLEANUP $WOOLDID_CLEANUP ppml
  qui est sto `wdfm'
  global WOOLDID_CLEANUP $WOOLDID_CLEANUP `=cond("`wdrepmode'" != "", "repmodel", "fullmodel")'
  if "`verbose'" == "verbose" dis as text "Wooldid has finished regression estimation and is preparing to begin margins calculations."

	* check on dropped obs
	* BS-37 (author ruling 2026-09-07): `dropcount' is the SAME exact-one expression the
	* sample rule uses, so it cannot see the rows that rule alone removed - a user who
	* typed "if flag" with flag in {1,2} specified 72 rows, got 60, and was told 0 were
	* dropped, contradicting the e(N_dropped) definition in the help. Add those rows. No
	* double count: a row with `ifvar' outside {0,1} has (`if') != 1 in double as well, so
	* `dropcount' is false there. wdvifexcl is 0 whenever the if() expression is a
	* comparison - i.e. on every ordinary if(), which is why no blessed baseline moves
	qui count if `dropcount' & e(sample) == 0
	local droppedobs = r(N) + `wdvifexcl'

  * grab some items from the model - sample indicator variable, residual degrees of freedom, etc
  tempvar esamp
  qui gen byte `esamp' = e(sample)
  local df_r = e(df_r)
  if "`poisson'" != "" {
    qui count if `esamp' == 0 & `targetflag' == 1 & `w' > 0 & `ifstatement'
    if r(N) > 0 dis as text "Caution: ppmlhdfe excluded `r(N)' treated observation(s) from the estimation sample (separation, or a singleton fixed-effect group), so their cohort-period cells contribute nothing to the reported effects."
  }
  * one cluster leaves e(df_r) = 0 and the run then died r(125) inside margins under a
  * caution blaming continuous treatment; stop here with the reason (bug sweep BS-35).
  * fast mode asks margins for no variance and may proceed
  if "`fastnose'" == "" & `df_r' < 1 {
    if e(N_clust) < 2 dis as error "Error: cluster(`cluster') leaves `=e(N_clust)' cluster(s) in the estimation sample; cluster-robust inference needs at least 2. Use a finer cluster variable, robust, or fastnose."
    else dis as error "Error: the model leaves no residual degrees of freedom (e(df_r) = `df_r'), so no variance can be estimated."
    error 459
  }

  * AMESE R1 (standing-issues family A): the structural cell-containment predicate.
  * Fires only with contreat + a cluster-type VCE, never under fastnose (no variances),
  * impvalues (IM dispersion inference is valid under containment), or unconditionalse
  * (different estimand; excluded from suppression by design). When it fires, every
  * downstream ame-family capture nulls its se/t/p/ci cells (the existing sentinel
  * machinery renders them as missing SEs), the ame-side V_es tests are skipped, and
  * e(V_es_ame*) is withheld while e(b_es_ame*) still posts.
  local amescrub 0
  if "`contreat'" != "" & "`impvalues'" == "" & "`unconditionalse'" == "" & "`fastnose'" == "" {
    local amescrubclus
    if "`cluster'" != "" local amescrubclus `cluster'
    local amecvw = word(`"`customvce'"', 1)
    if "`amescrubclus'" == "" & strlen("`amecvw'") >= 2 & strpos("cluster", "`amecvw'") == 1 {
      local amescrubclus : subinstr local customvce "`amecvw'" "", word
    }
    if "`amescrubclus'" != "" {
      wooldid_amescreen , mode(detect) icell(`_i') tcell(`_t') tflag(`targetflag') esamp(`esamp') clusvars(`amescrubclus') contreat(`contreat')
      local amescrub = r(fire)
    }
  }

  * AMESE R4: catch Stata's ereturn-post variance screen zeroing the whole VCE (the
  * quiet estimation call swallows the only warning) and SAY SO. Detection and the
  * caution are unconditional outside fast mode.
  *
  * The att-rescue refit that used to sit here - a re-estimation of the same
  * specification through regress with the fixed effects written out as explicit
  * dummies - was REMOVED on 2026-09-08 by author ruling. It was ruled in on
  * 2026-08-24 inside the larger AMESE package, where its runtime cost was never put
  * separately. What it bought was the subset of standard errors that survive a
  * specification whose cell-level containment has already made the ame inference
  * invalid; what it cost was a solve over hundreds of collinear dummy columns on
  * every such run, which delays the moment the user learns the specification is
  * broken, and which is numerically far noisier than the partialled fit it replaced
  * (measured 2026-09-08: reordering the input rows moves the explicit-dummies answers
  * by up to 1.6e-8 relative, against 1e-13 for the partialled fit; research record
  * review_process/research_2026-09-08/fastnose_precision.md). The run now reports the
  * all-zero covariance and continues with the standard errors it actually has.
  if "`fastnose'" == "" {
    wooldid_amescreen , mode(vzero) r1fired(`amescrub') contreat(`contreat')
    if r(vzero) == 1 {
      dis as text "Caution: the estimation command returned an all-zero variance-covariance matrix, so no standard error in this run is identified. The reported estimates stand; every standard error, test statistic, p-value and interval is missing. This is a property of the specification, not a numerical accident: `=cond("`contreat'" != "", "see the cell-containment caution above for the mechanism and for the alternatives", "the alternatives are")' (impvalues dispersion inference, a finer-than-cell clustering, or conditional robust standard errors)."
    }
  }

  *****************************************************************************************************************
  ******************************  PRODUCE THE WEIGHTS FOR CREATING VARIOUS OUTCOMES  ******************************
  *****************************************************************************************************************
  * special weights are required to go beyond a simple ATT; this generates them
  * note: constructed after estimation, on the realized estimation sample (esamp), so that cell/cohort
  * normalizations only count observations that actually entered the regression (estimators can drop
  * observations - e.g. reghdfe singletons, ppmlhdfe separation - after the initial sample screens)

  * handle custom weight multipliers that the user specifies
  local customweighttext
  if "`customweightmultiplier'" != "" local customweighttext * `customweightmultiplier'

  ***** begin manufacturing weights to produce averages of interest
  * 4 weight concepts to target
  * att is the "do nothing but use sample weights" weight
  * att_i says give the average of cohort-specific effects, obtaining said effects using sample weights
  * att_it says "calculate cohort-period specific effects using weights, then average those"
  * att_itime says "calculate cohort-period specific effects using weights, then average those by cohort, then average the cohort averages"

  foreach item in tag posttreatment   w_per_it w_per_i  t_per_i   w_att  w_att_it w_att_i w_att_itime   {
    tempvar `item'
  }

  * create necessary supporting object for computing att_i and att_itime weights
  * this object is 1 if in the treated period for treated units,
  * 2 if in the pre-treatment period for treated units and we are reporting effects for that time period
  * 0 if in the pre-treatment period for treated units and we are not reporting effects for that period (i.e., the reference period)
  * 3 if untreated and in sample
  * 4 if out of sample
  if "`att_i'" == "att_i" | "`att_itime'" == "att_itime" {
    qui gen `posttreatment' = `t' >= `ttre' if `esamp' == 1
    qui replace `posttreatment' = 2 if `t' < `ttre' & `ttre' != . & `_i' != 0 & `_t' != 0 & `esamp' == 1
    qui replace `posttreatment' = 3 if `ttre' == .  & `esamp' == 1
    qui replace `posttreatment' = 4 if `posttreatment' == .
  }

  * this variable helps count the number of i-t cells that exist in various subsets - useful for att_itime weights
  if "`att_itime'" == "att_itime" {
    qui gegen `tag' = tag(`i' `t') if  `esamp' == 1
    qui replace `tag' = 0 if `tag' == .
  }

  * generate some inputs for weights calculation
  * w_per_it contains the amount of weight falling into each i-t cell
  * w per i gives the total amount of sample weight associated for each i falling into the categories specified by posttreatment
  * t per i is the total number of periods appearing for i in each posttreatment-defined period
  * in each case, we only make the input if needed for a weight type the user actually requested
  * that said, we always produce w_att_it whenever a non att weight is requested, since they are needed for the event horizon estimates and sometimes for histogram
  if "`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" qui gegen `w_per_it' = total(`esamp' * `w')  , by(`i' `t')
  if "`att_i'" == "att_i"  qui gegen `w_per_i' = total(`w'*`esamp')      , by(`posttreatment' `i')
  if "`att_itime'" == "att_itime" qui gegen `t_per_i' = total(`tag')    , by(`posttreatment' `i')

  * att requires no special weights - it gets to stay as the given weight (w), albeit with the customweight text multiplier coming through
  * dividing w by w_per_it weights all i-t cells equally, but still lets w be influential within i-t cells (att_it)
  * dividing w by total w in i weights all i equally, but still lets w be influential within i-t cells and across i-t cells within i (att_i)
  * dividing w by total w in i-t cells * number of t in the i, yields w_att_itime. this lets w be influential inside i-t cells, but then
  *        weights every i-t cell equally when contributing to estimates of i-specific effects, which are themselves then averaged

  if "`att_i'" == "att_i" {
    qui gen double `w_att_i' = `w' `customweighttext'/`w_per_i'
    qui replace `w_att_i' = 0 if `w_att_i' == .
    qui drop `w_per_i'
    qui drop  `posttreatment'
  }

  if "`att_itime'" == "att_itime" {
    qui gen double `w_att_itime' = `w' `customweighttext' /(`w_per_it'*`t_per_i')
    qui replace `w_att_itime' = 0 if `w_att_itime' == .
    qui drop `tag' `t_per_i'
    qui cap drop  `posttreatment'
  }

  if "`att_it'" == "att_it" | "`att_itime'" == "att_itime"  | "`att_i'" == "att_i" {
    qui gen double `w_att_it' = `w' `customweighttext'/`w_per_it'
    qui replace `w_att_it' = 0 if `w_att_it' == .
    qui drop `w_per_it'
  }

  * generate w_att whenever att is requested
  if "`att'" == "att" qui gen double `w_att' = `w' `customweighttext'

  * some comments on the weights:
  * provided each i is always wholly contained in a given sg, we have no issues using these weights for sg purposes
  * since we don't allow fw type weights, weight scale here actually doesn't matter - no need to multiply weights by scalar adjustment factors
  * for event study estimates, note that att_i and att_itime are not conceptually all that sensible --
  *       we substitute att_it for them since it retains a spirit of 'all same across i'

  * save the name of the current frame  - will be useful going forward
  qui frame
  local defaultframe `r(currentframe)'

  * create the matrix relyrcolumn containing all relative periods of interest appearing in the regression - needed for constructing the results
  * do this by processing the column names of e(b)
  if `espre' + `espl' > 0 {
    local relyrlist
    tempname relyrcolumn
    local bvars: colnames e(b)
    local horizlist
    local horizonend = `horizonstarter' + `espre' + `espl' -1
    forval horiz = `horizonstarter'(1)`horizonend' {
      if regexm("`bvars'","[0-9].`_i'#`horiz'.`_t'#1.`targetflag'") {
        * the above regular expression resolves to 1 only if you see the following:
        * (anynumber not followed by a o or a b)._it#(horiz)._t#1.targetflag   <- so a valid, included in the model estimate at the given horizon
        local relyr = `horiz'-`min_t'-1
        local horizlist `horizlist' `horiz'
        if "`relyrlist'" != "" local relyrlist `relyrlist' , `relyr'
        if "`relyrlist'" == "" local relyrlist `relyrlist' `relyr'
      }
    }
    if "`relyrlist'" == "" {
      * no treated unit realizes any requested relative period: name the window and what the
      * sample actually contains instead of letting matrix () die with a bare r(198) (BS-33)
      * two causes, and the message must not blame the wrong one: either the sample has no
      * treated observation at any requested relative period, or it has them and the fitted
      * model carries no matching coefficient (the treatment terms were absorbed or dropped)
      local wdreal
      local wdhit 0
      qui levelsof `_t' if `_t' > 0 & `_t' < ., local(wdtl)
      foreach v of local wdtl {
        local wdreal `wdreal' `=`v' - `min_t' - 1'
        if `v' >= `horizonstarter' & `v' <= `horizonend' local wdhit 1
      }
      if `wdhit' {
        dis as error "Error: the fitted model carries no estimable treatment-effect coefficient at any of the requested event-study relative periods (`=`horizonstarter' - `min_t' - 1' through `=`horizonend' - `min_t' - 1'), although the estimation sample does contain those periods."
        dis as error "This means the cohort-by-period treatment terms were absorbed or dropped by the specification rather than estimated: check absorbed fixed effects (ccc_absorb() in particular, which can absorb the treatment interactions themselves), fe(), and any control set collinear with the treatment cells."
        error 459
      }
      dis as error "Error: no treated unit is observed at any of the requested event-study relative periods (`=`horizonstarter' - `min_t' - 1' through `=`horizonend' - `min_t' - 1'). Relative periods realized by treated units in the estimation sample: `wdreal'."
      dis as error "This usually means the onset period is missing for every treated cohort (an onset-period gap in the panel, or an onset that falls off the panel's time grid, such as an even-year onset in a biennial panel). Widen esprelength()/espostlength() to reach a realized period, or place ttre on the panel's grid."
      error 459
    }
    matrix `relyrcolumn' = (`relyrlist')'
  }


  * here, we check and see if the main treatment effects have been dropped from the regression for multicollinearity - my tests for this are imperfect, but still useful
  * amevoid (SI-A4/R6) and the scrub column floor are consumed downstream
  * unconditionally, so they are initialized outside the safety gate
  local amevoid 0
  local amescrublo 2
  * the collinearity screens live in wooldid_ccscreen (worker compiled-size headroom,
  * 2026-09-07): identical behaviour, the gate kept here so safety(off) never enters
  if "`safety'" != "off" wooldid_ccscreen , iv(`_i') tv(`_t') tf(`targetflag') wt(`w') es(`esamp') ifs(`ifstatement') mint(`min_t') ui(`i') ut(`t') lat(`lattice') ct(`contreat') ics(`interlist') icn(`interactivecontrols') `verbose' `poisson' `fastnose'




  **** extract important parameters from the model  - N, R2, etc
  * first, extract general parameters of interest from the regression
  local N = e(N)
  local Nclus = e(N_clust)
  local ll = e(ll)
  local ll_0 = e(ll_0)
  if "`regtype'" != "poisson" local rmse = e(rmse)
  if "`regtype'" != "poisson" local rss = e(rss)
  if "`poisson'" == "poisson" {
    local r2pseud = e(r2_p)
    local chi2 = e(chi2)
  }
  else {
    local r2 = e(r2)
    local r2adj = e(r2_a)
    local r2_within = e(r2_within)
    local r2_withinadj = e(r2_a_within)
    local F = e(F)
  }






  ***************************************************************************************************
  ***********************  F3: ATT(g,t) BUILDING BLOCKS (gteffects)  ********************************
  ***************************************************************************************************
  * runs inside the "clean e() window": the stage-1 model stored above is still the ACTIVE e()
  * here, so this reads e(b)/e(V)/e(sample) directly - no restore, no e-class command, and no
  * effect on anything downstream. binary designs only (wooldid_validate refuses contreat):
  * the <g>._i#<t>._t#1.targetflag coefficients ARE the cohort-by-relative-period effects on
  * the regression scale. the name parse and table assembly live in the wooldid_gtx mata
  * function at the bottom of this file to conserve worker compiled-size headroom.
  local gtn 0
  if "`gteffects'" != "" {
    tempvar gtin gttag gtnob gtw1 gtw2 gtw3 gtw4
    qui gen byte `gtin' = e(sample) & `targetflag' == 1
    qui gegen double `gtnob' = total(`gtin') if `gtin', by(`_i' `_t')
    local gtj 0
    foreach gtwv in `w_att' `w_att_it' `w_att_i' `w_att_itime' {
      local ++gtj
      cap confirm variable `gtwv'
      if _rc == 0 qui gegen double `gtw`gtj'' = total(`gtwv') if `gtin', by(`_i' `_t')
      else qui gen double `gtw`gtj'' = .
    }
    qui gegen byte `gttag' = tag(`_i' `_t') if `gtin'
    tempname gtraw gtsed gtm
    cap mkmat `_i' `_t' `ttre' `t' `gtnob' `gtw1' `gtw2' `gtw3' `gtw4' if `gttag' == 1, matrix(`gtraw')
    if _rc == 0 {
      cap matrix `gtsed' = vecdiag(e(V))
      if _rc != 0 matrix `gtsed' = J(1, colsof(e(b)), .)
      local gthstart = 1 + `min_t' - `espre'
      if "`esfixedbaseperiod'" == "" local gthstart 1
      local gthend = `gthstart' + `espre' + `espl' - 1
      mata: wooldid_gtx("`_i'", "`_t'", "`targetflag'", "`gtraw'", "`gtsed'", `min_t', `gthstart', `gthend', "`df_r'", "`gtm'", `wdlvl')
    }
    qui drop `gtin' `gtnob' `gtw1' `gtw2' `gtw3' `gtw4' `gttag'
  }

  **************************************************************************************************
  *****************************         PRODUCE SUMMARY STATS       ********************************
  **************************************************************************************************
  if "`summarystats'" == "summarystats" {

    * include among the summary stats information on size of reference period where relevant
    if "`esfixedbaseperiod'" == "" &  `espre' > 0  {
      qui count if e(sample) & `_t' == 0 & `ttre' < .
      local Npre_trebase = r(N)
      qui count if e(sample) & `_t' != 0 & `ttre' < . & `t' < `ttre'
      local Npre_treused = r(N)
      local Npre_rat = round(100*`Npre_trebase' / (`Npre_trebase' + `Npre_treused'),.01)
      tempvar overt
      qui gen `overt' = `_t' == 0
      qui total `w' if e(sample) & `ttre' < . & `t' < `ttre',over(`overt')
      drop `overt'
      local Wpre_trebase = e(b)[1,2]
      local Wpre_treused = e(b)[1,1]
      local Wpre_rat = round(100*`Wpre_trebase' / (`Wpre_trebase' + `Wpre_treused'),.01)
      qui est restore `wdfm'
    }


    * get some summary statistics on the outcome variable and the continuous treatment variable (if existent) in the estimation sample
    * this code block is inefficient / a result of some legacy choices
    * the summary-statistics table lives in wooldid_sumstat (worker compiled-size headroom,
    * 2026-09-08): identical behaviour; the matrix name is minted here so it survives the
    * call and can be posted as e(summarystats) below
    tempname summarystats
    wooldid_sumstat , ss(`summarystats') y(`y') tvar(`t') ttrevar(`ttre') ct(`contreat') ctc(`contreatcontrols') ctw(`contreatwithin') wi(`weightinserter')
    * NOTE the row names of e(summarystats) are set INSIDE wooldid_sumstat, not here. A
    * `matrix rownames `summarystats' = `rownamessum'' used to sit on this line and was
    * DELETED as inert (D-5, author ruling 2026-09-11): `rownamessum' is built inside that
    * program and is not among the locals it returns, so the line expanded to
    * `matrix rownames X = ' - measured a silent no-op in Stata 19 MP, neither an error nor
    * a clearing of the existing names. Do not delete wooldid_sumstat's own rownames call on
    * the belief that the caller does it: e(summarystats) would post with r1..r8 instead of
    * the y_*/c_* names that are the only thing identifying which sample each row describes.
    * (The `tempname summarystats' above also overwrites the OPTION macro of the same name,
    * which is why the post gate below tests != "" rather than == "summarystats")
  }


  ****** set some parameters that margins will rely on later
  * Poisson reporting scale. Both forms name the SAME prediction - the linear index -
  * so the documented xb-vs-predict(n) scale choice (see poisexpresults) is unaffected by
  * which one is used. They differ only in how margins differentiates it:
  *   expression() implies nochainrule, i.e. a NUMERICAL Jacobian (one interpreted pass
  *     over the estimation sample per coefficient). StataCorp recommend that posture
  *     after community-contributed estimators such as ppmlhdfe.
  *   predict()    reaches the ANALYTIC Jacobian.
  * For a plain dydx contrast of xb the margin is exactly linear in b, so the analytic
  * Jacobian is the EXACT one: measured against the exact identity J*b == margin, the
  * predict() form agrees to 2.3e-16 while the expression() form - a finite-difference
  * approximation - agrees only to 3.4e-13. So on dydx sites predict() is both faster
  * (2-10x end to end) and more accurate, and the change moves reported variances toward
  * the exact value rather than away from it.
  * The ey* operators (eyex, eydx) are a different matter and KEEP expression(): both code
  * paths take ln of the prediction, the poisson prediction is the sign-unconstrained
  * linear index, and which form aborts r(459) is data-dependent and goes in BOTH
  * directions. Preserving today's rc behaviour there - including its failures - is
  * deliberate; see the efficiency review and SI-20. dyex (dy/dln x) logs the REGRESSOR,
  * not the prediction, so it takes the predict(xb) form like dydx (research record
  * review_process/research_2026-09-06/F032.md: J*b == margin to 5.2e-16 under this form
  * against 2.7e-9 under expression(); 1.9-28x faster; author ruling 2026-09-06)
  local poismarginstext
  local poismarginstextC
  if "`poisson'" == "poisson" & "`poisexpresults'" == "" {
    local poismarginstext  expression(predict(xb))
    local poismarginstextC expression(predict(xb))
    if "`semielasticity'"         == "dydx" local poismarginstext  predict(xb)
    if "`contreatelasticitytype'" == "dydx" | "`contreatelasticitytype'" == "dyex" local poismarginstextC predict(xb)
  }




  ******************************************************************************************
  *****************************         PRODUCE HISTOGRAMS       ***************************
  ******************************************************************************************

  * F-collapse: margins evaluates over every row, but inside the collapse's gate every
  * quantity it computes is a cell-level weighted average and the derivative is constant
  * within an (_i, _t, targetflag) cell. So the margins calls below can be handed one row
  * per cell carrying the cell's SUMMED weight. `fcon' is wooldid_validate's realized
  * decision (on by default, off under fullprecisionslow or any ineligible option).
  * `preserve' is inline here on purpose - issued inside wooldid_fcol it would unwind at
  * that program's return and every margins call would silently run on full data. The
  * restore twin is likewise inline.
  local fcnoes
  local fconn
  * subgroup runs additionally screen each subgroup margins call for a numerically
  * singular VCE direction (wooldid_fcdeg): `fcsg' arms the screen, `fcguard' records
  * its outcome for e(wdid_fcguard), `fcrerun' is its per-call verdict
  local fcsg 0
  local fcguard
  local fcrerun 0
  if "`fcon'" == "1" {
    preserve
    * with a continuous treatment the surviving row also carries the DOSE, which varies
    * inside a cell, so two things change. (1) the user's cohort variable joins the key: `_i'
    * is `i' only where `_t' > 0 and 0 elsewhere (see its construction above), so without `i'
    * the reference-period rows pool across every cohort sharing an onset, and a read weighted
    * by w_att_it then consumes a mean built under w_att across cohorts whose weights are not
    * proportional (measured: a pre-trend joint p moving .0956 -> .1112). (2) the dose on the
    * surviving row is replaced by the cell's estimand-weighted mean, which is what makes every
    * margins read exact rather than a one-row draw of a within-cell-varying regressor.
    local fckeyx
    local fccont
    * F062 (quality review 2026-09-05, author ruling 2026-09-06): impvalues reads per-cluster
    * effects over(cluster), which the collapse reproduces exactly once the cluster variable
    * is part of the cell key (measured: a cluster crossed with i is zero at 1e-12 with the
    * key change, broken without it; a cluster nesting i is exact either way). The cell
    * count grows only when the cluster is crossed with i
    if "`impvalues'" != "" local fckeyx `fckeyx' `cluster'
    if "`contreat'" != "" {
      local fckeyx `fckeyx' `i'
      * the ringer (wooldid's double copy of the dose) is averaged the same way: wherever a
      * ringer set is in the model at all - absorbed as `contreatcontrolinserterABS' or
      * carried as an explicit regressor in `contreatcontrolinserterFE' - its level shifts
      * the round-off of every prediction. Unconditional on purpose: the ringer is built
      * whenever contreat() is given, and averaging it costs nothing where it is unused
      local fccont contreat(`contreat' `ringer')
    }
    wooldid_fcol `esamp' , key(`_i' `_t' `targetflag' `ifvar' `t' `ttre' `subgroups' `fckeyx') wvars(`w_att' `w_att_it' `w_att_i' `w_att_itime') `fccont'
    local fconn = r(ncell)
    * the collapsed rows carry no e(sample) marker, so every call in the window needs
    * noesample; rows in memory ARE the estimation sample (wooldid_fcol keeps on it)
    local fcnoes noesample
    if "`subgroups'" != "" local fcsg 1
    if "`verbose'" == "verbose" dis as text "Wooldid is evaluating the margins calls on one row per cell (`fconn' cells; see fastcollapse / fullprecisionslow in the help file)."
  }
  else if "`verbose'" == "verbose" {
    * F037 (2026-09-06): the user's own switch is the reason whenever it is given, with the
    * ineligible option(s) noted beside it; the option-naming text covers the silent skip
    if "`fullprecisionslow'" != "" {
      local fcwhy (fullprecisionslow)
      if "`fcskip'" != "" local fcwhy (fullprecisionslow; the default cell collapse would not apply with `fcskip' either)
      dis as text "Wooldid is evaluating the margins calls on the full data `fcwhy'."
    }
    else if "`fcskip'" != "" dis as text "Wooldid is evaluating the margins calls on the full data: the default cell collapse does not apply with `fcskip' (see fastcollapse in the help file)."
  }

  * since semielasticities are allowed, we need text for inserting that info into plots
  local effecttype Effects
  local Ceffecttype Effects
  * BS-6 (author ruling 2026-09-07): label what the user ASKED for. Under the poisson
  * xb-scale remap the operator locals hold dydx/dyex, but the quantity they produce IS the
  * requested semi-elasticity/elasticity, so the raw request is the true label. Off the
  * poisson path raw and operator are identical and nothing moves.
  if "`semielasraw'" == "eydx" local effecttype Semi-elasticities
  if "`cetyperaw'" == "eyex" local Ceffecttype Elasticities
  if "`cetyperaw'" == "eydx" local Ceffecttype Semi-elasticities (eydx)
  if "`cetyperaw'" == "dyex" local Ceffecttype Semi-elasticities (dyex)

  * for histograms of i-level effects, we do the following:
  * ask for a margins estimation to give us treated cohort specific effects without SEs
  * stick them into a matrix
  * use the matrix to populate a dataset that we will use to make the histogram
  if "`histogramcohorteffects'" != "" {
    * the histogram diagnostic lives in wooldid_hist (worker compiled-size headroom,
    * 2026-09-07): identical behaviour; the matrix name is minted here so it survives the
    * call and can be posted as e(histogramestimates) below
    tempname ___hist
    wooldid_hist , h(`___hist') iv(`_i') tv(`t') tf(`targetflag') ttrevar(`ttre') wtype(`weight') semi(`semielasticity') ifs(`ifstatement') dfr(`defaultframe') wa(`w_att') wai(`w_att_i') wait(`w_att_it') ct(`contreat') cet(`contreatelasticitytype') pm(`poismarginstext') pmc(`poismarginstextC') fcn(`fcnoes') sp(`saveplots') elbl(`effecttype') celbl(`Ceffecttype') `att' `att_it' `att_i' `att_itime' `verbose' `replace'
  }




  ***************************************************************************************************
  *****************************         GATHER MAIN EFFECTS AND HORIZON EFFECTS       ****************************
  ******************************************************************************************
  * the text below ensures margins produces p-values using t-tests with same dof as reghdfe; default is a z test w/ dof = infinity
  local marginsdof df(`df_r')

  * ensure that we request unconditional SEs where sought
  if "`unconditionalse'" == "unconditionalse" local marginsdof `marginsdof' vce(unconditional)

  * E1b: the two whole-sample pre-treatment horizon joint tests below hard-coded
  * df(`df_r') where every other margins call threads `marginsdof', so under
  * unconditionalse they silently used the conditional (delta) VCE. they cannot take
  * `marginsdof' itself: under impvalues `marginsdof' becomes `nose', which would leave
  * the joint test with no VCE at all (unreachable today - jointtests and impvalues are
  * mutually exclusive - but a latent trap), so the dof/vce pair is built separately
  local jointmarginsdof df(`df_r')
  if "`unconditionalse'" == "unconditionalse" local jointmarginsdof `jointmarginsdof' vce(unconditional)

  * F9: fast (no-SE) mode - every margins call runs with nose, so margins skips the
  * delta-method VCE (the dominant cost of a margins call in wooldid's design). r(table)
  * then comes back with ONE row instead of nine, which is why every six-row slice below
  * is padded through `f9tab' first. (the jointmarginsdof sites above are unreachable
  * under fastnose: jointtests + fastnose is refused at parse time)
  if "`fastnose'" != "" local marginsdof nose
  * the subgroup-by-horizon margins call keeps its VCE even in fast mode: the SI-7
  * estimability drop reads margins' omission marker, which nose suppresses, so under
  * fastnose a structural-zero horizon was reported (and ranked by permute()) and the
  * observed and replication column sets diverged, aborting bootstrap()/permute()
  * permanalyticse with r(459) at rep 0. the F2 b/V capture stays keyed on `marginsdof',
  * so fast mode still posts no e(V_es_*) (bug sweep BS-26)
  local sgmarginsdof `marginsdof'
  if "`fastnose'" != "" & "`subgroups'" != "" & `espre' + `espl' > 0 local sgmarginsdof df(`df_r')
  * 2026-09-13 author ruling: fast mode posts the event-study COEFFICIENT vector and
  * e(numpre_*) - sthlp:1295 has always said it does - while e(V_es_*) stays withheld
  * (BS-26 above). The F2 sites below therefore capture b alone wherever margins was asked
  * for nose, and the posting block reads an empty f2v_ as "no covariance to post".
  * permute()'s fast-mode default borrows `fastnose' from the shell (~5274) and is carved
  * out here: its stored-result contract is the permutation one, and its replications are
  * nested fastnose calls that consume e(b) alone. impvalues sets nose on its own account
  * below and is untouched (it is refused with fastnose in any case)
  local f2bonly = ("`fastnose'" != "" & "`permute'" == "")
  tempname f9tab

  * initialize the joint-test p-value locals
  local joint_horiz_i
  local joint_Choriz_i

  * do not calculate traditional SEs if producing Ibragimov and Muller style p-values
  if "`impvalues'" != "" {
    local marginsdof nose
    * --- IM (Ibragimov-Muller) degeneracy guard: shared tolerance + one-time caution flag ---
    * consumed by all six IM ttest sites. a collapsed cross-cluster spread (per-cluster estimates
    * numerically identical, e.g. noiseless/deterministic data) makes the IM SE not identified;
    * detect on the reported SE relative to the estimate (ttest returns r(se)/r(sd_1), NOT r(sd)),
    * strict relative only -- an absolute floor would false-fire on structural-zero estimands
    local imdegtol_rel = 1e-10
    * --- BS-10 (author ruling 2026-09-07): the IM reference distribution, for e(df_r) ---
    * every printed t, p and interval on this path comes from ttest, i.e. t(q-1) where q is
    * the number of clusters that contributed a non-missing per-cluster margin to THAT row,
    * and q varies row by row. The df that ereturn post's dof() carries is the stage-1
    * regression's residual df, which is always LARGER, so ereturn display / lincom / test /
    * esttab / coefplot reported p-values too small and intervals too narrow against the
    * table wooldid printed. Each surviving IM row appends its r(df_t) here; the MINIMUM is
    * written over e(df_r) after the post (conservative where the rows disagree). Posting it
    * missing was rejected in the ruling: the normal fallback is anti-conservative against
    * every row's t.
    local imdflist
    local imdegwarned
    local imdegmsg "Caution: for one or more Ibragimov-Muller estimates the per-cluster estimates were numerically identical, so the cross-cluster (IM) standard error is not identified; wooldid has set the standard error, test statistic, p-value, and 95% confidence interval to missing for the affected estimate(s) while retaining the point estimate. A collapsed IM standard error usually reflects deterministic or noiseless data - please verify your input if this was not expected."
  }


  * loop through requested estimate types - basically just different weights
  foreach item in `att' `att_it' `att_i' `att_itime' {
    if "`verbose'" == "verbose" dis as text "Wooldid is now computing `item' estimates for the full sample, including event study estimates."

    tempname ___`item'
    * use margins to produce the post treatment average effect within the treated group. captured via targetflag = 1 & t >= treatment period
    * do this for the non-impvalues case
    if "`impvalues'" == "" {
      qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' `marginsdof' `fcnoes'
      matrix `f9tab' = r(table)
      if rowsof(`f9tab') < 6 {
        matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
        matrix rownames `f9tab' = b se t pvalue ll ul
      }
      matrix `___`item'' = `f9tab'[1..6,2]'
    }
    * this block below produces IM p-values
    * basic principle: ignore normal SEs, get treatment effects per cluster. then, t-test the set of cluster-specific effects
    if "`impvalues'" == "impvalues" {
      *semielasticity injector here is constrained to dydx and eyex b/c x is binary
      qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' nose over(`cluster') `fcnoes'
      tempname ibmucoefs
      tempname ibmuframe
      * this column count / 2 object will appear frequently - it exists in the non-contreat margins commands b/c it is required to skip over a large slate of null estimates
      *     specifically, estimates for targetflag==0 have to be jumped over
      local matrixstartgrab = floor(`=colsof(r(b))'/2)+1
      matrix `ibmucoefs' =  r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
      frame create `ibmuframe'
      frame change `ibmuframe'
      qui svmat double `ibmucoefs'
      qui gstats sum `ibmucoefs'1, meanonly
      local b = r(mean)
      qui ttest `ibmucoefs'1 == 0
      local se = r(se)
      local imt = r(t)
      local imp = r(p)
      local im_ll = `b' + invt(r(df_t), (100 - `wdlvl')/200)*`se'
      local im_ul = `b' - invt(r(df_t), (100 - `wdlvl')/200)*`se'
      * BS-10 (ruling 2026-09-07): record this row's IM reference df. Only rows that SURVIVE
      * the degeneracy scrub immediately below count - the predicate here is the exact
      * negation of the scrub's, so the two can never disagree - and a row with a single
      * contributing cluster is excluded by !missing(`se') (ttest leaves r(se) missing and
      * r(df_t) = 0 there). r(df_t) is still live: the line above just read it, and nothing
      * between the ttest and here is r-class.
      if !missing(`se') & `se' >= `imdegtol_rel'*abs(`b') local imdflist `imdflist' `=r(df_t)'
      if !missing(`se') & `se' < `imdegtol_rel'*abs(`b') {
        local se .
        local imt .
        local imp .
        local im_ll .
        local im_ul .
        if "`imdegwarned'" == "" {
          dis as text "`imdegmsg'"
          local imdegwarned 1
        }
      }
      matrix `___`item'' = (`b',`se', `imt', `imp', `im_ll', `im_ul')
      frame change `defaultframe'
      frame drop `ibmuframe'
      matrix drop `ibmucoefs'
    }

    * average pre-treatment effect
    if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
      tempname ___`item'_pre

      if "`impvalues'" == "" {
        qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' `marginsdof' `fcnoes'
        matrix `f9tab' = r(table)
        if rowsof(`f9tab') < 6 {
          matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
          matrix rownames `f9tab' = b se t pvalue ll ul
        }
        matrix `___`item'_pre' = `f9tab'[1..6,2]'
      }

      * this produces an IM p-values
      if "`impvalues'" == "impvalues" {
        qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' nose over(`cluster') `fcnoes'
        tempname ibmucoefs
        tempname ibmuframe
        local matrixstartgrab floor(`=colsof(r(b))'/2)+1
        matrix `ibmucoefs' =  r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
        frame create `ibmuframe'
        frame change `ibmuframe'
        qui svmat double `ibmucoefs'
        qui gstats sum `ibmucoefs'1, meanonly
        local b = r(mean)
        qui ttest `ibmucoefs'1 == 0
        local se = r(se)
        local imt = r(t)
        local imp = r(p)
        local im_ll = `b' + invt(r(df_t), (100 - `wdlvl')/200)*`se'
        local im_ul = `b' - invt(r(df_t), (100 - `wdlvl')/200)*`se'
        * BS-10 (2026-09-07): this row's IM df for e(df_r); guard explained at the first IM site
        if !missing(`se') & `se' >= `imdegtol_rel'*abs(`b') local imdflist `imdflist' `=r(df_t)'
        if !missing(`se') & `se' < `imdegtol_rel'*abs(`b') {
          local se .
          local imt .
          local imp .
          local im_ll .
          local im_ul .
          if "`imdegwarned'" == "" {
            dis as text "`imdegmsg'"
            local imdegwarned 1
          }
        }
        matrix `___`item'_pre' = (`b',`se', `imt', `imp', `im_ll', `im_ul')
        frame change `defaultframe'
        frame drop `ibmuframe'
        matrix drop `ibmucoefs'
      }
      local mainpre 1
    }

    * effect by horizons
    if `espre' + `espl' > 0   {
      * att_i and att_itime are not particularly valuable concepts to bring to a horizons estimate
      * the reason: each i will only present once in a given horizon anyway
      * so the only way to get estimates where each i has an equal representation at each horizon
      * is to use att_it, which is how I address this
      local horizitem att
      if "`item'" != "att" local horizitem att_it
      cap confirm  matrix `___`horizitem'_horiz'
      * this prevents estimating horizons if done already for a non att - relevant b/c some weight types will always yield the same horizon specific estimate
      if _rc != 0 {
        * the idea behind this margins statement is to use the information in horizonstarter and horizonend to figure out which relative periods are requested
        * and then calculate all relative period effects in that range
        tempname ___`horizitem'_horiz

        if "`impvalues'" == "" {
          qui margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_t' >= `horizonstarter' & `_t' <= `horizonend' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof' `fcnoes'
          local matrixstartgrab = floor(`=colsof(r(b))'/2)+1
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `___`horizitem'_horiz' = ( `f9tab'[1..6, `matrixstartgrab'..`=colsof(r(b))']',`relyrcolumn')
          * F2: capture event-study b/V (relyr-labeled) for honestdid/lincom/test [pure-add]
          if "`marginsdof'" != "nose" | `f2bonly' {
            tempname f2bmat f2vmat
            matrix `f2bmat' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']
            local f2k `horizitem'
            local f2b_`f2k' `f2bmat'
            * nose leaves margins no r(V) to slice, so f2v_ stays empty and this family
            * posts its coefficients without a covariance
            if "`marginsdof'" != "nose" {
              matrix `f2vmat' = r(V)[`matrixstartgrab'..`=colsof(r(b))', `matrixstartgrab'..`=colsof(r(b))']
              local f2v_`f2k' `f2vmat'
            }
            local f2keys `f2keys' `f2k'
            local f2rl_`f2k'
            forval f2ri = 1/`=rowsof(`relyrcolumn')' {
              local f2rv = `relyrcolumn'[`f2ri',1]
              if `f2rv' < 0  local f2rl_`f2k' `f2rl_`f2k'' pre`=abs(`f2rv')'
              else          local f2rl_`f2k' `f2rl_`f2k'' post`f2rv'
            }
          }
          * NFADV-01/02/05: the adversarial statistics and display line live in
          * wooldid_advstat (worker size headroom); the finaladvtest estimability guard
          * stays here - mainly an issue w/ SGs later: sometimes we have enough coefs to
          * do the test overall, but not in a particular sg. since this is the overall
          * section, if it triggers badly here, call off adversarial as a whole
          *     threshold of 2 here b/c that corresponds w/ 1 coef (targetflag==0 will get dummy columns )
          local finaladvtest
          cap confirm matrix r(b)
          if _rc == 0 {
            local finaladvtest  `=colsof(r(b))'
          }
          if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
          if `finaladvtest' <= 2 {
            * this branch also runs when adversarial was never requested (or was already
            * cancelled upstream) - the clearing is a no-op then, so only Caution when
            * an actually-pending adversarial test is being called off here
            if "`adversarial'" == "adversarial" dis as text "Caution: adversarial pre-trend testing was cancelled: too few event-study coefficients were estimated to fit the reference line."
            local adversarial
          }
          if "`adversarial'" == "adversarial" & "`marginsdof'" != "nose"     {
            wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') fm(`wdfm')
            local advstatus `r(status)'
            * a gap-cancelled altmethod fit calls off adversarial as a whole (overall section)
            if "`advstatus'" == "gapskip" local adversarial
            else {
              if "`advstatus'" == "ok" local adversarialP_`horizitem' = r(p)
              tempname adversarialB_`horizitem'
              matrix `adversarialB_`horizitem'' = r(line)
            }
          }
          * end adversarial pre-trend handling
        }

        *constraints on margins setup requires looping through horizons for doing the impvalues approach. reason is you aren't allowed to do over (two sets of variables) without extra setup
        if "`impvalues'" == "impvalues" {
          local currenthoriz 1
          tempname horizhold
          foreach horizon in `horizlist' {
            qui  margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `_t' == `horizon' & `targetflag'==1  & `_t' != 0) noestimcheck `poismarginstext' nose over(`cluster') `fcnoes'
            tempname ibmucoefs
            tempname ibmuframe
            local matrixstartgrab floor(`=colsof(r(b))'/2)+1
            matrix `ibmucoefs' =  r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
            frame create `ibmuframe'
            frame change `ibmuframe'
            qui svmat double `ibmucoefs'
            qui gstats sum `ibmucoefs'1, meanonly
            local b = r(mean)
            qui ttest `ibmucoefs'1 == 0
            local se = r(se)
            local imt = r(t)
            local imp = r(p)
            local im_ll = `b' + invt(r(df_t), (100 - `wdlvl')/200)*`se'
            local im_ul = `b' - invt(r(df_t), (100 - `wdlvl')/200)*`se'
            * BS-10 (2026-09-07): this row's IM df for e(df_r); guard explained at the first IM site
            if !missing(`se') & `se' >= `imdegtol_rel'*abs(`b') local imdflist `imdflist' `=r(df_t)'
            if !missing(`se') & `se' < `imdegtol_rel'*abs(`b') {
              local se .
              local imt .
              local imp .
              local im_ll .
              local im_ul .
              if "`imdegwarned'" == "" {
                dis as text "`imdegmsg'"
                local imdegwarned 1
              }
            }
            if `currenthoriz' == 1 matrix `horizhold' = (`b',`se', `imt', `imp', `im_ll', `im_ul')
            if `currenthoriz' >1  matrix `horizhold' = `horizhold' \ (`b',`se', `imt', `imp', `im_ll', `im_ul')
            frame change `defaultframe'
            frame drop `ibmuframe'
            local ++currenthoriz
            matrix drop `ibmucoefs'
          }
          matrix `___`horizitem'_horiz' = ( `horizhold',`relyrcolumn')
        }
      }
    }
  }


  **** produce joint test of significance for pre-treatment horizons using the contrast option
  if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'`permute'"=="" {
    if "`verbose'" == "verbose" dis as text "Wooldid is now computing pre-treatment horizon effect joint tests."
    if "`att'" == "att" {
      * contrast(overjoint) joint-tests the hypothesis that each estimate by group specified in over() is the same
      * since we have specified our subpop to include _t == 0, which is the omitted case, this is a joint test vs 0
      qui margins [`weight'=`w_att'] ,`semielasticity'(`targetflag') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & (`_t' == 0 | `_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' `jointmarginsdof' `fcnoes'
      local joint_horiz = r(p)[1,2]
    }
    if "`joint_horiz_i'" == "" & (  "`att_i'" != "" | "`att_itime'" != "" | "`att_it'" != "" ){
      qui margins [`weight'=`w_att_it'],`semielasticity'(`targetflag') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . &  (`_t' == 0 | `_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' `jointmarginsdof' `fcnoes'
      local joint_horiz_i = r(p)[1,2]
    }
  }


  if "`verbose'" == "verbose" dis as text "Wooldid has completed production of all main average treatment effects for the full sample."









  **********************************************************************************************************************************************
  *****************************         MAKE CONTINUOUS TREATMENT EFFECTS PLOTS AND CALCULATE CONTINUOUS TREATMENT EFFECTS      ****************************
  **********************************************************************************************************************************************
  * graph window management: graphs are created with name(..., replace) so re-runs overwrite wooldid's own prior graphs; wooldid no longer clears the user's unrelated graphs from memory
  local cfxmatricestoreturn
  local cfxnp 0

  * continuous treatment estimates are constructed borderline identically to the main estimates until we get to cfx plots
  if "`contreat'" != "" {

    foreach item in `att' `att_it' `att_i' `att_itime' {
      if "`verbose'" == "verbose" dis as text "Wooldid is now working on `item'-weighted average marginal effect estimation and production of any requested continuous treatment effects plots, including plots by subgroup."
      tempname ___C`item'
      * posttreatment average marginal effect of the continuous treatment variable generated here - just a matter of specifying contreat vs targetflag to get an AME of some sort
      * we also no longer need the skip-half-the-matrix adjustment, because w/ a continuous variable things are not computed at 0 and 1
      if "`impvalues'" == "" {
        qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstextC' `marginsdof' `fcnoes'
        matrix `f9tab' = r(table)
        if rowsof(`f9tab') < 6 {
          matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
          matrix rownames `f9tab' = b se t pvalue ll ul
        }
        matrix `___C`item'' = `f9tab'[1..6,1]'
        * AMESE R1: under cell containment the CRVE ame variance is not the ame's
        if `amescrub' {
          forval amej = `amescrublo'/6 {
            matrix `___C`item''[1, `amej'] = .
          }
        }
      }

      if "`impvalues'" == "impvalues" {
        qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstextC' nose over(`cluster') `fcnoes'
        tempname ibmucoefs
        tempname ibmuframe
        matrix `ibmucoefs' =  r(b)[1, 1..`=colsof(r(b))']'
        frame create `ibmuframe'
        frame change `ibmuframe'
        qui svmat double `ibmucoefs'
        qui gstats sum `ibmucoefs'1, meanonly
        local b = r(mean)
        qui ttest `ibmucoefs'1 == 0
        local se = r(se)
        local imt = r(t)
        local imp = r(p)
        local im_ll = `b' + invt(r(df_t), (100 - `wdlvl')/200)*`se'
        local im_ul = `b' - invt(r(df_t), (100 - `wdlvl')/200)*`se'
        * BS-10 (2026-09-07): this row's IM df for e(df_r); guard explained at the first IM site
        if !missing(`se') & `se' >= `imdegtol_rel'*abs(`b') local imdflist `imdflist' `=r(df_t)'
        if !missing(`se') & `se' < `imdegtol_rel'*abs(`b') {
          local se .
          local imt .
          local imp .
          local im_ll .
          local im_ul .
          if "`imdegwarned'" == "" {
            dis as text "`imdegmsg'"
            local imdegwarned 1
          }
        }
        matrix `___C`item'' = (`b',`se', `imt', `imp', `im_ll', `im_ul')
        frame change `defaultframe'
        frame drop `ibmuframe'
        matrix drop `ibmucoefs'
      }

      * pretreatment average marginal effect of the continuous treatment variable
      if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
        tempname ___C`item'_pre

        if "`impvalues'" == "" {
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstextC' `marginsdof' `fcnoes'
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `___C`item'_pre' = `f9tab'[1..6,1]'
          * AMESE R1 scrub (see the main-ame twin)
          if `amescrub' {
            forval amej = `amescrublo'/6 {
              matrix `___C`item'_pre'[1, `amej'] = .
            }
          }
        }

        if "`impvalues'" == "impvalues" {
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstextC' nose over(`cluster') `fcnoes'
          tempname ibmucoefs
          tempname ibmuframe
          matrix `ibmucoefs' =  r(b)[1, 1..`=colsof(r(b))']'
          frame create `ibmuframe'
          frame change `ibmuframe'
          qui svmat double `ibmucoefs'
          qui gstats sum `ibmucoefs'1, meanonly
          local b = r(mean)
          qui ttest `ibmucoefs'1 == 0
          local se = r(se)
          local imt = r(t)
          local imp = r(p)
          local im_ll = `b' + invt(r(df_t), (100 - `wdlvl')/200)*`se'
          local im_ul = `b' - invt(r(df_t), (100 - `wdlvl')/200)*`se'
          * BS-10 (2026-09-07): this row's IM df for e(df_r); guard explained at the first IM site
          if !missing(`se') & `se' >= `imdegtol_rel'*abs(`b') local imdflist `imdflist' `=r(df_t)'
          if !missing(`se') & `se' < `imdegtol_rel'*abs(`b') {
            local se .
            local imt .
            local imp .
            local im_ll .
            local im_ul .
            if "`imdegwarned'" == "" {
              dis as text "`imdegmsg'"
              local imdegwarned 1
            }
          }
          matrix `___C`item'_pre' = (`b',`se', `imt', `imp', `im_ll', `im_ul')
          frame change `defaultframe'
          frame drop `ibmuframe'
          matrix drop `ibmucoefs'
        }
      }



      ************** MAKE CFX PLOTS - effects of treatment (ATT type and AME type) across the contreat distribution (across the lattice)
      if "`makecfxplots'" != "" & ("`item'" == "att" | "`safety'" == "off") {

        local cfxmatendpoint 6
        * D-1 (author ruling 2026-09-11): the permute() carve-out below applies to
        * cfxplotnose too. It used to truncate unconditionally, so under permute() every cfx
        * matrix came back one column wide, wooldid_permcellsx skipped all of them on its
        * six-column test, ndone stayed 0 - and so the drift caution never printed and
        * e(perm_cfx_blank) was never posted. The user got analytic plot cells under a
        * permutation header with nothing said. The 6-column shape is missing-padded here
        * (r(table) is short under nose and the padding below fills it), which is exactly
        * what the cell rewiring needs to fill p/ll/ul from the replications
        if "`cfxplotnose'" != "" & "`permute'" == "" local cfxmatendpoint 1
        * F9: fast mode keeps the cfx matrices in the established 1-column (point
        * estimates only) shape that cfxplotnose already produces, rather than minting
        * a new 6-column-with-missings shape
        * under permute() the 6-column shape is kept (missing-padded) so the shell's cell
        * rewiring can fill p/ll/ul from the replications - the same carve-out cfxplotnose
        * now takes above (D-1)
        if "`fastnose'" != "" & "`permute'" == "" local cfxmatendpoint 1
        * build the margins se/df option in a per-iteration copy rather than overwriting the user's
        * cfxplotnose flag: the overwrite leaked across the estimand loop, so every estimand after
        * the first was treated as if cfxplotnose had been specified, silently truncating its
        * returned cfx matrices to the estimate row only
        local cfxseopt `marginsdof'
        if "`cfxplotnose'" == "cfxplotnose" local cfxseopt nose

        * specify title for cfx plots x axis
        if "`cfxlattice'" != ""   local cfx_xtitle xtitle("estimation lattice: `cfxlattice'")
        if "`cfxlattice'" == "" & wordcount("`lattice'")==2 {
          local cfx_xtitle xtitle("`latticenum' quantiles (type: `latticetype')")
        }

        * use the lattice var for the plots if cfxlattice is not specified; a per-iteration copy
        * keeps the user's cfxlattice local pristine, so later iterations' x-axis titles do not
        * pick up the auto-lattice tempvar's internal name
        local cfxlatticeplot `cfxlattice'
        if "`cfxlatticeplot'" == "" local cfxlatticeplot `latticevar'

        * loop through plot types to make
        * note that we never apply impvalues here - except in rare cases, it probably (even more severely than normal) violates IM assumptions to do this at specific lattice points like this
        local cfxplotcaselist
        if `cfxpost' == 1 local cfxplotcaselist cfxpost
        if `cfxpre' == 1 local cfxplotcaselist `cfxplotcaselist'  cfxpre
        foreach cfxplotcase in `cfxplotcaselist'  `cfxrelyrlist' {
          local cfxif
          local nametext `cfxplotcase'
          if "`cfxplotcase'" == "cfxpost" & `cfxpost' ==1 {
            local cfxif & `t' >= `ttre' & `ttre' != .
            local titletext "Post-treatment `effecttype' Across the `contreat' distribution"
          }
          if "`cfxplotcase'" == "cfxpre" & `cfxpre' ==1 {
            local cfxif & `t' < `ttre' & `ttre' != .
            local titletext "Pre-treatment `effecttype' Across the `contreat' distribution"
          }
          if "`cfxplotcase'" != "cfxpre" & "`cfxplotcase'" != "cfxpost" {
            local cfxif & `_t' - `min_t' - 1 == `cfxplotcase' & `ttre' != .
            local nametext = "relyr" + subinstr("`cfxplotcase'","-","n",.)
            local titletext "`effecttype' Across the `contreat' distribution at Relative Time `cfxplotcase'"
          }

          if "`item'" == "att" local ytext "ATT on `y'"
          if "`item'" == "att_it" local ytext "ATT_it on `y'"
          if "`item'" == "att_itime" local ytext "ATT_itime on `y'"
          if "`item'" == "att_i" local ytext "ATT_i on `y'"

          * conceptually, all we are doing is asking for treatment effects to be generated using over(cfxlattice) - but this tends to be a pretty taxing request, actually
          * we then ask margins plot to handle plot production for us here
          qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstext' `fcnoes'
          tempname cfx`item'_`nametext'
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `cfx`item'_`nametext'' = `f9tab'[1..`cfxmatendpoint',floor(`=colsof(r(table))'/2)+1..`=colsof(r(table))']'
          * under bootstrap() the plot is drawn by wooldid_cfxplots from the REWIRED matrix
          * after the bootstrap (author ruling 2026-09-05), so neither the analytic plot nor
          * its file is produced here; the title and y-title it needs are kept per matrix
          if "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" {
            qui marginsplot, recastci(rarea) plotopts(color(navy)) ciopts(color(navy%50)) `cfx_xtitle'  scale(.85) ytitle("`ytext'") title("`titletext'") name("cfx`item'_`nametext'", replace)
          }
          if "`saveplots'" != "" & "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" qui graph save "`saveplots'_cfx`item'_`nametext'.gph", `replace'
          local cfxatt_ttl `"`titletext'"'
          local cfxatt_ytx `"`ytext'"'
          local cfxxt_meta `"`cfx_xtitle'"'

          * this code here just handles text to switch to contreat
          if "`item'" == "att" local ytext "AME on `y'"
          if "`item'" == "att_it" local ytext "AME_it on `y'"
          if "`item'" == "att_itime" local ytext "AME_itime on `y'"
          if "`item'" == "att_i" local ytext "AME_i on `y'"

          if "`cfxplotcase'" == "cfxpost" & `cfxpost' ==1 {
            local titletext "Post-treatment `Ceffecttype' Across the `contreat' distribution"
          }
          if "`cfxplotcase'" == "cfxpre" & `cfxpre' ==1 {
            local titletext "Pre-treatment `Ceffecttype' Across the `contreat' distribution"
          }
          if "`cfxplotcase'" != "cfxpre" & "`cfxplotcase'" != "cfxpost" {
            local titletext "`Ceffecttype' Across the `contreat' distribution at Relative Time `cfxplotcase'"
          }

          * proceed for continuous treatment
          local itemtext = subinstr("`item'","att","ame",.)
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstextC' `fcnoes'
          tempname cfx`itemtext'_`nametext'
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `cfx`itemtext'_`nametext'' = `f9tab'[1..`cfxmatendpoint',1..`=colsof(r(table))']'
          * AMESE R1: the ame cfx inference columns come from the same refused CRVE
          if `amescrub' & `cfxmatendpoint' > 1 {
            forval amer = 1/`=rowsof(`cfx`itemtext'_`nametext'')' {
              forval amej = `amescrublo'/`cfxmatendpoint' {
                matrix `cfx`itemtext'_`nametext''[`amer', `amej'] = .
              }
            }
          }
          local amecfxci recastci(rarea) ciopts(color(navy%50))
          if `amescrub' local amecfxci noci
          if "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" {
            qui marginsplot, `amecfxci' plotopts(color(navy)) `cfx_xtitle' ytitle("`ytext'")  scale(.85) title("`titletext'") name("cfx`itemtext'_`nametext'", replace)
          }
          if "`saveplots'" != "" & "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" qui graph save "`saveplots'_cfx`itemtext'_`nametext'.gph", `replace'
          local cfxmatricestoreturn `cfxmatricestoreturn' cfx`itemtext'_`nametext' cfx`item'_`nametext'
          * plot metadata by position in cfxmatricestoreturn (ame first, then att)
          local ++cfxnp
          local cfxttl`cfxnp' `"`titletext'"'
          local cfxytx`cfxnp' `"`ytext'"'
          local ++cfxnp
          local cfxttl`cfxnp' `"`cfxatt_ttl'"'
          local cfxytx`cfxnp' `"`cfxatt_ytx'"'
        }

        *********** MAKE CFX PLOTS BY SG (by looping through sgs)
        if "`makecfxplotsbysg'" == "makecfxplotsbysg" {
          foreach sgcase in `sglist' {
            local cfxplotcaselist
            if `cfxpost' == 1 local cfxplotcaselist cfxpost
            if `cfxpre' == 1 local cfxplotcaselist `cfxplotcaselist'  cfxpre
            foreach cfxplotcase in `cfxplotcaselist'  `cfxrelyrlist' {
              local cfxif
              local nametext `cfxplotcase'
              if "`cfxplotcase'" == "cfxpost" & `cfxpost' ==1 {
                local cfxif & `t' >= `ttre' & `ttre' != .
                local titletext "Post-treatment `effecttype' Across the `contreat' distribution for `subgroups' = `sgcase'"
              }
              if "`cfxplotcase'" == "cfxpre" & `cfxpre' ==1 {
                local cfxif & `t' < `ttre' & `ttre' != .
                local titletext "Pre-treatment `effecttype' Across the `contreat' distribution for `subgroups' = `sgcase'"
              }
              if "`cfxplotcase'" != "cfxpre" & "`cfxplotcase'" != "cfxpost" {
                local cfxif & `_t' - `min_t' - 1 == `cfxplotcase' & `ttre' != .
                local nametext = "relyr" + subinstr("`cfxplotcase'","-","n",.)
                local titletext "`effecttype' Across the `contreat' distribution at Relative Time `cfxplotcase' for `subgroups' = `sgcase'"
              }
              if "`item'" == "att" local ytext "ATT on `y'"
              if "`item'" == "att_it" local ytext "ATT_it on `y'"
              if "`item'" == "att_itime" local ytext "ATT_itime on `y'"
              if "`item'" == "att_i" local ytext "ATT_i on `y'"

              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `subgroups' == `sgcase' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstext' `fcnoes'
              tempname cfx`item'_`nametext'_sg`sgcase'
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `cfx`item'_`nametext'_sg`sgcase''  = `f9tab'[1..`cfxmatendpoint',floor(`=colsof(r(b))'/2)+1..`=colsof(r(b))']'
              if "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" {
                qui marginsplot, recastci(rarea) plotopts(color(navy)) ciopts(color(navy%50)) `cfx_xtitle' ytitle("`ytext'") title("`titletext'") scale(.85) name("cfx`item'_`nametext'_sg`sgcase'", replace)
              }
              if "`saveplots'" != "" & "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" qui graph save "`saveplots'_cfx`item'_`nametext'_sg`sgcase'.gph", `replace'
              local cfxatt_ttl `"`titletext'"'
              local cfxatt_ytx `"`ytext'"'
              if "`item'" == "att" local ytext "AME on `y'"
              if "`item'" == "att_it" local ytext "AME_it on `y'"
              if "`item'" == "att_itime" local ytext "AME_itime on `y'"
              if "`item'" == "att_i" local ytext "AME_i on `y'"

              if "`cfxplotcase'" == "cfxpost" & `cfxpost' ==1 {
                local titletext "Post-treatment `Ceffecttype' Across the `contreat' distribution for `subgroups' = `sgcase'"
              }
              if "`cfxplotcase'" == "cfxpre" & `cfxpre' ==1 {
                local titletext "Pre-treatment `Ceffecttype' Across the `contreat' distribution for `subgroups' = `sgcase'"
              }
              if "`cfxplotcase'" != "cfxpre" & "`cfxplotcase'" != "cfxpost" {
                local titletext "`Ceffecttype' Across the `contreat' distribution at Relative Time `cfxplotcase' for `subgroups' = `sgcase'"
              }

              local itemtext = subinstr("`item'","att","ame",.)
              qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `subgroups' == `sgcase' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstextC' `fcnoes'
              tempname cfx`itemtext'_`nametext'_sg`sgcase'
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `cfx`itemtext'_`nametext'_sg`sgcase'' = `f9tab'[1..`cfxmatendpoint',1..`=colsof(r(b))']'
              * AMESE R1 scrub (see the main cfx ame twin)
              if `amescrub' & `cfxmatendpoint' > 1 {
                forval amer = 1/`=rowsof(`cfx`itemtext'_`nametext'_sg`sgcase'')' {
                  forval amej = `amescrublo'/`cfxmatendpoint' {
                    matrix `cfx`itemtext'_`nametext'_sg`sgcase''[`amer', `amej'] = .
                  }
                }
              }
              local amecfxci recastci(rarea) ciopts(color(navy%50))
              if `amescrub' local amecfxci noci
              if "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" {
                qui marginsplot, plotopts(color(navy)) `amecfxci' `cfx_xtitle' ytitle("`ytext'") title("`titletext'")  scale(.85) name("cfx`itemtext'_`nametext'_sg`sgcase'", replace)
              }
              if "`saveplots'" != "" & "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" qui graph save "`saveplots'_cfx`itemtext'_`nametext'_sg`sgcase'.gph", `replace'
              local cfxmatricestoreturn `cfxmatricestoreturn'   cfx`itemtext'_`nametext'_sg`sgcase' cfx`item'_`nametext'_sg`sgcase'
              local ++cfxnp
              local cfxttl`cfxnp' `"`titletext'"'
              local cfxytx`cfxnp' `"`ytext'"'
              local ++cfxnp
              local cfxttl`cfxnp' `"`cfxatt_ttl'"'
              local cfxytx`cfxnp' `"`cfxatt_ytx'"'
            }
          }
        }
      }


      ********** average marginal effects by horizon  horizons
      if `espre' + `espl' > 0  {
        local horizitem att
        if "`item'" != "att" local horizitem att_it
        cap confirm  matrix `___C`horizitem'_horiz'
        * this prevents estimating horizons if done already for a non att
        if _rc != 0 {
          tempname ___C`horizitem'_horiz

          if "`impvalues'" == "" {
            qui  margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' ) noestimcheck   over(`_t') `poismarginstextC' `marginsdof' `fcnoes'
            matrix `f9tab' = r(table)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `___C`horizitem'_horiz' = (`f9tab'[1..6, 1..`=colsof(r(b))']',`relyrcolumn')
            * AMESE R1 scrub: null the se/t/p/ci columns of every ame ES horizon
            if `amescrub' {
              forval amer = 1/`=rowsof(`___C`horizitem'_horiz')' {
                forval amej = `amescrublo'/6 {
                  matrix `___C`horizitem'_horiz'[`amer', `amej'] = .
                }
              }
            }
            * F2: capture event-study b/V (relyr-labeled) for honestdid/lincom/test [pure-add]
            if "`marginsdof'" != "nose" | `f2bonly' {
              tempname f2bmat f2vmat
              matrix `f2bmat' = r(b)[1, 1..`=colsof(r(b))']
              local f2k = cond("`horizitem'"=="att","ame","ame_it")
              local f2b_`f2k' `f2bmat'
              if "`marginsdof'" != "nose" {
                matrix `f2vmat' = r(V)[1..`=colsof(r(b))', 1..`=colsof(r(b))']
                local f2v_`f2k' `f2vmat'
              }
              local f2keys `f2keys' `f2k'
              local f2rl_`f2k'
              forval f2ri = 1/`=rowsof(`relyrcolumn')' {
                local f2rv = `relyrcolumn'[`f2ri',1]
                if `f2rv' < 0  local f2rl_`f2k' `f2rl_`f2k'' pre`=abs(`f2rv')'
                else          local f2rl_`f2k' `f2rl_`f2k'' post`f2rv'
              }
            }
            * NFADV-01/02/05: adversarial statistics + display line via wooldid_advstat
            * (see the binary twin above; only the estimability guard stays in-worker)
            local finaladvtest
            cap confirm matrix r(b)
            if _rc == 0 {
              local finaladvtest  `=colsof(r(b))'
            }
            if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
            if `finaladvtest' <= 1 {
              * see the binary twin above: Caution only when a pending test is called off
              if "`adversarial'" == "adversarial" dis as text "Caution: adversarial pre-trend testing was cancelled: too few event-study coefficients were estimated to fit the reference line."
              local adversarial
            }
            if "`adversarial'" == "adversarial" & "`marginsdof'" != "nose" & `amescrub' == 0   {
              wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') fm(`wdfm')
              local advstatus `r(status)'
              if "`advstatus'" == "gapskip" local adversarial
              else {
                if "`advstatus'" == "ok" local adversarialP_C`horizitem' = r(p)
                tempname adversarialB_C`horizitem'
                matrix `adversarialB_C`horizitem'' = r(line)
              }
            }
            * end adversarial pre-trend handling
          }

          *constraints on margins setup requires looping through horizons
          if "`impvalues'" == "impvalues" {
            local currenthoriz 1
            tempname horizhold
            foreach horizon in `horizlist' {
              qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `_t' == `horizon' & `targetflag'==1  ) noestimcheck   over(`cluster') `poismarginstextC' nose `fcnoes'
              tempname ibmucoefs
              tempname ibmuframe
              matrix `ibmucoefs' =  r(b)[1, 1..`=colsof(r(b))']'
              frame create `ibmuframe'
              frame change `ibmuframe'
              qui svmat double `ibmucoefs'
              qui gstats sum `ibmucoefs'1, meanonly
              local b = r(mean)
              qui ttest `ibmucoefs'1 == 0
              local se = r(se)
              local imt = r(t)
              local imp = r(p)
              local im_ll = `b' + invt(r(df_t), (100 - `wdlvl')/200)*`se'
              local im_ul = `b' - invt(r(df_t), (100 - `wdlvl')/200)*`se'
              * BS-10 (2026-09-07): this row's IM df for e(df_r); guard explained at the first IM site
              if !missing(`se') & `se' >= `imdegtol_rel'*abs(`b') local imdflist `imdflist' `=r(df_t)'
              if !missing(`se') & `se' < `imdegtol_rel'*abs(`b') {
                local se .
                local imt .
                local imp .
                local im_ll .
                local im_ul .
                if "`imdegwarned'" == "" {
                  dis as text "`imdegmsg'"
                  local imdegwarned 1
                }
              }
              if `currenthoriz' == 1 matrix `horizhold' = (`b',`se', `imt', `imp', `im_ll', `im_ul')
              if `currenthoriz' >1  matrix `horizhold' = `horizhold' \ (`b',`se', `imt', `imp', `im_ll', `im_ul')
              frame change `defaultframe'
              frame drop `ibmuframe'
              local ++currenthoriz
              matrix drop `ibmucoefs'
            }
            matrix `___C`horizitem'_horiz' = ( `horizhold',`relyrcolumn')
          }
        }
      }


      ****** produce joint tests of pre-treatment AMEs by horizon - same method as for non continuous treatment
      * note - pret 1 case will suppress these already from att level processing
      * AMESE R3: under cell containment the joint pre-trend tests inherit the same
      * invalid ame covariance, so they are skipped (see the R1 caution)
      if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'`permute'"=="" & `amescrub' == 0 {
        if "`att'" == "att"  {
          qui margins [`weight'=`w_att'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != .  &  (`_t' == 0 |`_t' >= `horizonstarter')     ) contrast(overjoint) `poismarginstextC' `marginsdof' `fcnoes'
          local joint_Choriz = r(p)[1,1]
        }

        if "`joint_Choriz_i'" == "" & ( "`att_it'" != "" | "`att_itime'" != "" | "`att_i'" != "" ) {
          qui margins [`weight'=`w_att_it'],`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != .  &  (`_t' == 0 |`_t' >= `horizonstarter')  ) contrast(overjoint) `poismarginstextC' `marginsdof' `fcnoes'
          local joint_Choriz_i = r(p)[1,1]
        }
      }




    }
  }
  if "`verbose'" == "verbose" & "`contreat'" != "" dis as text "Wooldid has finished production of average continuous treatment effects plots and all full treatment sample average marginal effects."








  ***************************************************************************************************
  *****************************         MAKE MAIN RESULTS BY SG  ***********************************
  ***************************************************************************************************

  **************** repeat the production of main results by sg
  local sgpre
  if "`subgroups'" != "" {
    if "`verbose'" == "verbose" dis as text "Wooldid is now producing the final set of required margins calculations - those specific to particular subgroups."

    * the suppress command will lead us to only produce estimates comparing the sgs, not generate sg-specific estimates
    if "`suppresssgprimaryeffects'" == "" {
      foreach item in `att' `att_it' `att_i' `att_itime' {
        * posttreatment effect - same production method as before, but with an over(sg)  component
        tempname ___`item'_s
        if "`impvalues'" == "" {
          * EFF3: with more than one realized subgroup the sg-comparison pass below
          * issues this identical over(`subgroups') margins with pwcompare(effects)
          * bolted on, and Stata returns the subgroup LEVEL margins in r(table) beside
          * the pairwise rows in r(table_vs). defer the extraction to that pass rather
          * than paying for a second identical margins here. with a single realized
          * subgroup the comparison pass does not run at all, so the call stays.
          if wordcount("`sglist'") > 1 local sgeff3_`item' 1
          else {
            * BS-28 (author-ruled 2026-09-07): this dedicated call writes the SAME displayed
            * cells the screened pwcompare pass writes when more than one subgroup is realized,
            * so it takes the same F-collapse screen (see wooldid_fcdeg): a flagged call is
            * re-issued on the full data, and every later call runs there too (`fcsg' disarms,
            * `fcnoes' clears). `levels' is right here - nothing contrasts these rows, since the
            * pairwise and joint subgroup machinery lives on the pwcompare pass, which does not
            * run at one realized subgroup - so at n = 1 the predicate reduces to the omission
            * marker and a missing / non-positive level variance. The call is not posted, so no
            * un-post step is needed.
            forval fcpass = 1/2 {
              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups') `poismarginstext' `marginsdof' `fcnoes'
              if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') levels `verbose'
              if `fcsg' == 0 | `fcrerun' == 0 continue, break
              restore, preserve
              local fcnoes
              local fcsg 0
            }
            local matrixstartgrab floor(`=colsof(r(b))'/2)+1
            matrix `f9tab' = r(table)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `___`item'_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
          }
        }


        * pretreatment effects
        if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
          tempname ___`item'_pre_s
          local sgpre 1
          * F1 (2026-09-14): the subgroups the PRE-PERIOD subpopulation actually realizes.
          * when it realizes exactly one, margins drops the over() stripe altogether and names
          * the single row for the targetflag tempvar (1.__000005), so no returned name carries
          * the subgroup token and the name match below cannot say WHICH subgroup the row is.
          * that match used to rescue the shape off wordcount(`sglist'), the POST-period count,
          * which is a different quantity: at two post-period subgroups and one pre-period one
          * it never fired and every subgroup lost its pre-treatment row, estimable or not
          * domain = the subpop of the level call that produces the collapsed stripe (the
          * pwcompare twin adds & `ttre' != ., a subset, so this list is a superset of what
          * either call can realize and the rescue errs toward NOT firing). computed over the
          * data rather than e(sample) for the same reason. glevelsof is rclass: this site is
          * chosen because the posttreatment r() above has already been read into ___`item'_s
          * (or deferred to a later pass that re-issues its own margins), and nothing between
          * here and the pre-period margins call reads r()
          if "`sgprelistdone'" == "" {
            local sgprelist
            cap qui glevelsof(`subgroups') if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 , clean local(sgprelist)
            if _rc != 0 local sgprelist
            local sgprelistdone 1
          }

          if "`impvalues'" == "" {
            * EFF3: see the posttreatment twin. the pre-period comparison pass is
            * capture-wrapped (pwcompare hard-errors when the pre-period subpopulation
            * realizes a single subgroup), so its failure branch runs this call as a
            * fallback - the deferral must never be able to lose the level margins.
            if wordcount("`sglist'") > 1 local sgeff3pre_`item' 1
            else {
              * BS-28 (author-ruled 2026-09-07): screened like the posttreatment twin above.
              * `levels' additionally keeps the contrast test off a pre-period subpopulation
              * that can realize MORE than one subgroup even where `sglist' (glevelsof on
              * t >= ttre) holds one word - an unbalanced panel whose pre-period units have no
              * post rows - and nothing contrasts these level rows either way.
              forval fcpass = 1/2 {
                qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck over(`subgroups') `poismarginstext' `marginsdof' `fcnoes'
                if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') levels `verbose'
                if `fcsg' == 0 | `fcrerun' == 0 continue, break
                restore, preserve
                local fcnoes
                local fcsg 0
              }
              local matrixstartgrab floor(`=colsof(r(b))'/2)+1
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `___`item'_pre_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
            }
          }

        }

        * effects by horizon for sgs
        if `espre' + `espl' > 0 {
          foreach sgcase in `sglist' {
            local horizitem att
            if "`item'" != "att" local horizitem att_it
            cap confirm  matrix `___`horizitem'_horiz_s`sgcase''
            * this prevents estimating horizons if done already for a non att
            if _rc != 0 {
              * re: estimability by sg. given each sg is composed of cohorts, if the sg exists in sample, then
              *     then this margins should be estimable since we would have seen the sg's cohorts dropped earlier
              *     if there was no post period or an insufficient pre-period
              * problem: this doesn't guarantee we see the same relative periods, however
              * the per-estimand relyr list starts empty so a subgroup whose event study is
              * skipped entirely (below) names no ES rows for this estimand
              local relyrlist_s`sgcase'_`horizitem'
              * F-collapse screen (see wooldid_fcdeg): a flagged call is re-issued on the full
              * data, and every later call runs there too (`fcsg' disarms, `fcnoes' clears)
              forval fcpass = 1/2 {
                qui margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstext' `sgmarginsdof' `fcnoes'
                if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') levels `verbose'
                if `fcsg' == 0 | `fcrerun' == 0 continue, break
                restore, preserve
                local fcnoes
                local fcsg 0
              }
              local matrixstartgrab = floor(`=colsof(r(b))'/2)+1
              * capture the adversarial estimability count NOW, while r(b) is authoritative:
              * the single-relyr-level glevelsof fallback below is rclass and wipes r(b)
              local finaladvtest  `=colsof(r(b))'
              tempname ___`horizitem'_horiz_s`sgcase'
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              * F2: capture event-study b/V (relyr-labeled) for honestdid/lincom/test [pure-add]
              if "`marginsdof'" != "nose" | `f2bonly' {
                tempname f2bmat f2vmat
                matrix `f2bmat' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']
                local f2k `horizitem'_s`sgcase'
                local f2b_`f2k' `f2bmat'
                * sgmarginsdof keeps this ONE call's VCE in fast mode (BS-26, for the
                * estimability markers), so r(V) exists here where it does not elsewhere;
                * the grab stays keyed on `marginsdof' so fast mode still posts none
                if "`marginsdof'" != "nose" {
                  matrix `f2vmat' = r(V)[`matrixstartgrab'..`=colsof(r(b))', `matrixstartgrab'..`=colsof(r(b))']
                  local f2v_`f2k' `f2vmat'
                }
                local f2keys `f2keys' `f2k'
              }
              *the above will be updated w/ rownames / relyrs by the below

              * SGs need not have the same number of relative periods appear for each
              * this code block below figures out what relative periods appear for them
              tempname nameget relyrcolumn_s`sgcase'
              matrix `nameget' = r(table)[1, `matrixstartgrab'..`=colsof(r(b))']
              local names: colnames `nameget'
              matrix drop `nameget'
              * SI-7 (author-ruled 2026-08-27): margins marks a subgroup-by-horizon margin it
              * deems NOT ESTIMABLE with an omission marker (e.g. 5bno.<tvar>). Because this call
              * passes noestimcheck, margins computes such a margin anyway - so the marker is the
              * ONLY signal separating a real estimate from one produced under protest, and it must
              * be READ rather than stripped: stripping the o leaves a token that parses as a plain
              * number, and the non-estimable value is then reported as if it were an estimate
              * (the pre-fix behavior on builds where the token happened to parse). Instead drop
              * those horizons from THIS subgroup and say so; ragged subgroups are already
              * supported below ("not all sgs have all relyrs, necessarily"). The drop must be
              * positional: names, the estimates matrix and the F2 b/V all index the same margins
              * columns, and the matrix concatenation below requires them to stay conformable.
              * the scan, the drop and its caution, the relative-period labels and the positional
              * filter of the table and the F2 b/V live in wooldid_sglbl, shared with the
              * continuous-treatment twin below (moved out of this program on 2026-09-06 to
              * recover compiled-size headroom; the numbers are unchanged)
              local f2opt
              if "`marginsdof'" != "nose" local f2opt b(`f2bmat') v(`f2vmat')
              else if `f2bonly' local f2opt b(`f2bmat')
              wooldid_sglbl , names(`names') tab(`f9tab') first(`matrixstartgrab') mint(`min_t') old(`oldsyntax') sg(`sgcase') tvar(`_t') strip(1.`targetflag') fb(`ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' & `_t' != 0) rel(`relyrcolumn_s`sgcase'') lbl(relyrlist_s`sgcase'_`horizitem') f2l(f2rl_`horizitem'_s`sgcase') `f2opt'
              if `wdsg_skip' continue
              local sgndrop `wdsg_ndrop'
              local finaladvtest = `finaladvtest' - `sgndrop'
              local relyrlist_s`sgcase' `relyrlist_s`sgcase'_`horizitem''
              * BS-26 keeps this call's variance so the estimability markers exist in fast
              * mode; the variance itself is not part of fast mode's output, so discard it
              * here and leave the se/t/p/interval cells missing as before
              if "`fastnose'" != "" {
                forval fnj = 2/6 {
                  forval fni = 1/`=colsof(`f9tab')' {
                    matrix `f9tab'[`fnj', `fni'] = .
                  }
                }
              }
              matrix `___`horizitem'_horiz_s`sgcase'' = ( `f9tab'', `relyrcolumn_s`sgcase'')

              * we construct relyrcolumn_s`sgcase' because not all sgs have all relyrs, necessarily

              * handle the adversarial pre-trends test
              * NFADV-01/02/05: statistics + display line via wooldid_advstat (see the
              * overall binary site); a gap-cancelled fit only skips THIS subgroup
              * (finaladvtest was captured beside matrixstartgrab above: a late
              * cap-confirm read of r(b) here found glevelsof's r() instead whenever the
              * single-relyr-level fallback had fired - harmless only because that case
              * is <= 2 columns and fails the > 2 gate either way)
              if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
              if "`adversarial'" == "adversarial" & "`marginsdof'" != "nose"   {
                * SI-7 (author-ruled 2026-08-27): if any horizon was dropped as not estimable,
                * refuse the line test for THIS subgroup. The adversarial null is that the
                * coefficients lie on a single line in relative time; with a horizon missing the
                * fit is run on an incomplete - and possibly interior-holed - grid, so the
                * statistic would not mean what its name says. Other subgroups are unaffected.
                if `sgndrop' > 0 {
                  dis as text "Caution: the adversarial pre-trend test was skipped for subgroup `sgcase': `sgndrop' relative period effect(s) were dropped as not estimable, so its coefficients no longer span a complete relative-time grid and a line fit through them would not be interpretable."
                }
                else if  `finaladvtest' > 2  {
                  wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn_s`sgcase'') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') sg(`sgcase') fm(`wdfm')
                  local advstatus `r(status)'
                  if "`advstatus'" != "gapskip" {
                    if "`advstatus'" == "ok" local adversarialP_`horizitem'_s`sgcase' = r(p)
                    tempname adversarialB_`horizitem'_s`sgcase'
                    matrix `adversarialB_`horizitem'_s`sgcase'' = r(line)
                  }
                }
                else dis as text "Caution: the adversarial pre-trend test was skipped for subgroup `sgcase': too few event-study coefficients."
              }
            }
          }
        }
      }

      * produce joint tests of pre-treatment horizons by sg - accomplish by looping through sgs
      if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'`permute'"=="" {

        foreach sgcase in `sglist' {
          if "`att'" == "att"   {
             forval fcpass = 1/2 {
               qui margins  [`weight'=`w_att'] ,`semielasticity'(`targetflag') noestimcheck  subpop(if `ifstatement'  & `subgroups' == `sgcase' & `t' < `ttre' & `ttre' != .   &  (`_t' == 0 |`_t' >= `horizonstarter')    ) contrast(overjoint) over(`_t') `poismarginstext' `marginsdof' `fcnoes'
               if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') joint `verbose'
               if `fcsg' == 0 | `fcrerun' == 0 continue, break
               restore, preserve
               local fcnoes
               local fcsg 0
             }
             *this confirm matrix machinery handles the case where the joint test fails b/c a given sg has no estimable pre-treatment relative period effects
             cap confirm matrix r(p)
             if _rc == 0 {
               local joint_horiz_s`sgcase' = r(p)[1,2]
             }
             else{
               local joint_horiz_s`sgcase'
             }
          }
          if "`joint_horiz_i_s`sgcase''" == "" & ("`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" ){
            forval fcpass = 1/2 {
              qui margins  [`weight'=`w_att_it'] ,`semielasticity'(`targetflag') noestimcheck  subpop(if `ifstatement' & `subgroups' == `sgcase' & `t' < `ttre' & `ttre' != .   &  (`_t' == 0 |`_t' >= `horizonstarter')   ) contrast(overjoint) over(`_t') `poismarginstext' `marginsdof' `fcnoes'
              if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') joint `verbose'
              if `fcsg' == 0 | `fcrerun' == 0 continue, break
              restore, preserve
              local fcnoes
              local fcsg 0
            }
            cap confirm matrix r(p)
            if _rc == 0 {
              local joint_horiz_i_s`sgcase' = r(p)[1,2]
            }
            else{
              local joint_horiz_i_s`sgcase'
            }
          }
          * these parameters give you the 2nd half of the sg-specific joint test matrix, minus the joint test of all sgs
        }
      }

    }


    * joint test the sg effect estimates against each other - are they the same or different?
    * do not suppress these like normal joint tests; this is a key point of using SGs
    * in theory this could be expanded further as well
    * E1: ONE posted margins pass per (estimand x pre/post) now serves BOTH the pairwise
    * comparison table and the all-subgroups-equal joint test. the joint p comes from a
    * chain `test` on the posted margins, whose coefficient names are read back off the
    * posted e(b) (never built from subgroup values); contrast(overjoint) is not allowed
    * after a posted pwcompare (r322), so the separate contrast pass is gone. every post
    * is followed by est restore, and the pre-period pass is capture-wrapped because
    * pwcompare hard-errors (r111) where contrast(overjoint) merely no-ops when only one
    * subgroup level is realized in the pre-period subpopulation
    if wordcount("`sglist'") > 1 {
      * F9: under fast mode the pairwise table keeps its point estimates ("." SE cells)
      * but no joint test runs: the pass is not posted (a chain test after a nose post
      * would be r321), so the sgchain/test/restore machinery is skipped entirely
      local sgpost post
      if "`fastnose'" != "" local sgpost
      foreach item in `att' `att_it' `att_i' `att_itime' {
        * pairwise comparisons among the sgs AND the joint test that all sgs are the same
        * (F-collapse screen: a posted call is un-posted before the re-issue)
        forval fcpass = 1/2 {
          qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof' `sgpost' `fcnoes'
          if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') `verbose'
          if `fcsg' == 0 | `fcrerun' == 0 continue, break
          if "`sgpost'" != "" qui est restore `wdfm'
          restore, preserve
          local fcnoes
          local fcsg 0
        }
        tempname pwsgs_`item'
        matrix `f9tab' = r(table_vs)
        if rowsof(`f9tab') < 6 {
          matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
          matrix rownames `f9tab' = b se t pvalue ll ul
        }
        matrix `pwsgs_`item''  =  `f9tab'[1..6,floor(`=colsof(r(table_vs))'/2)+1..`=colsof(r(table_vs))']'
        * EFF3: the SAME pass returned the subgroup level margins in r(table); read
        * them here, before wooldid_sgchain (rclass) wipes r()
        if "`sgeff3_`item''" != "" {
          local matrixstartgrab floor(`=colsof(r(b))'/2)+1
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `___`item'_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
        }
        if "`sgpost'" != "" {
          wooldid_sgchain, half
          * rank guard (verifier e1_v1): on a rank-deficient posted VCE, -test- silently
          * drops redundant constraints and reports a residual-rank p where
          * contrast(overjoint) reported missing (a silent "." -> 0.0000 flip). when the
          * chain test is not provably full rank, defer to the old contrast engine so
          * the reported p is identical to the pre-E1 build by construction
          local sgfallback
          if `sgchainn' > 1 {
            qui test `sgchain'
            if r(df) == `sgchainn' - 1 local joint_`item'_sdiff = r(p)
            else {
              local joint_`item'_sdiff
              local sgfallback 1
            }
          }
          else local joint_`item'_sdiff
          qui est restore `wdfm'
          if "`sgfallback'" != "" {
            qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof' `fcnoes'
            local joint_`item'_sdiff = r(p)[1,2]
          }
        }
        else local joint_`item'_sdiff

        * repeat for pre-treatment
        if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
          forval fcpass = 1/2 {
            cap qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof' `sgpost' `fcnoes'
            local fcprerc = _rc
            if `fcsg' & `fcprerc' == 0 wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') `verbose'
            if `fcsg' == 0 | `fcprerc' != 0 | `fcrerun' == 0 continue, break
            if "`sgpost'" != "" qui est restore `wdfm'
            restore, preserve
            local fcnoes
            local fcsg 0
          }
          if `fcprerc' != 0 {
            * the pre-ame only exists for 1 sg (or margins otherwise refused): no joint
            * test and no pairwise rows, exactly as under the previous contrast-first order
            local prejoint_`item'_sdiff
            if "`sgpost'" != "" qui est restore `wdfm'
            * EFF3 fallback: the level margins were deferred to this pass, which just
            * refused. the plain over(`subgroups') call does NOT refuse here (it merely
            * returns the single realized subgroup), so run it now, after the restore
            if "`sgeff3pre_`item''" != "" {
              * BS-28 (author-ruled 2026-09-07): screened like the two twins above. This is a
              * LIVE site, not a no-op: the pwcompare loop just above breaks on `fcprerc' != 0
              * WITHOUT restoring, so `fcsg' can still be 1 here. `fcpass' is free to reuse -
              * that enclosing forval has already exited. Not posted (no `sgpost'), and the
              * estimates were restored on the line above, so no un-post step is needed.
              forval fcpass = 1/2 {
                qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck over(`subgroups') `poismarginstext' `marginsdof' `fcnoes'
                if `fcsg' wooldid_fcdeg , safety(`safetyinput') guard(`fcguard') levels `verbose'
                if `fcsg' == 0 | `fcrerun' == 0 continue, break
                restore, preserve
                local fcnoes
                local fcsg 0
              }
              local matrixstartgrab floor(`=colsof(r(b))'/2)+1
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `___`item'_pre_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
            }
          }
          else {
            tempname prepwsgs_`item'
            matrix `f9tab' = r(table_vs)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `prepwsgs_`item''  =  `f9tab'[1..6,floor(`=colsof(r(table_vs))'/2)+1..`=colsof(r(table_vs))']'
            * EFF3: level margins off this same pass (before wooldid_sgchain)
            if "`sgeff3pre_`item''" != "" {
              local matrixstartgrab floor(`=colsof(r(b))'/2)+1
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `___`item'_pre_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
            }
            * the pairwise rows are stacked into e(pwsgs) iff THIS flag is set - not iff
            * the joint p exists - so a chain-test failure after a successful pwcompare
            * cannot silently shorten e(pwsgs)
            local prepwok_`item' 1
            if "`sgpost'" != "" {
              wooldid_sgchain, half
              * rank guard: see the post-period twin above; a failed or not-full-rank
              * chain test defers to the old contrast engine (whose r(p) read keeps the
              * old cap semantics)
              local sgfallback
              if `sgchainn' > 1 {
                cap qui test `sgchain'
                if _rc == 0 & r(df) == `sgchainn' - 1 local prejoint_`item'_sdiff = r(p)
                else {
                  local prejoint_`item'_sdiff
                  local sgfallback 1
                }
              }
              else local prejoint_`item'_sdiff
              qui est restore `wdfm'
              if "`sgfallback'" != "" {
                qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof' `fcnoes'
                cap local prejoint_`item'_sdiff = r(p)[1,2]
                if _rc != 0 local prejoint_`item'_sdiff
              }
            }
            else local prejoint_`item'_sdiff
          }
        }
      }
    }


    * handle the continuous treatments by sg
    if "`contreat'" != "" {
      if "`suppresssgprimaryeffects'" == "" {
        foreach item in `att' `att_it' `att_i' `att_itime' {
          * posttreatment average marginal effect of the continuous treatment variable
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre'  ) noestimcheck over(`subgroups') `poismarginstextC' `marginsdof' `fcnoes'
          tempname ___C`item'_s
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `___C`item'_s' = `f9tab'[1..6,1..`=colsof(r(b))']'
          * AMESE R1 scrub (rows 2..6 = se/t/p/ci across all subgroup columns)
          if `amescrub' {
            forval amer = 1/`=rowsof(`___C`item'_s')' {
              forval amej = `amescrublo'/6 {
                matrix `___C`item'_s'[`amer', `amej'] = .
              }
            }
          }

          * pretreatment average marginal effect of the continuous treatment variable
          if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
            qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0  ) noestimcheck over(`subgroups') `poismarginstextC' `marginsdof' `fcnoes'
            tempname ___C`item'_pre_s
            matrix `f9tab' = r(table)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `___C`item'_pre_s' = `f9tab'[1..6,1..`=colsof(r(b))']'
            * AMESE R1 scrub (see the posttreatment twin)
            if `amescrub' {
              forval amer = 1/`=rowsof(`___C`item'_pre_s')' {
                forval amej = `amescrublo'/6 {
                  matrix `___C`item'_pre_s'[`amer', `amej'] = .
                }
              }
            }
          }

          * average marginal effects by horizon
          if `espre' + `espl' > 0  {
            local horizitem att
            if "`item'" != "att" local horizitem att_it

            foreach sgcase in `sglist' {
              cap confirm  matrix `___`horizitem'_Choriz_s`sgcase''
              * this prevents estimating horizons if done already for a non att
              if _rc != 0 {
                qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstextC' `marginsdof' `fcnoes'
                tempname ___`horizitem'_Choriz_s`sgcase'
                matrix `f9tab' = r(table)
                if rowsof(`f9tab') < 6 {
                  matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                  matrix rownames `f9tab' = b se t pvalue ll ul
                }
                * label this call from ITS OWN colnames (bug sweep BS-8): the binary call's relyr
                * column has already been SI-7-filtered on THAT call's omission markers, and margins
                * flags the dose slope independently of the level effect, so the two column sets can
                * differ - concatenating them was a conformability error r(503) on every SI-7 drop
                * under contreat(). r() is consumed here, before wooldid_sglbl (whose glevelsof
                * fallback is rclass) scans, labels and filters the call
                tempname namegetC relyrcolumn_Cs`sgcase'
                matrix `namegetC' = r(table)[1, 1..`=colsof(r(b))']
                local namesC: colnames `namegetC'
                matrix drop `namegetC'
                local finaladvtest = colsof(r(b))
                local f2kC = cond("`horizitem'" == "att", "ame", "ame_it")
                local f2optC
                if "`marginsdof'" != "nose" | `f2bonly' {
                  tempname f2bmat f2vmat
                  matrix `f2bmat' = r(b)[1, 1..`=colsof(r(b))']
                  local f2optC b(`f2bmat')
                  if "`marginsdof'" != "nose" {
                    matrix `f2vmat' = r(V)[1..`=colsof(r(b))', 1..`=colsof(r(b))']
                    local f2optC b(`f2bmat') v(`f2vmat')
                  }
                }
                wooldid_sglbl , names(`namesC') tab(`f9tab') mint(`min_t') old(`oldsyntax') sg(`sgcase') tvar(`_t') kind(continuous-treatment) fb(`ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' & `_t' != 0) rel(`relyrcolumn_Cs`sgcase'') lbl(relyrlist_s`sgcase'_`f2kC') f2l(f2rl_`f2kC'_s`sgcase') `f2optC'
                if `wdsg_skip' continue
                local finaladvtest = `finaladvtest' - `wdsg_ndrop'
                matrix `___`horizitem'_Choriz_s`sgcase'' = ( `f9tab'', `relyrcolumn_Cs`sgcase'')
                * AMESE R1 scrub: null se/t/p/ci columns of the per-sg ame ES rows
                if `amescrub' {
                  forval amer = 1/`=rowsof(`___`horizitem'_Choriz_s`sgcase'')' {
                    forval amej = `amescrublo'/6 {
                      matrix `___`horizitem'_Choriz_s`sgcase''[`amer', `amej'] = .
                    }
                  }
                }
                * F2: capture event-study b/V (relyr-labeled) for honestdid/lincom/test [pure-add]
                if "`marginsdof'" != "nose" | `f2bonly' {
                  local f2k `f2kC'_s`sgcase'
                  local f2b_`f2k' `f2bmat'
                  if "`marginsdof'" != "nose" local f2v_`f2k' `f2vmat'
                  local f2keys `f2keys' `f2k'
                }
                * handle the adversarial pre-trends test
                * NFADV-01/02/05: statistics + display line via wooldid_advstat (see the
                * overall binary site); a gap-cancelled fit only skips THIS subgroup
                if "`adversarial'" == "adversarial" & "`marginsdof'" != "nose" & `amescrub' == 0  {
                  if  `finaladvtest' > 1  {
                    wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn_Cs`sgcase'') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') sg(`sgcase') fm(`wdfm')
                    local advstatus `r(status)'
                    if "`advstatus'" != "gapskip" {
                      if "`advstatus'" == "ok" local adversarialP_C`horizitem'_s`sgcase' = r(p)
                      tempname adversarialB_C`horizitem'_s`sgcase'
                      matrix `adversarialB_C`horizitem'_s`sgcase'' = r(line)
                    }
                  }
                  else dis as text "Caution: the adversarial pre-trend test was skipped for subgroup `sgcase': too few event-study coefficients."
                }
              }
            }
          }
          * end item loop
        }

        * produce joint tests of pre-treatment AMES by horizon by sg
        * (AMESE R3: skipped under cell containment — same invalid covariance)
        if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'`permute'"=="" & `amescrub' == 0 {
          foreach sgcase in `sglist' {
            if "`att'" == "att" {
              qui margins  [`weight'=`w_att'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & `subgroups' == `sgcase'   &  (`_t' == 0 |`_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstextC' `marginsdof' `fcnoes'
              cap confirm matrix r(p)
              if _rc == 0 {
                * the continuous overjoint r(p) is 1x1 (no base column), so the real p
                * sits at [1,1] - unlike the binary twins' 1x2 geometry read at [1,2]
                local joint_Choriz_s`sgcase' = r(p)[1,1]
              }
              else{
                local joint_Choriz_s`sgcase'
              }
            }

            if ("`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" ) & "`joint_Choriz_i_s`sgcase''" == "" {
              qui margins [`weight'=`w_att_it'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & `subgroups' == `sgcase'  &  (`_t' == 0 |`_t' >= `horizonstarter')    ) contrast(overjoint) `poismarginstextC' `marginsdof' `fcnoes'
              cap confirm matrix r(p)
              if _rc == 0 {
                local joint_Choriz_i_s`sgcase' = r(p)[1,1]
              }
              else{
                local joint_Choriz_i_s`sgcase'
              }
            }
            * these parameters give you the 2nd half of the sg-specific joint test matrix, minus the joint test of all sgs
          }
        }

      }



      *************** joint test the sg AMEs  against each other -
      * E1: see the binary twin above - one posted pwcompare pass per (estimand x pre/post)
      * serves both the pairwise table and the all-subgroups-equal joint test
      if wordcount("`sglist'") > 1 {
        * F9: see the binary twin - no post, no joint test under fast mode
        * AMESE R3: same disablement under cell containment (the SG-equality joint
        * tests run on the invalid ame covariance); pairwise point estimates keep
        * flowing, with their se/t/p/ci cells scrubbed below
        local sgpost post
        if "`fastnose'" != "" local sgpost
        if `amescrub' == 1 local sgpost

        foreach item in `att' `att_it' `att_i' `att_itime' {

          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstextC' noestimcheck over(`subgroups') pwcompare(effects) `marginsdof' `sgpost' `fcnoes'
          tempname pwsgs_C`item'
          matrix `f9tab' = r(table_vs)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `pwsgs_C`item''  = `f9tab'[1..6,1..`=colsof(r(table_vs))']'
          * AMESE R1 scrub
          if `amescrub' {
            forval amer = 1/`=rowsof(`pwsgs_C`item'')' {
              forval amej = `amescrublo'/6 {
                matrix `pwsgs_C`item''[`amer', `amej'] = .
              }
            }
          }
          if "`sgpost'" != "" {
            wooldid_sgchain
            * rank guard: see the binary twin
            local sgfallback
            if `sgchainn' > 1 {
              qui test `sgchain'
              if r(df) == `sgchainn' - 1 local joint_C`item'_sdiff = r(p)
              else {
                local joint_C`item'_sdiff
                local sgfallback 1
              }
            }
            else local joint_C`item'_sdiff
            qui est restore `wdfm'
            if "`sgfallback'" != "" {
              qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstextC' noestimcheck over(`subgroups') contrast(overjoint) `marginsdof' `fcnoes'
              local joint_C`item'_sdiff = r(p)[1,1]
            }
          }
          else local joint_C`item'_sdiff

          * repeat for pre-treatment
          if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
            cap qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstextC' `marginsdof' `sgpost' `fcnoes'
            if _rc != 0 {
              local prejoint_C`item'_sdiff
              if "`sgpost'" != "" qui est restore `wdfm'
            }
            else {
              tempname prepwsgs_C`item'
              matrix `f9tab' = r(table_vs)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `prepwsgs_C`item''  =  `f9tab'[1..6,1..`=colsof(r(table_vs))']'
              * AMESE R1 scrub
              if `amescrub' {
                forval amer = 1/`=rowsof(`prepwsgs_C`item'')' {
                  forval amej = `amescrublo'/6 {
                    matrix `prepwsgs_C`item''[`amer', `amej'] = .
                  }
                }
              }
              * see the binary twin: stacking is gated on this flag, not the joint p
              local prepwok_C`item' 1
              if "`sgpost'" != "" {
                wooldid_sgchain
                * rank guard: see the binary twin
                local sgfallback
                if `sgchainn' > 1 {
                  cap qui test `sgchain'
                  if _rc == 0 & r(df) == `sgchainn' - 1 local prejoint_C`item'_sdiff = r(p)
                  else {
                    local prejoint_C`item'_sdiff
                    local sgfallback 1
                  }
                }
                else local prejoint_C`item'_sdiff
                qui est restore `wdfm'
                if "`sgfallback'" != "" {
                  qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstextC' `marginsdof' `fcnoes'
                  cap local prejoint_C`item'_sdiff = r(p)[1,1]
                  if _rc != 0 local prejoint_C`item'_sdiff
                }
              }
              else local prejoint_C`item'_sdiff
            }
          }
        }
      }


    }
  }

  * F-collapse: every margins call is behind us - put the full data back
  if "`fcon'" == "1" restore

  if "`verbose'" == "verbose" dis as text "Wooldid has completed estimation of all main effects and is now preparing to present them, as well as preparing to produce event study estimates."




  ***************************************************************************************************************************
  *****************************         CONSTRUCT RESULTS MATRIX TO BE RETURNED & CLEAN UP       ****************************
  ***************************************************************************************************************************


  * begin constructing results matrix
  * initialize column names
  local rownamespool

  * add the main estimates, pre-treatment estimates, and continuous treatment only average marginal effects (and pre-treatment analogues) + spacers
  local firsttimethrough 1
  local roundlist
  if   "`mainpre'" == "1" local roundlist pretreatment
  if "`contreat'" != "" & "`roundlist'" != "" local roundlist pretreatment cont contpre
  if "`contreat'" != "" & "`roundlist'" == "" local roundlist cont

  tempname wdidpooled

  * this processing is to make sure everything is stacked into one matrix and has the right row names and such
  foreach round in main `roundlist' {
    foreach object in `att'  `att_it' `att_i' `att_itime' {
      local processedobject `object'
      if "`round'" == "pretreatment" local processedobject `object'_pre
      if "`round'" == "cont" local processedobject C`object'
      if "`round'" == "contpre" local processedobject C`object'_pre
      if `firsttimethrough' == 1 matrix `wdidpooled' = ( `___`processedobject'',.)
      if `firsttimethrough' == 0 matrix `wdidpooled' = `wdidpooled' \  ( `___`processedobject'',.)
      matrix drop  `___`processedobject''
      local processedobject = subinstr("`processedobject'","Catt","ame",.)
      local processedobject = subinstr("`processedobject'","_pre","",.)
      local eqtext `round'
      if "`round'" == "cont" local eqtext main
      if "`round'" == "contpre" local eqtext pretreatment
      local rownamespool `rownamespool' `eqtext':`processedobject'
      * this rownamespool object will be built up successively across cases
      local firsttimethrough 0
    }
  }

  * if horizon estimates are to be produced, present them here
  if `espre' + `espl'  > 0 {
    if "`att'" == "att" {
      matrix `wdidpooled' = `wdidpooled' \ ( `___att_horiz')
      matrix drop `___att_horiz'
      foreach horiz in `horizlist' {
        local relyr = `horiz' - `min_t' - 1
        local arelyr = abs(`relyr')
        if `relyr' >= 0 local name "post`relyr'"
        if `relyr' == 0 & `oldsyntax' == 1 local name "contemp"
        if `relyr' < 0 local name "pre`arelyr'"
        local rownamespool `rownamespool' ES_att:`name'
      }
    }
    if "`att_i'" == "att_i"  |  "`att_it'" == "att_it" |  "`att_itime'" == "att_itime" {
      matrix `wdidpooled' = `wdidpooled' \ ( `___att_it_horiz')
      matrix drop `___att_it_horiz'
      foreach horiz in `horizlist' {
        local relyr = `horiz' - `min_t' - 1
        local arelyr = abs(`relyr')
        if `relyr' >= 0 local name "post`relyr'"
        if `relyr' == 0 & `oldsyntax' == 1 local name "contemp"
        if `relyr' < 0 local name "pre`arelyr'"
        local rownamespool `rownamespool' ES_att_it:`name'
      }
    }
    if "`contreat'" != "" {
      if "`att'" == "att" {
        matrix `wdidpooled' = `wdidpooled' \ ( `___Catt_horiz')
        matrix drop `___Catt_horiz'
        foreach horiz in `horizlist' {
          local relyr = `horiz' - `min_t' - 1
          local arelyr = abs(`relyr')
          if `relyr' >= 0 local name "post`relyr'"
          if `relyr' == 0 & `oldsyntax' == 1 local name "contemp"
          if `relyr' < 0 local name "pre`arelyr'"
          local rownamespool `rownamespool' ES_ame:`name'
        }
      }
      if "`att_i'" == "att_i"  |  "`att_it'" == "att_it" |  "`att_itime'" == "att_itime" {
        matrix `wdidpooled' = `wdidpooled' \ ( `___Catt_it_horiz')
        matrix drop `___Catt_it_horiz'
        foreach horiz in `horizlist' {
          local relyr = `horiz' - `min_t' - 1
          local arelyr = abs(`relyr')
          if `relyr' >= 0 local name "post`relyr'"
          if `relyr' == 0 & `oldsyntax' == 1 local name "contemp"
          if `relyr' < 0 local name "pre`arelyr'"
          local rownamespool `rownamespool' ES_ame_it:`name'
        }
      }
    }
  }

  matrix rownames `wdidpooled' = `rownamespool'
  if "`poisson'" == "" matrix colnames `wdidpooled' =  estimate se t p lb_`wdlvltxt'ci ub_`wdlvltxt'ci   relyr
  if "`poisson'" == "poisson" matrix colnames `wdidpooled' =  estimate se z p lb_`wdlvltxt'ci ub_`wdlvltxt'ci   relyr


  * with the success of the above, we should now repeat it by sg

  if "`subgroups'" != "" {
    if "`suppresssgprimaryeffects'" == "" {
      local sgnum 0
      local sgmatricestodrop
      foreach sgcase in `sglist' {
        *sgpre is for any sg; sghaspre is the sg being looped through. sghaspre yes is a default addressed later and updated only if necessary
        local sghaspre yes
        local ++sgnum
        * begin constructing results matrix
        * initialize column names
        * add the main estimates, pre-treatment estimates, and continuous treatment only average marginal effects (and pre-treatment analogues) + spacers
        local firsttimethrough 1
        local sgrowkeep
        local roundlist
        if  "`sgpre'" == "1" local roundlist pre
        if "`contreat'" != "" & "`roundlist'" != "" local roundlist pre cont contpre
        if "`contreat'" != "" & "`roundlist'" == "" local roundlist cont
        foreach round in main `roundlist' {
          foreach object in `att'  `att_it' `att_i' `att_itime' {
            local processedobject `object'_s
            if "`round'" == "pre" local processedobject `object'_pre_s
            if "`round'" == "cont" local processedobject C`object'_s
            if "`round'" == "contpre" local processedobject C`object'_pre_s

            * we require this since not all SGs have pre-treatment averages manufactured
            local truerow `sgnum'
            if ("`round'" == "pre" | "`round'" == "contpre")  {
              if "`sghaspre'" == "yes" {
                local sgsafetytest: rownames  `___`processedobject''
                local wordcounter 1
                local truerow
                foreach word in `sgsafetytest' {
                  local word = subinstr("`word'","bn.`subgroups'","",.)
                  local word = subinstr("`word'",".`subgroups'","",.)
                  if "`word'" == "`sgcase'" local truerow `wordcounter'
                  * include an adjustment in case only 1 sg realizes in the pre-period
                  * subpopulation: margins then returns the targetflag level rather than an
                  * over() stripe, so `word' carries no subgroup token at all. `sgprelist'
                  * (F1, built beside the pre-period margins call) is the one thing on hand
                  * that names the realizing subgroup, and the equality below is both tests at
                  * once - exactly one realized level, and it is this sgcase's
                  if strpos("`word'","`subgroups'") == 0 & "`sgprelist'" == "`sgcase'" local truerow `wordcounter'
                  local ++wordcounter
                }
                * F2 (2026-09-14): record the withholding rather than performing it silently.
                * the gate above skips every later round and object once sghaspre is no, so a
                * subgroup enters these three parallel lists exactly once. the reason is the
                * discriminator the site already holds: a subgroup missing from the realized
                * pre-period levels has no pre-period cells to average, while one that IS
                * there reached this line because margins returned its level carrying the
                * omission marker, which the strip above cannot turn back into the sg value
                * the round says WHICH rows go: flipping at pre costs the subgroup every
                * pre-treatment row it would have had (the gate suppresses contpre too),
                * while flipping at contpre costs only the continuous-treatment pre-treatment
                * marginal effect - its pre-treatment att row is already appended and stays.
                * (pre precedes contpre in roundlist, so a nocells subgroup always flips at
                * pre and the contpre x nocells pair is unreachable by construction)
                if "`truerow'" == "" {
                  local sghaspre no
                  local sgnopre `sgnopre' `sgcase'
                  if strpos(" `sgprelist' ", " `sgcase' ") == 0 local sgnoprewhy `sgnoprewhy' nocells
                  else local sgnoprewhy `sgnoprewhy' omitted
                  if "`round'" == "contpre" local sgnoprerow `sgnoprerow' ame
                  else local sgnoprerow `sgnoprerow' all
                }
              }
              else{
                local truerow
              }
            }

            if `firsttimethrough' == 1 tempname wdid_sg`sgcase'
            if `firsttimethrough' == 1 & ("`sghaspre'" == "yes" | "`round'" == "cont") matrix `wdid_sg`sgcase'' = ( `___`processedobject''[`truerow',1..6],.)
            if `firsttimethrough' == 0 & ("`sghaspre'" == "yes" | "`round'" == "cont") matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \  ( `___`processedobject''[`truerow',1..6],.)
            * record, in (round x object) order, whether this non-ES row was appended, so the
            * rowname filter below can name exactly the rows that exist (see the note there)
            if ("`sghaspre'" == "yes" | "`round'" == "cont") local sgrowkeep `sgrowkeep' 1
            else local sgrowkeep `sgrowkeep' 0
            if  "`truerow'" != "" local firsttimethrough 0
            local sgmatricestodrop `sgmatricestodrop' `___`processedobject''
          }
        }

        * if horizon estimates are to be produced, present them here
        if `espre' + `espl'  > 0 {
          * a subgroup whose event study was skipped entirely (no estimable relative
          * period) has no horizon matrix to append; its ES rownames are empty to match
          cap confirm matrix `___att_horiz_s`sgcase''
          if "`att'" == "att" & _rc == 0 {
            matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \ `___att_horiz_s`sgcase''
            matrix drop `___att_horiz_s`sgcase''
          }
          cap confirm matrix `___att_it_horiz_s`sgcase''
          if ("`att_i'" == "att_i"  |  "`att_it'" == "att_it" |  "`att_itime'" == "att_itime") & _rc == 0 {
            matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \ `___att_it_horiz_s`sgcase''
            matrix drop `___att_it_horiz_s`sgcase''
          }
          if "`contreat'" != "" {
            cap confirm matrix `___att_Choriz_s`sgcase''
            if "`att'" == "att" & _rc == 0 {
              matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \  `___att_Choriz_s`sgcase''
              matrix drop `___att_Choriz_s`sgcase''
            }
            cap confirm matrix `___att_it_Choriz_s`sgcase''
            if ("`att_i'" == "att_i"  |  "`att_it'" == "att_it" |  "`att_itime'" == "att_itime") & _rc == 0 {
              matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \  `___att_it_Choriz_s`sgcase''
              matrix drop `___att_it_Choriz_s`sgcase''
            }
          }
        }

        * ensure sg rownames are stripped of irrelevant es levels
        local rownamessg
        local sgnonesidx 0
        foreach rownametest in `rownamespool' {
          if strpos("`rownametest'","ES_")==0 {
            * the pooled non-ES rownames were generated by the same (round x object) loop that
            * fills `wdid_sg`sgcase'' above, in the same order, so the k-th pooled non-ES name is
            * the k-th sg non-ES row. keep it iff that row was actually appended: naming these
            * rows from the FINAL value of sghaspre instead lets an early pre round append a row
            * whose name a later pre round then suppresses, leaving the namelist shorter than the
            * matrix -- which stata silently pads by repeating the LAST name (r(0), no warning),
            * shifting every later row onto the label of the NEXT quantity and duplicating a colname
            local ++sgnonesidx
            local sgkeeprow : word `sgnonesidx' of `sgrowkeep'
            if "`sgkeeprow'" == "1" local rownamessg `rownamessg' `rownametest'
          }
          else{
            local testtext = subinstr("`rownametest'","ES_att_it:","",.)
            local testtext = subinstr("`testtext'","ES_ame_it:","",.)
            local testtext = subinstr("`testtext'","ES_ame:","",.)
            local testtext = subinstr("`testtext'","ES_att:","",.)
            local testtext = subinstr("`testtext'"," ","",.)
            * each ES block is named from ITS estimand's horizon list: att and att_it can
            * keep different horizon sets for a subgroup (SI-7 drops), and naming both
            * from one list either mislabels rows (stata pads a short namelist by
            * repeating its last name) or errors r(503). the contreat (ES_ame*) rows are
            * named from the AME call's own kept horizons (BS-8)
            local sgeslist `relyrlist_s`sgcase''
            if strpos("`rownametest'","ES_att:") == 1    local sgeslist `relyrlist_s`sgcase'_att'
            if strpos("`rownametest'","ES_att_it:") == 1 local sgeslist `relyrlist_s`sgcase'_att_it'
            if strpos("`rownametest'","ES_ame:") == 1    local sgeslist `relyrlist_s`sgcase'_ame'
            if strpos("`rownametest'","ES_ame_it:") == 1 local sgeslist `relyrlist_s`sgcase'_ame_it'
            if strpos("`sgeslist' ","`testtext' ") != 0 local rownamessg `rownamessg' `rownametest'
            * note: the spaces above are required to avoid triggering "1 is present" when presented with "11"
          }
        }
        matrix rownames `wdid_sg`sgcase'' = `rownamessg'
        local rownamessg`sgcase' `rownamessg'
        if "`poisson'" == "" matrix colnames `wdid_sg`sgcase'' =  estimate se t p lb_`wdlvltxt'ci ub_`wdlvltxt'ci   relyr
        if "`poisson'" == "poisson" matrix colnames `wdid_sg`sgcase'' =  estimate se z p lb_`wdlvltxt'ci ub_`wdlvltxt'ci   relyr


      }

      foreach matdrop in `sgmatricestodrop' {
        * do matrix clean up
        cap matrix drop `matdrop'
      }

    }

    ******* give column names to the pwsgs matrices; their default rownames are ok i think with eq name modifications
    ******* also, stack them into one
    if wordcount("`sglist'") > 1 {
      tempname pwsgs
      local itemround 0
      foreach item in `att'  `att_it' `att_i' `att_itime' {
        matrix rownames `pwsgs_`item'' = `item':
        if `itemround' == 0 matrix `pwsgs' = `pwsgs_`item''
        if `itemround' == 1 matrix `pwsgs' = `pwsgs' \ `pwsgs_`item''
        matrix drop `pwsgs_`item''
        local itemround 1
        if ("`prepwok_`item''" !=  "") & ( (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) ){
            matrix rownames `prepwsgs_`item'' = pre_`item':
            matrix `pwsgs' = `pwsgs' \ `prepwsgs_`item''
            matrix drop `prepwsgs_`item''
        }
      }
      if "`contreat'" != ""  {
        foreach item in `att'  `att_it' `att_i' `att_itime' {
          local itemtext = subinstr("`item'","att","ame",.)
          matrix rownames `pwsgs_C`item'' = `itemtext':
          matrix `pwsgs' = `pwsgs' \ `pwsgs_C`item''
          matrix drop `pwsgs_C`item''
          if ("`prepwok_C`item''" != "") & ( (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == ""))))) {
                matrix rownames `prepwsgs_C`item'' = pre_`itemtext':
                matrix `pwsgs' = `pwsgs' \ `prepwsgs_C`item''
                matrix drop `prepwsgs_C`item''
          }
        }
      }
      if "`poisson'" == "" matrix colnames `pwsgs' =  estimate se t p lb_`wdlvltxt'ci ub_`wdlvltxt'ci
      if "`poisson'" == "poisson" matrix colnames `pwsgs' =  estimate se z p lb_`wdlvltxt'ci ub_`wdlvltxt'ci

    }
  }







  ***************************************************************************************************************************
  ****************************         EVENT-STUDY COMPOSITION DIAGNOSTICS (F8)       **************************************
  ***************************************************************************************************************************
  * per reported horizon: distinct treatment-onset cohorts, distinct units of i (= the number of
  * unit-by-period coefficients averaged at that horizon), distinct clusters, estimation-sample
  * observations, and the horizon's share of the reported-window treated weight under both weight
  * concepts the event studies use (w_att = w x customweightmultiplier; w_att_it = that divided by
  * the i-by-t cell weight total, replicating the weight construction above - the tempvars
  * themselves are estimand-gated and w_per_it is already dropped, so the shares are recomputed
  * here so the matrix has fixed width regardless of which estimands were requested). rows come
  * from the SAME relyrcolumn / relyrcolumn_s<sg> sources as the reported ES matrices, so they
  * align one-for-one with the ES_* blocks (horizons without a surviving coefficient get no row).
  * composition is a descriptive counts-and-totals screen, not a test: flat columns rule out
  * changes in these REPORTED summaries (a cohort entering/leaving, a unit-count or weight-mass
  * shift); they do not rule out turnover of WHICH units contribute at unchanged counts, since
  * identities are not stored. moving columns do not by themselves establish composition bias.
  * this block runs NO margins call and NO e-class command, so the active e() is untouched and no
  * est restore is needed. it DOES clear r(); nothing between here and the ereturn post consumes
  * r(). gtools by() is hash-based: the user's sort order is not disturbed.
  local f8keys
  * skipped inside replication calls: this block only fills e(escomp*), which the
  * bootstrap and permutation engines never read (they consume coefficients), so it is
  * pure cost there - measured at 8-16% of rep-loop time on a realistic panel
  if `espre' + `espl' > 0 & "`wdrepmode'" != "1" {
    tempvar f8wa f8wpit f8wit
    qui gen double `f8wa' = `w' `customweighttext'
    qui gegen double `f8wpit' = total(`esamp' * `w'), by(`i' `t')
    qui gen double `f8wit' = `f8wa'/`f8wpit'
    qui replace `f8wit' = 0 if `f8wit' == .
    local f8cluster
    if "`cluster'" != "" & wordcount("`cluster'") == 1 local f8cluster `cluster'
    tempvar f8in f8tagc f8tagi f8tagcl f8tagt f8nc f8nu f8ncl f8n f8w f8wi
    qui gen byte `f8in' = (`esamp' == 1 & `targetflag' == 1 & `_t' >= `horizonstarter' & `_t' <= `horizonend')
    qui gegen byte   `f8tagc'  = tag(`_t' `ttre') if `f8in'
    qui gegen byte   `f8tagi'  = tag(`_t' `_i')   if `f8in'
    qui gegen byte   `f8tagt'  = tag(`_t')        if `f8in'
    qui gegen long   `f8nc'    = total(`f8tagc') if `f8in', by(`_t')
    qui gegen long   `f8nu'    = total(`f8tagi') if `f8in', by(`_t')
    if "`f8cluster'" != "" {
      qui gegen byte `f8tagcl' = tag(`_t' `f8cluster') if `f8in'
      qui gegen long `f8ncl'   = total(`f8tagcl') if `f8in', by(`_t')
    }
    else {
      qui gen byte `f8tagcl' = .
      qui gen double `f8ncl' = .
    }
    qui gegen long   `f8n'     = count(1)       if `f8in', by(`_t')
    qui gegen double `f8w'     = total(`f8wa')  if `f8in', by(`_t')
    qui gegen double `f8wi'    = total(`f8wit') if `f8in', by(`_t')
    tempname f8raw
    cap mkmat `_t' `f8nc' `f8nu' `f8ncl' `f8n' `f8w' `f8wi' if `f8tagt' == 1, matrix(`f8raw')
    local f8rawrc = _rc
    cap drop `f8in' `f8tagc' `f8tagi' `f8tagcl' `f8tagt' `f8nc' `f8nu' `f8ncl' `f8n' `f8w' `f8wi'
    if `f8rawrc' == 0 {
      tempname f8m
      local f8nr = rowsof(`relyrcolumn')
      matrix `f8m' = J(`f8nr', 8, 0)
      local f8wtot = 0
      local f8witot = 0
      local f8rn
      forval f8i = 1/`f8nr' {
        local f8rv = `relyrcolumn'[`f8i',1]
        local f8h  = `f8rv' + `min_t' + 1
        if `f8rv' < 0 local f8rn `f8rn' pre`=abs(`f8rv')'
        else if `f8rv' == 0 & `oldsyntax' == 1 local f8rn `f8rn' contemp
        else local f8rn `f8rn' post`f8rv'
        matrix `f8m'[`f8i',8] = `f8rv'
        forval f8r = 1/`=rowsof(`f8raw')' {
          if `f8raw'[`f8r',1] == `f8h' {
            matrix `f8m'[`f8i',1] = `f8raw'[`f8r',2]
            matrix `f8m'[`f8i',2] = `f8raw'[`f8r',3]
            matrix `f8m'[`f8i',3] = `f8raw'[`f8r',4]
            matrix `f8m'[`f8i',4] = `f8raw'[`f8r',5]
            matrix `f8m'[`f8i',5] = `f8raw'[`f8r',6]
            matrix `f8m'[`f8i',6] = `f8raw'[`f8r',6]
            matrix `f8m'[`f8i',7] = `f8raw'[`f8r',7]
            local f8wtot = `f8wtot' + `f8raw'[`f8r',6]
            local f8witot = `f8witot' + `f8raw'[`f8r',7]
          }
        }
      }
      * share denominators accumulate over the matched (reported) horizons only, so each
      * share column sums to 1 by construction
      forval f8i = 1/`f8nr' {
        if `f8wtot' > 0 & `f8wtot' < .   matrix `f8m'[`f8i',6] = `f8m'[`f8i',6]/`f8wtot'
        if `f8witot' > 0 & `f8witot' < . matrix `f8m'[`f8i',7] = `f8m'[`f8i',7]/`f8witot'
      }
      matrix rownames `f8m' = `f8rn'
      matrix colnames `f8m' = ncohorts nunits nclusters nobs sumw_att wshare_att wshare_att_it relyr
      local f8mat_main `f8m'
      local f8keys `f8keys' main
      * reference-pool and whole-treated-sample coverage scalars (the reported window need not
      * exhaust the treated sample; under esfixedbaseperiod no reference stats exist elsewhere)
      tempvar f8refc f8refu
      qui count if `esamp' == 1 & `_t' == 0 & `ttre' != .
      local f8refN = r(N)
      qui gegen byte `f8refc' = tag(`ttre') if `esamp' == 1 & `_t' == 0 & `ttre' != .
      qui gegen byte `f8refu' = tag(`i')    if `esamp' == 1 & `_t' == 0 & `ttre' != .
      qui count if `f8refc' == 1
      local f8refncoh = r(N)
      qui count if `f8refu' == 1
      local f8refnunits = r(N)
      qui drop `f8refc' `f8refu'
      qui sum `f8wa' if `esamp' == 1 & `_t' == 0 & `ttre' != ., meanonly
      local f8refW = cond(r(N) > 0, r(sum), 0)
      qui count if `esamp' == 1 & `ttre' != .
      local f8trtN = r(N)
      qui sum `f8wa' if `esamp' == 1 & `ttre' != ., meanonly
      local f8trtW = cond(r(N) > 0, r(sum), 0)
    }
    cap matrix drop `f8raw'

    * ---- subgroup scopes: one extra pass keyed on (_t, sg) covers every subgroup at once ----
    if "`subgroups'" != "" & "`suppresssgprimaryeffects'" == "" {
      tempvar f8sin f8stagc f8stagi f8stagcl f8stagt f8snc f8snu f8sncl f8sn f8sw f8swi
      qui gen byte `f8sin' = (`esamp' == 1 & `targetflag' == 1 & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' < .)
      qui gegen byte   `f8stagc'  = tag(`_t' `subgroups' `ttre') if `f8sin'
      qui gegen byte   `f8stagi'  = tag(`_t' `subgroups' `_i')   if `f8sin'
      qui gegen byte   `f8stagt'  = tag(`_t' `subgroups')        if `f8sin'
      qui gegen long   `f8snc'    = total(`f8stagc') if `f8sin', by(`_t' `subgroups')
      qui gegen long   `f8snu'    = total(`f8stagi') if `f8sin', by(`_t' `subgroups')
      if "`f8cluster'" != "" {
        qui gegen byte `f8stagcl' = tag(`_t' `subgroups' `f8cluster') if `f8sin'
        qui gegen long `f8sncl'   = total(`f8stagcl') if `f8sin', by(`_t' `subgroups')
      }
      else {
        qui gen byte `f8stagcl' = .
        qui gen double `f8sncl' = .
      }
      qui gegen long   `f8sn'     = count(1)        if `f8sin', by(`_t' `subgroups')
      qui gegen double `f8sw'     = total(`f8wa')   if `f8sin', by(`_t' `subgroups')
      qui gegen double `f8swi'    = total(`f8wit')  if `f8sin', by(`_t' `subgroups')
      tempname f8sraw
      cap mkmat `_t' `subgroups' `f8snc' `f8snu' `f8sncl' `f8sn' `f8sw' `f8swi' if `f8stagt' == 1, matrix(`f8sraw')
      local f8srawrc = _rc
      cap drop `f8sin' `f8stagc' `f8stagi' `f8stagcl' `f8stagt' `f8snc' `f8snu' `f8sncl' `f8sn' `f8sw' `f8swi'
      if `f8srawrc' == 0 {
        foreach sgcase in `sglist' {
          cap confirm matrix `relyrcolumn_s`sgcase''
          if _rc == 0 {
            tempname f8sm
            local f8snr = rowsof(`relyrcolumn_s`sgcase'')
            matrix `f8sm' = J(`f8snr', 8, 0)
            local f8swtot = 0
            local f8switot = 0
            local f8srn
            forval f8i = 1/`f8snr' {
              local f8rv = `relyrcolumn_s`sgcase''[`f8i',1]
              local f8h  = `f8rv' + `min_t' + 1
              if `f8rv' < 0 local f8srn `f8srn' pre`=abs(`f8rv')'
              else if `f8rv' == 0 & `oldsyntax' == 1 local f8srn `f8srn' contemp
              else local f8srn `f8srn' post`f8rv'
              matrix `f8sm'[`f8i',8] = `f8rv'
              forval f8r = 1/`=rowsof(`f8sraw')' {
                if `f8sraw'[`f8r',1] == `f8h' & `f8sraw'[`f8r',2] == `sgcase' {
                  matrix `f8sm'[`f8i',1] = `f8sraw'[`f8r',3]
                  matrix `f8sm'[`f8i',2] = `f8sraw'[`f8r',4]
                  matrix `f8sm'[`f8i',3] = `f8sraw'[`f8r',5]
                  matrix `f8sm'[`f8i',4] = `f8sraw'[`f8r',6]
                  matrix `f8sm'[`f8i',5] = `f8sraw'[`f8r',7]
                  matrix `f8sm'[`f8i',6] = `f8sraw'[`f8r',7]
                  matrix `f8sm'[`f8i',7] = `f8sraw'[`f8r',8]
                  local f8swtot = `f8swtot' + `f8sraw'[`f8r',7]
                  local f8switot = `f8switot' + `f8sraw'[`f8r',8]
                }
              }
            }
            forval f8i = 1/`f8snr' {
              if `f8swtot' > 0 & `f8swtot' < .   matrix `f8sm'[`f8i',6] = `f8sm'[`f8i',6]/`f8swtot'
              if `f8switot' > 0 & `f8switot' < . matrix `f8sm'[`f8i',7] = `f8sm'[`f8i',7]/`f8switot'
            }
            matrix rownames `f8sm' = `f8srn'
            matrix colnames `f8sm' = ncohorts nunits nclusters nobs sumw_att wshare_att wshare_att_it relyr
            local f8mat_s`sgcase' `f8sm'
            local f8keys `f8keys' s`sgcase'
          }
        }
      }
      cap matrix drop `f8sraw'
    }
    cap drop `f8wa' `f8wpit' `f8wit'
  }

  ***************************************************************************************************************************
  ***************************************         PRESENT AND RETURN THE RESULTS       **************************************
  ***************************************************************************************************************************

  ******************* In the first block, we signal that results are coming,
  * list the sample N and number of clusters
  * present diagnostics on sample size in the pre-treatment group if not using esfixed base
  * present some R^2 measures
  * and then present some basic summary statistics

  * result presentation must start with posting items to e - allows saving as .ster and using esttab
  * this also clears e
  tempname eb
  tempname ev
  tempname ebformer
  matrix `ebformer' = `wdidpooled'
  local rownamespool = subinstr("`rownamespool'",":","_",.)

  * the plan is just to push through result estimates and SEs to e(b) and e(V)
  * without giving any information about covariance (which we did not estimate anyway)
  * this means we will get the correct coef and SE when using esttab... but you won't be able to test estimates against each other, use lincom to combine estimates, etc.
  if "`subgroups'" != "" {
    if "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        matrix `ebformer' = `ebformer' \ `wdid_sg`sgcase''
        local rownamespool_sginj = subinstr("`rownamessg`sgcase''",":","_SG`sgcase'_",.)
        local rownamespool `rownamespool' `rownamespool_sginj'
      }
    }
  }

  * SI-A4/R6: ereturn post cannot accept missing coefficients, so the voided ame
  * entries re-enter e(b) as the regression's structural zeros (that one raw surface
  * is out of reach - the named caution says so); every wooldid-native surface
  * (display, wdidmainresults, sg matrices, plots) keeps them missing
  if `amevoid' {
    forval r6r = 1/`=rowsof(`ebformer')' {
      if `ebformer'[`r6r', 1] >= . matrix `ebformer'[`r6r', 1] = 0
    }
  }

  * AMESE generic guard: identical-SE collapse signature scan (caution only; the
  * structural predicate above is the enforcement layer, and already-scrubbed runs
  * have nothing left to scan)
  local amecvw2 = word(`"`customvce'"', 1)
  if ("`cluster'" != "" | (strlen("`amecvw2'") >= 2 & strpos("cluster", "`amecvw2'") == 1)) & `amescrub' == 0 & "`fastnose'" == "" {
    wooldid_amescreen , mode(sescan) semat(`ebformer')
  }

  matrix `eb' = `ebformer'[1..`=rowsof(`ebformer')',1]'
  matrix `ev' = diag(hadamard(`ebformer'[1..`=rowsof(`ebformer')',2],`ebformer'[1..`=rowsof(`ebformer')',2]))
  matrix rownames `ev' = :
  matrix rownames `eb' = :
  matrix colnames `ev' = :
  matrix colnames `eb' = :
  matrix rownames `ev' = `rownamespool'
  matrix colnames `ev' = `rownamespool'
  matrix colnames `eb' = `rownamespool'
  matrix rownames `eb' = `y'

  * the return statement also passes N and dof count and the function flagging the estimation sample
  cap ereturn post `eb' `ev', esample(`esamp') depname(`y') obs(`N') dof(`df_r')
  local evpostrc = _rc
  if `evpostrc' != 0 {
    * partial-identification rescue: ereturn post rejects a VCE carrying ANY missing diagonal cell,
    * so one estimand with an unidentified SE (an IM ttest with a single contributing cluster, or a
    * row nulled by the IM degeneracy guard) would otherwise cost every OTHER estimand its posted
    * variance. When at least one estimand IS identified, replace the missing diagonal cells with a
    * negative sentinel: Stata accepts it and every display/inference surface renders it as a
    * missing standard error (ereturn display/esttab show ., test drops the constraint, lincom
    * returns ., sqrt() of it is missing), so nothing is fabricated. Runs with no identified
    * estimand at all keep the b-only fallback below.
    local evpostgood 0
    local evpostmiss 0
    * capture the count now: a successful ereturn post consumes `eb'/`ev', so their
    * dimensions are no longer available when the caution below is displayed
    local evposttotal = colsof(`ev')
    forval evpostd = 1/`evposttotal' {
      if missing(`ev'[`evpostd',`evpostd']) local ++evpostmiss
      else if `ev'[`evpostd',`evpostd'] > 0 & `ev'[`evpostd',`evpostd'] < . local ++evpostgood
    }
    if `evpostgood' > 0 & `evpostmiss' > 0 {
      forval evpostd = 1/`=colsof(`ev')' {
        if missing(`ev'[`evpostd',`evpostd']) matrix `ev'[`evpostd',`evpostd'] = -1e100
      }
      cap ereturn post `eb' `ev', esample(`esamp') depname(`y') obs(`N') dof(`df_r')
      local evpostrc = _rc
      if `evpostrc' == 0 {
        dis as text "Caution: standard errors were not identified for `evpostmiss' of `evposttotal' reported estimates. e(V) contains placeholder entries for the affected estimates that display as missing standard errors; all other estimates retain their variances in e(V). See e(wdidmainresults) and the displayed tables for per-estimate detail."
      }
    }
    if `evpostrc' != 0 {
      cap ereturn post `eb' , esample(`esamp') depname(`y') obs(`N') dof(`df_r')
      if _rc != 0 ereturn post   , esample(`esamp') depname(`y') obs(`N') dof(`df_r')
      if `evpostmiss' > 0 & "`fastnose'" == "" dis as text "Caution: standard errors were not identified for `evpostmiss' of `evposttotal' reported estimates, so no variance-covariance matrix e(V) was posted. Point estimates remain in e(b); per-estimate results remain in e(wdidmainresults) and the displayed tables."
      if `evpostmiss' == 0 dis as text "Caution: wooldid could not post a variance-covariance matrix e(V) for this run. Per-estimate results remain in e(wdidmainresults) and the displayed tables."
    }
  }


  * BS-10 (author ruling 2026-09-07): under impvalues the reference distribution the table
  * actually used is IM's row-specific t(q-1), not the stage-1 residual df that dof() just
  * posted. Overwrite e(df_r) with the MINIMUM IM df across the reported rows so a hand-run
  * ereturn display / lincom / test / esttab is conservative against the printed table
  * instead of narrower than it, and keep the analytic value as the hidden twin - exactly
  * the convention wooldid_bootwire uses for the bootstrap reference distribution. Never
  * missing: the normal fallback is anti-conservative against every row's t. The reduction
  * and the ereturns live in a subprogram (wooldid_worker is at the compiled-size ceiling);
  * an eclass subprogram's ereturns persist into the caller's e(), the way wooldid_ptp's do.
  * NOTE: the e(gteffects) cells are stage-1 regression coefficients, computed before this
  * point on `df_r', so they remain on e(df_r_analytic) - see the help sentence added here.
  if "`impvalues'" != "" & "`imdflist'" != "" wooldid_imdf , dfl(`imdflist') ana(`df_r')

  * return the command line input and estimate_type
  * F2: expose captured event-study b/V (honestdid/lincom/test compatible) [pure-add, after post]
  foreach f2ky of local f2keys {
    local f2bm `f2b_`f2ky''
    local f2vm `f2v_`f2ky''
    local f2rl `f2rl_`f2ky''
    if "`f2bm'" != "" & "`f2rl'" != "" {
      if `=colsof(`f2bm')' == `:list sizeof f2rl' & `=colsof(`f2bm')' > 0 {
        cap mata: st_matrix("`f2bm'", st_matrix("`f2bm'"))
        matrix colnames `f2bm' = `f2rl'
        * fast mode captures b alone (`f2bonly' at the nose switch), so there is nothing
        * to stripe here and nothing to score below: f2ok stays 0 and the family posts its
        * coefficients under the Caution, which is what sthlp:1295 describes
        if "`f2vm'" != "" {
          cap mata: st_matrix("`f2vm'", st_matrix("`f2vm'"))
          matrix colnames `f2vm' = `f2rl'
          matrix rownames `f2vm' = `f2rl'
        }
        local f2np 0
        foreach f2t of local f2rl {
          if substr("`f2t'",1,3) == "pre" local ++f2np
        }
        * skip a fully non-identified VCV (e.g. timetrends): no positive-finite diagonal.
        * NFES-02: also count the degenerate rows (zero/missing/negative diagonals) -
        * partially degenerate families still POST (adjudicated: structural zeros are
        * genuine margins output), but silently shipping them under a joint-inference-
        * advertised name was the gap: caution + machine-readable marker below
        local f2ok 0
        local f2nbad 0
        local f2k 0
        if "`f2vm'" != "" local f2k = colsof(`f2vm')
        forval f2d = 1/`f2k' {
          if `f2vm'[`f2d',`f2d'] > 0 & `f2vm'[`f2d',`f2d'] < . local f2ok 1
          else local ++f2nbad
        }
        * 2026-09-08 author ruling: the POINT ESTIMATES post even when no diagonal is
        * usable. A family whose covariance is entirely non-identified - the all-zero-VCE
        * case that the removed explicit-dummies refit used to paper over, or timetrends -
        * still has meaningful estimates, and withholding them took the honestdid and
        * pretrendpower inputs down with the inference. The covariance still posts only
        * when some diagonal entry is positive and finite.
        * SI-A4/R6: under a voided ame block the captured points are structural
        * zeros, so b_es_ame* is withheld along with V_es_ame* (do not trust
        * f2ok=0 to do this - it was observed on the audited d3 runs, not
        * guaranteed across the class)
        if !(`amevoid' & strpos("`f2ky'", "ame") == 1) cap ereturn matrix b_es_`f2ky' = `f2bm'
        cap ereturn scalar numpre_`f2ky' = `f2np'
        if `f2ok' == 0 & !(`amevoid' & strpos("`f2ky'", "ame") == 1) dis as text "Caution: no horizon of the `f2ky' event study has an identified variance, so e(b_es_`f2ky') posts without e(V_es_`f2ky'). The estimates stand; anything needing their covariance - honestdid, pretrendpower(), a joint test - has no input on this run."
        if `f2ok' {
          * AMESE R1: withhold e(V_es_ame*) under cell containment — the captured V is
          * not the sampling covariance of the ame estimates (honestdid/lincom/test
          * consumers would ingest it silently); the points still post above
          local f2vrc 111
          if !(`amescrub' & strpos("`f2ky'", "ame") == 1) {
            cap ereturn matrix V_es_`f2ky' = `f2vm'
            local f2vrc = _rc
          }
          if `f2vrc' == 0 & `f2nbad' > 0 {
            dis as text "Caution: e(V_es_`f2ky') contains `f2nbad' of `f2k' horizons with zero or missing variance (structurally degenerate rows); inspect the diagonal before joint inference - honestdid/lincom on those rows is not meaningful."
            cap ereturn scalar es_ndegen_`f2ky' = `f2nbad'
          }
        }
      }
    }
  }
  * F3: expose the ATT(g,t) building blocks [pure-add, after post]
  if "`gteffects'" != "" {
    if `gtn' > 0 {
      cap ereturn matrix gteffects = `gtm'
      if _rc != 0 dis as text "Caution: e(gteffects) could not be posted (rc `=_rc')."
    }
    else dis as text "Caution: gteffects found no cohort-by-relative-period main-effect coefficients to report."
  }
  * F8: expose event-study composition diagnostics [pure-add, after post]
  foreach f8ky of local f8keys {
    local f8mm `f8mat_`f8ky''
    if "`f8mm'" != "" {
      if "`f8ky'" == "main" cap ereturn matrix escomp = `f8mm', copy
      else cap ereturn matrix escomp_`f8ky' = `f8mm', copy
    }
  }
  if strpos("`f8keys'", "main") != 0 {
    cap ereturn scalar escomp_refN = `f8refN'
    cap ereturn scalar escomp_refnunits = `f8refnunits'
    cap ereturn scalar escomp_refncoh = `f8refncoh'
    cap ereturn scalar escomp_refW = `f8refW'
    cap ereturn scalar escomp_trtN = `f8trtN'
    cap ereturn scalar escomp_trtW = `f8trtW'
  }
  * CORE-10: prefix the command name so e(cmdline) matches the Stata convention (the
  * worker's `0' is never rewritten, so this single posting site is authoritative)
  ereturn local cmdline `"wooldid `0'"'
  * posted here rather than as the worker's last statement: a failure in the display or
  * plot stages after the post left e() without e(cmd), so the shell's replay refusal was
  * unreachable and a bare wooldid fell through to a syntax error (bug sweep BS-32)
  ereturn local cmd "wooldid"
  local estimate_type "pooled_reference_period"
  if "`esfixedbaseperiod'" != "" local estimate_type "fixed_reference_period"
  ereturn local estimate_type "`estimate_type'"
  if `anticipation' > 0 ereturn scalar anticipation = `anticipation'





  * F11: ship every remaining result and the display-state macros to e() BEFORE any
  * display runs. the presentation itself now lives in wooldid_display / wooldid_esplots
  * (the worker is near stata's per-program compiled-size ceiling), and those
  * subprograms read e() only - which is also what makes replay (typing wooldid alone
  * after estimation) possible. the posting order below preserves the pre-F11 posting
  * order exactly; the hidden wdid_* display-state macros are appended at the end.
  if "`summarystats'" != "" {
    ereturn matrix summarystats = `summarystats'

    * most summary stats will live in the summary stats matrix, but these may be worth returning on their own
    ereturn scalar ypre_mean = `ypre_mean'
    ereturn scalar ypre_med =  `ypre_med'
    ereturn scalar ypre_sd =  `ypre_sd'
    ereturn scalar ypre_iqr =  `ypre_iqr'

    if "`contreat'" != "" {
      * under untreatedZero, pre-treatment contreat stats are never computed (see the matching guard in the summarystats block), so do not ereturn them
      if "`contreatcontrols'" != "untreatedZero" {
        ereturn scalar contreatpre_mean = `contreatpre_mean'
        ereturn scalar contreatpre_med =  `contreatpre_med'
        ereturn scalar contreatpre_sd =  `contreatpre_sd'
      }
      ereturn scalar contreatpost_mean = `contreatpost_mean'
      ereturn scalar contreatpost_med =  `contreatpost_med'
      ereturn scalar contreatpost_sd =  `contreatpost_sd'
      ereturn scalar contreat_mean = `contreat_mean'
      ereturn scalar contreat_med =  `contreat_med'
      ereturn scalar contreat_sd =  `contreat_sd'
    }
  }

  * return various scalars and summary stats
  local r2list r2 r2adj r2_within r2_withinadj F
  if "`poisson'" == "poisson" local r2list r2pseud chi2
  if "`regtype'" != "poisson" local r2list `r2list' rmse rss
  foreach scalar in N Nclus ll ll_0 `r2list'   {
    if "``scalar''" != "." ereturn scalar `scalar' = ``scalar''
  }
  ereturn scalar N_dropped = `droppedobs'

  * CS-4 (author ruling 2026-09-08): ONE disclosure contract for what e(V) is on this run.
  * The four states the package can post are: "analytic-diagonal" (variances only, no
  * covariance - a cross-estimate test on it is not valid), "analytic-blocks" (as
  * diagonal, plus real within-family event-study blocks where they exist), "bootstrap"
  * (the full replication covariance, with e(df_r) aligned so a hand test reproduces the
  * Wald display), and "absent" (permute's fast-mode default posts none). The marker is
  * machine-readable so downstream code can check rather than infer; the help states what
  * each means and what survives a replay.
  cap confirm matrix e(V)
  local wdhasv = (_rc == 0)
  ereturn hidden local wdid_vcvstate = cond(`wdhasv', "analytic-diagonal", "absent")

  * PS-62 (author rulings 2026-09-07): splice the within-family event-study covariance
  * blocks into e(V) - the one set of off-diagonals wooldid already has and never posted.
  * Runs AFTER the F2 block above (which posts the e(V_es_*) it reads) and after the
  * vcvstate line (which it overwrites with "analytic-blocks" when a block actually lands).
  * Body in a subprogram: wooldid_worker is at the per-program compiled-size ceiling.
  * Skipped inside a replication (wdrepmode): a nested bootstrap/permute call is read for
  * its coefficients, its e(V) DIAGONAL and its e(V_es_*) - never for e(V)'s off-diagonals -
  * so splicing there is pure cost, one Mata call per family per replication (measured ~9%
  * of a reps(150) boottstat fit, under concurrent load)
  if `wdhasv' & "`wdrepmode'" == "" wooldid_vsplice , sgl(`sglist')

  * PS-45 (author ruling 2026-09-07): the simulated analytic simultaneous band, opt-in and
  * defaulted off. Skipped under bootstrap(), where wooldid_bootwire builds the band from
  * the replications instead and would otherwise overwrite this one; refused outright under
  * permute() in wooldid_validate. Body in a subprogram (worker size ceiling)
  if "`simultaneousbands'" != "" & "`bootstrap'" == "" wooldid_esband , mode(analytic) level(`wdlvl') sgl(`sglist')

  * PS-1 (author ruling 2026-09-07): post the conventional names alongside wooldid's own,
  * so estout, esttab and anything else keying on Stata's usual cluster and variance
  * entries finds them. Pure addition - e(Nclus) and e(wdid_vcespec) are unchanged, and
  * nothing that already reads them moves.
  if "`Nclus'" != "." & "`Nclus'" != "" ereturn scalar N_clust = `Nclus'
  * e(clustvar) names the design and stands whether or not a variance was posted; e(vce)
  * and e(vcetype) name what e(V) IS, so a run that posted none leaves them empty rather
  * than labelling an absent matrix (the permutation wire does the same thing)
  if "`cluster'" != "" ereturn local clustvar `cluster'
  if `wdhasv' {
    if "`cluster'" != "" {
      ereturn local vce cluster
      ereturn local vcetype Robust
    }
    else if "`robust'" != "" {
      ereturn local vce robust
      ereturn local vcetype Robust
    }
    else if `"`customvce'"' != "" {
      local wdvce1 = word(`"`customvce'"', 1)
      ereturn local vce `wdvce1'
      if !inlist("`wdvce1'", "ols", "oim", "opg", "unadjusted") ereturn local vcetype Robust
    }
    else if "`unconditionalse'" != "" ereturn local vce unconditional
    else ereturn local vce = cond("`poisson'" == "poisson", "oim", "ols")
  }

  * return the joint  tests
  if "`suppressjointtests'" == "" {
    if "`joint_horiz'" != "" ereturn scalar joint_horiz = `joint_horiz'
    if "`joint_horiz_i'" != "" ereturn scalar joint_horiz_i = `joint_horiz_i'
    if "`joint_Choriz'" != "" ereturn scalar joint_Choriz = `joint_Choriz'
    if "`joint_Choriz_i'" != "" ereturn scalar joint_Choriz_i = `joint_Choriz_i'
  }

 if "`adversarial'" != "" {
   foreach scenario in _ _C {
     foreach type in att att_it {
       * NFADV-01/02: independent gates - the display line survives when inference is
       * unavailable (pooled reference regime, degenerate covariance)
       cap confirm matrix `adversarialB`scenario'`type''
       if _rc == 0 {
         local typetext `type'
         if "`scenario'" == "_C" local typetext = subinstr("`typetext'","att","ame",.)
         if "`adversarialP`scenario'`type''" != "" ereturn scalar adversarial_pval_`typetext' = `adversarialP`scenario'`type''
         ereturn matrix adversarial_line_`typetext' =   `adversarialB`scenario'`type''
       }
     }
   }
 }

  if "`subgroups'" != "" {
    if "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        if "`suppressjointtests'" == "" {
          * CS-60 (author ruling 2026-09-08): post these as numeric SCALARS, matching the
          * overall joint families. They were the only joint p-values stored as string
          * macros, so a single loop over subgroups got numbers for one family and strings
          * for the other. A p-value's conventional type is a scalar.
          if "`joint_horiz_s`sgcase''" != "" ereturn scalar joint_horiz_s`sgcase' = `joint_horiz_s`sgcase''
          if "`joint_horiz_i_s`sgcase''" != "" ereturn scalar joint_horiz_i_s`sgcase' = `joint_horiz_i_s`sgcase''
          if "`joint_Choriz_s`sgcase''" != "" ereturn scalar joint_Choriz_s`sgcase' = `joint_Choriz_s`sgcase''
          if "`joint_Choriz_i_s`sgcase''" != "" ereturn scalar joint_Choriz_i_s`sgcase' = `joint_Choriz_i_s`sgcase''
        }
        if "`adversarial'" != "" {
          foreach scenario in _ _C {
            foreach type in att att_it {
              * NFADV-01/02: independent gates (see the overall twin above)
              cap confirm matrix `adversarialB`scenario'`type'_s`sgcase''
              if _rc == 0 {
                local typetext `type'
                if "`scenario'" == "_C" local typetext = subinstr("`typetext'","att","ame",.)
                if "`adversarialP`scenario'`type'_s`sgcase''" != "" ereturn scalar adversarial_pval_`typetext'_s`sgcase' = `adversarialP`scenario'`type'_s`sgcase''
                ereturn matrix adversarial_line_`typetext'_s`sgcase' =   `adversarialB`scenario'`type'_s`sgcase''
              }
            }
          }
        }
      }
    }
    if wordcount("`sglist'") > 1 {
      ereturn matrix pwsgs = `pwsgs'
      foreach item in `att'  `att_it' `att_i' `att_itime' {
        if "`joint_`item'_sdiff'" != "" ereturn scalar joint_`item'_allsgssame =   `joint_`item'_sdiff'
        if "`prejoint_`item'_sdiff'" != "" ereturn scalar joint_pre_`item'_allsgssame =   `prejoint_`item'_sdiff'
      }
      if "`contreat'" != "" {
        foreach item in `att'  `att_it' `att_i' `att_itime' {
          local itemtext = subinstr("`item'","att","ame",.)
          if "`joint_C`item'_sdiff'" != "" ereturn scalar joint_`itemtext'_allsgssame =   `joint_C`item'_sdiff'
          if "`prejoint_C`item'_sdiff'" != "" ereturn scalar joint_pre_`itemtext'_allsgssame =   `prejoint_C`item'_sdiff'
        }
      }
    }
  }

  *return the histogram matrix
  if "`histogramcohorteffects'" != "" ereturn matrix histogramestimates = `___hist'

  * ship remaining matrices into ereturn
  if "`subgroups'" != "" & "`suppresssgprimaryeffects'" == "" {
    foreach sgcase in `sglist' {
      ereturn matrix wdidsgresults_`sgcase'  = `wdid_sg`sgcase''
    }
  }

  ereturn matrix wdidmainresults = `wdidpooled'

  *** finally, return the CFX plots
  * NFFAST-02: the margins-inherited row stripes embed wooldid tempvar names (volatile
  * across runs - two identical calls returned identical numbers under different row
  * identifiers, breaking documented downstream row-name joins); rename to stable
  * lattice_# keys before posting
  ereturn hidden local wdid_cfxmats `cfxmatricestoreturn'
  * P7 accident 2 (2026-09-08): on the ANALYTIC path the cfx figures come from marginsplot and
  * saveplots() wrote only the .gph files - the <name>_cfx.dta of plotted points was produced
  * under an engine and nowhere else. Record the stem here; the shell writes the file after
  * e() is posted, since the writer reads e(). Posted only when there is something to write, so
  * no existing run's e() surface changes
  if "`saveplots'" != "" & "`makecfxplots'" != "" & "`wdrepmode'" == "" & "`bootstrap'`permute'" == "" & "`cfxmatricestoreturn'" != "" ereturn hidden local wdid_cfxsave "`saveplots'"
  if "`e(wdid_cfxsave)'" != "" & "`replace'" != "" ereturn hidden local wdid_cfxrepl replace
  foreach matrix in `cfxmatricestoreturn' {
    wooldid_cfxstripe ``matrix''
    ereturn matrix `matrix' = ``matrix''
  }
  * per-plot metadata for wooldid_cfxplots, indexed by position in e(wdid_cfxmats): the
  * titles the worker composed, the x-axis title option verbatim, and the user's lattice
  * variable (empty under an auto lattice, whose tempvar does not outlive the worker)
  if `cfxnp' > 0 {
    forval wdk = 1/`cfxnp' {
      ereturn hidden local wdid_cfxttl`wdk' `"`cfxttl`wdk''"'
      ereturn hidden local wdid_cfxytx`wdk' `"`cfxytx`wdk''"'
    }
    ereturn hidden local wdid_cfxxt `"`cfxxt_meta'"'
    ereturn hidden local wdid_cfxlatvar "`cfxlattice'"
  }

  * display-state macros consumed by wooldid_display / wooldid_esplots (and by replay).
  * hidden: internal plumbing, not part of the documented e() surface (they still appear
  * to "local m : e(macros)" and so to the test battery). an ereturn local with an empty
  * value creates no macro, so the unconditional lines below are self-gating.
  * PS-20 (author ruling 2026-09-07): the Ibragimov-Muller path names its method in the header
  * the way bootstrap() and permute() name theirs. Under impvalues the reported p-values and
  * intervals are IM quantities on G-1 degrees of freedom; the header used to print the analytic
  * cluster specification, naming a variance the reported inference did not use. The cluster
  * specification is kept in the string, because it is what defines the IM groups
  ereturn hidden local wdid_vcespec `"`=cond("`impvalues'" != "", "Ibragimov-Muller p-values over " + cond("`cluster'" != "", "cluster(`cluster')", "the resampling groups") + " (t on G-1 df)", `"`vcespecifier'"')'"'
  ereturn hidden local wdid_poisson "`poisson'"
  ereturn hidden local wdid_poisexp "`poisexpresults'"
  ereturn hidden local wdid_contreat "`contreat'"
  ereturn hidden local wdid_sgvar "`subgroups'"
  ereturn hidden local wdid_sglist "`sglist'"
  * F2 (2026-09-14): the subgroups whose pre-treatment row was withheld, why, and which rows
  * went, written as the sibling of wdid_sgdegen/wdid_sgdegenr/wdid_sgdegenact - read back by
  * wooldid_display so replay reproduces the caution. the three lists are parallel (the k-th
  * reason and k-th row token belong to the k-th subgroup). reasons: "nocells" (the subgroup
  * realizes no pre-period cells outside the reference period) and "omitted" (it does, but
  * margins marked its level omitted). rows: "all" (every pre-treatment row that subgroup
  * would have had) and "ame" (only the continuous-treatment pre-treatment marginal effect;
  * its pre-treatment att row is present). an ereturn local with an empty value creates no
  * macro, so a run that withheld nothing posts none of the three and the caution cannot fire
  ereturn hidden local wdid_sgnopre "`sgnopre'"
  ereturn hidden local wdid_sgnoprewhy "`sgnoprewhy'"
  ereturn hidden local wdid_sgnoprerow "`sgnoprerow'"
  ereturn hidden local wdid_estimands "`att' `att_it' `att_i' `att_itime'"
  ereturn hidden local wdid_cwm "`customweightmultiplier'"
  ereturn hidden local wdid_semielas "`semielasticity'"
  * BS-6 (2026-09-07): the macro above holds the OPERATOR margins used; this one holds what
  * the user asked for, and is posted ONLY when the remap moved it - so no existing run's
  * e() surface changes and no baseline macro dump moves
  if "`semielasraw'" != "`semielasticity'" ereturn hidden local wdid_semielasraw "`semielasraw'"
  ereturn scalar level = `wdlvl'
  ereturn hidden local wdid_fcol "`fcon'"
  * F061 (2026-09-06): the realised cell count (missing when the collapse did not run, so
  * "no collapse" and "collapsed to zero cells" cannot be confused) and the skip reason(s)
  if "`fconn'" != "" ereturn hidden scalar wdid_fcncell = `fconn'
  else ereturn hidden scalar wdid_fcncell = .
  ereturn hidden local wdid_fcskip "`fcskip'"
  ereturn hidden local wdid_fcguard "`fcguard'"
  ereturn hidden local wdid_cetype "`contreatelasticitytype'"
  * BS-6 (2026-09-07): the user's request, posted only when the xb-scale remap moved it
  if "`cetyperaw'" != "`contreatelasticitytype'" ereturn hidden local wdid_cetyperaw "`cetyperaw'"
  ereturn hidden local wdid_suppressjt "`suppressjointtests'"
  ereturn hidden local wdid_suppresssgpe "`suppresssgprimaryeffects'"
  ereturn hidden local wdid_adversarial "`adversarial'"
  if "`summarystats'" != "" ereturn hidden local wdid_sumstats "summarystats"
  ereturn hidden local wdid_espre "`espre'"
  ereturn hidden local wdid_nptb "`Npre_trebase'"
  ereturn hidden local wdid_nprat "`Npre_rat'"
  ereturn hidden local wdid_wprat "`Wpre_rat'"
  ereturn hidden local wdid_esrel "`esrelativeto'"
  ereturn hidden local wdid_timetrends "`timetrends'"
  ereturn hidden local wdid_fastnose "`fastnose'"
  * replay contract version: bump this when the display contract changes (anything
  * wooldid_display / wooldid_esplots read from e(), or how they interpret it);
  * their gates exact-match this value, so a bump refuses cross-version replays
  ereturn hidden local wdid_replayv "2"


  ********* begin actual result presentation *************
  * F11: the presentation code lives in wooldid_display (non-eclass, reads e() only).
  * it is split in two parts so the eclass wooldid_ptp call below stays at worker level
  * (its ereturns would be discarded inside a non-eclass subprogram)
  * F-boot deferred display: under bootstrap the shell rewires e() to the bootstrap
  * covariance first and then calls the (e()-driven) display itself - printing here
  * would show analytic inference that is about to be replaced. pretrendpower is
  * refused with bootstrap at parse (plan A13), so the ptp call needs no gate
  * BS-15 (author ruling 2026-09-07): the realized-degeneracy screen on the subgroup pairwise
  * covariance. It rewrites e(pwsgs) and the allsgssame scalars, so it must run AFTER they are
  * posted and BEFORE anything displays them. Analytic path only: from here the bootstrap and
  * permutation shells own those cells and read the analytic pairwise standard error as their
  * own input (wooldid_bootjnt's `pw0', wooldid_permwire's se_analytic column). fastnose is
  * named rather than left to the missing-standard-error test: F9's contract is "." cells
  if "`bootstrap'`permute'" == "" & "`fastnose'" == "" & wordcount("`sglist'") > 1 wooldid_sgdeg , si(`safetyinput')
  if "`bootstrap'`permute'" == "" wooldid_display , part(1)
  * F7: Roth (2022) pre-trend power diagnostics (compute + display + ereturns live in a
  * subprogram; cap noi so a failure here cannot abort an otherwise-complete run)
  * F7: Roth (2022) pre-trend power diagnostics. Without bootstrap the whole diagnostic
  * runs here. WITH bootstrap only the data-dependent half runs here - it needs `ttre'
  * and `w', which are worker tempvars - and the shell re-enters this program after the
  * bootstrap rewiring so the power algebra reads the BOOTSTRAP e(V_es_*)
  if "`pretrendpower'" != "" & "`bootstrap'" == "" cap noi wooldid_ptp , levels(`pretrendpower') esrel(`esrelativeto') dfr(`df_r') pois(`poisson') poisexp(`poisexpresults') tvar(`t') ttrevar(`ttre') ivar(`i') wvar(`w') cwm(`customweightmultiplier') level(`wdlvl')
  if "`pretrendpower'" != "" & "`bootstrap'" != "" cap noi wooldid_ptp , levels(`pretrendpower') esrel(`esrelativeto') dfr(`df_r') pois(`poisson') poisexp(`poisexpresults') tvar(`t') ttrevar(`ttre') ivar(`i') wvar(`w') cwm(`customweightmultiplier') hbaronly level(`wdlvl')
  if "`bootstrap'`permute'" == "" wooldid_display , part(2)




  ***************************************************************************************************************************
  ********************************************         PRODUCE EVENT STUDY PLOTS       **************************************
  ***************************************************************************************************************************
  * present event study plots if requested; a little complicated because of naming practices in wdidpooled
  * this is also somewhat awkward placement b/c I have this coded to use the wdidpooled matrix
  * this choice is vestigial to a prior program design and is kept here on the principle of not redoing what already works
  * though I should probably redo this someday

  * Under bootstrap this call is DEFERRED to the shell, which runs it after
  * wooldid_bootwire. wooldid_esplots reads e(wdidmainresults) and e(wdidsgresults_*),
  * and running it here would read them BEFORE the rewiring - so the figure's bands were
  * analytic while the table's were bootstrap. Measured on one fixture, the two disagreed
  * about significance (ES_att pre2: table [-0.9951, 0.0113] includes zero, plot
  * [-1.1353, -0.1398] excludes it). Author ruling 2026-09-04. Because esplots reads e()
  * rather than worker tempvars, deferring needs no stash - unlike wooldid_ptp
  if `espre' +`espl' > 0 & "`makeplots'" == "makeplots" & "`bootstrap'`permute'" == "" {
    local wdesopts
    if "`saveplots'" != "" local wdesopts saveplots(`saveplots')
    wooldid_esplots , `wdesopts' `replace'
  }



  ***************************************************************************************************************************
  ***************************************         CLEAN UP AND EXIT       **************************************
  ***************************************************************************************************************************
  * (F11) all result matrices and scalars were shipped to e() above, BEFORE the display
  * and plot calls; only cleanup remains here

  if "`cleanmatrices'" != "" {
    clear mata
    matrix drop _all
  }

  if "`verbose'" != "verbose" | "`wdrepmode'" != "" qui est drop `wdfm'
  * clean up after ppmlhdfe. under verbose, retain _ppmlhdfe_d alongside the retained
  * ___wooldidfullmodel (author directive, 2026-07). note this alone does not make
  * predict/margins replay work after est restore: the model's regressors (_i, _t,
  * targetflag, and the internal weight) are tempvars that no longer exist at that
  * point - retention just stops closing that door for a future replay feature
  if "`poisson'" == "poisson" & "`unconditionalse'" == "" & "`customvce'" == "" & "`verbose'" != "verbose" qui cap drop _ppmlhdfe_d


  * program ends!
  dis as result  " "




  end


* ============================================================================
* S5: the complete early input gate - every parse-time check and option
* normalization wooldid performs before touching the dataset structure. moved out
* of wooldid_worker (per-program compiled-size ceiling) verbatim and in the original
* order, so multi-error calls report the same first error as before. re-parses `0'
* with a byte-identical copy of the worker's syntax line (test/check_validate_sync.py
* guards the duplication) and returns everything it parses or mutates via c_local.
* NOT in here (hard constraints, all measured): set emptycells (program-scoped, would
* be restored the instant this returns), tempvar creation (dies on return), and the
* update branch (its exit would only leave this subprogram).
cap prog drop wooldid_validate
prog def wooldid_validate
version 16.0

* S5: this line must remain BYTE-IDENTICAL to the worker's syntax line
syntax [varlist(numeric default=none)] [if/]  [aw fw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist numeric) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose fastcollapse fullprecisionslow GTEffects PERMute(string) PERMREPS(integer 0) PERMEXact PERMTStat PERMTci PERMSTRata(varname numeric) PERMPRobs(string asis) PERMSAve(string) PERMANalyticse PERMABsp PERMSide(string) PERMRAWp PERMSTEPdown SIMULTaneousbands BOOTstrap(string) reps(integer 0) ties Level(real 0) dots(integer -1) rseed(string) BOOTALPha(string) BOOTDf(string) BOOTDISPlay(string) BOOTTStat BOOTSAve(string)]

* Verify that the packages wooldid requires (reghdfe, ftools, gtools) are installed before doing anything else
* ppmlhdfe is deliberately not checked here: it is only required for poisson estimation and remains optional
  local missingpackages
  foreach requiredpackage in reghdfe ftools gtools {
    cap which `requiredpackage'
    if _rc != 0 {
      dis as error "Error: wooldid requires the package `requiredpackage', which does not appear to be installed."
      dis as text "To install it, run: ssc install `requiredpackage', replace"
      local missingpackages `missingpackages' `requiredpackage'
    }
  }
  if "`missingpackages'" != "" error 198

  if "`verbose'" == "verbose" dis as text "Wooldid is processing input in preparation for regression estimation."

 * This controls the names we use for ES plots; the new syntax is more compatible with the program event_plot
 if "`oldsyntax'" == "oldsyntax" {
    local oldsyntax 1
 }
 else {
    local oldsyntax 0
 }

  **** Begin a block of checks that command info is correctly specified; error out if not. handle all syntax errors here

  * check that sufficient input is given - note that the error will occur automatically if > 4 is given
  if wordcount("`varlist'") < 4 {
    dis as error "Error: command line input should feature all of y i t ttre; too few variables specified."
    error 198
  }

  if wordcount("`varlist'") > 4 {
    dis as error "Error: command line input should feature only y i t ttre; too many variables specified."
    error 198
  }

  *Confirm correct weight usage if using poisson
  if "`weight'" == "aweight" & "`poisson'" == "poisson" {
    dis as error "Error: poisson only accepts pweights, not aweights."
    error 198
  }

  *Issue warnings related to unconditional se reliability
  if "`safety'" != "off" & "`unconditionalse'" == "unconditionalse" & ("att_i" == "`att_i'" | "att_it" == "`att_it'" | "att_itime" == "`att_itime'" ) {
    dis as error "Error: unconditional standard errors may be unreliable with non-att estimand types."
    dis as text "Set option safety(off) to disregard this error."
    error 198
  }

  if "`safety'" != "off" & "`unconditionalse'" == "unconditionalse" & ("`weight'" == "aweight" ) {
    dis " "
    dis as text "Caution: unconditional standard errors may be unreliable with aweights."
    dis as text "Set option safety(off) to suppress this caution."
    dis " "
  }

  if "`adversarial'" == "adversarial" & "`impvalues'" == "impvalues" {
    dis as error "Error: Adversarial Pre-Trends Testing not allowed with Ibragimov and Muller p-values."
    error 198
  }

  * object to certain specification types when forced to switch to reg or poisson over an hdfe estimator
  if ("`customvce'" != "" | "`unconditionalse'" != "") & ( "`fe'" != "" | "`ccc_absorb'" != "")  {
    dis as text "Caution: unconditional SEs and custom VCE specifications rely on the plain regress command and so are not compatible with absorbed fixed effects."
    dis as text "You must specify anything in fe as controls, manually specifying any needed i. prefixes."
  }

  * check conditions for using im p-values
  if "`impvalues'" != "" & "`cluster'" == "" {
    dis as error "Error: Ibragimov and Muller p-values only allowed when clustering; please specify a clustering variable."
    error 198
  }

  * unconditionalse estimates through regress/poisson and takes its variance from margins,
  * vce(unconditional), which accepts one cluster variable; a second one reached regress as
  * a raw "too many variables specified" (bug sweep BS-14)
  if "`unconditionalse'" != "" & wordcount("`cluster'") > 1 {
    dis as error "Error: unconditionalse supports a single cluster() variable; multiway clustering is not available with unconditional standard errors."
    dis as text  "unconditionalse estimates through regress/poisson and obtains the variance from margins, vce(unconditional), which accepts only one cluster variable."
    error 198
  }
  if "`impvalues'" != "" & wordcount("`cluster'") > 1 {
    dis as error "Error: Ibragimov and Muller p-values only allowed with a single cluster variable."
    error 198
  }

  * check conditions for using im p-values
  if "`impvalues'" != "" & ( "`customvce'" != "" | "`unconditionalse'" != "") {
    dis as error "Error: Ibragimov and Muller p-values are not compatible with unconditional SE computation and custom VCE specifications."
    error 198
  }

  * complain if you specify customvce and something else with it
  if "`customvce'" != "" & ("`robust'" != "" | "`cluster'" != "") {
    dis as error "Error: customvce() is not compatible with robust and cluster(). You must manually specify the contents of reg's vce() statement in this case."
    error 198
  }

  * unconditional SEs are computed by margins via vce(unconditional), which requires the underlying
  * regression to have been fit with robust or clustered standard errors; catch this early rather than
  * failing with a bare r(322) at the first margins call after estimation has completed.
  * pweights are exempt: they force a robust vce in reg/poisson, satisfying the requirement on their own
  if "`unconditionalse'" == "unconditionalse" & "`robust'" == "" & "`cluster'" == "" & "`customvce'" == "" & "`weight'" != "pweight" {
    dis as error "Error: the option unconditionalse requires specifying either robust or clustered standard errors."
    dis as text "Specify robust or cluster(), or provide a robust/clustered vce via customvce()."
    error 198
  }


  *check conditions for using im p-values - in principle, SG estimates could be computed provided SGs never cross cut clusters, but I have not implemented this
  * NOTE (t10_v3): this gate is LOAD-BEARING for the E1 subgroup restructure - under
  * impvalues, marginsdof is nose, and a chain `test' after a nose `post' is r(321);
  * do not relax this pairing without also emptying `sgpost' in the SG-comparison blocks
  if "`impvalues'" != "" & "`subgroups'" != "" {
    dis as error "Error: Ibragimov and Muller p-values not compatible with subgroups."
    error 198
  }

  * F3: under a continuous treatment the cohort-by-period interaction coefficients are
  * intercepts at dose zero, not treatment effects, so the building-block table would be
  * silently meaningless there
  if "`gteffects'" != "" & "`contreat'" != "" {
    dis as error "Error: gteffects is not available with contreat(): the cohort-by-period interaction coefficients are dose-zero intercepts under a continuous treatment, not treatment effects. Aggregate e(b_es_ame) or the cfx materials instead."
    error 198
  }


  *check conditions for using im p-values - for CFX plots, we have problems similar to with SGs, but even worse
  if "`impvalues'" != "" & "`makecfxplots'" != "" {
    dis as error "Error: Ibragimov and Muller p-values not compatible with continuous treatment effects plots."
    error 198
  }

  * F7 pre-trend power diagnostics: validate the request early (checks live in a
  * subprogram - the worker is close to Stata's per-program size ceiling)
  if "`pretrendpower'" != "" wooldid_ptpcheck , levels(`pretrendpower') esf(`esfixedbaseperiod') impv(`impvalues') espre(`esprelength')


  * check conditions for using im p-values --
  * for the joint tests, I have not seen IM characterize multi-way tests. implementing a test here of some sort would be easy,
  * but I am not sure if it would technically be valid.
  if "`impvalues'" != "" & ( "`jointtests'" != "" ){
    dis as error "Error: Ibragimov and Muller p-values are not compatible with the jointtests option."
    error 198
  }

  * verify that if continuous treatment effects plots are requested, we have a lattice over which to produce them
  if "`makecfxplots'" != "" & "`cfxlattice'" == "" & ("`lattice'"=="" ) {
    dis as error  "Error: You requested continuous treatment effects plots, but either didn't specify a cfxlattice variable or a lattice option for the continuous treatment."
    error 198
  }

  * oostreatedcontrols determines whether or not units treated after the end of the sample window are allowed to contribute to pre-treatment horizons;
  * cancel estimation if you specify a negative value here (corresponding with a within window treatment)
  if `oostreatedcontrols' < 0 {
    dis as error  "Error: oostreatedcontrols specified as < 0, indicating that treatments occurring within sample should be ignored."
    dis as error  "This is not allowed; you must specially set up cases of this sort by modifying the ttre variable and the if statement."
    error 198
  }

  * anticipation() shifts the effective onset EARLIER (ttre - #); a negative value would
  * pretend onset is later, pushing genuinely treated periods into the control/reference
  * pool - a silent identification failure. delayed-onset designs need their own setup
  if `anticipation' < 0 {
    dis as error  "Error: anticipation specified as < 0, indicating that treatment onset should be treated as occurring later than recorded."
    dis as error  "This is not allowed; if you need a delayed-onset design, modify the ttre variable and the if statement directly."
    error 198
  }


  * check that controls are specified alright - error out if c. and i. are not specified for each variable
  if "`coarsecohortcontrols'" != "" {
    foreach case in `coarsecohortcontrols' {
      if strpos("`case'","c.") == 0 & strpos("`case'","i.") == 0 {
        dis as error  "Error: Your coarse cohort controls do not all specify whether they are continuous or categorical. You must do this using i. and c. style notation."
        error 198
      }
    }
  }

  * extract information about requested continuous treatment and error out if the user misspecified the information
  if "`lattice'" != "" & wordcount("`lattice'") != 1 {
    if wordcount("`lattice'") >= 3  {
      dis as error  "Error: lattice incorrectly specified; requires 2 arguments or a pre-made lattice variable"
      error 198
    }

    local latticenum = word("`lattice'",2)
    * validate by VALUE, not spelling: every downstream quantile engine accepts
    * integer-valued spellings like 02, +2, 2.0, and 1e2, so a confirm-integer-number
    * spelling check would reject inputs that work; conversely gquantiles silently
    * remaps nquantiles <= 0 to its default (and forval builds a degenerate no-op
    * lattice), so 0/negatives must be caught here by value. missing-evaluating
    * tokens fail confirm number or the mod() test (mod(., 1) != 0 is true)
    cap confirm number `latticenum'
    if _rc != 0 {
      dis as error "Error: lattice incorrectly specified; `latticenum' is not a valid number of lattice quantiles. lattice() requires a valid lattice type and an integer number of quantiles, or a pre-made lattice variable."
      error 198
    }
    if mod(`latticenum',1) != 0 {
      dis as error "Error: lattice incorrectly specified; the number of lattice quantiles (`latticenum') must be an integer."
      error 198
    }
    if `latticenum' <= 1 {
      dis as error "Error: lattice incorrectly specified; latticenum must be > 1 or else the lattice does not permit additional variation relative to not specifying a lattice at all."
      error 198
    }

    local latticetype = word("`lattice'",1)
    if inlist("`latticetype'","all","baseTreated","baseTreatedPost","baseTreatedPre", "baseTreatedRef") == 0 & inlist("`latticetype'","byCohort","byTime","byCohortByTime","byTreated","byTreatedByTime","byTtre","byTtreByTime") == 0{
      dis as error  "Error: lattice type is misspecified; `latticetype' is not a valid lattice type."
      error 198
    }
  }

  * extract contreatcontrol type and error out if misspecified
  if "`contreatcontrols'" != "" {
    foreach word in `contreatcontrols' {
      if inlist("`word'",  "byTreated", "byTtre", "byCohort", "byTime") == 0 & inlist("`word'","basic","untreatedZero") == 0 {
        dis as error  "Error: continuous treatment control type misspecified; `word' is not a valid continuous treatment control type."
        error 198
      }
    }
    if strpos("`contreatcontrols'","untreatedZero") != 0 & wordcount("`contreatcontrols'") > 1 {
      dis as error "Error: untreatedZero cannot be paired with other control types"
      error 198
    }
  }

  * error out if no contreatcontrol type is specified in a continuous treatment scenario
  if "`contreatcontrols'" == "" & "`contreat'" != "" {
    dis as error "Error: contreatcontrols type not specified; you must specify this - substantive judgement is required."
    error 198
  }

  * error out if you request continuous treatment controls w/o a continuous treatment variable
  if "`contreatcontrols'" != "" & "`contreat'" == "" {
    dis as error "Error: you cannot request continuous treatment controls without specifying a continuous treatment variable."
    error 198
  }

  * error out if a lattice is requested without a continuous treatment - the lattice-construction
  * code all sits inside the contreat block, so the option would otherwise be silently ignored
  if "`lattice'" != "" & "`contreat'" == "" {
    dis as error "Error: lattice() has no effect without contreat(); specify a continuous treatment variable or remove the lattice() option."
    error 198
  }

 * issue warnings pertaining to cfxplot types
  * CS-98 (author ruling 2026-09-08): resolve the EFFECTIVE safety() value before any
  * message keyed on it prints. contreatwithin with contreat() relaxes safety further down;
  * doing it there let the caution below announce a refusal that the relaxation then
  * overturned, so the run printed "will not produce these plots" and produced them. The
  * user's own value is kept as `safetyuser' for the checks that must see what was typed.
  local safetyuser `safety'
  if "`contreatwithin'" == "contreatwithin" & "`contreat'" != "" local safety off
  if "`makecfxplots'" != "" & ("`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" )   {
    if "`safety'" == "off" {
      dis as text "Caution: interpretation of continuous treatment effects plots for estimands other than the ATT is difficult unless i-t cells (and sometimes cohorts i) are wholly contained within lattice bins."
    }
    else{
      dis as text "Caution: interpretation of continuous treatment effects plots for estimands other than the ATT is difficult unless i-t cells (and sometimes cohorts i) are wholly contained within lattice bins."
      dis as text "As such, wooldid will not produce these plots unless you specify safety(off)."
    }
  }

  if `esprelength' < 0 {
    dis as error "Error: cannot accept negative values of esprelength."
    error 198
  }

  * this implements a hack for estimability check failures in certain special cases
  if `espostlength' < 0  {
    if `espostlength' != -999  {
      dis as error "Error: cannot accept negative values of espostlength."
      error 198
    }
    local espostlength 0
    local esposttrip no
  }

  if "`subgroups'" != "" &   `esprelength' > 0  & `espostlength' == 0 & "`esposttrip'" == "" {
    local espostlength 1
    * P6 (author ruling 2026-09-07): the line above RAISES the estimand window - the run reports
    * a post-treatment horizon the user did not ask for - so the fact of the change is announced
    * whatever safety() says. safety(off) suppresses the advice, never the disclosure: a silent
    * change to WHAT IS BEING ESTIMATED is not a safety check, and the option has never meant
    * "stop telling me what you did".
    dis as text "Caution: espostlength set to 1 due to use of subgroups in conjunction with esprelength > 0."
    if "`safety'" != "off" dis as text "Without this modification, subgroups may crash if pre-treatment effects are not estimable for them."
    if "`safety'" != "off" dis as text "Specify espostlength(-999) to override this issue. Specify safety(off) to suppress the explanation; the line above prints either way."
  }



  ******** parsing for specification errors complete, begin processing input on other dimensions, setting defaults, etc. handle errors that aren't purely syntax here

  ******* process the raw inputs

  * extract y, i, t, and ttre from the main variable list
  local y = word("`varlist'",1)
  local i = word("`varlist'",2)
  local t = word("`varlist'",3)
  local ttre_asgiven = word("`varlist'",4)

  * F6: permutation-inference screens. deep checks live in the subprogram
  * (wooldid_permctl) to conserve this program's compiled-size headroom
  if "`permute'" == "" {
    local permstray
    if `permreps' != 0 local permstray permreps()
    if "`permexact'" != "" local permstray `permstray' permexact
    if "`permtstat'" != "" local permstray `permstray' permtstat
    if "`permtci'" != "" local permstray `permstray' permtci
    if "`permstrata'" != "" local permstray `permstray' permstrata()
    if `"`permprobs'"' != "" local permstray `permstray' permprobs()
    if "`permsave'" != "" local permstray `permstray' permsave()
    if "`permanalyticse'" != "" local permstray `permstray' permanalyticse
    if "`permabsp'" != "" local permstray `permstray' permabsp
    if "`permside'" != "" local permstray `permstray' permside()
    if "`permrawp'" != "" local permstray `permstray' permrawp
    if "`permstepdown'" != "" local permstray `permstray' permstepdown
    if "`permstray'" != "" {
      dis as error "Error: `permstray' require(s) the permute() option."
      error 198
    }
  }
  * PS-45 (author ruling 2026-09-07): a simultaneous band is a CONFIDENCE BAND, and the
  * permutation engine's multiplicity object is the max-|studentized| step-down instead
  * (same ruling: "stays permutation-only"). Posting an analytic band under a permutation
  * header would also break the 2026-09-05 rule that only permutation quantities are shown
  * under permute(), so this combination is refused rather than silently dropped
  if "`permute'" != "" & "`simultaneousbands'" != "" {
    dis as error "Error: simultaneousbands is not available under permute(). A simultaneous band is a confidence-band object; the permutation engine's multiplicity tool is the max-|studentized| step-down, permstepdown. Use bootstrap() for a bootstrap band, or drop permute() for the simulated analytic one."
    error 198
  }

  * permtci's critical value comes from the permuted |t| distribution, so it implies
  * permtstat. Resolved HERE (it used to live in the shell alone) so that the fast-mode
  * default further down sees the implied value: permtci alone must not send the observed
  * run down the no-SE path it is about to be studentized against
  if "`permtci'" != "" & "`permtstat'" == "" local permtstat permtstat
  * NOTE which safety() value each screen gets, because the two calls differ. permctl is
  * handed si(`safety') - the EFFECTIVE value, already rewritten to off at the CS-98 site
  * above when contreatwithin is combined with contreat() - while bootctl below is handed
  * si(`safetyinput'), the user's own. They coincide on every run that reaches permctl's
  * subgroup screens only because permctl refuses contreatwithin outright at its top. If
  * that refusal is ever lifted (the contreat() refusal beside it was already adjudicated
  * the other way), a contreatwithin run would arrive here with si(off) and the subgroup
  * cell-containment and unit-splitting screens would silently downgrade from error 459 to
  * a caution - the exact failure the safetyinput split exists to prevent. Pass
  * si(`safetyuser') at that point, as bootctl already does
  if "`permute'" != "" {
    cap noi wooldid_permctl , permute(`permute') permreps(`permreps') `permexact' clus(`cluster') ttrevar(`ttre_asgiven') absp(`permabsp') side(`permside') rawp(`permrawp') permstrata(`permstrata') permprobs(`permprobs') impv(`impvalues') uncond(`unconditionalse') cw(`contreatwithin') os(`oldsyntax') ts(`permtstat') tc(`permtci') fn(`fastnose') cc(`contreatcontrols') ifr(`if') ptp(`pretrendpower') ssg(`suppresssgprimaryeffects') sg(`subgroups') ivar(`i') tvar(`t') si(`safety') psave(`permsave') level(`level') `replace'
    if _rc != 0 exit _rc   // RF-1: exit, not error - the cap noi above already printed it (rule at :61)
  }

  * capture the user's own safety() input, unaffected by the contreatwithin relaxation that
  * is now resolved earlier (CS-98); the identifying-assumption checks gate on this copy so
  * an internal relaxation does not silently disable them
  local safetyinput `safetyuser'

  * F-boot: bootstrap-SE screens; deep checks live in the subprogram (wooldid_bootctl)
  * to conserve this program's compiled-size headroom. called UNCONDITIONALLY: with
  * bootstrap() absent it enforces the stray-companion screen and exits. it sits AFTER
  * the safetyinput capture above deliberately - the bayes cell-weighted gate keys on
  * the user's own safety() value (contreatwithin rewrites the working copy later)
  cap noi wooldid_bootctl , bootstrap(`bootstrap') reps(`reps') `ties' level(`level') dots(`dots') rseed(`rseed') bootalpha(`bootalpha') bootdf(`bootdf') bdisp(`bootdisplay') bts(`boottstat') bsave(`bootsave') pois(`poisson') permute(`permute') impv(`impvalues') uncond(`unconditionalse') cvce(`customvce') fn(`fastnose') os(`oldsyntax') ssg(`suppresssgprimaryeffects') attit(`att_it') atti(`att_i') attitime(`att_itime') si(`safetyinput') clus(`cluster') ivar(`i') ifr(`if') jt(`jointtests') ptp(`pretrendpower') ttrevar(`ttre_asgiven') `replace'
  if _rc != 0 exit _rc   // RF-1: exit, not error - the cap noi above already printed it (rule at :61)

  * register verbose retention with the shell's error-path cleanup epilogue (top of this
  * file): under verbose, ___wooldidfullmodel (and _ppmlhdfe_d under poisson) are
  * documented as retained for inspection - on failed runs as well
  if "`verbose'" == "verbose" global WOOLDID_CLEANUP $WOOLDID_CLEANUP keep

  * the input-validity asserts below respect the user's if() restriction, mirroring the
  * estimation sample construction exactly: the sample keeps rows where the if expression
  * evaluates to 1 (gen ifvar = `if' further down), so the asserts test (`if') == 1 --
  * identical row selection for any expression shape, including bare flag variables with
  * nonzero or missing values, which a plain if-qualifier would wrongly include. empty when
  * no if() is given, leaving the asserts running over the full dataset exactly as before
  local ifassert
  if `"`if'"' != "" local ifassert if (`if') == 1

  * check that i is specified acceptably
  if "`safety'" != "off" {
    cap assert `i' > 0  & int(`i') == `i' `ifassert'
    if _rc != 0 {
      dis as error  "Error: `i' must consist of positive, non-zero integers, but does not"
      error 411
    }
  }

  * float storage cannot represent consecutive integers above 2^24 = 16,777,216, so an
  * 8-digit calendar code, a %tc datetime or a large ID held in float was rounded BEFORE
  * wooldid saw it: adjacent periods, cohorts or IDs merge with no other symptom and the run
  * completes with a wrong answer. no exact detector exists (the rounding happened
  * upstream) and coarse codes above the threshold are representable, so this is a caution
  * on the only available signal - storage type plus magnitude (bug sweep BS-34)
  if "`safety'" != "off" {
    foreach wdfv in `i' `t' `ttre_asgiven' {
      if "`: type `wdfv''" == "float" {
        qui gstats sum `wdfv' `ifassert', meanonly
        if r(N) > 0 & max(abs(r(min)), abs(r(max))) > 16777216 & max(abs(r(min)), abs(r(max))) < . {
          dis as text "Caution: `wdfv' is stored as float and holds values above 16,777,216 (2^24), beyond which float cannot represent consecutive integers. Its values were rounded before wooldid saw them, which can silently merge adjacent periods, cohorts or IDs. If it is an 8-digit calendar code or a %tc datetime, regenerate it as double from its source (recasting cannot recover the lost values); coarse or even-valued codes above the threshold are unaffected."
        }
      }
    }
  }

  * if poisson, check that y is not negative
  if "`poisson'" != "" {
    cap assert `y' >= 0 `ifassert'
    if _rc != 0 {
      dis as error "Error: the outcome variable contains negative values; this is not allowed with poisson."
      error 411
    }
  }


  ****** process other aspects of the command syntax

  * if no elasticity type specified, assume regular non-elasticity estimates; a supplied
  * but unrecognized value is rejected by the domain screen at the end of this program
  * (formerly any unrecognized value collapsed silently to dydx)
  if "`contreatelasticitytype'" == "" local contreatelasticitytype dydx
  if "`semielasticity'" == "" {
    local semielasticity dydx
  }
  else {
    local semielasticity eydx
  }

  * ---- BS-6 (author ruling 2026-09-07): the poisson xb-scale remap --------------------
  * Under poisson WITHOUT poisexpresults the reporting scale is the linear index
  * xb = ln E[y|x], on which margins' PLAIN operators are already the elasticity family:
  *     dydx(xb) = d ln E[y] / dx      = the semi-elasticity  (what semielasticity promises)
  *     dyex(xb) = d ln E[y] / d ln x  = the elasticity       (what eyex promises)
  * The ey* operators instead take a FURTHER log of the prediction - margins computes
  * d ln(xb)/dx, a log of a log - so the shipped number was scale-dependent in y and did
  * not match its own "Semi-elasticities"/"Elasticities" header. Measured on
  * test/data/dta_binary.dta (att, cluster(i)): semielasticity reported .1422311 (se .0494,
  * p .00400) where the plain poisson run reports .1240047 (se .02527, p 9.23e-07), and
  * multiplying ycount by 10 moved the reported "semi-elasticity" to .0399387 - a
  * semi-elasticity of E[y] cannot depend on the units of y. Remap to the operator that
  * DOES deliver the header's quantity; the raw request is carried for the labels.
  * Consequences, all intended: the reported number becomes the plain poisson one, the
  * reads become linear in the coefficients so the cell-collapse ban stops catching them,
  * and no ey* operator reaches margins under poisson any more, so the data-dependent
  * r(459)/missing-with-rc-0 behaviour frozen at the poismarginstext block cannot arise
  * there. poisexpresults is deliberately untouched: its scale is predict(n) = E[y] > 0,
  * where eydx/eyex are the genuinely correct operators.
  * ORDERING (two traps): test eyex/eydx by EQUALITY, never "anything but dydx/dyex" - the
  * contreatelasticitytype() domain screen at the end of this program must still see an
  * unrecognized value. And gate the contreat side on contreat(): contreatelasticitytype()
  * without contreat() has no requires-parent screen and is inert, so a caution there would
  * name a dose that does not exist.
  local semielasraw `semielasticity'
  local cetyperaw `contreatelasticitytype'
  if "`poisson'" == "poisson" & "`poisexpresults'" == "" {
    if "`semielasticity'" == "eydx" local semielasticity dydx
    if "`contreat'" != "" {
      if "`contreatelasticitytype'" == "eydx" local contreatelasticitytype dydx
      if "`contreatelasticitytype'" == "eyex" local contreatelasticitytype dyex
    }
    if "`semielasticity'" != "`semielasraw'" | "`contreatelasticitytype'" != "`cetyperaw'" {
      dis as text "Caution: under poisson without poisexpresults, results are on the linear-index scale xb = ln E[y|x], where a plain marginal effect already IS the semi-elasticity of E[y] and dyex already IS the elasticity. wooldid therefore computes the requested quantity with the xb-scale operator (semielasticity and contreatelasticitytype(eydx) as dydx; contreatelasticitytype(eyex) as dyex) rather than with eydx/eyex, which would take a further log of the linear index and report a number that depends on the units of y."
      dis as text "  The option changes the table header, not the estimate: these are the estimates a plain poisson run reports. Specify poisexpresults to obtain the same quantities from the exponentiated prediction instead."
    }
  }

  * ---- F-collapse gate ----------------------------------------------------------------
  * The collapse hands margins one row per (_i,_t,targetflag) cell. That is exact only
  * while the treatment derivative is constant within the cell, i.e. while no model term
  * interacts targetflag with something that varies inside a cell and the reported
  * quantity is a plain marginal effect on a scale linear in the coefficients.
  *
  * IT IS ON BY DEFAULT (author-ruled 2026-09-01) whenever the request permits it:
  *   fullprecisionslow          -> never collapse
  *   fastcollapse, ineligible   -> ERROR: an explicit request for something that would
  *                                 return a measured-wrong answer must not pass silently
  *   otherwise                  -> collapse iff eligible; SILENTLY skip if not, because
  *                                 an ordinary contreat() run must not start erroring
  * The refusals carry no safety(off) escape: safety(off) means "I accept a risk", never
  * "give me numbers I have been shown are wrong". `fcon' is the realized decision; the
  * worker consumes it and records it in e(wdid_fcol). `fcskip' names the skip reason(s)
  * for the verbose line.
  * subgroups() is NOT on this list (author-ruled 2026-09-01, after a 2026-09-01 refusal):
  * the subgroup reads are exact under the collapse except in a numerically singular VCE
  * direction, where the exact zero they branch on is cancellation residue. The worker
  * screens each subgroup margins call for that (wooldid_fcdeg) and re-issues the call on
  * the full data when it fires - or, under safety(off), proceeds with a caution.
  local fcbad
  * contreat() itself is eligible: every collapsible read is linear in the dose inside a cell
  * (dydx(targetflag) = delta_g + beta_g*x, dydx(contreat) = beta_g, dyex = beta_g*x), so a
  * summed weight and a weighted mean reproduce it exactly. what is refused is the reads that
  * are NOT linear in the dose, and the requests that split a cell on a variable outside the key
  if "`contreat'" != "" {
    if `contreatpoly' > 1                                                       local fcbad `fcbad' contreatpoly(>1)
    if `"`lattice'"' != ""                                                      local fcbad `fcbad' lattice()
    * only makecfxplots reads the bins (from lattice() or cfxlattice()); cfxlattice() on its
    * own is announced as ignored and must not decide the collapse (bug sweep BS-27)
    if "`makecfxplots'" != ""                                                   local fcbad `fcbad' makecfxplots
    * BS-6 (author ruling 2026-09-07): under poisson without poisexpresults these two were
    * remapped to dyex|dydx above, so they never reach this line and the collapse ban lifts
    * on exactly the F032 argument - the remapped reads are linear in the coefficients.
    * Measured on test/suite/data/continuous_control_dense.dta: e(wdid_fcol) 0 -> 1, and the
    * collapsed answer agrees with fullprecisionslow to 9.1e-14, the same figure a plain
    * poisson contreat() run posts. Off the poisson path, and under poisexpresults, nothing
    * was remapped and the ban stands (test_fastcollapse F11b/F11c legs 2, 3 and 5).
    if "`contreatelasticitytype'" == "eyex" | "`contreatelasticitytype'" == "eydx" local fcbad `fcbad' contreatelasticitytype(`contreatelasticitytype')
    * poisson + contreatelasticitytype(dyex) left this list on 2026-09-06 (F032): dyex now
    * takes the analytic predict(xb) form under poisson too, and the collapse reproduces the
    * full-data ame rows to 6.3e-16 there (research record F032.md)
    * the treatment variable inside a user control term (controls(c.dose#c.x) and the like)
    * makes the DOSE derivative vary within a cell on a variable the collapse does not average
    local fcctrl `"`controls' `coarsecohortcontrols' `fe' `ccc_absorb'"'
    foreach tok in "#" "(" ")" "," {
      local fcctrl = subinstr(`"`fcctrl'"', "`tok'", " ", .)
    }
    local fcctrlhit 0
    foreach w of local fcctrl {
      * a Stata variable name cannot contain a dot, so whatever follows the LAST dot is the
      * name and everything before it is operator notation - c.dose, 2.dose, L.dose and the
      * stacked c.L.dose all reduce to dose, where stripping a single prefix left L.dose and
      * missed the gate
      local w2 = substr("`w'", strrpos("`w'", ".") + 1, .)
      if "`w2'" == "" continue
      * controls(), fe() and ccc_absorb() are varlists, so Stata has already unabbreviated
      * them by the time we see them; coarsecohortcontrols() is declared string and has NOT
      * been, so an abbreviated or wildcarded name would otherwise walk past this check and
      * collapse a request that must not be collapsed. expand each token and look for the
      * treatment variable in what it expands to
      cap unab fcw : `w2'
      if _rc == 0 {
        if `: list posof "`contreat'" in fcw' local fcctrlhit 1
      }
      else if "`w2'" == "`contreat'" local fcctrlhit 1
    }
    if `fcctrlhit' local fcbad `fcbad' contreat-in-controls
  }
  if "`interactivecontrols'" != ""          local fcbad `fcbad' interactivecontrols()
  if "`unconditionalse'" != ""              local fcbad `fcbad' unconditionalse
  * BS-6: as the contreat twin above - remapped to dydx under poisson without
  * poisexpresults, so the ban lifts there and stands everywhere else
  if "`semielasticity'" != "dydx"           local fcbad `fcbad' semielasticity
  if "`poisexpresults'" != ""               local fcbad `fcbad' poisexpresults
  if `"`customvce'"' != ""                  local fcbad `fcbad' customvce()
  if "`fastcollapse'" != "" & "`fullprecisionslow'" != "" {
    dis as error "Error: fastcollapse and fullprecisionslow are contradictory - the first asks for the cell collapse and the second suppresses it. Specify at most one. (The collapse is on by default wherever the request permits it, so neither is needed to get it.)"
    error 198
  }
  if "`fastcollapse'" != "" {
    if "`fcbad'" != "" {
      dis as error "Error: fastcollapse is incompatible with: `fcbad'."
      dis as error "fastcollapse evaluates the margins calls on one row per cohort-by-period-by-treatment cell, which is exact only when the reported quantity is linear in the coefficients and, within a cell, linear in the treatment variable. Each option above breaks one of those conditions:"
      dis as error "  interactivecontrols() interacts the treatment dummy with a variable that varies inside a cell, so the effect is no longer cell-constant. contreat() itself is eligible - the collapse carries the cell's weighted mean of the dose - but the reads listed below are not."
      dis as error "  contreatpoly() above 1 makes the effect a polynomial in the dose, which the cell mean does not reproduce; lattice(), makecfxplots and cfxlattice() split a cell on a binning of the dose, which is not part of the cell key."
      dis as error "  contreat-in-controls means the contreat variable also appears inside controls(), coarsecohortcontrols(), fe() or ccc_absorb(); the dose derivative then varies inside the cell on a variable the collapse does not average (measured: the ame rows move by up to 0.6 percent while the att rows stay exact)."
      dis as error "  semielasticity (and the contreat elasticity forms eyex and eydx) divide by the prediction, which re-introduces the within-cell variation that a plain marginal effect cancels. contreatelasticitytype(dydx) and (dyex) are linear in the dose and stay eligible."
      dis as error "  poisexpresults reports on the exponentiated scale, which is not linear in the coefficients (measured: estimates moved in the 7th significant digit)."
      dis as error "  unconditionalse asks margins for a variance built from observation-level influence functions; margins verifies the visible row count against e(N) and errors when it is collapsed."
      dis as error "  customvce() is not verified for fastcollapse rather than known-broken; it is refused pending measurement."
      dis as error "Drop fastcollapse, or drop the conflicting option. (No safety(off) escape: these are measured wrong answers, not risks.) Without an explicit fastcollapse the run proceeds normally on the full data - the collapse simply does not engage."
      error 198
    }
  }
  local fcon 0
  local fcskip `fcbad'
  if "`fullprecisionslow'" == "" & "`fcbad'" == "" local fcon 1

  * do some light processing of the "number of pre-treatment coefs" input - it needs some adjustment when esfixedbaseperiod is used since the period normalized to 0 could fall within your range of requested coefs
  * F-boot A8: preserve the user's raw esprelength for rep-line reconstruction. the
  * normalization below is RE-APPLIED whenever a nested rep re-parses, so embedding the
  * mutated value in a rep command line widens the window once per nesting level
  * (measured: user 3 -> c_local 4 -> rep re-parse 5). rep builders emit this raw copy.
  local esprelengthraw = `esprelength'
  if "`esfixedbaseperiod'" == "esfixedbaseperiod" &   `esprelength' > 0  {
    if abs(`esrelativeto') <= `esprelength' local ++esprelength
    * this is needed b/c we arent counting the reference period coef as a coef
  }

  * process the joint tests suppression flag  - it is on by default unless jointtests is specified
  local suppressjointtests suppressjointtests
  if "jointtests" == "`jointtests'" local suppressjointtests

  * what type of results do we want estimated  - if nothing specified, default estimand to att
  if "`att'" == "" & "" == "`att_it'" & "" == "`att_itime'" & "" == "`att_i'" {
    local att att
  }

  * handle contreatwithin - indicator to use only within variation
  * crash if probably illegitimate contreatcontrols are specified
  * also, set safety to off
  if "`contreatwithin'" == "contreatwithin" & "`contreat'" != "" {
    if strpos("`contreatcontrols'","byTime") != 0  {
      if "`safetyuser'" == "off" dis as text "Caution: contreatwithin is specified in conjunction with the byTime contreatcontrol type."
      if "`safetyuser'" != "off" dis as error "Error: contreatwithin is specified in conjunction with the byTime contreatcontrol type."
      if "`verbose'" == "verbose" dis as text "Typically, this partials out some treatment effects and yields unidentified estimates."
      if "`safetyuser'" != "off" {
        dis as error "To proceed anyway in a new estimation run, manually set safety(off) before proceeding."
        error 198
      }
    }
    * safety was already relaxed above, where the effective value is resolved before any
    * message keys on it (CS-98). The relaxation is required because contreatwithin
    * intentionally makes the base effects irrelevant
  }


  * process cfxplottypes
  * note that the tests here on what cfxplottypes are allowed are kinda weak in that you can request impractical relative periods (relyrs) and such, but better than no checks
  if "`makecfxplots'" != "" {
    local cfxpost 0
    local cfxpre 0
    local cfxrelyrlist
    if "`cfxplottypes'" == "" local cfxplottypes post
    * VAL-G13A: exact-token parse. the old substring matching accepted garbage tokens
    * (postpost counted as post, 2post2 became 2 2), and fractional numeric entries fell
    * through both cap'd window filters silently (rc 0, no plot, no message). keywords are
    * matched exactly per token; numeric tokens must be non-missing integers, loudly.
    foreach item in `cfxplottypes' {
      if "`item'" == "post" {
        local cfxpost 1
        continue
      }
      if "`item'" == "pre" {
        * without esf or pre-periods, pre is dropped (documented precondition in the sthlp)
        if "`esfixedbaseperiod'" == "esfixedbaseperiod" | `esprelength' > 0 local cfxpre 1
        continue
      }
      * remaining entries must be numeric: a stray token would otherwise be evaluated as a
      * variable reference below -- crashing with an unexplained r(111) if no such variable
      * exists, or silently injecting the variable into the plot list if one does (stata
      * does not short-circuit |, so the reference-period comparison always runs)
      cap confirm number `item'
      if _rc != 0 {
        dis as error "Error: cfxplottypes() contains an unrecognized entry: `item'"
        dis as error "Valid cfxplottypes() entries are post, pre, and/or integers specifying relative time periods."
        error 198
      }
      * a leading + is a natural way to write positive relative years (+2, matching -2), but
      * the raw token is later embedded in matrix names, where + is an invalid character;
      * strip it so +2 means 2 (and the post-loop dedup then also collapses +2 with 2)
      if substr("`item'",1,1) == "+" local item = substr("`item'",2,.)
      if `item' >= . | int(`item') != `item' {
        dis as error "Error: cfxplottypes() contains a missing or non-integer entry: `item'"
        dis as error "Valid cfxplottypes() entries are post, pre, and/or integers specifying relative time periods."
        error 198
      }
      if "`esfixedbaseperiod'" == "" | `item' != `esrelativeto' {
        if `item' >= 0 & (`item' <= `espostlength' | "`esfixedbaseperiod'" == "esfixedbaseperiod") local cfxrelyrlist `cfxrelyrlist' `item'
        if `item' < 0 & (abs(`item') <= `esprelength' | "`esfixedbaseperiod'" == "esfixedbaseperiod") local cfxrelyrlist `cfxrelyrlist' `item'
      }
    }
    * a repeated relyr entry would build the same result matrix twice; the final ereturn loop
    * moves (not copies) each matrix, so the duplicate's second pass dies r(111). list uniq
    * preserves first-occurrence order, so duplicate-free input is untouched
    local cfxrelyrlist : list uniq cfxrelyrlist
  }


  * weight-sign screen: the worker's `w' is its own tempvar when no weights are
  * given (identically 1, making the assert vacuous there), so screen the user's
  * weight expression directly instead
  local w `exp'
  * gate on safetyinput, not safety: contreatwithin's internal safety-off overwrite sits
  * upstream and previously skipped this screen entirely - contreatwithin + a negative
  * weight ran rc 0 with no message (verifier s5_v2, 2026-07-26). an explicit user
  * safety(off) still bypasses
  if  "`safetyinput'" != "off" & "`w'" != "" {
    cap assert `w' >= 0 `ifassert'
    if _rc != 0 {
      dis as error  "Error: weights must not be negative."
      error 411
    }
  }

  * VAL-G15A: the negative-subgroup screen moved into the worker (beside the integer
  * screen) so it runs on the estimation-aware domain (ifstatement folds markout;
  * internal ttre < . spans exactly the rows any over(subgroups) subpop can consume).
  * the raw-if-domain check that lived here falsely rejected runs whose only negative
  * values sat on deterministically-discarded rows (e.g. missing outcomes).

  * ---- S5: additional early screens for inputs that previously failed late with ----
  * ---- opaque messages or were silently mis-handled. placed AFTER all pre-existing ----
  * ---- checks so no existing multi-error precedence changes. ----

  * makecfxplots without a continuous treatment was a silent no-op: the whole cfx
  * production block is gated on contreat(), so makecfxplots + cfxlattice() alone
  * returned rc 0 with no plots and no message
  if "`makecfxplots'" != "" & "`contreat'" == "" {
    dis as error "Error: makecfxplots requires a continuous treatment; specify contreat()."
    error 198
  }

  * one-token lattice() was previously unvalidated: a forgotten quantile count surfaced
  * as a bare "variable all not found" r(111), and bad variable contents survived to the
  * regression stage. variable lookup runs first, so a genuine variable that happens to
  * share a lattice-type name still wins, exactly as before
  if "`lattice'" != "" & wordcount("`lattice'") == 1 {
    cap confirm numeric variable `lattice'
    if _rc != 0 {
      local wdvlatone = strlower("`lattice'")
      if inlist("`wdvlatone'","all","basetreated","basetreatedpost","basetreatedpre","basetreatedref") | inlist("`wdvlatone'","bycohort","bytime","bycohortbytime","bytreated","bytreatedbytime","byttre","byttrebytime") {
        dis as error "Error: lattice(`lattice') specifies a lattice type without a quantile count; use, for example, lattice(`lattice' 4)."
        error 198
      }
      dis as error "Error: lattice() must contain either an existing numeric variable or a lattice type followed by a quantile count."
      error 198
    }
    * gate on safetyinput (the user's own safety() setting, captured before the
    * contreatwithin block relaxes `safety'), matching the F13 checks: an internal
    * relaxation must not silently disable a data-integrity screen
    if "`safetyinput'" != "off" {
      cap assert (`lattice' >= 0 & int(`lattice') == `lattice') | `lattice' >= . `ifassert'
      if _rc != 0 {
        dis as error "Error: the lattice() variable must contain only non-negative integers (it is used as a factor variable)."
        error 411
      }
    }
  }

  * cfxlattice() content screens: margins over() accepts only non-negative integers, and
  * previously bad contents died late at the margins/regression stage with generic
  * messages. worse, margins silently DROPS missing-cfxlattice rows and wooldid calls
  * margins quietly, so the omitted-observations note never reached the user (measured:
  * rc 0 with materially different cfx results) - hence the missing-values Caution
  if "`cfxlattice'" != "" & "`makecfxplots'" == "" {
    * inert option: without makecfxplots the variable is never passed to margins, and
    * such calls ran fine before these screens existed - caution, do not error
    dis as text "Caution: cfxlattice() has no effect without makecfxplots and is being ignored."
  }
  if "`cfxlattice'" != "" & "`makecfxplots'" != "" {
    cap confirm numeric variable `cfxlattice'
    if _rc != 0 {
      dis as error "Error: cfxlattice() must specify a numeric variable."
      error 198
    }
    if "`safetyinput'" != "off" {
      * VAL-G12A/G12B: both screens tolerate rows with a missing recorded onset - every
      * executed CFX subpop requires internal ttre != . (never-treated rows feed no plot),
      * so values there can neither crash margins nor silently shift a plotted average.
      * mirrors the ttre_asgiven >= . tolerance of the commensurability screen below
      cap assert (`cfxlattice' >= 0 & int(`cfxlattice') == `cfxlattice') | `cfxlattice' >= . | `ttre_asgiven' >= . `ifassert'
      if _rc != 0 {
        dis as error "Error: the cfxlattice() variable must contain only non-negative integers on observations with a nonmissing treatment onset (margins allows only non-negative integers in over())."
        error 411
      }
      local wdvcfxmiss `cfxlattice' >= . & `ttre_asgiven' < .
      if `"`if'"' != "" local wdvcfxmiss (`wdvcfxmiss') & ((`if') == 1)
      qui count if `wdvcfxmiss'
      if r(N) > 0 dis as text "Caution: the cfxlattice() variable is missing on `r(N)' observation(s) with a nonmissing treatment onset; margins omits those observations when computing the continuous treatment effect plot points, so the plotted averages reflect only the nonmissing-cfxlattice subsample."
    }
  }

  * non-commensurate ttre: fractional t - ttre previously either died late at the
  * regression stage with opaque errors (a bare "invalid syntax", a tempvar-named
  * factor-variable error, or a misleading insufficient-pre-periods r(459)) or - worse -
  * ran to completion with the onset period silently mislabeled as pre-treatment.
  * integer SPACING is not required - biennial or 5-year panels pass - only that onsets
  * sit on the same grid as t. tiered like the F13 checks: offending onsets that fall
  * wholly outside the usable window (dropped as pre-sample or treated as never-treated
  * downstream) downgrade to a caution rather than an error
  if "`safetyinput'" != "off" {
    cap assert (`t' - `ttre_asgiven') == int(`t' - `ttre_asgiven') | `ttre_asgiven' >= . `ifassert'
    if _rc != 0 {
      local wdvtrng `t' < .
      if `"`if'"' != "" local wdvtrng (`wdvtrng') & ((`if') == 1)
      qui gstats sum `t' if `wdvtrng', meanonly
      local wdvtmin = r(min)
      local wdvtlo = `wdvtmin' + `anticipation'
      local wdvthi = r(max) + `oostreatedcontrols' + `anticipation'
      * under timetrends the worker needs a two-period reference block, so it drops treated
      * units with onsets at or below the SECOND observed period (mintp1, per-unit); lift the
      * tier's lower bound to the second-smallest distinct t so fractional onsets in that
      * band are screened by the worker's own drop rule rather than erroring here. exact for
      * any period spacing (mintp1_i >= min{t : t > min(t)}); falls back to the first-period
      * bound if the panel has a single distinct period
      if "`timetrends'" != "" {
        qui gstats sum `t' if (`wdvtrng') & `t' > `wdvtmin', meanonly
        if r(N) > 0 local wdvtlo = r(min) + `anticipation'
      }
      cap assert (`t' - `ttre_asgiven') == int(`t' - `ttre_asgiven') | `ttre_asgiven' >= . | `ttre_asgiven' <= `wdvtlo' | `ttre_asgiven' > `wdvthi' `ifassert'
      if _rc != 0 {
        dis as error "Error: t minus ttre takes non-integer values, so relative-period horizons cannot be constructed. Treatment onset times must sit on the same grid as the time variable; if the time variable itself is on a fractional grid (e.g. year + quarter/4), recode both to consecutive integers."
        error 459
      }
      dis as text "Caution: some treatment onset times are not on the same grid as t (t minus ttre is non-integer), but all such onsets fall outside the usable estimation window and will be dropped or treated as never-treated downstream."
    }
  }

  * esrelativeto() without esfixedbaseperiod was silently ignored (output is identical
  * with and without it)
  if `esrelativeto' != -1 & "`esfixedbaseperiod'" == "" {
    dis as text "Caution: esrelativeto() has no effect without esfixedbaseperiod and is being ignored."
  }

  * ---- F9: fast (no-SE) mode gates ----
  * placed last so no pre-existing check's position (and so no multi-error call's first
  * message) moves. fastnose runs every margins call with nose: point estimates only.
  * everything downstream that needs a variance is therefore unavailable and is refused
  * here rather than silently degraded - margins itself rejects nose+vce(unconditional)
  * (r198), test after a nose post is r321, and contrast(overjoint) under nose returns a
  * MISSING p rather than erroring, which would otherwise be displayed and ereturn'd.
  if "`fastnose'" != "" {
    if "`unconditionalse'" != "" {
      dis as error "Error: fastnose is not compatible with unconditionalse; unconditional standard errors are produced by margins, which fast mode does not ask for."
      error 198
    }
    if "`impvalues'" != "" {
      dis as error "Error: fastnose is not compatible with impvalues; Ibragimov-Muller inference reports standard errors and p-values, which fast mode does not estimate."
      error 198
    }
    * CS-23 (author ruling 2026-09-08): under permute() without permtstat/permanalyticse,
    * fast mode is what wooldid itself selects below, and the adversarial and joint-test
    * surfaces are permutation quantities the shell builds after the replications. Typing
    * the word that is already in force was an error while omitting it ran fine; these two
    * refusals now stand down in exactly that case. Every other fastnose incompatibility
    * (unconditionalse, impvalues, pretrendpower) is unaffected.
    local fnperm = ("`permute'" != "" & "`permtstat'" == "" & "`permanalyticse'" == "")
    if "`adversarial'" != "" & !`fnperm' {
      dis as error "Error: fastnose is not compatible with adversarial; the adversarial pre-trends test needs a variance-covariance matrix, which fast mode does not estimate."
      error 198
    }
    if "`jointtests'" != "" & !`fnperm' {
      dis as error "Error: fastnose is not compatible with jointtests; fast mode reports point estimates only."
      error 198
    }
    if "`pretrendpower'" != "" {
      dis as error "Error: fastnose is not compatible with pretrendpower(); the Roth (2022) diagnostics need the event-study variance-covariance matrix, which fast mode does not estimate."
      error 198
    }
    * belt-and-braces: joint tests are already off whenever jointtests is absent, and the
    * gate above rejects jointtests, so this is redundant today - it is kept explicit so a
    * future relaxation of that gate cannot silently re-enable a joint test under nose
    local suppressjointtests suppressjointtests
  }
  * ---- permute() fast-mode default (author ruling 2026-09-05) ----
  * without permtstat the observed run is a fast (no-SE) run as well: nothing the
  * permutation inference reports consumes the analytic standard errors, so the margins
  * variance work they cost is wasted. permanalyticse restores them (the tables then carry
  * them labelled se_analytic). Placed AFTER the fast-mode refusals above, which key on the
  * placed AFTER the fastnose incompatibility gates on purpose: clearing adversarial here
  * would otherwise mask the documented "fastnose is not compatible with adversarial"
  * refusal, turning it into a caution about the window instead (test_fastnose T4).
  * CS-38 (author rulings 2026-09-08): the defect was that a requested adversarial test
  * could be dropped in SILENCE. The first ruling screened it like its sibling
  * pretrendpower() and refused; measured against the suite that refused 121 of 501
  * configurations, every one carrying the option across a grid rather than asking for the
  * test deliberately, which is how a user setting options from a macro would hit it too.
  * The second ruling is therefore to CAUTION rather than refuse: the user still learns the
  * test was not performed, which is the whole defect, and a decorative option no longer
  * breaks a working script. The caution for a window the DATA fails to realize lives in
  * the worker; this one covers a window that could never have supported the test.
  if "`adversarial'" != "" {
    if `esprelength' <= 0 | `espostlength' <= 0 {
      dis as text "Caution: adversarial pre-trend testing was dropped: the test fits a line through the reported relative-period coefficients and needs at least one pre-treatment and one post-treatment horizon, but this window requests `esprelength' pre and `espostlength' post. Specify esprelength() and espostlength() of at least 1 to run it."
      local adversarial
    }
    else if "`timetrends'" != "" & `esprelength' + `espostlength' + 1 < 3 {
      dis as text "Caution: adversarial pre-trend testing was dropped: with timetrends the test needs at least three event-study coefficients in the requested window (the reference period counts as one), and this window requests `esprelength' pre and `espostlength' post. Widen the window or drop timetrends to run it."
      local adversarial
    }
  }

  * USER's fastnose: under permute() the adversarial and joint-test surfaces are
  * permutation quantities built by the shell after the replications, not analytic ones
  * that fast mode cannot produce, and pretrendpower() is refused by wooldid_permctl.
  * suppressjointtests is deliberately NOT set here: the worker's analytic joint-test
  * blocks are gated on permute() themselves, and the display must still print the
  * permutation joint tests the shell posts under the analytic names
  if "`permute'" != "" & "`permtstat'" == "" & "`permanalyticse'" == "" local fastnose fastnose

  * ---- option-domain screens (issue-file batch VAL-G19/G21A/G21B/G21C/G21E) ----
  * appended after all pre-existing checks so no existing multi-error call's first message
  * moves. each closes a silent-accept path: these inputs previously ran rc 0 with
  * silently different semantics than requested.
  * contreatwithin without contreat() was silently ignored AND disabled the worker's
  * binary-design estimability checks (which gate on the bare flag)
  if "`contreatwithin'" != "" & "`contreat'" == "" {
    dis as error "Error: you cannot request contreatwithin without specifying a continuous treatment variable via contreat()."
    error 198
  }
  * a supplied but unrecognized elasticity type formerly collapsed silently to dydx
  if !inlist("`contreatelasticitytype'", "dydx", "eyex", "eydx", "dyex") {
    dis as error "Error: contreatelasticitytype() must be dydx, eyex, eydx, or dyex."
    error 198
  }
  * syntax's int type rejects fractions/missing but admits 0 and negatives, which
  * formerly ran identically to order 1
  if `contreatpoly' < 1 {
    dis as error "Error: contreatpoly() must be a positive integer."
    error 198
  }
  * tested on the raw user value (safetyinput): contreat+contreatwithin rewrites the
  * working copy to off above, which would mask an invalid user token here
  if "`safetyinput'" != "" & !inlist("`safetyinput'", "on", "off") {
    dis as error "Error: safety() must be on or off."
    error 198
  }
  * SI-9 (author-ruled 2026-08-25): formerly accepted silently with zero effect
  if "`poisexpresults'" != "" & "`poisson'" == "" {
    dis as error "Error: poisexpresults requires poisson to have any effect; specify poisson, or drop poisexpresults."
    error 198
  }
  * under esfixedbaseperiod the worker formerly clamped esrelativeto() values above -1
  * to -1 silently; the reference period must lie strictly before onset
  if "`esfixedbaseperiod'" != "" & `esrelativeto' > -1 {
    dis as error "Error: esrelativeto() must be a negative integer (a pre-onset reference period) when esfixedbaseperiod is specified."
    error 198
  }
  * NFADV-07: altmethod_adversarial is a modifier of adversarial; every gate and work
  * branch keys on the adversarial flag, so the standalone spelling ran rc 0 as a silent
  * no-op. error rather than imply, matching the modifier-without-parent convention
  * (lattice()-without-contreat, contreatwithin-without-contreat)
  if "`altmethod_adversarial'" != "" & "`adversarial'" == "" {
    dis as error "Error: altmethod_adversarial modifies the adversarial pre-trends test; specify the adversarial option as well (or remove altmethod_adversarial)."
    error 198
  }

  * ---- return the parsed options and every derived/normalized local to the worker ----
  * blanket-by-construction: every option the syntax line defines is listed (plus the
  * derived locals below), so a future edit that mutates one more local here cannot
  * silently strand the worker on the stale syntax-parsed value. varlist/if/weight/exp
  * are excluded: the worker's own syntax call already parsed them identically and this
  * subprogram never alters them.
  * CS-2 (consistency sweep 2026-09-06): an unspecified level() means c(level) - the
  * analytic tables already used it (F054) while every engine rewiring fell back to a
  * literal 95, so `set level 90' produced lb_90ci columns holding 95% intervals under
  * bootstrap() / permute(). Resolve it once here; every consumer sees a number
  if `level' == 0 local level = c(level)
  foreach wdv_ret in att att_it att_itime att_i cluster customvce unconditionalse robust impvalues fe coarsecohortcontrols ccc_absorb controls interactivecontrols timetrends subgroups customweightmultiplier esprelength espostlength saveplots poisson poisexpresults oostreatedcontrols anticipation esfixedbaseperiod esrelativeto safety makeplots makecfxplots summarystats histogramcohorteffects regtol jointtests suppresssgprimaryeffects contreat contreatpoly lattice latticeignoreweights ccnoabsorb contreatcontrols verbose cfxplottypes cfxlattice makecfxplotsbysg cfxplotnose semielasticity contreatelasticitytype contreatwithin cleanmatrices emptycellsoverride replace altmethod_adversarial adversarial update oldsyntax pretrendpower fastnose fastcollapse fullprecisionslow fcon fcskip gteffects permute permreps permexact permtstat permtci permstrata permprobs permsave permanalyticse permabsp permside permrawp permstepdown simultaneousbands y i t ttre_asgiven safetyinput suppressjointtests latticenum latticetype esposttrip cfxpost cfxpre cfxrelyrlist bootstrap reps ties level dots rseed bootalpha bootdf bootdisplay boottstat bootsave esprelengthraw semielasraw cetyperaw {
    c_local `wdv_ret' `"``wdv_ret''"'
  }

end

* ============================================================================
* F4 helper subprograms (kept out of the worker to conserve its size headroom)
* ============================================================================

* warn when an effective onset lands on a period not observed in the data - on panels
* whose t is not unit-spaced (e.g. biennial), an off-grid anticipation shift can merely
* relabel horizons (pooled reference), move the fixed reference period in unintended
* ways (esfixedbaseperiod), or make estimation fail
cap prog drop wooldid_anticgrid
prog def wooldid_anticgrid
version 16.0
  syntax , tvar(string) ttrevar(string) ifv(string)
  qui glevelsof `ttrevar' if `ifv' == 1 & `ttrevar' != ., local(anticonsets)
  local anticoffgrid
  foreach anticv of local anticonsets {
    qui count if `tvar' == `anticv' & `ifv' == 1
    if r(N) == 0 local anticoffgrid `anticoffgrid' `anticv'
  }
  if "`anticoffgrid'" != "" {
    dis as text "Caution: anticipation() is expressed in the units of `tvar', and the effective onset value(s) `anticoffgrid' fall on period(s) not observed in the data. Depending on the specification, this can merely relabel horizons, move the reference period in unintended ways, or cause estimation to fail. Verify that the anticipation length matches the spacing of `tvar'."
  }
end

* under esfixedbaseperiod, report units dropped only BECAUSE the anticipation shift
* moved their reference period onto a period they are not observed in (called after the
* reference-period screen zeroes ifvar, before ttre is nulled for the dropped units)
cap prog drop wooldid_anticesf
prog def wooldid_anticesf
version 16.0
  syntax , ivar(string) tvar(string) rawttre(string) refobs(string) esrel(string) antic(string) ttrevar(string) ifv(string)
  * the raw-reference check must use the SAME domain the reference screen itself used:
  * rows in-sample before this screen ran. the screen has already zeroed ifv for its own
  * victims, so that domain = still-in-sample rows plus the just-zeroed victims. without
  * this restriction, a raw reference observed only on rows an EARLIER screen (e.g. a
  * user if()) had excluded would mis-attribute an unrelated drop to the shift
  tempvar rawref rawrefobs
  qui gen byte `rawref' = `tvar' == `rawttre' + `esrel' if `ifv' == 1 | (`refobs' == 0 & `ttrevar' != .)
  qui gegen `rawrefobs' = max(`rawref'), by(`ivar')
  qui count if `refobs' == 0 & `ttrevar' != . & `rawrefobs' == 1
  if r(N) > 0 {
    local anticdropobs = r(N)
    tempvar antictag
    qui gegen `antictag' = tag(`ivar') if `refobs' == 0 & `ttrevar' != . & `rawrefobs' == 1
    qui count if `antictag' == 1
    dis as text "Caution: anticipation(`antic') places the effective reference period of `r(N)' unit(s) on a period unusable for them (not observed, or dropped for lack of a comparison unit); those unit(s) and their `anticdropobs' observations are dropped (their reference period would have been usable without the anticipation shift)."
    qui drop `antictag'
  }
  qui drop `rawref' `rawrefobs'
end

* BS-37 (bug sweep 2026-09-06, author ruling 2026-09-07): wooldid keeps the rows where
* the if() expression equals exactly 1; Stata's own if-qualifier keeps every nonzero
* value, missing included. The rule is deliberate - it is the reason a bare flag
* variable cannot drag its missing and out-of-range rows into the sample - and every
* forwarded copy of the expression uses it, so the estimation sample, the permutation
* design and the bootstrap frame agree. The defect was the silence: the rows the rule
* alone removed drew no message and were counted nowhere. Name them, with the count,
* and hand the count back so the dropped-observation total can include it. A
* subprogram because wooldid_worker sits at Stata's per-program compiled-size ceiling.
* `ifv' is the CALLER's tempvar, passed by name (a tempvar made here would die on
* return); the count comes back by c_local into the caller's `wdvifexcl'
cap prog drop wooldid_ifrule
prog def wooldid_ifrule
version 16.0
  syntax , ifv(string)
  * missing is excluded from the count on purpose: excluding missing-valued rows is the
  * POINT of the rule, not an accident of it, and the help now states it under Syntax
  qui count if `ifv' != 0 & `ifv' != 1 & `ifv' < .
  local wdvifn = r(N)
  if `wdvifn' > 0 dis as text "Caution: the if() expression is nonzero but not exactly 1 on `wdvifn' observation(s). wooldid keeps only the rows where it equals exactly 1, where Stata's own if-qualifier keeps any nonzero value (missing included), so those observations are excluded here and are counted in the dropped-observation total. Write the restriction as a comparison - for example, if <var> != 0 - for Stata's usual behaviour."
  c_local wdvifexcl `wdvifn'
end

* ---------------------------------------------------------------------------
* F-collapse helper (on by default; fastcollapse / fullprecisionslow). Replaces the
* estimation sample with ONE row per (_i, _t, targetflag [, subgroup]) cell, carrying that
* cell's SUMMED weight and keeping the cell's FIRST row.
*
* WHY THE RETAINED ROW'S OTHER VALUES DO NOT MATTER. margins forms dydx(targetflag) as
* xb(tf=1) - xb(tf=0). Every regressor that does not interact with targetflag contributes
* identically to both terms and cancels, so the derivative is the cell coefficient and is
* constant within the cell. That is exactly the condition wooldid_validate's F-collapse
* gate enforces; outside it (interactivecontrols(), a non-dydx operator, the exponentiated
* poisson scale) the derivative varies observation by observation and the retained row would
* BE the estimand rather than an irrelevant carrier - which is why those combinations skip
* the collapse by default and are refused when it is asked for explicitly, rather than
* cautioned. contreat() is inside the condition rather than outside it: the retained row
* carries the cell's weighted mean of the dose (below), which reproduces every read that is
* linear in the dose.
*
* The cancellation is algebraic, not bitwise: (S+c1)-(S+c0) keeps a residue of order
* ULP(|xb|), where S is the retained row's shared partial sum. Results therefore agree
* with the uncollapsed path to machine precision but are not bit-identical - measured to
* move no more than permuting the input rows moves the uncollapsed path itself.
*
* preserve/restore are NOT here on purpose: preserve is scoped to the program that issues
* it, so a preserve here would unwind the instant this program returned and every margins
* call would silently run on full data. Both live inline in the worker body.
cap prog drop wooldid_fcol
prog def wooldid_fcol, rclass
  syntax varname(numeric) , KEY(varlist) [WVARS(string) CONTreat(varlist numeric)]
  * rows in memory become exactly the estimation sample, so the margins calls' noesample
  * is equivalent to their usual e(sample) restriction, and any user if() is already folded in
  qui keep if `varlist' == 1
  * only the weight variables for the requested estimands exist
  local wv
  foreach v of local wvars {
    cap confirm variable `v'
    if _rc == 0 local wv `wv' `v'
  }
  tempvar ord pick
  qui gen long `ord' = _n
  * with contreat(), the surviving row carries the cell's estimand-weighted mean of the
  * continuous treatment. every estimand weight is w x customweightmultiplier divided by a
  * quantity that is constant inside an (i,t) cell (w_per_it, w_per_i, t_per_i - see the weight
  * construction), and `i' is in the key here, so the first weight variable that exists serves
  * as the analytic weight and the four estimands give the same mean.
  * the variable is recast to double first: a float store of the mean injects ~2.4e-8.
  * the mean is formed as sum(w*x)/sum(w) and written onto EVERY row of the cell, not by
  * gcollapse (mean) [aw=w] merge: that form writes MISSING onto zero-weight rows, and the
  * surviving row is the cell's first IN-SAMPLE row whether or not its weight is zero. a zero
  * sampling weight cannot produce that (Stata drops those rows from e(sample)), but a
  * customweightmultiplier() of exactly 0 can, and did (measured on the earlier form: a
  * pre-trend horizon moving p .167 -> .0022). the `!= 0' guard rather than `> 0' costs nothing
  * and keeps a signed customweightmultiplier() under safety(off) exact wherever the cell's
  * weights do not cancel to exactly zero. two cases are left alone deliberately: a cell whose
  * weights are ALL zero keeps its raw values and carries no weight in any read of either arm;
  * and a cell whose non-zero signed weights cancel to exactly zero enters every read with
  * weight zero, losing the dose level it carried (measured at 2.7e-4 relative on the main att
  * family) - unrepresentable in a (summed weight, weighted mean) pair, and reachable only
  * through a signed customweightmultiplier() under safety(off), which is caveat emptor.
  local cmv
  if "`contreat'" != "" {
    local wsrc
    if "`wv'" != "" local wsrc `: word 1 of `wv''
    * the per-variable temporaries are indexed by POSITION, not by the variable's name:
    * a name-built macro (cm_<varname>) would exceed Stata's 31-character local-name limit
    * for a contreat variable of 29 characters or more.
    * each variable carries its OWN weight sum, restricted to the rows where that variable
    * is non-missing. with a non-missing dose - which the caller guarantees, by zeroing the
    * if-marker on any row with a missing contreat value - this is the same denominator for
    * every variable and the same number as a single shared sum, bit for bit. it matters
    * only if a missing value ever reaches here: a shared denominator would divide a
    * numerator taken over the non-missing rows by a weight total taken over all of them,
    * and a cell whose values were ALL missing would sum to 0 over a positive weight total
    * and have its missing dose silently rewritten to 0. with the per-variable form the
    * partial case gives the correct weighted mean over the rows that have a value, and the
    * all-missing case leaves the weight total at 0, which the guard below declines.
    * note the scope: this makes the MEAN missing-aware, not the worker. the cell's summed
    * estimand weight, built further down with the row pick, still totals every row that
    * survived the e(sample) keep, so a row with a missing dose would contribute weight
    * without contributing a dose. that pairing is unreachable for the same reason the rest
    * of this is defensive - such a row is not in e(sample) - and is left alone deliberately
    * so the weight arithmetic stays bit-identical to the binary case
    local ck 0
    foreach cv of local contreat {
      local ctype : type `cv'
      if "`ctype'" != "double" qui recast double `cv'
      local ++ck
      tempvar cm`ck' cw`ck'
      if "`wsrc'" != "" {
        qui gen double `cm`ck'' = `cv' * `wsrc' if !missing(`cv')
        qui gen double `cw`ck'' = `wsrc'        if !missing(`cv')
      }
      else {
        qui gen double `cm`ck'' = `cv' if !missing(`cv')
        qui gen double `cw`ck'' = 1    if !missing(`cv')
      }
      local cmv `cmv' `cm`ck'' `cw`ck''
    }
    qui gcollapse (sum) `cmv', by(`key') merge replace
    local ck 0
    foreach cv of local contreat {
      local ++ck
      qui replace `cv' = `cm`ck'' / `cw`ck'' if `cw`ck'' != 0 & `cw`ck'' < .
    }
    qui drop `cmv'
  }
  * `replace' is required because each weight target is also its own source (the cell
  * sum is written back over the row-level weight in place)
  if "`wv'" != "" qui gcollapse (sum) `wv' (min) `pick' = `ord', by(`key') merge replace
  else            qui gcollapse (min) `pick' = `ord', by(`key') merge
  qui keep if `ord' == `pick'
  qui drop `ord' `pick'
  * deterministic row order: a tag-and-keep build returns panel order where a full
  * collapse returns key order, and row order alone perturbs a near-zero coefficient by
  * more than the collapse does - so the order is pinned rather than left to chance
  qui sort `key'
  return scalar ncell = _N
end

* ============================================================================
* E1 subprogram: build the "every subgroup margin is equal" constraint chain for
* -test- out of the coefficient names Stata itself posted in e(b) after
*   margins ..., over(<sg>) pwcompare(effects) <marginsdof> post
* The names are READ BACK (colnames/coleq), never constructed from subgroup
* values, so they are valid by construction whatever the level values, level
* labels or variable name length; omitted (o.) columns are skipped. -half-
* selects the binary layout, where the posted e(b) carries the base-outcome half
* followed by the 1.<targetflag> half and only the second half is wanted;
* without -half- (continuous treatment) every column is used.
* Returns: sgchain (the chain, empty if fewer than 2 usable columns)
*          sgchainn (number of usable columns; the caller must require > 1,
*          since a one-term chain would test that margin against 0 instead)
* A subprogram because wooldid_worker sits near Stata's per-program compiled-size
* ceiling and this is used at four call sites.
* ============================================================================
cap prog drop wooldid_sgchain
prog def wooldid_sgchain
version 16.0
  syntax [, HALF]
  tempname sgcb
  matrix `sgcb' = e(b)
  local sgck = colsof(`sgcb')
  local sgccn : colnames `sgcb'
  local sgcce : coleq `sgcb'
  local sgcj0 = 1
  if "`half'" == "half" local sgcj0 = floor(`sgck'/2) + 1
  local sgcchain
  local sgcn 0
  forval sgcj = `sgcj0'/`sgck' {
    local sgcnm : word `sgcj' of `sgccn'
    local sgceq : word `sgcj' of `sgcce'
    if strmatch("`sgcnm'","*o.*") continue
    local sgcn = `sgcn' + 1
    if `sgcn' == 1 local sgcchain [`sgceq']`sgcnm'
    else local sgcchain `sgcchain' = [`sgceq']`sgcnm'
  }
  if `sgcn' < 2 local sgcchain
  c_local sgchain `"`sgcchain'"'
  c_local sgchainn `sgcn'
end

* F-collapse degeneracy screen for the subgroup reads (author-ruled 2026-09-01). The
* subgroup section is the one place wooldid BRANCHES on a margins result being exactly
* zero: the omission marker margins writes into a returned name when b and V are both
* zero, the rank -test- assigns the subgroup-equality chain, a missing pairwise se. In a
* numerically singular direction of the VCE that zero is cancellation residue of order
* 1e-16, and the collapse - which changes the summation order - can land it on either
* side, silently changing which horizons and joint tests are reported (measured on 8 of
* 260 subgroup configs, every one in such a direction; none elsewhere). So each subgroup
* margins call is screened on its OWN r() the moment it returns, before anything reads
* it. A flag under the default safety re-issues the call, and every call after it, on the
* full data (the worker does the restore inline: preserve is scoped to the program that
* issued it, and a posted call is un-posted first); under safety(off) the run proceeds
* on the collapsed rows with a caution, as ruled. Reads r() only through Mata, which
* never disturbs it. c_locals: fcrerun (1 = re-issue), fcguard (fallback / caution, for
* e(wdid_fcguard)). Kept out of the worker for compiled-size headroom.
* `joint' marks the overjoint pre-trend calls: when margins returns no r(p) there (a
* subgroup whose pre-treatment window holds only the reference period - an early cohort
* used as a subgroup - has one structurally zero level: b = 0, V = 0, no test in either
* arm) nothing is reported from the call, so there is nothing for the collapse to move
* and the zero-variance predicate would only cost the speedup (2026-09-01 probe, A7).
cap prog drop wooldid_fcdeg
prog def wooldid_fcdeg
version 16.0
  syntax , [SAFety(string) GUARD(string) LEVels VERbose JOINT]
  c_local fcrerun 0
  if "`joint'" != "" {
    cap confirm matrix r(p)
    if _rc exit
  }
  * the tolerance and the `levels' switch are both argued in wooldid_fcdegx's own Mata
  * header and only cited here. 1e-9: residue sits near 1e-16 and the smallest real subgroup
  * contrast measured in the assessment was 3.5e-5, so there are about seven orders of
  * headroom either side. `levels' suppresses the adjacent-difference contrast test, because
  * event-study LEVELS are contrasted by the adversarial line test and by whatever consumes
  * e(V_es_*), both of which make their own rank decision at 1e-12 of the largest eigenvalue
  * - three orders tighter - and a 1e-9 contrast test there fires on every healthy run whose
  * cluster count is below its horizon count (measured 2026-09-06 at 4, 5 and 6 clusters
  * against 8 horizons; BS-29). Anyone moving the number must move it knowing it is paired
  * with that other tolerance, not chosen alone
  mata: wooldid_fcdegx(1e-9, "`levels'" != "")
  if "`fcdeg'" != "1" exit
  if "`safety'" != "off" {
    c_local fcrerun 1
    c_local fcguard fallback
    if "`verbose'" == "verbose" dis as text "Wooldid found `fcdegwhy' at a subgroup margins call, so it is setting the cell collapse aside and evaluating the remaining subgroup results on the full data."
    exit
  }
  c_local fcguard caution
  if "`guard'" == "caution" exit
  dis as text "Caution: wooldid found `fcdegwhy' at a subgroup margins call. Wooldid treats that as a numerically singular direction of the VCE, where an exact zero margins returns can be cancellation residue rather than an estimate, so which subgroup horizons and joint tests are reported can depend on the summation order - which the default cell collapse (see fastcollapse in the help file) changes. Under the default safety setting wooldid would evaluate the subgroup results on the full data from this point; with safety(off) it proceeds on the collapsed rows as requested. Specify fullprecisionslow to evaluate everything on the full data."
end



*** F6: parse-time screens for permutation inference (author's brief: cluster nesting
*** in treatment-onset levels; refuse impvalues/unconditionalse/contreatwithin; reps
*** xor exact; strata subset checks). kept out of wooldid_validate for size headroom.
cap prog drop wooldid_permctl
prog def wooldid_permctl
version 16.0
  syntax , [permute(string) permreps(integer 0) PERMEXact clus(varlist) ttrevar(varname) permstrata(varname numeric) permprobs(string asis) impv(string) uncond(string) cw(string) os(string) ts(string) tc(string) fn(string) cc(string asis) ifr(string asis) ptp(string asis) ssg(string) sg(varname numeric) ivar(varname) tvar(varname) si(string) psave(string) level(real 0) absp(string) side(string) rawp(string) replace]

  local permmethod `permute'
  if "`permmethod'" != "shuffle" & "`permmethod'" != "probabilities" {
    dis as error "Error: permute() must be permute(shuffle) or permute(probabilities). Strata and assignment probabilities have their own options: permstrata() and permprobs()."
    error 198
  }
  if "`impv'" != "" {
    dis as error "Error: permute() is not compatible with impvalues."
    error 198
  }
  * F052/F064 (2026-09-06): the p-value conventions are mutually exclusive
  if "`side'" != "" & !inlist("`side'", "left", "right") {
    dis as error "Error: permside() must be permside(left) - the share of replications at or below the observed statistic - or permside(right) - the share at or above it."
    error 198
  }
  if ("`side'" != "") + ("`rawp'" != "") + ("`absp'" != "") > 1 {
    dis as error "Error: permside(), permrawp and permabsp are alternative p-value conventions; specify at most one."
    error 198
  }
  if "`uncond'" != "" {
    dis as error "Error: permute() is not compatible with unconditionalse."
    error 198
  }
  if "`cw'" != "" {
    dis as error "Error: permute() is not compatible with contreatwithin: randomization inference reassigns onset timing across clusters, which requires a comparison group."
    error 198
  }
  * NFPERM-15: untreatedZero's model formulation requires the dose to be zero outside
  * treated cohort-periods; replications reassign onset while the dose path stays frozen,
  * so every rep would violate that structural premise and silently contaminate the null
  * distribution. plain equality suffices: validate already forces untreatedZero to be
  * the sole contreatcontrols token. contreat in general stays ALLOWED (adjudicated:
  * the null is timing-as-if-random given dose)
  if `"`cc'"' == "untreatedZero" {
    dis as error "Error: permute() is not compatible with contreatcontrols(untreatedZero): replications reassign treatment onset while the dose path stays fixed, violating untreatedZero's requirement that the dose be zero outside treated cohort-periods."
    error 198
  }
  if "`os'" == "1" {
    dis as error "Error: permute() is not compatible with oldsyntax."
    error 198
  }
  * chunk-1 refusals (author ruling 2026-09-05). pretrendpower: the Roth diagnostics are
  * computed on an event-study covariance, and randomization inference produces none - the
  * author allowed the refusal rather than a permutation-covariance substitute.
  * suppresssgprimaryeffects: mirrors bootctl - the suppressed subgroup entries are what
  * the replication vectors need in order to permute the subgroup statistics
  if `"`ptp'"' != "" {
    dis as error "Error: permute() is not compatible with pretrendpower(): the Roth (2022) diagnostics are computed on an event-study variance-covariance matrix, and randomization inference produces none (its p-values and intervals come from the permutation distribution). Run pretrendpower() in a separate call without permute()."
    error 198
  }
  if "`ssg'" != "" {
    dis as error "Error: permute() may not be combined with suppresssgprimaryeffects; the suppressed subgroup entries are exactly what the replication vectors need in order to permute the subgroup statistics."
    error 198
  }
  if `level' != 0 & (`level' < 10 | `level' > 99.99) {
    dis as error "Error: level() must be between 10 and 99.99."
    error 198
  }
  * permsave(): the name rules are bootsave()'s (plain file name, cwd, replace to
  * overwrite), checked here BEFORE the replication loop runs rather than after it
  if "`psave'" != "" {
    local wdpfn `psave'
    if lower(substr("`wdpfn'", -4, 4)) == ".dta" local wdpfn = substr("`wdpfn'", 1, strlen("`wdpfn'") - 4)
    if strpos("`wdpfn'", "/") | strpos("`wdpfn'", "\") | "`wdpfn'" == "" {
      dis as error "Error: permsave() takes a plain file name (no directory part); the file is written to the current working directory."
      error 198
    }
    cap confirm file "`wdpfn'.dta"
    if _rc == 0 & "`replace'" == "" {
      dis as error "Error: permsave(`wdpfn') - `wdpfn'.dta already exists in the working directory; add the replace option to overwrite it."
      error 602
    }
  }
  if ("`ts'" != "" | "`tc'" != "") & "`fn'" != "" {
    dis as error "Error: permtstat/permtci randomize on t-statistics, which fastnose does not compute."
    error 198
  }
  if `permreps' < 0 {
    dis as error "Error: permreps() must be a positive integer."
    error 198
  }
  if `permreps' > 0 & `permreps' < 20 {
    dis as error "Error: permreps() must be at least 20: the engine refuses inference on fewer than 20 completed replications, so a smaller request can only ever burn the run."
    error 198
  }
  * F036 (quality review 2026-09-05, author ruling 2026-09-06: caution, floor stays 20). The
  * sampled-draw p-value has floor 2/(permreps + 1) under the default equal-tailed
  * convention; 1/(permreps + 1) under permabsp and permside(); no floor under permrawp,
  * whose plain c/n share can be 0 (bug sweep BS-20: the caution follows the convention in
  * force, the five-way code wooldid_permwire builds), so below 39 replications a 5% rejection
  * is unreachable by construction and every row prints the floor as if it were evidence
  if `permreps' > 0 & "`permexact'" == "" {
    * the p-value convention code, byte-identical to the `wdpabs' expression in
    * wooldid_permwire and duplicated here only to pick the floor this caution prints:
    * 0 equal-tailed (default), 1 absolute add-one (permabsp), 2 raw c/n share without
    * add-one (permrawp), 3 permside(left), 4 permside(right). KEEP THE TWO IN SYNC - if
    * they drift, this caution announces a smallest attainable p-value the engine does not
    * use, i.e. a wrong number told to the user about their own resolution limit. The
    * cond() ladder's precedence is never exercised: the mutual-exclusivity screen above
    * guarantees at most one of permside()/permrawp/permabsp is set, so the order is an
    * artefact of the nesting and not a ruling
    local wdpcv = cond("`rawp'" != "", 2, cond("`side'" == "left", 3, cond("`side'" == "right", 4, ("`absp'" != ""))))
    local wdp_floor = cond(`wdpcv' == 0, 2/(`permreps' + 1), cond(`wdpcv' == 2, 0, 1/(`permreps' + 1)))
    if `wdp_floor' > .05 dis as text "Caution: with permreps(`permreps') the smallest attainable per-coefficient permutation p-value is " %5.3f `wdp_floor' " (" cond(`wdpcv' == 0, "2", "1") "/(permreps + 1)), above .05, so no coefficient p-value can reject at the 5% level whatever the data say, and a p-value printed at that floor is the resolution limit rather than evidence. Use at least " cond("`absp'" != "", "20", "40") " replications for a 5% test (roughly 100 for .02, 200 for .01 and 1000 for .002). The joint, equality and step-down tests are one-sided add-one shares with floor 1/(permreps + 1)."
  }
  if "`permmethod'" == "probabilities" {
    cap confirm numeric variable `clus'
    if _rc != 0 {
      dis as error "Error: permute(probabilities) requires a numeric cluster variable (the permprobs() grammar lists clusters by numeric value)."
      error 198
    }
    if `"`permprobs'"' == "" {
      dis as error "Error: permute(probabilities) requires permprobs()."
      error 198
    }
    if "`permexact'" != "" {
      dis as error "Error: permute(probabilities) supports permreps() only, not permexact."
      error 198
    }
    if `permreps' <= 0 {
      dis as error "Error: permute(probabilities) requires permreps()."
      error 198
    }
  }
  if "`permmethod'" == "shuffle" {
    if `"`permprobs'"' != "" {
      dis as error "Error: permprobs() is only for permute(probabilities)."
      error 198
    }
    if "`permexact'" != "" & `permreps' > 0 {
      dis as error "Error: specify permreps() or permexact, not both."
      error 198
    }
    if "`permexact'" == "" & `permreps' <= 0 {
      dis as error "Error: permute(shuffle) requires permreps() or permexact."
      error 198
    }
  }
  if wordcount("`clus'") != 1 {
    dis as error "Error: permute() requires cluster() with a single variable - the unit whose treatment-onset assignment is randomized."
    error 198
  }

  * sample for the data checks: the user's if-restriction, rows with a cluster value
  tempvar insmp
  * NFPERM-03: exact-one rule ((exp)==1), matching estimation's if-semantics (8d7bd01):
  * native truth admitted rows with if-values like 2 or missing into the randomization
  * design that never enter estimation (false r459 nesting errors on bare-flag ifs;
  * design-population drift). byte-parallel with permrun's smp below
  if `"`ifr'"' != "" qui gen byte `insmp' = ((`ifr') == 1) & !missing(`clus')
  else qui gen byte `insmp' = !missing(`clus')
  qui replace `insmp' = 0 if `insmp' >= .

  * the level accounting treats all missing codes as the single never-treated group,
  * but enumeration/probability matching distinguishes extended missings; require the
  * plain-missing convention outright
  qui count if `insmp' & `ttrevar' > .
  if r(N) > 0 {
    dis as error "Error: permute() requires never-treated observations coded as plain missing (.) in the timing variable; `r(N)' in-sample observation(s) carry extended missing codes (.a, .b, ...)."
    error 459
  }

  * nesting (author's brief): each cluster must sit inside exactly one treatment-onset
  * level - onset constant within cluster, and treated/never-treated rows never mixed
  tempvar tmin tmax ntr nnv
  qui gegen double `tmin' = min(`ttrevar') if `insmp', by(`clus')
  qui gegen double `tmax' = max(`ttrevar') if `insmp', by(`clus')
  qui gegen double `ntr' = total(!missing(`ttrevar')) if `insmp', by(`clus')
  qui gegen double `nnv' = total(missing(`ttrevar')) if `insmp', by(`clus')
  qui count if `insmp' & ((`ntr' > 0 & `nnv' > 0) | (`ntr' > 0 & `tmin' != `tmax'))
  if r(N) > 0 {
    dis as error "Error: permute() requires clusters to be subsets of treatment-onset groups (levels of the timing variable, with never-treated its own group), but `r(N)' in-sample observation(s) sit in clusters that mix onset values or mix treated and never-treated rows."
    error 459
  }

  * each replication draws ONE onset per cluster and writes it to every observation of
  * the cluster, so a unit (i) straddling two clusters would carry two onsets in a
  * replication - undefined for the estimator (every replication fails, or under
  * permexact the enumeration is over the wrong design). The cluster() variable must
  * nest i: clustering at i or coarser passes, finer is refused (the bootstrap's re-id
  * leg has the same requirement). Added 2026-09-05 (review finding econ_perm F4)
  * F049 (quality review 2026-09-05, author ruling 2026-09-06): the studentized statistic is
  * MacKinnon-Webb (2020)'s remedy for atypical treated clusters, but it presumes a
  * cluster-robust se that sees each treated cluster's own variance. With cluster() at i
  * and one observation per i-period the treated cohort-period cells are absorbed by their
  * own terms, so the se cannot see them and permtstat does not restore size (measured:
  * 7 of 30 null datasets rejected at 5% with treated-cluster error sd 3 vs 1, against 6
  * of 30 for the raw coefficient). Say so; the run proceeds
  if "`ts'" != "" & "`clus'" == "`ivar'" dis as text "Caution: permtstat with cluster() at the i level. The studentized permutation statistic protects against treated clusters being atypical in variance only when the cluster-robust standard error is informative about them - many observations per treated cohort-period cell, with cluster() coarser than the cell. With one observation per i-period the treated cells are absorbed by their own cohort-period terms, the standard error cannot see their variance, and permtstat is not protected there either (MacKinnon and Webb 2020); the test presumes exchangeable clusters in both forms."
  * P7 accident 1 (author ruling 2026-09-08; research_2026-09-08/engine_parity.md): permute()
  * gets its OWN exchangeability caution rather than the bootstrap's few-cluster ones ported
  * across. The failure modes differ in KIND. The pairs cluster bootstrap fails when clusters are
  * not exchangeable as SAMPLING draws. A randomization test fails on two other axes: onset not
  * being exchangeable in ASSIGNMENT, which is the design assumption and cannot be diagnosed
  * here; and, diagnosably, too few ever-treated clusters, which makes the set of distinct
  * assignments small - a limit no replication count can lift, because the draws repeat.
  if "`ts'" == "" & "`clus'" == "`ivar'" dis as text "Caution: cluster() is at the i level and the permuted statistic is the coefficient itself. Randomization inference presumes clusters are exchangeable in ASSIGNMENT; it does not presume they are alike in variance, and the coefficient statistic carries no protection against a treated cluster whose errors are unusually dispersed. permtstat is MacKinnon and Webb's (2020) device for that, with the proviso that it protects only when the cluster-robust standard error can see the treated clusters' own variance - which at cluster() = i with one observation per i-period it cannot."
  tempvar wdpttag
  qui gegen byte `wdpttag' = tag(`clus') if `insmp' & `ntr' > 0
  qui count if `wdpttag' == 1
  local wdpntc = r(N)
  cap drop `wdpttag'
  if `wdpntc' > 0 & `wdpntc' < 6 dis as text "Caution: only `wdpntc' of the resampling clusters are ever treated, so the number of DISTINCT onset assignments is small. Sampled draws will repeat, and the reference set - hence the set of p-values the run can produce - is coarse however large permreps() is; the limit is the design, not the replication count. permexact enumerates the assignments and reports the exact floor instead of sampling them. This is the randomization-inference analogue of the bootstrap's few-treated-cluster caution, not the same condition: that one is about resampling variance, this one about how many assignments exist."
  if "`clus'" != "`ivar'" {
    tempvar cgrp cmin cmax
    qui gegen double `cgrp' = group(`clus') if `insmp'
    qui gegen double `cmin' = min(`cgrp') if `insmp', by(`ivar')
    qui gegen double `cmax' = max(`cgrp') if `insmp', by(`ivar')
    qui count if `insmp' & `cmin' != `cmax'
    if r(N) > 0 {
      dis as error "Error: permute() requires the cluster() variable to nest i (every i wholly inside one cluster): each replication assigns one onset per cluster to all of its observations, which is undefined when a unit straddles clusters (`r(N)' in-sample observation(s) do). Cluster at or above the i level."
      error 459
    }
  }

  * subgroups() under permute() (chunk 1, 2026-09-05): every cluster can become treated in
  * a replication, so the worker's subgroup screens - which run on the OBSERVED treated
  * units only - must hold on the whole design population, or a replication fails by a
  * rule that is correlated with the assignment it drew, which silently biases the
  * reference distribution. The four screens below are the worker's (VAL-G15A negative
  * values, the integer screen, F13 cell containment, F13 unit splitting) with their
  * onset clauses dropped; the sample is the design sample (the user's if and a
  * nonmissing cluster value), slightly wider than the estimation sample
  if "`sg'" != "" {
    cap assert `sg' >= 0 | `sg' >= . if `insmp'
    if _rc != 0 {
      dis as error "Error: permute() with subgroups() requires the subgroup variable to be non-negative (or missing) on every in-sample row, not only on treated units' rows: any cluster can be assigned treatment in a replication, and the replications refuse negative subgroup values on treated rows."
      error 411
    }
    cap assert `sg' == int(`sg') if `insmp'
    if _rc != 0 {
      dis as error "Error: permute() with subgroups() requires integer-valued (or missing) subgroup values on every in-sample row, because any cluster can be assigned treatment in a replication."
      error 411
    }
    tempvar sgmin sgmax
    qui gegen double `sgmin' = min(cond(`sg' >= ., ., `sg')) if `insmp', by(`ivar' `tvar')
    qui gegen double `sgmax' = max(cond(`sg' >= ., ., `sg')) if `insmp', by(`ivar' `tvar')
    cap assert `sgmin' == `sgmax' if `insmp'
    if _rc != 0 & "`si'" != "off" {
      dis as error "Error: permute() with subgroups(): the subgroup variable `sg' varies within `ivar'-by-`tvar' cells of the sample (over ALL units, because any cluster can be assigned treatment in a replication). Define `ivar' at the finer level over which the subgroup varies, or drop subgroups(). To ignore this error, specify safety(off)."
      error 459
    }
    qui drop `sgmin' `sgmax'
    qui gegen double `sgmin' = min(cond(`sg' >= ., ., `sg')) if `insmp', by(`ivar')
    qui gegen double `sgmax' = max(cond(`sg' >= ., ., `sg')) if `insmp', by(`ivar')
    cap assert `sgmin' == `sgmax' if `insmp'
    if _rc != 0 {
      if "`si'" != "off" {
        dis as error "Error: permute() with subgroups(): the subgroup variable `sg' splits some units of `ivar' across multiple subgroups over time (over ALL units, because any cluster can be assigned treatment in a replication). A unit cannot belong to more than one subgroup. To proceed anyway, specify safety(off); share weighting then applies inside the replications as it does in the observed run."
        error 459
      }
      else dis as text "Caution: the subgroup variable `sg' splits some units of `ivar' across multiple subgroups over time; under safety(off) the replications proceed with share weighting, as the observed run does."
    }
    qui drop `sgmin' `sgmax'
    if "`permstrata'" == "" & "`permmethod'" == "shuffle" dis as text "Caution: subgroups() with permute(shuffle) and no permstrata(): the shuffled assignments also re-randomize which subgroups the treated clusters fall in, so the reference distribution of the subgroup-specific statistics mixes subgroup composition with onset timing. When subgroup membership is fixed at the cluster level, permstrata(`sg') holds each subgroup's onset composition fixed across replications."
  }

  * strata: clusters must be strict subsets of strata
  if "`permstrata'" != "" {
    if "`permmethod'" != "shuffle" {
      dis as error "Error: permstrata() is only for permute(shuffle)."
      error 198
    }
    tempvar smin smax
    qui gegen double `smin' = min(`permstrata') if `insmp', by(`clus')
    qui gegen double `smax' = max(`permstrata') if `insmp', by(`clus')
    qui count if `insmp' & (`smin' != `smax' | missing(`permstrata'))
    if r(N) > 0 {
      dis as error "Error: permstrata() requires every cluster to sit inside a single stratum with no missing stratum values (`r(N)' in-sample observation(s) violate this)."
      error 459
    }
  }

  * enough freedom to randomize: at least 2 clusters and 2 distinct onset groups
  tempvar ctag consv
  qui gegen byte `ctag' = tag(`clus') if `insmp'
  qui count if `ctag' == 1
  local permnclus = r(N)
  qui gen double `consv' = cond(`ntr' > 0, `tmin', .) if `ctag' == 1
  qui glevelsof `consv' if `ctag' == 1, missing local(permlvls)
  local permnlvl : word count `permlvls'
  if `permnclus' < 2 | `permnlvl' < 2 {
    dis as error "Error: permute() found `permnclus' cluster(s) across `permnlvl' treatment-onset group(s); there is nothing to randomize."
    error 459
  }

  * permute(probabilities): validate the permprobs() grammar against the data's actual
  * onset levels and cluster values (the same mata parser the engine uses, run in
  * check-only mode, so the two can never drift)
  if "`permmethod'" == "probabilities" {
    qui glevelsof `ttrevar' if `insmp', missing local(permonslvls)
    qui glevelsof `clus' if `insmp', local(permclusvals)
    local permperr
    mata: wooldid_permp(st_local("permprobs"), "`ttrevar'", "`permonslvls'", "`permclusvals'", "", "")
    if "`permperr'" != "" {
      dis as error "`permperr'"
      error 198
    }
  }
end


cap prog drop wooldid_bootctl
prog def wooldid_bootctl
version 16.0
*** F-boot (plan: review_process/bootstrap_plan.md): the bootstrap option SCREEN. It runs
*** no loop, calls no Mata and returns nothing - the engine is wooldid_bootrun and the
*** rewiring is wooldid_bootwire. The deep checks live here rather than in wooldid_validate
*** for compiled-size headroom (wooldid_validate is within ~1 KB of the per-program
*** ceiling). jt() and ptp() are retained-but-UNREAD option slots: they are parsed and
*** passed by wooldid_validate and are kept so the byte-identical-syntax contract between
*** the two programs - checked by test/check_validate_sync.py - and the call site need not
*** change. Do not wire a gate onto them expecting one to already exist, and do not delete
*** them while the call site still passes them.
  syntax , [bootstrap(string) reps(integer 0) ties level(real 0) dots(integer -1) rseed(string) bootalpha(string) bootdf(string) bdisp(string) bts(string) bsave(string) pois(string) permute(string) impv(string) uncond(string) cvce(string asis) fn(string) os(string) ssg(string) attit(string) atti(string) attitime(string) si(string) clus(varlist) ivar(varname) ttrevar(varname) ifr(string asis) jt(string) ptp(string asis) replace]

  * rseed() is validated before the stray screen: it is shared with permute() (since
  * 2026-09-05), whose runs leave this program at the screen's exit below
  if "`rseed'" != "" {
    cap confirm number `rseed'
    local badseed = _rc != 0
    if `badseed' == 0 {
      if mod(`rseed',1) != 0 | `rseed' < 0 | `rseed' > 2147483647 local badseed 1
    }
    if `badseed' {
      dis as error "Error: rseed() must be an integer between 0 and 2,147,483,647 (it is passed to set seed)."
      error 198
    }
  }
  * stray-companion screen: every companion requires the master option; level(),
  * dots() and rseed() are shared with permute()
  if "`bootstrap'" == "" {
    local bootstray
    if `reps' != 0 local bootstray reps()
    if "`ties'" != "" local bootstray `bootstray' ties
    if `dots' != -1 & "`permute'" == "" local bootstray `bootstray' dots()
    if "`rseed'" != "" & "`permute'" == "" local bootstray `bootstray' rseed()
    if "`bootalpha'" != "" local bootstray `bootstray' bootalpha()
    if "`bootdf'" != "" local bootstray `bootstray' bootdf()
    if "`bdisp'" != "" local bootstray `bootstray' bootdisplay()
    if "`bts'" != "" local bootstray `bootstray' boottstat
    if "`bsave'" != "" local bootstray `bootstray' bootsave()
    if "`bootstray'" != "" {
      dis as error "Error: `bootstray' require(s) the bootstrap() option (dots() and rseed() are also accepted with permute(); level() applies to every run)."
      error 198
    }
    exit 0
  }

  * master-option token
  if !inlist("`bootstrap'", "standard", "bayes") {
    dis as error "Error: bootstrap() must be either bootstrap(standard) - cluster resampling with replacement - or bootstrap(bayes) - Bayesian cluster reweighting. Received: bootstrap(`bootstrap')."
    error 198
  }

  * reps() is mandatory (author ruling: no floor - low-rep runs are legitimate; see the
  * help file on how many replications inference needs)
  if `reps' <= 0 {
    dis as error "Error: bootstrap() requires reps(#) with a positive number of replications."
    error 198
  }

  * CS-20 (consistency sweep 2026-09-06, author ruling 2026-09-07: "caution, no floor").
  * permreps() has carried a resolution caution since F036 and a hard floor of 20; reps()
  * carried neither, so reps(20) printed percentile p-values on a .10 grid in silence.
  * The footer that would have said so is conditional, not routine: wooldid_pfloor prints
  * only when a displayed p is exactly 0, and exits under bootdisplay(wald) without
  * boottstat - so a low-rep run whose p-values all happen to be non-zero was told nothing
  * at all. The FLOOR deliberately stays absent: a handful of replications is a legitimate
  * way to time a specification, which is the asymmetry the author ruled defensible (a
  * permutation p has a hard combinatorial floor, a bootstrap p only a Monte-Carlo one).
  * Which number binds depends on the family the table carries. boottstat's studentized p
  * is the add-one rank with floor 2/(reps + 1) (F033, same ruling date). The
  * percentile/bc/bca coefficient p is a plain two-sided share on the 2/reps grid, so it
  * CAN still reject - but only by printing exactly 0, which means no completed
  * replication reached the far side of the null. bootdisplay(wald) without boottstat
  * prints analytic normal/t p-values with no Monte-Carlo grid at all, leaving only the
  * one-sided add-one joint, equality and adversarial family, floor 1/(reps + 1).
  * reps() is the REQUEST: failed replications are excluded from the completed count
  * (e(boot_failed)), so the realized grid can be coarser than the one named here.

  * companion-value validation
  if `level' != 0 & (`level' < 10 | `level' > 99.99) {
    dis as error "Error: level() must be between 10 and 99.99."
    error 198
  }
  * the "not given" sentinels differ by program and are worth stating once: dots() uses -1
  * because dots(0) is a legal user value meaning "suppress" (wooldid_permrun then remaps
  * -1 back to 1), while permctl's permreps() uses 0, which cannot be a legal count.
  * CAUTION: unlike rseed(), this range check was not hoisted above the permute() exit
  * further up, so under permute() a negative dots() passes unvalidated, fails every
  * `if `dots' > 0' test in the engine and silently behaves as dots(0)
  if `dots' != -1 & `dots' < 0 {
    dis as error "Error: dots() must be a nonnegative integer (0 suppresses the progress display)."
    error 198
  }
  * bootsave(): the name rules (plain file name, cwd, replace to overwrite) are checked
  * here BEFORE the replication loop runs. Until 2026-09-05 only wooldid_bootsave checked
  * them, after the loop, so a name clash cost the whole run (permsave() checks in permctl)
  if "`bsave'" != "" {
    local wdbfn `bsave'
    if lower(substr("`wdbfn'", -4, 4)) == ".dta" local wdbfn = substr("`wdbfn'", 1, strlen("`wdbfn'") - 4)
    if strpos("`wdbfn'", ",") {
      * RF-6 (2026-09-13): the check below rejected only a directory part, so
      * bootsave(name, replace) - the obvious analogy to Stata's own saving(f, replace),
      * which this option deliberately does NOT follow - was swallowed whole and wrote a
      * file literally called "name, replace.dta" at rc 0.
      dis as error "Error: bootsave() takes a plain file name and no comma; replace is a separate option of wooldid, so write bootsave(`=subinstr("`wdbfn'", ",", " ", .)') replace."
      error 198
    }
    if strpos("`wdbfn'", "/") | strpos("`wdbfn'", "\") | "`wdbfn'" == "" {
      dis as error "Error: bootsave() takes a plain file name (no directory part); the file is written to the current working directory."
      error 198
    }
    cap confirm file "`wdbfn'.dta"
    if _rc == 0 & "`replace'" == "" {
      dis as error "Error: bootsave(`wdbfn') - `wdbfn'.dta already exists in the working directory; add the replace option to overwrite it."
      error 602
    }
  }
  if "`bootalpha'" != "" {
    if !inlist("`bootalpha'", "size", "equal") {
      dis as error "Error: bootalpha() must be bootalpha(equal) - Dirichlet(1,...,1) across clusters, the Rubin convention and the default - or bootalpha(size) - per-cluster alpha equal to cluster size, the bayesboot id() convention, whose standard errors understate cluster-level dispersion on multi-observation clusters."
      error 198
    }
    if "`bootstrap'" == "standard" {
      dis as error "Error: bootalpha() applies only to bootstrap(bayes)."
      error 198
    }
  }

  * F067 follow-up (bug sweep BS-16): the acceleration jackknife deletes one resampling
  * unit per fit; without cluster() that unit is the observation, so it would be N
  * delete-one-observation fits rather than G. refuse, as bootdf(clusters) does
  if "`bdisp'" == "bca" & "`clus'" == "" {
    dis as error "Error: bootdisplay(bca) requires cluster(): without it the resampling unit is the observation, so the acceleration jackknife would be one delete-one-observation fit per observation (N fits, not G). Give cluster(), or use bootdisplay(bc)."
    error 198
  }

  * bootdf(): the reference distribution for the reported statistics. normal (the
  * default) matches Stata's own bootstrap prefix; clusters uses t with G-1 degrees of
  * freedom, G being the number of resampling units, which is the usual small-sample
  * correction when the number of clusters is modest; a positive integer sets the
  * degrees of freedom directly
  if "`bootdf'" != "" {
    * bootdf(clusters) without cluster() used to resolve silently to t(N-1) - the
    * observation count is the resampling-unit count there (quality review 2026-09-05
    * F040, author ruling 2026-09-06: refuse rather than substitute)
    if "`bootdf'" == "clusters" & "`clus'" == "" {
      dis as error "Error: bootdf(clusters) requires cluster(): without it the resampling unit is the observation, so the cluster-count degrees of freedom would silently be t(N-1), indistinguishable from the normal. Use bootdf(normal) - the default without cluster() - or give the degrees of freedom as a positive integer."
      error 198
    }
    if !inlist("`bootdf'", "normal", "clusters") {
      cap confirm number `bootdf'
      local baddf = _rc != 0
      if `baddf' == 0 {
        if mod(`bootdf',1) != 0 | `bootdf' < 1 local baddf 1
      }
      if `baddf' {
        dis as error "Error: bootdf() must be bootdf(normal) - standard normal critical values, the default and the convention Stata's own bootstrap prefix uses - or bootdf(clusters) - t with (number of resampling clusters - 1) degrees of freedom - or a positive integer setting the degrees of freedom directly."
        error 198
      }
    }
  }

  * bootdisplay(): which inference family the printed tables show. Every family is
  * computed and stored in e() on every bootstrap run regardless (author ruling
  * 2026-09-03); this option chooses only what is displayed. percentile (the default)
  * shows the classic asymmetric percentile interval, the percentile p dual to it, and
  * the standard error implied by inverting that interval - which deliberately does NOT
  * satisfy b +/- c*se, the author having ruled that the interval is the primary object.
  * sd shows the older internally-consistent Wald package on the bootstrap standard
  * deviation, where se, statistic, p and interval do all reproduce one another
  * boottstat needs each replication to produce a standard error that is the same
  * functional of the replicate data as the observed one. Two paths cannot (author
  * ruling 2026-09-04, measured E[se_r^2]/se_obs^2 of 1.44-1.83): poisson runs its
  * replications on pweights, whose robust meat squares the duplication weight; and
  * bootstrap(bayes) reweights by continuous Dirichlet mass, which is neither an
  * integer count for an fweight nor a set of rows for a physical expansion. Refuse
  * boottstat ALONE - the percentile, bias-corrected and Wald families read only the
  * coefficient draws and are unaffected
  if "`bts'" != "" & "`pois'" != "" {
    dis as error "Error: boottstat may not be combined with poisson. The Poisson replications are estimated on sampling weights, so a replication's own standard error is not the same functional of the resampled data as the observed standard error, and the studentized bootstrap would be anti-conservative. Drop boottstat; the percentile, bias-corrected and Wald bootstrap families remain available with poisson."
    error 198
  }
  if "`bts'" != "" & "`bootstrap'" == "bayes" {
    dis as error "Error: boottstat may not be combined with bootstrap(bayes). The Bayesian bootstrap reweights clusters by continuous Dirichlet mass, which cannot be expressed as a duplication count, so a replication's own standard error is not the same functional of the reweighted data as the observed standard error. Use bootstrap(standard) for a studentized bootstrap, or drop boottstat; the percentile, bias-corrected and Wald families remain available with bootstrap(bayes)."
    error 198
  }
  if "`bdisp'" != "" {
    if !inlist("`bdisp'", "percentile", "bc", "bca", "wald") {
      dis as error "Error: bootdisplay() must be bootdisplay(percentile) - the default: the p-value and confidence interval from the percentile ranks of the bootstrap distribution - or bootdisplay(bc) - the bias-corrected percentile method - or bootdisplay(bca) - the bias-corrected and accelerated method, whose acceleration comes from a delete-one-cluster jackknife of G extra fast fits - or bootdisplay(wald) - the p-value and interval from a normal or t approximation, the choice between them governed by bootdf(). The reported standard error is the standard deviation of the bootstrap distribution under all three; every family is computed and stored in e() whichever is displayed. Under boottstat only bootdisplay(percentile) is available: the studentized draws have no Wald or bias-corrected family."
      error 198
    }
    * bootdisplay(wald) is meaningless after studentizing: boottstat reports the ANALYTIC
    * standard error (there is no standard deviation of a studentized distribution to put
    * in that column), so a normal/t approximation built on a bootstrap standard deviation
    * has no standard error to be built from. Refused rather than silently reinterpreted
    if "`bdisp'" == "wald" & "`bts'" != "" {
      dis as error "Error: bootdisplay(wald) may not be combined with boottstat. A studentized bootstrap reports the analytic standard error and takes its p-value and confidence interval from the percentile ranks of the studentized distribution; there is no bootstrap standard deviation of a studentized statistic for a normal or t approximation to be built on. Drop boottstat for the Wald package, or drop bootdisplay(wald) for the studentized one."
      error 198
    }
    * bootdisplay(bc) under boottstat is refused as well (quality review 2026-09-05 F003,
    * author ruling 2026-09-06). The former BC-t twin applied Efron's z0 shift to the
    * quantiles of the centred t-draws and then reflected them through the percentile-t
    * construction, which moves the interval AGAINST the median bias - its upper endpoint
    * coincided with the naive percentile-of-b_r endpoint that BC was invented to repair
    * (measured on t_r = z - 1: percentile-t [.52, 2.48], BC-t [-.15, 1.48]). No BC-t
    * construction exists in Efron-Tibshirani or Davison-Hinkley; percentile-t is already
    * second-order accurate. Refused rather than silently shown as percentile-t
    if inlist("`bdisp'", "bc", "bca") & "`bts'" != "" {
      dis as error "Error: bootdisplay(`bdisp') may not be combined with boottstat. The studentized (percentile-t) bootstrap is already second-order accurate, and applying the bias-corrected shift to the studentized quantiles and then reflecting them moves the interval against the median bias, so there is no bias-corrected (or accelerated) studentized family. Drop boottstat for the bias-corrected percentile interval, or drop bootdisplay(`bdisp') for the studentized one."
      error 198
    }
  }

  * interaction refusals
  if "`permute'" != "" {
    dis as error "Error: bootstrap() may not be combined with permute(); each orchestrates its own replication loop. Run them in separate wooldid calls."
    error 198
  }
  if "`impv'" != "" {
    dis as error "Error: bootstrap() may not be combined with impvalues; the imputation-method point estimates are cross-cluster means whose replication vector would not correspond to the reported results."
    error 198
  }
  if "`fn'" != "" {
    dis as error "Error: bootstrap() may not be combined with fastnose; fastnose requests point estimates only, leaving nothing for the bootstrap to attach inference to (the replications use it internally)."
    error 198
  }
  if "`uncond'" != "" {
    dis as error "Error: bootstrap() may not be combined with unconditionalse; the reported standard errors would be bootstrap-based, contradicting the unconditional-VCE request."
    error 198
  }
  if `"`cvce'"' != "" {
    dis as error "Error: bootstrap() may not be combined with customvce(); the reported standard errors would be bootstrap-based, contradicting the custom-VCE request."
    error 198
  }
  * oldsyntax arrives NORMALIZED to 0/1 (validate rewrites it early) - gate on "1",
  * never on nonempty (the permctl precedent; a nonempty gate fires on every run)
  if "`os'" == "1" {
    dis as error "Error: bootstrap() is not available under oldsyntax."
    error 198
  }
  if "`ssg'" != "" {
    dis as error "Error: bootstrap() may not be combined with suppresssgprimaryeffects; the suppressed subgroup entries are exactly what the replication vectors need in order to bootstrap the subgroup statistics."
    error 198
  }
  if wordcount("`clus'") > 1 {
    dis as error "Error: bootstrap() supports a single (numeric) cluster() variable; multiway clustering is not available with bootstrap-based standard errors."
    error 198
  }

  * cell-weighted estimand gates (author ruling 2026-08-30). att_it/att_i/att_itime
  * aggregate with equal weight per treatment cell, so a within-cluster weight
  * multiplier cancels out of them: reweighting-based replications hold the cell set
  * fixed and omit the cell-sampling variance component.
  local cellwtd `attit' `atti' `attitime'
  if "`bootstrap'" == "bayes" & "`cellwtd'" != "" {
    if "`si'" != "off" {
      dis as error "Error: bootstrap(bayes) with `cellwtd' reweights clusters but cannot resample treatment cells, so those bootstrap standard errors omit the cell-sampling variance component and can be too small when treatment effects vary across clusters. To proceed anyway with a caution, specify safety(off); bootstrap(standard) resamples cells and does not have this limitation."
      error 198
    }
    dis as text "Caution: bootstrap(bayes) standard errors for `cellwtd' omit the cell-sampling variance component (cluster reweighting cannot resample treatment cells) and can be too small when treatment effects vary across clusters. Everything computed from that covariance inherits the understatement - the joint tests, the pairwise subgroup contrasts and the confidence intervals for those estimands, where p-values can be extreme - so read them alongside a bootstrap(standard) run, which does not have this limitation."
  }
  * the physical re-id leg - taken for the cell-weighted estimands, and since 2026-09-04
  * also for boottstat under cluster() - requires the cluster() variable to nest i. The
  * two routes reach it for different reasons and get different messages
  local wdb_btsclus = ("`bts'" != "" & "`clus'" != "" & "`bootstrap'" == "standard")
  if "`bootstrap'" == "standard" & ("`cellwtd'" != "" | `wdb_btsclus') {
    * the resampling copies of a drawn cluster must be re-identifiable as distinct
    * units, which is defined only when every i lies wholly inside one resampling
    * cluster. observation-level resampling (no cluster()) can never satisfy this.
    if "`clus'" == "" & "`cellwtd'" != "" {
      dis as error "Error: bootstrap(standard) with `cellwtd' requires cluster(): observation-level resampling cannot re-identify duplicated treatment cells. Specify cluster() at or above the i level, or request att/ame estimands only."
      error 198
    }
    if "`clus'" != "" {
      tempvar smp cmin cmax
      if `"`ifr'"' != "" qui gen byte `smp' = ((`ifr')==1) & !missing(`clus')
      else qui gen byte `smp' = !missing(`clus')
      qui gegen double `cmin' = min(`clus') if `smp', by(`ivar')
      qui gegen double `cmax' = max(`clus') if `smp', by(`ivar')
      cap assert `cmin' == `cmax' if `smp'
      if _rc != 0 {
        if "`cellwtd'" != "" dis as error "Error: bootstrap(standard) with `cellwtd' requires the cluster() variable to nest i (every i wholly inside one cluster): resampled clusters are duplicated as whole units, which is undefined when a unit straddles clusters. Cluster at or above the i level, or request att/ame estimands only."
        else dis as error "Error: boottstat with cluster() requires the cluster() variable to nest i (every i wholly inside one cluster). A studentized bootstrap has to estimate each replication on a physical expansion of the resampled clusters so that its standard error clusters on the resampled cluster identity, and that expansion is undefined when a unit straddles clusters. Cluster at or above the i level, or drop boottstat (the percentile, bias-corrected and Wald families do not need it)."
        error 459
      }
    }
  }

  * cluster-count sanity (resampling/reweighting needs at least two resampling units)
  if "`clus'" != "" {
    tempvar smp2 cgrp
    if `"`ifr'"' != "" qui gen byte `smp2' = ((`ifr')==1) & !missing(`clus')
    else qui gen byte `smp2' = !missing(`clus')
    qui gegen double `cgrp' = group(`clus') if `smp2'
    qui gstats sum `cgrp' if `smp2', nod
    if r(max) < 2 {
      dis as error "Error: bootstrap() requires at least 2 clusters in the estimation sample."
      error 459
    }
  // RF-2 (2026-09-13): this caution was printed BEFORE the boottstat nesting check and the
  // two-cluster check, both of which abort with error 459 - so a user whose run never
  // executed was advised about its p-value resolution. Moved below every aborting
  // validation, above the other bootstrap cautions.
  if `reps' > 0 {
    local wdbgrid = cond("`bts'" != "", 2/(`reps' + 1), cond("`bdisp'" == "wald", 1/(`reps' + 1), 2/`reps'))
    if `wdbgrid' > .05 & "`bts'" != "" dis as text "Caution: with reps(`reps') the smallest studentized coefficient p-value this run can produce is " %5.3f `wdbgrid' " (2/(reps + 1)), above .05, so no coefficient p-value can reject at the 5% level whatever the data say, and a p-value printed at that value is the resolution limit rather than evidence. reps() has no minimum - small values are a legitimate way to time a specification - but use at least 39 replications for a 5% test (roughly 100 for .02, 200 for .01 and 1000 for .002). Failed replications are excluded, so the realized resolution can be coarser than this."
    else if `wdbgrid' > .05 & "`bdisp'" == "wald" dis as text "Caution: with reps(`reps') the joint, equality and adversarial p-values - one-sided add-one shares over their own completed replication counts - cannot fall below " %5.3f `wdbgrid' " (1/(reps + 1)), above .05, so none of them can reject at the 5% level whatever the data say. The coefficient p-values under bootdisplay(wald) are analytic and carry no such limit. reps() has no minimum - small values are a legitimate way to time a specification - but use at least 19 replications for a 5% joint test (roughly 50 for .02, 100 for .01 and 500 for .002). Failed replications are excluded, so the realized resolution can be coarser than this."
    else if `wdbgrid' > .05 dis as text "Caution: with reps(`reps') the coefficient p-values fall on a grid of " %5.3f `wdbgrid' " (2/reps), so the only value below .05 this run can produce is exactly 0 - which records that no completed replication reached the far side of the null, not that the evidence is that strong. reps() has no minimum - small values are a legitimate way to time a specification - but use at least 40 replications for a 5% test (roughly 100 for .02, 200 for .01 and 1000 for .002). Failed replications are excluded, so the realized grid can be coarser than this. The joint, equality and adversarial p-values are one-sided add-one shares with floor 1/(reps + 1)."
  }
    if r(max) < 10 dis as text "Caution: only " r(max) " clusters are available for the bootstrap; bootstrap inference is unreliable with very few clusters."
    * F050 (quality review 2026-09-05, author ruling 2026-09-06): the total-cluster caution
    * above cannot see the case the few-cluster literature is about - few TREATED clusters.
    * The pairs cluster bootstrap keeps only the replications that draw at least one
    * treated cluster (the others fail and are excluded), so the retained distribution is
    * conditional on that event, and MacKinnon-Webb (2017, 2018) find it badly sized there;
    * the Bayesian bootstrap shares its first two moments. The wild cluster bootstrap-t is
    * the literature's remedy and is NOT offered (a self-contained wild bootstrap of the
    * margins functional would be a new engine with a refit per replication - feasibility
    * verdict F1, 2026-09), so the caution names permute() and boottest on the stage-1
    * regression. Keyed on the treated resampling clusters; the unclustered leg resamples
    * observations and is not covered
    if "`ttrevar'" != "" {
      tempvar trt ctr ctag
      qui gen byte `trt' = (`ttrevar' < .) if `smp2'
      qui gegen byte `ctr' = max(`trt') if `smp2', by(`cgrp')
      qui gegen byte `ctag' = tag(`cgrp') if `smp2' & `ctr' == 1
      qui count if `ctag' == 1
      if r(N) < 6 {
        dis as text "Caution: only " r(N) " treated cluster(s) are available for resampling. What is implemented is the " cond("`bootstrap'" == "bayes", "Bayesian bootstrap, which shares the pairs cluster bootstrap's first two moments", "pairs cluster bootstrap, which keeps only the replications that draw at least one treated cluster, so the retained distribution is conditional on that event") "; with few treated clusters the literature finds this bootstrap badly sized (MacKinnon and Webb 2017, 2018), and boottstat is not a cure for it. Consider permute() - randomization inference over the treatment assignment - or the wild cluster bootstrap-t (boottest, on the stage-1 regression outside wooldid), which wooldid does not offer."
      }
    }
  }

  * RECORD OF TWO LIFTED REFUSALS. Nothing in this block describes a live gate - both
  * refusals it records were removed, and the entries are kept so that neither is
  * reinstated by someone who rediscovers the original reasoning. (The staged refusals this
  * block once also carried - of jointtests and of bootstrap(standard) itself, each "until
  * the engine lands in a later commit" - are gone with their commits: the engine is
  * wooldid_bootrun and the bootstrapped joint tests are wooldid_bootjnt, which the shell
  * calls with jt(`jointtests').)
  * cleanmatrices: a refusal was shipped here and in wooldid_permctl and has been REMOVED.
  * The r(3499) it guarded only reproduces when wooldid is loaded with `do <path>`, as this
  * project's own harness does; an adopath-resolved ado recovers, because Stata re-reads it
  * when a Mata function is missing. The refusal was therefore blocking a combination that
  * works for every installed user - see SI-17/SI-18.
  * pretrendpower: the plan-A13 refusal was LIFTED 2026-09-04 on the author's ruling. The
  * diagnostic now runs from the shell AFTER wooldid_bootwire, so Roth's own formula is
  * fed the bootstrap event-study covariance instead of the analytic one. The refusal had
  * rested on an ordering constraint - the ptp subprogram ran inside the worker, before
  * any rewiring - not on an impossibility. The author declined a simulation-based power
  * calculation as going beyond what Roth (2022) specifies for a feature that cites him
end

cap prog drop wooldid_bootsave
prog def wooldid_bootsave
version 16.0
  * bootsave(): write the replication accumulators to <file>.dta in the working directory.
  * Plain file name only (it is the user's stated design: the file lands in the cwd), an
  * existing file needs wooldid's replace option, and the build runs in its own frame so
  * the user's data are untouched
  syntax , file(string) acc(string) obsb(string) [seacc(string) aseobs(string) pwseacc(string) wmat(string) wobs(real 1) jwacc(string) jobs(string) awacc(string) aobs(string) seed(string) rng(string) replace]
  local fn `file'
  if lower(substr("`fn'", -4, 4)) == ".dta" local fn = substr("`fn'", 1, strlen("`fn'") - 4)
  if strpos("`fn'", "/") | strpos("`fn'", "\") | "`fn'" == "" {
    dis as error "Error: bootsave() takes a plain file name (no directory part); the file is written to the current working directory."
    error 198
  }
  cap confirm file "`fn'.dta"
  if _rc == 0 & "`replace'" == "" {
    dis as error "Error: bootsave(`fn') - `fn'.dta already exists in the working directory; add the replace option to overwrite it."
    error 602
  }
  tempname bsf
  cap frame drop `bsf'
  frame create `bsf'
  cap noi frame `bsf' {
    mata: wooldid_bootsavemk("`acc'", "`obsb'", "`seacc'", "`aseobs'", "`pwseacc'", "`wmat'", `wobs', "`jwacc'", "`jobs'", "`awacc'", "`aobs'")
    * F059: the reproducibility record travels with the file
    char _dta[wooldid_engine] "bootstrap"
    char _dta[wooldid_seed] "`seed'"
    char _dta[wooldid_rngstate] "`rng'"
    label var _rep "bootstrap replication (0 = observed estimates)"
    save "`fn'.dta", replace
    * `_N' here is the NEW frame's row count and the local survives the closing brace,
    * because `frame X { }' switches the current DATASET, not the macro scope - the local
    * belongs to this program either way. Tidying this line out of the frame block would
    * silently count the USER's observations instead and report them as replication rows
    local wdbsn = _N
  }
  * `_rc' must be captured on the line immediately after the cap block and the `cap frame
  * drop' must stay BELOW it: _rc is reset by the very next command, so swapping the two
  * turns a failed save into rc 0 - the error below is skipped and the program prints "the
  * replication draws were written to ..." for a file that was never created. wooldid_permsave
  * has the identical shape and the same constraint
  local wdbsrc = _rc
  cap frame drop `bsf'
  if `wdbsrc' != 0 {
    dis as error "Error: bootsave(`fn') could not build the replication dataset (rc=`wdbsrc')."
    error `wdbsrc'
  }
  dis as text "bootsave(): the replication draws were written to `fn'.dta (`=`wdbsn'-1' replication rows plus the _rep = 0 observed row)."
end

cap prog drop wooldid_bootcellsx
prog def wooldid_bootcellsx, eclass
version 16.0
  * rewires e(gteffects) and every e(cfx*) matrix to bootstrap inference from the cell
  * accumulators filled in wooldid_bootrun. Same families and the same conventions as the
  * coefficient rows; se is the bootstrap SD (non-studentized) or the analytic cell se
  * (studentized). Twins are posted per family so every family is retrievable
  syntax , cell(string) [ncfx(integer 0) names(string) bdisp(string) bts(string) level(real 0) bootdf(string) nclus(integer 0) ties]
  local bdisp = cond("`bdisp'" == "", "percentile", "`bdisp'")
  * F067: the cells have no jackknife (the acceleration would need a delete-one-cluster
  * jackknife of every cell), so under bootdisplay(bca) they carry the bias-corrected
  * family, said so by e(boot_gteffects_family) and a caution
  if "`bdisp'" == "bca" {
    local bdisp bc
    dis as text "Caution: the gteffects/cfx cells carry the bias-corrected (bc) family under bootdisplay(bca): the BCa acceleration is not computed for the cells (it would need a delete-one-cluster jackknife of every cell)."
  }
  * the accumulator index k pairs with the k-th NAME wooldid_bootrun recorded; a count
  * mismatch would misalign the cell tables silently, so it is an error, not a fallback
  if `ncfx' > 0 & `: word count `names'' != `ncfx' {
    dis as error "Error: wooldid_bootcellsx received `ncfx' continuous-treatment cell accumulators but `: word count `names'' matrix names; the bootstrap cell rewiring cannot pair them. This is an internal inconsistency - please report it."
    error 498
  }
  * the shell's level() defaults to 0 (meaning c(level)); resolve it here exactly as the
  * coefficient-row rewiring does, else al = 0.5 collapses every interval onto the median
  if `level' == 0 local level = c(level)
  local wdbdf = .
  if "`bootdf'" != "" & "`bootdf'" != "normal" {
    if "`bootdf'" == "clusters" local wdbdf = `nclus' - 1
    else local wdbdf = `bootdf'
    if `wdbdf' < 1 local wdbdf = .
  }
  local zc = cond(`wdbdf' < ., invt(`wdbdf', 1 - (100 - `level')/200), invnormal(1 - (100 - `level')/200))
  local famlist percentile bc wald
  * studentized: percentile-t only (bootdisplay(bc) is refused with boottstat, F003)
  if "`bts'" != "" local famlist percentile
  local wdbcellmissall 0
  * ---- gteffects: cols estimate(5) se(6) p(7) lb(8) ub(9)
  cap confirm matrix e(gteffects)
  if _rc == 0 {
    cap confirm matrix `cell'_gt
    if _rc == 0 {
      tempname gtsel
      foreach fam of local famlist {
        local ftag = cond("`fam'" == "percentile", "pctl", "`fam'")
        tempname gt`fam'
        matrix `gt`fam'' = e(gteffects)
        * NOT dead code: `wdbcellmiss' is an OUT parameter, written back by
        * wooldid_bootcellwire through st_local, so the 0 here is a reset rather than a
        * value and deleting these lines deletes the user-facing caution at the foot of
        * this program. The `fam' guard on the accumulation below is there because
        * bootcellwire is called once per family over the same accumulator, so counting
        * every family would multiply the reported cell count by the family count;
        * percentile is always in famlist, which is why it is the one counted. And the
        * number is narrow: it counts ONLY studentized cells whose observed analytic se is
        * missing or non-positive. Cells blanked because their key is absent from the rep-0
        * accumulator, because the observed estimate is missing, because fewer than two
        * usable draws survived, or because the bootstrap SD came back non-positive are
        * blanked silently and are not in it
        local wdbcellmiss 0
        if "`bts'" != "" mata: wooldid_bootcellwire("`gt`fam''", 1, "`cell'_gt", "`cell'_gtse", 5, 6, 7, 8, 9, "`fam'", `level', `zc', `wdbdf', ("`ties'" != ""))
        else            mata: wooldid_bootcellwire("`gt`fam''", 1, "`cell'_gt", "", 5, 6, 7, 8, 9, "`fam'", `level', `zc', `wdbdf', ("`ties'" != ""))
        if "`fam'" == "percentile" local wdbcellmissall = `wdbcellmissall' + `wdbcellmiss'
        matrix colnames `gt`fam'' = `: colnames e(gteffects)'
        if "`fam'" == "`bdisp'" matrix `gtsel' = `gt`fam''
        * `ereturn matrix name = mat' MOVES its source: `gt`fam'' does not exist after
        * this line. That is why the colnames and the selected-family copy above must run
        * first. The cfx loop below has the same ordering requirement with a nastier
        * failure mode - see the note at its own move
        ereturn hidden matrix boot_`ftag'_gteffects = `gt`fam''
      }
      * `gtsel' is set iff `bdisp' is a member of `famlist', which is an invariant
      * assembled three places away: the empty->percentile default and the bca->bc remap
      * above, plus wooldid_bootctl's refusals of bootdisplay(wald), bootdisplay(bc) and
      * bootdisplay(bca) under boottstat (F003 rulings 2026-09-04/06), which is what makes
      * the `famlist percentile' cut safe. The `cap confirm' is a belt, NOT a fallback: if
      * the invariant ever breaks - a new bootdisplay() token, a new famlist branch, one of
      * those refusals relaxed - e(gteffects) and the cfx matrices keep their ANALYTIC
      * se/p/ll/ul while the family label posted on the next line announces a bootstrap
      * family, i.e. analytic inference under a bootstrap header with no message. The two
      * lines are meant to stand or fall together
      cap confirm matrix `gtsel'
      if _rc == 0 ereturn matrix gteffects = `gtsel'
      ereturn local boot_gteffects_family `bdisp'
    }
  }
  * ---- cfx matrices: cols b(1) se(2) t(3) pvalue(4) ll(5) ul(6) when the full shape is present
  forval c = 1/`ncfx' {
    local nm : word `c' of `names'
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    cap confirm matrix `cell'_cfx`c'
    if _rc != 0 continue
    * a cfx matrix has exactly two legal shapes: 6 columns (rows 1..6 of r(table)
    * transposed) or 1 column (b alone), the latter under cfxplotnose, where there are no
    * inference columns to fill - nothing in between is ever produced, which is what makes
    * the hard-coded 1, 2, 4, 5, 6 below safe. wooldid_bootrun still built an accumulator
    * for the skipped matrices (it counts every e(wdid_cfxmats) entry), so this leaves an
    * unused accumulator rather than a mismatch - which is why the ncfx/names pairing check
    * above is done on the full list and not on the surviving one
    if colsof(e(`nm')) < 6 continue
    tempname cfsel
    foreach fam of local famlist {
      local ftag = cond("`fam'" == "percentile", "pctl", "`fam'")
      tempname cf`fam'
      matrix `cf`fam'' = e(`nm')
      local wdbcellmiss 0
      if "`bts'" != "" mata: wooldid_bootcellwire("`cf`fam''", 0, "`cell'_cfx`c'", "`cell'_cfx`c'se", 1, 2, 4, 5, 6, "`fam'", `level', `zc', `wdbdf', ("`ties'" != ""))
      else            mata: wooldid_bootcellwire("`cf`fam''", 0, "`cell'_cfx`c'", "", 1, 2, 4, 5, 6, "`fam'", `level', `zc', `wdbdf', ("`ties'" != ""))
      if "`fam'" == "percentile" local wdbcellmissall = `wdbcellmissall' + `wdbcellmiss'
      * the t column is b/se under the reported se, whichever family
      * NAMESPACED (author ruling 2026-09-11, documentation review D-3). A bare `M =' at a
      * `mata:' prompt creates a GLOBAL Mata object that outlives the statement - Mata has no
      * statement-local scope here - so the old `M' silently overwrote any Mata matrix the
      * user held under that name, once per family per cfx matrix, with no error raised. The
      * three sibling one-liners (wooldid_permwire's gteffects blanking and both statements
      * in wooldid_permcellsx) were renamed with it
      mata: wooldid_Mtmp = st_matrix("`cf`fam''"); wooldid_Mtmp[., 3] = wooldid_Mtmp[., 1] :/ wooldid_Mtmp[., 2]; st_matrix("`cf`fam''", wooldid_Mtmp)
      matrix colnames `cf`fam'' = `: colnames e(`nm')'
      matrix rownames `cf`fam'' = `: rowfullnames e(`nm')'
      if "`fam'" == "`bdisp'" matrix `cfsel' = `cf`fam''
      * e() names are capped at 32 characters; a cfx matrix name too long for the twin
      * prefix keeps only its selected family in the live matrix
      * this move (`ereturn matrix' erases its source) is CONDITIONAL on the 32-character
      * e()-name cap, so an ordering mistake above - reading `cf`fam'' after this line -
      * would abort r(111) only when the cfx matrix name happens to be short, and pass when
      * it is long. That is a failure whose presence depends on the length of the user's
      * continuous-treatment variable name and which reads as a data-dependent bug
      if strlen("boot_`ftag'_`nm'") <= 32 ereturn hidden matrix boot_`ftag'_`nm' = `cf`fam''
    }
    cap confirm matrix `cfsel'
    if _rc == 0 ereturn matrix `nm' = `cfsel'
  }
  if `wdbcellmissall' > 0 dis as text "Caution: `wdbcellmissall' gteffects/cfx cell(s) have no analytic standard error to studentize by (for example the ame covariance is withheld under cell containment), so under boottstat their bootstrap inference is left blank rather than filled from the cell draws; the non-studentized bootstrap (without boottstat) reports those cells from the draws."
end

cap prog drop wooldid_bootjnt
prog def wooldid_bootjnt, eclass
version 16.0
*** F-boot commit E: bootstrapped joint tests, subgroup contrast inference and CI
*** flavors, all from the name-aligned replication matrix. Every inverted quadratic
*** form uses a LISTWISE complete-case covariance with a PSD/inversion guard
*** (amendment A11; a dropped constraint posts MISSING plus a caution, never a
*** residual-rank p - the ed780ac precedent). Joint-test e() names and the
*** per-subgroup macro-vs-scalar split mirror the analytic surface exactly, so the
*** existing display prints them unchanged (amendment A14).
  syntax , [level(real 0) acc(string) obsb(string) comp(string) jt(string) ties done(string) bootdf(string) bts(string) bdisp(string) nclus(integer 0) pwseacc(string) jwacc(string) jobs(string) aseobs(string) seacc(string) awacc(string) aobs(string) adv(string) advalt(string) advesrel(string) advtt(string) clus(varlist) jkacc(string)]
  * a pre-rewrite copy of the observed pairwise table: column 2 is the analytic margins se,
  * which the polarisation route needs as its observed comparator. Taken before anything
  * below touches e(pwsgs)
  tempname pw0
  * NOTE the deliberately missing guard: with e(pwsgs) absent this `cap matrix' leaves a
  * FABRICATED 1x1 missing matrix (the F2-class trap the walk below guards against with
  * `cap confirm'), so the failure would surface far away as a conformability error on
  * `pw0'[1..., 2]' and only after joint-test scalars had already been posted - a
  * half-rewired e(). What keeps that safe is a cross-program chain: the read runs only
  * under `wdb_sgok', which requires `pwseacc', which wooldid_bootrun builds only when the
  * observed run posted a pairwise table. Relaxing the `wdb_sgok' test removes it
  cap matrix `pw0' = e(pwsgs)
  local seacc0 `seacc'
  local bdisp = cond("`bdisp'" == "", "percentile", "`bdisp'")
  if `level' == 0 local level = c(level)
  * reference distribution (bootdf()): missing = normal/chi-squared, else t/F with
  * these denominator degrees of freedom
  local wdbdf = .
  if "`bootdf'" != "" & "`bootdf'" != "normal" {
    if "`bootdf'" == "clusters" local wdbdf = `nclus' - 1
    else local wdbdf = `bootdf'
    if `wdbdf' < 1 local wdbdf = .
  }
  tempname wm
  matrix `wm' = e(wdidmainresults)
  local nmain = rowsof(`wm')
  local sgl `e(wdid_sglist)'
  local wdbncaut 0

  * -------- per-column (eq, name, relyr, scope) metadata via the positional walk
  * ALIGNMENT CONTRACT: k is a COLUMN INDEX into `acc' and `obsb', not merely a counter.
  * `obsb' is e(b); `acc' is J(reps, colsof(obsb)) carrying obsb's colnames, and every
  * replication is grabbed BY NAME into those columns (wooldid_bootrun). So "row k of the
  * stacked results tables = column k of e(b)" is load-bearing, and it holds for exactly
  * one reason: e(b), e(wdidmainresults) and every e(wdidsgresults_*) are cut from the
  * single rownamespool the worker builds - main rows first, then each subgroup block in
  * e(wdid_sglist) order, which is the order walked here. The `continue' below (subgroup
  * table absent) would shift every later k while `acc', built from e(b), would not skip:
  * silently mispairing every subsequent estimand. It cannot fire today because the only
  * configuration that drops a subgroup table is suppresssgprimaryeffects, which
  * bootstrap() refuses. If they ever do diverge the fix is a name check of
  * `: colfullnames `acc'' against the walked row names, not a re-derivation
  local kk 0
  local rnames : rowfullnames `wm'
  foreach rn of local rnames {
    local ++kk
    gettoken eqn nmn : rn, parse(":")
    local nmn = subinstr("`nmn'", ":", "", 1)
    local wdbeq`kk' `eqn'
    local wdbnm`kk' `nmn'
    local wdbsg`kk' .
    * column 7 of the results tables is the RELATIVE YEAR. `wdbrl' is that value, and the
    * `< 0' test in the joint blocks below is the only thing that makes the fixed-metric
    * test a PRE-trend test. Change the table layout and this silently selects a different
    * quantity while still producing a well-formed p-value. wooldid_permjnt reads the same
    * magic 7, so the two engines have to move together
    local wdbrl`kk' = `wm'[`kk', 7]
  }
  foreach sg of local sgl {
    tempname wsg
    * confirm before reading: `cap matrix X = e(absent)' yields a fabricated 1x1
    * missing scalar rather than failing (the F2-class trap bootwire already guards)
    cap confirm matrix e(wdidsgresults_`sg')
    if _rc != 0 continue
    matrix `wsg' = e(wdidsgresults_`sg')
    local rnames : rowfullnames `wsg'
    local rr 0
    foreach rn of local rnames {
      local ++kk
      local ++rr
      gettoken eqn nmn : rn, parse(":")
      local nmn = subinstr("`nmn'", ":", "", 1)
      local wdbeq`kk' `eqn'
      local wdbnm`kk' `nmn'
      local wdbsg`kk' `sg'
      local wdbrl`kk' = `wsg'[`rr', 7]
    }
  }
  local K `kk'

  * -------- joint pre-trend tests (only when the user requested jointtests).
  * Author ruling 2026-09-03: a joint test that cannot be BOOTSTRAPPED is ABSENT, not
  * approximated - "the joint test is absent if not bootstrappable. NEVER show an
  * analytic SE when bootstrap is requested." Bootstrapping a joint test means
  * bootstrapping its test STATISTIC, which requires each replication's own covariance
  * among the estimands and therefore requires boottstat (the replications run without
  * fastnose only in that mode). So without boottstat the whole block is skipped and
  * nothing is posted: these names are created ONLY here under bootstrap, since the
  * worker suppresses the analytic joint tests whenever bootstrap() is on, so skipping
  * leaves them genuinely absent rather than stale
  if "`jt'" != "" & "`bts'" == "" {
    * F034 (quality review 2026-09-05, author ruling 2026-09-06): without boottstat the
    * joint pre-trend test is the FIXED-METRIC rank test the permutation engine already
    * uses (wooldid_permqf's construction): the replication draws recentred at the
    * observed estimates (the rep-0 comparator, never the draw mean - the F018 lesson),
    * W_r = (b_r - b)' S^-1 (b_r - b) on the listwise bootstrap covariance S of the
    * family's pre-horizon columns, against W_obs = b' S^-1 b, p = add-one share of
    * W_r >= W_obs. It is a bootstrap quantity (the ruling that blanked these tests asked
    * for one), not the Wald-on-V_boot test a hand `test' reproduces, and is labelled
    * e(boot_jt_family) = fixedmetric; under boottstat the studentized statistic is kept
    * and labelled studentized
    local wdbjw
    foreach fam in att att_it ame ame_it {
      local ename = cond(strpos("`fam'", "ame"), "joint_Choriz", "joint_horiz")
      if "`fam'" == "att_it" | "`fam'" == "ame_it" local ename `ename'_i
      foreach sgtag in "" `sgl' {
        local sgw = cond("`sgtag'" == "", ".", "`sgtag'")
        local cols
        forval k = 1/`K' {
          * the event-study row pool gives a SINGLE equation name, ES_att_it, to att_it,
          * att_i and att_itime alike (the worker builds one shared horizon block for the
          * three) - which is why the family list above has four entries while the
          * allsgssame loop below has eight: LEVEL rows do keep the estimands distinct
          * (main:att_i, main:att_itime). Adding ES_att_i / ES_att_itime to the family
          * list finds no columns and falls into the silent `continue' below; the walk is
          * not broken, the two row-naming regimes simply differ
          if "`wdbeq`k''" == "ES_`fam'" & "`wdbsg`k''" == "`sgw'" & `wdbrl`k'' < 0 local cols `cols' `k'
        }
        local q : word count `cols'
        if `q' < 2 continue
        tempname Cm Wr
        matrix `Cm' = I(`q')
        local jp .
        local jn 0
        local jw .
        mata: wooldid_bootqf("`acc'", "`obsb'", "`cols'", "`Cm'", "jp", "jn", "jw", "`Wr'")
        local en = cond("`sgtag'" == "", "`ename'", "`ename'_s`sgtag'")
        if `jp' < . {
          * CS-60 (author ruling 2026-09-07): a scalar at ALL posting sites, so a loop over subgroups gets one type
          ereturn scalar `en' = `jp'
          matrix colnames `Wr' = `en'
          tempname Wo
          matrix `Wo' = J(1, 1, `jw')
          matrix colnames `Wo' = `en'
          if "`wdbjw'" == "" {
            tempname jwall jwobs
            matrix `jwall' = `Wr'
            matrix `jwobs' = `Wo'
            local wdbjw 1
          }
          else {
            matrix `jwall' = `jwall', `Wr'
            matrix `jwobs' = `jwobs', `Wo'
          }
        }
        else {
          local ++wdbncaut
          local wdbc`wdbncaut' the bootstrap joint pre-trend test `en' is unavailable (fewer than 2 complete replications carry every pre-treatment horizon, or the bootstrap covariance of those horizons is rank-deficient); it is reported as absent rather than as an analytic test
        }
      }
    }
    if "`wdbjw'" != "" {
      ereturn matrix boot_jointW = `jwall'
      ereturn matrix boot_jointWobs = `jwobs'
      ereturn local boot_jt_family fixedmetric
    }
    else {
      * a flag, not just a caution, so a REPLAYED run still carries the disclosure
      ereturn local boot_jt_absent 1
    }
  }
  if "`jt'" != "" & "`bts'" != "" {
    * BOOTSTRAPPED joint statistics: p = add-one share of recentred W_r >= W_obs, the W_r
    * having been computed inside each replication from its own e(V_es_*) (wooldid_bootjwr)
    * and W_obs here from the observed analytic covariance (wooldid_bootjwobs). This
    * replaces the earlier draw-covariance Wald (wooldid_bootjw), which the author ruled
    * is not a bootstrapped test. A family whose replications gave no usable W_r posts
    * nothing, with a caution - never an analytic p under a bootstrap header
    local wdb_jtok = ("`jwacc'" != "" & "`jobs'" != "")
    if `wdb_jtok' {
      cap confirm matrix `jwacc'
      if _rc != 0 local wdb_jtok 0
    }
    if `wdb_jtok' local wdb_jfk : colnames `jwacc'
    foreach fam in att att_it ame ame_it {
      local ename = cond(strpos("`fam'", "ame"), "joint_Choriz", "joint_horiz")
      if "`fam'" == "att_it" | "`fam'" == "ame_it" local ename `ename'_i
      * overall
      local fk `fam'
      local jp .
      if `wdb_jtok' {
        local jk : list posof "`fk'" in wdb_jfk
        if `jk' > 0 {
          local wobs .
          mata: wooldid_bootjwobs("`fk'", "`jobs'", "wobs")
          if `wobs' < . {
            mata: st_local("jp", strofreal((1 + sum(select(st_matrix("`jwacc'")[., `jk'], st_matrix("`jwacc'")[., `jk'] :< .) :>= strtoreal(st_local("wobs")))) / (sum(st_matrix("`jwacc'")[., `jk'] :< .) + 1), "%21x"))
            mata: st_local("jn", strofreal(sum(st_matrix("`jwacc'")[., `jk'] :< .)))
            if `jn' < 2 local jp .
          }
        }
      }
      if `jp' < . ereturn scalar `ename' = `jp'
      else if `: list posof "`fk'" in wdb_jfk' > 0 {
        local ++wdbncaut
        local wdbc`wdbncaut' the bootstrapped `ename' joint test is unavailable (no usable replication statistic - rank-deficient replication covariances, or a horizon missing from the replications)
      }
      * per subgroup. CS-60 (consistency sweep 2026-09-06, author ruling 2026-09-07): this
      * block used to post an e() MACRO while the no-subgroups case posted a scalar, so a
      * loop over subgroups got a number for one family and a string for the other. All
      * three posting sites - the worker's analytic one, this one and wooldid_permjnt's -
      * are scalars now. Do not put the macro back, and do not write display code that
      * dereferences these as macros (the sweep found such a block too)
      foreach sg of local sgl {
        local fk `fam'_s`sg'
        local jp .
        if `wdb_jtok' {
          local jk : list posof "`fk'" in wdb_jfk
          if `jk' > 0 {
            local wobs .
            mata: wooldid_bootjwobs("`fk'", "`jobs'", "wobs")
            if `wobs' < . {
              mata: st_local("jp", strofreal((1 + sum(select(st_matrix("`jwacc'")[., `jk'], st_matrix("`jwacc'")[., `jk'] :< .) :>= strtoreal(st_local("wobs")))) / (sum(st_matrix("`jwacc'")[., `jk'] :< .) + 1), "%21x"))
              mata: st_local("jn", strofreal(sum(st_matrix("`jwacc'")[., `jk'] :< .)))
              if `jn' < 2 local jp .
            }
          }
        }
        if `jp' < . ereturn scalar `ename'_s`sg' = `jp'   // CS-60: scalar, not macro
        else if `: list posof "`fk'" in wdb_jfk' > 0 {
          local ++wdbncaut
          local wdbc`wdbncaut' the bootstrapped `ename' joint test for subgroup `sg' is unavailable (no usable replication statistic)
        }
      }
    }
    if `wdb_jtok' {
      matrix colnames `jwacc' = `wdb_jfk'
      cap ereturn matrix boot_jointW = `jwacc'
      ereturn local boot_jt_family studentized
    }
  }

  * -------- adversarial line-test p-values (author ruling 2026-09-04: bootstrapped under
  * boottstat, ABSENT otherwise). The worker posted the ANALYTIC p at estimation time; that
  * is never left standing under a bootstrap header. Names mirror the worker's posting:
  * adversarial_pval_<typetext>[_s<sg>], typetext in att att_it ame ame_it
  if "`adv'`advalt'" != "" {
    local wdb_advok = ("`bts'" != "" & "`awacc'" != "" & "`aobs'" != "")
    if `wdb_advok' {
      cap confirm matrix `awacc'
      if _rc != 0 local wdb_advok 0
    }
    local wdb_alt = ("`advalt'" != "")
    if `wdb_advok' local wdb_acols : colnames `awacc'
    local wdb_advblank 0
    local wdb_advmiss 0
    foreach ty in att att_it ame ame_it {
      foreach sgtag in "" `sgl' {
        local fk `ty'
        local ename adversarial_pval_`ty'
        if "`sgtag'" != "" {
          local fk `ty'_s`sgtag'
          local ename `ename'_s`sgtag'
        }
        if "`e(`ename')'" == "" continue
        local jp .
        if `wdb_advok' {
          local acol `fk'_a`wdb_alt'
          local ak : list posof "`acol'" in wdb_acols
          if `ak' > 0 {
            local wobs .
            mata: wooldid_bootadvobs("`fk'", `wdb_alt', "`aobs'", `advesrel', ("`advtt'" != ""), "wobs")
            if `wobs' < . {
              mata: st_local("jp", strofreal((1 + sum(select(st_matrix("`awacc'")[., `ak'], st_matrix("`awacc'")[., `ak'] :< .) :>= strtoreal(st_local("wobs")))) / (sum(st_matrix("`awacc'")[., `ak'] :< .) + 1), "%21x"))
              mata: st_local("jn", strofreal(sum(st_matrix("`awacc'")[., `ak'] :< .)))
              if `jn' < 2 local jp .
            }
          }
        }
        if `jp' < . ereturn scalar `ename' = `jp'
        else {
          ereturn scalar `ename' = .
          if `wdb_advok' local ++wdb_advmiss
          else local ++wdb_advblank
        }
      }
    }
    if `wdb_advblank' > 0 {
      ereturn local boot_adv_absent 1
      local ++wdbncaut
      local wdbc`wdbncaut' the adversarial line-test p-value(s) were not computed: bootstrapping them means bootstrapping the test statistic, which needs boottstat. They are absent rather than shown as analytic tests, which would not be bootstrap quantities
    }
    if `wdb_advmiss' > 0 {
      ereturn local boot_adv_absent 1
      local ++wdbncaut
      local wdbc`wdbncaut' `wdb_advmiss' bootstrapped adversarial line-test p-value(s) are unavailable (no usable replication statistic - degenerate replication covariances, or a horizon missing from the replications)
    }
    if `wdb_advok' {
      cap ereturn matrix boot_advW = `awacc'
    }
  }

  * -------- allsgssame equality chains (mirror the analytic surface: computed
  * whenever the analytic run posted the scalar, i.e. >1 realized sg; blanked-then-
  * filled is now compute-directly)
  if `: word count `sgl'' > 1 {
    * The item list MUST cover every name the worker can post, which is the att family
    * PLUS its contreat twin under `itemtext = subinstr("`item'","att","ame",.)' (worker
    * 4218-4219). That mapping turns att_i into ame_i and att_itime into ame_itime, and
    * an earlier six-name list here omitted exactly those two - so under
    * bootstrap(bayes) + contreat() + subgroups() their ANALYTIC equality p-values were
    * never blanked and printed alongside the footer asserting they were absent, with
    * e(boot_sgjt_absent) set. Keep all eight in step with the worker's posting sites
    foreach it in att att_it att_i att_itime ame ame_it ame_i ame_itime {
      foreach pre in "" "pre_" {
        local eqwant = cond("`pre'" == "", "main", "pretreatment")
        local cols
        local ok 1
        foreach sg of local sgl {
          local found
          forval k = 1/`K' {
            if "`wdbeq`k''" == "`eqwant'" & "`wdbnm`k''" == "`it'" & "`wdbsg`k''" == "`sg'" local found `k'
          }
          if "`found'" == "" local ok 0
          else local cols `cols' `found'
        }
        if "`bts'" == "" {
          * F034 (2026-09-06): the fixed-metric rank test on the chain contrast of the
          * subgroup level columns (recentred draws, bootstrap-covariance metric), the
          * permutation engine's construction. The worker's ANALYTIC scalar is always
          * overwritten or blanked first, so no analytic p can stand under a bootstrap
          * header (the blank is not gated on `ok', which reports whether every
          * subgroup's column could be mapped)
          if "`e(joint_`pre'`it'_allsgssame)'" != "" ereturn scalar joint_`pre'`it'_allsgssame = .
          local jp .
          if `ok' {
            local mm : word count `cols'
            tempname Cm Wr
            matrix `Cm' = J(`mm'-1, `mm', 0)
            forval j = 1/`=`mm'-1' {
              matrix `Cm'[`j', `j'] = 1
              matrix `Cm'[`j', `j'+1] = -1
            }
            local jn 0
            local jw .
            mata: wooldid_bootqf("`acc'", "`obsb'", "`cols'", "`Cm'", "jp", "jn", "jw", "`Wr'")
          }
          if `jp' < . {
            ereturn scalar joint_`pre'`it'_allsgssame = `jp'
            matrix colnames `Wr' = joint_`pre'`it'_allsgssame
            tempname Wo
            matrix `Wo' = J(1, 1, `jw')
            matrix colnames `Wo' = joint_`pre'`it'_allsgssame
            if "`wdbsw'" == "" {
              tempname swall swobs
              matrix `swall' = `Wr'
              matrix `swobs' = `Wo'
              local wdbsw 1
            }
            else {
              matrix `swall' = `swall', `Wr'
              matrix `swobs' = `swobs', `Wo'
            }
          }
          else if `ok' | "`e(joint_`pre'`it'_allsgssame)'" != "" {
            ereturn scalar joint_`pre'`it'_allsgssame = .
            local wdbsgjtblank 1
            * per-item caution carrying the real reason (bug sweep BS-19), mirroring the
            * joint-test twin above; the lumped wording blaming a missing boottstat predated F034.
            * The caution says why at run time; `wdbsgwhy' carries the same reason into e() so
            * the WHY survives a replay, which the bare absent flag never did (bug sweep 2.4)
            if `ok' local wdbsgwhy = trim("`wdbsgwhy' " + cond(strpos(" `wdbsgwhy' ", " nostatistic ") == 0, "nostatistic", ""))
            else local wdbsgwhy = trim("`wdbsgwhy' " + cond(strpos(" `wdbsgwhy' ", " unmapped ") == 0, "unmapped", ""))
            local ++wdbncaut
            if `ok' local wdbc`wdbncaut' the bootstrap `pre'`it' subgroup-equality test is unavailable (fewer than 2 complete replications carry every subgroup level, or the bootstrap covariance of those levels is rank-deficient); it is reported as absent rather than as an analytic test
            else local wdbc`wdbncaut' the `pre'`it' subgroup-equality test could not be mapped to replication columns and is reported as absent
          }
        }
        else if `ok' {
          * NOT gated on the worker having posted the analytic scalar: under boottstat this
          * is a new bootstrap quantity, and the worker withholds its analytic twin in exactly
          * the configurations (clustering at or above i) where the bootstrap one is wanted
          * BOOTSTRAPPED subgroup-equality statistic by polarisation (wooldid_bootsgw): the
          * per-replication level covariance is rebuilt from the replication's own level
          * standard errors (seacc) and pairwise standard errors (pwseacc), W_r is the chain
          * Wald of (b_r - b_obs), and p is the add-one tail against W_obs from the observed
          * analytic level and pairwise standard errors. Replaces the draw-covariance Wald
          local wdbjp .
          local wdbjn 0
          local wdb_sgok = ("`pwseacc'" != "" & "`aseobs'" != "")
          if `wdb_sgok' {
            cap confirm matrix `pwseacc'
            if _rc != 0 local wdb_sgok 0
            cap confirm matrix `aseobs'
            if _rc != 0 local wdb_sgok 0
          }
          if `wdb_sgok' {
            tempname pairidx pwobs0
            * the OBSERVED pairwise se by row, read from the pre-rewrite pw copy taken at the
            * top of this program (column 2 is still the analytic margins se there)
            * ORDERING CONSTRAINT: wooldid_bootpairidx does NOT see `pw0'. It reads the
            * row stripe of the LIVE e(pwsgs) and returns e(pwsgs) ROW NUMBERS, which are
            * then used as column positions in this transposed SNAPSHOT. The two align
            * only because nothing between the capture at the top of this program and here
            * reposts e(pwsgs) - the pairwise rewiring block that does repost it starts
            * below. Move that block above these tests, or repost e(pwsgs) with a changed
            * row order, and every off-diagonal entry of W_obs is built from the wrong
            * pair's standard error: a wrong subgroup-equality p-value, no error, no
            * caution. If this ever has to move, key pwobs0 by name rather than position
            matrix `pwobs0' = `pw0'[1..., 2]'
            * TWO EQUATION VOCABULARIES meet in this block three lines apart. `eqwant'
            * (set above) is "main"/"pretreatment", the equation stripe of the level
            * tables. `eqwant0' is "att"/"pre_att", the equation stripe of e(pwsgs), which
            * wooldid_bootpairidx matches with a plain string compare. Collapsing one into
            * the other raises no error: the index comes back all-missing, wooldid_bootsgw
            * returns on its first off-diagonal pair, and the test is blanked with the
            * caution below blaming an unavailable pairwise margins SE - a wholly
            * misleading diagnosis for what is a name mismatch
            local eqwant0 = cond("`pre'" == "", "`it'", "pre_`it'")
            mata: wooldid_bootpairidx("`eqwant0'", "`sgl'", "`pairidx'")
            mata: wooldid_bootsgw("`acc'", "`seacc0'", "`pwseacc'", "`obsb'", "`aseobs'", "`pwobs0'", "`cols'", "`pairidx'", "wdbjp", "wdbjn")
          }
          if `wdbjp' < . ereturn scalar joint_`pre'`it'_allsgssame = `wdbjp'
          else {
            ereturn scalar joint_`pre'`it'_allsgssame = .
            local wdbsgjtblank 1
            if strpos(" `wdbsgwhy' ", " nostudentizedinputs ") == 0 local wdbsgwhy `wdbsgwhy' nostudentizedinputs
            local ++wdbncaut
            local wdbc`wdbncaut' the bootstrapped `pre'`it' subgroup-equality test is unavailable: forming its studentized statistic needs the observed and per-replication pairwise-contrast standard errors, and at least one of these is missing (the pairwise margins standard error is unavailable in this specification, e.g. under clustering at or above the i level)
          }
        }
        else if "`e(joint_`pre'`it'_allsgssame)'" != "" {
          * boottstat requested but the columns could not be mapped: still never leave
          * the analytic value standing
          ereturn scalar joint_`pre'`it'_allsgssame = .
          local wdbsgjtblank 1
          if strpos(" `wdbsgwhy' ", " unmapped ") == 0 local wdbsgwhy `wdbsgwhy' unmapped
          local ++wdbncaut
          local wdbc`wdbncaut' the `pre'`it' subgroup-equality test could not be mapped to replication columns and is reported as absent
        }
      }
    }
    if "`wdbsw'" != "" {
      ereturn matrix boot_sgW = `swall'
      ereturn matrix boot_sgWobs = `swobs'
    }
    if "`wdbsgjtblank'" != "" {
      ereturn local boot_sgjt_absent 1
      * bug sweep 2.4: the bare flag says THAT an across-subgroup equality test was withheld and
      * never WHY, and the per-item cautions that do say why are gone by the time anyone replays
      * the estimate. One or more of: nostatistic (the test was attempted on the fixed metric and
      * no statistic survived - fewer than 2 complete replications carry every subgroup level, or
      * the bootstrap covariance of those levels is rank-deficient); nostudentizedinputs (under
      * boottstat, the observed or per-replication pairwise standard errors the studentization
      * needs were unavailable); unmapped (the tested columns could not be matched to replication
      * columns). Absence of the flag means no test was withheld, not that none was attempted.
      if "`wdbsgwhy'" != "" ereturn local boot_sgjt_reason `wdbsgwhy'
    }
  }

  * -------- pairwise subgroup contrasts: keep the analytic point differences,
  * rewire se/z/p/ci from the bootstrap covariance of the matching estimand columns
  cap confirm matrix e(pwsgs)
  if _rc == 0 {
    tempname pw pwsd pwpct pwbc pwbca pwtp pwts pwt1
    local wdbpbcamiss 0
    matrix `pw' = e(pwsgs)
    * under boottstat the contrasts are STUDENTIZED (author ruling 2026-09-03: boottstat
    * bootstraps t-statistics for the contrasts too). Column 2 of `pw' is still the
    * observed ANALYTIC pairwise se at this point - bootwire does not touch e(pwsgs) - and
    * it is the studentization's denominator; it is read here, BEFORE any rewrite
    local wdbpwts = ("`bts'" != "" & "`pwseacc'" != "")
    if `wdbpwts' {
      cap confirm matrix `pwseacc'
      if _rc != 0 local wdbpwts 0
    }
    if `wdbpwts' {
      matrix `pwtp' = `pw'
      * DELIBERATELY TRANSPOSED relative to every other table here (statistics in rows,
      * contrasts in columns), because it is posted as e(boot_tstud_pwsgs) with one named
      * column per contrast, matching e(boot_tstud) for coefficients. wooldid_bootpwt
      * returns a 6x1 column laid out (n, t_obs, p_pctl, qlo_pctl, qhi_pctl, live), so the
      * `matrix `pwts'[1, `rr'] = `pwt1'' below drops a whole COLUMN in at (1, rr) rather
      * than setting one cell. Row 6 is a three-state flag - missing = never called
      * (`pwj' == 0), 0 = called and declined, 1 = usable - which is why the guard below
      * tests `== 1' and why the never-called case is safe. The rownames that decode all
      * of this are only set at the end of the block
      matrix `pwts' = J(6, rowsof(`pw'), .)
      local pwtmiss 0
    }
    * both families for the contrast rows, exactly as for the coefficient tables:
    * always computed, always stored, bootdisplay() selects which one is displayed
    matrix `pwsd' = `pw'
    matrix `pwpct' = `pw'
    matrix `pwbc' = `pw'
    matrix `pwbca' = `pw'
    local prn : rowfullnames `pw'
    local rr 0
    local pwunmapped 0
    foreach rn of local prn {
      local ++rr
      gettoken eqn nmn : rn, parse(":")
      local nmn = subinstr("`nmn'", ":", "", 1)
      * eq = item (att/att_it/.../ame; pre-rows carry a pre-marked eq), name = JvsK
      local ispre = strpos("`eqn'", "pre") > 0
      local item = subinstr("`eqn'", "pre", "", 1)
      local item = subinstr("`item'", "_", "", 1)
      if !`ispre' local item `eqn'
      local eqwant = cond(`ispre', "pretreatment", "main")
      * with three or more subgroups Stata marks the base level in the pairwise row
      * name (1vs0bn.sg, not 1vs0.sg), so the marker must be tolerated on BOTH sides
      if regexm("`nmn'", "^([0-9]+)(bn?)?vs([0-9]+)(bn?)?([.].*)?$") {
        local sga = regexs(1)
        local sgb = regexs(3)
        local ka
        local kb
        forval k = 1/`K' {
          if "`wdbeq`k''" == "`eqwant'" & "`wdbnm`k''" == "`item'" & "`wdbsg`k''" == "`sga'" local ka `k'
          if "`wdbeq`k''" == "`eqwant'" & "`wdbnm`k''" == "`item'" & "`wdbsg`k''" == "`sgb'" local kb `k'
        }
        if "`ka'" != "" & "`kb'" != "" {
          local zc = cond(`wdbdf' < ., invt(`wdbdf', 1 - (100 - `level')/200), invnormal(1 - (100 - `level')/200))
          * the sd family for this contrast (the historical behaviour)
          local wdbse .
          * the only Mata call here with no output-name argument: wooldid_bootpw writes
          * the caller's local `wdbse' through a name HARDCODED inside the Mata function
          * (st_local("wdbse", ...)), so the local's name is part of that contract.
          * Renaming it here - the ordinary, safe-looking refactor - lands the Mata write
          * in a local nobody reads: `wdbse' stays missing, every row fails the test
          * below, and columns 2-6 of all four contrast families are blanked with no error
          * and no caution (pwunmapped and pwtmiss never increment)
          mata: wooldid_bootpw("`acc'", `ka', `kb')
          if `wdbse' < . & `wdbse' > 0 {
            matrix `pwsd'[`rr', 2] = `wdbse'
            matrix `pwsd'[`rr', 3] = `pwsd'[`rr', 1]/`wdbse'
            matrix `pwsd'[`rr', 4] = cond(`wdbdf' < ., 2*ttail(`wdbdf', abs(`pwsd'[`rr', 3])), 2*normal(-abs(`pwsd'[`rr', 3])))
            matrix `pwsd'[`rr', 5] = `pwsd'[`rr', 1] - `zc'*`wdbse'
            matrix `pwsd'[`rr', 6] = `pwsd'[`rr', 1] + `zc'*`wdbse'
          }
          else {
            forval cc = 2/6 {
              matrix `pwsd'[`rr', `cc'] = .
            }
          }
          * the percentile family for the same contrast: author ruling 2026-09-03, a
          * 1-vs-1 contrast without boottstat is bootstrapped as a COEFFICIENT - the
          * difference of the two subgroup coefficients - so it gets exactly the package
          * a coefficient row gets, interval first and the standard error inverted from it
          * the contrast's standard error is the standard deviation of the contrast draws
          * (wooldid_bootpw, computed just above) under every non-studentized display -
          * never the interval inverted. Author ruling 2026-09-04
          mata: wooldid_bootpp("`acc'", `ka', `kb', `level', "`comp'", ("`ties'" != ""), "wdplb", "wdpub", "wdpp", "wdpn", "wdpblb", "wdpbub", "wdpbp", "wdpz0")
          if `wdbse' < . & `wdbse' > 0 & `wdplb' < . & `wdpub' < . {
            matrix `pwpct'[`rr', 2] = `wdbse'
            matrix `pwpct'[`rr', 3] = `pwpct'[`rr', 1]/`wdbse'
            matrix `pwpct'[`rr', 4] = `wdpp'
            matrix `pwpct'[`rr', 5] = `wdplb'
            matrix `pwpct'[`rr', 6] = `wdpub'
          }
          else {
            forval cc = 2/6 {
              matrix `pwpct'[`rr', `cc'] = .
            }
          }
          if `wdbse' < . & `wdbse' > 0 & `wdpblb' < . & `wdpbub' < . {
            matrix `pwbc'[`rr', 2] = `wdbse'
            matrix `pwbc'[`rr', 3] = `pwbc'[`rr', 1]/`wdbse'
            matrix `pwbc'[`rr', 4] = `wdpbp'
            matrix `pwbc'[`rr', 5] = `wdpblb'
            matrix `pwbc'[`rr', 6] = `wdpbub'
          }
          else {
            forval cc = 2/6 {
              matrix `pwbc'[`rr', `cc'] = .
            }
          }
          * F067: the BCa twin for the contrast, acceleration from the jackknife contrasts
          local wdpalb .
          local wdpaub .
          local wdpap .
          if "`jkacc'" != "" mata: wooldid_bootppa("`acc'", `ka', `kb', `level', "`comp'", ("`ties'" != ""), "`jkacc'", "wdpalb", "wdpaub", "wdpap")
          if `wdbse' < . & `wdbse' > 0 {
            matrix `pwbca'[`rr', 2] = `wdbse'
            matrix `pwbca'[`rr', 3] = `pwbca'[`rr', 1]/`wdbse'
            * BS-17: an unavailable BCa family (no acceleration, or an adjusted endpoint level
            * outside (0,1)) blanks the p-value and interval only; the bootstrap SD stands
            if `wdpalb' < . & `wdpaub' < . {
              matrix `pwbca'[`rr', 4] = `wdpap'
              matrix `pwbca'[`rr', 5] = `wdpalb'
              matrix `pwbca'[`rr', 6] = `wdpaub'
            }
            else {
              forval cc = 4/6 {
                matrix `pwbca'[`rr', `cc'] = .
              }
              local ++wdbpbcamiss
            }
          }
          else {
            forval cc = 2/6 {
              matrix `pwbca'[`rr', `cc'] = .
            }
          }
          * the studentized (percentile-t) family for this contrast (boottstat): se_analytic,
          * t_obs, rank p, equal-tailed interval. No BC-t twin: bootdisplay(bc) is refused
          * under boottstat (F003, 2026-09-06)
          if `wdbpwts' {
            local pwj = 0
            local pwacn : colfullnames `pwseacc'
            local pwj : list posof "`rn'" in pwacn
            if `pwj' > 0 & `pw'[`rr', 2] < . & `pw'[`rr', 2] > 0 {
              mata: wooldid_bootpwt("`acc'", `ka', `kb', "`pwseacc'", `pwj', st_matrix("`pw'")[`rr', 1], st_matrix("`pw'")[`rr', 2], `level', ("`ties'" != ""), "`pwt1'")
              matrix `pwts'[1, `rr'] = `pwt1'
            }
            if `pwts'[6, `rr'] == 1 {
              matrix `pwtp'[`rr', 2] = `pw'[`rr', 2]
              matrix `pwtp'[`rr', 3] = `pwts'[2, `rr']
              matrix `pwtp'[`rr', 4] = `pwts'[3, `rr']
              * the endpoints are CROSSED on purpose and both subtract: the equal-tailed
              * percentile-t interval is [b - q_hi*se, b - q_lo*se] on the signed
              * quantiles of t(r) = (b_r - b)/se_r, and that inversion is what makes it
              * second-order accurate and deliberately asymmetric about b (the same
              * construction wooldid_boottp uses for coefficients; boottstat in the help).
              * "Fixing" the apparent row swap ships every studentized contrast interval
              * mirrored about the point estimate - the right width, the wrong interval,
              * and nothing in the suite would flag it
              matrix `pwtp'[`rr', 5] = `pw'[`rr', 1] - `pwts'[5, `rr']*`pw'[`rr', 2]
              matrix `pwtp'[`rr', 6] = `pw'[`rr', 1] - `pwts'[4, `rr']*`pw'[`rr', 2]
            }
            else {
              forval cc = 2/6 {
                matrix `pwtp'[`rr', `cc'] = .
              }
              local ++pwtmiss
            }
          }
        }
        else {
          * unmapped: blank rather than leave the analytic se/t/p/CI in a table the
          * display attributes to the bootstrap covariance (the commit-C blanking
          * invariant, restored here for every future label surprise as well)
          forval cc = 2/6 {
            matrix `pwsd'[`rr', `cc'] = .
            matrix `pwpct'[`rr', `cc'] = .
            matrix `pwbc'[`rr', `cc'] = .
            matrix `pwbca'[`rr', `cc'] = .
            if `wdbpwts' matrix `pwtp'[`rr', `cc'] = .
          }
          local ++pwunmapped
        }
      }
      else {
        forval cc = 2/6 {
          matrix `pwsd'[`rr', `cc'] = .
          matrix `pwpct'[`rr', `cc'] = .
          matrix `pwbc'[`rr', `cc'] = .
          matrix `pwbca'[`rr', `cc'] = .
          if `wdbpwts' matrix `pwtp'[`rr', `cc'] = .
        }
        local ++pwunmapped
      }
    }
    local pwc : colnames `pw'
    local pwc : subinstr local pwc "t" "`=cond(`wdbdf' < ., "t", "z")'", word
    local pwc : subinstr local pwc "z" "`=cond(`wdbdf' < ., "t", "z")'", word
    * the CI columns carry the requested level here too (F5): bootwire renames the
    * main and subgroup tables, and this table must not keep saying 95
    if `level' != 95 {
      local wdbplvl = cond(mod(`level',1)==0, strofreal(`level'), subinstr(strofreal(`level',"%9.0g"), ".", "p", .))
      local pwc : subinstr local pwc "lb_95ci" "lb_`wdbplvl'ci", word
      local pwc : subinstr local pwc "ub_95ci" "ub_`wdbplvl'ci", word
    }
    * select the displayed family, then post both twins. `ereturn matrix' MOVES its
    * source, so the selection has to happen first
    if `wdbpwts' {
      * studentized display: percentile-t (bootdisplay(bc) and bootdisplay(wald) are
      * refused with boottstat)
      matrix `pw' = `pwtp'
      local pwc : subinstr local pwc "t" "t_pt", word
      local pwc : subinstr local pwc "z" "t_pt", word
      local pwc : subinstr local pwc "se" "se_analytic", word
    }
    else if "`bdisp'" == "percentile"  matrix `pw' = `pwpct'
    else if "`bdisp'" == "bc"           matrix `pw' = `pwbc'
    else if "`bdisp'" == "bca"          matrix `pw' = `pwbca'
    else                                matrix `pw' = `pwsd'
    local pwc0 : subinstr local pwc "t_pt" "`=cond(`wdbdf' < ., "t", "z")'", word
    local pwc0 : subinstr local pwc0 "se_analytic" "se", word
    * F053: the rank-based displays label the Wald statistic column t_wald / z_wald
    local pwcw : subinstr local pwc0 "t" "t_wald", word
    local pwcw : subinstr local pwcw "z" "z_wald", word
    if !`wdbpwts' & "`bdisp'" != "wald" local pwc `pwcw'
    matrix colnames `pw' = `pwc'
    matrix colnames `pwsd' = `pwc0'
    matrix colnames `pwpct' = `pwcw'
    matrix colnames `pwbc' = `pwcw'
    matrix colnames `pwbca' = `pwcw'
    ereturn matrix pwsgs = `pw'
    ereturn hidden matrix boot_wald_pwsgs = `pwsd'
    ereturn hidden matrix boot_pctl_pwsgs = `pwpct'
    ereturn hidden matrix boot_bc_pwsgs = `pwbc'
    if "`jkacc'" != "" ereturn hidden matrix boot_bca_pwsgs = `pwbca'
    if `wdbpwts' {
      matrix colnames `pwtp' = `pwc'
      ereturn hidden matrix boot_tpctl_pwsgs = `pwtp'
      matrix rownames `pwts' = n t_obs p_pctl qlo_pctl qhi_pctl live
      matrix colnames `pwts' = `prn'
      ereturn matrix boot_tstud_pwsgs = `pwts'
      * the per-replication contrast standard errors go out compacted on the SAME row mask
      * as e(bootdist), so d_r = bootdist[ka]-bootdist[kb] pairs row-by-row with these
      tempname pwsed
      mata: wooldid_bootdist("`acc'", "`pwsed'", "`pwseacc'", "`pwsed'")
      cap matrix colnames `pwsed' = `prn'
      cap ereturn matrix bootdist_pwse = `pwsed'
      if `pwtmiss' > 0 {
        local ++wdbncaut
        local wdbc`wdbncaut' `pwtmiss' pairwise subgroup contrast row(s) had no usable studentized bootstrap distribution; their statistic, p-value and interval are left blank rather than filled from the coefficient-draw families
      }
    }
    if `wdbpbcamiss' > 0 & "`bdisp'" == "bca" {
      local ++wdbncaut
      local wdbc`wdbncaut' `wdbpbcamiss' pairwise subgroup contrast row(s) had no usable BCa family; their p-value and interval are left blank rather than filled from another family (the standard-error column is unaffected)
    }
    if `pwunmapped' > 0 {
      local ++wdbncaut
      local wdbc`wdbncaut' `pwunmapped' pairwise subgroup contrast row(s) could not be mapped to replication columns; their bootstrap inference is left missing
    }
    * disclosure, not a guard: under clustering at or above i the analytic pairwise
    * margins se is degenerate on some designs (missing, or ~1e-9 with p = 0, seen at HEAD
    * on a saturated i-as-cohort fixture), and boottstat studentizes by exactly that se
    * in every replication. boottstat + cluster() always nests i (the re-id leg requires it).
    * Gated on cluster() itself: nclus() is the observation count on an unclustered run
    if "`bts'" != "" & "`clus'" != "" {
      local ++wdbncaut
      local wdbc`wdbncaut' under clustering at or above the i level the analytic pairwise-contrast standard errors can be degenerate (the same limitation behind the refusal of jointtests with analytic clustered SEs), and boottstat studentizes each contrast by that standard error in every replication, so the studentized contrasts and the across-subgroup equality test can be absent or numerically unstable here; the non-studentized bootstrap (bootstrap() without boottstat) forms the contrasts from the coefficient draws instead
    }
  }

  * -------- the CI flavors, the percentile p and the BC z0 were computed and posted by
  * wooldid_bootwire, which runs before this program. They moved there because the
  * percentile package became a display default and the table rewiring needs the
  * percentile interval and its dual p in hand; see the comment at that call site

  forval cc = 1/`wdbncaut' {
    dis as text "Caution: `wdbc`cc''."
  }
end

*** F6: the permutation engine. builds the cluster->onset map, generates assignments
*** (mata: wooldid_permx), runs the rep loop as nested quiet wooldid calls on a
*** permuted onset column, and hands the summary statistics back to the shell via
*** c_local (the shell restores the observed results and posts them to e())
cap prog drop wooldid_permrun
prog def wooldid_permrun
version 16.0
  syntax , y(varname) i(varname) t(varname) ttrevar(varname) clus(varlist) permute(string) estimands(string) distname(string) [permreps(integer 0) PERMEXact permstrata(varname numeric) permprobs(string asis) PERMTStat repopts(string asis) wtype(string) wexp(string) ifr(string asis) compname(string) obsb(string) accname(string) seaccname(string) cmpbname(string) cmpsename(string) cellacc(string) cell0(string) dots(integer 1) asgname(string) rseed(string)]
  if `dots' == -1 local dots 1
  * rseed(): seeds the assignment draws (Mata's jumble/runiform share Stata's RNG
  * state), set before anything random runs; accepted under permute() since 2026-09-05
  if "`rseed'" != "" set seed `rseed'
  * F059 (2026-09-06): the RNG state before the first draw (e(perm_rngstate); permsave())
  local wdp_rng0 = c(rngstate)

  * comparison-statistic scalars, filled by the mandatory rep-0 run below (NFPERM-10)
  tempname wdpo1 wdpo2 wdpo3 wdpo4

  * one row per cluster: its onset level and (optional) stratum
  tempvar smp
  if `"`ifr'"' != "" {
    * NFPERM-03: exact-one rule, byte-parallel with permctl's insmp
    qui gen byte `smp' = ((`ifr') == 1)
    qui replace `smp' = 0 if `smp' >= .
  }
  else qui gen byte `smp' = 1
  qui replace `smp' = 0 if missing(`clus')
  tempname mpfr
  frame put `clus' `ttrevar' `permstrata' if `smp', into(`mpfr')
  frame `mpfr' {
    qui bysort `clus': keep if _n == 1
    if "`permstrata'" == "" {
      tempvar onestr
      qui gen byte `onestr' = 1
      local permstrata `onestr'
    }
    qui sort `permstrata' `clus'
    qui gen double __wdpassign = .
    * NFPERM-02: freeze the generation-time row identity. the assignment matrix is
    * generated in the (strata, cluster) sort order above, but frlink below RE-SORTS
    * this frame by the cluster key (frlink is sortpreserve for the master frame only),
    * so a positional _n read at application time attaches assignments to the WRONG
    * clusters whenever ascending cluster order crosses stratum order - silently wrong
    * stratified inference (probed: crossing-strata exact designs reported p=0 with the
    * observed assignment absent from their own enumeration). __wdprow keys each row to
    * its generation-time matrix row, making any later reorder harmless
    qui gen long __wdprow = _n
    * permsave(): the cluster values and observed onsets in generation-row order, so the
    * assignment block of the saved file can be labelled by cluster
    if "`asgname'" != "" {
      mata: st_local("wdp_clusvals", invtokens(strofreal(st_data(., "`clus'")', "%20.0g")))
      mata: st_local("wdp_obsons", invtokens(strofreal(st_data(., "`ttrevar'")', "%20.0g")))
    }
  }

  * assignments + counts (mata): fills the C x R matrix, sets wdp_R, the possible-
  * assignment count (double + string), and caution flags
  tempname wdpassign wdpP wdpL
  if "`asgname'" != "" local wdpassign `asgname'
  if "`permute'" == "probabilities" {
    * re-run the shared parser in engine mode: fills the probability matrix (rows in
    * ascending cluster order = the map frame's order) and the level vector
    frame `mpfr' {
      qui glevelsof `ttrevar', missing local(permonslvls)
      qui glevelsof `clus', local(permclusvals)
    }
    local permperr
    mata: wooldid_permp(st_local("permprobs"), "`ttrevar'", "`permonslvls'", "`permclusvals'", "`wdpP'", "`wdpL'")
    if "`permperr'" != "" {
      dis as error "`permperr'"
      error 198
    }
  }
  mata: wooldid_permx("`mpfr'", "`ttrevar'", "`permstrata'", `permreps', ("`permexact'" != ""), "`wdpassign'", "`wdpP'", "`wdpL'")
  if "`wdp_caut0'" != "" dis as text "Caution: `wdp_caut0'"
  if "`wdp_err'" != "" {
    dis as error "`wdp_err'"
    error 198
  }
  * NFPERM-09: conditionally degenerate designs (every stratum onset-homogeneous, or every
  * permprobs() row deterministic) pass the global freedom precheck but admit exactly one
  * possible assignment; they formerly ran to completion reporting a useless p (p=1, or
  * p=1/(R+1) when a deterministic table targets a non-observed assignment). both
  * generators report the support count in wdp_npd, and M==1 iff lnM==0 exactly (the
  * per-stratum lnfactorial terms cancel exactly; smallest nondegenerate M is 2), so the
  * numeric compare cannot false-fire. placed BEFORE the reps-vs-support cautions so a
  * degenerate run is refused rather than first advised to consider permexact
  if real("`wdp_npd'") == 1 {
    dis as error "Error: the requested randomization design admits only one possible assignment (every stratum is onset-homogeneous, or every permprobs() row is deterministic); there is nothing to randomize."
    error 459
  }
  if "`wdp_caut1'" != "" dis as text "Caution: `wdp_caut1'"
  if "`wdp_caut2'" != "" dis as text "Caution: `wdp_caut2'"

  * link the map into the data once; each rep refreshes the onset column from it
  tempvar lnk repons
  qui frlink m:1 `clus', frame(`mpfr') generate(`lnk')
  local ndone 0
  local nfail 0
  local nblank 0
  tempname res
  matrix `res' = J(`wdp_R', `: word count `estimands'', .)
  local wtif
  if "`wtype'" != "" local wtif [`wtype'=`wexp']
  if `"`ifr'"' != "" local wtif `wtif' if (`ifr')
  local repcall `y' `i' `t'
  * statistic-path consistency: draws are estimated under fastnose, whose point
  * estimates differ from the full run's in the last digits (the F9 contract), so the
  * comparison statistic comes from a rep-0 fastnose run of the TRUE assignment - this
  * keeps the exact test's observed-is-a-member property bitwise. the display still
  * shows the full run's estimate. the statistics are held in tempname SCALARS: a
  * local's decimal round-trip re-parses one ulp off, which breaks the
  * observed-is-a-member property (measured: b >= local-copy-of-b evaluates FALSE).
  * Since the 2026-09-05 fast-mode default the observed run itself takes the fast path
  * unless permtstat or permanalyticse, and the comparator then equals it bitwise; the
  * run is kept unconditional (one extra fit) because under permanalyticse it is still
  * load-bearing
  * ---- THE NESTED-REPLICATION IDIOM (this shape recurs at both sites here and at six in
  * wooldid_bootrun; it is explained once, here). Three parts, each load-bearing:
  *  (1) `nobreak' is what makes WOOLDID_REPMODE un-leakable. It is a GLOBAL, and an
  *      interrupt landing between the set and the clear would carry it into the user's
  *      session, where wooldid_worker consumes it read-once and would silently skip the
  *      diagnostic e() returns of the user's NEXT run.
  *  (2) `_rc' is stashed into a local on the very next line because any program called in
  *      between resets it. MEASURED: with `prog def pp' whose body is `cap error 0',
  *      `cap error 459' followed by `pp' leaves the caller at _rc == 0. (Clearing a
  *      global does not reset it - also measured - so the stash is merely defensive here,
  *      and genuinely load-bearing in wooldid_bootrun, where a `dis' and a frame block
  *      sit between the capture and the read.) `cap error `wdrepdone'' re-raises it so
  *      `if _rc == 1' can pass a user break through while other failures become `badrep'.
  *  (3) What a missed rc costs: nothing in wooldid clears e() on failure, so a failed
  *      replication leaves the PREVIOUS replication's e(b) in place. Without `badrep' the
  *      grab would succeed and write that draw into row r a second time - a duplicated
  *      point in the reference distribution, with no error and no `x' printed
  nobreak {
    global WOOLDID_REPMODE 1
  cap qui wooldid `repcall' `ttrevar' `wtif' , `repopts'
    local wdrepdone = _rc
    global WOOLDID_REPMODE
  }
  cap error `wdrepdone'
  * NFPERM-10: the comparator fit and its scalars must be valid. formerly a rep-0
  * failure fell through silently (retaining full-run-provenance scalars, breaking the
  * exact test's bitwise observed-is-a-member property), and a missing comparison
  * statistic (e.g. a degenerate observed SE under permtstat) silently zeroed every
  * tail count - exact p = 0 with no message. the shell restores the observed
  * estimation results bit-intact on these error paths
  if _rc == 1 error 1
  if _rc != 0 {
    dis as error "Error: the rep-0 comparison run of the observed assignment failed (rc `=_rc'); permutation inference is impossible. The observed estimation results are unaffected."
    error 459
  }
  local ec 0
  foreach ee in `estimands' {
    local ++ec
    if "`permtstat'" == "" cap scalar `wdpo`ec'' = _b[main_`ee']
    else cap scalar `wdpo`ec'' = _b[main_`ee']/_se[main_`ee']
    if _rc != 0 {
      dis as error "Error: the observed comparison statistic for `ee' could not be read from the rep-0 run; permutation inference is impossible."
      error 459
    }
    if missing(`wdpo`ec'') {
      dis as error "Error: the observed comparison statistic for `ee' is missing (unidentified observed estimate or standard error); permutation inference is impossible."
      error 459
    }
  }
  * NFPERM-12 (lean kernel): export the comparison statistics actually used for the
  * tail counts, guarded valid above; the caller posts them as e(perm_comparator)
  if "`compname'" != "" {
    matrix `compname' = J(1, `: word count `estimands'', .)
    local ec 0
    foreach ee in `estimands' {
      local ++ec
      matrix `compname'[1, `ec'] = `wdpo`ec''
    }
    matrix colnames `compname' = `estimands'
  }
  * chunk 1 (2026-09-05): the FULL comparator vector and, name-aligned to the observed
  * e(b) column set, the R x K accumulator every reported surface is built from. Under
  * permtstat a parallel R x K of replication standard errors (the rep's own e(V)
  * diagonal) and the comparator's own standard errors
  local wdp_K 0
  if "`accname'" != "" {
    matrix `cmpbname' = e(b)
    local wdpnames0 : colfullnames `obsb'
    local wdpnames1 : colfullnames e(b)
    if "`wdpnames0'" != "`wdpnames1'" {
      dis as error "Error: the rep-0 coefficient vector does not match the observed run's column set; permutation accumulation would be misaligned. This is unexpected - please report it."
      error 459
    }
    local wdp_K = colsof(`obsb')
    matrix `accname' = J(`wdp_R', `wdp_K', .)
    matrix colnames `accname' = `wdpnames0'
    if "`seaccname'" != "" {
      matrix `seaccname' = J(`wdp_R', `wdp_K', .)
      matrix colnames `seaccname' = `wdpnames0'
      matrix `cmpsename' = J(1, `wdp_K', .)
      matrix colnames `cmpsename' = `wdpnames0'
      cap confirm matrix e(V)
      if _rc == 0 mata: wooldid_permcse("`cmpsename'", `wdp_K')
    }
  }
  * cell accumulators for the continuous-treatment (cfx) matrices, keyed by row name from
  * the rep-0 run (the same lattice levels the observed run posted), plus a one-row
  * comparator harvest of the rep-0 cells themselves
  local wdp_ncfx 0
  if "`cellacc'" != "" {
    foreach cm in `e(wdid_cfxmats)' {
      cap confirm matrix e(`cm')
      if _rc != 0 continue
      local ++wdp_ncfx
      mata: wooldid_bootcellkeys("e(`cm')", "`cellacc'_cfx`wdp_ncfx'", `wdp_R', 0)
      mata: wooldid_bootcellkeys("e(`cm')", "`cell0'_cfx`wdp_ncfx'", 1, 0)
      if "`seaccname'" != "" {
        mata: wooldid_bootcellkeys("e(`cm')", "`cellacc'_cfx`wdp_ncfx'se", `wdp_R', 0)
        mata: wooldid_bootcellkeys("e(`cm')", "`cell0'_cfx`wdp_ncfx'se", 1, 0)
      }
      local wdb_cfxname`wdp_ncfx' `cm'
    }
    if `wdp_ncfx' > 0 mata: wooldid_bootcells("`cell0'", 1, `wdp_ncfx', ("`seaccname'" != ""))
  }
  c_local wdp_ncfx `wdp_ncfx'
  forval kk = 1/`wdp_ncfx' {
    c_local wdb_cfxname`kk' `wdb_cfxname`kk''
  }
  if `dots' > 0 dis as text "Permutation replications (`wdp_R'): " _continue
  forval r = 1/`wdp_R' {
    * F066 (2026-09-06): time replications 1-5 and print a projection (Stata timer 100)
    if `r' == 1 & `wdp_R' > 5 {
      timer clear 100
      timer on 100
    }
    if `r' == 6 {
      timer off 100
      qui timer list 100
      local wdp_proj = r(t100)/5*`wdp_R'
      timer clear 100
      local wdp_ptxt = cond(`wdp_proj' < 90, strofreal(round(`wdp_proj')) + " s", strofreal(round(`wdp_proj'/60)) + " min")
      if `dots' > 0 dis as text " (~`wdp_ptxt' projected for `wdp_R' replications)" _continue
    }
    frame `mpfr': qui replace __wdpassign = el("`wdpassign'", __wdprow, `r')
    cap drop `repons'
    qui frget `repons' = __wdpassign, from(`lnk')
    nobreak {
      global WOOLDID_REPMODE 1
    cap qui wooldid `repcall' `repons' `wtif' , `repopts'
      local wdrepdone = _rc
      global WOOLDID_REPMODE
    }
    cap error `wdrepdone'
    if _rc == 1 {
      if `dots' > 0 dis as result " "
      error 1
    }
    local badrep 0
    if _rc != 0 local badrep 1
    if `badrep' == 0 {
      local ec 0
      foreach ee in `estimands' {
        local ++ec
        if "`permtstat'" == "" {
          cap matrix `res'[`r', `ec'] = _b[main_`ee']
          * NFPERM-10: mirror the t-branch missing check so whole-row blanking keeps a
          * single completed-row set (nuse == ndone for every column)
          if _rc == 0 {
            if missing(`res'[`r', `ec']) local badrep 1
          }
        }
        else {
          cap matrix `res'[`r', `ec'] = _b[main_`ee']/_se[main_`ee']
          if _rc == 0 {
            if missing(`res'[`r', `ec']) local badrep 1
          }
        }
        if _rc != 0 local badrep 1
      }
    }
    if `badrep' == 0 {
      local ++ndone
      if "`accname'" != "" {
        * name-aligned grab of the full rep e(b) into row r (ame-family structural zeros
        * blanked, A12(b) - a voided ame replication is |0| < |obs| and would otherwise
        * count as a non-exceedance); under permtstat also this rep's own standard errors
        local wdb_gblank 0
        mata: wooldid_bootgrab("`accname'", `r', "`obsb'", "`seaccname'", "")
        local nblank = `nblank' + `wdb_gblank'
        if `wdp_ncfx' > 0 mata: wooldid_bootcells("`cellacc'", `r', `wdp_ncfx', ("`seaccname'" != ""))
      }
    }
    else {
      local ++nfail
      * blank the whole row: a partially-failed rep (reachable under permtstat with a
      * degenerate SE on one of several estimands) must not leak its good columns into
      * the per-estimand summaries, and the permdist packing keys on column 1
      local ec 0
      foreach ee in `estimands' {
        local ++ec
        cap matrix `res'[`r', `ec'] = .
      }
    }
    if `dots' > 0 {
      * one mark every dots() replications - the documented cadence (a failed
      * replication always prints its x) - with the count after every 50 marks
      if `badrep' dis as error "x" _continue
      else if mod(`r', `dots') == 0 dis as text "." _continue
      if mod(`r', 50*`dots') == 0 dis as text " " `r' _newline "  " _continue
    }
  }
  if `dots' > 0 dis as text " done"
  cap drop `lnk'
  cap drop `repons'
  frame drop `mpfr'

  * the 20 is the same floor wooldid_permctl enforces on permreps() itself ("the engine
  * refuses inference on fewer than 20 completed replications"), so the two numbers must
  * move together; the second floor is half the replications the user asked for.
  * permexact is exempt BY RULING, not by omission: under an enumeration the observed
  * assignment is still a member of its own reference set, so the shares stay exact over
  * whatever completed and what incompleteness costs is RESOLUTION, recorded per element
  * in e(perm_n) - see the BS-23 note below. A sampled run that lost half its draws has
  * no such guarantee, which is why it is refused instead
  if `ndone' < max(20, `wdp_R'/2) & "`permexact'" == "" {
    dis as error "Error: only `ndone' of `wdp_R' permutation replications estimated successfully; too few for inference. Check whether reassigned onset configurations are frequently inestimable on these data."
    error 459
  }
  * NFPERM-01 (lean residue): under permexact ANY failure means the enumeration is
  * incomplete, so the caution fires at nfail > 0 with wording that states the
  * substance; sampled mode keeps the 10% threshold (adjudicated accounting design)
  * BS-23 (author ruling 2026-09-07): what incompleteness costs is RESOLUTION, not the
  * convention. The observed assignment is one of the enumerated columns, so the shares stay
  * exact - over the assignments that completed. The old wording ("shares over completed
  * assignments only ... the observed-is-a-member and 2/M-floor readings no longer apply")
  * contradicted the per-element caution wooldid_permwire prints, and described a rule the
  * code did not use.
  if `nfail' > 0 & "`permexact'" != "" dis as text "Caution: `nfail' of `wdp_R' enumerated assignments failed to estimate, so the enumeration is incomplete. The p-values stay exact shares - the observed assignment is still a member of its own reference set - but over the `=`wdp_R' - `nfail'' assignments that completed: the smallest attainable p is 2/n under the default convention (1/n under the absolute-value and one-sided conventions) rather than 2/`wdp_R', with each element's n in e(perm_n)."
  * (a failed assignment blanks its WHOLE accumulator row, so every element loses it and
  * n = M - nfail across the board; per-element loss on an assignment that DID estimate is
  * what wooldid_permwire's e(perm_incomplete) caution covers.)
  else if `nfail' > 0 & `nfail' > 0.1*`wdp_R' dis as text "Caution: `nfail' of `wdp_R' permutation replications failed to estimate and were excluded."

  * dense legacy permdist matrix (completed rows only, keyed on the legacy estimands),
  * built once from res. The per-element statistics now come from the name-aligned
  * accumulator (wooldid_permwire); this matrix keeps the documented estimand-keyed shape
  local rr 0
  matrix `distname' = J(max(`ndone', 1), `: word count `estimands'', .)
  forval r = 1/`wdp_R' {
    if !missing(el("`res'", `r', 1)) {
      local ++rr
      local ec 0
      foreach ee in `estimands' {
        local ++ec
        matrix `distname'[`rr', `ec'] = el("`res'", `r', `ec')
      }
    }
  }
  matrix colnames `distname' = `estimands'

  c_local wdp_r_reps `wdp_R'
  c_local wdp_r_done `ndone'
  c_local wdp_r_fail `nfail'
  c_local wdp_r_blank `nblank'
  c_local wdp_r_npd `wdp_npd'
  c_local wdp_r_nps `wdp_nps'
  c_local wdp_r_rng0 `wdp_rng0'
  * BS-23 (author ruling 2026-09-07): the accumulator row carrying the observed
  * assignment (0 outside an enumeration); wooldid_permx posts it
  c_local wdp_r_obsrow `wdp_obsrow'
  if "`asgname'" != "" {
    c_local wdp_clusvals `wdp_clusvals'
    c_local wdp_obsons `wdp_obsons'
  }
end

cap prog drop wooldid_permwire
prog def wooldid_permwire, eclass
version 16.0
*** chunk 1 (2026-09-05): rewires the reported inference surfaces to the permutation
*** distribution. Point estimates never change. Every p-value column becomes the
*** equal-tailed permutation p (coverage by subtraction: an element with no usable
*** permutation distribution is BLANK, never analytic), every interval the shift
*** interval, the statistic column the permuted t-ratio under permtstat and blank
*** otherwise; the standard-error column is left as the worker posted it - the analytic
*** standard error (labelled se_analytic) under permtstat or permanalyticse, blank
*** under the fast-mode default. The per-element kernel and the stripe-aligned vectors
*** are posted by name; the analytic tables are kept as hidden twins
  syntax , acc(string) obsb(string) cmpb(string) [seacc(string) cmpse(string) level(real 0) m(string) orow(integer 0) done(integer 0) fail(integer 0) blank(integer 0) ts(integer 0) method(string) reps(integer 0) ase(string) clus(varlist) gt(string) absp(string) side(string) rawp(string)]
  if `level' == 0 local level = c(level)
  * p-value convention code (F052/F064, 2026-09-06): 0 equal-tailed (default), 1 absolute
  * add-one (permabsp), 2 raw |s| share c/n without add-one (permrawp: the ritest / Stata
  * permute number), 3 left one-sided, 4 right one-sided (permside())
  * this expression is duplicated byte-for-byte in wooldid_permctl, where it picks the
  * floor the resolution caution prints; the two must be kept in sync
  local wdpabs = cond("`rawp'" != "", 2, cond("`side'" == "left", 3, cond("`side'" == "right", 4, ("`absp'" != ""))))
  local wdpM = cond("`m'" == "" | "`m'" == ".", ".", "`m'")
  if `ts' == 0 {
    local seacc
    local cmpse
  }

  * -------- column-position walk: e(b) col k <-> results-matrix row (positional: e(b)
  * is built by stacking wdidpooled + the SG blocks in sglist order). excluded columns =
  * rows whose reported ESTIMATE is missing. The analytic se column is captured in the
  * same walk (hidden e(perm_ase_obs); permsave() writes it as the observed se_* row)
  tempname wm excl aseobs
  matrix `wm' = e(wdidmainresults)
  local nmain = rowsof(`wm')
  local K = colsof(`obsb')
  matrix `excl' = J(1, `K', 0)
  matrix `aseobs' = J(1, `K', .)
  local sgl `e(wdid_sglist)'
  local kk 0
  forval rr = 1/`nmain' {
    local ++kk
    if `wm'[`rr', 1] >= . matrix `excl'[1, `kk'] = 1
    matrix `aseobs'[1, `kk'] = `wm'[`rr', 2]
  }
  foreach sg of local sgl {
    tempname wsg`sg'
    cap confirm matrix e(wdidsgresults_`sg')
    if _rc != 0 continue
    matrix `wsg`sg'' = e(wdidsgresults_`sg')
    forval rr = 1/`=rowsof(`wsg`sg'')' {
      local ++kk
      if `wsg`sg''[`rr', 1] >= . matrix `excl'[1, `kk'] = 1
      matrix `aseobs'[1, `kk'] = `wsg`sg''[`rr', 2]
    }
  }
  if `kk' != `K' {
    dis as error "Error: permutation rewiring found `kk' results-matrix rows against `K' coefficient columns; the inference surfaces cannot be aligned. This is unexpected - please report it."
    error 459
  }

  * -------- the per-element kernel (mata): n, s_obs, p, quantiles, sd, shift interval,
  * shift se, null band at level(), live; permabsp switches the p to the absolute-value
  * convention and the interval to its symmetric dual
  tempname kern
  local wdp_incomplete
  local wdp_addone
  mata: wooldid_permstat("`acc'", "`seacc'", "`obsb'", "`cmpb'", "`cmpse'", "`excl'", `level', `wdpM', `orow', `wdpabs', "`kern'")
  local obn : colfullnames `obsb'

  * -------- rewire the results matrices: col 3 = permuted t-ratio (permtstat) or blank,
  * col 4 = permutation p, cols 5-6 = shift interval. col 2 stays as posted (analytic se
  * or blank). Blanked rows are counted for the cautions
  local wdpmiss 0
  local wdplow 0
  local wdpdeg
  * NOTE `wdptsmiss' is currently incremented (once per coefficient blanked under
  * permtstat for want of a usable observed analytic se) and never read: the caution at
  * the foot of this program states both causes in one sentence without being able to say
  * which applies. The two are different fixes for the user - too few completed
  * replications says raise permreps(), no usable observed se says the element has no
  * identified standard error to studentize by - so the counter is kept for a caution
  * split, not deleted
  local wdptsmiss 0
  local wdpoff
  tempname wma
  matrix `wma' = `wm'
  local kk 0
  * KERNEL ROW MAP, used by both rewiring passes below (defined by wooldid_permstat and
  * extended by wooldid_permtail): 1 n, 2 s_obs, 3 p, 4 qlo, 5 qhi, 6 sd, 7 lb, 8 ub,
  * 9 se_shift, 10 offcentre flag, 11 bandlo, 12 bandhi, 13 live, 14-18 the tail block.
  * Row 13 gates the whole rewrite (the element has a usable distribution at all); 2, 3
  * and 7/8 fill the table; 1, 6 and 10 drive the cautions. These subscripts are a
  * POSITIONAL contract with `out = J(18, K, .)' in wooldid_permtail: insert a row there
  * and every subscript here shifts, and because rows 1-13 all hold finite numbers the
  * result is a silently mislabelled table (the p column filled from qlo, say), not an
  * error. The rownames that would decode it are only attached far below, at
  * e(perm_kernel)
  forval rr = 1/`nmain' {
    local ++kk
    matrix `wm'[`rr', 3] = .
    if `kern'[13, `kk'] == 1 {
      if `ts' matrix `wm'[`rr', 3] = `kern'[2, `kk']
      matrix `wm'[`rr', 4] = `kern'[3, `kk']
      matrix `wm'[`rr', 5] = `kern'[7, `kk']
      matrix `wm'[`rr', 6] = `kern'[8, `kk']
      if `kern'[1, `kk'] < 0.5*`done' local ++wdplow
      if `kern'[6, `kk'] == 0 local wdpdeg `wdpdeg' `: word `kk' of `obn''
      if `kern'[10, `kk'] == 1 local wdpoff `wdpoff' `: word `kk' of `obn''
    }
    else {
      forval cc = 4/6 {
        matrix `wm'[`rr', `cc'] = .
      }
      if `excl'[1, `kk'] == 0 {
        local ++wdpmiss
        if `ts' & (`aseobs'[1, `kk'] >= . | `aseobs'[1, `kk'] <= 0) local ++wdptsmiss
      }
    }
  }
  local wdplvltxt = cond(mod(`level',1)==0, strofreal(`level'), subinstr(strofreal(`level',"%9.0g"), ".", "p", .))
  local wmc : colnames `wm'
  local wmc : subinstr local wmc "se" "se_analytic", word
  local wmc : subinstr local wmc "t" "t_perm", word
  local wmc : subinstr local wmc "z" "t_perm", word
  local wmc : subinstr local wmc "p" "p_perm", word
  * NOTE the renames above are the load-bearing ones: se -> se_analytic and t/z -> t_perm
  * are the only signal to a downstream reader that column 2 is the analytic standard error
  * and column 3 the permuted t-ratio. A pair of lb_95ci/ub_95ci renames used to sit here
  * and was DELETED as unreachable (D-5, author ruling 2026-09-11): the worker names the
  * interval columns lb_<level>ci from the run's own level and wooldid_permctl hands this
  * program that same level, so at 95 the names are already lb_95ci and at any other level
  * they already carry it - the subinstr could never match
  matrix colnames `wm' = `wmc'
  ereturn hidden matrix perm_analytic_mainresults = `wma'
  ereturn matrix wdidmainresults = `wm'
  foreach sg of local sgl {
    * PASS 2 of the e(b) stripe walk. `kk' is deliberately NOT reset here: it indexes
    * e(b)'s columns in the order they were stacked (the pooled block, then one block per
    * value of e(wdid_sglist)), so it must run monotonically across the main-table loop
    * and this one in BOTH passes. This pass skips a subgroup on the tempname matrix
    * existing, which is a PROXY for pass 1's test on e(wdidsgresults_`sg') - the two
    * agree by construction, not by contract, and the `kk' != `K'' guard that would catch
    * a mis-walk runs after pass 1 only. If they ever desync, every element from the first
    * skipped subgroup onwards takes another element's kernel column: the p-value, shift
    * interval and t-ratio are attached to the wrong coefficient, with no error and no
    * caution
    cap confirm matrix `wsg`sg''
    if _rc != 0 continue
    tempname wsa`sg'
    matrix `wsa`sg'' = `wsg`sg''
    forval rr = 1/`=rowsof(`wsg`sg'')' {
      local ++kk
      matrix `wsg`sg''[`rr', 3] = .
      if `kern'[13, `kk'] == 1 {
        if `ts' matrix `wsg`sg''[`rr', 3] = `kern'[2, `kk']
        matrix `wsg`sg''[`rr', 4] = `kern'[3, `kk']
        matrix `wsg`sg''[`rr', 5] = `kern'[7, `kk']
        matrix `wsg`sg''[`rr', 6] = `kern'[8, `kk']
        if `kern'[1, `kk'] < 0.5*`done' local ++wdplow
        if `kern'[6, `kk'] == 0 local wdpdeg `wdpdeg' `: word `kk' of `obn''
        if `kern'[10, `kk'] == 1 local wdpoff `wdpoff' `: word `kk' of `obn''
      }
      else {
        forval cc = 4/6 {
          matrix `wsg`sg''[`rr', `cc'] = .
        }
        if `excl'[1, `kk'] == 0 {
          local ++wdpmiss
          if `ts' & (`aseobs'[1, `kk'] >= . | `aseobs'[1, `kk'] <= 0) local ++wdptsmiss
        }
      }
    }
    local wsc : colnames `wsg`sg''
    local wsc : subinstr local wsc "se" "se_analytic", word
    local wsc : subinstr local wsc "t" "t_perm", word
    local wsc : subinstr local wsc "z" "t_perm", word
    local wsc : subinstr local wsc "p" "p_perm", word
    matrix colnames `wsg`sg'' = `wsc'
    * the hidden twin is skipped when its name would exceed 32 characters (a 9+ digit
    * subgroup value); the visible table is unaffected
    if strlen("perm_analytic_sgresults_`sg'") <= 32 ereturn hidden matrix perm_analytic_sgresults_`sg' = `wsa`sg''
    ereturn matrix wdidsgresults_`sg' = `wsg`sg''
  }

  * -------- gteffects: the cell keys are cohort-by-onset-by-relative-period, and a
  * cohort's onset is exactly what a replication reassigns, so no cell has a permutation
  * distribution; the stage-1 inference columns are blanked (never analytic under a
  * permutation header) and the display note says so
  if "`gt'" != "" {
    cap confirm matrix e(gteffects)
    if _rc == 0 {
      tempname gtm
      matrix `gtm' = e(gteffects)
      local gtcn : colnames `gtm'
      local gtrn : rowfullnames `gtm'
      * columns 6-9 are se, p, lb, ub; column 5 is the estimate and is deliberately kept
      * (point estimates do not change under permute()). The layout is defined in
      * wooldid_bootcellsx ("gteffects: cols estimate(5) se(6) p(7) lb(8) ub(9)") and the
      * two must move together: a layout change silently blanks the estimate column or,
      * worse, leaves analytic se/p/ci standing under a permutation header - the exact
      * outcome the block comment above says must never happen. `gtcn'/`gtrn' are captured
      * and re-applied because st_matrix() writes a bare matrix with no stripes. NOTE this
      * statement leaves a global Mata object behind in the session; it is named
      * wooldid_Mtmp rather than `M' so it cannot collide with the user's own (D-3)
      mata: wooldid_Mtmp = st_matrix("`gtm'"); wooldid_Mtmp[., 6..9] = J(rows(wooldid_Mtmp), 4, .); st_matrix("`gtm'", wooldid_Mtmp)
      matrix colnames `gtm' = `gtcn'
      matrix rownames `gtm' = `gtrn'
      ereturn matrix gteffects = `gtm'
      ereturn scalar perm_gt_blank = 1
    }
  }

  * -------- the stripe-aligned surfaces, all colnamed like e(b)
  tempname pv ci sesh sdv nv nb cmp2 ase2 kern2
  matrix `pv' = `kern'[3, 1...]
  matrix colnames `pv' = `obn'
  ereturn matrix perm_p = `pv'
  matrix `ci' = `kern'[7..8, 1...]
  matrix colnames `ci' = `obn'
  matrix rownames `ci' = ll ul
  ereturn matrix perm_ci = `ci'
  matrix `sesh' = `kern'[9, 1...]
  matrix colnames `sesh' = `obn'
  ereturn matrix perm_se_shift = `sesh'
  matrix `sdv' = `kern'[6, 1...]
  matrix colnames `sdv' = `obn'
  ereturn matrix perm_sd = `sdv'
  matrix `nv' = `kern'[1, 1...]
  matrix colnames `nv' = `obn'
  ereturn matrix perm_n = `nv'
  matrix `nb' = `kern'[11..12, 1...]
  matrix colnames `nb' = `obn'
  matrix rownames `nb' = lo hi
  ereturn matrix perm_nullband = `nb'
  * F064 (2026-09-06): the one-sided shares, the Monte-Carlo se of the reported p and its
  * level() Clopper-Pearson interval, by name
  tempname ptl pth pse pci
  matrix `ptl' = `kern'[14, 1...]
  matrix colnames `ptl' = `obn'
  ereturn matrix perm_p_lo = `ptl'
  matrix `pth' = `kern'[15, 1...]
  matrix colnames `pth' = `obn'
  ereturn matrix perm_p_hi = `pth'
  matrix `pse' = `kern'[16, 1...]
  matrix colnames `pse' = `obn'
  ereturn matrix perm_p_se = `pse'
  matrix `pci' = `kern'[17..18, 1...]
  matrix colnames `pci' = `obn'
  matrix rownames `pci' = lo hi
  ereturn matrix perm_p_ci = `pci'
  matrix `cmp2' = `cmpb'
  matrix colnames `cmp2' = `obn'
  ereturn matrix perm_comparator_b = `cmp2'
  if `ts' {
    tempname cmps2
    matrix `cmps2' = `cmpse'
    matrix colnames `cmps2' = `obn'
    ereturn matrix perm_comparator_se = `cmps2'
  }
  matrix `kern2' = `kern'
  matrix rownames `kern2' = n s_obs p qlo qhi sd lb ub se_shift offcentre bandlo bandhi live p_lo p_hi se_p cp_lo cp_hi
  matrix colnames `kern2' = `obn'
  ereturn hidden matrix perm_kernel = `kern2'
  matrix `ase2' = `aseobs'
  matrix colnames `ase2' = `obn'
  ereturn hidden matrix perm_ase_obs = `ase2'
  * above 10,000 completed replications the draw matrices are not posted AT ALL (not
  * truncated) - the same cap wooldid_bootwire applies to e(bootdist), and documented in
  * the help. permsave() writes from the accumulators and is NOT subject to it, so a user
  * who needs the draws from a 50,000-replication run has that route. The `cap's below
  * mean a post that fails for any other reason also passes without a message
  if `done' <= 10000 {
    tempname pd pds
    if `ts' {
      mata: wooldid_bootdist("`acc'", "`pd'", "`seacc'", "`pds'")
      cap matrix colnames `pds' = `obn'
      cap ereturn matrix permdist_se = `pds'
    }
    else mata: wooldid_bootdist("`acc'", "`pd'")
    cap matrix colnames `pd' = `obn'
    cap ereturn matrix permdist_b = `pd'
  }
  ereturn scalar perm_level = `level'
  ereturn scalar perm_exact = (`wdpM' < .)
  ereturn local perm_incomplete `wdp_incomplete'
  * BS-23 (author ruling 2026-09-07): the elements that actually take the add-one form -
  * their reference set does not contain the observed statistic
  ereturn local perm_addone `wdp_addone'
  ereturn local perm_offcentre `wdpoff'
  ereturn scalar perm_blank = `blank'
  if "`ase'" != "" ereturn scalar perm_analyticse = 1
  ereturn local perm_pconv = cond(`wdpabs' == 1, "absolute", cond(`wdpabs' == 2, "raw", cond(`wdpabs' == 3, "left", cond(`wdpabs' == 4, "right", "equaltailed"))))
  ereturn scalar perm_absp = (`wdpabs' == 1)

  * -------- header: the display reads e(wdid_vcespec) and e(wdid_fastnose)
  local wdpcl = cond("`clus'" != "", "cluster(`clus')", "unclustered")
  if `ts' local wdphdr randomization inference on t-ratios (permtstat; permute(`method'), `reps' replications); se_analytic = the analytic `wdpcl' standard error each t-ratio is formed with
  else if "`ase'" != "" local wdphdr randomization inference (permute(`method'), `reps' replications); se_analytic = the analytic `wdpcl' standard error, shown for reference only
  else local wdphdr randomization inference (permute(`method'), `reps' replications); analytic standard errors not computed (permanalyticse restores them)
  ereturn hidden local wdid_vcespec "`wdphdr'"
  * CS-4: the permutation engine's covariance state. Under the fast-mode default no e(V)
  * is posted at all; under permanalyticse the posted matrix is the ANALYTIC one, shown for
  * reference, and is not a permutation quantity.
  cap confirm matrix e(V)
  local wdphasv = (_rc == 0)
  * PS-62 (2026-09-07): under permute() the posted covariance is still the ANALYTIC one the
  * worker built - this engine reposts none of its own - so if the worker's event-study
  * splice landed, e(V) really does carry blocks and the label must keep saying so
  local wdpvs0 "`e(wdid_vcvstate)'"
  ereturn hidden local wdid_vcvstate = cond(`wdphasv', cond("`wdpvs0'" == "analytic-blocks", "analytic-blocks", "analytic-diagonal"), "absent")
  * PS-1: with no covariance posted the conventional variance-type names describe nothing,
  * so clear them rather than leave the analytic run's labels beside an absent e(V)
  if !`wdphasv' {
    ereturn local vce
    ereturn local vcetype
  }
  * PS-46 (author ruling 2026-09-10, on the measurement at
  * review_process/research_2026-09-10/permute_fastmode_esttab.md). With no e(V) posted,
  * esttab and estout do NOT refuse: they print an EMPTY inference row beneath a footer
  * asserting the numbers are there, and write that assertion into .tex and .csv. estout
  * reads posted e(p)/e(ci_l)/e(ci_u) AHEAD of anything derived from e(V) (its own bVs2
  * list), so posting them under those derivable names makes a coefficient table carry THIS
  * run's real permutation p-values and shift intervals - which is the point: a user should
  * be able to table the coefficients with the standard tools even where there is no
  * standard error to report.
  * Nothing is invented. No covariance is posted, so test, lincom, _se[] and estat vce go on
  * refusing exactly as before; and no e(se) is posted, because a permutation run HAS no
  * standard error - the shift se is a different object and naming it se would be the
  * fabrication this route exists to avoid.
  * Fast mode only: under permanalyticse an e(V) exists and the tools already work from it.
  * Read back from e() rather than reused: `pv' and `ci' were MOVED by their own posts.
  if !`wdphasv' {
    tempname psp psc psl psu
    cap confirm matrix e(perm_p)
    if _rc == 0 {
      matrix `psp' = e(perm_p)
      * e(perm_p) is row 3 of the kernel and carries its "r3" rowname; a 1 x K matrix posted
      * under the name p should say p, since a user reads it with matrix list
      matrix rownames `psp' = p
      ereturn matrix p = `psp'
    }
    cap confirm matrix e(perm_ci)
    if _rc == 0 {
      matrix `psc' = e(perm_ci)
      if rowsof(`psc') >= 2 {
        matrix `psl' = `psc'[1, 1...]
        matrix `psu' = `psc'[2, 1...]
        ereturn matrix ci_l = `psl'
        ereturn matrix ci_u = `psu'
      }
    }
  }
  ereturn hidden local wdid_fastnose ""

  * -------- cautions (estimation-time; the footer carries the standing disclosures)
  * D-6 (author ruling 2026-09-11): the two causes are SPLIT, because they are different
  * fixes for the user - too few completed replications says raise permreps(), while no
  * usable observed analytic standard error says that element has no identified SE to
  * studentize by, and no number of replications will change it. `wdptsmiss' counts the
  * second; it was incremented and never read before this
  local wdpmissrep = `wdpmiss' - `wdptsmiss'
  if `wdpmissrep' < 0 local wdpmissrep 0
  if `wdpmissrep' > 0 dis as text "Caution: `wdpmissrep' reported coefficient(s) have no usable permutation distribution (fewer than 2 completed replications carry them); their p-value and interval are left blank rather than filled from any other source. More replications would raise this."
  if `wdptsmiss' > 0 dis as text "Caution: `wdptsmiss' reported coefficient(s) could not be studentized under permtstat - their observed analytic standard error is missing or non-positive, so there is no t-ratio to permute; their p-value and interval are left blank rather than filled from any other source. This is a property of the specification, not of permreps(): more replications will not change it."
  if `wdplow' > 0 dis as text "Caution: `wdplow' reported coefficient(s) were present in fewer than half of the completed replications (a reassigned onset can make a horizon or a subgroup level inestimable in a replication); their permutation p-values rest on the counts in e(perm_n)."
  if "`wdpdeg'" != "" dis as text "Caution: the permutation distribution is degenerate (every completed draw equal) for: `wdpdeg'; those permutation p-values are uninformative."
  * The worked example this caution used to carry - the ame rows under lattice() with
  * contreatcontrols(basic), "where the baseline dose slope is not partialled out", plus a
  * pointer to "the standing issue on the reference cell" - was removed on 2026-09-11. SI-20
  * FIXED that condition on 2026-09-05 (1df8613: basic now always absorbs i.lattice##c.(ringer),
  * and the measured reference margin went from 0.27..2.00 to exactly 0 on all 59 affected
  * configs), and the ticket moved to wooldid_solved_issues.md, so the caution was sending a
  * user after a cause that no longer exists and a section that no longer stands. The mechanism
  * is what matters here and is kept; resist re-attaching a worked example unless it is one you
  * have just measured, since an example in a caution outlives the defect that motivated it
  if "`wdpoff'" != "" dis as text "Caution: the permutation distribution is not centred at zero (its central `level'% band excludes 0) for: `wdpoff'. The shift interval takes the estimate's sampling error to have that null distribution, so it is not interpretable there and is left blank; the permutation p-value stands. This happens when a quantity is not zero under reassigned onsets, or under a permute(probabilities) design far from the observed assignment."
  if `blank' > 0 dis as text "Caution: `blank' replication cell(s) in ame-family columns returned exact structural zeros (reassignments without dose variation in some cells) and were excluded from those columns' permutation distributions."
  * BS-23 (author ruling 2026-09-07): incompleteness costs resolution, not the convention;
  * the add-one wording now sits on the case that actually uses the add-one form
  if "`wdp_incomplete'" != "" dis as text "Caution: under permexact the following elements were not estimable in every enumerated assignment, so their exact shares are taken over the assignments that completed - a reduced reference set whose size is in e(perm_n): `wdp_incomplete'"
  if "`wdp_addone'" != "" dis as text "Caution: under permexact the observed assignment's own replication did not carry the following elements, so their reference set does not contain the observed statistic and their p-values use the add-one form (1 + count)/(n + 1) rather than an exact share: `wdp_addone'"
end

cap prog drop wooldid_permjnt
prog def wooldid_permjnt, eclass
version 16.0
*** chunk 1 (2026-09-05): the permutation joint pre-trend tests, across-subgroup
*** equality tests, adversarial line tests and pairwise subgroup contrasts. The three
*** quadratic-form families are one-sided upper-tail permutation tests of W = d' S^-1 d
*** on a FIXED metric S, the permutation covariance of the relevant coefficient draws
*** (replications plus the comparator row, listwise complete), so they exist in both
*** the coefficient and the t-ratio modes and consume no per-replication covariance.
*** Contrasts are permuted as coefficient differences in both modes. Names mirror the
*** analytic surface so the display prints them unchanged; an analytic value that the
*** observed run posted (permtstat or permanalyticse) is always overwritten or blanked
*** obsb() is accepted for signature symmetry with wooldid_permwire and is deliberately
*** NOT read: every statistic on this surface (W_obs in permqf, t_obs in permstep, W_obs
*** in permadv, s_obs in permpw) comes from cmpb(), the rep-0 comparator run of the TRUE
*** assignment taken down the same fastnose path as the replications (wooldid_permrun,
*** see its header), so the exact test's "the observed value is a member of the reference
*** set" property holds BITWISE rather than to within a numerical tolerance. Two
*** deliberate consequences, invisible at the sites below: e(adversarial_line_*) is
*** fitted on the comparator's coefficients, not on the ones the event-study table
*** prints, so under permtstat/permanalyticse it differs from them in the last digits;
*** and the pairwise block MIXES the two on purpose - p ranks the comparator difference
*** while the shift interval is anchored on the displayed analytic difference read out of
*** e(pwsgs), which is still analytic here because wooldid_permwire does not touch it
*** (wooldid_bootjnt states the same for its own read of e(pwsgs))
  syntax , acc(string) obsb(string) cmpb(string) [level(real 0) m(string) orow(integer 0) ts(integer 0) jt(string) adv(string) advalt(string) advesrel(string) advtt(string) esf(string) absp(string) side(string) rawp(string) sd(string)]
  if `level' == 0 local level = c(level)
  * p-value convention code (F052/F064, 2026-09-06): 0 equal-tailed (default), 1 absolute
  * add-one (permabsp), 2 raw |s| share c/n without add-one (permrawp: the ritest / Stata
  * permute number), 3 left one-sided, 4 right one-sided (permside())
  local wdpabs = cond("`rawp'" != "", 2, cond("`side'" == "left", 3, cond("`side'" == "right", 4, ("`absp'" != ""))))
  local wdpM = cond("`m'" == "" | "`m'" == ".", ".", "`m'")
  if "`advesrel'" == "" local advesrel -1
  tempname wm
  matrix `wm' = e(wdidmainresults)
  local nmain = rowsof(`wm')
  local sgl `e(wdid_sglist)'
  local ncaut 0
  local wdpprec 0

  * -------- per-column (eq, name, sg, relyr) metadata via the positional walk
  * ALIGNMENT CONTRACT: slot k of this walk IS column k of `acc'. That holds because the
  * accumulator is built with colnames = colfullnames of e(b) (wooldid_permrun) and e(b)
  * stacks wdidpooled then the subgroup blocks in e(wdid_sglist) order - the order walked
  * here. Every Mata call below (permqf, permstep, permadv, permadvline, permpw) indexes
  * `acc' BY POSITION with these k's, so a walk that slips by one is a wrong number, not
  * an error. This program carries no check of its own on purpose: wooldid_permwire runs
  * the identical walk first from the same caller (the shell calls permwire, then permjnt)
  * and hard-aborts with error 459 on any mismatch. Never call permjnt without permwire
  * ahead of it, and never reorder those two calls in the shell
  local kk 0
  local rnames : rowfullnames `wm'
  foreach rn of local rnames {
    local ++kk
    gettoken eqn nmn : rn, parse(":")
    local nmn = subinstr("`nmn'", ":", "", 1)
    local wdpeq`kk' `eqn'
    local wdpnm`kk' `nmn'
    local wdpsg`kk' .
    local wdprl`kk' = `wm'[`kk', 7]
    if "`eqn'" == "pretreatment" | `wdprl`kk'' < 0 local wdpprec 1
  }
  foreach sg of local sgl {
    tempname wsg
    * the `cap confirm' guards the trap wooldid_bootjnt documents at its own walk: a bare
    * `matrix X = e(namethatdoesnotexist)' fabricates a 1x1 missing matrix instead of
    * failing. The `continue' below is NOT a benign skip - the accumulator still carries
    * that subgroup's columns, so every later k would be shifted onto the wrong ones - it
    * is unreachable in practice only because of permwire's error 459 (see above)
    cap confirm matrix e(wdidsgresults_`sg')
    if _rc != 0 continue
    matrix `wsg' = e(wdidsgresults_`sg')
    local rnames : rowfullnames `wsg'
    local rr 0
    foreach rn of local rnames {
      local ++kk
      local ++rr
      gettoken eqn nmn : rn, parse(":")
      local nmn = subinstr("`nmn'", ":", "", 1)
      local wdpeq`kk' `eqn'
      local wdpnm`kk' `nmn'
      local wdpsg`kk' `sg'
      local wdprl`kk' = `wsg'[`rr', 7]
    }
  }
  local K `kk'
  ereturn scalar perm_precaveat = `wdpprec'

  * -------- joint pre-trend tests (jointtests): identity constraint over the family's
  * pre-horizon columns (at least two, as the analytic test requires)
  local wdpjw
  if "`jt'" != "" {
    foreach fam in att att_it ame ame_it {
      local ename = cond(strpos("`fam'", "ame"), "joint_Choriz", "joint_horiz")
      if "`fam'" == "att_it" | "`fam'" == "ame_it" local ename `ename'_i
      foreach sgtag in "" `sgl' {
        local sgw = cond("`sgtag'" == "", ".", "`sgtag'")
        local cols
        forval k = 1/`K' {
          if "`wdpeq`k''" == "ES_`fam'" & "`wdpsg`k''" == "`sgw'" & `wdprl`k'' < 0 local cols `cols' `k'
        }
        local q : word count `cols'
        if `q' < 2 continue
        tempname Cm Wr
        matrix `Cm' = I(`q')
        local jp .
        local jn 0
        local jw .
        mata: wooldid_permqf("`acc'", "`cmpb'", "`cols'", "`Cm'", `wdpM', `orow', "jp", "jn", "jw", "`Wr'")
        local en = cond("`sgtag'" == "", "`ename'", "`ename'_s`sgtag'")
        if `jp' < . {
          * CS-60 (author ruling 2026-09-07): a scalar at ALL posting sites, so a loop over subgroups gets one type
          ereturn scalar `en' = `jp'
          matrix colnames `Wr' = `en'
          tempname Wo
          matrix `Wo' = J(1, 1, `jw')
          matrix colnames `Wo' = `en'
          if "`wdpjw'" == "" {
            tempname jwall jwobs
            matrix `jwall' = `Wr'
            matrix `jwobs' = `Wo'
            local wdpjw 1
          }
          else {
            matrix `jwall' = `jwall', `Wr'
            matrix `jwobs' = `jwobs', `Wo'
          }
        }
        else {
          local ++ncaut
          * no blanking step here, unlike the equality and adversarial blocks below. The
          * worker computes the analytic joint pre-trend tests only when
          * "`suppressjointtests'`bootstrap'`permute'" == "", so under permute() these
          * names were never posted and "absent rather than analytic" is literally true.
          * The allsgssame scalars and the adversarial p-values are NOT gated on permute()
          * and therefore do survive, which is why those two blocks blank explicitly.
          * Dropping the `permute' term from that worker guard would silently put analytic
          * joint p-values back under a permutation header here
          local wdpc`ncaut' the permutation joint pre-trend test `en' is unavailable (fewer than 2 complete replications carry every pre-treatment horizon, or the permutation-covariance metric is degenerate); it is absent rather than analytic
        }
      }
    }
    if "`wdpjw'" != "" {
      ereturn matrix perm_jointW = `jwall'
      ereturn matrix perm_jointWobs = `jwobs'
    }
  }

  * -------- F065 (quality review 2026-09-05, author ruling 2026-09-06): multiplicity across
  * the event-study horizons - the max-|studentized| step-down (Westfall-Young / Romano-Wolf)
  * from the same accumulator, exact under the sharp null: every horizon's draws are
  * studentized on the fixed metric (the SD of its null draws, as the kernel's row 6), the
  * horizons are ordered by |t_obs|, and the j-th adjusted p is the share of replications
  * whose max |t_r| over the horizons not yet rejected reaches |t_obs(j)|, made monotone.
  * Posted per family as e(perm_stepdown[_i][_pre][_s<sg>]) with rows t_obs, p_raw (the
  * per-horizon absolute-value share), p_singlestep (max over all horizons of the family)
  * and p_stepdown, one named column per horizon.
  * PS-45 (author ruling 2026-09-07), two changes. (a) OPT-IN: the procedure used to run
  * whenever at least two post horizons existed, so it shipped output on runs that never
  * asked for it; it now runs only under permstepdown. Nothing else about it moved, so a
  * permstepdown run reproduces the numbers this block printed before. (b) The PRE-treatment
  * horizons are a SEPARATE family, run through the identical procedure and posted and
  * reported apart from the post-treatment one rather than pooled with it - pooling would
  * change the post-treatment adjusted p-values by enlarging their family, and the two ask
  * different questions. The pre-period family inherits the pre-period caveat the footer
  * already carries (e(perm_precaveat)): in staggered designs the reassigned onsets widen
  * the null for pre-period statistics, so its adjusted p-values are conservative in a way
  * the post-treatment family's are not.
  * The engine stays permutation-only (same ruling): the bootstrap path gets simultaneous
  * bands instead, which are a different object - a uniform confidence band rather than a
  * family-wise-error-controlling test.
  if "`sd'" != "" {
    foreach fam in att att_it ame ame_it {
      local sname = cond(strpos("`fam'", "ame"), "perm_Cstepdown", "perm_stepdown")
      if "`fam'" == "att_it" | "`fam'" == "ame_it" local sname `sname'_i
      foreach wdpsw in post pre {
        local sname2 = cond("`wdpsw'" == "pre", "`sname'_pre", "`sname'")
        foreach sgtag in "" `sgl' {
          local sgw = cond("`sgtag'" == "", ".", "`sgtag'")
          local cols
          forval k = 1/`K' {
            * relyr missing is +infinity in stata, so "< 0" already excludes it
            if "`wdpeq`k''" == "ES_`fam'" & "`wdpsg`k''" == "`sgw'" {
              if "`wdpsw'" == "post" & `wdprl`k'' >= 0 & `wdprl`k'' < . local cols `cols' `k'
              if "`wdpsw'" == "pre"  & `wdprl`k'' < 0                  local cols `cols' `k'
            }
          }
          if `: word count `cols'' < 2 continue
          tempname SD
          mata: wooldid_permstep("`acc'", "`cmpb'", "`cols'", `wdpM', `orow', "`SD'")
          * wooldid_permstep signals "no usable step-down" with st_matrix(out, J(0,0,.)).
          * MEASURED (Stata 19): st_matrix() with a 0x0 argument DROPS the Stata matrix
          * rather than creating an empty one, so `confirm' and `rowsof' both give r(111).
          * The `cap confirm' is load-bearing, not defensive. Keep `tempname SD' INSIDE
          * this loop: hoisting it (which the surrounding style invites) would let a
          * failed family inherit the previous family's SD, pass both guards, and be
          * posted under this family's e() name
          cap confirm matrix `SD'
          if _rc != 0 continue
          if rowsof(`SD') < 4 continue
          local en = cond("`sgtag'" == "", "`sname2'", "`sname2'_s`sgtag'")
          if strlen("`en'") <= 32 ereturn matrix `en' = `SD'
        }
      }
    }
  }

  * -------- across-subgroup equality (allsgssame): chain contrast over the subgroup
  * level columns of each estimand. Every name the worker can post is walked so that no
  * analytic value survives under a permutation header
  if `: word count `sgl'' > 1 {
    foreach it in att att_it att_i att_itime ame ame_it ame_i ame_itime {
      foreach pre in "" "pre_" {
        local eqwant = cond("`pre'" == "", "main", "pretreatment")
        local cols
        local ok 1
        foreach sg of local sgl {
          local found
          forval k = 1/`K' {
            if "`wdpeq`k''" == "`eqwant'" & "`wdpnm`k''" == "`it'" & "`wdpsg`k''" == "`sg'" local found `k'
          }
          if "`found'" == "" local ok 0
          else local cols `cols' `found'
        }
        local jp .
        if `ok' {
          local mm : word count `cols'
          tempname Cm Wr
          matrix `Cm' = J(`mm'-1, `mm', 0)
          forval j = 1/`=`mm'-1' {
            matrix `Cm'[`j', `j'] = 1
            matrix `Cm'[`j', `j'+1] = -1
          }
          local jn 0
          local jw .
          mata: wooldid_permqf("`acc'", "`cmpb'", "`cols'", "`Cm'", `wdpM', `orow', "jp", "jn", "jw", "`Wr'")
          if `jp' < . {
            ereturn scalar joint_`pre'`it'_allsgssame = `jp'
            matrix colnames `Wr' = joint_`pre'`it'_allsgssame
            tempname Wo
            matrix `Wo' = J(1, 1, `jw')
            matrix colnames `Wo' = joint_`pre'`it'_allsgssame
            * CAUTION: `wdpsw' carries two unrelated meanings in this program - the
            * foreach loop macro over post/pre in the step-down block above, and from
            * here the "accumulator already started" flag. The test is correct only
            * because Stata clears a foreach macro when its loop exits (measured). The
            * sibling flags wdpjw and wdpaw are initialised explicitly just above their
            * blocks; this one is not, so rewriting that loop as a forvalues, or setting
            * `wdpsw' anywhere above, sends the first post down the else branch into
            * `matrix  = , <tmp>' - r(198) after a full permutation run has been paid for
            if "`wdpsw'" == "" {
              tempname swall swobs
              matrix `swall' = `Wr'
              matrix `swobs' = `Wo'
              local wdpsw 1
            }
            else {
              matrix `swall' = `swall', `Wr'
              matrix `swobs' = `swobs', `Wo'
            }
          }
          else {
            local ++ncaut
            local wdpc`ncaut' the permutation `pre'`it' subgroup-equality test is unavailable (fewer than 2 complete replications carry every subgroup level, or the permutation-covariance metric is degenerate); it is absent rather than analytic
          }
        }
        if `jp' >= . & "`e(joint_`pre'`it'_allsgssame)'" != "" ereturn scalar joint_`pre'`it'_allsgssame = .
      }
    }
    if "`wdpsw'" != "" {
      ereturn matrix perm_sgW = `swall'
      ereturn matrix perm_sgWobs = `swobs'
    }
  }

  * -------- adversarial line tests: wooldid_advW on the permutation-covariance metric
  * over the family's horizon columns; p under the sharp null of no effect at any
  * horizon (esfixedbaseperiod only, the analytic test's own domain); the V-free display
  * line is posted where the fast-mode observed run left none
  local wdpaw
  if "`adv'" != "" {
    local alt = ("`advalt'" != "")
    local esfflag = ("`esf'" != "")
    foreach fk in att att_it ame ame_it {
      foreach sgtag in "" `sgl' {
        local sgw = cond("`sgtag'" == "", ".", "`sgtag'")
        local cols
        local rels
        forval k = 1/`K' {
          if "`wdpeq`k''" == "ES_`fk'" & "`wdpsg`k''" == "`sgw'" & `wdprl`k'' < . local cols `cols' `k'
          if "`wdpeq`k''" == "ES_`fk'" & "`wdpsg`k''" == "`sgw'" & `wdprl`k'' < . local rels `rels' `wdprl`k''
        }
        if `: word count `cols'' < 1 continue
        local ename adversarial_pval_`fk'`=cond("`sgtag'" == "", "", "_s`sgtag'")'
        local lname adversarial_line_`fk'`=cond("`sgtag'" == "", "", "_s`sgtag'")'
        if `esfflag' {
          tempname Wr
          local jp .
          local jn 0
          local jw .
          mata: wooldid_permadv("`acc'", "`cmpb'", "`cols'", "`rels'", `advesrel', ("`advtt'" != ""), `alt', `wdpM', `orow', "jp", "jn", "jw", "`Wr'")
          if `jp' < . {
            ereturn scalar `ename' = `jp'
            matrix colnames `Wr' = `ename'
            tempname Wo
            matrix `Wo' = J(1, 1, `jw')
            matrix colnames `Wo' = `ename'
            if "`wdpaw'" == "" {
              tempname awall awobs
              matrix `awall' = `Wr'
              matrix `awobs' = `Wo'
              local wdpaw 1
            }
            else {
              matrix `awall' = `awall', `Wr'
              matrix `awobs' = `awobs', `Wo'
            }
          }
          else {
            if "`e(`ename')'" != "" ereturn scalar `ename' = .
            local ++ncaut
            local wdpc`ncaut' the permutation adversarial line test `ename' is unavailable (fewer than 2 complete replications carry every horizon, or the permutation-covariance metric leaves no testable residual dimension); it is absent rather than analytic
          }
        }
        else if "`e(`ename')'" != "" ereturn scalar `ename' = .
        cap confirm matrix e(`lname')
        if _rc != 0 {
          tempname ln
          mata: wooldid_permadvline("`cmpb'", "`cols'", "`rels'", `advesrel', ("`advtt'" != ""), `esfflag', `alt', "`ln'")
          if `ln'[1, 1] < . {
            matrix colnames `ln' = slope intercept_rt0
            matrix rownames `ln' = y1
            ereturn matrix `lname' = `ln'
          }
        }
      }
    }
    if "`wdpaw'" != "" {
      ereturn matrix perm_advW = `awall'
      ereturn matrix perm_advWobs = `awobs'
    }
    ereturn scalar perm_adv = 1
  }

  * -------- pairwise subgroup contrasts: coefficient differences in both modes
  cap confirm matrix e(pwsgs)
  if _rc == 0 {
    tempname pw pwa pk o1
    matrix `pw' = e(pwsgs)
    matrix `pwa' = `pw'
    local prn : rowfullnames `pw'
    local npw : word count `prn'
    matrix `pk' = J(18, `npw', .)
    local rr 0
    local pwunmapped 0
    local pwmiss 0
    foreach rn of local prn {
      local ++rr
      gettoken eqn nmn : rn, parse(":")
      local nmn = subinstr("`nmn'", ":", "", 1)
      local ispre = strpos("`eqn'", "pre") > 0
      local item = subinstr("`eqn'", "pre", "", 1)
      local item = subinstr("`item'", "_", "", 1)
      if !`ispre' local item `eqn'
      local eqwant = cond(`ispre', "pretreatment", "main")
      local ok 0
      if regexm("`nmn'", "^([0-9]+)(bn?)?vs([0-9]+)(bn?)?([.].*)?$") {
        local sga = regexs(1)
        local sgb = regexs(3)
        local ka
        local kb
        forval k = 1/`K' {
          if "`wdpeq`k''" == "`eqwant'" & "`wdpnm`k''" == "`item'" & "`wdpsg`k''" == "`sga'" local ka `k'
          if "`wdpeq`k''" == "`eqwant'" & "`wdpnm`k''" == "`item'" & "`wdpsg`k''" == "`sgb'" local kb `k'
        }
        if "`ka'" != "" & "`kb'" != "" {
          mata: wooldid_permpw("`acc'", "`cmpb'", `ka', `kb', st_matrix("`pw'")[`rr', 1], `level', `wdpM', `orow', `wdpabs', "`o1'")
          matrix `pk'[1, `rr'] = `o1'
          local ok 1
        }
      }
      * column 3 (t/z) is blanked unconditionally and renamed t_perm below: contrasts are
      * permuted as COEFFICIENT DIFFERENCES in both modes (see the program header), so
      * there is no permuted t to report even under permtstat - unlike the main table,
      * whose column 3 does carry one (wooldid_permwire). The name is kept so the display
      * prints the analytic header layout unchanged
      matrix `pw'[`rr', 3] = .
      * `pk' is read by RAW ROW INDEX: 13 is the liveness flag, 3 is p, 7 and 8 are the
      * interval endpoints. The layout is defined by wooldid_permtail (rows 1-13 in its
      * header, rows 14-18 where they are appended); the J(18, ...) allocation above and
      * the eighteen-name `matrix rownames `pk'' below must track that kernel
      if `ok' & `pk'[13, `rr'] == 1 {
        matrix `pw'[`rr', 4] = `pk'[3, `rr']
        matrix `pw'[`rr', 5] = `pk'[7, `rr']
        matrix `pw'[`rr', 6] = `pk'[8, `rr']
      }
      else {
        forval cc = 4/6 {
          matrix `pw'[`rr', `cc'] = .
        }
        if !`ok' local ++pwunmapped
        else local ++pwmiss
      }
    }
    local pwc : colnames `pw'
    local pwc : subinstr local pwc "se" "se_analytic", word
    local pwc : subinstr local pwc "t" "t_perm", word
    local pwc : subinstr local pwc "z" "t_perm", word
    local pwc : subinstr local pwc "p" "p_perm", word
    if `level' != 95 {
      local wdplvltxt = cond(mod(`level',1)==0, strofreal(`level'), subinstr(strofreal(`level',"%9.0g"), ".", "p", .))
      local pwc : subinstr local pwc "lb_95ci" "lb_`wdplvltxt'ci", word
      local pwc : subinstr local pwc "ub_95ci" "ub_`wdplvltxt'ci", word
    }
    matrix colnames `pw' = `pwc'
    ereturn hidden matrix perm_analytic_pwsgs = `pwa'
    ereturn matrix pwsgs = `pw'
    matrix rownames `pk' = n s_obs p qlo qhi sd lb ub se_shift offcentre bandlo bandhi live p_lo p_hi se_p cp_lo cp_hi
    matrix colnames `pk' = `prn'
    ereturn matrix perm_pwsgs_kernel = `pk'
    if `pwunmapped' > 0 {
      local ++ncaut
      local wdpc`ncaut' `pwunmapped' pairwise subgroup contrast row(s) could not be mapped to replication columns; their inference is left blank
    }
    if `pwmiss' > 0 {
      local ++ncaut
      local wdpc`ncaut' `pwmiss' pairwise subgroup contrast row(s) have no usable permutation distribution (fewer than 2 completed replications carry both levels); their inference is left blank
    }
  }

  forval cc = 1/`ncaut' {
    dis as text "Caution: `wdpc`cc''."
  }
end

cap prog drop wooldid_permcellsx
prog def wooldid_permcellsx, eclass
version 16.0
  * rewires every e(cfx*) matrix to permutation inference from the cell accumulators
  * (chunk 1, 2026-09-05): p = the comparator cell ranked among the replications' cells,
  * interval = the shift interval about the displayed cell; studentized by the cell
  * standard errors under permtstat. Lattice types whose bins are functions of the
  * assignment (baseTreated*, byTreated*, byTtre*) are not permutation statistics -
  * lattice_k is a different dose range in every replication - so those matrices keep
  * their points and have their inference columns blanked, with a caution
  syntax , cell(string) cell0(string) [ncfx(integer 0) names(string) level(real 0) m(string) orow(integer 0) ts(integer 0) lattype(string) cfxlat(string) absp(string) side(string) rawp(string)]
  if `level' == 0 local level = c(level)
  * p-value convention code (F052/F064, 2026-09-06): 0 equal-tailed (default), 1 absolute
  * add-one (permabsp), 2 raw |s| share c/n without add-one (permrawp: the ritest / Stata
  * permute number), 3 left one-sided, 4 right one-sided (permside())
  local wdpabs = cond("`rawp'" != "", 2, cond("`side'" == "left", 3, cond("`side'" == "right", 4, ("`absp'" != ""))))
  local wdpM = cond("`m'" == "" | "`m'" == ".", ".", "`m'")
  * this whitelist is the exact COMPLEMENT of the latticetype list the validator accepts:
  * the four assignment-invariant types (all, byCohort, byTime, byCohortByTime) are the
  * ones left out. The test is "do this type's bins depend on which units are treated, or
  * on ttre" - a new latticetype that bins the dose on treated units' rows must be added
  * here, or its cells get a permutation p and interval computed across replications whose
  * bins are different objects: a wrong statistic, silently, with no caution. The
  * cfxlat clause is not an exception to that test but an instance of it: a user-supplied
  * cfxlattice() variable replaces latticevar as the plot's over() variable, and a data
  * variable cannot move under reassignment
  * both `mata:' one-liners in the loop below leave a GLOBAL Mata object behind; it is named
  * wooldid_Mtmp rather than `M' so it cannot collide with a user's own (D-3, ruling
  * 2026-09-11) - see the same note in wooldid_bootcellsx
  local drift = (inlist("`lattype'", "baseTreated", "baseTreatedPost", "baseTreatedPre", "baseTreatedRef") | inlist("`lattype'", "byTreated", "byTreatedByTime", "byTtre", "byTtreByTime")) & "`cfxlat'" == ""
  local ndone 0
  forval c = 1/`ncfx' {
    local nm : word `c' of `names'
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    * six columns is the transposed r(table) layout (1 b, 2 se, 3 t, 4 pvalue, 5 ll,
    * 6 ul) that the 1, 2, 4, 5, 6 arguments to wooldid_permcellwire below index. Fewer
    * than six means cfxmatendpoint was 1, which under permute() means cfxplotnose -
    * fastnose is carved out of that truncation precisely so this rewiring has columns to
    * fill. CAUTION: this skip is silent. permute() + cfxplotnose leaves ndone at 0, so
    * the caution at the foot of this program never prints and e(perm_cfx_blank) is never
    * posted, and the cfx cells stay analytic under a permutation header with nothing said
    if colsof(e(`nm')) < 6 continue
    cap confirm matrix `cell'_cfx`c'
    if _rc != 0 continue
    tempname cf cfa
    matrix `cf' = e(`nm')
    matrix `cfa' = `cf'
    if `drift' {
      * what is actually blanked: columns 4-6 (p, ll, ul) always, and column 3 (t) only
      * when the statistic is not the t-ratio. Column 2 keeps the ANALYTIC standard error,
      * as every permutation surface does, and under permtstat column 3 keeps the analytic
      * b/se. The caution printed at the foot of this program says "inference columns are
      * blank" - which, since the D-2 ruling of 2026-09-11, says exactly this: the p-value
      * and interval columns are blank, the se column is analytic, and the t is blank unless
      * permtstat is in force
      mata: wooldid_Mtmp = st_matrix("`cf'"); wooldid_Mtmp[., 4..6] = J(rows(wooldid_Mtmp), 3, .); if (!`ts') wooldid_Mtmp[., 3] = J(rows(wooldid_Mtmp), 1, .); st_matrix("`cf'", wooldid_Mtmp)
    }
    else {
      mata: wooldid_permcellwire("`cf'", "`cell'_cfx`c'", "`cell0'_cfx`c'", "`=cond(`ts', "`cell'_cfx`c'se", "")'", "`=cond(`ts', "`cell0'_cfx`c'se", "")'", 1, 2, 4, 5, 6, `level', `wdpM', `orow', `wdpabs')
      mata: wooldid_Mtmp = st_matrix("`cf'"); wooldid_Mtmp[., 3] = (`ts' ? wooldid_Mtmp[., 1] :/ wooldid_Mtmp[., 2] : J(rows(wooldid_Mtmp), 1, .)); st_matrix("`cf'", wooldid_Mtmp)
    }
    * the worker's original column names are restored, so e(cfx*) column 2 ships under the
    * plain name "se" while holding the ANALYTIC standard error beside a permutation p and
    * a permutation shift interval (wooldid_permcellwire leaves that column as posted).
    * This is the one permutation surface that does not rename it to se_analytic - the
    * pairwise block above and wooldid_permwire both do - because these column names are
    * consumed by wooldid_cfxplots and by marginsplot-shaped downstream code that expects
    * the r(table) names b/se/t/pvalue/ll/ul. The provenance is disclosed in the footer
    * instead of in the name
    matrix colnames `cf' = `: colnames e(`nm')'
    matrix rownames `cf' = `: rowfullnames e(`nm')'
    ereturn matrix `nm' = `cf'
    if strlen("perm_analytic_`nm'") <= 32 ereturn hidden matrix perm_analytic_`nm' = `cfa'
    * the cell draws and the comparator cells, keyed by cell, so every cell p and
    * interval is reconstructable from e() (row 1 of permcmp_* is the comparator)
    * each twin pair is guarded by its own prefix length (coefficient twins at 9
    * characters, se twins at 11): a 22-23 character cfx name keeps permdist_/permcmp_
    * and loses only the se twins, which is what the help says
    if strlen("permdist_`nm'") <= 32 & !`drift' {
      tempname cd cd0
      matrix `cd' = `cell'_cfx`c'
      matrix `cd0' = `cell0'_cfx`c'
      cap ereturn hidden matrix permdist_`nm' = `cd'
      cap ereturn hidden matrix permcmp_`nm' = `cd0'
    }
    if `ts' & strlen("permdistse_`nm'") <= 32 & !`drift' {
      tempname cds cds0
      matrix `cds' = `cell'_cfx`c'se
      matrix `cds0' = `cell0'_cfx`c'se
      cap ereturn hidden matrix permdistse_`nm' = `cds'
      cap ereturn hidden matrix permcmpse_`nm' = `cds0'
    }
    local ++ndone
  }
  if `drift' & `ndone' > 0 {
    ereturn scalar perm_cfx_blank = 1
    dis as text "Caution: the lattice type `lattype' bins the dose on treated units' rows, so its bins are functions of the treatment assignment and differ from replication to replication; the continuous-treatment plot cells therefore have no permutation distribution and carry no permutation inference: their p-value and interval columns are blank, the standard error column keeps the analytic value the estimation posted (as it does on every permutation surface), and the t-ratio is blank unless permtstat is in force. latticetype all, byCohort, byTime, byCohortByTime, or a user-supplied cfxlattice() variable, give assignment-invariant bins that are permuted."
  }
end

cap prog drop wooldid_permsave
prog def wooldid_permsave
version 16.0
  * permsave(): write the replication accumulator to <file>.dta in the working directory
  * (the bootsave() layout: _rep 0 = observed estimates, 1..R = replications, failed ones
  * all-missing; se_* under permtstat) plus one onset_<cluster> column per cluster holding
  * the onset that cluster carried in each replication (row 0: the observed onset)
  syntax , file(string) acc(string) obsb(string) [seacc(string) aseobs(string) asg(string) clusvals(string) obsons(string) seed(string) rng(string) replace]
  local fn `file'
  if lower(substr("`fn'", -4, 4)) == ".dta" local fn = substr("`fn'", 1, strlen("`fn'") - 4)
  if strpos("`fn'", "/") | strpos("`fn'", "\") | "`fn'" == "" {
    dis as error "Error: permsave() takes a plain file name (no directory part); the file is written to the current working directory."
    error 198
  }
  cap confirm file "`fn'.dta"
  if _rc == 0 & "`replace'" == "" {
    dis as error "Error: permsave(`fn') - `fn'.dta already exists in the working directory; add the replace option to overwrite it."
    error 602
  }
  tempname psf
  cap frame drop `psf'
  frame create `psf'
  cap noi frame `psf' {
    mata: wooldid_bootsavemk("`acc'", "`obsb'", "`seacc'", "`aseobs'", "")
    local wdpsasgerr
    local wdpsasgdrop 0
    if "`asg'" != "" mata: wooldid_permsaveasg("`asg'", "`clusvals'", "`obsons'")
    char _dta[wooldid_engine] "permute"
    char _dta[wooldid_seed] "`seed'"
    char _dta[wooldid_rngstate] "`rng'"
    label var _rep "permutation replication (0 = observed estimates)"
    save "`fn'.dta", replace
    local wdpsn = _N
  }
  local wdpsrc = _rc
  cap frame drop `psf'
  if `wdpsrc' != 0 {
    dis as error "Error: permsave(`fn') could not build the replication dataset (rc=`wdpsrc')."
    error `wdpsrc'
  }
  dis as text "permsave(): the permutation draws were written to `fn'.dta (`=`wdpsn'-1' replication rows plus the _rep = 0 observed row`=cond("`asg'" != "", "; the onset_* columns hold the onset assigned to each cluster", "")')."
  * D-10 (author ruling 2026-09-11): the assignment block's two failure modes are announced
  * rather than left to be discovered in the saved file
  if "`wdpsasgerr'" == "1" dis as text "Caution: permsave(`fn') could not write the per-cluster assignment block (the cluster-value and observed-onset lists do not conform to the assignment matrix), so `fn'.dta carries the replication draws but NO onset_* columns."
  if `=cond("`wdpsasgdrop'" == "", "0", "`wdpsasgdrop'")' > 0 dis as text "Caution: permsave(`fn') omitted `wdpsasgdrop' cluster assignment column(s) from `fn'.dta because their onset_<id> names were already taken (a name collision with another column, or two cluster ids that map to the same Stata name). The file is that many clusters short and nothing in it marks the gap."
end




cap prog drop wooldid_bootrun
prog def wooldid_bootrun
version 16.0
*** F-boot engine (plan sections 2-4 + amendments): builds the cluster map frame,
*** generates the Dirichlet draw matrix (mata, rwgen-bitwise arithmetic), runs the
*** rep loop as nested fastnose wooldid calls under reweighted [aw/pw], and
*** accumulates the FULL rep e(b) vectors name-aligned into the caller's matrix.
*** bayes semantics (2026-08-30; default flipped to bootalpha(equal) on the MC evidence
*** in review_process/bootstrap_report_2026-08-30.md section 4; the equal-alpha row
*** multiplier corrected 2026-09-06, quality-review F001): cluster c draws Dirichlet mass
*** draw_c. Under bootalpha(equal) every row of the cluster is multiplied by draw_c (cluster
*** total = W_c * draw_c, Rubin's cluster Bayesian bootstrap, centred on the observation-
*** weighted estimand); under bootalpha(size) the mass is spread EQUALLY over the cluster's
*** rows (row multiplier = draw_c / n_c, n_c the row count and also the Gamma shape), which
*** matches bayesboot BITWISE on unweighted data. The divisor was the cluster's base-weight
*** total W_c until 2026-09-13: identical unweighted, since W_c == n_c on a total() of ones,
*** and a displaced replication distribution under weights (quality-review F001).
*** The rwgen priorpowers(-n_c) oracle for equal is its weights times n_c.
  syntax , [mode(string) reid phys reps(integer 0) obsb(string) accname(string) y(varname) i(varname) t(varname) ttrevar(varname) clus(varlist) dots(integer 1) rseed(string) bootalpha(string) repopts(string asis) wtype(string) wexp(string) ifr(string asis) pois(string) compname(string) seaccname(string) pwseaccname(string) jwaccname(string) jobsname(string) awaccname(string) aobsname(string) advesrel(string) advtt(string) advalt(string) cellacc(string) wmat(string) jkacc(string)]

  if `dots' == -1 local dots 1
  if "`bootalpha'" == "" local bootalpha equal

  * sample concept (A7): user-if sample + nonmissing cluster (permrun byte-parallel)
  tempvar smp
  if `"`ifr'"' != "" {
    qui gen byte `smp' = ((`ifr') == 1)
    qui replace `smp' = 0 if `smp' >= .
  }
  else qui gen byte `smp' = 1
  * observation-level bootstrap: cluster := observation (a unique row id)
  tempvar obsclus
  local clusv `clus'
  if "`clusv'" == "" {
    qui gen double `obsclus' = _n
    local clusv `obsclus'
  }
  qui replace `smp' = 0 if missing(`clusv')

  * base weight, materialized once (reghdfe rejects weight expressions - measured)
  tempvar basew
  if "`wexp'" != "" qui gen double `basew' = (`wexp') if `smp'
  else qui gen double `basew' = 1 if `smp'

  * map frame: one row per cluster with n_c, the cluster's ROW COUNT. n_c is both the
  * Gamma shape under bootalpha(size) and the divisor that spreads the cluster's drawn
  * mass over its rows - the two have to be the same quantity or the replications centre
  * somewhere the point estimate does not (see wooldid_bootdw). The base-weight total W_c
  * that used to be built here as __wdbW was that divisor until 2026-09-13 and is now
  * nothing's input, so it is gone rather than left standing: a live weight-total variable
  * beside a shape-proportional divisor is the trap that produced the defect.
  * frozen generation-time row ids (NFPERM-02: frlink re-sorts the linked frame)
  tempname mbfr
  frame put `clusv' `basew' if `smp', into(`mbfr')
  frame `mbfr' {
    qui bysort `clusv': gen double __wdbn = _N
    * in-sample ROW COUNT (the map frame holds the in-sample rows only, one per
    * observation until the collapse below): the observed row of a bootsave() file under
    * bootalpha(size) is written on the replication scale 1/N (BS-18). The multipliers a
    * replication actually applies are p_c/n_c, which sum to 1 over ROWS - sum_c n_c *
    * (p_c/n_c) = sum_c p_c = 1 - so 1/N is the flat multiplier on that same scale. It was
    * 1/W, sum(basew), while the divisor was W_c; on unweighted data sum(basew) == N, so
    * only weighted bootsave() output moves
    local wdb_nrows = _N
    qui bysort `clusv': keep if _n == 1
    qui sort `clusv'
    qui gen double __wdbdraw = .
    qui gen long __wdbrow = _n
    local wdb_C = _N
  }
  if `wdb_C' < 2 {
    dis as error "Error: bootstrap() requires at least 2 clusters in the estimation sample."
    error 459
  }
  * F063 (2026-09-06): the resample itself, one row per replication and one column per
  * cluster (map-frame order) - Dirichlet multipliers under bayes, multiplicities under
  * standard - kept for bootsave(). Only under cluster(): observation resampling would
  * mean one column per observation
  local wdb_dow = ("`wmat'" != "" & "`clus'" != "")
  if `wdb_dow' {
    matrix `wmat' = J(`reps', `wdb_C', .)
    frame `mbfr': mata: wooldid_bootwnames("`wmat'", "`clusv'")
  }

  * draws. bayes: C x R Dirichlet matrix up front - per-rep (outer) per-cluster-
  * ascending (inner) Gamma draws from the seeded mt64 stream, the exact rwgen bayes
  * consumption order; row multiplier = (g/colsum)/n_c. standard: draws happen
  * INSIDE the rep loop (one bsample per rep - the official bootstrap prefix's
  * consumption sequence, so same-seed runs are comparable; nested wooldid calls
  * leave c(rngstate) untouched, measured)
  if "`rseed'" != "" set seed `rseed'
  * F059 (2026-09-06): the RNG state before the first draw, posted as e(boot_rngstate)
  * and written into bootsave() files, so a run is replayable from e() alone
  local wdb_rng0 = c(rngstate)
  tempname wdbdrawm
  if "`mode'" == "bayes" mata: wooldid_bootdw("`mbfr'", "`bootalpha'", `reps', "`wdbdrawm'")

  * rep command line
  local wtif
  local wtype2 `wtype'
  if "`wtype2'" == "" local wtype2 = cond("`pois'" != "", "pweight", "aweight")
  * UNCLUSTERED STUDENTIZED FIX (2026-09-04). On this one path the replication multiplicity
  * is a nonnegative INTEGER count, so passing it as an fweight rather than an aweight makes
  * each replication's standard error the same functional of the replicate data that the
  * observed standard error is of the original - which is what percentile-t requires and what
  * the aweight form violated. An aweight loses the duplicated rows from N (wrong residual
  * DOF) and squares the weight in the robust meat; an fweight is a literal duplication count
  * and does neither. Measured arithmetically IDENTICAL to physically expanding the rows
  * (se_fw^2/se_phys^2 = 1.0000000 on all 400 pairs), and the coefficient draws are BITWISE
  * unchanged, so no reported estimate or interval family moves.
  * Every gate is necessary: `wtype' empty means the base weight is 1 so the multiplicity is
  * the whole weight and is integral; `seaccname' non-empty is boottstat, the only mode that
  * consumes a replication's own standard error; `clus' empty is the unclustered state, the
  * clustered defect being a different one that no weight type repairs; and mode standard
  * excludes bayes, whose Dirichlet mass is continuous and cannot be an fweight
  * the fifth condition, `pois' == "", is a belt rather than a gate: wooldid_bootctl
  * already refuses boottstat with poisson at parse, and `seaccname' non-empty IS
  * boottstat, so it can never bind. And note why an fweight is legal at all here - the
  * shell refuses fweights from users and accepts one only while ${WOOLDID_REPMODE} is
  * set, then demotes it back to an aweight immediately, so the fweight reaches the
  * REGRESSION and nothing else (margins, gstats and gquantiles keep aweight semantics).
  * The companion case, an unclustered boottstat run WITH a user weight, is handled by the
  * physical-expansion leg below; together the two cover the partition
  if "`wtype'" == "" & "`pois'" == "" & "`seaccname'" != "" & "`clus'" == "" & "`mode'" == "standard" local wtype2 fweight
  local repcall `y' `i' `t' `ttrevar'
  tempvar repw lnk drawv

  * rep-0 comparator: fastnose run of the OBSERVED weights (A3/NFPERM-10 precedent;
  * BC/ties in a later commit compare replicates against THESE values, and the
  * master-name sanity check below guards the name-aligned accumulation)
  local wtif0
  if "`wtype'" != "" local wtif0 [`wtype'=`wexp']
  if `"`ifr'"' != "" local wtif0 `wtif0' if (`ifr')
  * the nested-replication idiom (nobreak / stash _rc / re-raise) recurs six times in this
  * program and is explained in full at its first occurrence, in wooldid_permrun. The two
  * parts that matter most here: the nobreak is what guarantees WOOLDID_REPMODE is cleared
  * even on a user Break, and the `if _rc == 1' re-raise is what stops a Ctrl-Break being
  * tallied as just another failed replication and the loop running on to completion,
  * producing inference the user believed they had cancelled
  nobreak {
    global WOOLDID_REPMODE 1
  cap qui wooldid `repcall' `wtif0' , `repopts'
    local wdrepdone = _rc
    global WOOLDID_REPMODE
  }
  cap error `wdrepdone'
  if _rc == 1 error 1
  if _rc != 0 {
    dis as error "Error: the rep-0 comparison run of the observed data failed (rc `=_rc'); bootstrap inference is impossible. The observed estimation results are unaffected."
    error 459
  }
  if "`compname'" != "" matrix `compname' = e(b)
  * master-name sanity: the rep e(b) column set must equal the observed run's
  local wdbnames0 : colfullnames `obsb'
  local wdbnames1 : colfullnames e(b)
  if "`wdbnames0'" != "`wdbnames1'" {
    dis as error "Error: the rep-0 coefficient vector does not match the observed run's column set; bootstrap accumulation would be misaligned. This is unexpected - please report it."
    error 459
  }

  * accumulation matrix: R x K over the OBSERVED column set (A17)
  local wdb_K = colsof(`obsb')
  matrix `accname' = J(`reps', `wdb_K', .)
  matrix colnames `accname' = `wdbnames0'
  * under boottstat a parallel R x K of per-replication STANDARD ERRORS, filled by the
  * same name-aligned grab; the reps are not run under fastnose in that mode, so each
  * one has an e(V) whose diagonal supplies them
  if "`seaccname'" != "" {
    matrix `seaccname' = J(`reps', `wdb_K', .)
    matrix colnames `seaccname' = `wdbnames0'
  }
  * under boottstat, ALSO each replication's own pairwise-contrast standard errors: one
  * column per row of the observed e(pwsgs), keyed by that row's full name (eq:JvsK), and
  * filled by a name-anchored grab of the replication's e(pwsgs) column 2. That column is
  * the analytic margins pairwise se, which already folds in the covariance between the two
  * subgroup levels - the posted e(V) is diagonal, so this is the only place a per-rep
  * contrast se exists. Absent when the observed run posted no pairwise table
  * under boottstat, the recentred ES joint statistic W_r per family key, computed INSIDE
  * each replication (wooldid_bootjwr) from that replication's own e(V_es_*), against the
  * observed pre-horizon coefficients stacked here from the LIVE e() - which at this point
  * belongs to the rep-0 COMPARATOR, not to the observed run (see the F018 note below).
  * Keys mirror wooldid_bootwire's family walk: att att_it ame ame_it and their _s<sg>
  * twins. A key with no observed e(b_es_*) simply contributes no column
  * WHOSE e() IS LIVE HERE: from the rep-0 comparator fit above down to the frlink that
  * opens the replication loop, the live e() is the COMPARATOR's, not the observed run's.
  * Every key, name list and observed statistic harvested in this accumulator block is
  * therefore the comparator's. That is safe because the comparator is the observed data
  * at the observed weights and its column set is checked against the observed run's
  * above - F018 measured the two key sets to coincide on every configuration tested - but
  * anyone adding an accumulator here must not assume that e() entries which only the
  * observed run posts are available
  if "`jwaccname'" != "" {
    local wdb_jkeys att att_it ame ame_it
    foreach sg in `e(wdid_sglist)' {
      local wdb_jkeys `wdb_jkeys' att_s`sg' att_it_s`sg' ame_s`sg' ame_it_s`sg'
    }
    local wdb_jobsn 0
    mata: wooldid_bootjobs("`wdb_jkeys'", "`jobsname'")
    if `wdb_jobsn' > 0 {
      local wdb_jfk
      local wdb_jcs : colfullnames `jobsname'
      foreach cf of local wdb_jcs {
        gettoken fk rest : cf, parse(":")
        if `: list posof "`fk'" in wdb_jfk' == 0 local wdb_jfk `wdb_jfk' `fk'
      }
      matrix `jwaccname' = J(`reps', `: word count `wdb_jfk'', .)
      matrix colnames `jwaccname' = `wdb_jfk'
    }
    else {
      local jwaccname
      local jobsname
    }
  }
  * under boottstat + adversarial: the recentred adversarial statistic per (family, alt),
  * computed inside each replication by wooldid_bootadvr against the observed coefficients
  * over ALL horizons stacked here
  * cell-level accumulators (gteffects rows keyed by i|onset|relyr; each cfx matrix keyed by
  * its row names). `cellacc' is a prefix; the matrices are `cellacc'_gt, `cellacc'_gtse and
  * `cellacc'_cfx<k>, `cellacc'_cfx<k>se with a name map in the locals returned to the shell.
  * NB (F018 research record, 2026-09-06): these keys are minted from the REP-0 comparator
  * run's e(), which is live here (nothing restores the observed e() between the comparator
  * assignment above and this point), not from the observed run's. Measured on every
  * configuration tested the two key sets coincide (no observed cell lacked a rep-0 key and
  * the cfx matrix name lists match); a keyless observed cell would be left with missing
  * inference by wooldid_bootcellwire rather than mis-keyed
  local wdb_ncfx 0
  if "`cellacc'" != "" {
    cap confirm matrix e(gteffects)
    if _rc == 0 {
      mata: wooldid_bootcellkeys("e(gteffects)", "`cellacc'_gt", `reps', 1)
      if "`seaccname'" != "" mata: wooldid_bootcellkeys("e(gteffects)", "`cellacc'_gtse", `reps', 1)
      if "`wdb_cellkeydup'" != "" dis as text "Caution: the cohort-by-period cells could not be keyed uniquely for the bootstrap; e(gteffects) inference is left analytic and blanked."
    }
    foreach cm in `e(wdid_cfxmats)' {
      cap confirm matrix e(`cm')
      if _rc != 0 continue
      local ++wdb_ncfx
      mata: wooldid_bootcellkeys("e(`cm')", "`cellacc'_cfx`wdb_ncfx'", `reps', 0)
      if "`seaccname'" != "" mata: wooldid_bootcellkeys("e(`cm')", "`cellacc'_cfx`wdb_ncfx'se", `reps', 0)
      local wdb_cfxname`wdb_ncfx' `cm'
    }
  }
  c_local wdb_ncfx `wdb_ncfx'
  forval kk = 1/`wdb_ncfx' {
    c_local wdb_cfxname`kk' `wdb_cfxname`kk''
  }
  if "`awaccname'" != "" {
    local wdb_akeys att att_it ame ame_it
    foreach sg in `e(wdid_sglist)' {
      local wdb_akeys `wdb_akeys' att_s`sg' att_it_s`sg' ame_s`sg' ame_it_s`sg'
    }
    local wdb_jobsn 0
    mata: wooldid_bootjobs("`wdb_akeys'", "`aobsname'", 1)
    if `wdb_jobsn' > 0 {
      local wdb_afk
      local wdb_acs : colfullnames `aobsname'
      foreach cf of local wdb_acs {
        gettoken fk rest : cf, parse(":")
        if `: list posof "`fk'" in wdb_afk' == 0 local wdb_afk `wdb_afk' `fk'
      }
      local wdb_acols
      foreach fk of local wdb_afk {
        local wdb_acols `wdb_acols' `fk'_a0
        if "`advalt'" != "" local wdb_acols `wdb_acols' `fk'_a1
      }
      matrix `awaccname' = J(`reps', `: word count `wdb_acols'', .)
      matrix colnames `awaccname' = `wdb_acols'
    }
    else {
      local awaccname
      local aobsname
    }
  }
  if "`pwseaccname'" != "" {
    cap confirm matrix e(pwsgs)
    if _rc == 0 {
      local wdb_pwn : rowfullnames e(pwsgs)
      local wdb_npw : word count `wdb_pwn'
      matrix `pwseaccname' = J(`reps', `wdb_npw', .)
      matrix colnames `pwseaccname' = `wdb_pwn'
    }
    else local pwseaccname
  }

  qui frlink m:1 `clusv', frame(`mbfr') generate(`lnk')
  local ndone 0
  local nfail 0
  local nblank 0
  if `dots' > 0 dis as text "Bootstrap replications (`reps'): " _continue
  tempname wdbcf wdbrf
  forval r = 1/`reps' {
    * F066 (2026-09-06): time replications 1-5 and print a projection for the run (uses
    * Stata timer 100 for the duration of those five replications)
    if `r' == 1 & `reps' > 5 {
      timer clear 100
      timer on 100
    }
    if `r' == 6 {
      timer off 100
      qui timer list 100
      local wdb_proj = r(t100)/5*(`reps' + cond("`jkacc'" != "", `wdb_C', 0))
      timer clear 100
      * no "." in the text: the dots line is counted by tests and by eye
      local wdb_ptxt = cond(`wdb_proj' < 90, strofreal(round(`wdb_proj')) + " s", strofreal(round(`wdb_proj'/60)) + " min")
      if `dots' > 0 dis as text " (~`wdb_ptxt' projected for `reps' replications`=cond("`jkacc'" != "", " + `wdb_C' jackknife fits", "")')" _continue
    }
    if "`mode'" == "bayes" {
      frame `mbfr': qui replace __wdbdraw = el("`wdbdrawm'", __wdbrow, `r')
      if `wdb_dow' frame `mbfr': mata: wooldid_bootwrow("`wmat'", `r', "__wdbdraw")
    }
    else if "`phys'" != "" {
      * physical-expansion leg WITHOUT re-identification (unclustered + user weight +
      * boottstat): duplicate the drawn observations, keep the original i, carry the user's
      * weight expression onto every copy via `wtif0'. The replication's own VCE is then the
      * ordinary robust/homoskedastic one on the expanded rows, the same functional as the
      * observed - which a multiplicity aweight could not deliver (it loses duplicated rows
      * from N and squares the weight in the robust meat). Same-seed draws equal the count
      * leg's (one bsample per rep, same form)
      cap frame drop `wdbrf'
      qui frame
      frame copy `r(currentframe)' `wdbrf'
      frame `wdbrf' {
        qui keep if `smp'
        bsample, cluster(`clusv') idcluster(__wdbns)
        nobreak {
          global WOOLDID_REPMODE 1
        if "${WOOLDID_BOOTDEBUG}" != "" cap noi wooldid `y' `i' `t' `ttrevar' `wtif0' , `repopts'
        else cap qui wooldid `y' `i' `t' `ttrevar' `wtif0' , `repopts'
          local wdrepdone = _rc
          if "${WOOLDID_BOOTDEBUG}" != "" dis as error "WOOLDID_BOOTDEBUG: replication `r' returned rc=`wdrepdone'"
          global WOOLDID_REPMODE
        }
        cap error `wdrepdone'
        local wdb_rrc = _rc
      }
      if `wdb_rrc' == 1 {
        cap frame drop `wdbrf'
        if `dots' > 0 dis as result " "
        error 1
      }
      local badrep = `wdb_rrc' != 0
      if `badrep' == 0 {
        local wdb_gblank 0
        mata: wooldid_bootgrab("`accname'", `r', "`obsb'", "`seaccname'", "`pwseaccname'")
        if "`jwaccname'" != "" mata: wooldid_bootjwr("`jwaccname'", `r', "`jobsname'")
        if "`awaccname'" != "" mata: wooldid_bootadvr("`awaccname'", `r', "`aobsname'", `advesrel', ("`advtt'" != ""))
        if "`cellacc'" != "" mata: wooldid_bootcells("`cellacc'", `r', `wdb_ncfx', ("`seaccname'" != ""))
        local nblank = `nblank' + `wdb_gblank'
        local ++ndone
      }
      else local ++nfail
      cap frame drop `wdbrf'
      * D-4 (author ruling 2026-09-11): dots() is honoured here as it is on the count path,
      * the re-identification leg and the jackknife loop. This leg used to print one dot per
      * replication whatever dots() said, and never the 50-mark running count - so the option
      * meant one thing on three legs and another on this one, and a test counting marks here
      * had to expect `reps' of them rather than reps/dots
      if `dots' > 0 {
        * the documented cadence, identical to the count path, the re-identification leg and
        * wooldid_permrun: one mark every dots() replications, a failed replication ALWAYS
        * prints its x, and the count after every 50 marks
        if `badrep' dis as error "x" _continue
        else if mod(`r', `dots') == 0 dis as text "." _continue
        if mod(`r', 50*`dots') == 0 dis as text " " `r' _newline "  " _continue
      }
      * this leg is complete: without this, the count path's common tail below would run
      * the same replication a second time on frequency weights and overwrite row r
      continue
    }
    else if "`reid'" == "" {
      * standard count path: cluster multiplicities from bsample on a COPY of the
      * map frame (map-frame draws == full-data draws under the cluster()+
      * idcluster() form, measured; plain bsample diverges - A17)
      cap frame drop `wdbcf'
      frame copy `mbfr' `wdbcf'
      frame `wdbcf' {
        bsample, cluster(`clusv') idcluster(__wdbns)
        if `wdb_dow' mata: wooldid_bootwrowx("`wmat'", `r', "`clusv'", "__wdbns")
        qui gegen double __wdbcnt = count(__wdbns), by(`clusv')
        qui bysort `clusv': keep if _n == 1
      }
      frame `mbfr' {
        cap drop __wdblnk
        cap drop __wdbcnt
        qui frlink m:1 `clusv', frame(`wdbcf') generate(__wdblnk)
        qui frget __wdbcnt = __wdbcnt, from(__wdblnk)
        qui replace __wdbdraw = cond(missing(__wdbcnt), 0, __wdbcnt)
        qui drop __wdbcnt __wdblnk
        if `wdb_dow' mata: wooldid_bootwrow("`wmat'", `r', "__wdbdraw")
      }
      frame drop `wdbcf'
    }
    else {
      * physical re-id leg (cell-weighted estimands; plan section 7): resample the
      * FULL estimation universe and re-identify duplicated clusters' units so a
      * twice-drawn cluster contributes its cells twice - the textbook cluster
      * bootstrap of the equal-cell-weighted estimands. cluster nests i (bootctl
      * gate A1), so group(drawid, i) assigns each copy fresh unit codes. same-seed
      * draws equal the count path's (one bsample per rep, same form)
      cap frame drop `wdbrf'
      qui frame
      frame copy `r(currentframe)' `wdbrf'
      frame `wdbrf' {
        qui keep if `smp'
        bsample, cluster(`clusv') idcluster(__wdbns)
        if `wdb_dow' mata: wooldid_bootwrowx("`wmat'", `r', "`clusv'", "__wdbns")
        qui gegen double __wdbnewi = group(__wdbns `i')
        * under boottstat the replication's own VCE must cluster on the RESAMPLED cluster
        * identity __wdbns (which bsample just generated), not on the user's original
        * cluster variable that repopts carries: on this expanded frame a twice-drawn
        * cluster is two physical copies, and clustering on the original id would pool
        * them back into one score - the exact m^2-for-m inflation this leg exists to
        * remove. Without boottstat the replication runs under fastnose and has no VCE,
        * so the substitution is inert there and is skipped to keep those runs bitwise
        local wdb_ropts `repopts'
        * the mechanism is a literal STRING substitution into an option string built far
        * upstream, so it can silently do nothing. It hits only because that string is
        * rebuilt canonically as `cluster(`cluster')' - no spaces, no abbreviation, never
        * the user's own typing. Change that emission (a space, an abbreviation, a
        * vce(cluster ...) form, a reordered cluster list) and subinstr no-ops, the
        * replication clusters on the ORIGINAL ids again, and the per-replication standard
        * error silently returns to the E[m^2]/E[m] = 2 - 1/G (about sqrt(2)) inflation
        * this entire leg exists to remove. Nothing checks it and there is no r()
        if "`seaccname'" != "" & "`clus'" != "" local wdb_ropts : subinstr local wdb_ropts "cluster(`clus')" "cluster(__wdbns)"
        nobreak {
          global WOOLDID_REPMODE 1
        cap qui wooldid `y' __wdbnewi `t' `ttrevar' `wtif0' , `wdb_ropts'
          local wdrepdone = _rc
          global WOOLDID_REPMODE
        }
        cap error `wdrepdone'
        local wdb_rrc = _rc
      }
      if `wdb_rrc' == 1 {
        cap frame drop `wdbrf'
        if `dots' > 0 dis as result " "
        error 1
      }
      local badrep = `wdb_rrc' != 0
      if `badrep' == 0 {
        local wdb_gblank 0
        mata: wooldid_bootgrab("`accname'", `r', "`obsb'", "`seaccname'", "`pwseaccname'")
        if "`jwaccname'" != "" mata: wooldid_bootjwr("`jwaccname'", `r', "`jobsname'")
        if "`awaccname'" != "" mata: wooldid_bootadvr("`awaccname'", `r', "`aobsname'", `advesrel', ("`advtt'" != ""))
        if "`cellacc'" != "" mata: wooldid_bootcells("`cellacc'", `r', `wdb_ncfx', ("`seaccname'" != ""), "`wdbrf'", "`i'")
        local nblank = `nblank' + `wdb_gblank'
        local ++ndone
      }
      else local ++nfail
      cap frame drop `wdbrf'
      if `dots' > 0 {
        if `badrep' dis as error "x" _continue
        else if mod(`r', `dots') == 0 dis as text "." _continue
        if mod(`r', 50*`dots') == 0 dis as text " " `r' _newline "  " _continue
      }
      continue
    }
    cap drop `drawv'
    qui frget `drawv' = __wdbdraw, from(`lnk')
    cap drop `repw'
    qui gen double `repw' = `basew' * `drawv'
    nobreak {
      global WOOLDID_REPMODE 1
    cap qui wooldid `repcall' [`wtype2'=`repw'] , `repopts'
      local wdrepdone = _rc
      global WOOLDID_REPMODE
    }
    cap error `wdrepdone'
    if _rc == 1 {
      if `dots' > 0 dis as result " "
      error 1
    }
    local badrep = _rc != 0
    if `badrep' == 0 {
      * name-aligned grab of the full rep e(b) into row r; ame-family structural
      * zeros (rep 0.0 where observed b != 0) blanked per A12(b); wdb_gblank returns
      * the count blanked this rep
      local wdb_gblank 0
      mata: wooldid_bootgrab("`accname'", `r', "`obsb'", "`seaccname'", "`pwseaccname'")
        if "`jwaccname'" != "" mata: wooldid_bootjwr("`jwaccname'", `r', "`jobsname'")
        if "`awaccname'" != "" mata: wooldid_bootadvr("`awaccname'", `r', "`aobsname'", `advesrel', ("`advtt'" != ""))
        if "`cellacc'" != "" mata: wooldid_bootcells("`cellacc'", `r', `wdb_ncfx', ("`seaccname'" != ""))
      local nblank = `nblank' + `wdb_gblank'
    }
    if `badrep' == 0 local ++ndone
    else local ++nfail
    if `dots' > 0 {
      * one mark every dots() replications - the documented cadence (a failed
      * replication always prints its x) - with the count after every 50 marks
      if `badrep' dis as error "x" _continue
      else if mod(`r', `dots') == 0 dis as text "." _continue
      if mod(`r', 50*`dots') == 0 dis as text " " `r' _newline "  " _continue
    }
  }
  if `dots' > 0 dis as text " done"
  * F067 (quality review 2026-09-05, author ruling 2026-09-06): the BCa acceleration
  * needs a delete-one-cluster jackknife - G extra fast fits, each dropping one cluster
  * by a zero weight through the same map-frame link the count leg uses (a zero weight
  * excludes the observation under every weight type). Rows of jkacc = clusters in
  * map-frame order; a fit that fails leaves its row missing (the acceleration then uses
  * the clusters that fitted, with a caution downstream). Only under bootdisplay(bca)
  if "`jkacc'" != "" {
    matrix `jkacc' = J(`wdb_C', `wdb_K', .)
    matrix colnames `jkacc' = `wdbnames0'
    if `dots' > 0 dis as text "Jackknife for the BCa acceleration (`wdb_C' clusters): " _continue
    forval g = 1/`wdb_C' {
      frame `mbfr': qui replace __wdbdraw = (__wdbrow != `g')
      cap drop `drawv'
      qui frget `drawv' = __wdbdraw, from(`lnk')
      cap drop `repw'
      qui gen double `repw' = `basew' * `drawv'
      nobreak {
        global WOOLDID_REPMODE 1
      cap qui wooldid `repcall' [`wtype2'=`repw'] , `repopts'
        local wdrepdone = _rc
        global WOOLDID_REPMODE
      }
      cap error `wdrepdone'
      if _rc == 1 {
        if `dots' > 0 dis as result " "
        error 1
      }
      if _rc == 0 {
        local wdb_gblank 0
        mata: wooldid_bootgrab("`jkacc'", `g', "`obsb'")
      }
      if `dots' > 0 {
        if `wdrepdone' != 0 dis as error "x" _continue
        else if mod(`g', `dots') == 0 dis as text "." _continue
      }
    }
    if `dots' > 0 dis as text " done"
  }
  cap drop `lnk'
  cap drop `drawv'
  cap drop `repw'
  frame drop `mbfr'

  if `ndone' == 0 {
    dis as error "Error: none of the `reps' bootstrap replications estimated successfully; bootstrap inference is impossible. The observed estimation results are unaffected."
    error 459
  }
  if `nfail' > 0 dis as text "Caution: `nfail' of `reps' bootstrap replications failed to estimate and were excluded (official-bootstrap convention: complete replications only)."

  c_local wdb_r_done `ndone'
  c_local wdb_r_fail `nfail'
  c_local wdb_r_blank `nblank'
  * the resampling-unit count: clusters, or observations under the observation-level
  * branch. bootdf(clusters) turns this into the t degrees of freedom
  c_local wdb_r_nclus `wdb_C'
  c_local wdb_r_nrows `wdb_nrows'
  c_local wdb_r_rng0 `wdb_rng0'
end

cap prog drop wooldid_bootwire
prog def wooldid_bootwire, eclass
version 16.0
*** F-boot rewiring (plan A6/A11/A12/A17): replaces the reported inference surfaces
*** with the bootstrap covariance. Point estimates NEVER change. Surfaces rewired
*** here: e(wdidmainresults) + e(wdidsgresults_*) se/z/p/ci columns, e(V) (listwise
*** complete-case bootstrap covariance; analytic V preserved hidden as
*** e(V_analytic)), e(V_es_*) (name-aligned listwise subblocks; analytic kept hidden
*** as e(V_es_analytic_*)), the VCE header text, and the e(boot_*) metadata.
*** pwsgs/prepwsgs inference columns and the allsgssame scalars are BLANKED in this
*** build (their bootstrap Wald versions land with the joint-test commit).
  syntax , [mode(string) reps(integer 0) acc(string) obsb(string) done(string) fail(string) blank(string) level(real 0) bootalpha(string) clus(varlist) comp(string) TIES rseed(string) bootdf(string) bdisp(string) nclus(integer 0) seacc(string) aseout(string) jkacc(string) sb(string)]
  * OUTPUT CHANNEL to sibling programs, not a local: the caller declares the tempname and
  * passes it in as aseout(); this program fills it near the end with the OBSERVED analytic
  * standard errors, and it survives this program's exit because a tempname matrix is
  * dropped when the program that DECLARED it ends, not when the program that filled it
  * ends. wooldid_bootjnt and the simultaneous band read it from there. Outside boottstat
  * nothing is ever written into the name, so both consumers must treat a nonexistent
  * matrix as "no studentized inputs" and fall back - the same silent contract governs the
  * else-branches that blank jwacc/jobs and awacc/aobs locally while the caller still
  * passes those tempnames on
  local wdbaseobs_out `aseout'
  if "`wdbaseobs_out'" == "" tempname wdbaseobs_out
  if `level' == 0 local level = c(level)
  if "`bootalpha'" == "" local bootalpha equal
  * reference distribution for the reported statistics (bootdf()): missing = normal
  local wdbdf = .
  if "`bootdf'" != "" & "`bootdf'" != "normal" {
    if "`bootdf'" == "clusters" local wdbdf = `nclus' - 1
    else local wdbdf = `bootdf'
    if `wdbdf' < 1 {
      dis as text "Caution: bootdf(`bootdf') resolves to `wdbdf' degrees of freedom; falling back to normal critical values."
      local wdbdf = .
    }
  }
  local wdblab = cond(`wdbdf' < ., "t", "z")

  * -------- column-position walk: e(b) col k <-> results-matrix row (positional:
  * e(b) is built by stacking wdidpooled + the SG blocks in sglist order, so the
  * correspondence is structural). excluded columns (A12a) = rows whose reported
  * ESTIMATE is missing (amevoid and friends re-enter e(b) as structural zeros, so
  * e(b) itself is never missing - the results matrices carry the truth)
  tempname wm excl
  matrix `wm' = e(wdidmainresults)
  local nmain = rowsof(`wm')
  local K = colsof(`obsb')
  matrix `excl' = J(1, `K', 0)
  * boottstat needs the OBSERVED run's own standard errors - the same formula the
  * replications use - to studentize with, so they are captured here, in the same walk,
  * before the rewiring below replaces column 2 with the bootstrap ones
  tempname aseobs
  matrix `aseobs' = J(1, `K', .)
  local sgl `e(wdid_sglist)'
  * THE e(b) COLUMN CURSOR. e(b) is built by stacking wdidpooled then the subgroup blocks
  * in e(wdid_sglist) order, so column k of e(b) is row k of the concatenated results
  * matrices. `kk' walks that concatenation and is deliberately NOT reset at the subgroup
  * boundary. Four separate main-table walks in this program each restart it and run to
  * `nmain', and the subgroup loop then continues from whichever ran last - only THIS
  * first walk is validated (the `kk' != `K'' check below). The later walks, and the
  * `kk0`sg''/`kkP'/`kkS' re-walks inside the subgroup loop, are unguarded: a new walk
  * that resets `kk' and does not end at exactly `nmain', or a stray `local kk 0' between
  * the last main-table walk and the subgroup loop, shifts every subgroup standard error
  * one column with no error and no caution - the coefficients stay right and the
  * inference belongs to the neighbouring estimate
  local kk 0
  forval rr = 1/`nmain' {
    local ++kk
    if `wm'[`rr', 1] >= . matrix `excl'[1, `kk'] = 1
    matrix `aseobs'[1, `kk'] = `wm'[`rr', 2]
  }
  foreach sg of local sgl {
    tempname wsg`sg'
    cap matrix `wsg`sg'' = e(wdidsgresults_`sg')
    if _rc == 0 {
      forval rr = 1/`=rowsof(`wsg`sg'')' {
        local ++kk
        if `wsg`sg''[`rr', 1] >= . matrix `excl'[1, `kk'] = 1
        matrix `aseobs'[1, `kk'] = `wsg`sg''[`rr', 2]
      }
    }
  }
  if `kk' != `K' {
    * structural mismatch between the stacked results matrices and e(b) - should be
    * impossible (suppresssgprimaryeffects is refused with bootstrap); fail loudly
    * rather than rewire misaligned inference
    dis as error "Error: bootstrap rewiring found `kk' results-matrix rows against `K' coefficient columns; the inference surfaces cannot be aligned. This is unexpected - please report it."
    error 459
  }

  * -------- moments (mata): per-column se over available reps (R-1, dev from rep
  * mean); LISTWISE complete-case covariance over non-excluded columns (A11: PSD by
  * construction; pairwise-complete is measurably indefinite); sentinel -1e100
  * diagonal + zeroed off-diagonals for excluded / <2-available columns; low-
  * availability disclosure list
  tempname bse bmean bV bcnt
  local wdb_nlist 0
  local wdb_lowav
  * CONSEQUENCE of the two different denominators, which a user will hit: se[k] is computed
  * over the reps available for column k, while V uses only the rows complete over ALL
  * non-excluded columns. So on any run with missing replication cells
  * sqrt(diag(e(V_boot))) does NOT equal the printed se column, and a hand-run test or
  * lincom on e(V) will not reproduce the table's Wald numbers exactly. e(boot_N_V) is the
  * listwise count behind V; e(boot_navail) the per-column counts behind the se.
  * SIDE CHANNEL: this helper also returns through st_local, setting `wdb_nlw' (the
  * listwise row count posted as e(boot_N_V)), `wdb_exclnames' and `wdb_lowav'. Several
  * other Mata helpers in this program do the same (wdb_vesok, wdb_gblank, wdb_cellkeydup,
  * wdb_jobsn), which is why locals appear below with no assignment anywhere in the ado. A
  * helper that returns early without setting its local leaves the ado reading an empty
  * macro: `ereturn scalar boot_N_V = ' is r(198), and `if "`wdb_vesok'" == "1"' silently
  * takes the failure branch
  mata: wooldid_bootfin("`acc'", "`excl'", "`bse'", "`bmean'", "`bV'", "`bcnt'")

  * -------- the critical value for the reported level. HOISTED to here (it used to be
  * computed just below, at the head of the rewiring loop) because the percentile
  * package now needs it before the loop runs: the displayed standard error under that
  * package is the interval inverted, (ub-lb)/(2*zc). It depends only on `wdbdf', fixed
  * above, and `level', so moving it up is a pure relocation
  local zc = cond(`wdbdf' < ., invt(`wdbdf', 1 - (100 - `level')/200), invnormal(1 - (100 - `level')/200))

  * -------- CI flavors: normal / percentile / bias-corrected at level(), the classic
  * percentile p dual to the percentile interval, and the BC median-bias z0 (compared
  * against the REP-0 fastnose comparator per amendment A3; `ties' adjusts z0 for
  * replicate-comparator ties, the built-in bootstrap's ties semantics).
  * This call USED TO LIVE IN wooldid_bootjnt, which runs after this program. It was
  * hoisted here because the percentile package is now a display default and the
  * rewiring loop below has to read `cip'/`cipp' to build it; bootjnt cannot supply
  * them in time, and it can see neither `excl' nor `seacc', both of which the
  * blanking and boottstat rules need
  tempname cin cip cib cipp cbp cz0
  mata: wooldid_bootci("`acc'", "`obsb'", "`comp'", `level', ("`ties'" != ""), "`cin'", "`cip'", "`cib'", "`cipp'", "`cbp'", "`cz0'", `wdbdf')
  * F067: the BCa family (only under bootdisplay(bca), when the jackknife ran): the
  * acceleration a_k from the delete-one-cluster jackknife, z0 as the BC family's, the
  * endpoints at the BCa-adjusted levels by the same order-statistic rule, and the p
  * that inverts them at the null
  local wdbbca = ("`jkacc'" != "")
  if `wdbbca' {
    cap confirm matrix `jkacc'
    if _rc != 0 local wdbbca 0
  }
  tempname ciba cbap cacc
  if `wdbbca' mata: wooldid_bootbca("`acc'", "`obsb'", "`comp'", "`jkacc'", `level', ("`ties'" != ""), "`ciba'", "`cbap'", "`cacc'")
  else {
    matrix `ciba' = J(2, colsof(`obsb'), .)
    matrix `cbap' = J(1, colsof(`obsb'), .)
    matrix `cacc' = J(1, colsof(`obsb'), .)
  }

  * -------- boottstat: the studentized (percentile-t) tail p and critical value per
  * column, from t_r = (b_r - b_obs)/se_r against t_obs = b_obs/se_obs
  local wdbts = ("`seacc'" != "")
  tempname bts
  if `wdbts' {
    cap confirm matrix `seacc'
    if _rc != 0 {
      dis as text "Caution: boottstat was requested but no per-replication standard errors were harvested; the reported inference falls back to the coefficient-draw families."
      local wdbts 0
    }
    else mata: wooldid_boottp("`acc'", "`seacc'", "`obsb'", "`aseobs'", "`excl'", `level', ("`ties'" != ""), "`bts'")
  }

  * -------- rewire the results matrices: cols 2-6 = se t/z p lb ub at level()
  * (`zc' is computed above, before the CI flavors, which consume it)
  local kk 0
  forval rr = 1/`nmain' {
    local ++kk
    if `excl'[1, `kk'] == 0 & `bse'[1, `kk'] < . & `bse'[1, `kk'] > 0 {
      matrix `wm'[`rr', 2] = `bse'[1, `kk']
      matrix `wm'[`rr', 3] = `wm'[`rr', 1]/`bse'[1, `kk']
      matrix `wm'[`rr', 4] = cond(`wdbdf' < ., 2*ttail(`wdbdf', abs(`wm'[`rr', 3])), 2*normal(-abs(`wm'[`rr', 3])))
      matrix `wm'[`rr', 5] = `wm'[`rr', 1] - `zc'*`bse'[1, `kk']
      matrix `wm'[`rr', 6] = `wm'[`rr', 1] + `zc'*`bse'[1, `kk']
    }
    else {
      forval cc = 2/6 {
        matrix `wm'[`rr', `cc'] = .
      }
    }
  }
  * column labels: the test statistic is z-based at level() (bootstrap convention)
  * a non-integer level still needs a legal matrix column name: 92.5 -> lb_92p5ci
  local wdblvltxt = cond(mod(`level',1)==0, strofreal(`level'), subinstr(strofreal(`level',"%9.0g"), ".", "p", .))
  local wmc : colnames `wm'
  local wmc : subinstr local wmc "t" "`wdblab'", word
  local wmc : subinstr local wmc "z" "`wdblab'", word
  if `level' != 95 {
    local wmc : subinstr local wmc "lb_95ci" "lb_`wdblvltxt'ci", word
    local wmc : subinstr local wmc "ub_95ci" "ub_`wdblvltxt'ci", word
  }
  matrix colnames `wm' = `wmc'

  * -------- both inference families, always. Author ruling 2026-09-03: every flavour is
  * computed and stored in e() on every bootstrap run, and bootdisplay() chooses ONLY
  * which one the printed table carries. The sd family is what the loop above just
  * wrote; the percentile family is built here from the percentile interval and its
  * dual p. The percentile standard error is the interval INVERTED, (ub-lb)/(2*zc),
  * which deliberately does not satisfy b +/- zc*se when the interval is asymmetric -
  * the author ruled the interval the primary object and accepted that mismatch. The
  * statistic column keeps its z/t label and holds b/se_inv; its tail probability is
  * NOT the reported p, which is a count, and the help says so
  local bdisp = cond("`bdisp'" == "", "percentile", "`bdisp'")
  tempname wmsd wmpct wmbc wmbca
  local wdbbcamiss 0
  matrix `wmsd' = `wm'
  matrix `wmpct' = `wm'
  matrix `wmbc' = `wm'
  matrix `wmbca' = `wm'
  * The standard error column is the bootstrap STANDARD DEVIATION in every non-studentized
  * display, percentile and bias-corrected included. Author ruling 2026-09-04, overruling
  * his own earlier preference for a standard error implied by inverting the reported
  * interval: that construction is not what the literature or Stata's own bootstrap prefix
  * does, and wooldid follows the convention. So b +/- c*se still does not reproduce a
  * percentile or BC interval - but that is the ordinary, expected reading of a bootstrap
  * table, not a wooldid peculiarity.
  * The endpoints are copied by MATRIX reference, never through a local: `local x =
  * matrix[i,j]' round-trips the double through a decimal string, which loses the last
  * bits and would leave the displayed interval merely close to e(boot_ci_percentile)
  * rather than identical to it
  local kk 0
  forval rr = 1/`nmain' {
    local ++kk
    if `excl'[1, `kk'] == 0 & `bse'[1, `kk'] < . & `bse'[1, `kk'] > 0 & `cip'[1, `kk'] < . & `cip'[2, `kk'] < . {
      matrix `wmpct'[`rr', 2] = `bse'[1, `kk']
      matrix `wmpct'[`rr', 3] = `wmpct'[`rr', 1]/`bse'[1, `kk']
      matrix `wmpct'[`rr', 4] = `cipp'[1, `kk']
      matrix `wmpct'[`rr', 5] = `cip'[1, `kk']
      matrix `wmpct'[`rr', 6] = `cip'[2, `kk']
    }
    else {
      forval cc = 2/6 {
        matrix `wmpct'[`rr', `cc'] = .
      }
    }
    if `excl'[1, `kk'] == 0 & `bse'[1, `kk'] < . & `bse'[1, `kk'] > 0 & `cib'[1, `kk'] < . & `cib'[2, `kk'] < . {
      matrix `wmbc'[`rr', 2] = `bse'[1, `kk']
      matrix `wmbc'[`rr', 3] = `wmbc'[`rr', 1]/`bse'[1, `kk']
      matrix `wmbc'[`rr', 4] = `cbp'[1, `kk']
      matrix `wmbc'[`rr', 5] = `cib'[1, `kk']
      matrix `wmbc'[`rr', 6] = `cib'[2, `kk']
    }
    else {
      forval cc = 2/6 {
        matrix `wmbc'[`rr', `cc'] = .
      }
    }
    if `excl'[1, `kk'] == 0 & `bse'[1, `kk'] < . & `bse'[1, `kk'] > 0 {
      matrix `wmbca'[`rr', 2] = `bse'[1, `kk']
      matrix `wmbca'[`rr', 3] = `wmbca'[`rr', 1]/`bse'[1, `kk']
      * BS-17: an unavailable BCa family (no acceleration, or an adjusted endpoint level
      * outside (0,1)) blanks the p-value and interval only; the bootstrap SD stands
      if `ciba'[1, `kk'] < . & `ciba'[2, `kk'] < . {
        matrix `wmbca'[`rr', 4] = `cbap'[1, `kk']
        matrix `wmbca'[`rr', 5] = `ciba'[1, `kk']
        matrix `wmbca'[`rr', 6] = `ciba'[2, `kk']
      }
      else {
        forval cc = 4/6 {
          matrix `wmbca'[`rr', `cc'] = .
        }
        local ++wdbbcamiss
      }
    }
    else {
      forval cc = 2/6 {
        matrix `wmbca'[`rr', `cc'] = .
      }
    }
  }
  * NB: `ereturn matrix' MOVES its source (the tempname is dropped), so the selected
  * family must be copied into `wm' BEFORE the two posts, not after
  * F053 (quality review 2026-09-05, author ruling 2026-09-06): under the rank-based
  * displays the statistic column is b over the bootstrap standard deviation - a Wald
  * object whose tail probability is NOT the printed rank-based p - so it is labelled
  * t_wald / z_wald there (and in the percentile and bc twins); the wald twin keeps t / z
  local wmcw : colnames `wm'
  local wmcw : subinstr local wmcw "`wdblab'" "`wdblab'_wald", word
  matrix colnames `wmpct' = `wmcw'
  matrix colnames `wmbc' = `wmcw'
  matrix colnames `wmbca' = `wmcw'
  if "`bdisp'" == "percentile" matrix `wm' = `wmpct'
  else if "`bdisp'" == "bc"       matrix `wm' = `wmbc'
  else if "`bdisp'" == "bca"      matrix `wm' = `wmbca'
  ereturn hidden matrix boot_wald_mainresults = `wmsd'
  * CS-54 (author ruling 2026-09-08): every other display family's table name derives from
  * e(boot_display) - bc, bca, wald - but the DEFAULT posts "percentile" against
  * boot_pctl_mainresults, so the derivation fails on exactly the case a user hits first.
  * Post the derivable name as an alias: pure addition, and the short name still works.
  * NB the `, copy' - see the note above: `ereturn matrix' MOVES its source, so the first
  * of the two posts has to copy or the second finds no matrix (r(111) on every bootstrap)
  ereturn hidden matrix boot_pctl_mainresults = `wmpct', copy
  ereturn hidden matrix boot_percentile_mainresults = `wmpct'
  ereturn hidden matrix boot_bc_mainresults = `wmbc'
  if `wdbbca' ereturn hidden matrix boot_bca_mainresults = `wmbca'

  * boottstat: the Wald table is preserved as a hidden twin and columns 3-6 are replaced
  * by the studentized statistic, its percentile-t tail probability and the interval that
  * inverts it. column 2 keeps the bootstrap standard deviation, which is the quantity a
  * bootstrap run is asked for; note that the reported statistic is therefore b/se_obs
  * rather than b/se_reported - the studentization has to use the same standard error the
  * replications use, and the help says so
  if `wdbts' {
    * The studentized (percentile-t) family. boottstat chooses WHAT is bootstrapped (the
    * t-statistics); how that distribution is summarised is fixed at its equal-tailed
    * percentile ranks - bootdisplay(wald) is refused at parse (no standard deviation of
    * a studentized statistic) and so is bootdisplay(bc) (F003, author ruling 2026-09-06:
    * the BC shift on reflected studentized quantiles moved the interval against the
    * median bias; the former BC-t twin is withdrawn). Column 2 is the
    * ANALYTIC standard error, labelled se_analytic (ruling 2026-09-04): it is the
    * studentization's own denominator, and a studentized distribution has no standard
    * deviation in coefficient units. Column 3 is t_obs = b/se_analytic. The interval is
    * EQUAL-TAILED from the signed t-quantiles, so it is not symmetric about b.
    * A column whose studentized family is unavailable has columns 3-6 BLANKED rather
    * than left holding the coefficient-family numbers under a studentized header (that
    * unmarked mixture was a recorded defect), and is counted for the caution below
    tempname wmtp
    matrix `wmtp' = `wm'
    local wdbtmiss 0
    local wdbtlow 0
    local kk 0
    forval rr = 1/`nmain' {
      local ++kk
      * `bts' ROW MAP (rownames are only attached far below): 1 n, 2 t_obs, 3 p_pctl,
      * 4 qlo_pctl, 5 qhi_pctl, 6 live. Row 6 is the liveness flag and `aseobs' must be
      * non-missing independently, which is why both are tested here. The interval lines
      * below look transposed - the LOWER bound uses row 5 (qhi) and the UPPER row 4
      * (qlo) - and are not: that is the percentile-t reflection, [b - q_hi*se_obs,
      * b - q_lo*se_obs], which is equal-tailed rather than symmetric about b (author
      * rulings 2026-09-03/04, recorded in wooldid_boottp)
      local tlive = (`excl'[1, `kk'] == 0 & `bts'[6, `kk'] == 1 & `aseobs'[1, `kk'] < .)
      if `tlive' {
        matrix `wmtp'[`rr', 2] = `aseobs'[1, `kk']
        matrix `wmtp'[`rr', 3] = `bts'[2, `kk']
        matrix `wmtp'[`rr', 4] = `bts'[3, `kk']
        matrix `wmtp'[`rr', 5] = `wm'[`rr', 1] - `bts'[5, `kk']*`aseobs'[1, `kk']
        matrix `wmtp'[`rr', 6] = `wm'[`rr', 1] - `bts'[4, `kk']*`aseobs'[1, `kk']
        if `bts'[1, `kk'] < 0.5*`done' local ++wdbtlow
      }
      else {
        forval cc = 2/6 {
          matrix `wmtp'[`rr', `cc'] = .
        }
        if `excl'[1, `kk'] == 0 local ++wdbtmiss
      }
    }
    local wmc2 : colnames `wm'
    local wmc2 : subinstr local wmc2 "`wdblab'" "t_pt", word
    local wmc2 : subinstr local wmc2 "`wdblab'_wald" "t_pt", word
    local wmc2 : subinstr local wmc2 "se" "se_analytic", word
    matrix colnames `wmtp' = `wmc2'
    matrix `wm' = `wmtp'
    ereturn hidden matrix boot_tpctl_mainresults = `wmtp'
  }
  ereturn matrix wdidmainresults = `wm'
  foreach sg of local sgl {
    cap confirm matrix `wsg`sg''
    if _rc != 0 continue
    local kk0`sg' = `kk'
    forval rr = 1/`=rowsof(`wsg`sg'')' {
      local ++kk
      if `excl'[1, `kk'] == 0 & `bse'[1, `kk'] < . & `bse'[1, `kk'] > 0 {
        matrix `wsg`sg''[`rr', 2] = `bse'[1, `kk']
        matrix `wsg`sg''[`rr', 3] = `wsg`sg''[`rr', 1]/`bse'[1, `kk']
        matrix `wsg`sg''[`rr', 4] = cond(`wdbdf' < ., 2*ttail(`wdbdf', abs(`wsg`sg''[`rr', 3])), 2*normal(-abs(`wsg`sg''[`rr', 3])))
        matrix `wsg`sg''[`rr', 5] = `wsg`sg''[`rr', 1] - `zc'*`bse'[1, `kk']
        matrix `wsg`sg''[`rr', 6] = `wsg`sg''[`rr', 1] + `zc'*`bse'[1, `kk']
      }
      else {
        forval cc = 2/6 {
          matrix `wsg`sg''[`rr', `cc'] = .
        }
      }
    }
    local wsc : colnames `wsg`sg''
    local wsc : subinstr local wsc "t" "`wdblab'", word
    local wsc : subinstr local wsc "z" "`wdblab'", word
    if `level' != 95 {
      local wsc : subinstr local wsc "lb_95ci" "lb_`wdblvltxt'ci", word
      local wsc : subinstr local wsc "ub_95ci" "ub_`wdblvltxt'ci", word
    }
    matrix colnames `wsg`sg'' = `wsc'

    * both families for the subgroup coefficient rows too, same rule as the main table
    * above (and the same `ereturn matrix' move-not-copy caveat: select into the live
    * matrix BEFORE posting the twins)
    tempname wssd`sg' wspct`sg' wsbc`sg' wsbca`sg'
    matrix `wssd`sg'' = `wsg`sg''
    matrix `wspct`sg'' = `wsg`sg''
    matrix `wsbc`sg'' = `wsg`sg''
    matrix `wsbca`sg'' = `wsg`sg''
    local kkP = `kk0`sg''
    forval rr = 1/`=rowsof(`wsg`sg'')' {
      local ++kkP
      if `excl'[1, `kkP'] == 0 & `bse'[1, `kkP'] < . & `bse'[1, `kkP'] > 0 & `cip'[1, `kkP'] < . & `cip'[2, `kkP'] < . {
        matrix `wspct`sg''[`rr', 2] = `bse'[1, `kkP']
        matrix `wspct`sg''[`rr', 3] = `wspct`sg''[`rr', 1]/`bse'[1, `kkP']
        matrix `wspct`sg''[`rr', 4] = `cipp'[1, `kkP']
        matrix `wspct`sg''[`rr', 5] = `cip'[1, `kkP']
        matrix `wspct`sg''[`rr', 6] = `cip'[2, `kkP']
      }
      else {
        forval cc = 2/6 {
          matrix `wspct`sg''[`rr', `cc'] = .
        }
      }
      if `excl'[1, `kkP'] == 0 & `bse'[1, `kkP'] < . & `bse'[1, `kkP'] > 0 & `cib'[1, `kkP'] < . & `cib'[2, `kkP'] < . {
        matrix `wsbc`sg''[`rr', 2] = `bse'[1, `kkP']
        matrix `wsbc`sg''[`rr', 3] = `wsbc`sg''[`rr', 1]/`bse'[1, `kkP']
        matrix `wsbc`sg''[`rr', 4] = `cbp'[1, `kkP']
        matrix `wsbc`sg''[`rr', 5] = `cib'[1, `kkP']
        matrix `wsbc`sg''[`rr', 6] = `cib'[2, `kkP']
      }
      else {
        forval cc = 2/6 {
          matrix `wsbc`sg''[`rr', `cc'] = .
        }
      }
      if `excl'[1, `kkP'] == 0 & `bse'[1, `kkP'] < . & `bse'[1, `kkP'] > 0 {
        matrix `wsbca`sg''[`rr', 2] = `bse'[1, `kkP']
        matrix `wsbca`sg''[`rr', 3] = `wsbca`sg''[`rr', 1]/`bse'[1, `kkP']
        * BS-17: an unavailable BCa family (no acceleration, or an adjusted endpoint level
        * outside (0,1)) blanks the p-value and interval only; the bootstrap SD stands
        if `ciba'[1, `kkP'] < . & `ciba'[2, `kkP'] < . {
          matrix `wsbca`sg''[`rr', 4] = `cbap'[1, `kkP']
          matrix `wsbca`sg''[`rr', 5] = `ciba'[1, `kkP']
          matrix `wsbca`sg''[`rr', 6] = `ciba'[2, `kkP']
        }
        else {
          forval cc = 4/6 {
            matrix `wsbca`sg''[`rr', `cc'] = .
          }
          local ++wdbbcamiss
        }
      }
      else {
        forval cc = 2/6 {
          matrix `wsbca`sg''[`rr', `cc'] = .
        }
      }
    }
    local wscw : colnames `wsg`sg''
    local wscw : subinstr local wscw "`wdblab'" "`wdblab'_wald", word
    matrix colnames `wspct`sg'' = `wscw'
    matrix colnames `wsbc`sg'' = `wscw'
    matrix colnames `wsbca`sg'' = `wscw'
    if "`bdisp'" == "percentile"    matrix `wsg`sg'' = `wspct`sg''
    else if "`bdisp'" == "bc"       matrix `wsg`sg'' = `wsbc`sg''
    else if "`bdisp'" == "bca"      matrix `wsg`sg'' = `wsbca`sg''
    * the hidden twins are skipped when a name would exceed 32 characters (a 13+ digit
    * subgroup value); the visible table is unaffected
    if strlen("boot_pctl_sgresults_`sg'") <= 32 {
      ereturn hidden matrix boot_wald_sgresults_`sg' = `wssd`sg''
      ereturn hidden matrix boot_pctl_sgresults_`sg' = `wspct`sg''
      ereturn hidden matrix boot_bc_sgresults_`sg' = `wsbc`sg''
      if `wdbbca' ereturn hidden matrix boot_bca_sgresults_`sg' = `wsbca`sg''
    }

    * boottstat: same treatment as the main table - Wald kept hidden, columns 3-6 replaced
    if `wdbts' {
      tempname wstp`sg'
      matrix `wstp`sg'' = `wsg`sg''
      local kkS = `kk0`sg''
      forval rr = 1/`=rowsof(`wsg`sg'')' {
        local ++kkS
        local tlive = (`excl'[1, `kkS'] == 0 & `bts'[6, `kkS'] == 1 & `aseobs'[1, `kkS'] < .)
        if `tlive' {
          matrix `wstp`sg''[`rr', 2] = `aseobs'[1, `kkS']
          matrix `wstp`sg''[`rr', 3] = `bts'[2, `kkS']
          matrix `wstp`sg''[`rr', 4] = `bts'[3, `kkS']
          matrix `wstp`sg''[`rr', 5] = `wsg`sg''[`rr', 1] - `bts'[5, `kkS']*`aseobs'[1, `kkS']
          matrix `wstp`sg''[`rr', 6] = `wsg`sg''[`rr', 1] - `bts'[4, `kkS']*`aseobs'[1, `kkS']
          if `bts'[1, `kkS'] < 0.5*`done' local ++wdbtlow
        }
        else {
          forval cc = 2/6 {
            matrix `wstp`sg''[`rr', `cc'] = .
          }
          if `excl'[1, `kkS'] == 0 local ++wdbtmiss
        }
      }
      local wsc2 : colnames `wsg`sg''
      local wsc2 : subinstr local wsc2 "`wdblab'" "t_pt", word
      local wsc2 : subinstr local wsc2 "`wdblab'_wald" "t_pt", word
      local wsc2 : subinstr local wsc2 "se" "se_analytic", word
      matrix colnames `wstp`sg'' = `wsc2'
      matrix `wsg`sg'' = `wstp`sg''
      if strlen("boot_tpctl_sgresults_`sg'") <= 32 ereturn hidden matrix boot_tpctl_sgresults_`sg' = `wstp`sg''
      matrix colnames `wsg`sg'' = `wsc2'
    }
    ereturn matrix wdidsgresults_`sg' = `wsg`sg''
  }

  * -------- e(V): preserve the analytic V hidden, repost the bootstrap covariance
  tempname vold vnew
  * F2 guard: reading an absent e(V) via cap matrix yields a fabricated 1x1 missing
  * scalar - confirm the matrix exists before preserving it
  cap confirm matrix e(V)
  local wdb_hadv = _rc == 0
  if `wdb_hadv' {
    matrix `vold' = e(V)
    ereturn hidden matrix V_analytic = `vold'
  }
  matrix `vnew' = `bV'
  local obn : colfullnames `obsb'
  matrix rownames `vnew' = `obn'
  matrix colnames `vnew' = `obn'
  * repost, not post, ON PURPOSE: `ereturn repost V=' does NOT re-run `ereturn post''s
  * positive-definiteness screen - MEASURED, a 2x2 with V[1,1] = -1e100 reposts rc 0 and
  * comes back verbatim - and the -1e100 sentinel rows wooldid_bootfin manufactures for
  * excluded or under-available columns would otherwise trip it. That is the PS-62 ruling:
  * an unidentified variance prints as missing rather than as a fabricated 0. The cost is
  * that e(V) can leave this program INDEFINITE, so any downstream consumer that inverts
  * it or takes the square root of its diagonal (test, lincom, honestdid, the
  * pretrendpower path) must screen for a non-positive-finite diagonal first
  cap ereturn repost V = `vnew'
  local wdb_vposted = _rc == 0
  if !`wdb_vposted' dis as text "Caution: the bootstrap covariance could not be posted to e(V) (rc `=_rc'); it remains available as e(V_boot). Reported tables carry the bootstrap standard errors either way."
  * NOT redundant: `ereturn repost V = `vnew'' MOVES its source exactly as `ereturn matrix
  * name = mat' does - measured, `confirm matrix' on the source gives rc 0 before and
  * r(111) after - so `vnew' no longer exists here and has to be rebuilt from `bV' (the
  * retained master copy, kept for precisely this reason) with its stripes re-applied
  * before e(V_boot) can be posted. Delete these three lines and every bootstrap run dies
  * r(111) on the ereturn below
  cap matrix `vnew' = `bV'
  matrix rownames `vnew' = `obn'
  matrix colnames `vnew' = `obn'
  ereturn matrix V_boot = `vnew'

  * -------- e(V_es_*): name-aligned listwise subblocks replace the analytic F2
  * surface (honestdid entry point); analytic kept hidden. b_es is untouched (point
  * estimates). subblock columns are matched by relyr-anchored e(b) names.
  local wdbf2 att att_it ame ame_it
  foreach sg of local sgl {
    local wdbf2 `wdbf2' att_s`sg' att_it_s`sg' ame_s`sg' ame_it_s`sg'
  }
  foreach fk of local wdbf2 {
    cap confirm matrix e(V_es_`fk')
    if _rc != 0 continue
    tempname evold evnew
    matrix `evold' = e(V_es_`fk')
    mata: wooldid_bootves("`acc'", "`excl'", "`fk'", "`evold'", "`evnew'")
    if "`wdb_vesok'" == "1" {
      ereturn hidden matrix V_es_analytic_`fk' = `evold'
      ereturn matrix V_es_`fk' = `evnew'
    }
    else {
      * the swap failed, so e(V_es_`fk') KEEPS ITS ANALYTIC VALUE. That used to be
      * inert; since pretrendpower() became compatible with bootstrap (2026-09-04) a
      * consumer reads this matrix, so a silent fallback would feed Roth's formula an
      * analytic covariance under a bootstrap header. Record it. Note the marker used
      * downstream is the ABSENCE of e(V_es_analytic_`fk'), which is created only on the
      * success branch above - a signal that already exists rather than a new one
      local wdb_vesbad `wdb_vesbad' `fk'
    }
  }
  if "`wdb_vesbad'" != "" {
    ereturn local boot_ves_failed `wdb_vesbad'
    dis as text "Caution: the bootstrap event-study covariance could not be built for: `wdb_vesbad'. Those e(V_es_*) matrices remain ANALYTIC. Any diagnostic that consumes them - pretrendpower() in particular - is skipped for the affected families rather than computed on an analytic covariance under a bootstrap header."
  }


  * -------- header + metadata
  local wdbmodetext = cond("`mode'" == "bayes", "Bayesian cluster bootstrap", "cluster bootstrap")
  local wdbcl `clus'
  if "`wdbcl'" == "" {
    local wdbmodetext = cond("`mode'" == "bayes", "Bayesian bootstrap", "bootstrap")
    local wdbcl observation-level
  }
  else local wdbcl cluster(`wdbcl')
  * -------- the CI flavors and the percentile p (computed above, before the rewiring
  * loop that consumes them). Posted here, after the repost, alongside the other
  * bootstrap surfaces. `ereturn repost' preserves existing e() entries, so the order
  * relative to 6681 is not load-bearing; keeping them together is for legibility
  local obn : colfullnames `obsb'
  foreach cm in cin cip cib {
    matrix colnames ``cm'' = `obn'
    matrix rownames ``cm'' = ll ul
  }
  matrix colnames `cipp' = `obn'
  matrix colnames `cbp' = `obn'
  matrix colnames `cz0' = `obn'
  ereturn matrix boot_ci_normal = `cin'
  ereturn matrix boot_ci_percentile = `cip'
  ereturn matrix boot_ci_bc = `cib'
  ereturn matrix boot_p_percentile = `cipp'
  * F064 (2026-09-06): the bootstrap analogue of the permutation tail record - the
  * one-sided percentile shares about 0, the Monte-Carlo se of the percentile p and its
  * level() Clopper-Pearson interval, hidden, by name
  tempname wdbpt
  mata: wooldid_bootptails("`acc'", "`obsb'", `level', "`wdbpt'")
  cap matrix colnames `wdbpt' = `: colnames e(boot_p_percentile)'
  matrix rownames `wdbpt' = p_lo p_hi se_p cp_lo cp_hi
  ereturn hidden matrix boot_ptails = `wdbpt'
  ereturn matrix boot_p_bc = `cbp'
  ereturn matrix boot_z0 = `cz0'
  * 2026-09-13 author ruling: disclose the withheld bias-corrected family. A column whose
  * completed draws all sit on one side of the rep-0 comparator has z0 = invnormal(0) or
  * invnormal(1), which is undefined, so wooldid_bootci leaves that column's z0, BC interval
  * and BC p missing - and the BCa twin with them, since it starts from the same z0 - while
  * the percentile and Wald families stand. Read the POSTED matrices, not `cz0'/`cip':
  * `ereturn matrix' moved those. The discriminator is a live percentile endpoint; a column
  * with fewer than 2 available reps or no dispersion carries no bootstrap inference at all
  * and is not a BC withholding. Named columns rather than a bare flag, on the
  * boot_ves_failed precedent: the caution scrolls away, the macro survives a replay.
  * 2026-09-13: the remedy clause is QUALIFIED rather than flat. It used to end "a larger
  * reps() makes a one-sided draw set less likely", which is true only when the cause is
  * sampling noise about the estimate. When the replication distribution is DISPLACED from
  * the estimate the draws are one-sided at every reps() - the weighted bootalpha(size)
  * divisor defect measured at 3.5-9.4 bootstrap deviations was exactly that - so the flat
  * advice sent a user to spend hours of replications on an interval that would still
  * exclude their own estimate. The same sentence in the sthlp bootdisplay() entry is
  * qualified with it; keep the two in step
  tempname wdbz0m wdbpcm
  matrix `wdbz0m' = e(boot_z0)
  matrix `wdbpcm' = e(boot_ci_percentile)
  local wdb_bcn : colnames `wdbz0m'
  local wdb_bcbad
  forval wdbk = 1/`=colsof(`wdbz0m')' {
    if missing(`wdbz0m'[1,`wdbk']) & !missing(`wdbpcm'[1,`wdbk']) {
      local wdb_bcw : word `wdbk' of `wdb_bcn'
      local wdb_bcbad `wdb_bcbad' `wdb_bcw'
    }
  }
  if "`wdb_bcbad'" != "" {
    ereturn local boot_bc_absent `wdb_bcbad'
    dis as text "Caution: the bias-corrected family is withheld for: `wdb_bcbad'. Every completed replication of those columns fell on one side of the rep-0 comparator, so the median-bias correction z0 is invnormal(0) or invnormal(1) and is undefined; e(boot_z0), e(boot_ci_bc) and e(boot_p_bc) are missing there, and the BCa twin with them where it was requested. The percentile and Wald families are unaffected; under bootdisplay(bc) the affected rows print a blank p-value and interval rather than one borrowed from another family. Whether a larger reps() helps depends on the cause: it makes a one-sided draw set less likely when the draws are simply sparse around the estimate, and does nothing when the replication distribution sits away from the estimate altogether, which is a property of the specification rather than of the replication count."
  }
  if `wdbbca' {
    matrix colnames `ciba' = `obn'
    matrix rownames `ciba' = ll ul
    matrix colnames `cbap' = `obn'
    matrix colnames `cacc' = `obn'
    ereturn matrix boot_ci_bca = `ciba'
    ereturn matrix boot_p_bca = `cbap'
    ereturn matrix boot_accel = `cacc'
    mata: st_local("wdbjkn", strofreal(sum(rowmissing(st_matrix("`jkacc'")) :== 0)))
    ereturn scalar boot_jk_completed = `wdbjkn'
    if `wdbjkn' < rowsof(`jkacc') dis as text "Caution: `=rowsof(`jkacc') - `wdbjkn'' of `=rowsof(`jkacc')' delete-one-cluster jackknife fits failed; the BCa acceleration uses the clusters that fitted."
    if `wdbbcamiss' > 0 & "`bdisp'" == "bca" dis as text "Caution: `wdbbcamiss' reported coefficient(s) had no usable BCa family (no acceleration could be computed, or an acceleration-adjusted endpoint level fell outside (0,1)); their p-value and interval are left blank rather than filled from another family. The standard-error column is unaffected."
  }
  if "`ties'" != "" ereturn local boot_ties ties

  * PS-45 (author ruling 2026-09-07): the bootstrap simultaneous band, opt-in. Built from
  * the same accumulator the joint test uses, after the tables are rewired, so the band and
  * the printed pointwise intervals come from one set of replications

  ereturn hidden local wdid_vcespec "`wdbmodetext', `wdbcl', reps(`reps')"
  * CS-4: the bootstrap engine posts the full replication covariance, unless the repost
  * failed (the caution above), in which case whatever e(V) held before still stands
  * PS-62 follow-through (2026-09-10): on that fallback the matrix left standing may be the
  * SPLICED one the worker built, so read the label back rather than demoting it to
  * analytic-diagonal - the wooldid_permwire twin does exactly the same
  local wdb_vs0 "`e(wdid_vcvstate)'"
  ereturn hidden local wdid_vcvstate = cond(`wdb_vposted', "bootstrap", cond(`wdb_hadv', cond("`wdb_vs0'" == "analytic-blocks", "analytic-blocks", "analytic-diagonal"), "absent"))
  * PS-1: Stata's own bootstrap prefix labels the variance this way
  if `wdb_vposted' {
    ereturn local vce bootstrap
    ereturn local vcetype Bootstrap
  }
  ereturn local boot_method `mode'
  ereturn scalar boot_reps = `reps'
  ereturn scalar boot_completed = `done'
  ereturn scalar boot_failed = `fail'
  ereturn scalar boot_level = `level'
  ereturn local boot_stat = cond(`wdbts', "t_pt", cond("`bdisp'" == "wald", "`wdblab'", "`wdblab'_wald"))
  ereturn local boot_display `bdisp'
  if `wdbts' {
    ereturn scalar boot_tstat = 1
    * the whole studentized kernel goes out with row and column names, so every reported
    * studentized number is reconstructable from e() by name (not by position)
    matrix rownames `bts' = n t_obs p_pctl qlo_pctl qhi_pctl live
    matrix colnames `bts' = `obn'
    ereturn matrix boot_tstud = `bts'
    matrix colnames `aseobs' = `obn'
    matrix `wdbaseobs_out' = `aseobs'
    ereturn hidden matrix boot_ase_obs = `aseobs'
    if `wdbtmiss' > 0 dis as text "Caution: `wdbtmiss' reported coefficient(s) had no usable studentized bootstrap distribution (fewer than 2 replications with a positive standard error, or no observed standard error); their statistic, p-value and interval are left blank rather than filled from the coefficient-draw families."
    if `wdbtlow' > 0 dis as text "Caution: `wdbtlow' reported coefficient(s) used fewer than half of the completed replications in their studentized distribution, because the remaining replications produced no usable standard error for that coefficient. Replications are dropped by a rule that correlates with the resample being degenerate, so the studentized tail is trimmed in the anti-conservative direction there; inspect e(boot_tstud) row n."
    dis as text "Caution: boottstat reports studentized (percentile-t) inference: the statistic t_pt is b/se_analytic, and the p-value and interval come from the ranks of the replications' own t-statistics (b_r - b)/se_r, so the standard-error column is the ANALYTIC standard error (se_analytic), the studentization's denominator, not a bootstrap standard deviation; the interval is equal-tailed from the signed t-quantiles and is not symmetric about the estimate. Under boottstat the pairwise subgroup contrasts are studentized the same way from each replication's own pairwise standard error, and the joint pre-trend tests, the across-subgroup equality tests and the adversarial line tests are bootstrapped as test statistics. bootdisplay(bc) and bootdisplay(wald) are refused under boottstat: the studentized draws have no bias-corrected or Wald (normal or t approximation) family. The studentized kernel is in e(boot_tstud) and the coefficient and standard-error draws in e(bootdist) and e(bootdist_se)."
  }
  * PS-45: the bootstrap simultaneous band. Placed HERE, after the studentized block above,
  * because under boottstat the band is a percentile-t and needs the observed analytic
  * standard errors that block builds - at its earlier position `wdbaseobs_out' did not yet
  * exist and the band silently fell back to the coefficient-draw SD.
  if "`sb'" != "" wooldid_esband , mode(bootstrap) level(`level') sgl(`sgl') acc(`acc') obsb(`obsb') seacc(`seacc') aseobs(`wdbaseobs_out')
  local wdbdfr0 = e(df_r)
  * O-3 (2026-09-10): the analytic residual df is preserved WHATEVER bootdf() says. It used
  * to be kept only on the branch below - the bootdf(normal) case that blanks e(df_r) - so
  * the common clustered run silently lost it, while three help sites promised it
  * unconditionally "alongside the analytic covariances". The covariances are kept
  * unconditionally; now so is this. Pure addition to e(); no reported number moves.
  if `wdbdfr0' < . ereturn hidden scalar df_r_analytic = `wdbdfr0'
  if `wdbdf' < . {
    ereturn scalar boot_df = `wdbdf'
    * align e(df_r) with the reference distribution actually used, so that test and
    * lincom run by hand on the bootstrap covariance reproduce the reported p-values
    ereturn scalar df_r = `wdbdf'
  }
  else {
    * the reported statistic is z / chi-squared here, so a finite e(df_r) inherited from
    * the analytic fit would send a hand-run test or lincom to a t/F metric on the
    * BOOTSTRAP covariance and disagree with the printed p (measured on a 12-cluster run:
    * test reported 1.781e-11 against a printed p of 0 to 14 decimal places). blanking it
    * sends test and lincom to chi-squared/z, which reproduces the printed values; the
    * analytic value is kept as a hidden twin, the way e(V_analytic) preserves the
    * analytic covariance
    ereturn scalar df_r = .
  }
  ereturn scalar boot_N_V = `wdb_nlw'
  ereturn scalar boot_v_posted = `wdb_vposted'
  if "`rseed'" != "" ereturn local boot_seed `rseed'
  * bias = replicate mean - observed (official convention)
  tempname wdbbias
  matrix `wdbbias' = `bmean' - `obsb'
  matrix colnames `wdbbias' = `obn'
  ereturn matrix bias_boot = `wdbbias'

  if "`mode'" == "bayes" ereturn local boot_alpha `bootalpha'
  tempname bmn
  matrix `bmn' = `bmean'
  matrix colnames `bmn' = `obn'
  ereturn matrix b_boot = `bmn'
  tempname bcn
  matrix `bcn' = `bcnt'
  matrix colnames `bcn' = `obn'
  ereturn matrix boot_navail = `bcn'
  if "`comp'" != "" {
    tempname bcm
    matrix `bcm' = `comp'
    cap ereturn hidden matrix boot_comparator = `bcm'
  }
  * above 10,000 COMPLETED replications e(bootdist) is not posted at all (not truncated),
  * so reps(10005) with six failures still posts while reps(20000) leaves it absent. The
  * constant is a policy cap on the size of e(), not a Stata ceiling - J(12000, 3, .) is
  * legal under Stata 19 MP - so if it is ever revisited, revisit it on a measurement of
  * e() size rather than on an assumed matrix limit. bootsave() writes from the
  * accumulators and is NOT subject to the cap, which is the route for a user who needs
  * the draws from a very large run
  if `done' <= 10000 {
    tempname bdist
    if `wdbts' {
      * both accumulators compacted on ONE mask so the two matrices stay row-aligned
      tempname bdists
      mata: wooldid_bootdist("`acc'", "`bdist'", "`seacc'", "`bdists'")
      cap matrix colnames `bdist' = `obn'
      cap ereturn matrix bootdist = `bdist'
      cap matrix colnames `bdists' = `obn'
      cap ereturn matrix bootdist_se = `bdists'
    }
    else {
      mata: wooldid_bootdist("`acc'", "`bdist'")
      cap matrix colnames `bdist' = `obn'
      cap ereturn matrix bootdist = `bdist'
    }
  }
  if `blank' > 0 dis as text "Caution: `blank' replication cell(s) in ame-family columns returned exact structural zeros (resamples without dose variation in some cells) and were excluded from those columns' bootstrap moments."
  if "`wdb_lowav'" != "" dis as text "Caution: the following estimates were available in fewer than half of the completed replications; their bootstrap standard errors rest on the listed counts: `wdb_lowav'"
  if "`wdb_exclnames'" != "" dis as text "Caution: the following estimates have no identified observed value and are excluded from bootstrap inference: `wdb_exclnames'"
end


* ============================================================================
* F7 subprograms: Roth (2022) pre-trend power diagnostics.
* Implemented as subprograms because wooldid_worker sits close to Stata's
* per-program compiled-size ceiling (verified: sizable inline additions make
* the ado fail to load with r(1000)).
* ============================================================================

* parse-time validation for pretrendpower()
cap prog drop wooldid_ptpcheck
prog def wooldid_ptpcheck
version 16.0
  syntax , [levels(numlist) esf(string) impv(string) espre(string)]
  if "`impv'" != "" {
    dis as error "Error: pretrendpower() is not compatible with impvalues (no event-study variance-covariance matrix exists under Ibragimov-Muller inference)."
    error 198
  }
  if "`esf'" == "" {
    dis as error "Error: pretrendpower() requires esfixedbaseperiod. Under the pooled reference the comparison period is a cohort-specific pool of pre-periods, so the single omitted reference period that the Roth (2022) linear-trend normalization requires is not defined."
    dis as error "Specify esfixedbaseperiod (with esrelativeto() if you want a reference other than -1)."
    error 198
  }
  if `espre' <= 0 {
    dis as error "Error: pretrendpower() requires esprelength() of at least 1; without pre-treatment event-study coefficients there is nothing to diagnose."
    error 198
  }
  foreach ptpl of local levels {
    if `ptpl' > 0 & `ptpl' < 1 {
      dis as error "Error: pretrendpower() takes power levels in percent (for example 80, not .8)."
      error 198
    }
    if `ptpl' != int(`ptpl') | `ptpl' <= 5 | `ptpl' >= 100 {
      dis as error "Error: pretrendpower() power levels must be integers strictly between 5 and 100 (for example 50, 80, 90, 95)."
      error 198
    }
  }
end

* compute + display + ereturn the diagnostics. Runs AFTER the ereturn post, so it may
* only use r-class data commands and matrix algebra (an estimation command here would
* wipe the posted e()). All inputs arrive as passed tempvar names or e() objects.
cap prog drop wooldid_ptp
prog def wooldid_ptp, eclass
version 16.0
  syntax , levels(numlist) [esrel(string) dfr(string) pois(string) poisexp(string) tvar(string) ttrevar(string) ivar(string) wvar(string) cwm(string) hba(string) hbi(string) HBARonly level(real 95)]

  * duplicate power levels (e.g. pretrendpower(80 80)) would otherwise be processed
  * twice, yielding a duplicated row in e(roth_power_*) and in the displayed table;
  * keep the first occurrence of each level, preserving the order given
  local levels : list uniq levels

  local ptpref = `esrel'
  * degrees of freedom: use the regression's residual df when available (t and F
  * distributions, matching the critical values wooldid's own tables and joint tests
  * use); fall back to normal / chi-squared otherwise (e.g. some poisson paths)
  local ptpdfok 0
  if "`dfr'" != "" & "`dfr'" != "." {
    if real("`dfr'") > 0 & real("`dfr'") < . local ptpdfok 1
  }
  if `ptpdfok' local ptpcrit = invttail(`dfr', (100 - `level')/200)
  else local ptpcrit = invnormal(1 - (100 - `level')/200)
  * BS-11 (bug sweep 2026-09-06, author ruling 2026-09-07, "pin the pre-test at five
  * percent"): the NIS (any-pre-period-individually-significant) pre-test gets its own
  * critical value, level()-independent. ptpcrit above is a CONFIDENCE level and belongs
  * to the slope interval (F044); the NIS threshold is a TEST size, and the mdt_wald leg
  * below is hard-coded at .95, so letting level() move one and not the other powered the
  * two detectable-trend columns of one table at different sizes. Written as (100 - 95)/200
  * rather than .025 so that at c(level) = 95 the expression is character-identical to
  * ptpcrit's and the default-level baselines cannot shift by an ulp.
  if `ptpdfok' local ptpniscrit = invttail(`dfr', (100 - 95)/200)
  else local ptpniscrit = invnormal(1 - (100 - 95)/200)

  * hbar: estimand-weighted mean post-treatment relative period over the FULL post
  * support of the estimation sample - the pooled estimands aggregate over all post
  * cells (subpop t >= onset), not only the reported event-study window, and using the
  * reported window instead understates the implied bias whenever espostlength is
  * smaller than the largest post relative period
  * BOOTSTRAP PATH: hba()/hbi() carry the two hbar scalars already computed at worker
  * level, so the whole data block below is skipped. It has to be, because `ttrevar' and
  * `wvar' are worker TEMPVARS (ado ~437 and ~348) that Stata drops when the worker
  * exits - and under bootstrap this program is called from the shell, afterwards, so
  * that it reads the bootstrap e(V_es_*) rather than the analytic one it would have
  * seen inside the worker. Everything below this block is pure algebra on
  * e(b_es_*)/e(V_es_*) plus these two scalars
  if "`hba'" != "" | "`hbi'" != "" {
    local ptphbar_att = cond("`hba'" == "" | "`hba'" == ".", ., real("`hba'"))
    local ptphbar_it  = cond("`hbi'" == "" | "`hbi'" == ".", ., real("`hbi'"))
  }
  else {
  tempvar ptprely ptpwa ptppost ptpnum
  qui gen double `ptprely' = `tvar' - `ttrevar'
  qui gen double `ptpwa' = `wvar'
  if "`cwm'" != "" qui replace `ptpwa' = `ptpwa' * `cwm'
  qui gen byte `ptppost' = e(sample) & `ttrevar' != . & `ptprely' >= 0 & `ptprely' < .
  qui gen double `ptpnum' = `ptpwa' * `ptprely' if `ptppost'
  qui sum `ptpnum' if `ptppost', meanonly
  local ptphnum = r(sum)
  qui sum `ptpwa' if `ptppost', meanonly
  local ptphbar_att = cond(r(sum) > 0 & r(sum) < ., `ptphnum'/r(sum), .)
  * att_it-family hbar: weight divided by the (i,t) cell weight total, replicating the
  * weight-construction block (cell totals use the raw weight, without any custom
  * weight multiplier - the documented asymmetry)
  tempvar ptpwpit ptpwit ptpnum2
  qui gegen double `ptpwpit' = total(cond(e(sample), `wvar', 0)), by(`ivar' `tvar')
  qui gen double `ptpwit' = `ptpwa'/`ptpwpit'
  qui replace `ptpwit' = 0 if `ptpwit' == .
  qui gen double `ptpnum2' = `ptpwit' * `ptprely' if `ptppost'
  qui sum `ptpnum2' if `ptppost', meanonly
  local ptphnum2 = r(sum)
  qui sum `ptpwit' if `ptppost', meanonly
  local ptphbar_it = cond(r(sum) > 0 & r(sum) < ., `ptphnum2'/r(sum), .)
  cap drop `ptprely' `ptpwa' `ptppost' `ptpnum' `ptpwpit' `ptpwit' `ptpnum2'
  }
  * stash them so the bootstrap path can hand them back in (hidden: an implementation
  * detail, not part of the documented surface). `ereturn repost' preserves e(), so these
  * survive the bootstrap rewiring that happens between the worker and the shell call
  cap ereturn hidden scalar wdid_ptphbar_att = `ptphbar_att'
  cap ereturn hidden scalar wdid_ptphbar_it  = `ptphbar_it'
  * hbaronly: the worker's pass under bootstrap. The two scalars need the worker's
  * tempvars, so they must be computed HERE; the power algebra needs the bootstrap
  * covariance, which does not exist yet, so it must NOT be computed here
  if "`hbaronly'" != "" exit

  local ptpheader 0
  foreach ptpty in att att_it ame ame_it {
    * on a bootstrap run, refuse to run on a family whose bootstrap covariance swap
    * failed: e(V_es_`ptpty') would still be the ANALYTIC matrix. The success marker is
    * the hidden analytic twin, which wooldid_bootwire creates only when the swap worked
    * (families that were never estimated are skipped silently first: the caution is for
    * a family that exists but whose bootstrap swap failed - review 2026-09-05)
    cap confirm matrix e(b_es_`ptpty')
    if _rc != 0 continue
    if "`e(boot_method)'" != "" {
      cap confirm matrix e(V_es_analytic_`ptpty')
      if _rc != 0 {
        dis as text "Caution: the bootstrap event-study covariance is unavailable for `ptpty'; pretrendpower() diagnostics skipped for `ptpty' rather than computed on the analytic covariance."
        continue
      }
    }
    cap confirm matrix e(b_es_`ptpty')
    if _rc != 0 continue
    cap confirm matrix e(V_es_`ptpty')
    if _rc != 0 continue
    tempname ptpb ptpV
    matrix `ptpb' = e(b_es_`ptpty')
    matrix `ptpV' = e(V_es_`ptpty')
    local ptplbl : colnames `ptpb'
    * pre-period indices and the trend vector m = relyr - reference. m is built from
    * the actual labels: under esrelativeto(< -1) the pre relyrs straddle the
    * reference, so m carries BOTH signs - no ordering or sign assumption is safe
    local ptpk 0
    local ptppreidx
    local ptpmvals
    local ptpj 0
    foreach ptpl of local ptplbl {
      local ++ptpj
      if substr("`ptpl'", 1, 3) == "pre" {
        local ++ptpk
        local ptppreidx `ptppreidx' `ptpj'
        local ptpmvals `ptpmvals' `=-real(substr("`ptpl'", 4, .)) - `ptpref''
      }
    }
    if `ptpk' == 0 {
      dis as text "Caution: pretrendpower() requested, but the `ptpty' event study has no pre-treatment coefficients; diagnostics skipped for `ptpty'."
      continue
    }
    tempname ptpS ptpbp ptpm
    matrix `ptpS' = J(`ptpk', `ptpk', 0)
    matrix `ptpbp' = J(`ptpk', 1, 0)
    matrix `ptpm' = J(`ptpk', 1, 0)
    local ptpr 0
    foreach ptpi of local ptppreidx {
      local ++ptpr
      matrix `ptpbp'[`ptpr', 1] = `ptpb'[1, `ptpi']
      local ptpmv : word `ptpr' of `ptpmvals'
      matrix `ptpm'[`ptpr', 1] = `ptpmv'
      local ptpc 0
      foreach ptpi2 of local ptppreidx {
        local ++ptpc
        matrix `ptpS'[`ptpr', `ptpc'] = `ptpV'[`ptpi', `ptpi2']
      }
    }
    * screens, in hazard order: a missing element would make invsym abort the run; a
    * non-positive diagonal (degeneracy sentinel, zero-variance row) would otherwise be
    * silently g-inversed into a plausible-looking wrong number
    local ptpbad 0
    forval ptpr = 1/`ptpk' {
      if missing(`ptpbp'[`ptpr', 1]) local ptpbad 1
      if missing(`ptpS'[`ptpr', `ptpr']) | `ptpS'[`ptpr', `ptpr'] <= 0 local ptpbad 1
      forval ptpc = 1/`ptpk' {
        if missing(`ptpS'[`ptpr', `ptpc']) local ptpbad 1
      }
    }
    if `ptpbad' {
      dis as text "Caution: the `ptpty' pre-period variance-covariance block contains missing or non-positive variance entries; pretrendpower() diagnostics skipped for `ptpty'."
      continue
    }
    matrix `ptpS' = (`ptpS' + `ptpS'')/2
    tempname ptpSi ptpwork
    matrix `ptpSi' = invsym(`ptpS')
    local ptpkeff = `ptpk' - diag0cnt(`ptpSi')
    if `ptpkeff' <= 0 {
      dis as text "Caution: the `ptpty' pre-period variance-covariance block is not invertible; pretrendpower() diagnostics skipped for `ptpty'."
      continue
    }
    if `ptpkeff' < `ptpk' dis as text "Caution: the `ptpty' pre-period variance-covariance block is rank-deficient (rank `ptpkeff' of `ptpk'). The slope, its CI, the joint test (now `ptpkeff' degrees of freedom), and the detectable-trend columns then reflect only the linearly independent subset of pre-period coefficients a generalized inverse retains - NOT a summary of all pre-period evidence - and can differ arbitrarily from statistics built on the discarded coefficients; interpret with care."
    matrix `ptpwork' = `ptpm''*`ptpSi'*`ptpm'
    local ptpc2 = `ptpwork'[1, 1]
    if missing(`ptpc2') | `ptpc2' <= 0 {
      dis as text "Caution: the linear-trend information m'V^(-1)m for `ptpty' is not positive; pretrendpower() diagnostics skipped for `ptpty'."
      continue
    }
    matrix `ptpwork' = `ptpm''*`ptpSi'*`ptpbp'
    local ptpslope = `ptpwork'[1, 1]/`ptpc2'
    local ptpse = 1/sqrt(`ptpc2')
    local ptpll = `ptpslope' - `ptpcrit'*`ptpse'
    local ptpul = `ptpslope' + `ptpcrit'*`ptpse'
    matrix `ptpwork' = `ptpbp''*`ptpSi'*`ptpbp'
    local ptpW = `ptpwork'[1, 1]
    if `ptpdfok' local ptpprep = Ftail(`ptpkeff', `dfr', `ptpW'/`ptpkeff')
    else local ptpprep = chi2tail(`ptpkeff', `ptpW')
    local ptphbar `ptphbar_att'
    if "`ptpty'" == "att_it" | "`ptpty'" == "ame_it" local ptphbar `ptphbar_it'

    * NIS setup: the pre-test Roth's pretrends package powers - reject when ANY pre
    * coefficient is individually significant at the 5% level. BS-11 (author ruling
    * 2026-09-07): that 5% is deliberate and level()-independent - it matches the .95
    * hard-coded in the mdt_wald leg below and Roth's own pretrends convention, so the
    * two detectable-trend columns power tests of the SAME size. level() reaches the
    * slope interval (ptpcrit) and nothing here; the threshold below is ptpniscrit
    * exact quadrature cost explodes with dimension (~2.9s per evaluation at k=6, ~85
    * evaluations per power level per type); cap the NIS column at k <= 4 and point
    * larger designs at the pretrends package, which consumes e(b_es_*)/e(V_es_*)
    local ptpnisok = (`ptpk' <= 4)
    local ptpnissize .
    if `ptpnisok' {
      tempname ptpsd ptpmrow
      matrix `ptpsd' = J(1, `ptpk', 0)
      forval ptpr = 1/`ptpk' {
        matrix `ptpsd'[1, `ptpr'] = sqrt(`ptpS'[`ptpr', `ptpr'])
      }
      matrix `ptpmrow' = `ptpm''
      cap mata: st_local("ptpnissize", strofreal(1 - mvnormalcvqp(-`ptpniscrit'*st_matrix("`ptpsd'"), `ptpniscrit'*st_matrix("`ptpsd'"), J(1, `ptpk', 0), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
      if _rc != 0 local ptpnisok 0
      else if missing(real("`ptpnissize'")) local ptpnisok 0
    }

    local ptpnlev : word count `levels'
    tempname ptpmat
    matrix `ptpmat' = J(`ptpnlev', 5, .)
    local ptprn
    local ptplr 0
    foreach ptplev of local levels {
      local ++ptplr
      local ptprn `ptprn' p`ptplev'
      matrix `ptpmat'[`ptplr', 1] = `ptplev'
      * Wald-based minimum detectable trend: powers the joint pre-trend test above,
      * using the F distribution with the regression's residual df (a chi-squared
      * shortcut is optimistic by ~25 percent at 10 clusters)
      if `ptpdfok' local ptplam = npnF(`ptpkeff', `dfr', invF(`ptpkeff', `dfr', .95), 1 - `ptplev'/100)
      else local ptplam = npnchi2(`ptpkeff', invchi2(`ptpkeff', .95), 1 - `ptplev'/100)
      if !missing(`ptplam') {
        matrix `ptpmat'[`ptplr', 2] = sqrt(`ptplam'/`ptpc2')
        if !missing(`ptphbar') matrix `ptpmat'[`ptplr', 4] = sqrt(`ptplam'/`ptpc2')*(`ptphbar' - `ptpref')
      }
      * NIS-based minimum detectable trend, by bisection on the exact rejection
      * probability (deterministic quadrature; self-checked against the target)
      if `ptpnisok' {
        * NFPTP-04: a requested power level at or below the NIS test's null rejection
        * rate (its size) is achieved by a zero trend - the mathematically exact boundary
        * answer. formerly these cells were left missing, indistinguishable from a
        * numerical failure
        if `ptplev'/100 <= real("`ptpnissize'") {
          matrix `ptpmat'[`ptplr', 3] = 0
          if !missing(`ptphbar') matrix `ptpmat'[`ptplr', 5] = 0
        }
        if `ptplev'/100 > real("`ptpnissize'") {
          local ptpglo = 0
          local ptpghi = 8*`ptpse'
          local ptpiter 0
          local ptppw 0
          while `ptppw' < `ptplev'/100 & `ptpiter' < 60 {
            local ptpghi = 2*`ptpghi'
            cap mata: st_local("ptppw2", strofreal(1 - mvnormalcvqp(-`ptpniscrit'*st_matrix("`ptpsd'"), `ptpniscrit'*st_matrix("`ptpsd'"), `ptpghi'*st_matrix("`ptpmrow'"), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
            if _rc != 0 local ptpiter 999
            else local ptppw = real("`ptppw2'")
            local ++ptpiter
          }
          if `ptpiter' < 999 {
            forval ptpbi = 1/80 {
              local ptpgmid = (`ptpglo' + `ptpghi')/2
              cap mata: st_local("ptppw2", strofreal(1 - mvnormalcvqp(-`ptpniscrit'*st_matrix("`ptpsd'"), `ptpniscrit'*st_matrix("`ptpsd'"), `ptpgmid'*st_matrix("`ptpmrow'"), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
              if _rc == 0 {
                if real("`ptppw2'") < `ptplev'/100 local ptpglo = `ptpgmid'
                else local ptpghi = `ptpgmid'
              }
            }
            local ptpgn = (`ptpglo' + `ptpghi')/2
            cap mata: st_local("ptppw2", strofreal(1 - mvnormalcvqp(-`ptpniscrit'*st_matrix("`ptpsd'"), `ptpniscrit'*st_matrix("`ptpsd'"), `ptpgn'*st_matrix("`ptpmrow'"), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
            if _rc == 0 {
              if abs(real("`ptppw2'") - `ptplev'/100) < 1e-4 {
                matrix `ptpmat'[`ptplr', 3] = `ptpgn'
                if !missing(`ptphbar') matrix `ptpmat'[`ptplr', 5] = `ptpgn'*(`ptphbar' - `ptpref')
              }
            }
          }
        }
      }
    }
    matrix rownames `ptpmat' = `ptprn'
    matrix colnames `ptpmat' = power mdt_wald mdt_nis bias_wald bias_nis

    * display
    if `ptpheader' == 0 {
      dis as result " "
      dis as result "Pre-Trend Power Diagnostics (Roth 2022):"
      local ptpheader 1
    }
    dis as result "  `ptpty': GLS linear pre-trend fit on `ptpk' pre-treatment coefficient(s); reference relative period = `ptpref'"
    dis as result "    slope = " %-9.0g `ptpslope' " (se " %-9.0g `ptpse' "; `level'% CI " %-9.0g `ptpll' " to " %-9.0g `ptpul' ")"
    if `ptpdfok' dis as result "    joint pre-trend test: F(`ptpkeff', `dfr') = " %-8.0g `=`ptpW'/`ptpkeff'' ", p = " %-6.4f `ptpprep'
    else dis as result "    joint pre-trend test: chi2(`ptpkeff') = " %-8.0g `ptpW' ", p = " %-6.4f `ptpprep'
    if `ptpnisok' dis as result "    size of the any-pre-period-individually-significant pre-test: " %-6.4f real("`ptpnissize'")
    if !`ptpnisok' & `ptpk' > 4 dis as text "    note: mdt_nis and the NIS-test size are omitted (`ptpk' pre-periods exceed the exact-quadrature budget of 4); mdt_wald is unaffected. For the NIS convention at larger k, feed e(b_es_`ptpty')/e(V_es_`ptpty')/e(numpre_`ptpty') to Roth's pretrends package."
    dis as result "    weighted mean post-treatment relative period (hbar) = " %-9.0g `ptphbar'
    dis as result "    minimum detectable linear trend, and the bias such a trend would imply for the"
    dis as result "    `ptpty' pooled estimand if it continued into the post-treatment period:"
    * NFPTP-03: the att_it/ame_it block is the shared event-study curve; when the user
    * requested att_i/att_itime the block stands in for them but its hbar/bias columns
    * use att_it-family post-period weights - say so rather than implying the requested
    * pooled estimand's own weighting. exact-token match (att_i is a prefix of att_it
    * and att_itime; cf. the cfxplottypes exact-token precedent)
    if inlist("`ptpty'", "att_it", "ame_it") {
      local ptpsurr
      local ptpests `e(wdid_estimands)'
      if `: list posof "att_i" in ptpests' > 0 local ptpsurr att_i
      if `: list posof "att_itime" in ptpests' > 0 local ptpsurr `ptpsurr' att_itime
      if "`ptpsurr'" != "" {
        if "`ptpty'" == "ame_it" local ptpsurr = subinstr("`ptpsurr'", "att", "ame", .)
        dis as text "    note: this `ptpty'-weighted diagnostic stands in for the requested `ptpsurr' estimand(s); its hbar and bias columns use `ptpty' post-period weights."
      }
    }
    matlist `ptpmat', format(%12.4g)
    dis as text "    note: mdt_wald powers the joint pre-trend test above; mdt_nis powers the test that any"
    dis as text "    pre-period coefficient is individually significant (the convention of Roth's pretrends package)."
    dis as text "    Both pre-tests are at the 5% level, whatever level() is in force; level() moves the slope CI only."
    if "`pois'" == "poisson" & "`poisexp'" == "" dis as text "    note: poisson event-study coefficients are on the linear-index (xb) scale, and so are these diagnostics."

    * ereturns (matrix uses copy so the object survives if this block is ever resited)
    cap ereturn matrix roth_power_`ptpty' = `ptpmat', copy
    cap ereturn scalar roth_slope_`ptpty' = `ptpslope'
    cap ereturn scalar roth_slope_se_`ptpty' = `ptpse'
    cap ereturn scalar roth_slope_ll_`ptpty' = `ptpll'
    cap ereturn scalar roth_slope_ul_`ptpty' = `ptpul'
    cap ereturn scalar roth_prep_`ptpty' = `ptpprep'
    cap ereturn scalar roth_hbar_`ptpty' = `ptphbar'
    cap ereturn scalar roth_npre_`ptpty' = `ptpk'
    cap ereturn scalar roth_keff_`ptpty' = `ptpkeff'
    if `ptpnisok' cap ereturn scalar roth_nissize_`ptpty' = real("`ptpnissize'")
    * F11: the joint-test statistic, so replay can re-render the F/chi2 line from e()
    cap ereturn scalar roth_W_`ptpty' = `ptpW'
  }
  if `ptpheader' == 0 dis as text "Caution: pretrendpower() was requested, but no event-study coefficient vector with a valid variance-covariance matrix was available (e.g. under timetrends); no diagnostics were computed."
  * completion marker: the caller runs this program under cap noi and every ereturn
  * above is individually cap'd, so an abort partway through the type loop would leave
  * the earlier types' roth_* returns in e() with nothing marking the run incomplete.
  * e(roth_ok) is set only here, after the loop has finished; any roth_* return present
  * without e(roth_ok) therefore identifies a partially populated diagnostics pass.
  cap ereturn scalar roth_ok = 1
end


* F11: all result presentation lives here, out of wooldid_worker (which is near stata's
* per-program compiled-size ceiling), reading ONLY e() - the live estimation path and
* replay (wooldid typed alone after estimation) both call this, so their output is
* identical by construction. non-eclass on purpose: ereturn inside a non-eclass program
* is a hard error, so this program is structurally incapable of altering the posted
* results. part(1) = header through the adversarial lines; part(2) = subgroup tables
* through the pairwise comparisons (split so the eclass wooldid_ptp call can stay in
* the caller between them); part(roth) redisplays the Roth diagnostics from e(roth_*)
* on replay without re-running the quadrature.
* NFFAST-02: rename CFX matrix row stripes from margins' tempvar-embedded subnames
* (e.g. "1bn.__000ABC" - volatile every run) to stable lattice_# keys, and null the row
* equations. lives outside wooldid_worker (compiled-size ceiling). the lattice level is
* read from the matrix's own row subnames - no r() capture (margins' r(by#) is not a
* documented stored result). validation guarantees nonnegative integer levels, so the
* generated names are always valid; a parse failure is loud, never silent
cap prog drop wooldid_cfxstripe
prog def wooldid_cfxstripe
version 16.0
  args cfxm
  local rn : rowfullnames `cfxm'
  local newnames
  foreach rtok of local rn {
    local rsub `rtok'
    if strpos("`rtok'", ":") > 0 local rsub = substr("`rtok'", strpos("`rtok'", ":") + 1, .)
    if !regexm("`rsub'", "^([0-9]+)(bn|b|o)?\.") {
      dis as error "Error: wooldid could not parse the CFX matrix row name '`rtok'' into a lattice level; please report this."
      error 498
    }
    local newnames `newnames' lattice_`=regexs(1)'
  }
  matrix rownames `cfxm' = `newnames'
  matrix roweq `cfxm' = _
end

* PS-3 (author ruling 2026-09-07, built 2026-09-10): "keep the table structure, give each
* column an explicit format". matlist's default is ONE general format for the whole matrix,
* which is why a p of 1.84e-21 used to sit beside a t of 13.07929 where every comparator
* prints 0.000 and 13.08. This builds matlist's cspec() and rspec() from the matrix's OWN
* names, so no call site hard-codes a column count and a table that gains a column cannot
* break: statistics get two decimals, p-values three, the relative-year column prints as the
* integer it is, and everything else keeps a general format so no estimate loses precision.
* Note the exact-match rule on p-values - the summary-statistics table has columns p25 and
* p75, which are percentiles and must stay on the general format.
* rspec() is required with cspec() and takes one character per horizontal rule: a "-" draws
* one, a "&" does not, and the count is rows + 2 (measured on all three tables, and
* independent of how many row-equation groups they carry). The rules are placed where the
* default put them - under the header, at every equation-group boundary, and at the foot -
* so the structure is unchanged and only the number formatting moves.
cap prog drop wooldid_mlfmt
prog def wooldid_mlfmt
version 16.0
  syntax , mat(string)
  local wdmrn : rowfullnames `mat'
  * THE RESULTS COME BACK AS TWO FIXED CALLER LOCALS, wdmlfmt and wdmlfmtr, which the
  * caller cannot rename - so every `matlist X, cspec(`wdmlfmt') rspec(`wdmlfmtr')' must be
  * immediately preceded by this program's call FOR THAT MATRIX. Reusing a previous table's
  * macros is not a cosmetic slip: MEASURED on Stata 19 MP, a cspec whose column count is
  * short gives r(198) "& or | expected", too long gives "extraneous input in cspec()", and
  * a short rspec gives "too few specifications in rspec()" - so a mispaired format string
  * aborts a COMPLETED estimation mid-display, and the shell's WOOLDID_LASTRC marker then
  * refuses the replay. The widths below: 12 is matlist's own default minimum row-label
  * width, 26 the cap that keeps the label column from pushing the statistics off the line
  * (matlist abbreviates a longer name with a tilde rather than truncating). Only the text
  * after ":" is measured on equation-prefixed rownames, because matlist prints the equation
  * name on a line of its own
  local wdmw 12
  foreach wdmn of local wdmrn {
    local wdmnm = cond(strpos("`wdmn'", ":") > 0, substr("`wdmn'", strpos("`wdmn'", ":") + 1, .), "`wdmn'")
    if strlen("`wdmnm'") > `wdmw' local wdmw = strlen("`wdmnm'")
    if strlen("`wdmn'") > `wdmw' & strpos("`wdmn'", ":") == 0 local wdmw = strlen("`wdmn'")
  }
  if `wdmw' > 26 local wdmw 26
  local wdmf "& %-`wdmw's |"
  foreach wdmc in `: colnames `mat'' {
    local wdmc1 = lower("`wdmc'")
    if "`wdmc1'" == "relyr" | substr("`wdmc1'", 1, 6) == "numpre" local wdmf `wdmf' %7.0f &
    else if "`wdmc1'" == "p" | "`wdmc1'" == "pvalue" | substr("`wdmc1'", 1, 2) == "p_" local wdmf `wdmf' %8.3f &
    else if "`wdmc1'" == "t" | "`wdmc1'" == "z" | substr("`wdmc1'", 1, 2) == "t_" | substr("`wdmc1'", 1, 2) == "z_" local wdmf `wdmf' %8.2f &
    else local wdmf `wdmf' %10.0g &
  }
  * character 1 governs the rule ABOVE the header, which the default does not draw; the next
  * one sits under the header, which it does; then one per row, ruled at each equation-group
  * boundary; and the last is the foot, which the default leaves open
  local wdmr "&-"
  local wdmprev "~~none~~"
  local wdmfirst 1
  foreach wdmn of local wdmrn {
    local wdmeq = cond(strpos("`wdmn'", ":") > 0, substr("`wdmn'", 1, strpos("`wdmn'", ":") - 1), "")
    if `wdmfirst' local wdmfirst 0
    else if "`wdmeq'" != "`wdmprev'" local wdmr "`wdmr'-"
    else local wdmr "`wdmr'&"
    local wdmprev "`wdmeq'"
  }
  local wdmr "`wdmr'&"
  c_local wdmlfmt `"`wdmf'"'
  c_local wdmlfmtr `"`wdmr'"'
end

* F038 (quality review 2026-09-05, author ruling 2026-09-06): a rank-based bootstrap p of
* exactly 0 means no completed replication fell on the far side of the null - the
* resolution limit 2/n, not a probability of 0. matlist renders a numeric matrix, so the
* cell keeps the 0 (e() keeps it too, by the ruling); this line under the table says what
* the 0 means. Under permute() only the plain c/n share (permrawp) can print 0 from a
* sampled enumeration - the add-one conventions never do, and an exact enumeration's 0
* would be a true share (bug sweep BS-25); the Wald display's p is a tail probability of the t/normal
* reference rather than a rank, so it has no floor and exits at the top of the program.
* The boottstat pass-through in that same gate is vacuous by F033: the studentized family
* is add-one with floor 2/(n+1) and can never print 0, so the 2/n message would be wrong
* if it ever fired.
cap prog drop wooldid_pfloor
prog def wooldid_pfloor
  args mat
  if "`e(boot_completed)'" != "" {
    if "`e(boot_display)'" == "wald" & "`e(boot_tstat)'" != "1" exit
    local wdpn = e(boot_completed)
    local wdpres = 2/`wdpn'
    local wdpmsg (2/`wdpn' completed replications, the resolution of the rank-based p): no completed replication fell on the far side of the null. The cell and e() keep the 0.
  }
  else if "`e(perm_completed)'" != "" {
    if "`e(perm_pconv)'" != "raw" | "`e(perm_exact)'" == "1" exit
    local wdpn = e(perm_completed)
    local wdpres = 1/`wdpn'
    local wdpmsg (1/`wdpn' completed replications, the resolution of the plain c/n share): no completed replication reached the observed absolute statistic. The cell and e() keep the 0; e(perm_p_ci) holds the Clopper-Pearson upper bound.
  }
  else exit
  cap confirm matrix `mat'
  if _rc != 0 exit
  tempname pm
  matrix `pm' = `mat'
  local pc = colnumb(`pm', "p")
  * the fallback to column 4 is the LAYOUT CONSTANT of the three tables that reach this
  * program - estimate, se, t|z, p, lb_<lvl>ci, ub_<lvl>ci [, relyr] - in which p is column
  * 4 and stays column 4 under the bootstrap and permute renames (those touch only
  * se -> se_analytic and t|z -> t_pt / t_perm). It is kept for a store whose colnames were
  * lost. Any NEW caller must carry that layout or pass a matrix with a real "p" colname:
  * on a mislabelled matrix this counts zeros in whatever column 4 holds and prints the
  * "p = 0 means p < ..." annotation under a table whose zeros are not p-values, with no
  * error. The `pc' > colsof' exit below catches a short matrix, never a mislabelled one
  if `pc' >= . local pc 4
  if `pc' > colsof(`pm') exit
  mata: st_local("wdp0", strofreal(sum(st_matrix("`pm'")[., `pc'] :== 0)))
  if `wdp0' > 0 dis as text "  p = 0 above means p < " %5.3f `wdpres' " `wdpmsg'"
end

* ============================================================================
* BS-10 (author ruling 2026-09-07): reduce the per-row Ibragimov-Muller degrees of
* freedom collected at the six IM ttest sites and align e(df_r) with them. eclass so
* the ereturns survive the return (same reason wooldid_ptp is eclass); called only
* after wooldid_worker's ereturn post, so it adds to an existing e() rather than
* replacing it. Kept out of wooldid_worker for compiled-size headroom.
* ana() is a string, not an integer option, so a poisson fit with no e(df_r) can pass "."
cap prog drop wooldid_imdf
prog def wooldid_imdf, eclass
version 16.0
  syntax , dfl(string) [ana(string)]
  local wdimin .
  local wdimax .
  mata: st_local("wdimin", strofreal(min(strtoreal(tokens(st_local("dfl"))))))
  mata: st_local("wdimax", strofreal(max(strtoreal(tokens(st_local("dfl"))))))
  * a degenerate reduction leaves e(df_r) exactly as posted - the ruling's "do not post
  * missing" applies to the fallback as much as to the IM value
  if `wdimin' <= 0 | `wdimin' >= . exit
  if "`ana'" != "" & "`ana'" != "." ereturn hidden scalar df_r_analytic = `ana'
  ereturn scalar df_r = `wdimin'
  ereturn scalar imdf_min = `wdimin'
  ereturn scalar imdf_max = `wdimax'
end

* ============================================================================
* PS-62 (six-lens feasibility review; author rulings 2026-09-07). e(V) has always been
* diag(se^2): wooldid estimates no cross-estimate covariance, so a test or lincom spanning
* two of its rows was never valid. The review found that for ONE family of rows that is an
* accounting gap rather than a missing estimate. Every reported relative-period effect of an
* event-study family comes from a SINGLE margins call, so that call's r(V) IS their complete
* joint covariance; wooldid already captures it and posts it as e(V_es_*), where only
* honestdid, pretrendpower() and the joint tests read it. This program copies those genuine
* within-family covariances into the corresponding diagonal BLOCKS of e(V), so that a
* cross-horizon test or lincom becomes valid on the posted matrix.
* Scope is exactly the ruling's: the event-study blocks and nothing else. Every off-diagonal
* between two families, or between a family and a pooled estimate, stays zero - those are the
* cells the review found genuinely blocked (three modes have no delta-method object at all,
* and the pooled effect is a linear combination of the horizons, so the general matrix would
* be rank-deficient by construction rather than by accident). No new estimation happens here.
* Labelling is the ruling's too: help text plus a stored marker - e(wdid_vcvstate) becomes
* "analytic-blocks" and e(wdid_vblocks) names the families actually spliced - and
* deliberately NO footer under the printed table.
* eclass because its ereturns must survive the return, the way wooldid_imdf's do; kept out
* of wooldid_worker for compiled-size headroom.
* ============================================================================
* PS-45 (author ruling 2026-09-07): simultaneous (sup-t / uniform) confidence bands over
* an event study's horizons, DEFAULTED OFF, behind simultaneousbands. The printed table
* keeps its per-row pointwise intervals - a simultaneous band is a property of the family,
* not of a row (same ruling) - so this posts the endpoints for export and the plot layer
* draws them; no shipped figure changes until asked.
* Both engines the ruling names are here, chosen by mode():
*   bootstrap - the upper quantile over replications of the largest studentized deviation,
*               one pass over the accumulator that already feeds the joint test. Called
*               from wooldid_bootwire, after the coefficient tables are rewired.
*   analytic  - simulated from the family's own covariance e(V_es_*). Called from
*               wooldid_worker. Deliberately reads e(V_es_*) and NOT the e(V) blocks the
*               PS-62 splice writes: the two hold the same numbers, but reading the family
*               matrix directly means the band does not depend on whether the splice fired.
* The analytic draw runs on a FIXED internal seed with the caller's RNG state saved and
* restored around it, so the band is reproducible, does not depend on where the user's
* stream happens to be, and does not disturb it.
* Under permute() neither engine runs: validate refuses the combination.
cap prog drop wooldid_esband
prog def wooldid_esband, eclass
version 16.0
  syntax , mode(string) [level(real 95) sgl(string) acc(string) obsb(string) nd(integer 10000) seacc(string) aseobs(string)]
  if `level' <= 0 | `level' >= 100 local level = 95
  * bootwire's family list, so the two engines agree on what a family is
  local wdsfl att att_it ame ame_it
  foreach wdss of local sgl {
    local wdsfl `wdsfl' att_s`wdss' att_it_s`wdss' ame_s`wdss' ame_it_s`wdss'
  }
  * the caller's RNG state is saved here and restored after the loop, so the simulated
  * analytic band never disturbs the user's stream. The SEED itself is set per family inside
  * the loop - see the D-9 note there
  local wdsrng
  if "`mode'" != "bootstrap" local wdsrng `c(rngstate)'
  local wdsdone
  * D-7 (author ruling 2026-09-11): families that were ATTEMPTED and produced no band are
  * collected here and named in a caution below. Before this, a family whose band could not
  * be formed simply vanished from its figure whenever some other family got one - so a user
  * could see a dashed band on the ATT figure, none on the AME_it figure, and be told
  * nothing. An attempted family is one whose inputs are present: the kernels return an
  * EMPTY reason for a family merely absent from this run, and those are not named.
  local wdsnoband
  foreach wdsk of local wdsfl {
    * the e(b)/e(V) row-name prefix for this family (worker: rownamespool, then the _SG
    * injection) - the same rule wooldid_vsplice uses
    local wdsb `wdsk'
    local wdsg
    if regexm("`wdsk'", "^(.+)_s([0-9]+)$") {
      local wdsb = regexs(1)
      local wdsg = regexs(2)
    }
    local wdspx = cond("`wdsg'" == "", "ES_`wdsb'_", "ES_`wdsb'_SG`wdsg'_")
    if strlen("esband_`wdsk'") > 32 continue
    tempname SB
    local wdscrit .
    local wdsn .
    local wdsm 0
    local wdsstud 0
    local wdsbwhy
    if "`mode'" == "bootstrap" {
      * cap: see the note in wooldid_vsplice (cleanmatrices clears wooldid's own Mata
      * kernels; SI-17/SI-18). CONTINUE and not exit - the c(rngstate) restore at the end of
      * this loop must run, or a run that asked for a band would leak the fixed internal seed
      * into the caller's stream. With every family skipped the named caution below fires.
      cap mata: st_local("wdsm", strofreal(wooldid_bootband("`acc'", "`obsb'", "`wdspx'", `level', "`SB'", "wdscrit", "wdsn", "`seacc'", "`aseobs'", "wdsstud", "wdsbwhy")))
      if _rc != 0 continue
    }
    else {
      * D-9 (author ruling 2026-09-11): the seed is set PER FAMILY, derived from the fixed
      * internal seed and the family's own name. It used to be set once before the loop, so
      * each family consumed successive draws and a family's critical value depended on how
      * many families had been banded before it - adding subgroups() or an estimand moved an
      * unrelated family's band in its low-order digits. Deriving the offset from the NAME
      * rather than the loop position is what makes a family's band independent of which
      * other families this run happens to band. Reproducibility is unchanged: the seed is
      * still fixed, and the caller's stream is still saved and restored around the loop
      local wdsseed 20260907
      cap mata: st_local("wdsseed", strtrim(strofreal(20260907 + mod(sum(ascii(st_local("wdsk")) :* (1..strlen(st_local("wdsk")))), 1000000), "%12.0f")))
      set seed `wdsseed'
      cap confirm matrix e(V_es_`wdsk')
      local wdsrc1 = _rc
      cap confirm matrix e(b_es_`wdsk')
      if `wdsrc1' != 0 | _rc != 0 continue
      tempname BV BB
      matrix `BV' = e(V_es_`wdsk')
      matrix `BB' = e(b_es_`wdsk')
      * cap + continue: see the bootstrap twin above
      cap mata: st_local("wdsm", strofreal(wooldid_esbanda("`BV'", "`BB'", `level', `nd', "`SB'", "wdscrit")))
      if _rc != 0 continue
      * the analytic kernel DROPS a degenerate horizon and bands the rest, so it comes back
      * short only when fewer than two bandable horizons are left. Its two input matrices
      * were confirmed just above, so this family really was attempted
      if `wdsm' < 2 local wdsbwhy fewer than two horizons of this family have a usable variance
    }
    if `wdsm' < 2 {
      if `"`wdsbwhy'"' != "" local wdsnoband `"`wdsnoband'`wdsk' (`wdsbwhy'); "'
      continue
    }
    ereturn matrix esband_`wdsk' = `SB'
    if strlen("esband_crit_`wdsk'") <= 32 ereturn scalar esband_crit_`wdsk' = `wdscrit'
    if "`mode'" == "bootstrap" & strlen("esband_n_`wdsk'") <= 32 ereturn scalar esband_n_`wdsk' = `wdsn'
    local wdsdone `wdsdone' `wdsk'
    if `wdsstud' local wdsanystud 1
  }
  if "`wdsrng'" != "" set rngstate `wdsrng'
  if "`wdsdone'" == "" {
    dis as text "Caution: simultaneousbands was requested but no simultaneous band could be built. A band needs at least two horizons of one event-study family with a usable variance - so it is unavailable under fastnose, on a run with no event study, and on a family whose covariance is entirely unidentified. The pointwise intervals in the table are unaffected."
    if `"`wdsnoband'"' != "" dis as text "  Families attempted and declined: `wdsnoband'"
    exit
  }
  * D-7/D-8: some families banded and some did not - say which, and why. The bootstrap
  * kernel refuses a whole family for one unusable horizon (its critical value is a
  * maximum taken replication by replication, so a horizon that cannot supply a
  * denominator at every retained replication changes what the maximum is over), while
  * the analytic kernel drops the horizon and bands the rest. That asymmetry is
  * deliberate (author ruling 2026-09-11), which is why the reason names the HORIZON
  * rather than only the family
  if `"`wdsnoband'"' != "" {
    dis as text "Caution: a simultaneous band was drawn for some event-study families and not others. No band for: `wdsnoband'The figures for those families show their pointwise intervals only, which are unaffected."
  }
  * PS-45 (author ruling 2026-09-10): under boottstat the band is a genuine percentile-t -
  * each replication's own standard error in the denominator, the observed analytic ones in
  * the endpoints - so it is on the SAME scale as the pointwise intervals printed beside it.
  * It used to be scaled by the coefficient-draw standard deviation whatever the display,
  * which meant a simultaneous band could come out NARROWER than the pointwise interval it
  * accompanies. Three engine values now, so a store says which construction produced it.
  ereturn local esband_engine = cond("`mode'" != "bootstrap", "analytic-simulated", cond("`wdsanystud'" == "1", "bootstrap-studentized", "bootstrap"))
  ereturn local esband_families `wdsdone'
  ereturn scalar esband_level = `level'
  if "`mode'" != "bootstrap" ereturn scalar esband_ndraw = `nd'
end

* ============================================================================
* PS-45: the plot layer. Fills __wdsblo/__wdsbhi on the plot dataset from the stored
* band of ONE family and hands back the twoway layer that draws them, empty when the
* family has no band - so the figure is byte-identical to today's on every run that did
* not ask for a band, and on replay of a store written before 2026-09-07. Reading the
* endpoints from e() rather than recomputing them is what makes replay work.
cap prog drop wooldid_esbandv
prog def wooldid_esbandv
version 16.0
  syntax , fam(string) [sfx(string)]
  cap drop __wdsblo
  cap drop __wdsbhi
  c_local wdsblayer
  cap confirm matrix e(esband_`fam'`sfx')
  if _rc != 0 exit
  tempname WSB
  matrix `WSB' = e(esband_`fam'`sfx')
  * THE STORED BAND IS READ ENTIRELY BY ROW POSITION below: row 1 lower endpoint, row 2
  * upper, row 3 the relative year to match on; the columns are the surviving horizons.
  * That layout is produced by wooldid_esbanda and its bootstrap twin wooldid_bootband,
  * which do attach a real rowstripe - this program never reads it, so inserting a row (a
  * critical value, a count) ahead of those three draws a band from the wrong series with
  * no error. This guard and the `cap confirm' above drop the band for THIS FAMILY with no
  * message, while wooldid_esband's caution fires only when NO family got a band at all:
  * under simultaneousbands a user can get a dashed band on the ATT figure, none on the
  * AME_it figure, and be told nothing
  if rowsof(`WSB') < 3 exit
  qui gen double __wdsblo = .
  qui gen double __wdsbhi = .
  forval wsbc = 1/`=colsof(`WSB')' {
  * `sfx' is in the e() KEY above and deliberately not in this row match, because the two
  * namespaces really do differ: wooldid_esband stores bands keyed by family-plus-subgroup
  * (att_s3, att_it_s3), while the plot frame's rownamesvar comes from the row equations of
  * e(wdidsgresults_3), which are NOT suffixed (the _SG3_ injection touches only the stacked
  * e(b)/e(V) names). "Harmonizing" this to ES_`fam'`sfx' matches nothing, leaves both
  * variables all-missing, and makes every subgroup figure's band vanish with no error -
  * the r(N) > 0 test below just returns an empty layer
    qui replace __wdsblo = `WSB'[1, `wsbc'] if rownamesvar == "ES_`fam'" & relyr == `WSB'[3, `wsbc']
    qui replace __wdsbhi = `WSB'[2, `wsbc'] if rownamesvar == "ES_`fam'" & relyr == `WSB'[3, `wsbc']
  }
  qui count if __wdsblo < . & __wdsbhi < .
  if r(N) > 0 c_local wdsblayer (rline __wdsblo __wdsbhi relyr, lpattern(dash) lcolor(navy))
end

cap prog drop wooldid_vsplice
prog def wooldid_vsplice, eclass
version 16.0
  syntax , [sgl(string)]
  cap confirm matrix e(V)
  if _rc != 0 exit
  * the family list is bootwire's, so the two engines agree on what a family is
  local wdvfl att att_it ame ame_it
  foreach wdvs of local sgl {
    local wdvfl `wdvfl' att_s`wdvs' att_it_s`wdvs' ame_s`wdvs' ame_it_s`wdvs'
  }
  tempname WV
  matrix `WV' = e(V)
  local wdvdone
  local wdvnp 0
  foreach wdvk of local wdvfl {
    cap confirm matrix e(V_es_`wdvk')
    if _rc != 0 continue
    * e(V) names an overall family's rows ES_<key>_<horizon> and a subgroup's
    * ES_<base>_SG<sg>_<horizon> (worker: rownamespool, then the _SG injection)
    local wdvb `wdvk'
    local wdvg
    if regexm("`wdvk'", "^(.+)_s([0-9]+)$") {
      local wdvb = regexs(1)
      local wdvg = regexs(2)
    }
    local wdvpx = cond("`wdvg'" == "", "ES_`wdvb'_", "ES_`wdvb'_SG`wdvg'_")
    tempname WE
    matrix `WE' = e(V_es_`wdvk')
    local wdvn 0
    * cleanmatrices runs `clear mata' at the top of the worker, which drops wooldid's own
    * Mata kernels for the rest of the session. SI-17/SI-18 measured that an adopath-resolved
    * ado recovers (Stata re-reads the file when a function is missing) and only a
    * `do <path>' load - this project's harness - does not, which is why the refusal that
    * once guarded this was removed. So the missing kernel is caught rather than raised:
    * the run keeps its analytic-diagonal e(V), exactly the surface it had before PS-62,
    * instead of dying r(3499). Exiting is right, not continuing: if the kernel is gone it
    * is gone for every family, and there is no state to restore here.
    cap mata: st_local("wdvn", strofreal(wooldid_vspl("`WV'", "`WE'", "`wdvpx'")))
    if _rc != 0 exit
    if `wdvn' > 0 {
      local wdvdone `wdvdone' `wdvk'
      local wdvnp = `wdvnp' + `wdvn'
    }
  }
  if "`wdvdone'" == "" exit
  * repost, not post: `ereturn repost V=' replaces the matrix without re-running the
  * screen that `ereturn post' applies (measured: repost accepts an indefinite matrix and
  * a -1e100 diagonal verbatim), so the sentinel rows the partial-identification rescue
  * wrote survive the splice exactly as they were posted
  cap ereturn repost V = `WV'
  if _rc != 0 exit
  ereturn hidden local wdid_vcvstate analytic-blocks
  ereturn hidden local wdid_vblocks `wdvdone'
  ereturn hidden scalar wdid_vblockn = `wdvnp'
end

cap prog drop wooldid_sgdeg
prog def wooldid_sgdeg, eclass
version 16.0
  * BS-15 (author ruling 2026-09-07): blank the subgroup pairwise standard errors and both
  * across-subgroup equality tests, with a named caution, when the analytic covariance of the
  * subgroup effects is numerically rank-deficient in the contrast direction. This is the
  * fourth site of the policy the bootstrap, boottstat and permtstat paths already carry; it
  * was the only one that computed and printed the residue in silence while jointtests on the
  * identical call named the same degeneracy.
  *
  * ANALYTIC PATH ONLY (see the gate at the call site). Under bootstrap() wooldid_bootjnt takes
  * `pw0' = e(pwsgs) and reads column 2 as boottstat's studentization denominator and as
  * wooldid_bootsgw's observed comparator; under permute() wooldid_permwire keeps column 2 and
  * renames it se_analytic. Blanking it here would delete THEIR inference, not wooldid's.
  *
  * DETECTION. The construction is the siblings' (wooldid_bootsgw's W_obs leg): rebuild the
  * level covariance of the subgroup effects by polarisation from the two standard-error
  * columns already in hand, C_ab = (V_a + V_b - V_ab)/2, and reject a numerically singular
  * correlation matrix. Per pairwise row that is the implied correlation
  *     rho = (V_a + V_b - V_ab) / (2*sqrt(V_a*V_b)) ,
  * and at two subgroups it is the only form carrying information: the chain-contrast
  * covariance is then 1x1 and its correlation matrix is [1], so the siblings' eigen guard
  * cannot fire there. The eigen form on the m x m level correlation matrix is equivalent at
  * m = 2 (min/max eigenvalue = (1-rho)/(1+rho)) and is the generalisation if a family-level
  * test is ever wanted at m > 2.
  * The TOLERANCE is deliberately NOT the siblings' 1e-8, which is a floating-point rank test
  * and is inert on this failure. Measured on test/data/dta_binary.dta, the fixture the bug
  * sweep used, with two disjoint sets of 25 treated units:
  *              main att row     pretreatment att row
  *   cluster(i)   rho .99940            rho .99715
  *   robust       rho .51141            rho .49779
  * and on the clustered run the contrast standard errors are 1/44 and 1/20 of their robust
  * counterparts - the same pathology on both rows. The cut is 1 - 1e-2, which catches both
  * and still sits half the range away from a healthy design. A cut at 1 - 1e-3, tried first,
  * caught the main row and left the pretreatment row reporting a contrast standard error
  * twenty times too small, which is the defect the ruling names. The cut is a judgment call
  * on a number the ruling did not fix; it lives in ONE local below so it can be moved.
  syntax , [si(string)]
  cap confirm matrix e(pwsgs)
  if _rc != 0 exit
  local sgl `e(wdid_sglist)'
  if wordcount("`sgl'") < 2 exit
  * the ruled tolerance: a row is flagged when 1 - rho < this
  local wdsgdtol 1e-2
  tempname pw
  matrix `pw' = e(pwsgs)
  * the per-subgroup level tables, read once. their rows are "main:att", "pretreatment:ame",
  * ... - the same (eq, item) keys wooldid_bootjnt's positional walk builds. ABSENT under
  * suppresssgprimaryeffects, where the level margins are never computed at all: every row is
  * then skipped and the screen is inert (documented in the help)
  foreach sg of local sgl {
    cap confirm matrix e(wdidsgresults_`sg')
    if _rc != 0 continue
    tempname wsg`sg'
    matrix `wsg`sg'' = e(wdidsgresults_`sg')
    local wsgn`sg' : rowfullnames `wsg`sg''
  }
  local wdsgdn 0
  local wdsgdrows
  local wdsgdfam
  local wdsgdmax 0
  local prn : rowfullnames `pw'
  local rr 0
  foreach rn of local prn {
    local ++rr
    * the eq/name split and the (eq, item) key, byte-for-byte the walk wooldid_bootjnt uses on
    * the same matrix; with three or more subgroups stata marks the base level in the row name
    gettoken eqn nmn : rn, parse(":")
    local nmn = subinstr("`nmn'", ":", "", 1)
    local ispre = strpos("`eqn'", "pre") > 0
    local item = subinstr("`eqn'", "pre", "", 1)
    local item = subinstr("`item'", "_", "", 1)
    if !`ispre' local item `eqn'
    local eqwant = cond(`ispre', "pretreatment", "main")
    if !regexm("`nmn'", "^([0-9]+)(bn?)?vs([0-9]+)(bn?)?([.].*)?$") continue
    local sga = regexs(1)
    local sgb = regexs(3)
    local sea .
    local seb .
    local kk 0
    foreach wn of local wsgn`sga' {
      local ++kk
      if "`wn'" == "`eqwant':`item'" local sea = `wsg`sga''[`kk', 2]
    }
    local kk 0
    foreach wn of local wsgn`sgb' {
      local ++kk
      if "`wn'" == "`eqwant':`item'" local seb = `wsg`sgb''[`kk', 2]
    }
    if `sea' >= . | `seb' >= . | `sea' <= 0 | `seb' <= 0 continue
    local sed = `pw'[`rr', 2]
    * a MISSING contrast standard error COUNTS AS rank-deficient, exactly as the sibling guard
    * does (wooldid_bootsgw returns on missing(sc)): it is the extreme of the same failure - the
    * boottstat disclosure names it "missing, or numerically zero with a p-value of 0" - and on
    * the binary path `sgpost' is still set there, so the chain test or its overjoint fallback
    * can still have posted an equality p that must be withheld. Blanking an already-blank row
    * is idempotent; withholding the family's equality test is the point. Fast mode and the
    * AMESE R1 scrub never reach this line: both blank columns 2-6 of the LEVEL matrices too, so
    * the level-standard-error check above has already skipped those rows. A contrast standard
    * error of exactly 0 gives rho >= 1 and is flagged.
    local rho = cond(`sed' >= ., 2, (`sea'^2 + `seb'^2 - `sed'^2)/(2*`sea'*`seb'))
    if `rho' <= 1 - `wdsgdtol' continue
    local ++wdsgdn
    local wdsgdrows `wdsgdrows' `rr'
    if `rho' > `wdsgdmax' local wdsgdmax = `rho'
    local wdsgds = cond(`ispre', "joint_pre_`item'_allsgssame", "joint_`item'_allsgssame")
    if strpos(" `wdsgdfam' ", " `wdsgds' ") == 0 local wdsgdfam `wdsgdfam' `wdsgds'
  }
  if `wdsgdn' == 0 exit
  * the disclosure state, read back by wooldid_display so replay reproduces the caution
  ereturn hidden scalar wdid_sgdegen = `wdsgdn'
  ereturn hidden local wdid_sgdegenr = strofreal(`wdsgdmax', "%6.4f")
  * safety(off): caution and report, do not blank. gates on the USER's safety() value
  * (safetyinput), never on `safety', which contreatwithin overwrites internally
  if "`si'" == "off" {
    ereturn hidden local wdid_sgdegenact report
    exit
  }
  ereturn hidden local wdid_sgdegenact blank
  foreach cc of local wdsgdrows {
    forval jj = 2/6 {
      matrix `pw'[`cc', `jj'] = .
    }
  }
  ereturn matrix pwsgs = `pw'
  * the equality tests are blanked to MISSING, not dropped - the convention the bootstrap and
  * permutation paths already use for a withheld allsgssame scalar
  foreach s of local wdsgdfam {
    if "`e(`s')'" != "" ereturn scalar `s' = .
  }
end

* THE PART MAP. This program's entire control flow is `part(string)', matched by string
* equality with NO else - so an unknown or mistyped part prints nothing and returns 0, the
* silent-option failure this codebase has been bitten by before. Five parts, and four
* callers take different subsets:
*   1     header, summary statistics, the main table, joint tests, adversarial
*   roth   REPLAY ONLY - on a live run the Roth pre-trend-power block is printed by
*          wooldid_ptp itself, called between part(1) and part(2) by the worker on a plain
*          run and by the shell on the bootstrap path
*   2     subgroup tables, pairwise contrasts, the gteffects block
*   perm   permutation footer, existence-gated on e(perm_method)
*   boot   bootstrap footer, existence-gated on e(boot_method)
* The shell's replay calls all five in order; the permute path calls 1, 2, perm (no roth -
* validate refuses pretrendpower under permute); the bootstrap path calls 1, wooldid_ptp,
* 2, boot; the worker calls 1 and 2 only, and only when neither engine is active. The two
* engine footers are mutually exclusive by their gates, so replaying the other kind prints
* nothing rather than erroring.
cap prog drop wooldid_display
prog def wooldid_display
version 16.0
  syntax , part(string)
  * replay contract gate: e(wdid_replayv) is the version of the display contract the
  * store was written under. bump the posted value (in wooldid_worker) whenever the
  * display contract changes; the exact match here refuses cross-version replays
  * instead of misrendering a store this code no longer understands
  if "`e(wdid_replayv)'" != "2" {
    dis as error "Error: the stored results were posted by a wooldid version without replay support or with a different replay contract (e(wdid_replayv)='`e(wdid_replayv)'', this version requires '2'); nothing can be redisplayed. Refit wooldid."
    error 301
  }
  * NFRPL-03: mirror wooldid_esplots' existence guard - on a corrupt/fabricated store the
  * matlist below would otherwise render a synthetic 1x1 "." table (an absent e() matrix
  * reads as a missing scalar). every genuine wooldid post carries the matrix, so this
  * cannot fire on a valid store. scope decision (deliberate, realistic-vs-pathological):
  * the central matrix only - the summarystats/wdidsgresults/pwsgs matlist sites remain
  * fabrication-exposed
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    dis as error "Error: the stored results no longer contain e(wdidmainresults); rerun wooldid."
    error 301
  }

  * rehydrate the display state from e(). each local below mirrors the worker local of
  * the same meaning byte-for-byte: option/label strings come from the hidden wdid_*
  * macros (string-exact by construction), data-derived values from posted e() scalars
  * (locals were themselves "= e(...)" evaluations upstream, so re-evaluating here
  * reproduces the same strings, including "." for a never-set value)
  local y `e(depvar)'
  local vcespecifier `"`e(wdid_vcespec)'"'
  local poisson `e(wdid_poisson)'
  local contreat `e(wdid_contreat)'
  local subgroups `e(wdid_sgvar)'
  local sglist `e(wdid_sglist)'
  * F2: the withheld-pre-treatment disclosure, posted by wooldid_worker. absent on every run
  * that withheld nothing, and on any store written before the macro existed
  local sgnopre `e(wdid_sgnopre)'
  local sgnoprewhy `e(wdid_sgnoprewhy)'
  local sgnoprerow `e(wdid_sgnoprerow)'
  local estimandlist `e(wdid_estimands)'
  local customweightmultiplier `e(wdid_cwm)'
  local semielasticity `e(wdid_semielas)'
  local contreatelasticitytype `e(wdid_cetype)'
  * BS-6 (2026-09-07): the two macros above hold the OPERATOR margins used. The headers must
  * name what the user asked for - which the remapped quantity now genuinely is - so prefer
  * the raw request where it was recorded. Absent on every run where no remap fired and on
  * pre-BS-6 stores, so replay of an older estimate is unchanged.
  if "`e(wdid_semielasraw)'" != "" local semielasticity `e(wdid_semielasraw)'
  if "`e(wdid_cetyperaw)'" != "" local contreatelasticitytype `e(wdid_cetyperaw)'
  local suppressjointtests `e(wdid_suppressjt)'
  local suppresssgprimaryeffects `e(wdid_suppresssgpe)'
  local adversarial `e(wdid_adversarial)'
  local fastnose `e(wdid_fastnose)'
  local summarystats `e(wdid_sumstats)'
  local espre `e(wdid_espre)'
  local Npre_trebase `e(wdid_nptb)'
  local Npre_rat `e(wdid_nprat)'
  local Wpre_rat `e(wdid_wprat)'
  local esfixedbaseperiod
  if "`e(estimate_type)'" == "fixed_reference_period" local esfixedbaseperiod esfixedbaseperiod
  local N = e(N)
  local Nclus = e(Nclus)
  local droppedobs = e(N_dropped)
  local r2 = e(r2)
  local r2adj = e(r2adj)
  local r2_within = e(r2_within)
  local r2_withinadj = e(r2_withinadj)
  local r2pseud = e(r2pseud)

  local poistext
  if "`poisson'" == "poisson" local poistext " Poisson Regression"
  local elasinj
  if "`semielasticity'" == "eydx" local elasinj Treatment Effects Presented as Semi-elasticities
  if "`contreatelasticitytype'" != "dydx" & "`elasinj'" != "" local elasinj `elasinj',
  if "`contreatelasticitytype'" == "eyex" local elasinj `elasinj' Marginal Treatment Effects Presented as Elasticities
  if "`contreatelasticitytype'" == "eydx" local elasinj `elasinj' Marginal Treatment Effects Presented as Semi-elasticities (eydx)
  if "`contreatelasticitytype'" == "dyex" local elasinj `elasinj' Marginal Treatment Effects Presented as Semi-elasticities (dyex)

  if "`part'" == "1" {
  dis as result  " "
  if "`fastnose'" == "" dis as result  "Wooldid Estimation for Outcome: `y'; Standard Errors: `vcespecifier'"
  else dis as result  "Wooldid Estimation for Outcome: `y'; Standard Errors: not computed (fast mode)"
  if "`fastnose'" != "" dis as text "SEs and related materials are not estimated in fast mode."
  dis as result  " "
  if "`Nclus'" != "." & "`Nclus'" != "" {
    dis as result  "N clusters = `Nclus'; N = `N' (`droppedobs' obs dropped from initial sample)"
  }
  else{
    dis as result  "N = `N' (`droppedobs' obs dropped from initial sample)"
  }
  if "`poisson'" == "" {
    if "`r2_within'" != "." {
      dis as result  "R2 = " %-6.4f `r2' ";  R2adj = " %-6.4f `r2adj' ";  R2-within = " %-6.4f `r2_within' ";  R2-withinadj = " %-6.4f `r2_withinadj'
    }
    else{
      dis as result  "R2 = " %-6.4f `r2' ";  R2adj = " %-6.4f `r2adj'
    }
  }
  else{
    dis as result  "pseudo-R2 = " %-6.4f `r2pseud'
  }

  dis as result  " "


  if "`summarystats'" != "" {

    if "`esfixedbaseperiod'" == "" &  `espre' > 0  {
      dis as result  " "
      dis as result  "Treatment group reference period size: N = `Npre_trebase'."
      dis as result  "`Npre_rat'% of pre-treatment treatment group observations are in the reference period."
      dis as result  "The reference period contains `Wpre_rat'% of the pre-treatment treatment group by sample weight."
      dis as result  " "
    }

    local contreatsummaryinjector
    if "`contreat'" != "" local contreatsummaryinjector " and Continuous Treatment Variable (c)"
    dis as result  "Summary Stats on Outcome Variable (y)`contreatsummaryinjector':"
    wooldid_mlfmt , mat(e(summarystats))
    matlist e(summarystats), cspec(`wdmlfmt') rspec(`wdmlfmtr')

    if "`poisson'" == "poisson" dis as result  "Note: summary stats above are for untransformed Y, while coefficients are from Y=exp(regression equation) and have a log-linear type interpretation."
  }


  ************* now show main results with matlist
  dis as result  " "
  dis as result  "Main`poistext' Results (Full Estimation Sample): `elasinj'"
  if "`customweightmultiplier'" != "" dis as result  "Note: all Estimates Shown Apply `customweightmultiplier' Weights Post-Estimation"
  wooldid_mlfmt , mat(e(wdidmainresults))
  matlist e(wdidmainresults), cspec(`wdmlfmt') rspec(`wdmlfmtr')
  wooldid_pfloor e(wdidmainresults)

  * BS-10 (ruling 2026-09-07): the IM reference distribution is row-specific - t(q-1) with q
  * the clusters contributing to that row - so when q is not common across the reported rows
  * no single posted df reproduces every printed p. e(df_r) carries the MINIMUM; this line
  * says so. Self-gating on a store written before BS-10: both scalars read missing and
  * ". > ." is false, which is why this read alone takes no wdid_replayv bump.
  if e(imdf_min) < . & e(imdf_max) > e(imdf_min) {
    dis as text "  Ibragimov-Muller degrees of freedom are not common across the rows above: they run from `=e(imdf_min)' to `=e(imdf_max)', one less than the number of clusters contributing to each row. e(df_r) carries the minimum, `=e(imdf_min)', so a test or lincom run by hand is conservative for every row rather than matching the widest one; each row's own p-value and interval in the table above use that row's own degrees of freedom. The stage-1 regression's residual degrees of freedom are kept as hidden e(df_r_analytic), and the e(gteffects) cells - stage-1 regression coefficients, not IM quantities - use that value."
  }

  * F8: composition table prints with the other descriptive material under summarystats
  cap confirm matrix e(escomp)
  if "`summarystats'" != "" & _rc == 0 {
    dis as result  " "
    dis as result  "Event-Study Composition by Horizon (Full Estimation Sample; also in e(escomp)):"
    matlist e(escomp), format(%12.4g)
    * DOC-03: the six coverage scalars were computed and returned but never displayed;
    * cap-confirm gate (they post under cap, so benign absence is possible - skip
    * silently, matching the surrounding style). neutral pool wording: "reference
    * period" would be accurate in at most one of the four reference regimes
    cap confirm scalar e(escomp_refN)
    if _rc == 0 {
      dis as text "  Reference-normalization pool: N = " %-9.0g e(escomp_refN) " ; treated units = " %-6.0g e(escomp_refnunits) " ; onset cohorts = " %-6.0g e(escomp_refncoh) " ; att-style weight = " %-9.0g e(escomp_refW)
      dis as text "  Entire treated estimation sample: N = " %-9.0g e(escomp_trtN) " ; att-style weight = " %-9.0g e(escomp_trtW)
    }
  }

  * show the joint  tests
  if "`suppressjointtests'" == "" {
    dis as result  " "
    if "`e(joint_horiz)'" != "" {
      dis as result  "p-value for joint test that all pre-treatment att horizon estimates are 0: " %-6.4f e(joint_horiz)
    }
    if "`e(joint_horiz_i)'" != "" {
      dis as result  "p-value for joint test that all pre-treatment att_it horizon estimates are 0: " %-6.4f e(joint_horiz_i)
    }
    if "`e(joint_Choriz)'" != "" {
      dis as result  "p-value for joint test that all pre-treatment ame horizon estimates are 0: " %-6.4f e(joint_Choriz)
    }
    if "`e(joint_Choriz_i)'" != "" {
      dis as result  "p-value for joint test that all pre-treatment ame_it horizon estimates are 0: " %-6.4f e(joint_Choriz_i)
    }
  }
  dis as result  " "

 if "`adversarial'" != "" {
   foreach typetext in att att_it ame ame_it {
     if e(adversarial_pval_`typetext') < . {
       if "`e(perm_method)'" == "" dis as result "p-value for joint test of null that the `typetext' coefficients and the reference period are on a single line: " %-6.4f e(adversarial_pval_`typetext')
       else dis as result "permutation p-value (sharp null of no effect at any horizon) of the adversarial line-fit statistic on the `typetext' coefficients: " %-6.4f e(adversarial_pval_`typetext')
     }
   }
 }
  }

  if "`part'" == "roth" {
  * replay-only: redisplay the Roth (2022) pre-trend power block from e(roth_*) without
  * re-running the quadrature (the NIS computation can cost minutes). cautions that are
  * reconstructible from posted e() ARE replayed (rank-deficiency from roth_keff vs
  * roth_npre; interrupted-pass and completed-but-unavailable states from roth_ok);
  * per-type estimation-time skip cautions remain live-only. on the live path the worker
  * calls wooldid_ptp itself between part(1) and part(2), so this part is never run live.
  * NOTE (deliberate replay-contract exception): the roth_ok reads below consult an e()
  * member without a wdid_replayv bump - the only affected window is same-day pre-marker
  * dev stores, and the degradation is a conservative caution + skip, versus a bump
  * refusing every older store
  * NFPTP-07: e(roth_ok) is written after a completed diagnostics pass; the worker runs
  * wooldid_ptp under cap noi with individually-cap'd ereturns, so an abort partway
  * through the type loop leaves earlier types' roth_* in e() with no marker. formerly
  * nothing read the marker, so replay presented an interrupted pass as complete
  local rothany 0
  foreach ptpty in att att_it ame ame_it {
    cap confirm matrix e(roth_power_`ptpty')
    if _rc == 0 local rothany 1
  }
  if `rothany' & e(roth_ok) != 1 {
    dis as text "Caution: the stored pretrendpower() diagnostics are incomplete (the diagnostics pass was interrupted); partial results are not replayed. Refit wooldid to recompute them."
  }
  else if `rothany' == 0 & e(roth_ok) == 1 {
    dis as text "Caution: pretrendpower() was requested, but no event-study coefficient vector with a valid variance-covariance matrix was available (e.g. under timetrends); no diagnostics were computed."
  }
  else if `rothany' & e(roth_ok) == 1 {
  local dfr = e(df_r)
  local ptpdfok 0
  if "`dfr'" != "" & "`dfr'" != "." {
    if real("`dfr'") > 0 & real("`dfr'") < . local ptpdfok 1
  }
  local ptpref `e(wdid_esrel)'
  local ptpheader 0
  foreach ptpty in att att_it ame ame_it {
    cap confirm matrix e(roth_power_`ptpty')
    if _rc != 0 continue
    local ptpk = e(roth_npre_`ptpty')
    local ptpkeff = e(roth_keff_`ptpty')
    local ptpslope = e(roth_slope_`ptpty')
    local ptpse = e(roth_slope_se_`ptpty')
    local ptpll = e(roth_slope_ll_`ptpty')
    local ptplvl = cond(e(level) < ., strofreal(e(level)), "95")
    local ptpul = e(roth_slope_ul_`ptpty')
    local ptpprep = e(roth_prep_`ptpty')
    local ptphbar = e(roth_hbar_`ptpty')
    local ptpW = e(roth_W_`ptpty')
    local ptpnisok = ("`e(roth_nissize_`ptpty')'" != "")
    local ptpnissize = e(roth_nissize_`ptpty')
    * NFPTP-02/07: the rank-deficiency disclosure is reconstructible from e() and is
    * now replayed (byte-identical text and position: live prints it pre-header)
    if `ptpkeff' < `ptpk' dis as text "Caution: the `ptpty' pre-period variance-covariance block is rank-deficient (rank `ptpkeff' of `ptpk'). The slope, its CI, the joint test (now `ptpkeff' degrees of freedom), and the detectable-trend columns then reflect only the linearly independent subset of pre-period coefficients a generalized inverse retains - NOT a summary of all pre-period evidence - and can differ arbitrarily from statistics built on the discarded coefficients; interpret with care."
    if `ptpheader' == 0 {
      dis as result " "
      dis as result "Pre-Trend Power Diagnostics (Roth 2022):"
      local ptpheader 1
    }
    dis as result "  `ptpty': GLS linear pre-trend fit on `ptpk' pre-treatment coefficient(s); reference relative period = `ptpref'"
    dis as result "    slope = " %-9.0g `ptpslope' " (se " %-9.0g `ptpse' "; `ptplvl'% CI " %-9.0g `ptpll' " to " %-9.0g `ptpul' ")"
    if `ptpdfok' dis as result "    joint pre-trend test: F(`ptpkeff', `dfr') = " %-8.0g `=`ptpW'/`ptpkeff'' ", p = " %-6.4f `ptpprep'
    else dis as result "    joint pre-trend test: chi2(`ptpkeff') = " %-8.0g `ptpW' ", p = " %-6.4f `ptpprep'
    if `ptpnisok' dis as result "    size of the any-pre-period-individually-significant pre-test: " %-6.4f real("`ptpnissize'")
    if !`ptpnisok' & `ptpk' > 4 dis as text "    note: mdt_nis and the NIS-test size are omitted (`ptpk' pre-periods exceed the exact-quadrature budget of 4); mdt_wald is unaffected. For the NIS convention at larger k, feed e(b_es_`ptpty')/e(V_es_`ptpty')/e(numpre_`ptpty') to Roth's pretrends package."
    dis as result "    weighted mean post-treatment relative period (hbar) = " %-9.0g `ptphbar'
    dis as result "    minimum detectable linear trend, and the bias such a trend would imply for the"
    dis as result "    `ptpty' pooled estimand if it continued into the post-treatment period:"
    * NFPTP-03: the att_it/ame_it block is the shared event-study curve; when the user
    * requested att_i/att_itime the block stands in for them but its hbar/bias columns
    * use att_it-family post-period weights - say so rather than implying the requested
    * pooled estimand's own weighting. exact-token match (att_i is a prefix of att_it
    * and att_itime; cf. the cfxplottypes exact-token precedent)
    if inlist("`ptpty'", "att_it", "ame_it") {
      local ptpsurr
      local ptpests `e(wdid_estimands)'
      if `: list posof "att_i" in ptpests' > 0 local ptpsurr att_i
      if `: list posof "att_itime" in ptpests' > 0 local ptpsurr `ptpsurr' att_itime
      if "`ptpsurr'" != "" {
        if "`ptpty'" == "ame_it" local ptpsurr = subinstr("`ptpsurr'", "att", "ame", .)
        dis as text "    note: this `ptpty'-weighted diagnostic stands in for the requested `ptpsurr' estimand(s); its hbar and bias columns use `ptpty' post-period weights."
      }
    }
    matlist e(roth_power_`ptpty'), format(%12.4g)
    dis as text "    note: mdt_wald powers the joint pre-trend test above; mdt_nis powers the test that any"
    dis as text "    pre-period coefficient is individually significant (the convention of Roth's pretrends package)."
    dis as text "    Both pre-tests are at the 5% level, whatever level() is in force; level() moves the slope CI only."
    if "`poisson'" == "poisson" & "`e(wdid_poisexp)'" == "" dis as text "    note: poisson event-study coefficients are on the linear-index (xb) scale, and so are these diagnostics."
  }
  }
  }

  if "`part'" == "2" {
  dis as result  " "

  * now show the sg results
  if "`subgroups'" != "" {
    if "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        dis as result  " "
        local poistext
        if "`poisson'" == "poisson" local poistext "Poisson Regression "
        dis as result  "`poistext'Results for Subgroup `subgroups' == `sgcase': `elasinj'"
        if "`customweightmultiplier'" != "" dis as result  "Note: all Estimates Shown Apply `customweightmultiplier' Weights Post-Estimation"
        wooldid_mlfmt , mat(e(wdidsgresults_`sgcase'))
        matlist e(wdidsgresults_`sgcase'), cspec(`wdmlfmt') rspec(`wdmlfmtr')
        wooldid_pfloor e(wdidsgresults_`sgcase')
        * F2 (2026-09-14): name the missing pre-treatment row. posted by wooldid_worker and
        * read from e() so replay reproduces it, as the wdid_sgdegen caution is; written in
        * the register of the SI-7 horizon caution, whose sibling it is - same mechanism on
        * the whole-pre-history row rather than on a single relative period
        local sgnopi : list posof "`sgcase'" in sgnopre
        if `sgnopi' > 0 {
          local sgnopr : word `sgnopi' of `sgnoprewhy'
          local sgnopq : word `sgnopi' of `sgnoprerow'
          local sgnopref The subgroups withheld this way are named in e(wdid_sgnopre), with the reason for each in e(wdid_sgnoprewhy) and the rows withheld in e(wdid_sgnoprerow).
          if "`sgnopr'" == "nocells" dis as text "Caution: subgroup `subgroups' == `sgcase' has no usable pre-treatment observations outside the reference period, so it has no pre-treatment average and no pre-treatment row appears in the table above; its remaining rows, and all other subgroups, are unaffected. `sgnopref'"
          else if "`sgnopq'" == "ame" dis as text "Caution: the pre-treatment average marginal effect of `contreat' for subgroup `subgroups' == `sgcase' is not estimable (margins reports its level as omitted), so its pre-treatment ame row is absent from the table above; the subgroup's pre-treatment att row, its remaining rows, and all other subgroups, are unaffected. This typically arises when the specified fixed effects absorb the subgroup-by-time variation in the pre-period. `sgnopref'"
          else dis as text "Caution: the pre-treatment average(s) for subgroup `subgroups' == `sgcase' are not estimable (margins reports its level as omitted), so no pre-treatment row appears in the table above; its remaining rows, and all other subgroups, are unaffected. This typically arises when the specified fixed effects absorb the subgroup-by-time variation in the pre-period. `sgnopref'"
        }
        dis as result  " "
        * F8: per-subgroup composition table under summarystats
        cap confirm matrix e(escomp_s`sgcase')
        if "`summarystats'" != "" & _rc == 0 {
          dis as result  "Event-Study Composition by Horizon, Subgroup `subgroups' == `sgcase' (also in e(escomp_s`sgcase')):"
          matlist e(escomp_s`sgcase'), format(%12.4g)
          dis as result  " "
        }
        * show the joint  tests
        if "`suppressjointtests'" == "" {

          if "`e(joint_horiz_s`sgcase')'" != "" {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment att horizon estimates are 0: " %-6.4f `e(joint_horiz_s`sgcase')'
          }
          if "`e(joint_horiz_i_s`sgcase')'" != ""  {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment att_it horizon estimates are 0: " %-6.4f `e(joint_horiz_i_s`sgcase')'
          }


          if "`e(joint_Choriz_s`sgcase')'" != "" {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment ame horizon estimates are 0: " %-6.4f `e(joint_Choriz_s`sgcase')'
          }

          if "`e(joint_Choriz_i_s`sgcase')'" != ""  {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment ame_it horizon estimates are 0: " %-6.4f `e(joint_Choriz_i_s`sgcase')'
          }

          dis as result  " "
        }
        if "`adversarial'" != "" {
          foreach typetext in att att_it ame ame_it {
            if e(adversarial_pval_`typetext'_s`sgcase') < . {
              if "`e(perm_method)'" == "" dis as result "p-value for joint test of null that the `typetext' coefficients for sg=`sgcase' and the reference period are on a single line: " %-6.4f e(adversarial_pval_`typetext'_s`sgcase')
              else dis as result "permutation p-value (sharp null of no effect at any horizon) of the adversarial line-fit statistic on the `typetext' coefficients for sg=`sgcase': " %-6.4f e(adversarial_pval_`typetext'_s`sgcase')
            }
          }
        }
        dis as result " "
      }
    }
    * having shown the main sg materials, show the pairwise comparisons
    if wordcount("`sglist'") > 1 {
      dis as result  "Pairwise Comparisons of Different Subgroups' Estimates: `elasinj'"
      wooldid_mlfmt , mat(e(pwsgs))
      matlist e(pwsgs), cspec(`wdmlfmt') rspec(`wdmlfmtr')
      wooldid_pfloor e(pwsgs)
      * BS-15 (author ruling 2026-09-07): the realized-degeneracy disclosure, posted by
      * wooldid_sgdeg on the analytic path and read from e() so replay reproduces it. Written as
      * the sibling of the jointtests caution: same mechanism, same remedies
      if "`e(wdid_sgdegenact)'" != "" {
        if "`e(wdid_sgdegenact)'" == "blank" dis as text "Caution: the analytic covariance of the subgroup effects is numerically rank-deficient in the contrast direction for `e(wdid_sgdegen)' pairwise row(s) (largest implied correlation between two subgroups' estimates `e(wdid_sgdegenr)'). The saturated cell dummies leave each cell's residuals with no within-cluster variation to identify, so the contrast's standard error is cancellation residue rather than an estimated variance - the same degeneracy named in the jointtests caution. The standard error, statistic, p-value and interval of those rows, and the across-subgroup equality test(s) covering them, are reported as absent rather than computed from that covariance; the point estimates are unaffected. wooldid's randomization inference (permute()) or bootstrap-based clustered SEs give subgroup contrasts and an equality test that are valid under this design - use bootstrap() WITHOUT boottstat, which studentizes each contrast by exactly this standard error - as does clustering strictly finer than the cohort-by-period cells."
        else dis as text "Caution: the analytic covariance of the subgroup effects is numerically rank-deficient in the contrast direction for `e(wdid_sgdegen)' pairwise row(s) (largest implied correlation between two subgroups' estimates `e(wdid_sgdegenr)'), so those contrast standard errors are cancellation residue rather than estimated variances. safety(off) was specified, so they and the across-subgroup equality tests are reported rather than withheld; read them as the residue they are. permute(), bootstrap() without boottstat, or clustering strictly finer than the cohort-by-period cells gives inference that is valid under this design."
      }
      dis as result " "
      foreach item in `estimandlist' {
        if "`e(joint_`item'_allsgssame)'" != "" dis as result  "p-value for joint test that all SGs have the same `item': " %-6.4f e(joint_`item'_allsgssame)
        if "`e(joint_pre_`item'_allsgssame)'" != "" dis as result  "p-value for joint test that all SGs have the same pre-treatment `item': " %-6.4f e(joint_pre_`item'_allsgssame)
      }
      dis as result  " "
      if "`contreat'" != "" {
        foreach item in `estimandlist' {
          local itemtext = subinstr("`item'","att","ame",.)
          if "`e(joint_`itemtext'_allsgssame)'" != "" dis as result  "p-value for joint test that all SGs have the same `itemtext': " %-6.4f e(joint_`itemtext'_allsgssame)
          if "`e(joint_pre_`itemtext'_allsgssame)'" != "" dis as result  "p-value for joint test that all SGs have the same pre-treatment `itemtext': " %-6.4f e(joint_pre_`itemtext'_allsgssame)
        }
        dis as result  " "
      }
    }
  }

  * F3: cohort-by-relative-period building blocks - existence-gated on e(gteffects),
  * so replay reproduces it and older stores simply print nothing
  cap confirm matrix e(gteffects)
  if _rc == 0 {
    tempname gtd
    matrix `gtd' = e(gteffects)
    local gtrows = rowsof(`gtd')
    dis as result " "
    dis as result "Cohort-by-Relative-Period Effects (regression scale; full table incl. weight totals in e(gteffects)):"
    if "`poisson'" == "poisson" dis as text "Note: under poisson these are on the linear-index (xb) scale; poisexpresults/semielasticity adjustments are not applied."
    if "`e(perm_gt_blank)'" == "1" dis as text "Note: under permute() these cells have no permutation distribution (their keys are cohort onsets, which the replications reassign), so their inference columns are blank; rows outside the reported event-study window have inwindow = 0."
    else if "`e(boot_method)'" != "" dis as text "Note: under bootstrap() the se, p, lb and ub columns are the `=cond("`e(boot_gteffects_family)'" == "", "percentile", "`e(boot_gteffects_family)'")' bootstrap family of each cell's own replications (every family is stored in e(boot_*_gteffects)); rows outside the reported event-study window have inwindow = 0."
    else dis as text "Note: standard errors are the stage-1 regression's coefficient-level VCE; rows outside the reported event-study window have inwindow = 0."
    if `gtrows' <= 200 {
      tempname gtds
      * e(gteffects) is 15 columns (the stripe is set in Mata: i, onset, caltime, relyr,
      * estimate, se, p, lb_ci, ub_ci, inwindow, nobs, w_att, w_att_it, w_att_i,
      * w_att_itime). 9 is the boundary between the printed block and the export-only
      * columns, not a count of anything, and it must move if that stripe does - the slice
      * is positional, so a restripe prints the wrong columns silently. NOTE the notes above
      * point at inwindow, which is column 10 and therefore NOT in this printed slice: they
      * are directing the reader to e(gteffects). The 200-row cap keeps a cohort-by-period
      * table on high-cardinality i from burying the rest of the output
      matrix `gtds' = `gtd'[1..., 1..9]
      matlist `gtds', format(%9.0g)
    }
    else dis as text "  (`gtrows' rows; table not printed - see e(gteffects))"
  }
  }
  * F6: permutation-inference footer - existence-gated on e(perm_method), so replay
  * reproduces it and older stores print nothing. Rewritten 2026-09-05 (chunk 1): the
  * tables above now carry the permutation inference, so this block is the provenance
  * statement plus the pooled-estimand summaries
  if "`part'" == "perm" {
    if "`e(perm_method)'" != "" {
      dis as result "Randomization Inference (permute(`e(perm_method)')): `e(perm_completed)' of `e(perm_reps)' replications completed; `e(perm_npossible_str)' possible assignments."
      local wdplv = cond("`e(perm_level)'" == "", "95", "`e(perm_level)'")
      * BS-23 (author ruling 2026-09-07): e(perm_exact) is posted as (M < .), i.e. 1 on ANY
      * permexact run, complete or not - so the phrase names the set the shares are actually
      * taken over. The membership claim holds either way. (e(perm_failed) can be missing in
      * an old store and . > 0 is TRUE in Stata, hence the explicit < . test.)
      local wdpexs = cond((e(perm_failed) < . & e(perm_failed) > 0) | "`e(perm_incomplete)'" != "", "completed", "enumerated")
      local wdpex = cond("`e(perm_exact)'" == "1", "exact shares over the `wdpexs' assignments (the observed assignment is a member of its own reference set)", "the add-one Monte-Carlo form (1 + count)/(replications + 1)")
      if "`e(perm_pconv)'" == "absolute" {
        local wdppc two-sided absolute-value permutation p-values (permabsp): the share of replications whose absolute statistic is at or above the observed absolute statistic, `wdpex'
        local wdpiv symmetric shift intervals dual to that p-value - the estimate plus or minus the `wdplv'th-percentile absolute
      }
      else if "`e(perm_pconv)'" == "raw" {
        local wdppc two-sided absolute-value permutation p-values WITHOUT the add-one adjustment (permrawp): c/n, the share of replications whose absolute statistic is at or above the observed absolute statistic - the number ritest and Stata's permute prefix report; a p of exactly 0 means no replication reached the observed absolute statistic (the note under the table gives the resolution)
        local wdpiv symmetric shift intervals - the estimate plus or minus the `wdplv'th-percentile absolute
      }
      else if "`e(perm_pconv)'" == "left" | "`e(perm_pconv)'" == "right" {
        local wdppc one-sided `e(perm_pconv)'-tail permutation p-values (permside(`e(perm_pconv)')): the share of replications whose statistic is at or `=cond("`e(perm_pconv)'" == "left", "below", "above")' the observed statistic, `wdpex'; the two-sided equal-tailed p-values are 2 min(e(perm_p_lo), e(perm_p_hi))
        local wdpiv two-sided shift intervals at level() - the estimate minus the equal-tailed quantiles of the
      }
      else {
        local wdppc two-sided equal-tailed permutation p-values, `wdpex'
        local wdpiv shift intervals - the estimate minus the equal-tailed quantiles of the
      }
      if "`e(perm_stat)'" == "t" dis as result "  Every p-value and interval in the tables above is permutation-based. The permuted statistic is the t-ratio (permtstat): each replication's estimate over its own analytic standard error, ranked against the observed t-ratio t_perm = estimate/se_analytic. p-values are `wdppc'; the `wdplv'% intervals are `wdpiv' permuted t-ratio`=cond("`e(perm_pconv)'" == "absolute", "", "s")', times se_analytic."
      else dis as result "  Every p-value and interval in the tables above is permutation-based: the observed estimate is ranked among the replications' estimates under reassigned treatment onsets. p-values are `wdppc'; the `wdplv'% intervals are `wdpiv' permut`=cond("`e(perm_pconv)'" == "absolute", "ed estimate", "ation distribution")'."
      dis as result "  A shift interval inverts the permutation test under the fixed-zero-null pivotality approximation (it takes the estimate's sampling error to have the null permutation distribution); it is not a sharp-null randomization-test inversion. Its width grows with the true effect (the relabelled outcomes carry it), so it is calibrated only near the null; permtstat removes that growth for the post-treatment rows but not for the pre-treatment rows or the contrasts. Its implied standard errors are in e(perm_se_shift); the per-coefficient p-values, intervals, null-distribution SDs and replication counts are in e(perm_p), e(perm_ci), e(perm_sd) and e(perm_n), and the draws in e(permdist_b)."
      if "`e(perm_stat)'" == "t" dis as result "  The standard-error column is the ANALYTIC standard error (se_analytic), the t-ratio's own denominator; e(V) and e(V_es_*) hold the analytic covariance, so test or lincom run by hand on them are analytic, not permutation, quantities."
      else if "`e(perm_analyticse)'" == "1" dis as result "  The standard-error column is the ANALYTIC standard error (se_analytic, permanalyticse), shown for reference only; the t_perm column is blank because the permuted statistic is the estimate itself. e(V) and e(V_es_*) hold the analytic covariance, so test or lincom run by hand on them are analytic, not permutation, quantities."
      else dis as result "  No analytic standard errors were computed (the permute() fast-mode default; permanalyticse restores them), so the se_analytic and t_perm columns are blank and no e(V) or e(V_es_*) is posted."
      if "`e(perm_exact)'" == "1" dis as result "  The assignments are enumerated rather than sampled, so the p-values are exact shares with no Monte-Carlo error; e(perm_p_se) and e(perm_p_ci) are left missing (an element named below as taking the add-one form carries both). The one-sided shares are in e(perm_p_lo) and e(perm_p_hi)."
      else dis as result "  Monte-Carlo precision of the p-values: e(perm_p_se) holds the standard error of each p-value's underlying share and e(perm_p_ci) its `wdplv'% Clopper-Pearson interval (doubled under the equal-tailed rule); the one-sided shares are in e(perm_p_lo) and e(perm_p_hi)."
      local wdpqf 0
      cap confirm matrix e(perm_jointW)
      if _rc == 0 local wdpqf 1
      cap confirm matrix e(perm_sgW)
      if _rc == 0 local wdpqf 1
      cap confirm matrix e(perm_advW)
      if _rc == 0 local wdpqf 1
      if `wdpqf' dis as result "  The joint pre-trend tests, across-subgroup equality tests and adversarial line tests shown are one-sided permutation tests of a quadratic-form statistic in the permutation-covariance metric (the covariance of the replications' coefficients), so they need no per-replication covariance; the adversarial statistic is wooldid's line-fit statistic on that metric, and its p-value tests the sharp null of no effect at any horizon rather than the composite on-a-line null of the analytic test. Replication statistics: e(perm_jointW), e(perm_sgW), e(perm_advW)."
      * F065: multiplicity-adjusted p-values across the horizons of a family. PS-45 (author
      * ruling 2026-09-07): opt-in under permstepdown, and the pre-treatment horizons are a
      * separate family reported apart from the post-treatment one
      foreach wdpsdn in perm_stepdown perm_stepdown_i perm_Cstepdown perm_Cstepdown_i perm_stepdown_pre perm_stepdown_i_pre perm_Cstepdown_pre perm_Cstepdown_i_pre {
        cap confirm matrix e(`wdpsdn')
        if _rc != 0 continue
        local wdpsdb = subinstr("`wdpsdn'", "_pre", "", 1)
        local wdpsdw = cond(strpos("`wdpsdn'", "_pre") > 0, "pre-treatment", "post-treatment")
        local wdpsdf = cond("`wdpsdb'" == "perm_stepdown", "att", cond("`wdpsdb'" == "perm_stepdown_i", "att_it", cond("`wdpsdb'" == "perm_Cstepdown", "ame", "ame_it")))
        dis as result "  Multiplicity across the `wdpsdw' `wdpsdf' horizons (max-|studentized| step-down, Romano-Wolf; the studentization is the SD of each horizon's null draws): horizon: p_raw -> step-down p, where p_raw is the two-sided |t| share (add-one) in the step-down's own metric, not the table's convention - compare the two columns with each other rather than with the table`=cond("`wdpsdw'" == "pre-treatment", ". The pre-treatment horizons are their own family, adjusted among themselves and not pooled with the post-treatment ones; read them with the pre-treatment caveat below", "")'"
        tempname wdpsdm
        matrix `wdpsdm' = e(`wdpsdn')
        local wdpsdc : colnames `wdpsdm'
        local wdpsdk 0
        local wdpsdl
        foreach wdpsdh of local wdpsdc {
          local ++wdpsdk
          * e(perm_stepdown*) is 4 x (one column per horizon) with rows in the order
          * wooldid_permstep sets: 1 t_obs, 2 p_raw (the two-sided |t| add-one share in the
          * step-down's OWN studentized metric, not the table's p), 3 p_singlestep, 4
          * p_stepdown. This prints rows 2 and 4; row 3 is computed and stored but
          * deliberately not displayed. The subscripts are positional - reorder or insert a
          * row in that Mata stripe and this prints p_singlestep under the "step-down p"
          * label, a materially larger number, with no error and no shape change
          local wdpsdl `wdpsdl' `=subinstr("`wdpsdh'", "ES_`wdpsdf':", "", 1)': `=string(`wdpsdm'[2, `wdpsdk'], "%5.3f")' -> `=string(`wdpsdm'[4, `wdpsdk'], "%5.3f")';
        }
        dis as result "    `wdpsdl'"
      }
      if wordcount("`e(wdid_sglist)'") > 1 dis as result "  Pairwise subgroup contrasts are permuted as coefficient differences (in both modes); their kernel is e(perm_pwsgs_kernel)."
      * direction corrected 2026-09-05 (review): the reassigned onsets place actual
  * post-treatment outcomes into the draws' pre-treatment horizons, so the null
  * distribution WIDENS with the true effect while the observed pre-period statistics
  * do not (measured: null SD x7, p .49 -> .69 on a staggered fixture at effect 30; no
  * movement at all with a single onset date). The tests lose power against pre-trends
  if "`e(perm_precaveat)'" == "1" dis as result "  Pre-treatment rows (pretreatment estimates, pre-onset event-study horizons, their contrasts, and the joint pre-trend and adversarial statistics) are tested against the sharp null of no effect at ANY horizon. In staggered designs the reassigned onsets move actual post-treatment outcomes into the draws' pre-treatment horizons, so the null distribution widens with the true post-treatment effect while the observed pre-period statistics do not: their permutation p-values rise toward 1 as the true effect grows, and a large pre-period p-value is not by itself evidence of no pre-trend (with a single onset date the reference distribution is unaffected)."
      if "`e(perm_gt_blank)'" == "1" dis as result "  The gteffects cells are keyed by cohort onset, which the replications reassign, so they have no permutation distribution: their inference columns are blank."
      if "`e(perm_cfx_blank)'" == "1" dis as result "  The continuous-treatment plot cells were not permuted (assignment-dependent lattice bins); their p-value and interval columns are blank, and the standard error column is analytic."
      if "`e(perm_incomplete)'" != "" dis as result "  Elements not estimable in every enumerated assignment (exact shares over a reduced reference set; their n is in e(perm_n)): `e(perm_incomplete)'"
      if "`e(perm_addone)'" != "" dis as result "  Elements whose reference set does not contain the observed statistic (add-one p-values rather than exact shares): `e(perm_addone)'"
      if "`e(perm_offcentre)'" != "" dis as result "  Shift intervals left blank because the permutation distribution is not centred at zero (the p-values stand): `e(perm_offcentre)'"
      dis as result "  Pooled estimand summaries (p, `wdplv'% null band and SD of the permutation distribution`=cond("`e(perm_stat)'" == "t", ", on the t-ratio scale", "")'):"
      local wdptcihedge 0
      foreach wdpe in `e(wdid_permest)' {
        dis as result "    `wdpe': estimate " %10.0g e(b)[1, colnumb(e(b), "main_`wdpe'")] "   p " %6.4f e(perm_p_`wdpe') "   null band [" %10.0g e(perm_lb_`wdpe') ", " %10.0g e(perm_ub_`wdpe') "]   SD " %10.0g e(perm_sd_`wdpe')
        if e(perm_tlb_`wdpe') < . {
          dis as result "        inverted-|t| `wdplv'% CI (permtci) [" %10.0g e(perm_tlb_`wdpe') ", " %10.0g e(perm_tub_`wdpe') "]"
          local wdptcihedge 1
        }
      }
      if `wdptcihedge' dis as text "The inverted-|t| CI is symmetric about the estimate (the `wdplv'th percentile of the permuted absolute t-ratios times se_analytic) and rests on the same pivotality approximation as the shift intervals; it is not a sharp-null randomization-test inversion`=cond("`e(perm_pconv)'" == "absolute", "; under permabsp it coincides with the interval the table carries for the pooled estimands", "")'."
      dis as result " "
    }
  }


 if "`part'" == "boot" {
   * F-boot summary block (existence-gated: replay of non-bootstrap stores prints
   * nothing). all content from e() only
   if "`e(boot_method)'" == "" exit
   * keyed on e(wdid_vcespec), which names the resampling unit (e(clustvar), added by
   * PS-1, names the analytic cluster variable and not this); both modes take the
   * observation-level test (bug sweep BS-31)
   local wdbobs = strpos("`e(wdid_vcespec)'", "observation-level") > 0
   local wdbmt = cond("`e(boot_method)'" == "bayes", cond(`wdbobs', "Bayesian bootstrap (Dirichlet reweighting over observations, alpha = `e(boot_alpha)')", "Bayesian cluster bootstrap (Dirichlet reweighting, alpha = `e(boot_alpha)')"), cond(`wdbobs', "observation-level bootstrap (resampling with replacement)", "cluster bootstrap (resampling with replacement)"))
   dis as result "Bootstrap Inference (`wdbmt'):"
   dis as result "  " e(boot_completed) " of " e(boot_reps) " replications completed" _continue
   if e(boot_failed) > 0 dis as result " (" e(boot_failed) " failed and excluded)" _continue
   * the reference distribution is resolved from e(boot_df), NOT from e(boot_stat):
   * under boottstat the statistic label is "t_pt", which is neither "t" nor "z", so the
   * old equality test silently took the normal branch and printed the wrong reference
   local wdbrefd = cond(e(boot_df) < ., "t(" + strofreal(e(boot_df)) + ")", "standard normal")
   local wdbdsp = cond("`e(boot_display)'" == "", "percentile", "`e(boot_display)'")
   if "`e(boot_tstat)'" == "1" {
     dis as result "; the reported " e(boot_level) "% CIs and p-values come from the percentile ranks of the STUDENTIZED bootstrap distribution, and the standard error column is the ANALYTIC standard error (labelled se_analytic), which is the denominator the studentization uses. A studentized distribution has no standard deviation in coefficient units, so no bootstrap standard error is reported here."
   }
   else if "`wdbdsp'" == "wald" {
     dis as result "; reported standard errors are the standard deviation of the bootstrap distribution, and the `e(boot_stat)' statistics, p-values and " e(boot_level) "% CIs are referred to the `wdbrefd' distribution."
   }
   else {
     dis as result "; reported standard errors are the standard deviation of the bootstrap distribution, and the p-values and " e(boot_level) "% CIs come from the " cond("`wdbdsp'" == "bc", "BIAS-CORRECTED percentile", cond("`wdbdsp'" == "bca", "BIAS-CORRECTED AND ACCELERATED (BCa) percentile", "PERCENTILE")) " ranks of that distribution."
     if "`wdbdsp'" == "bca" {
       local wdbacc0 0
       cap confirm matrix e(boot_accel)
       if _rc == 0 mata: st_local("wdbacc0", strofreal(sum(st_matrix("e(boot_accel)") :< .) > 0))
       if `wdbacc0' dis as result "  The BCa acceleration comes from a delete-one-cluster jackknife (" e(boot_jk_completed) " fits; e(boot_accel)); its intervals are in e(boot_ci_bca) and p-values in e(boot_p_bca). The gteffects/cfx cells carry the bias-corrected family."
       else dis as result "  No BCa acceleration could be computed (" e(boot_jk_completed) " jackknife fits completed; the acceleration needs at least 3 fits with dispersion), so the BCa p-values and intervals are blank; e(boot_pctl_mainresults) and e(boot_bc_mainresults) carry the percentile and bias-corrected families."
     }
     dis as result "  A percentile interval is not symmetric about the estimate, so the estimate plus or minus a critical value times the reported standard error will not in general reproduce the reported interval - this is the ordinary reading of a bootstrap table, not a wooldid peculiarity. Use bootdisplay(wald) for the normal/t package in which they do agree."
   }
   dis as result "  Every family is computed whichever is displayed: intervals in e(boot_ci_normal)/e(boot_ci_percentile)/e(boot_ci_bc), p-values in e(boot_p_percentile)/e(boot_p_bc), and the complete tables in the hidden e(boot_wald_mainresults)/e(boot_pctl_mainresults)/e(boot_bc_mainresults)."
   * F033 (author ruling 2026-09-07): the studentized family now uses the add-one rank, so
   * its resolution and its floor are both 2/(completed + 1) and it never prints 0; the
   * coefficient-draw families keep the plain share and can still print 0
   if "`e(boot_tstat)'" == "1" dis as result "  The studentized coefficient p-values are add-one ranks over " e(boot_completed) " completed replications, so their resolution and their smallest attainable value are both 2/(n + 1) = " %5.3f 2/(e(boot_completed) + 1) ": a studentized p is never exactly 0. The joint, equality and adversarial p-values are add-one shares over their own listwise-complete replication counts (their smallest value is 1/(n + 1), never 0)."
   else if "`wdbdsp'" != "wald" dis as result "  The rank-based coefficient p-values have resolution 2/" e(boot_completed) " = " %5.3f 2/e(boot_completed) ": a p printed as 0 means p below that (no completed replication on the far side of the null); it is stored as 0 in e(). The joint, equality and adversarial p-values are add-one shares over their own listwise-complete replication counts (their smallest value is 1/(n + 1), never 0)."
   if e(boot_v_posted) == 1 dis as result "  e(V) holds the bootstrap covariance (complete-replication listwise, e(boot_N_V) replications; the analytic covariance, when one existed, remains in the hidden e(V_analytic))."
   else dis as result "  The bootstrap covariance is stored in e(V_boot); no covariance is posted to e(V) for this run."
   * e(V) and e(df_r) carry the sd-family inference whichever family is displayed, so a
   * hand-run test/lincom reproduces the sd package and NOT a displayed percentile p
   * PS-2 (author ruling 2026-09-10): name the EXPORTERS here too. test and lincom are the
   * commands a user runs deliberately; esttab, coefplot and honestdid read the same e(V)
   * without being asked to, so a table or a figure built from this run carries the
   * standard-deviation family even though the printed table does not. That is the whole
   * content of the ruling - the stored family itself is unchanged, per the 2026-09-03/04
   * display rulings that an interval-implied standard error is not a standard error.
   if "`wdbdsp'" != "wald" | "`e(boot_tstat)'" == "1" dis as result "  e(V) and e(df_r) carry the standard-deviation family referred to the `wdbrefd' distribution, so test and lincom run by hand reproduce the bootdisplay(wald) numbers rather than the columns printed above - and so do the tools that read e(V) without being asked: esttab, coefplot and honestdid export the standard-deviation family, not the columns above."
   if "`e(boot_jt_family)'" == "fixedmetric" dis as result "  The joint pre-trend tests (and the across-subgroup equality tests, when subgroups are present) are rank tests of the replication draws recentred at the observed estimates on the fixed bootstrap-covariance metric, p = add-one share of replications at or beyond the observed statistic (e(boot_jt_family) = fixedmetric; the replication statistics are in e(boot_jointW) and e(boot_sgW)). Under boottstat they would instead be bootstrapped test statistics."
   if "`e(boot_jt_absent)'" != "" | "`e(boot_sgjt_absent)'" != "" | "`e(boot_adv_absent)'" != "" dis as result "  Some joint tests, across-subgroup equality tests and/or adversarial line tests are NOT reported for this run (no usable replication statistic - too few complete replications or a rank-deficient bootstrap covariance of the tested columns). They are absent rather than shown as analytic tests, which would not be bootstrap quantities."
   if "`e(boot_tstat)'" == "1" dis as result "  Subgroup contrasts, when requested, are studentized by their own analytic pairwise standard error in every replication (boottstat) and carry the same package as the coefficient rows."
   else dis as result "  Subgroup contrasts, when requested, are bootstrapped as coefficient differences and carry the same package as the coefficient rows."
 }

end

* ============================================================================
* The post-estimation collinearity screens: the dropped treated cohort-period dummy
* check (error 322, with the poisson separation branch of bug sweep BS-4), the
* interaction-drop caution, the SI-A4/R6 continuous-treatment dose-slope screen, and the
* BS-3 interactivecontrols twin. Moved verbatim out of wooldid_worker (per-program
* compiled-size ceiling) on 2026-09-07; the gate stays at the call site so safety(off)
* never enters, and the three locals the R6 branch sets come back by c_local.
* Deliberately PLAIN (not rclass, not eclass): the screens read colnames e(b) of the
* still-active stage-1 model, and an rclass return would blank the caller's r() on the
* path where the inline code left it untouched. `bvars' is recomputed here rather than
* passed, since colnames e(b) grows with the cohort-by-period cell count and could
* overflow a command line on high-cardinality i.
cap prog drop wooldid_ccscreen
prog def wooldid_ccscreen
version 16.0
  syntax , iv(name) tv(name) tf(name) wt(string) es(name) ifs(string asis) mint(real) ui(name) ut(name) [lat(string asis) ct(string) ics(string) icn(string) VERbose POISson FASTnose]
  * re-establish the worker's own names so the screen body below is byte-identical to
  * the code it replaced
  local _i `iv'
  local _t `tv'
  local targetflag `tf'
  local w `wt'
  local esamp `es'
  local ifstatement `ifs'
  local min_t `mint'
  local i `ui'
  local t `ut'
  local lattice `lat'
  local contreat `ct'
  local interlist `ics'
  local interactivecontrols `icn'
  * THIS IS A HARD CEILING ON THE CELL COUNT wooldid can run at, and it binds here and at
  * the worker's own earlier copy (which is reached first, so safety(off) does not avoid
  * it). MEASURED (Stata 19 MP, set maxvar 5000): `local x : colnames M' on a 1 x 15,000
  * matrix whose names total 543,895 characters returns rc 0; at 20,000 such names it fails
  * r(1000) and leaves the local EMPTY rather than truncating. The binding cap is the macro
  * length limit Stata reports as 645,200 characters (it grows 129 characters per unit of
  * set maxvar). At roughly 33 characters per expanded i#t#targetflag name that is about
  * 19,500 treated cell terms - reachable at, say, 1,000 cohorts x 20 periods. The failure
  * is a bare r(1000) raised AFTER the expensive regression, naming no wooldid concept
  local bvars: colnames e(b)
  local catastrophe 0
  local bigcatastrophe 0
  * "o." appears in e(b) colnames iff some term carries an omission flag (a dot only ever
  * follows a factor operator, and a token cannot straddle whitespace), so when it is
  * absent the loop below provably cannot fire -- skip the whole interpreted scan
  if strpos("`bvars'", "o.") != 0 {
  tokenize `bvars'
  local maxwordcount = wordcount("`bvars'")
  forval testbvar = 1(1)`maxwordcount' {
    if strpos("``testbvar''", "o.") != 0   & strpos("``testbvar''", "0b.") == 0 & strpos("``testbvar''", "`_i'") != 0 & strpos("``testbvar''", "`_t'") != 0  & ( strpos("``testbvar''", "1o.`targetflag'") != 0 | strpos("``testbvar''", "1.`targetflag'") != 0)  {
      * the above triggers if and only if (a) the variable has an omission flag (and thus was dropped), and
      * (b) it contains _i, _t, and targetflag (and thus is associated with a main treatment effect), and
      * (c) targetflag (omitted or not) has a value of 1 (meaning it is supposed to correspond with a treatment effect), and
      * (d) base values of 0 are not present anywhere (meaning it is supposed to correspond with a treatment effect)
      * note that these conditions are required in combination since some nonsense cells can form with values like "treatflag = 1 for _i == 0 ", which can't occur in the data
      * note also that this doesn't always work very well for continuous treatments - it is really only high value for the regular non contreat items

      * below, we will test to confirm that there isnt just something weird going on w/ duplicative entries --
      * check that the variable minus the omitteds isn't present somewhere.
      * note that the white space below actually matters!
      local checkifpresent = subinstr("``testbvar''","o.",".",.)
      local checkifpresent = subinstr("`checkifpresent'"," ","",.)
      if strpos("`bvars'"," `checkifpresent' ") == 0 {
        * passing into this if means that the omitted object isn't just a duplicate of any actually existing object (which shouldn't usually happen anyway)
        * next, figure out exactly where the problem occurred
        * proceed only if there is serious reason to believe the coef should've existed

        * strip out the _i _t targetflag etc from the word as well as the hashes
        * what's left behind should be the values of _i and _t associated with the missing object
        foreach word in b. . `_i' `_t' `targetflag' {
          local checkifpresent = subinstr("`checkifpresent'","`word'","",.)
        }
        local checkifpresent = subinstr("`checkifpresent'","#"," ",.)
        local problem_i = word("`checkifpresent'",1)
        local problem_t = word("`checkifpresent'",2)
        * 2 remaining cases. wordcount = 3, in which case what we have is probably a dropped i-t dummy
        * or wordcount > 3, in which case we have something else, like a contreat interaction
        local wordcountproblem = wordcount("`checkifpresent'")

        * are there observations falling in the _i _t cell associated with the missing coefficient? if not, this isn't a problem
        * note: this regrettably is a bit computationally intensive -- and its result is
        * consumed only through the r(N) gate below, whose branches ALL require
        * wordcount <= 3, no lattice, or verbose: when the flagged term is an interaction
        * (wordcount > 3) with a lattice active and verbose off, nothing can consume the
        * count, so skip the full-data pass entirely (lattice drops there are expected)
        if `wordcountproblem' <= 3 | "`lattice'" == "" | "`verbose'" != "" {
        qui count if `_i' == `problem_i' & `_t' == `problem_t' & `targetflag' == 1 & `w' > 0 & `ifstatement'
        if r(N) > 0 {
          * if we proceeded to this point, then the omitted object has associated observations, so this is an actual problem!

          if `wordcountproblem' <= 3 {
            * UNIT CHANGE on this line: `problem_t' arrives in internal _t units (parsed
            * out of a coefficient name; 1-based, 0 = reference/never-treated) and leaves in
            * RELATIVE-HORIZON units, which is what the user-facing message prints. The
            * conversion is the package's load-bearing identity relyr = _t - min_t - 1, set
            * where _t is shifted positive so it can serve as a factor variable. The poisson
            * separation count below writes `+ `min_t' + 1' to undo it, because that count
            * runs against the _t variable still in memory. Getting the direction wrong
            * flips the diagnosis without any error
            local problem_t = `problem_t' - `min_t' - 1
            * under poisson an omitted cell dummy whose rows ppmlhdfe dropped is perfect
            * separation, not multicollinearity (bug sweep BS-4)
            local wdsep 0
            if "`poisson'" != "" {
              qui count if `_i' == `problem_i' & `_t' == `problem_t' + `min_t' + 1 & `targetflag' == 1 & `esamp' & `ifstatement'
              local wdsep = r(N) == 0
            }
            if `wdsep' dis as error "Error: ppmlhdfe dropped every observation of the treated cohort-period cell `i' == `problem_i', relative-`t' == `problem_t' (perfect separation: the outcome is zero on all of them, so the cell's coefficient does not exist). This is not multicollinearity. Coarsen `i', drop the cohort, or specify safety(off) to proceed with the separated cells excluded."
            else dis as error  "Error: a treated cohort-period dummy coefficient was dropped for multicollinearity."
            if "`verbose'" != "" dis as error  "The specific problematic coefficient is: ``testbvar''"
            if "`verbose'" != "" dis as error  "For reference, _i is `_i', _t is `_t', and targetflag is `targetflag'"
            if !`wdsep' dis as error  "The problem seems to be with the treatment dummy for `i' == `problem_i', relative-`t' == `problem_t'"
            if !`wdsep' & "`contreat'" != "" dis as error  "Specify contreatwithin if you have a continuous treatment intensity variable and wish to partial out any effect of treatment when the treatment intensity is 0."
            dis as error  "You can ignore this check if you specify safety(off), but be careful."
            dis as error  " "
            * verbose not specified triggers error here (rather than later after looping through them all )
            if "`verbose'" == "" error 322
            local bigcatastrophe 1
          }
          if `wordcountproblem' > 3  & ("`lattice'" == "" | "`verbose'" == "verbose"){
            * > 3 implies we have more than _i, _t, and targetflag (the variable that is 1 for obs contributed to an estimate and 0 otherwise) forming the dropped variable
            * this implies it is not a main effect treatment dummy
            local problem_t = `problem_t' - `min_t' - 1
            if `catastrophe' == 0 dis as text  "Caution: a treatment period coefficient interacted with something (likely, a term implementing a continuous treatment variable or an interactive control variable) was dropped for multicollinearity."
            if `catastrophe' == 0 & "`lattice'" != "" dis as text  "If using a lattice, this is likely innocuous."
            if `catastrophe' == 0 dis as text  "This could imply your results are not identified, depending on exactly what was dropped."
            if `catastrophe' == 0 dis as text  "Re-estimate with option verbose for details on what was dropped."
            if "`verbose'" != "" dis as text " "
            if "`verbose'" != "" dis as text  "The specific problematic coefficient is: ``testbvar''"
            if "`verbose'" != "" dis as text  "For reference, _i is `_i', _t is `_t', and targetflag is `targetflag'"
            if "`verbose'" != "" dis as text  "Given we are looking at something with an _i _t (other items here) structure,"
            if "`verbose'" != "" dis as text  "the problem values of `i' and relative-`t' are `problem_i', `problem_t'."
            local catastrophe = `catastrophe' + 1
          }
        }
        }
      }
    }
  }
  }
  if `catastrophe' > 0 dis as text "Note: the multicollinearity caution above was triggered `catastrophe' times total."
  if `bigcatastrophe' == 1 error 322

  * SI-A4/R6 (author-ruled 2026-08-25): when EVERY treated-cell dose-slope term was
  * dropped for multicollinearity (zero live c.contreat treated terms while dropped
  * ones exist - the d3 class, e.g. dose constant within each cohort), the within-cell
  * ame estimand is unidentified and the regression's ame entries are exact structural
  * zeros. Name the refusal, void the ame block INCLUDING the point estimates (the R1
  * scrub keeps b; R6 must not - the zeros are fabricated), and continue: the att side
  * is untouched by construction. Scoping on targetflag==1 exempts untreatedZero's
  * expected zero-dose reference periods by construction. Skipped under fastnose (R1
  * convention) and under safety(off) (this whole screen block).
  if "`contreat'" != "" & "`fastnose'" == "" & strpos("`bvars'", "o.") != 0 {
    local amelive 0
    local amedrop 0
    foreach r6tok of local bvars {
      * omitted continuous interactions carry the co. marker (not a bare c. plus
      * o. elsewhere), so a dose term is c.<contreat> OR co.<contreat>
      if strpos("`r6tok'", "`_i'") != 0 & strpos("`r6tok'", "`_t'") != 0 & (strpos("`r6tok'", "1.`targetflag'") != 0 | strpos("`r6tok'", "1o.`targetflag'") != 0) & (strpos("`r6tok'#", "c.`contreat'#") != 0 | strpos("`r6tok'#", "co.`contreat'#") != 0) {
        if strpos("`r6tok'", "o.") == 0 local ++amelive
        else local ++amedrop
      }
    }
    if `amelive' == 0 & `amedrop' > 0 {
      dis as text "Caution: every treated cohort-period dose-slope term for the continuous treatment variable `contreat' was dropped for multicollinearity (see the caution counts above); this occurs, for example, when `contreat' is constant within each cohort. The within-cell average marginal effects (ame) that wooldid reports for continuous treatments are NOT IDENTIFIED in this design - the regression's ame entries are structural zeros, not estimates - so wooldid is reporting the ame results as missing and skipping the ame-based tests (only the raw e(b) coefficient vector retains the structural zeros: Stata cannot post missing coefficients). The att-side results are unaffected."
      dis as text "Consider whether a coarser definition of `i' (so that cohort-period cells contain dose variation) matches your design, or estimate binary att effects without contreat(). You can ignore this check if you specify safety(off), but be careful: with safety(off) the ame output will contain the structural zeros."
      c_local amevoid 1
      c_local amescrub 1
      c_local amescrublo 1
    }
  }
  * interactivecontrols() twin of the screen above (bug sweep BS-3): when no cohort-period
  * cell holds within-cell variation in a control, every treated-cell slope term drops
  * out and the run silently reduces to controls(i.t#c.x); the scan lives in
  * wooldid_icscreen (worker headroom)
  if "`interactivecontrols'" != "" & "`fastnose'" == "" & strpos("`bvars'", "o.") != 0 wooldid_icscreen , bvars(`bvars') iv(`_i') tv(`_t') tf(`targetflag') ics(`interlist') names(`interactivecontrols')
end

* ============================================================================
* Relative-time construction on the POOLED-reference branch: the anticipation shift, the
* usable-window screens, the estimation-sample narrowing, and the birth of the relative
* period variable and the pooled cohort identifier. Moved verbatim out of wooldid_worker
* (per-program compiled-size ceiling) on 2026-09-07. The esfixedbaseperiod twin stays
* inline. `_t', `_i', `ifvar' and `ttre' are the CALLER's tempvars, passed by name and
* mutated in place - a tempvar created here would die on return. Plain program: it reads
* no r() or e() and posts none.
cap prog drop wooldid_reltime
prog def wooldid_reltime
version 16.0
  * the three window integers sit in the OPTIONAL list: Stata rejects a default value on
  * a REQUIRED option (r(197) at the syntax line), and the caller always supplies them
  syntax , ivar(name) tvar(name) ttrevar(name) rawttre(string) ifv(name) relt(name) relid(name) ifs(string asis) [epre(integer 0) epost(integer 0) antic(integer 0) tt(string) sfty(string) ct(string) cw(string)]
  * re-establish the worker's own names so the body below is byte-identical to the code
  * it replaced
  local i `ivar'
  local t `tvar'
  local ttre `ttrevar'
  local ttre_asgiven `rawttre'
  local ifvar `ifv'
  local _t `relt'
  local _i `relid'
  local esprelength `epre'
  local espostlength `epost'
  local anticipation `antic'
  local ifstatement `ifs'
  local timetrends `tt'
  local safety `sfty'
  local contreat `ct'
  local contreatwithin `cw'

  *  calculate the minimum time unit observed w/in the sample for each unit i
  * make a variable to drop all i where we see no pre period
  tempvar mint

  qui gegen `mint' = min(`t') if `ifstatement', by(`i')


  * if requesting time trends, figuring the 2nd from the minimum time unit observed for each i
  if "`timetrends'" == "timetrends" {
    tempvar mintp1
    tempvar tp1
    qui gen double `tp1' = `t' if `t' != `mint'
    qui gegen `mintp1' = min(`tp1' ) if `ifstatement' , by(`i')
    qui drop `tp1'
  }

  * drop cohorts which do not have a pre-treatment period (or 2 such periods if using time trends)
  * F4: before applying the drops, report units that are dropped only BECAUSE of the
  * anticipation shift (shifted onset at/before their first usable period while the
  * recorded onset is not) - otherwise the only signal is a smaller N in the header
  if `anticipation' > 0 {
    local anticmintvar `mint'
    if "`timetrends'" == "timetrends" local anticmintvar `mintp1'
    qui count if `ttre' <= `anticmintvar' & `ttre_asgiven' > `anticmintvar' & `ttre_asgiven' < . & `ifstatement'
    if r(N) > 0 {
      local anticdroppedobs = r(N)
      tempvar antictag
      qui gegen `antictag' = tag(`i') if `ttre' <= `anticmintvar' & `ttre_asgiven' > `anticmintvar' & `ttre_asgiven' < . & `ifstatement'
      qui count if `antictag' == 1
      dis as text "Caution: anticipation(`anticipation') places the effective onset of `r(N)' unit(s) at or before their first usable pre-period; those unit(s) and their `anticdroppedobs' observations are dropped from the sample."
      qui drop `antictag'
    }
  }
  * the drop is UNIT-level by construction (gegen max by i): with a within-unit-varying
  * ttre (reachable only under safety(off), F13 forbids it otherwise) the old row-level
  * form could remove a surviving unit's t==mint rows, leaving the per-unit mint/mintp1
  * below stale and the reference-forcing clauses targeting a period that no longer
  * exists - a silent collinearity route (dropped treated dummies, biased estimates)
  tempvar nopre
  if "`timetrends'" == "" qui gegen byte `nopre' = max(`ttre' <= `mint') if `ifstatement', by(`i')
  if "`timetrends'" == "timetrends" qui gegen byte `nopre' = max(`ttre' <= `mintp1') if `ifstatement', by(`i')
  qui replace `ifvar' = 0 if `nopre' == 1
  qui drop `nopre'
  qui replace `ttre' = . if `ifvar' == 0
  * invariant guard: recompute mint/mintp1 on the post-drop sample. with the unit-level
  * drop above this is a provable no-op (drops remove whole units, so surviving units'
  * first periods are unchanged); it stays as insurance that the reference-forcing
  * clauses below can never target a period absent from the surviving sample
  qui drop `mint'
  qui gegen `mint' = min(`t') if `ifstatement', by(`i')
  if "`timetrends'" == "timetrends" {
    tempvar tp1g
    qui drop `mintp1'
    qui gen double `tp1g' = `t' if `t' != `mint'
    qui gegen `mintp1' = min(`tp1g') if `ifstatement', by(`i')
    qui drop `tp1g'
  }


  * construct _t, which contains dummies for treatment periods of interest - it will be specified in terms of relative periods
  * _t is born as a relative period variable, which we missing out for relyrs before the event study prelength and for the very first observation we have for everyone
  * note that if the first observation for a unit is a period for which treatment effects are requested, this observation still ends up being the reference period
  * this creates much mischief in interpreting event study results when it happens and using the fixed base period approach should be used when this is relevant
  * the third replace statement just insures that the minimum value is always 1; the 4th 0s the missings

  qui gen `_t' = `t' - `ttre' if `ifstatement'

  if "`timetrends'" == "" {
    * the caution branch runs whenever a pre-treatment horizon is requested at all
    * (bug sweep BS-2: <= 1 gated it off at esprelength(1), the smallest case it exists for)
    if "`safety'" == "off" | `esprelength' == 0 {
      qui replace `_t' = . if `_t' < -`esprelength'  | `t' == `mint'
    }
    else{
      qui replace `_t' = . if `_t' < -`esprelength'
      cap assert `_t' == . if  `t' == `mint' &  `ifstatement'
      if _rc != 0 {
        qui replace `_t' = . if `t' == `mint'
        dis as text " "
        dis as text "Caution: reference period for some units forced (for want of panel length) to include time period for which treatment effects were requested."
        dis as text "This does not break estimation, but event study estimates can carry counterintuitive meanings and reflect effect levels by cohort."
        dis as text "Recommendation: use option esfixedbaseperiod to change estimation approach to using a specific time period relative to treatment as your reference period."
        if `esprelength' > 1 dis as text "(Or request estimates for fewer pre-treatment horizons.)"
        dis as text " "
        * If this condition triggers, event study coefficients are very difficult to interpret
      }
    }
  }
  else{
    * if using timetrends, never use the 2nd observation either - dont just always skip the first!
    * the caution branch runs whenever a pre-treatment horizon is requested at all
    * (bug sweep BS-2: <= 1 gated it off at esprelength(1), the smallest case it exists for)
    if "`safety'" == "off" | `esprelength' == 0 {
      qui replace `_t' = . if `_t' < -`esprelength'  | `t' <= `mintp1'
    }
    else{
      qui replace `_t' = . if `_t' < -`esprelength'
      cap assert `_t' == . if  `t' <= `mintp1' &  `ifstatement'
      if _rc != 0 {
        qui replace `_t' = . if `t' <= `mintp1'
        dis as text " "
        dis as text "Caution: reference period for some units forced (for want of panel length) to include time period for which treatment effects were requested."
        dis as text "This does not break estimation, but event study estimates can carry counterintuitive meanings and reflect effect levels by cohort."
        dis as text "Recommendation: use option esfixedbaseperiod to change estimation approach to using a specific time period relative to treatment as your reference period."
        if `esprelength' > 1 dis as text "(Or: request estimates for fewer pre-treatment horizons, or drop the request to include timetrends.)"
        else dis as text "(Or drop the request to include timetrends.)"
        dis as text " "
      }
    }
    qui drop `mintp1'
  }


  * update _t and ttre
  qui replace `_t' = . if `ifvar' == 0
  qui replace `ttre'  = . if `ifvar' == 0

  * produce a check that each _t interaction will be estimable in that there is an untreated cohort for it
  * note that you can have periods with only treated units in them in the pre period provided they are periods pooled into an pooled reference period
  * allowing for the above is why we use _t instead of ttre here like in the ensuing section
  * do not apply test if not using a purely untreated control group (contreatwithin case)
  if "`contreatwithin'" == "" {
    tempvar controlcheck
    tempvar controlcheck2
    qui gen `controlcheck' = `_t' ==. & `ifstatement'
    qui gegen `controlcheck2' = max(`controlcheck') if `ifstatement' , by(`t')
    qui replace `ifvar'  = 0 if `controlcheck2' == 0

    * update _t and ttre
    qui replace `_t' = . if `ifvar' == 0
    qui replace `ttre'  = . if `ifvar' == 0

    if "`safety'" != "off"  {
      cap assert `controlcheck2' != 0
      if _rc != 0 {
        dis as text "Caution: estimability checks found time periods exist in sample without a control unit, causing them to be dropped."
        if "`contreat'" != "" dis as text "If intending to only use within treated-cohorts variation in your continuous treatment variable, specify option contreatwithin to prevent these time periods from being dropped."
      }
    }
    qui drop `controlcheck' `controlcheck2'
  }

  *get min viable _t and max; set up relative period variable to have no negative values
  qui gstats sum `_t' if `ifstatement', meanonly
  * min_t counts the realized pre-treatment horizons; when no treated unit is observed at or
  * before its onset (an onset-period gap, or an onset off the panel's time grid) r(min) is
  * positive and abs() would invent a pre-window that nothing realizes (bug sweep BS-33)
  local min_t = max(0, -r(min))
  local max_t = r(max) + 1
  qui replace `_t' = `_t' + `min_t' + 1
  qui replace `_t' = 0 if `_t' == . & `ifstatement'

  * check that some cohort has a post period and that some has a pre-period
  cap assert `max_t' > 0 & `max_t' < .
  if _rc != 0 {
    dis as error  "Error: cannot find either a treated time period for any treated units or sufficient untreated periods for treated units."
    if "`timetrends'" != "" dis as error  "Note that including cohort-specific time trends raises the number of untreated periods required."
    if `anticipation' > 0 dis as error  "Note: anticipation(`anticipation') shifted every effective onset earlier; a shift larger than the available pre-period window can leave no estimable treated periods."
    error 459
  }
  * produce modified identifiers _i
  * _i is unit identifiers that are positive integers if i is ever treated in sample and 0 otherwise
  qui gen double `_i' = `i' if  `_t' > 0 & `_t' < . & `ifstatement'
  qui replace `_i' = 0 if `_i' == . & `ifstatement'

  local espl = min(`max_t',`espostlength')
  local espre = `min_t'
  drop `mint'
  * end data setup work in the case not using relyr = -1 as a fixed base period
  c_local min_t `min_t'
  c_local max_t `max_t'
  c_local espl `espl'
  c_local espre `espre'
end

* ============================================================================
* Subgroup-domain validation and the realized subgroup list: the integrality and
* within-unit-constancy screens (F13 family), the split-unit ruling under safety(off),
* the name-length caution, and glevelsof over the treated estimation rows. Moved
* verbatim out of wooldid_worker (per-program compiled-size ceiling) on 2026-09-07;
* the caller keeps the gate, so a run without subgroups() never enters. Plain program,
* NOT rclass: glevelsof's r() residue is left exactly as the inline code left it. The
* only thing the worker needs back is `sglist'.
cap prog drop wooldid_sgdom
prog def wooldid_sgdom
version 16.0
  syntax , sg(name) ivar(name) tvar(name) ttrevar(name) ifs(string asis) [si(string) MAKECfxplotsbysg]
  * re-establish the worker's own names so the body below is byte-identical to the code
  * it replaced
  local subgroups `sg'
  local i `ivar'
  local t `tvar'
  local ttre `ttrevar'
  local ifstatement `ifs'
  local safetyinput `si'
  * VAL-G15A: negative-value screen, moved here from wooldid_validate. domain =
  * in-sample rows (ifstatement folds the deterministic markout, curing false
  * rejections on e.g. missing-outcome rows) of units with a nonmissing internal
  * onset - both pre- and post-onset rows, because the SG pre-trend margins consume
  * pre-onset rows too; strictly looser than the old raw-if domain (never-treated,
  * marked-out, and out-of-window-onset rows newly tolerated; nothing newly rejected).
  * placed before the sglist screens so negative values error first and never enter
  * sglist. the integer screen below keeps its own adjudicated t >= ttre domain
  cap assert `subgroups' >= 0 | `subgroups' >= . if `ifstatement' & `ttre' < .
  if _rc != 0 {
    dis as error "Error: the subgroups() variable contains negative values on in-sample rows of treated units (pre- or post-onset); subgroup values must be non-negative integers (or missing) there."
    error 411
  }
  qui glevelsof(`subgroups') if `ifstatement' & `t' >= `ttre' , clean local(sglist)
  if wordcount("`sglist'") == 0 {
    dis as error  "Error: subgroups() is specified, but the subgroup variable contains no values within sample during the treatment period."
    error 459
  }
  * subgroup values reaching the treated estimation sample must be integer-valued: a non-integer
  * value yields an sglist token like 1.5 that later crashes the naming machinery (locals,
  * matrices, ereturns named ..._s1.5) with an unexplained r(198). the check deliberately uses
  * the levelsof domain above, so values appearing only outside it (never-treated units, screened
  * observations) stay tolerated exactly as before, and missing values pass (int(.) == .)
  cap assert `subgroups' == int(`subgroups') if `ifstatement' & `t' >= `ttre'
  if _rc != 0 {
    dis as error "Error: the subgroup variable contains non-integer values within the treated estimation sample; it must contain only non-negative integers."
    error 411
  }
  * F13: subgroup containment. sg never enters the regression - the unit-by-period cell
  * coefficients are shared by every subgroup a unit touches. sg varying WITHIN a unit-period
  * cell therefore makes by-subgroup contrasts structurally unable to reflect within-cell
  * heterogeneity (each subgroup only reweights the same pooled cell effects; contrasts can
  * come back identical under strong true heterogeneity) -> error. sg varying within a unit
  * only ACROSS periods leaves by-sg att/att_it correct (row-weighted) but makes
  * att_i/att_itime weight a straddling unit by its weight share rather than equally ->
  * caution. domains: missing sg rows belong to no subgroup by design and never-treated
  * units never enter any sg subpop, so both are excluded (cond() folds extended missings)
  if "`safetyinput'" != "off" {
    tempvar f13sgitmin f13sgitmax
    qui gegen double `f13sgitmin' = min(cond(`subgroups' >= ., ., `subgroups')) if `ifstatement' & `ttre' != ., by(`i' `t')
    qui gegen double `f13sgitmax' = max(cond(`subgroups' >= ., ., `subgroups')) if `ifstatement' & `ttre' != ., by(`i' `t')
    cap assert `f13sgitmin' == `f13sgitmax' if `ifstatement' & `ttre' != .
    if _rc != 0 {
      dis as error "Error: the subgroup variable `subgroups' varies within `i'-by-`t' cells of the treated estimation sample."
      dis as error "wooldid estimates effects at the `i'-by-period level, so subgroup estimates built from cells that mix subgroups cannot reflect within-cell subgroup differences - each subgroup only reweights the same pooled cell effects, and subgroup contrasts can come back identical even when true effects differ sharply by subgroup."
      dis as error "Define `i' at the finer level over which the subgroup varies (e.g., egen newi = group(`i' `subgroups')), or drop the subgroups() option."
      dis as error "To ignore this error, specify safety(off)."
      error 459
    }
    qui drop `f13sgitmin' `f13sgitmax'
  }
  * a unit split across multiple subgroups over time (subgroup varies within i across
  * periods) is refused under default safety - a unit cannot belong to more than one
  * subgroup, and the share-weighted fallback deserves an explicit opt-in. under
  * safety(off) the run proceeds with share weighting and a caution (author directive,
  * 2026-07-28; this check runs regardless of the safety setting)
  tempvar f13sgimin f13sgimax
  qui gegen double `f13sgimin' = min(cond(`subgroups' >= ., ., `subgroups')) if `ifstatement' & `ttre' != ., by(`i')
  qui gegen double `f13sgimax' = max(cond(`subgroups' >= ., ., `subgroups')) if `ifstatement' & `ttre' != ., by(`i')
  cap assert `f13sgimin' == `f13sgimax' if `ifstatement' & `ttre' != .
  if _rc != 0 {
    if "`safetyinput'" != "off" {
      dis as error "Error: the subgroup variable `subgroups' splits some units of `i' across multiple subgroups over time. A unit cannot belong to more than one subgroup."
      dis as error "To proceed anyway, specify safety(off): by-subgroup att and att_it then weight each observation to its own subgroup as usual, while att_i and att_itime enter a split unit into each subgroup weighted by its share of sample weight there, rather than with the equal unit weighting those estimands normally imply."
      error 459
    }
    else dis as text "Caution: the subgroup variable `subgroups' splits some units of `i' across multiple subgroups over time; proceeding under safety(off) using share weighting. att_i and att_itime enter a split unit into each subgroup weighted by its share of sample weight there, rather than with the equal unit weighting those estimands normally imply; by-subgroup att and att_it results are unaffected."
  }
  qui drop `f13sgimin' `f13sgimax'
  * long subgroup values overflow stata's 31/32-character name limits when embedded in
  * constructed local/matrix/rowname/graph/e() names (the worst prefixes reach 24-25
  * characters: the makecfxplotsbysg tempnames, the adversarial-by-sg ereturns, and the
  * pooled pretreatment_SG<value>_att_itime rownames all crash at 8-digit values; the
  * bysg tempnames also embed the relative year, so with 3-digit relyr plot requests
  * even 7-digit values overflow), and values beyond 2,147,483,647 exceed the
  * factor-variable level cap that margins over() relies on. whether a given run
  * crashes depends on the options requested, so this is a caution rather than an
  * error. the sglist domain matches the checks above: only values reaching the
  * treated estimation sample are ever injected into names
  local sgmaxlen 0
  foreach sgcase in `sglist' {
    if strlen("`sgcase'") > `sgmaxlen' local sgmaxlen = strlen("`sgcase'")
  }
  if `sgmaxlen' >= 8 | (`sgmaxlen' == 7 & "`makecfxplotsbysg'" != "") {
    dis as text "Caution: the subgroup variable takes values up to `sgmaxlen' digits long. wooldid embeds subgroup values in the names of stored results, matrices, and plots, and Stata caps such names at 31-32 characters, so values 8 or more digits long (7, with makecfxplotsbysg and 3-digit relative-year plot requests) can sometimes cause unexplained r(198) invalid-name crashes; att_itime results tables, adversarial, and makecfxplotsbysg are most exposed. Values above 2,147,483,647 will generally fail outright at Stata's factor-variable level cap."
    dis as text "If wooldid crashes, recode ID-style subgroup values to small integers first (e.g., egen newsg = group(oldsg))."
  }
  * NOTE a loop initialising joint_horiz_i_s<sg> / joint_Choriz_i_s<sg> used to sit here and
  * was DELETED as inert (D-5, author ruling 2026-09-11): those locals are scoped to this
  * program and die at its `end', so the caller never received them, and it was already a
  * no-op before the 2026-09-07 move out of wooldid_worker. The caller's guards test for the
  * EMPTY string and an unset local reads as empty, which is what makes the absence correct.
  * Do not "repair" this by adding c_locals - they would overwrite values the caller sets
  * for itself
  c_local sglist `sglist'
end

* ============================================================================
* The histogramcohorteffects() diagnostic: per-cohort and per-cohort-period marginal
* effects, their histograms, and the matrix behind e(histogramestimates). Moved verbatim
* out of wooldid_worker (per-program compiled-size ceiling) on 2026-09-07. `___hist' is
* created by the CALLER and passed in by name: a tempname minted here would be dropped on
* return and the ereturn of e(histogramestimates) would fail. Plain program - it reads the
* fitted model from e() and must not post into it, and its margins calls consume their own
* r() immediately, so nothing outside depends on r() here.
cap prog drop wooldid_hist
prog def wooldid_hist
version 16.0
  syntax , h(name) iv(name) tv(name) tf(name) ttrevar(name) wtype(string) semi(string) ifs(string asis) dfr(name) [att att_it att_i att_itime VERbose replace wa(string) wai(string) wait(string) ct(string) cet(string) pm(string asis) pmc(string asis) fcn(string) sp(string) elbl(string asis) celbl(string asis)]
  * re-establish the worker's own names so the body below is byte-identical to the code
  * it replaced (`whichweight' resolves w_att / w_att_i by double dereference, so both
  * weight variables must exist under their original names here)
  local ___hist `h'
  local _i `iv'
  local t `tv'
  local targetflag `tf'
  local ttre `ttrevar'
  local weight `wtype'
  local semielasticity `semi'
  local ifstatement `ifs'
  local defaultframe `dfr'
  local w_att `wa'
  local w_att_i `wai'
  local w_att_it `wait'
  local contreat `ct'
  local contreatelasticitytype `cet'
  local poismarginstext `pm'
  local poismarginstextC `pmc'
  local fcnoes `fcn'
  local saveplots `sp'
  local effecttype `elbl'
  local Ceffecttype `celbl'
  if "`verbose'" == "verbose" dis as text "Wooldid is now producing histograms."
  * graph window management: graphs are created with name(..., replace) so re-runs overwrite wooldid's own prior graphs; wooldid no longer clears the user's unrelated graphs from memory

  * we pool att with att_i b/c they conceptually should match into the same histograms
  if "`att_i'" == "att_i" |  "`att'" == "att"  {
    * WHY A DIFFERENT WEIGHT GIVES THE SAME PICTURE - this is exact, not approximate.
    * w_per_i and t_per_i are totals by(posttreatment, i), and every margins call here
    * carries subpop(... & `t' >= `ttre'), which pins posttreatment to 1 - so within each
    * over(i) group the divisors are constant, making w_att_i proportional to w_att and
    * w_att_itime proportional to w_att_it. margins' over() margin is a within-group
    * weighted mean and is invariant to a positive within-group rescaling of the weights.
    * The substitution exists only so a weight that may not have been generated is never
    * referenced (w_att needs att; w_att_itime needs att_itime). THE TRAP: widen the subpop
    * to include pre-periods and the divisors stop being constant within the group, and the
    * two weights stop agreeing - silently, since nothing here checks
    local whichweight w_att
    if "`att'" == "" local whichweight w_att_i
    qui margins [`weight'=``whichweight''], `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext' `fcnoes'
    * the SECOND HALF of r(b): targetflag enters as a factor variable, so margins returns
    * 2 x (number of over groups) columns - the 0.targetflag base first (always the omitted
    * level), the 1.targetflag effects second. The contreat call below reads ALL of r(b)
    * because contreat is continuous and there is no half to skip; "harmonizing" the two
    * would silently fill the histogram, and e(histogramestimates), with base-level margins.
    * The row names of the result are margins' column subnames and therefore carry wooldid's
    * own tempvars - volatile every run. That is deliberate scope (the help says the
    * histogram itself is the intended output); the cfx matrices get stable keys from
    * wooldid_cfxstripe and this matrix pointedly does not
    local matrixstartgrab floor(`=colsof(r(b))'/2)+1
    matrix `___hist' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
    matrix colnames `___hist' = att_i
    if "`contreat'" != "" {
      qui margins [`weight'=``whichweight''], `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstextC' `fcnoes'
      tempname ___Chist
      matrix `___Chist' = r(b)[1, 1..`=colsof(r(b))']'
      matrix `___hist' = (`___hist',`___Chist')
      matrix drop `___Chist'
      matrix colnames `___hist' = att_i ame_i
    }
    tempname wooldidframe
    frame create `wooldidframe'
    frame change `wooldidframe'
    qui svmat double `___hist'
    rename `___hist'1 att_i
    if "`contreat'" != "" rename `___hist'2 ame_i
    local bininserter
    if _N < 100 local bininserter bin(20)
    qui hist att_i, freq   xtitle("Cohort-Specific Average Treatment `effecttype'") name(hist_att_i, replace) gap(5) `bininserter' color(navy)
    if "`saveplots'" != "" qui graph save "`saveplots'_hist_att_i.gph", `replace'
    if "`contreat'" != "" qui hist ame_i, freq   xtitle("Cohort-Specific Average Marginal Treatment `Ceffecttype'") name(hist_ame_i, replace) gap(5) `bininserter' color(navy)
    if "`contreat'" != "" & "`saveplots'" != ""  qui graph save "`saveplots'_hist_ame_i.gph", `replace'
    frame change `defaultframe'
    frame drop `wooldidframe'
  }

  * bundle these for conceptually both generating the same histograms
  if "`att_itime'" == "att_itime" |  "`att_it'" == "att_it"  {
    * ONE ROW PER TREATED COHORT, under every estimand: this call is over(`_i') exactly as
    * the att/att_i call above is, so the att_it/att_itime column is the COHORT-LEVEL
    * average of that cohort's i-t cell effects taken under the cell weights - not a
    * distribution over cohort-by-period cells, and there is no per-cell object anywhere in
    * this program. The equal row count is load-bearing: it is the only reason the
    * horizontal concatenation below is conformable rather than r(503)
    qui margins [`weight'=`w_att_it'], `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext' `fcnoes'
    local matrixstartgrab floor(`=colsof(r(b))'/2)+1
    tempname ___hist2
    matrix `___hist2' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
    matrix colnames `___hist2' = att_itime
    if "`contreat'" != "" {
      tempname ___Chist2
      qui margins [`weight'=`w_att_it'], `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstextC' `fcnoes'
      matrix `___Chist2' = r(b)[1, 1..`=colsof(r(b))']'
      matrix `___hist2' = (`___hist2',`___Chist2')
      matrix drop `___Chist2'
      matrix colnames `___hist2' = att_itime ame_itime
    }
    tempname wooldidframe
    frame create `wooldidframe'
    frame change `wooldidframe'
    qui svmat double `___hist2'
    rename `___hist2'1 att_itime
    if "`contreat'" != "" rename `___hist2'2 ame_itime
    local bininserter
    if _N < 100 local bininserter bin(20)
    qui hist att_itime, freq   xtitle("Cohort-Specific i-t Cell Average Treatment `effecttype'") name(hist_att_itime, replace) gap(5) `bininserter' color(navy)
    if "`saveplots'" != "" qui graph save "`saveplots'_hist_att_itime.gph", `replace'
    if "`contreat'" != "" qui hist ame_itime, freq   xtitle("Cohort-Specific i-t Cell Average Marginal Treatment `Ceffecttype'") name(hist_ame_itime, replace) gap(5) `bininserter' color(navy)
    if "`contreat'" != "" &   "`saveplots'" != ""   qui graph save "`saveplots'_hist_ame_itime.gph", `replace'
    frame change `defaultframe'
    frame drop `wooldidframe'
    if "`att_i'" == "att_i" |  "`att'" == "att"  {
      matrix `___hist' = (`___hist',`___hist2')
    }
    else{
      matrix `___hist' = `___hist2'
    }
    cap matrix drop `___hist2'
  }
  if "`verbose'" == "verbose" dis as text "Wooldid has finished producing histograms."

end

* ============================================================================
* The summarystats table: outcome and continuous-treatment distributions over the
* estimation sample, its pre and post halves, and the comparison rows. Moved verbatim out
* of wooldid_worker (per-program compiled-size ceiling) on 2026-09-08. The matrix name is
* minted by the CALLER and passed in: a tempname created here would be dropped on return
* and the ereturn of e(summarystats) would fail. The reference-period branch stays in the
* worker so the temp-name allocation order is untouched. Plain program, NOT rclass: the
* gstats calls consume their own r() line by line and the caller's r() is left as the
* inline code left it. Only the thirteen locals the e() posts read come back.
cap prog drop wooldid_sumstat
prog def wooldid_sumstat
version 16.0
  syntax , ss(name) y(name) tvar(name) ttrevar(name) [ct(string) ctc(string) ctw(string) wi(string asis)]
  * re-establish the worker's own names so the body below is byte-identical to the code
  * it replaced (the stat locals are built by double dereference off `sumvar', so `y' and
  * `contreat' must exist under exactly those names here)
  local summarystats `ss'
  local t `tvar'
  local ttre `ttrevar'
  local contreat `ct'
  local contreatcontrols `ctc'
  local contreatwithin `ctw'
  local weightinserter `wi'
local rownamessum
if "`contreat'" == "" local sumvarlist y
if "`contreat'" != "" local sumvarlist y contreat
* gstats sum accepts pweights directly, so the regression weights pass through to these descriptive computations unchanged
* Two consequences worth knowing before quoting this table. Official summarize does NOT
* accept pweights at all (measured: r(101) "pweights not allowed"), so a user cannot
* reproduce e(summarystats) with summarize under pweights. And gstats' pweighted dispersion
* is not its aweighted dispersion: measured on y = 1..10 with w = 1,1,1,1,1,100,...,100,
* [pw=w] gives sd 1.4998428 and [aw=w] gives 1.5794071 - the ratio is sqrt(N/(N-1)), i.e.
* pweights take the population divisor and aweights the N-1 one. Means and percentiles agree
* between the two, so only the sd column depends on which weight keyword the USER typed
* (the weight inserter carries their keyword through verbatim).
foreach sumvar in `sumvarlist' {
  local sumvar2 c
  if "`sumvar'" == "y" local sumvar2 y
  qui gstats sum ``sumvar'' `weightinserter' if e(sample), detail
  local `sumvar'_mean = r(mean)
  local `sumvar'_med = r(p50)
  local `sumvar'_sd = r(sd)
  local `sumvar'_p25 = r(p25)
  local `sumvar'_p75 = r(p75)
  local `sumvar'_iqr = ``sumvar'_p75' - ``sumvar'_p25'
  local `sumvar'_min = r(min)
  local `sumvar'_max = r(max)
  local rownamessum `rownamessum' `sumvar2'_insample

  * do the same in the pre-treatment treated sample
  if "`sumvar'" == "y" | ("`contreatcontrols'" != "untreatedZero" & "`sumvar'" == "contreat") {
    qui gstats sum ``sumvar'' `weightinserter' if e(sample) & `ttre' <. & `t' < `ttre', detail
    local `sumvar'pre_mean = r(mean)
    local `sumvar'pre_med = r(p50)
    local `sumvar'pre_sd = r(sd)
    local `sumvar'pre_p25 = r(p25)
    local `sumvar'pre_p75 = r(p75)
    local `sumvar'pre_iqr = ``sumvar'pre_p75' - ``sumvar'pre_p25'
    local `sumvar'pre_min = r(min)
    local `sumvar'pre_max = r(max)
    local rownamessum `rownamessum' `sumvar2'_pretreat
  }

  * do the same in the post-treatment treated sample
  qui gstats sum ``sumvar'' `weightinserter' if e(sample) & `ttre' <. & `t' >= `ttre', detail
  local `sumvar'post_mean = r(mean)
  local `sumvar'post_med = r(p50)
  local `sumvar'post_sd = r(sd)
  local `sumvar'post_p25 = r(p25)
  local `sumvar'post_p75 = r(p75)
  local `sumvar'post_iqr = ``sumvar'post_p75' - ``sumvar'post_p25'
  local `sumvar'post_min = r(min)
  local `sumvar'post_max = r(max)
  local rownamessum `rownamessum' `sumvar2'_posttreat

  * do the same in the control group (units never treated in-window: `ttre' == . covers never-treated and out-of-window-treated units)
  if "`contreatwithin'" == "" & ("`sumvar'" == "y" | ("`contreatcontrols'" != "untreatedZero" & "`sumvar'" == "contreat"))  {
    qui gstats sum ``sumvar'' `weightinserter' if e(sample) & `ttre' == ., detail
    local `sumvar'con_mean = r(mean)
    local `sumvar'con_med = r(p50)
    local `sumvar'con_sd = r(sd)
    local `sumvar'con_p25 = r(p25)
    local `sumvar'con_p75 = r(p75)
    local `sumvar'con_iqr = ``sumvar'con_p75' - ``sumvar'con_p25'
    local `sumvar'con_min = r(min)
    local `sumvar'con_max = r(max)
    local rownamessum `rownamessum' `sumvar2'_control

  }
}

* store these into a summary stats matrix
matrix `summarystats' = (`y_mean' , `y_med', `y_sd', `y_iqr', `y_p25', `y_p75', `y_min', `y_max')
matrix `summarystats' = `summarystats' \ (`ypre_mean' , `ypre_med', `ypre_sd', `ypre_iqr', `ypre_p25', `ypre_p75', `ypre_min', `ypre_max')
matrix `summarystats' = `summarystats' \ (`ypost_mean' , `ypost_med', `ypost_sd', `ypost_iqr', `ypost_p25', `ypost_p75', `ypost_min', `ypost_max')
if "`contreatwithin'" == "" matrix `summarystats' = `summarystats' \ (`ycon_mean' , `ycon_med', `ycon_sd', `ycon_iqr', `ycon_p25', `ycon_p75', `ycon_min', `ycon_max')

if "`contreat'" != "" {
  matrix `summarystats' = `summarystats' \ (`contreat_mean' , `contreat_med', `contreat_sd', `contreat_iqr', `contreat_p25', `contreat_p75', `contreat_min', `contreat_max')
  if "`contreatcontrols'" != "untreatedZero"  matrix `summarystats' = `summarystats' \ (`contreatpre_mean' , `contreatpre_med', `contreatpre_sd', `contreatpre_iqr', `contreatpre_p25', `contreatpre_p75', `contreatpre_min', `contreatpre_max')
  matrix `summarystats' = `summarystats' \ (`contreatpost_mean' , `contreatpost_med', `contreatpost_sd', `contreatpost_iqr', `contreatpost_p25', `contreatpost_p75', `contreatpost_min', `contreatpost_max')
  if "`contreatwithin'" == "" & "`contreatcontrols'" != "untreatedZero"  matrix `summarystats' = `summarystats' \ (`contreatcon_mean' , `contreatcon_med', `contreatcon_sd', `contreatcon_iqr', `contreatcon_p25', `contreatcon_p75', `contreatcon_min', `contreatcon_max')
}

matrix colnames `summarystats' = mean median sd iqr p25 p75 min max
matrix rownames `summarystats' = `rownamessum'
  c_local ypre_mean `ypre_mean'
  c_local ypre_med `ypre_med'
  c_local ypre_sd `ypre_sd'
  c_local ypre_iqr `ypre_iqr'
  c_local contreat_mean `contreat_mean'
  c_local contreat_med `contreat_med'
  c_local contreat_sd `contreat_sd'
  c_local contreatpre_mean `contreatpre_mean'
  c_local contreatpre_med `contreatpre_med'
  c_local contreatpre_sd `contreatpre_sd'
  c_local contreatpost_mean `contreatpost_mean'
  c_local contreatpost_med `contreatpost_med'
  c_local contreatpost_sd `contreatpost_sd'
end

* BS-8 (2026-09-06): label a subgroup-by-horizon margins call from its OWN colnames.
* positional omission scan (SI-7, author-ruled 2026-08-27): every column carrying margins'
* o. marker is dropped from the six-row table (and from b/V when given), a caution names
* the dropped horizons, and the relative periods of the kept columns come back as the
* column matrix rel(), the rowname membership list lbl() and the strict per-column label
* list f2l(). serves the binary call (first() = the 1.targetflag half, strip(1.targetflag))
* and the continuous-treatment call (kind(continuous-treatment)). lives outside
* wooldid_worker to conserve that program's compiled-size headroom. the caller consumes
* r() first: the single-level fallback below is rclass and wipes it
cap prog drop wooldid_sglbl
prog def wooldid_sglbl
  syntax , names(string) tab(name) mint(real) old(integer) sg(string) tvar(name) fb(string asis) rel(name) lbl(name) [f2l(name) b(name) v(name) first(integer 1) strip(string) kind(string)]
  local kk = cond("`kind'" == "", "", "`kind' ")
  local keep
  local droplbl
  local ndrop 0
  local idx 0
  foreach nm of local names {
    local ++idx
    if strpos("`nm'", "o.") != 0 {
      local ++ndrop
      if regexm("`nm'", "^([0-9]+)") {
        local drel = real(regexs(1)) - `mint' - 1
        if `drel' < 0 local droplbl `droplbl' pre`=abs(`drel')'
        else          local droplbl `droplbl' post`drel'
      }
    }
    else local keep `keep' `idx'
  }
  c_local wdsg_ndrop `ndrop'
  c_local wdsg_skip 0
  c_local `lbl'
  if "`f2l'" != "" c_local `f2l'
  if `ndrop' > 0 {
    dis as text "Caution: the `kk'relative period effect(s) `droplbl' are not estimable for subgroup `sg' (margins reports them as omitted), so they are dropped from that subgroup's `kk'event study; its remaining relative periods, and all other subgroups, are unaffected. This typically arises when the specified fixed effects absorb the subgroup-by-time variation (for example fe() interacting the subgroups() variable with time), or when a subgroup has no usable observations at that relative period."
    if "`keep'" == "" {
      dis as text "Caution: no `kk'relative period effect is estimable for subgroup `sg', so its `kk'event study is skipped entirely."
      c_local wdsg_skip 1
      exit
    }
  }
  * relative periods from the kept names. margins names a single-level over() result
  * without the level (it returns 1.targetflag, or nothing, as the relevant level, because
  * the subpop made over(_t) irrelevant), so fall back to the realized levels of the
  * subpop when a kept name does not parse - a levelsof pass costs more than the text
  * processing on large data, hence not by default
  local kn
  foreach k of local keep {
    local w : word `k' of `names'
    local kn `kn' `w'
  }
  local kn = subinstr("`kn'", "bn", "", .)
  local kn = subinstr("`kn'", "`tvar'", "", .)
  if "`strip'" != "" local kn = subinstr("`kn'", "`strip'", "", .)
  local ok 1
  foreach w of local kn {
    if !regexm("`w'", "^[0-9]+[.]?$") local ok 0
  }
  * TWO THINGS THIS FALLBACK COSTS. (a) Alignment: the recovered levels are paired with
  * `keep' BY POSITION in the loop below, so correctness needs glevelsof's ascending order
  * and a subpop that realizes exactly this margins call's estimable levels. The wordcount
  * check below errors 459 on a size mismatch and cannot see a same-size set mismatch -
  * that case mislabels every horizon silently. (b) r() clobber: this program is NOT rclass,
  * so glevelsof's r() overwrites the CALLER's r() whenever this branch fires. Callers must
  * consume r() before calling wooldid_sglbl
  if `ok' == 0 | wordcount("`kn'") == 0 qui glevelsof(`tvar') if `fb', clean local(kn)
  if wordcount("`kn'") != wordcount("`keep'") {
    dis as error "Error: the `kk'event-study margins for subgroup `sg' returned `=wordcount("`keep'")' estimable columns but `=wordcount("`kn'")' relative-period labels could be recovered for them. This is unexpected - please report it."
    error 459
  }
  local relcol
  local lbls
  local f2lbls
  foreach w of local kn {
    * the package-wide internal-time identity, inverted: _t = relyr + min_t + 1, so
    * relyr = level - min_t - 1. `mint' must be the SAME min_t the margins call was built
    * on, and `w' is the factor-variable token margins put in the colname. An off-by-one
    * here does not error - it relabels every horizon of this subgroup's event study (pre1
    * printed as post0), and those labels then drive rowname membership and the f2l column
    * labels. The drop caution's labels use the identical line, so the two must stay in step.
    * NOTE the r == 0 case under oldsyntax appends TWO tokens to `lbls' (post0 AND contemp)
    * on purpose: the ES row for relative period 0 is named contemp under oldsyntax, and
    * lbl() is a MEMBERSHIP list consumed by a strpos() test, so both spellings must pass.
    * Never index lbl() positionally and never take a column count from it - f2l() is the
    * strict one-token-per-kept-column list. (The worker records what positional use costs:
    * a namelist shorter than the matrix makes Stata pad by repeating the LAST name, r(0)
    * and no warning, shifting every later row onto the next quantity's label)
    local r = `w' - `mint' - 1
    if `r' >= 0 local lbls `lbls' post`r'
    if `r' == 0 & `old' == 1 local lbls `lbls' contemp
    if `r' < 0 local lbls `lbls' pre`=abs(`r')'
    local f2lbls `f2lbls' `=cond(`r' < 0, "pre" + strofreal(abs(`r')), "post" + strofreal(`r'))'
    if "`relcol'" == "" local relcol `r'
    else local relcol `relcol' , `r'
  }
  matrix `rel' = (`relcol')'
  * the six-row table keeps the estimable columns only (column by column, as the in-line
  * code did, so an undropped call is the same slice as before)
  tempname kt
  local firstk 1
  foreach k of local keep {
    local c = `first' + `k' - 1
    if `firstk' {
      matrix `kt' = `tab'[1..6, `c']
      local firstk 0
    }
    else matrix `kt' = (`kt', `tab'[1..6, `c'])
  }
  matrix `tab' = `kt'
  if "`b'" != "" & `ndrop' > 0 {
    local keepm = subinstr(trim("`keep'"), " ", ",", .)
    mata: st_matrix("`b'", st_matrix("`b'")[1, (`keepm')])
    * fast mode hands b alone (no covariance is captured or posted there), so v() is optional
    if "`v'" != "" mata: st_matrix("`v'", st_matrix("`v'")[(`keepm'), (`keepm')])
  }
  c_local `lbl' `lbls'
  if "`f2l'" != "" c_local `f2l' `f2lbls'
end

* bug sweep BS-3 (2026-09-06): interactivecontrols() is silently inert when no cohort-period
* cell holds within-cell variation in a control - every treated-cell slope term
* i._i#i._t#1.targetflag#c.<ic> is dropped for collinearity and the run reduces to
* controls(i.t#c.x). twin of the contreat R6 screen in the worker: a named caution, since
* the generic collinearity caution cannot say which control or what the run reduced to.
* ics() are the demeaned copies in the user's order, so word k maps back to names()
cap prog drop wooldid_icscreen
prog def wooldid_icscreen
  syntax , bvars(string asis) iv(name) tv(name) tf(name) ics(string) names(string)
  local k 0
  foreach ic of local ics {
    local ++k
    local live 0
    local drop 0
    foreach tok of local bvars {
      if strpos("`tok'", "`iv'") != 0 & strpos("`tok'", "`tv'") != 0 & (strpos("`tok'", "1.`tf'") != 0 | strpos("`tok'", "1o.`tf'") != 0) & (strpos("`tok'#", "c.`ic'#") != 0 | strpos("`tok'#", "co.`ic'#") != 0) {
        if strpos("`tok'", "o.") == 0 local ++live
        else local ++drop
      }
    }
    if `live' == 0 & `drop' > 0 {
      local nm : word `k' of `names'
      dis as text "Caution: every treated cohort-period slope term for the interactive control `nm' was dropped for multicollinearity: `nm' does not vary within the cohort-period cells (for example i is as fine as the unit and `nm' is a unit-level variable, or `nm' is constant across units within a period), so the treatment effects cannot vary with it and this run is identical to controls(i.t#c.`nm'). Coarsen i so that cells hold variation in `nm', or move it to controls()."
    }
  }
end

* NFADV-01/02/05: the adversarial pre-trend machinery. One caller replaces the four
* near-identical frame/svmat/`test' blocks in wooldid_worker (per-program size ceiling).
* Inputs are the F2 event-study capture of the SAME margins run (b 1xk, V kxk full
* covariance), the aligned relative-time column relyr (kx1), min_t, esrelativeto, and
* the reference-regime/timetrends/altmethod flags. It does two things:
*  (1) fits the DISPLAY LINE by the same unweighted OLS `reg b t' the old in-worker
*      frames ran - identical rows (real horizon coefficients in horizon order, then
*      the planted zero-valued reference row(s)) and the identical internal regressor
*      t = relyr + min_t + 1 (the worker's LOAD-BEARING identity), so the published
*      slope is bit-identical to the old one - and re-expresses the intercept at
*      relative time 0 (NFADV-05: the old _cons was the fitted value at relyr =
*      -min_t-1), publishing columns [slope, intercept_rt0];
*  (2) computes a CALIBRATED p-value on the event-study covariance, esfixedbaseperiod
*      only (the pooled reference has no single omitted period; the worker prints one
*      Caution per run and the line is returned without inference). Primary = a
*      minimum-distance test of the composite null that all horizon effects lie on
*      SOME line (NFADV-01: the old code froze the fitted line as text constants in
*      `test', a weighted chi-square whose null size depends on the margins VCV).
*      Altmethod = the post coefficients tested against the pre-period-line forecasts
*      using the FULL forecast covariance S = V_QQ + A V_PP A' - A V_PQ - V_QP A'
*      (NFADV-02: the old code used V_QQ alone, anti-conservative). The statistics
*      live in mata (wooldid_advx, end of file); F vs chi-squared reference follows
*      the wooldid_ptp df_r convention.
* Returns r(status): "ok" (line + p), "nop" (line only), "gapskip" (nothing - the B16
* gap backstop: the altmethod fit sample has no real pre-period coefficient).
* The rest of the return contract, which the call sites branch on: r(p) exists only on
* "ok"; r(line) exists on "ok" and "nop" but not on "gapskip". "nop" means no p-value was
* computable because esfixedbaseperiod was not in force (a pooled reference has no single
* omitted period), and only the display line comes back. "gapskip" obliges the caller to
* CANCEL adversarial - at the overall call site that clears the option for the whole run,
* at the per-subgroup call sites it skips that subgroup only. The display line is fit over
* ALL horizons by default and over the pre-period horizons only under
* altmethod_adversarial (the empty-vs-nonempty `advcut' inside the program).
cap prog drop wooldid_advstat
prog def wooldid_advstat, rclass
version 16.0
  syntax , b(string) v(string) relyr(string) mint(string) esrel(string) [dfr(string) esf(string) tt(string) alt(string) sg(string) fm(string)]
  if "`fm'" == "" local fm ___wooldidfullmodel

  * B16 gap backstop (altmethod): the planted reference row(s) always satisfy the
  * pre-period cut, so require at least one REAL pre-period coefficient beyond them
  * (both regimes' cuts reduce to relyr <= -1: esf uses t <= min_t and the pooled
  * regime uses t <= espre with espre = min_t)
  if "`alt'" != "" {
    mata: st_local("advnpre", strofreal(sum(st_matrix("`relyr'") :<= -1)))
    if `advnpre' < 1 {
      if "`sg'" == "" dis as text "Caution: adversarial pre-trend testing was cancelled: after accounting for gaps in the data, too few pre-treatment coefficients remain to fit the adversarial reference line."
      else dis as text "Caution: the adversarial pre-trend test was skipped for subgroup `sg': after accounting for gaps in the data, too few pre-treatment coefficients remain to fit the reference line."
      return local status gapskip
      exit
    }
  }

  * the display line. rows are assembled exactly as the old adversarial frames were:
  * real coefficients first (b transposed against t = relyr + min_t + 1), then the
  * planted reference row(s) - esf: relyr = esrelativeto (preceded by esrelativeto - 1
  * under timetrends); pooled: one horizon below the earliest real horizon present.
  * missing coefficients become 0, as the old frame's `replace b = 0 if b == .' did.
  local k = rowsof(`relyr')
  tempname advF advB advL
  matrix `advF' = `b'' , (`relyr' + J(`k',1,`mint'+1))
  forval advi = 1/`k' {
    if `advF'[`advi',1] >= . matrix `advF'[`advi',1] = 0
  }
  if "`esf'" != "" {
    if "`tt'" != "" matrix `advF' = `advF' \ (0, `esrel'+`mint')
    matrix `advF' = `advF' \ (0, `esrel'+`mint'+1)
  }
  else {
    mata: st_local("advminr", strofreal(min(st_matrix("`relyr'"))))
    matrix `advF' = `advF' \ (0, `advminr'+`mint')
  }
  qui frame
  local advfr0 `r(currentframe)'
  tempname advfr
  frame create `advfr'
  frame change `advfr'
  qui svmat double `advF'
  local advcut
  if "`alt'" != "" local advcut if `advF'2 <= `mint'
  qui reg `advF'1 `advF'2 `advcut'
  matrix `advB' = e(b)
  frame change `advfr0'
  frame drop `advfr'
  * the little reg clobbered e(); put the main model back (stored at worker start,
  * under the replication name inside a bootstrap/permutation replication)
  qui est restore `fm'
  matrix `advL' = (`advB'[1,1], `advB'[1,2] + (`mint'+1)*`advB'[1,1])
  matrix colnames `advL' = slope intercept_rt0
  matrix rownames `advL' = y1

  * inference (esf only; on the F2 covariance, in mata)
  local advstatus nop
  if "`esf'" != "" {
    local advdfr .
    if "`dfr'" != "" & "`dfr'" != "." {
      if real("`dfr'") > 0 & real("`dfr'") < . local advdfr `dfr'
    }
    tempname advP
    local advpstat
    mata: wooldid_advx("`b'", "`v'", "`relyr'", `esrel', ("`tt'"!=""), ("`alt'"!=""), `advdfr', "`advP'")
    if "`advpstat'" == "ok" {
      return scalar p = scalar(`advP')
      local advstatus ok
    }
    else {
      * covers rank-0/degenerate covariance (e.g. timetrends), missing entries, a
      * residual outside the covariance's column space, and no testable residual
      * dimensions (e.g. only two real horizons under the primary composite null)
      if "`sg'" == "" dis as text "Caution: adversarial: the composite-null covariance is degenerate; no p-value computed."
      else dis as text "Caution: adversarial: the composite-null covariance is degenerate for subgroup `sg'; no p-value computed."
    }
  }
  return local status `advstatus'
  return matrix line = `advL'
end


* AMESE (standing-issues family A, R1/R4): screens for cluster-robust variance
* pathologies around the continuous-treatment (ame) machinery. Lives outside
* wooldid_worker for per-program size headroom. FOUR modes, in the order the program
* implements them - detect and vzero are described below; the other two are:
*   bycoh: the per-treated-cohort POINT-identification screen for
*     contreatcontrols(byCohort) (2026-08-25 Group B ruling). Returns r(fail); the
*     caller raises error 459. It runs REGARDLESS of the fastnose/impvalues/unconditionalse
*     gate below, because point identification rather than variance is at stake.
*   sescan: the belt-and-braces duplicate-standard-error signature. Returns r(ndup);
*     caution only. This is the mode the worker's ame-family call asks for.
*   detect: the STRUCTURAL predicate — with a cluster-type VCE and contreat() active,
*     if ANY treated cohort-by-period cell of the realized estimation sample lies
*     wholly inside a single cluster, the plug-in CRVE cannot see within-cell dose
*     variation: the ame "standard errors" it would deliver collapse toward the
*     variance of the absorbed baseline dose control (in the basic/byTime cases they
*     ARE that variance, identically across every ame quantity), so they do not
*     measure the sampling variability of the ame estimates. ANY-cell semantics is
*     load-bearing: partial containment (a cluster variable nesting only some cells)
*     still contaminates the V while presenting a full-rank matrix. Rank/eigenvalue
*     guards are deliberately NOT used (the byTime case posts a full-rank V that is
*     still 100 percent blind to within-cell noise). Returns r(fire), r(ncells) and
*     prints the explanatory caution once when firing.
*   vzero: detects the estimation VCE arriving as an all-zero matrix — Stata's
*     `ereturn post` screen zeroes a nonsymmetric/highly-singular VCE WHOLESALE
*     (destroying valid entries in transit) and the quiet estimation call swallows
*     the only warning. Returns r(vzero) and prints its caution when firing.
* The predicate is never evaluated on the fastnose path (no variances exist there),
* under impvalues (IM inference is cross-cluster dispersion — the estimator class
* that remains valid under cell saturation), or under unconditionalse (a different
* estimand whose margins variance carries live sampling-of-X components; excluded
* from suppression by design — see the R1 design note).
cap prog drop wooldid_amescreen
prog def wooldid_amescreen, rclass
version 16.0
  syntax , mode(string) [icell(string) tcell(string) tflag(string) esamp(string) clusvars(string) contreat(string) sgnote(string) semat(string) iorig(string) tvar(string) r1fired(string)]
  if "`mode'" == "detect" {
    local fire 0
    local ncells 0
    tempvar amecl ametag
    foreach cv of local clusvars {
      cap drop `amecl'
      qui gegen double `amecl' = nunique(`cv') if `esamp' & `tflag' == 1, by(`icell' `tcell')
      cap drop `ametag'
      qui gegen byte `ametag' = tag(`icell' `tcell') if `esamp' & `tflag' == 1 & `amecl' == 1
      qui count if `ametag' == 1
      if r(N) > 0 {
        local fire 1
        if r(N) > `ncells' local ncells = r(N)
      }
    }
    if `fire' {
      dis as text "Caution: with standard errors clustered on `clusvars', `ncells' treated i-t cell(s) lie completely within a single cluster, preventing valid analytic estimation of SEs for AME estimates (and with some adverse implications for the ATTs as well). Consider switching to use of bootstrap SEs, IM p-values, or randomization inference if finer clustering levels are non-viable."
    }
    return scalar fire = `fire'
    return scalar ncells = `ncells'
    exit
  }
  if "`mode'" == "bycoh" {
    * GROUP B (author-ruled 2026-08-25): per-treated-cohort POINT-identification
    * screen for contreatcontrols(byCohort). The cohort-specific baseline dose slope
    * is identified only by dose variation across the cohort's untreated/reference-
    * pool observations (the omitted base cell, tcell == 0; every other cell carries
    * its own dose term). Zero variation there -> the slope is unestimable and that
    * cohort's AME decomposition is solver-arbitrary -> hard refusal (caller raises
    * r459, in the style of the pre-treatment-period sufficiency errors). Variation
    * below 1/1000 of the sample dose sd -> near-unidentified direction where
    * estimation tolerances bind -> computational-fragility caution. Runs regardless
    * of VCE/fastnose/impvalues: point identification is at stake, not variances.
    * NOTE the cell coding: `_i' pools ALL tcell==0 rows (reference + never-treated)
    * into cohort code 0, so treated cohorts must be grouped by the ORIGINAL cohort
    * variable among treated rows (tvar = the working onset variable; treated rows
    * carry a nonmissing onset)
    tempvar bctr bcref bchas bcsd bcn
    qui gen byte `bctr' = (`tvar' < .) & `esamp'
    qui gen byte `bcref' = (`tcell' == 0) if `bctr'
    qui gegen byte `bchas' = max(`bcref') if `bctr', by(`iorig')
    qui gegen double `bcsd' = sd(`contreat') if `bctr' & `tcell' == 0, by(`iorig')
    qui gegen double `bcn' = count(`contreat') if `bctr' & `tcell' == 0, by(`iorig')
    local fail 0
    tempname bcmin bcsamp
    scalar `bcmin' = .
    qui count if `bctr' & `bchas' == 0
    if r(N) > 0 {
      qui gstats sum `iorig' if `bctr' & `bchas' == 0, meanonly
      scalar `bcmin' = r(min)
    }
    qui count if `bctr' & `bcref' == 1 & (`bcn' < 2 | `bcsd' == 0 | `bcsd' >= .)
    if r(N) > 0 {
      qui gstats sum `iorig' if `bctr' & `bcref' == 1 & (`bcn' < 2 | `bcsd' == 0 | `bcsd' >= .), meanonly
      if r(min) < `bcmin' scalar `bcmin' = r(min)
    }
    if `bcmin' < . {
      local fail 1
      dis as error "Error: the continuous treatment variable `contreat' has no variation across the untreated/reference-period observations of treated cohort `=`bcmin'', so the cohort-specific baseline dose slope required by contreatcontrols(byCohort) cannot be estimated for that cohort."
      dis as error "Examine this cohort's `contreat' values in its untreated/reference periods, and consider dropping the cohort from the estimation sample before proceeding."
    }
    else {
      qui gstats sum `contreat' if `esamp'
      scalar `bcsamp' = r(sd)
      qui count if `bctr' & `bcref' == 1 & `bcsd' > 0 & `bcsd' < 1e-3*`bcsamp' & `bcsd' < .
      if r(N) > 0 {
        tempvar bctag
        qui gegen byte `bctag' = tag(`iorig') if `bctr' & `bcref' == 1 & `bcsd' > 0 & `bcsd' < 1e-3*`bcsamp' & `bcsd' < .
        qui count if `bctag' == 1
        local bck = r(N)
        qui gstats sum `iorig' if `bctr' & `bcref' == 1 & `bcsd' > 0 & `bcsd' < 1e-3*`bcsamp' & `bcsd' < ., meanonly
        dis as text "Caution: under contreatcontrols(byCohort), `bck' treated cohort(s) (smallest: `=r(min)') have untreated/reference-period variation in `contreat' below 1/1000 of its sample standard deviation. The cohort-specific baseline dose slopes are nearly unidentified there, so AME estimation is likely computationally fragile: consider increasing regtol(), or examine and drop the affected cohort(s)."
      }
    }
    return scalar fail = `fail'
    exit
  }
  if "`mode'" == "sescan" {
    * generalizable belt-and-braces guard: >=3 reported standard errors that are
    * pairwise identical to ~14 digits under a cluster-type VCE is the signature of a
    * degenerate variance being mis-attributed across quantities (legitimately
    * identical SEs arise on balanced UNCLUSTERED designs, hence the caller's gate).
    * Caution only - the structural predicate is the enforcement layer.
    tempname sesc
    matrix `sesc' = `semat'
    local sn = rowsof(`sesc')
    local sdup 0
    if `sn' >= 3 {
      forval sa = 1/`=`sn'-1' {
        local sva = `sesc'[`sa', 2]
        local sba = `sesc'[`sa', 1]
        if `sva' >= . | `sva' == 0 continue
        local smatch 1
        forval sb = `=`sa'+1'/`sn' {
          local svb = `sesc'[`sb', 2]
          local sbb = `sesc'[`sb', 1]
          * identical SE with a DIFFERING point estimate is the mis-attribution
          * signature; identical (b, se) pairs are legitimately coincident estimands
          if `svb' < . & `svb' != 0 & reldif(`sva', `svb') <= 1e-12 & reldif(`sba', `sbb') > 1e-9 local ++smatch
        }
        if `smatch' > `sdup' local sdup = `smatch'
      }
    }
    if `sdup' >= 3 {
      dis as text "Caution: `sdup' reported standard errors appear to be identical under a cluster-type VCE. This usually indicates a failure of analytic clustered standard error estimation; consider using a different approach to SE estimation (or consider clustering at a finer level of aggregation)."
    }
    return scalar ndup = `sdup'
    exit
  }
  if "`mode'" == "vzero" {
    local vzero 0
    tempname vchk
    cap matrix `vchk' = e(V)
    if _rc == 0 {
      mata: st_local("vzero", strofreal(max(abs(st_matrix(st_local("vchk")))) == 0))
    }
    * the caution defers to the more precise cell-containment caution when that one
    * already displayed (r1fired); detection is still returned for the refit logic
    if `vzero' & "`r1fired'" != "1" {
      if "`contreat'" != "" dis as text "Caution: estimation returned an all-zero variance-covariance matrix. The issue is likely related to continuous treatment effect estimation coupled with use of analytic clustered SEs, with clusters occuring at a high level of aggregation. Consider switching to use of bootstrap SEs, IM p-values, or randomization inference if finer clustering levels are non-viable."
      else dis as text "Caution: estimation returned an all-zero variance-covariance matrix. This usually reflects a degenerate clustering (very few clusters, or clusters at a high level of aggregation relative to the design). Consider switching to bootstrap SEs, IM p-values, or randomization inference."
    }
    return scalar vzero = `vzero'
    exit
  }
  dis as error "wooldid_amescreen: unknown mode `mode'"
  error 198
end


* F11: event-study plots, moved out of wooldid_worker (per-program size ceiling). reads
* e() plus the two per-call options only; runs inside its own tempname frame and returns
* to the caller's frame on the normal path. on an error inside this program, stata
* auto-drops the tempname frame and the wooldid shell's rc-gated epilogue restores the
* caller's frame - the same net behavior as the pre-F11 inline code.
cap prog drop wooldid_cfxplots
prog def wooldid_cfxplots
version 16.0
  * continuous-treatment (cfx) plots drawn from the POSTED e(cfx*) matrices - under
  * bootstrap() those are the rewired ones, so the bands are the displayed family's.
  * A twoway rendering (connected estimate, rarea band) of what marginsplot drew from the
  * margins results in an analytic run: same names, same titles, same saveplots() files,
  * plus one <saveplots>_cfx.dta holding every plotted point. Rows without both bounds
  * (cfxplotnose, fast mode, or a studentized cell with no analytic se) plot as points only
  * P7 accident 2 (author ruling 2026-09-08): dataonly builds the plot dataset and writes
  * <saveplots>_cfx.dta WITHOUT drawing or saving any graph. It exists for the ANALYTIC path,
  * where marginsplot has already drawn the figures and saved the .gph files, and only the
  * plot-data file was missing - an engine-only artifact for no reason anyone could state.
  syntax [, SAVEplots(name) replace dataonly]
  if "`e(wdid_replayv)'" != "2" {
    dis as error "Error: wooldid_cfxplots is internal to wooldid and needs a stored wooldid estimation in e(); run wooldid first."
    error 301
  }
  local mats `e(wdid_cfxmats)'
  if "`mats'" == "" exit 0
  qui frame
  local defaultframe `r(currentframe)'
  * the x-axis labels: the lattice variable's value labels when it has them (marginsplot
  * shows those), read here in the caller's frame before we leave it
  local latvar `e(wdid_cfxlatvar)'
  local vl
  if "`latvar'" != "" {
    cap confirm variable `latvar'
    if _rc == 0 local vl : value label `latvar'
  }
  local levels
  foreach nm of local mats {
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    foreach rt in `: rownames e(`nm')' {
      * regexs() must be read on a LATER line: a `=regexs(1)' on the regexm line itself is
      * macro-expanded before the match runs and returns the previous match
      if regexm("`rt'", "^lattice_([0-9]+)$") {
        local lv = regexs(1)
        local levels `levels' `lv'
      }
    }
  }
  local levels : list uniq levels
  local levels : list sort levels
  local xlab
  if "`levels'" != "" {
    local xlab xlabel(
    foreach v of local levels {
      if "`vl'" != "" {
        local lbl : label `vl' `v', strict
        if `"`lbl'"' == "" local lbl `v'
        local xlab `xlab' `v' `"`lbl'"'
      }
      else local xlab `xlab' `v'
    }
    local xlab `xlab')
  }
  local xt `"`e(wdid_cfxxt)'"'
  tempname pf af
  frame create `pf'
  frame create `af'
  frame change `pf'
  local k 0
  foreach nm of local mats {
    * `local ++k' MUST stay above the `cap confirm' / `continue' pair below. k indexes the
    * per-plot metadata e(wdid_cfxttl`k')/e(wdid_cfxytx`k'), which the worker posted BY
    * POSITION in e(wdid_cfxmats) (pairs laid down ame first, then att), so a matrix that
    * fails the confirm must still consume its index. Move the increment below the continue
    * - the obvious tidy-up - and every later figure inherits the previous entry's title and
    * y-title: a plot of the average marginal effects titled as the ATT, on a figure the
    * user will publish
    local ++k
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    tempname M
    matrix `M' = e(`nm')
    local rn : rownames `M'
    local nr = rowsof(`M')
    local nc = colsof(`M')
    clear
    qui set obs `nr'
    qui gen str32 matname = "`nm'"
    qui gen double lattice = .
    qui gen double b = .
    qui gen double se = .
    qui gen double p = .
    qui gen double ll = .
    qui gen double ul = .
    forval r = 1/`nr' {
      local rt : word `r' of `rn'
      if regexm("`rt'", "^lattice_([0-9]+)$") {
        local lv = regexs(1)
        qui replace lattice = `lv' in `r'
      }
      qui replace b = `M'[`r', 1] in `r'
      * the cfx matrices are transposed margins r(table) blocks, so the columns are
      * b(1) se(2) t(3) pvalue(4) ll(5) ul(6) - column 3 is skipped here because the plot
      * has no use for the t. This is a WHOLE-MATRIX shape test, not a per-row condition:
      * the worker sets cfxmatendpoint to 1 under cfxplotnose or fastnose, and then the
      * matrix is one column wide and se/p/ll/ul are missing for every row. The header's
      * "rows without both bounds plot as points only" describes the per-row path and does
      * not cover this one
      if `nc' >= 6 {
        qui replace se = `M'[`r', 2] in `r'
        qui replace p = `M'[`r', 4] in `r'
        qui replace ll = `M'[`r', 5] in `r'
        qui replace ul = `M'[`r', 6] in `r'
      }
    }
    * LOAD-BEARING, not tidying: twoway's connected/line/rarea join points in DATASET
    * ORDER, not in x order (see the same sort before each event-study twoway)
    sort lattice
    qui count if ll < . & ul < . & b < .
    local band = r(N) > 0
    local ttl `"`e(wdid_cfxttl`k')'"'
    local ytx `"`e(wdid_cfxytx`k')'"'
    if "`dataonly'" == "" {
      if `band' twoway (rarea ll ul lattice, color(navy%50)) (connected b lattice, color(navy)), `xt' `xlab' ytitle(`"`ytx'"') title(`"`ttl'"') scale(.85) legend(off) name("`nm'", replace)
      else twoway (connected b lattice, color(navy)), `xt' `xlab' ytitle(`"`ytx'"') title(`"`ttl'"') scale(.85) legend(off) name("`nm'", replace)
    }
    if "`saveplots'" != "" {
      if "`dataonly'" == "" qui graph save "`saveplots'_`nm'.gph", `replace'
      tempfile part
      qui save `part'
      frame `af' {
      * the `cap' here has no stated cause and no realistic one: every part file is built
      * by the same seven gen statements above (str32 matname plus six doubles), so the
      * parts are structurally identical. If an append ever did fail, this whole plot's rows
      * would be missing from <saveplots>_cfx.dta with no message and no rc - which is the
      * silent no-op the shell's analytic twin of this write explicitly refuses (it uses
      * `cap noi' and prints a caution). Flagged rather than changed: bringing it under the
      * same rule is a behaviour change
        cap append using `part'
      }
    }
  }
  if "`saveplots'" != "" {
    frame `af' {
      qui save "`saveplots'_cfx.dta", `replace'
    }
  }
  frame change `defaultframe'
  frame drop `pf'
  frame drop `af'
end

cap prog drop wooldid_esplots
prog def wooldid_esplots
version 16.0
  syntax [, SAVEplots(name) replace]
  * guard the hand-call surface: without a posted wooldid result, the body would die
  * later inside its own tempname frame with a misleading "relyr not found" (an absent
  * e() matrix silently reads as a 1x1 missing scalar), stranding a caller who typed
  * this directly from a non-default frame in the default frame
  * exact-match on the replay contract version, mirroring the wooldid_display gate
  if "`e(wdid_replayv)'" != "2" {
    dis as error "Error: wooldid_esplots is internal to wooldid and needs a stored wooldid estimation in e() (e(wdid_replayv)='`e(wdid_replayv)'', this version requires '2'); run wooldid first."
    error 301
  }
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    dis as error "Error: the stored results no longer contain e(wdidmainresults); rerun wooldid."
    error 301
  }
  qui frame
  local defaultframe `r(currentframe)'

  * rehydrate what the plot code reads from e()
  local y `e(depvar)'
  local esrelativeto `e(wdid_esrel)'
  local timetrends `e(wdid_timetrends)'
  local semielasticity `e(wdid_semielas)'
  local contreatelasticitytype `e(wdid_cetype)'
  * BS-6 (2026-09-07): the two macros above hold the OPERATOR margins used. The headers must
  * name what the user asked for - which the remapped quantity now genuinely is - so prefer
  * the raw request where it was recorded. Absent on every run where no remap fired and on
  * pre-BS-6 stores, so replay of an older estimate is unchanged.
  if "`e(wdid_semielasraw)'" != "" local semielasticity `e(wdid_semielasraw)'
  if "`e(wdid_cetyperaw)'" != "" local contreatelasticitytype `e(wdid_cetyperaw)'
  local subgroups `e(wdid_sgvar)'
  local sglist `e(wdid_sglist)'
  local suppresssgprimaryeffects `e(wdid_suppresssgpe)'
  local esfixedbaseperiod
  if "`e(estimate_type)'" == "fixed_reference_period" local esfixedbaseperiod esfixedbaseperiod
  tempname wdidpooled
  matrix `wdidpooled' = e(wdidmainresults)



    * this creates and enters a frame where we will produce graphs using wdidpooled and similar for sgs
    tempname __wooldidframe
    frame create `__wooldidframe'
    frame change `__wooldidframe'
    qui svmat double `wdidpooled', names(col)
    * the interval columns are read by NAME below; a bootstrap or permutation run at a
    * level other than 95 renames them lb_<level>ci/ub_<level>ci, so map them back
    * NOTE the `if _rc != 0' test below is separated from this `confirm' by two `local'
    * statements and is still correct: only `capture' sets _rc, and a plain local assignment
    * leaves it untouched. The surrounding convention is to stash _rc into a named local the
    * moment it is produced (wooldid_esband does), and a reader applying the codebase's own
    * "_rc is reset by the very next command" rule will read this as a live bug - it is not
    cap confirm variable lb_95ci
    local wdes_lbname
    local wdes_ubname
    if _rc != 0 {
      cap ds lb_*ci
      if _rc == 0 & wordcount("`r(varlist)'") == 1 {
        local wdes_lbname `r(varlist)'
        rename `r(varlist)' lb_95ci
      }
      cap ds ub_*ci
      if _rc == 0 & wordcount("`r(varlist)'") == 1 {
        local wdes_ubname `r(varlist)'
        rename `r(varlist)' ub_95ci
      }
    }
    qui gen rownamesvar = ""
    local names : roweq `wdidpooled'
    forvalues i=1/`: word count `names'' {
      qui replace rownamesvar=`"`: word `i' of `names''"' in `i'
    }
    * THE FIGURE LIST IS DECIDED HERE, and it must be decided BEFORE the reference row is
    * planted below. rownamesvar holds the row equations of e(wdidmainresults), whose ES_
    * families are exactly ES_att, ES_att_it, ES_ame, ES_ame_it; non-event-study rows carry
    * a missing relyr and are what this drop removes. The planted reference observation is
    * given rownamesvar == "reference" precisely so the restriction inside each twoway puts
    * the shared zero on every family's figure - plant it above this glevelsof and
    * "reference" becomes a graph type in its own right, and wooldid draws a fifth
    * one-point figure and saves it under the user's saveplots stem
    qui drop if relyr == .
    qui glevelsof(rownamesvar), clean local(graphtypes)
    qui gen postperiod =  relyr >= 0
    keep   postperiod relyr  estimate lb_95ci ub_95ci rownamesvar

    * add the pre treatment observation / reference period
    if "`esfixedbaseperiod'" == "" {
      local newN = _N + 1
      qui gstats sum relyr, meanonly
      local max = r(max)
      local min = r(min)
      local initialflagt = r(min) - 1
      qui set obs `newN'
      qui replace relyr = `initialflagt' if _n == _N
      qui replace lb_95ci = 0 if _n == _N
      qui replace ub_95ci = 0 if _n == _N
      qui replace estimate = 0 if _n == _N
      qui replace rownamesvar = "reference" if _n == _N
      * TWO IDENTIFICATION RESTRICTIONS ARE BEING DRAWN AS DATA. Under the pooled
      * reference there is no single omitted period at all - the comparison period is a
      * cohort-specific pool of pre-periods, which is why pretrendpower() is refused there -
      * so this branch parks one zero at an invented position (min - 1) and labels it with a
      * bare asterisk rather than a number. The esfixedbaseperiod branch below plants ONE
      * zero normally and TWO under timetrends (at esrelativeto and esrelativeto - 1),
      * because a linear trend costs a second normalization: those horizons are zero by
      * construction, not by estimation. Without this, a reader takes the zeros for
      * estimates that happened to come out at zero with zero width.
      * SYNTAX NOTE: `xlabel(' is opened on this line and not closed until the `local
      * xlabelstarter `xlabelstarter' )' several statements below - the indented forval
      * between them is a sequential append to the same macro, not a nested block. Stata
      * locals store raw text, so the unbalanced "(" is legal. Reindent, reorder or return
      * early out of that stretch and the twoway call ships with an unterminated xlabel(
      * and an r(198) that points nowhere near the edit
      local xlabelstarter xlabel(`initialflagt' "*"
        forval xlabeltick = `min'(1)`max' {
          local xlabelstarter `xlabelstarter' `xlabeltick' "`xlabeltick'"
        }
        local xlabelstarter `xlabelstarter' )

    }
    else{
      qui   gstats sum relyr, meanonly
      local max = r(max)
      local min = r(min)
      local newN = _N + 1
			if "`timetrends'" == "timetrends" local ++newN
      qui set obs `newN'
			qui replace relyr = `esrelativeto' if _n == _N
			if "`timetrends'" == "timetrends" qui replace relyr = `esrelativeto' - 1 if _n == _N - 1
      qui replace lb_95ci = 0 if _n == _N
      qui replace ub_95ci = 0 if _n == _N
      qui replace estimate = 0 if _n == _N
			qui replace rownamesvar = "reference" if _n == _N
			if "`timetrends'" == "timetrends" qui replace lb_95ci = 0 if _n == _N -1
			if "`timetrends'" == "timetrends" qui replace ub_95ci = 0 if _n == _N -1
			if "`timetrends'" == "timetrends" qui replace estimate = 0 if _n == _N -1
			if "`timetrends'" == "timetrends" qui replace rownamesvar = "reference" if _n == _N - 1
			if `esrelativeto' < `min' local min `esrelativeto'
			if "`timetrends'" == "timetrends" if `esrelativeto' -1  < `min' local min = `esrelativeto' - 1
      local xlabelstarter xlabel(
        forval xlabeltick = `min'(1)`max' {
          local asterisk
					if `xlabeltick' == `esrelativeto' local asterisk *
					if "`timetrends'" == "timetrends" & `xlabeltick' == `esrelativeto' -1 local asterisk *
          local xlabelstarter `xlabelstarter' `xlabeltick' "`xlabeltick'`asterisk'"
        }
        local xlabelstarter `xlabelstarter' )
    }
    sort  relyr


    local xlineinserter
    qui   gstats sum relyr, meanonly
    if r(min) < 0 local xlineinserter xline(-0.5,   lpattern(dash) lcolor(black))

    * inject elasticity info
    local elasinj Effects
    if "`semielasticity'" == "eydx" local elasinj "Semi-elasticities"
    local elasinj2 Effects
    if "`contreatelasticitytype'" == "eyex" local elasinj2 Elasticities
    if "`contreatelasticitytype'" == "eydx" local elasinj2 Semi-elasticities (eydx)
    if "`contreatelasticitytype'" == "dyex" local elasinj2 Semi-elasticities (dyex)

    * generate and save the actual graphs; then save the data
    foreach graph in `graphtypes' {
      local graphtext = subinstr("`graph'","ES_","",.)
      local att_it_text ATT_it `elasinj'
      local att_text ATT `elasinj'
      local ame_it_text AME_it `elasinj2'
      local ame_text AME `elasinj2'
      * PS-45 (author ruling 2026-09-07): the simultaneous band, drawn as a dashed pair
      * over the pointwise shading when simultaneousbands stored one for THIS family.
      * `wdsblayer' is empty otherwise, so the command is unchanged on every other run
      local wdsblayer
      wooldid_esbandv , fam(`graphtext')
      twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr,color(navy%50)) `wdsblayer' if (rownamesvar =="reference" | rownamesvar=="`graph'")  , `xlineinserter' xtitle("Time Relative to Time of Initial Treatment") ytitle("``graphtext'_text' on `y' by Horizon") legend(off)   name(`graphtext'_fullmodel, replace) `xlabelstarter'
      if "`saveplots'" != "" qui graph save "`saveplots'_`graph'.gph", `replace'
    }
    * the band columns are working variables: the endpoints are exported through
    * e(esband_*), and dropping them here keeps saveplots' saved data as it was
    cap drop __wdsblo
    cap drop __wdsbhi


    if "`saveplots'" != "" {
      * the saved data keep the table's own interval names (lb_<level>ci): the lb_95ci
      * alias above served the twoway calls only (saved as lb_95ci at every level until
      * 2026-09-05)
      if "`wdes_lbname'" != "" rename lb_95ci `wdes_lbname'
      if "`wdes_ubname'" != "" rename ub_95ci `wdes_ubname'
      qui save "`saveplots'_ES.dta", `replace'
    }
    clear

    * repeat by SG if appropriate
    if "`subgroups'" != "" & "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        tempname wdidsgmat
        matrix `wdidsgmat' = e(wdidsgresults_`sgcase')
        qui svmat double `wdidsgmat', names(col)
        cap confirm variable lb_95ci
        local wdes_lbname
        local wdes_ubname
        if _rc != 0 {
          cap ds lb_*ci
          if _rc == 0 & wordcount("`r(varlist)'") == 1 {
            local wdes_lbname `r(varlist)'
            rename `r(varlist)' lb_95ci
          }
          cap ds ub_*ci
          if _rc == 0 & wordcount("`r(varlist)'") == 1 {
            local wdes_ubname `r(varlist)'
            rename `r(varlist)' ub_95ci
          }
        }
        qui gen rownamesvar = ""
        local names : roweq `wdidsgmat'
        forvalues i=1/`: word count `names'' {
          qui replace rownamesvar=`"`: word `i' of `names''"' in `i'
        }
        qui drop if relyr == .
        qui glevelsof(rownamesvar), clean local(graphtypes)
        qui gen postperiod =  relyr >= 0
        keep  rownamesvar postperiod relyr  estimate lb_95ci ub_95ci

        * add the pre treatment observation
        if "`esfixedbaseperiod'" == "" {
          local newN = _N + 1
          qui   gstats sum relyr, meanonly
          local max = r(max)
          local min = r(min)
          local initialflagt = r(min) - 1
          qui set obs `newN'
          qui replace relyr = `initialflagt' if _n == _N
          qui replace lb_95ci = 0 if _n == _N
          qui replace ub_95ci = 0 if _n == _N
          qui replace estimate = 0 if _n == _N
          qui replace rownamesvar = "reference" if _n == _N

          local xlabelstarter xlabel(`initialflagt' "*"
            forval xlabeltick = `min'(1)`max' {
              local xlabelstarter `xlabelstarter' `xlabeltick' "`xlabeltick'"
            }
            local xlabelstarter `xlabelstarter' )
        }
        else{
          qui   gstats sum relyr, meanonly
          local max = r(max)
          local min = r(min)
          local newN = _N + 1
					if "`timetrends'" == "timetrends" local ++newN
          qui set obs `newN'
          qui replace relyr = `esrelativeto' if _n == _N
					if "`timetrends'" == "timetrends" qui replace relyr = `esrelativeto' - 1 if _n == _N - 1
					qui replace lb_95ci = 0 if _n == _N
          qui replace ub_95ci = 0 if _n == _N
          qui replace estimate = 0 if _n == _N
					qui replace rownamesvar = "reference" if _n == _N
					if "`timetrends'" == "timetrends" qui replace lb_95ci = 0 if _n == _N -1
          if "`timetrends'" == "timetrends" qui replace ub_95ci = 0 if _n == _N -1
          if "`timetrends'" == "timetrends" qui replace estimate = 0 if _n == _N -1
					if "`timetrends'" == "timetrends" qui replace rownamesvar = "reference" if _n == _N - 1

          if `esrelativeto' < `min' local min `esrelativeto'
          if "`timetrends'" == "timetrends" if `esrelativeto' -1  < `min' local min = `esrelativeto' - 1
          local xlabelstarter xlabel(
            forval xlabeltick = `min'(1)`max' {
              local asterisk
							if `xlabeltick' == `esrelativeto' local asterisk *
							if "`timetrends'" == "timetrends" & `xlabeltick' == `esrelativeto' -1 local asterisk *
              local xlabelstarter `xlabelstarter' `xlabeltick' "`xlabeltick'`asterisk'"
            }
            local xlabelstarter `xlabelstarter' )
        }

        * LOAD-BEARING, not tidying: twoway's connected/line/rarea join points in DATASET
        * ORDER, not in x order. The reference observation was appended as the LAST
        * observation but sits at the far LEFT of the axis, so without this sort the line
        * runs left to right across the estimates and then doubles back to the reference
        * point, and the rarea shading self-intersects - a wrong-looking, error-free figure
        sort relyr
        foreach graph in `graphtypes' {
          local graphtext = subinstr("`graph'","ES_","",.)
          local att_it_text ATT_it `elasinj'
          local att_text ATT `elasinj'
          local ame_it_text AME_it `elasinj2'
          local ame_text AME `elasinj2'
          * PS-45 (author ruling 2026-09-07): the subgroup family's own simultaneous band
          local wdsblayer
          wooldid_esbandv , fam(`graphtext') sfx(_s`sgcase')
          twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr,color(navy%50)) `wdsblayer' if (rownamesvar =="reference" | rownamesvar=="`graph'")  , `xlineinserter' xtitle("Time Relative to Time of Initial Treatment") ytitle("``graphtext'_text' on `y' by Horizon") legend(off)   name(`graphtext'_sg`sgcase', replace) `xlabelstarter'   title("`subgroups' == `sgcase'")
          if "`saveplots'" != "" qui graph save "`saveplots'_`graph'_sg`sgcase'.gph", `replace'
        }
        cap drop __wdsblo
        cap drop __wdsbhi

        if "`saveplots'" != "" {
          if "`wdes_lbname'" != "" rename lb_95ci `wdes_lbname'
          if "`wdes_ubname'" != "" rename ub_95ci `wdes_ubname'
          qui save "`saveplots'_ES_`subgroups'_`sgcase'.dta", `replace'
        }
        clear
      }
    }

    frame change `defaultframe'
    frame drop `__wooldidframe'
end


*** F3: mata worker for the gteffects extraction (name parse + table assembly; kept out
*** of wooldid_worker to conserve compiled-size headroom). the cap mata drop makes the
*** file safely re-loadable within a session (two-build harnesses, dev reloads)
cap mata: mata drop wooldid_gtx()
mata:
void wooldid_gtx(string scalar gtiv, string scalar gttv, string scalar gttf,
                 string scalar rawname, string scalar sedname, real scalar gtmint,
                 real scalar gthstart, real scalar gthend, string scalar gtdfr,
                 string scalar outname, | real scalar gtlvl)
{
    real matrix cells, raw, out
    real matrix b, sed
    string matrix names
    string scalar pat
    real scalar j, r, q, df, g, tt, est, v, se, z, p, lb, ub, crit, n
    names = st_matrixcolstripe("e(b)")[., 2]
    b = st_matrix("e(b)")
    sed = st_matrix(sedname)
    raw = st_matrix(rawname)
    pat = "^([0-9]+)\." + gtiv + "#([0-9]+)\." + gttv + "#1\." + gttf + "$"
    cells = J(0, 3, .)
    for (j = 1; j <= rows(names); j++) {
        if (regexm(names[j], pat)) cells = cells \ (strtoreal(regexs(1)), strtoreal(regexs(2)), j)
    }
    n = rows(cells)
    st_local("gtn", strofreal(n))
    if (n == 0) return
    df = strtoreal(gtdfr)
    // F054: the interval level (default 95 when the caller passes none)
    gtal = (args() >= 11 ? (100 - gtlvl)/200 : .025)
    gtlt = subinstr(strofreal(100 - 200*gtal, "%9.0g"), ".", "p", .)
    out = J(n, 15, .)
    for (r = 1; r <= n; r++) {
        g = cells[r, 1]
        tt = cells[r, 2]
        j = cells[r, 3]
        est = b[1, j]
        v = sed[1, j]
        se = (v < . & v > 0) ? sqrt(v) : .
        p = lb = ub = .
        if (se < .) {
            z = est / se
            if (df < .) {
                p = 2 * ttail(df, abs(z))
                crit = invttail(df, gtal)
            }
            else {
                p = 2 * normal(-abs(z))
                crit = invnormal(1 - gtal)
            }
            lb = est - crit * se
            ub = est + crit * se
        }
        out[r, 1] = g
        out[r, 4] = tt - gtmint - 1
        out[r, 5] = est
        out[r, 6] = se
        out[r, 7] = p
        out[r, 8] = lb
        out[r, 9] = ub
        out[r, 10] = (tt >= gthstart & tt <= gthend)
        for (q = 1; q <= rows(raw); q++) {
            if (raw[q, 1] == g & raw[q, 2] == tt) {
                out[r, 2] = raw[q, 3]
                out[r, 3] = raw[q, 4]
                out[r, 11] = raw[q, 5]
                out[r, 12] = raw[q, 6]
                out[r, 13] = raw[q, 7]
                out[r, 14] = raw[q, 8]
                out[r, 15] = raw[q, 9]
                break
            }
        }
    }
    out = sort(out, (1, 4))
    st_matrix(outname, out)
    st_matrixcolstripe(outname, (J(15, 1, ""), ("i", "onset", "caltime", "relyr", "estimate", "se", "p", "lb_" + gtlt + "ci", "ub_" + gtlt + "ci", "inwindow", "nobs", "w_att", "w_att_it", "w_att_i", "w_att_itime")'))
}
end


*** F-collapse degeneracy predicate (see wooldid_fcdeg) on the r() of the margins call
*** just issued: r(b)/r(V) hold the level margins, whose SECOND half of columns are the
*** 1.targetflag margins the worker reads (the first half is the 0.targetflag base, always
*** omitted); r(V_vs), present after pwcompare, holds the pairwise contrasts likewise. At
*** relative tolerance tol the call is flagged when a second-half name carries an omission
*** marker; a second-half level variance is missing, all are zero, or one is <= tol x the
*** largest; the adjacent-difference contrast VCE of the second-half levels (the object
*** the subgroup-equality chain -test- and the overjoint pre-trend tests act on; skipped
*** for the event-study levels under `levels': those levels ARE contrasted - by the
*** adversarial line test and by whatever consumes e(V_es_*) - but both make their own
*** rank decision at 1e-12 of the largest eigenvalue, three orders tighter than tol, so
*** a residue-level direction is dropped identically whether the call ran collapsed or
*** on the full data, while a 1e-9 contrast test here would fire on every run whose
*** cluster count is below its horizon count, which is healthy and common (measured
*** 2026-09-06: 4/5/6 clusters against 8 horizons; bug sweep BS-29)) has its largest
*** eigenvalue <= tol x the largest level variance or its smallest <= tol x the larger of
*** its largest and the largest level variance; or a second-half pairwise variance is
*** missing or <= tol x the largest level variance.
*** Residue sits near 1e-16 and the smallest real subgroup contrast measured in the
*** assessment at 3.5e-5, so tol = 1e-9 has about seven orders of headroom either side.
*** A nose call carries no r(V), and margins writes no omission marker under nose, so
*** there is nothing to screen: not flagged. Sets locals fcdeg (0/1) and fcdegwhy.
cap mata: mata drop wooldid_fcdegx()
mata:
void wooldid_fcdegx(real scalar tol, real scalar levels)
{
    real matrix V, C, L
    real colvector d, dv
    real rowvector ev
    string matrix cs
    real scalar K, j0, n, k, dmax
    string scalar why
    st_local("fcdeg", "0")
    st_local("fcdegwhy", "")
    V = st_matrix("r(V)")
    K = cols(V)
    if (K == 0) return
    j0 = floor(K/2) + 1
    n = K - j0 + 1
    // the marker lives in the r(b)/r(table) names (what the read-1 parse uses), never
    // in r(V)'s - 0 of 148 screened calls carried it there in the 2026-09-01 probe
    cs = st_matrixcolstripe("r(b)")
    why = ""
    for (k=j0; k<=K; k++) {
        if (strpos(cs[k,2], "o.")) {
            why = "a subgroup margin that margins marks as not estimable"
            break
        }
    }
    d = diagonal(V)[|j0 \ K|]
    dmax = (hasmissing(d) ? . : max(d))
    if (why == "" & (dmax == . | dmax <= 0)) why = "subgroup margins with zero variance"
    if (why == "" & min(d) <= tol*dmax) why = "a subgroup margin whose variance is below " + strofreal(tol) + " of the largest one"
    if (why == "" & n >= 2 & !levels) {
        L = (I(n-1), J(n-1, 1, 0)) - (J(n-1, 1, 0), I(n-1))
        C = L * V[|j0, j0 \ K, K|] * L'
        C = (C + C') / 2
        ev = symeigenvalues(C)
        // the singular-direction floor is anchored to the level scale as well as to the
        // contrast scale: a rank-1 contrast VCE whose largest eigenvalue is itself ~1e-8
        // of dmax passed the contrast-relative test alone (2026-09-01 probe, b_q1e3)
        if (max(ev) <= tol*dmax) why = "subgroup contrasts with zero variance"
        else if (min(ev) <= tol*max((dmax, max(ev)))) why = "a contrast among the subgroup margins lying in a numerically singular direction of the VCE"
    }
    if (why == "") {
        V = st_matrix("r(V_vs)")
        K = cols(V)
        if (K > 0) {
            j0 = floor(K/2) + 1
            dv = diagonal(V)[|j0 \ K|]
            if (hasmissing(dv)) why = "a pairwise subgroup contrast with missing variance"
            else if (min(dv) <= tol*dmax) why = "a pairwise subgroup contrast whose variance is below " + strofreal(tol) + " of the largest subgroup margin variance"
        }
    }
    if (why == "") return
    st_local("fcdeg", "1")
    st_local("fcdegwhy", why)
}
end


*** F6: assignment generation for permutation inference. fills <assignname> (C x R,
*** frame row order) and c_locals: wdp_R, wdp_npd/wdp_nps (possible-assignment count,
*** double and string), wdp_caut0/1/2 (cautions), wdp_err (fatal message, caller raises)
cap mata: mata drop wooldid_permx()
cap mata: mata drop wooldid_permx_enum()
cap mata: mata drop wooldid_permx_msp()
mata:
void wooldid_permx(string scalar fr, string scalar onsv, string scalar strv,
                   real scalar reps, real scalar exact, string scalar assignname,
                   string scalar pname, string scalar lname)
{
    real matrix A
    real colvector ons, str, idx, levels, v
    real scalar C, si, lnM, M, R, r, k, capexact, capwarn, oc
    string scalar nps
    st_local("wdp_err", "")
    st_local("wdp_caut1", "")
    st_local("wdp_caut2", "")
    // BS-23 (author ruling 2026-09-07): the accumulator ROW that holds the observed
    // assignment. 0 = "not an enumeration" - sampled draws, and the permute(probabilities)
    // generator, which returns before the exact branch below.
    st_local("wdp_obsrow", "0")
    string scalar curfr
    curfr = st_framecurrent()
    st_framecurrent(fr)
    ons = st_data(., onsv)
    str = st_data(., strv)
    st_framecurrent(curfr)
    C = rows(ons)
    st_local("wdp_caut0", "")
    if (pname != "" & rows(st_matrix(pname)) > 0) {
        wooldid_permxp(ons, st_matrix(pname), st_matrix(lname), reps, assignname)
        return
    }
    // possible-assignment count: product over strata of multiset-permutation counts
    lnM = 0
    for (si = 1; si <= rows(uniqrows(str)); si++) {
        v = select(ons, str :== uniqrows(str)[si])
        lnM = lnM + lnfactorial(rows(v))
        levels = uniqrows(v)
        for (k = 1; k <= rows(levels); k++) {
            lnM = lnM - lnfactorial(sum(v :== levels[k] :| (v :>= . :& levels[k] :>= .)))
        }
    }
    M = (lnM < 700 ? round(exp(lnM)) : .)
    if (lnM < 34) nps = strofreal(round(exp(lnM)), "%20.0f")
    else nps = strofreal(exp(lnM - floor(lnM/ln(10))*ln(10)), "%9.4f") + "e+" + strofreal(floor(lnM/ln(10)), "%3.0f")
    st_local("wdp_npd", strofreal(M, "%20.0f"))
    st_local("wdp_nps", strtrim(nps))
    // TWO OF THESE FIVE NUMBERS ARE ARITHMETIC, THREE ARE POLICY. Arithmetic: 700 is the
    // double-overflow guard (ln of the largest double is about 709.78), which is what makes
    // the count deliberately MISSING rather than +inf; 34 is the exact-integer boundary
    // (exp(34) is about 5.8e14, under 2^53), so below it round() and the %20.0f print are
    // lossless and above it the code must switch to the mantissa/exponent string. Policy,
    // and tunable: capexact is the permexact ceiling - it bounds the C x M enumeration
    // matrix and the length of the rep loop, not just the user's patience - capwarn the
    // slow-run warning, and the 0.1*M below the sampling-overlap caution
    capexact = 10000
    capwarn = 1000
    if (exact) {
        if (M >= . | M > capexact) {
            st_local("wdp_err", "Error: permexact requires enumerating every assignment, but there are " + strtrim(nps) + " (ceiling: " + strofreal(capexact) + "). Use permreps().")
            return
        }
        if (M > capwarn) st_local("wdp_caut1", "permexact is enumerating " + strtrim(nps) + " assignments (over " + strofreal(capwarn) + "); this may be slow.")
        R = M
    }
    else {
        R = reps
        if (M < . & reps > 0.1*M) st_local("wdp_caut2", "permreps(" + strofreal(reps) + ") exceeds 10% of the " + strtrim(nps) + " possible assignments; consider permexact.")
    }
    st_local("wdp_R", strofreal(R))
    // NOTE the possible-assignment count and the enumeration are produced by two
    // INDEPENDENT algorithms: R comes from the log-factorial formula above, the columns of
    // A from the recursive enumerator below (this J() allocation is thrown away in the
    // exact branch). Only their AGREEMENT makes the rep loop correct, and nothing
    // reconciles them. The loop is driven by wdp_R and reads the matrix with el(), which
    // returns MISSING out of range rather than erroring (measured, Stata 19), and a missing
    // onset is wooldid's never-treated code - so R > cols(A) yields replications in which
    // every cluster is silently a control, and R < cols(A) yields a non-exhaustive
    // enumeration whose completeness flag still claims exact shares. The two routines agree
    // because both bucket every missing onset as ONE never-treated level (the
    // `:| (v :>= . :& levels[k] :>= .)' clause in the counter, the same test in
    // wooldid_permx_msp): exact for `.', wrong for extended missings (.a, .b), which
    // uniqrows() keeps as distinct levels
    A = J(C, R, .)
    if (exact) {
        // enumerate per-stratum multiset permutations, then take the cartesian product
        A = wooldid_permx_enum(ons, str)
        // BS-23 (author ruling 2026-09-07): locate the observed assignment among the
        // enumerated columns. It is always one of them (the enumeration is every distinct
        // multiset permutation), but NOT necessarily the first: each stratum is enumerated
        // from its SORTED onset vector, so the identity column sits wherever the observed
        // order falls. Missing onsets (never-treated) are one level here, exactly as in
        // wooldid_permx_msp's own multiset test.
        for (oc = 1; oc <= cols(A); oc++) {
            if (all((A[., oc] :== ons) :| (A[., oc] :>= . :& ons :>= .))) {
                st_local("wdp_obsrow", strofreal(oc))
                oc = cols(A)
            }
        }
    }
    else {
        for (r = 1; r <= R; r++) {
            v = ons
            for (si = 1; si <= rows(uniqrows(str)); si++) {
                idx = selectindex(str :== uniqrows(str)[si])
                v[idx] = jumble(ons[idx])
            }
            A[., r] = v
        }
    }
    st_matrix(assignname, A)
}

real matrix wooldid_permx_enum(real colvector ons, real colvector str)
{
    real matrix A, B, S1
    real colvector idx
    real scalar si, n1, n2, a, b
    // per-stratum enumerations, combined by cartesian product
    A = J(rows(ons), 0, .)
    for (si = 1; si <= rows(uniqrows(str)); si++) {
        idx = selectindex(str :== uniqrows(str)[si])
        S1 = wooldid_permx_msp(sort(ons[idx], 1))
        if (cols(A) == 0) {
            A = J(rows(ons), cols(S1), .)
            for (b = 1; b <= cols(S1); b++) A[idx, b] = S1[., b]
        }
        else {
            n1 = cols(A)
            n2 = cols(S1)
            B = J(rows(ons), n1*n2, .)
            for (a = 1; a <= n1; a++) {
                for (b = 1; b <= n2; b++) {
                    B[., (a-1)*n2 + b] = A[., a]
                    B[idx, (a-1)*n2 + b] = S1[., b]
                }
            }
            A = B
        }
    }
    return(A)
}

real matrix wooldid_permx_msp(real colvector v)
{
    real matrix R, S
    real colvector rest, lv
    real scalar k, j
    // all distinct permutations of multiset v (v sorted), columns = permutations
    if (rows(v) == 1) return(v)
    lv = uniqrows(v)
    R = J(rows(v), 0, .)
    for (k = 1; k <= rows(lv); k++) {
        j = selectindex(v :== lv[k] :| (v :>= . :& lv[k] :>= .))[1]
        if (j == rows(v)) rest = v[|1 \ j-1|]
        else if (j == 1) rest = v[|2 \ .|]
        else rest = v[|1 \ j-1|] \ v[|j+1 \ .|]
        S = wooldid_permx_msp(rest)
        R = R, (J(1, cols(S), lv[k]) \ S)
    }
    return(R)
}
end


*** F6: probabilities-method draws (design distribution): per-cluster categorical
*** draws from user-specified probabilities; count = product of support sizes
cap mata: mata drop wooldid_permxp()
mata:
void wooldid_permxp(real colvector ons, real matrix P, real rowvector lvl,
                    real scalar reps, string scalar assignname)
{
    real matrix A
    real colvector cum
    real scalar C, K, c, r, k, u, lnS, zp, obsk
    string scalar nps
    C = rows(ons)
    K = cols(lvl)
    lnS = 0
    zp = 0
    for (c = 1; c <= C; c++) {
        lnS = lnS + ln(sum(P[c, .] :> 0))
        obsk = 0
        for (k = 1; k <= K; k++) {
            if (ons[c] == lvl[k]) obsk = k
        }
        if (obsk == 0 | (obsk > 0 ? P[c, obsk] <= 0 : 0)) zp = 1
    }
    if (lnS < 34) nps = strofreal(round(exp(lnS)), "%20.0f")
    else nps = strofreal(exp(lnS - floor(lnS/ln(10))*ln(10)), "%9.4f") + "e+" + strofreal(floor(lnS/ln(10)), "%3.0f")
    st_local("wdp_npd", strofreal((lnS < 700 ? round(exp(lnS)) : .), "%20.0f"))
    st_local("wdp_nps", strtrim(nps))
    st_local("wdp_R", strofreal(reps))
    if (exp(lnS) < . & reps > 0.1*exp(lnS)) st_local("wdp_caut2", "permreps(" + strofreal(reps) + ") exceeds 10% of the " + strtrim(nps) + " possible assignments in the design's support.")
    if (zp) st_local("wdp_caut0", "the OBSERVED assignment has zero probability under the specified design, so the reported p-value is a tail probability under that design distribution, not a randomization test that includes the observed assignment.")
    A = J(C, reps, .)
    for (r = 1; r <= reps; r++) {
        for (c = 1; c <= C; c++) {
            // this loop can leave A[c, r] MISSING if u exceeds the last cumulative bin,
            // and a missing onset means never-treated downstream - so the cluster silently
            // becomes a control in that replication, rc 0, no message. It cannot happen
            // because wooldid_permp (NFPERM-06) rescales every row so its cumulative sum
            // ends at 1 to within about K*eps, which the 1e-12 cushion below covers; a P
            // matrix arriving by any other route must be normalized first. `k = K + 1' is
            // the deliberate loop break, written as a counter assignment
            cum = runningsum(P[c, .])'
            u = runiform(1, 1)
            for (k = 1; k <= K; k++) {
                if (u <= cum[k] + 1e-12) {
                    A[c, r] = lvl[k]
                    k = K + 1
                }
            }
        }
    }
    st_matrix(assignname, A)
}
end

*** F6: the shared permprobs() parser: validates the grammar against the data's onset
*** levels and cluster values; in engine mode (pname != "") also emits the probability
*** matrix with rows in ascending cluster order and the level vector in user order.
*** on any violation, sets the local permperr (caller raises r(198))
cap mata: mata drop wooldid_permp()
mata:
void wooldid_permp(string scalar spec, string scalar onsname, string scalar lvlstr,
                   string scalar clusstr, string scalar pname, string scalar lname)
{
    string vector grp, groups
    string scalar cur
    string rowvector tok, dlv, dcl
    real rowvector lvl, seen, dl, dc
    real matrix P
    real colvector cval, ord
    real scalar G, g, K, C, c, k, j, ssum, v, ing, hit
    st_local("permperr", "")
    grp = tokens(subinstr(subinstr(spec, "(", " __GRP__ ", .), ")", " __END__ ", .))
    // regroup
    groups = J(0, 1, "")
    cur = ""
    ing = 0
    for (g = 1; g <= cols(grp); g++) {
        if (grp[g] == "__GRP__") {
            if (ing) {
                st_local("permperr", "Error: permprobs() has a '(' inside a group; check the parentheses.")
                return
            }
            ing = 1
            cur = ""
        }
        else if (grp[g] == "__END__") {
            if (!ing) {
                st_local("permperr", "Error: permprobs() has an unmatched ')'.")
                return
            }
            ing = 0
            groups = groups \ cur
        }
        else if (ing) cur = cur + " " + grp[g]
        else {
            st_local("permperr", "Error: permprobs() has text outside parenthesized groups: " + grp[g])
            return
        }
    }
    if (ing) {
        st_local("permperr", "Error: permprobs() has an unmatched '('.")
        return
    }
    G = rows(groups)
    if (G < 2) {
        st_local("permperr", "Error: permprobs() needs the onset-level group plus one group per cluster.")
        return
    }
    // group 1: onset varname + level order
    tok = tokens(groups[1])
    if (cols(tok) < 3) {
        st_local("permperr", "Error: permprobs() first group must be (onsetvar level level ...).")
        return
    }
    if (tok[1] != onsname) {
        st_local("permperr", "Error: permprobs() first group names '" + tok[1] + "' but the timing variable is '" + onsname + "'.")
        return
    }
    K = cols(tok) - 1
    lvl = J(1, K, .)
    for (k = 1; k <= K; k++) {
        // NFPERM-05: strict token rule. only the literal "." names the never-treated
        // level; any other token must parse to a nonmissing number. previously a garbage
        // token strtoreal'd to missing, silently aliasing the never-treated level, and a
        // duplicate listed level passed the one-way membership check - either way a
        // realized onset level silently vanished from the randomization design (rc 0)
        if (tok[k+1] == ".") lvl[k] = .
        else {
            lvl[k] = strtoreal(tok[k+1])
            if (lvl[k] >= .) {
                st_local("permperr", "Error: permprobs() onset level token '" + tok[k+1] + "' is not a number or the never-treated marker '.'.")
                return
            }
        }
        // duplicate rejection (mata missing == missing is true, so this also catches a
        // second '.'); together with the count and exactly-once membership checks below
        // this makes listed-vs-realized levels a bijection
        for (j = 1; j < k; j++) {
            if (lvl[j] == lvl[k]) {
                st_local("permperr", "Error: permprobs() lists onset level '" + tok[k+1] + "' more than once.")
                return
            }
        }
    }
    // level set must equal the data's level set
    dlv = tokens(lvlstr)
    if (cols(dlv) != K) {
        st_local("permperr", "Error: permprobs() lists " + strofreal(K) + " onset levels but the data (within sample) has " + strofreal(cols(dlv)) + ".")
        return
    }
    dl = strtoreal(dlv)
    for (k = 1; k <= K; k++) {
        if (sum(dl :== lvl[k]) != 1) {
            st_local("permperr", "Error: permprobs() onset level '" + tok[k+1] + "' does not match the data's onset levels (" + lvlstr + ").")
            return
        }
    }
    // cluster groups
    C = G - 1
    dcl = tokens(clusstr)
    if (cols(dcl) != C) {
        st_local("permperr", "Error: permprobs() has " + strofreal(C) + " cluster group(s) but the data (within sample) has " + strofreal(cols(dcl)) + " cluster value(s).")
        return
    }
    cval = J(C, 1, .)
    P = J(C, K, .)
    for (g = 2; g <= G; g++) {
        tok = tokens(groups[g])
        if (cols(tok) != K + 1) {
            st_local("permperr", "Error: permprobs() cluster group (" + groups[g] + " ) must have the cluster value plus exactly " + strofreal(K) + " probabilities.")
            return
        }
        cval[g-1] = strtoreal(tok[1])
        ssum = 0
        for (k = 1; k <= K; k++) {
            v = strtoreal(tok[k+1])
            if (v >= . | v < 0) {
                st_local("permperr", "Error: permprobs() has an invalid probability in group (" + groups[g] + " ).")
                return
            }
            P[g-1, k] = v
            ssum = ssum + v
        }
        if (abs(ssum - 1) > 1e-6) {
            st_local("permperr", "Error: permprobs() probabilities for cluster value " + tok[1] + " sum to " + strofreal(ssum) + ", not 1.")
            return
        }
        // NFPERM-06: normalize the row so its cumulative sum ends exactly at 1 (within
        // ~K*eps, absorbed by the draw loop's +1e-12 cushion). rows accepted a hair
        // under 1 by the tolerance formerly left a deficit tail where a uniform draw
        // matched no probability column, leaving the assignment cell missing - which
        // downstream is wooldid's never-treated code (silent). division by a positive
        // ssum preserves zero/positive sign patterns, so the observed-support check and
        // support counts are unaffected; exact-sum rows are bitwise unchanged (ssum==1)
        if (ssum != 1) P[g-1, .] = P[g-1, .] :/ ssum
    }
    // cluster set must equal the data's, no duplicates
    dc = strtoreal(dcl)
    seen = J(1, C, 0)
    for (c = 1; c <= C; c++) {
        hit = 0
        for (k = 1; k <= C; k++) {
            if (dc[k] == cval[c]) {
                if (seen[k]) {
                    st_local("permperr", "Error: permprobs() lists cluster value " + strofreal(cval[c]) + " more than once.")
                    return
                }
                seen[k] = 1
                hit = 1
            }
        }
        if (!hit) {
            st_local("permperr", "Error: permprobs() cluster value " + strofreal(cval[c]) + " is not a value of the cluster variable (within sample).")
            return
        }
    }
    if (pname != "") {
        // rows in ascending cluster order = the engine's map-frame order
        // THAT EQUALITY IS NOT SELF-EVIDENT and is enforced 7,000 lines away:
        // wooldid_permrun sorts the map frame by (permstrata, clus), and wooldid_permctl
        // REFUSES permstrata() with permute(probabilities) - so under probabilities the
        // stratum key is constant and the frame order collapses to ascending cluster.
        // wooldid_permxp then pairs P's row c with ons's row c positionally, with no key.
        // If stratified probability designs are ever allowed (the caution text beside that
        // refusal actively recommends permstrata() for the shuffle case), every cluster
        // silently receives another cluster's probability row wherever cluster order
        // crosses stratum order - a wrong randomization design with rc 0. Key P to the
        // frame by cluster value at that point, the way NFPERM-02's __wdprow keys the
        // assignment matrix; that hazard was caught there by measurement (crossing-strata
        // exact designs reported p = 0 with the observed assignment absent from its own
        // enumeration) and the probability matrix has the same exposure
        ord = order(cval, 1)
        st_matrix(pname, P[ord, .])
        st_matrix(lname, lvl)
    }
}
end


*** NFADV-01/02: the adversarial p-value statistics, called by wooldid_advstat (kept
*** out of wooldid_worker to conserve compiled-size headroom; the cap mata drop makes
*** the file safely re-loadable within a session). Inputs are the F2 event-study
*** capture (b 1xk, V kxk) and the aligned relative-time column vector; planted
*** reference rows carry b = 0 with zero variance. Sets st_local advpstat = "ok" and
*** the scalar `pn' on success, advpstat = "degen" otherwise (never errors).
***
*** alt == 0, primary (NFADV-01), REWRITTEN under BS-12 (author ruling 2026-09-08; research
*** record review_process/research_2026-09-08/adversarial_null.md). Minimum-distance test of
*** the ONE-parameter null that the horizon effects lie on a line THROUGH THE REFERENCE
*** PERIOD. Under esfixedbaseperiod b_h = theta_h - theta_r, so a linear underlying path
*** implies b_h = c*(h - esrel): the design is X = relyr :- esrel, one free parameter, and
*** NOT the free-intercept pair X = [1, relyr] this test used to fit. c is estimated by GLS
*** under a generalized inverse of V; the minimized distance W = min_c (b - Xc)' V^+ (b - Xc)
*** is computed in the whitened eigenspace of V (z = L^-1/2 U' b, M = L^-1/2 U' X,
*** W = ||z - M chat||^2), chi-squared with q = rank(V) - 1 under the null.
*** Under timetrends TWO periods are normalised (esrel and esrel - 1, both planted at 0), so
*** the only line through both is the zero line: the null collapses to b = 0, W = z'z and
*** q = k1. That branch is written out rather than handed a 0-column design because svdsv()
*** on a k1 x 0 matrix is not a shape to put inside a replication loop.
*** WHY THE FREE INTERCEPT WENT, since the old form is the one a reader will expect: the
*** constant vector lies in col([1, relyr]), so the old statistic was exactly invariant to
*** b :+ c and its power against a reference-period level shift - the signature of an
*** anomalous or anticipatory reference year - was exactly alpha at any magnitude.
*** col(relyr - esrel) is nested in col([1, relyr]), so the new statistic's noncentrality is
*** at least as large against EVERY alternative; the price is one degree of freedom, q rising
*** from rank(V) - 2 to rank(V) - 1.
*** NOTE ON THE REFERENCE PERIOD (comment corrected 2026-09-08; the previous text claimed
*** the planted zero-variance reference rows receive zero weight under V^+ and are
*** mathematically inert, so that omitting them is identical to including them - both halves
*** of that were wrong). The reference row is not in b/V at all: relyrcolumn is built only
*** from relative periods carrying an estimable coefficient, so under esfixedbaseperiod the
*** omitted base level never enters. wooldid_advstat plants the reference rows into the
*** DISPLAY frame only, and passes the unplanted b/V here.
***
*** alt == 1, altmethod (NFADV-02): the forecast test. P = the altmethod fit set
*** (real relyr <= -1 plus the planted reference row(s)), Q = the real post horizons.
*** A = Z_Q (Z_P' Z_P)^-1 Z_P' is the unweighted-OLS forecast operator (the same fit
*** as the display line's pre-period cut), d = b_Q - A b_P, and the test covariance
*** propagates the fitted line's uncertainty and the pre/post covariance:
*** S = V_QQ + A V_PP A' - A V_PQ - V_QP A'; W = d' S^+ d, q = rank(S). The S^+,
*** rank, and the residual-in-column-space check share one eigensystem and one
*** scale-aware tolerance (a residual component outside col(S) is a degeneracy, not
*** a testable direction).
***
*** Both statistics use Ftail(q, dfr, W/q) when dfr is a valid residual df and
*** chi2tail(q, W) otherwise (the wooldid_ptp convention).
cap mata: mata drop wooldid_advx()
cap mata: mata drop wooldid_advW()
cap mata: mata drop wooldid_bootadvr()
cap mata: mata drop wooldid_bootadvobs()
cap mata: mata drop wooldid_bootadvrel()
cap mata: mata drop wooldid_bootcellkey()
cap mata: mata drop wooldid_bootcellkeys()
cap mata: mata drop wooldid_bootcellkeym()
cap mata: mata drop wooldid_bootcells()
cap mata: mata drop wooldid_bootsavemk()
cap mata: mata drop wooldid_bootwnames()
cap mata: mata drop wooldid_bootwrow()
cap mata: mata drop wooldid_bootwrowx()
cap mata: mata drop wooldid_bootcellwire()
mata:
// the adversarial statistic W (minimum-distance line test, alt==0; forecast test, alt==1)
// as a PURE function of (b, V, rel) - factored out of wooldid_advx so the bootstrap can
// evaluate the SAME construction on a replication's recentred coefficients and its own
// covariance. Returns W, or missing when the covariance is degenerate / untestable; sets
// the second output to the residual degrees of freedom q
// "PURE" ABOVE IS NOT LITERAL: Mata passes by ADDRESS, so the `V = (V + V')/2'
// symmetrization below MUTATES THE CALLER'S MATRIX (verified, Stata 19 MP: a caller's
// X = (1,2\3,4) comes back as (1,2\2.5,4)). That is safe only because (V+V')/2 is bitwise
// idempotent on an already-symmetric V - transpose is exact, x+x is exact, /2 is exact - and
// wooldid_permadv leans on it: it builds ONE covariance, computes W_obs before the first
// call and every W_r after, and its membership shortcut `memb = (Wr[orow] < .)' is sound
// only because W_r at the observed row equals W_obs bitwise, which needs one identical
// metric across all those calls. Never add a non-idempotent in-place touch of V here; copy
// first. The same by-address rule is what makes qout an output argument - callers must pass
// an initialised variable, not an expression.
real scalar wooldid_advW(real colvector b, real matrix V, real colvector rel,
                         real scalar esrel, real scalar tt, real scalar alt, real scalar qout)
{
    real colvector z, g, r, sel, idxP, idxQ, bP, relP, d, c, resid
    real matrix X, U, T, M, ZP, ZQ, ZZi, Aop, VPP, VPQ, VQQ, S, US
    real rowvector lam, sv, lamS
    real scalar k1, np, npr, nq, scale, tol, rankM, q, W
    qout = .
    if (hasmissing(V) | hasmissing(b) | hasmissing(rel)) return(.)
    if (alt == 0) {
        // BS-12 (author ruling 2026-09-08; research record
        // review_process/research_2026-09-08/adversarial_null.md): the ONE-parameter null.
        // Under esfixedbaseperiod b_h = theta_h - theta_r, so a linear underlying path
        // implies b_h = c*(h - esrel): a line THROUGH the reference period, one free
        // parameter, not the free-intercept pair X = [1, relyr] this used to fit. The
        // constant vector lies in col([1, relyr]), so the old statistic was exactly
        // invariant to b :+ c and its power against a reference-period level shift - the
        // signature of an anomalous or anticipatory reference year - was exactly alpha at
        // any magnitude. col(relyr - esrel) is nested in col([1, relyr]), so the new
        // statistic's noncentrality is at least as large against EVERY alternative; the
        // price is one degree of freedom, q rising from rank(V) - 2 to rank(V) - 1.
        // timetrends normalises TWO periods (esrel and esrel - 1, both planted at 0), so
        // the only line through both is the zero line and the null collapses to b = 0:
        // no free parameter at all, W = z'z with q = k1. Branch on tt rather than build a
        // 0-column design - svdsv() on a k1 x 0 matrix is not a shape to hand a
        // replication loop.
        V = (V + V')/2
        U = .
        lam = .
        symeigensystem(V, U, lam)
        scale = max(abs(lam))
        if (scale <= 0 | scale >= .) return(.)
        tol = scale*1e-12
        sel = selectindex(lam' :> tol)
        k1 = rows(sel)
        if (k1 == 0) return(.)
        T = diag(1:/sqrt(lam[sel]')) * U[., sel]'
        z = T*b
        if (tt) {
            q = k1
            W = z'z
        }
        else {
            X = rel :- esrel
            M = T*X
            sv = svdsv(M)
            rankM = sum(sv :> max(sv)*1e-12)
            q = k1 - rankM
            if (q <= 0) return(.)
            g = svsolve(M, z)
            r = z - M*g
            W = r'r
        }
    }
    else {
        idxP = selectindex(rel :<= -1)
        idxQ = selectindex(rel :>= 0)
        npr = rows(idxP)
        nq = rows(idxQ)
        if (npr == 0 | nq == 0) return(.)
        bP = b[idxP] \ 0
        relP = rel[idxP] \ esrel
        if (tt) {
            bP = bP \ 0
            relP = relP \ (esrel - 1)
        }
        np = rows(relP)
        ZP = J(np, 1, 1), relP
        ZQ = J(nq, 1, 1), rel[idxQ]
        ZZi = invsym(ZP'ZP)
        if (ZZi[1,1] == 0 | ZZi[2,2] == 0) return(.)
        Aop = ZQ*ZZi*ZP'
        VPP = J(np, np, 0)
        VPP[|1,1 \ npr,npr|] = V[idxP, idxP]
        VPQ = J(np, nq, 0)
        VPQ[|1,1 \ npr,nq|] = V[idxP, idxQ]
        VQQ = V[idxQ, idxQ]
        d = b[idxQ] - Aop*bP
        S = VQQ + Aop*VPP*Aop' - Aop*VPQ - VPQ'*Aop'
        S = (S + S')/2
        US = .
        lamS = .
        symeigensystem(S, US, lamS)
        scale = max(abs(lamS))
        if (scale <= 0 | scale >= .) return(.)
        tol = scale*1e-12
        sel = selectindex(lamS' :> tol)
        q = rows(sel)
        if (q == 0) return(.)
        c = US[., sel]'d
        resid = d - US[., sel]*c
        if (norm(resid) > 1e-8*max((norm(d), sqrt(scale)))) return(.)
        W = sum(c:^2 :/ lamS[sel]')
    }
    if (W >= .) return(.)
    qout = q
    return(W)
}

void wooldid_advx(string scalar bn, string scalar vn, string scalar rn,
                  real scalar esrel, real scalar tt, real scalar alt,
                  real scalar dfr, string scalar pn)
{
    real colvector b, rel, z, g, r, sel, idxP, idxQ, bP, relP, d, c, resid
    real matrix V, X, U, T, M, ZP, ZQ, ZZi, Aop, VPP, VPQ, VQQ, S, US
    real rowvector lam, sv, lamS
    real scalar k1, np, npr, nq, scale, tol, rankM, q, W, p
    st_local("advpstat", "degen")
    b = st_matrix(bn)'
    V = st_matrix(vn)
    rel = st_matrix(rn)
    if (hasmissing(V) | hasmissing(b) | hasmissing(rel)) return
    q = .
    W = wooldid_advW(b, V, rel, esrel, tt, alt, q)
    if (W >= .) return
    if (dfr < .) p = Ftail(q, dfr, W/q)
    else p = chi2tail(q, W)
    if (p >= .) return
    st_numscalar(pn, p)
    st_local("advpstat", "ok")
}
end


* F-boot mata kernels. drop guards: the file must stay session-reloadable
cap mata: mata drop wooldid_bootdw()
cap mata: mata drop wooldid_boottp()
cap mata: mata drop wooldid_bootgrab()
cap mata: mata drop wooldid_bootfin()
cap mata: mata drop wooldid_bootves()
cap mata: mata drop wooldid_bootdist()
cap mata: mata drop wooldid_bootqf()
cap mata: mata drop wooldid_bootpw()
cap mata: mata drop wooldid_bootpp()
cap mata: mata drop wooldid_bootpwt()
cap mata: mata drop wooldid_bootjwr()
cap mata: mata drop wooldid_bootsgw()
cap mata: mata drop wooldid_bootjwobs()
cap mata: mata drop wooldid_bootjobs()
cap mata: mata drop wooldid_bootpairidx()
cap mata: mata drop wooldid_bootci()
cap mata: mata drop wooldid_bootbca()
cap mata: mata drop wooldid_bootppa()
cap mata: mata drop wooldid_bootqs()
cap mata: mata drop wooldid_bootveskey()
cap mata: mata drop wooldid_permq()
cap mata: mata drop wooldid_permmemb()
cap mata: mata drop wooldid_permtail()
cap mata: mata drop wooldid_bootptails()
cap mata: mata drop wooldid_permstat()
cap mata: mata drop wooldid_permpw()
cap mata: mata drop wooldid_permqf()
cap mata: mata drop wooldid_permstep()
cap mata: mata drop wooldid_permadv()
cap mata: mata drop wooldid_permcellwire()
cap mata: mata drop wooldid_permsaveasg()
cap mata: mata drop wooldid_permcse()
cap mata: mata drop wooldid_permadvline()
cap mata: mata drop wooldid_vspl()
cap mata: mata drop wooldid_esrelx()
cap mata: mata drop wooldid_esbanda()
cap mata: mata drop wooldid_bootband()
mata:
// draw matrix: per-rep (outer) per-cluster-ascending (inner) Gamma(alpha,1) draws -
// the exact rwgen bayes stream-consumption order; alpha = n_c (size) or 1 (equal);
// row multiplier = (g/colsum)/n_c. writes a C x R matrix (transposed for el() reads)
void wooldid_bootdw(string scalar fr, string scalar amode, real scalar R, string scalar outname)
{
  string scalar cur
  real matrix N, E, P
  real scalar C, r, k, a
  cur = st_framecurrent()
  st_framecurrent(fr)
  N = st_data(., "__wdbn")
  st_framecurrent(cur)
  C = rows(N)
  E = J(C, R, .)
  for (r = 1; r <= R; r++) {
    for (k = 1; k <= C; k++) {
      a = (amode == "size" ? N[k] : 1)
      E[k, r] = rgamma(1, 1, a, 1)
    }
  }
  P = E :/ colsum(E)
  // bootalpha(size) only: spread the cluster's Dirichlet mass EQUALLY over its rows, so
  // that every row carries p_c / n_c. The divisor has to be proportional to the Gamma
  // SHAPE, which under size is the cluster's row count n_c: E[p_c] = n_c/N, so a row of
  // cluster c carries expected multiplier 1/N whatever its cluster's size or weight, and
  // the replications centre on the same weighted estimand the point estimate targets.
  // Dividing by the base-weight TOTAL W_c instead (the form shipped until 2026-09-13) put
  // expected multiplier n_c/(N*W_c) on each row, so a cluster's expected replicate mass
  // was its ROW-COUNT share n_c/N rather than its weight share W_c/W - the replications
  // bootstrapped the row-count-weighted estimand while the point estimate was the user's.
  // Unweighted it is invisible because W_c == n_c identically (a total() of ones), and
  // weighted it displaced the replication distribution from the estimate by a measured
  // 3.5 to 9.4 bootstrap deviations with every draw on one side, which also killed the
  // whole bias-corrected family (measured 2026-09-13 on dta_s_stateyear, cluster(st),
  // [aw=pop]: max |replicate mean - estimate| / bootstrap sd 9.36 -> 0.10, and
  // e(boot_bc_absent) six columns -> empty). The comment that stood here cited "p_c / W_c -
  // the bayesboot / rwgen priorpowers convention": that convention was MEASURED and is
  // p_c / n_c - rwgen bayes puts p_c/n_c on each row with the weight total appearing
  // nowhere, the recovered Dirichlet parameters match the row counts, and bayesboot
  // refuses user weights outright. See test_bootstrap_bayes.do O2/O2d and
  // review_process/quality_review_2026-09-05.md F001.
  // Under bootalpha(equal) (the default) every row carries p_c itself, so the cluster's
  // TOTAL replicate mass is W_c * p_c: Rubin's Bayesian bootstrap with clusters as the
  // exchangeable units, whose replications stay centred on the observation-weighted
  // estimand. Until 2026-09-06 the split applied under equal too, which centred the
  // replications on the equal-cluster-weighted estimand (rep mean .80 vs att 1.69 on a
  // 12-state unbalanced fixture)
  if (amode == "size") P = P :/ N
  st_matrix(outname, P)
}

// name-aligned grab of the current rep e(b) into row r of the accumulation matrix;
// ame-family structural zeros (rep value exactly 0 where the observed b is nonzero)
// are blanked and counted (A12b)
// boottstat: the studentized (percentile-t) tail probability and critical value, per
// column. t_r = (b_r - b_obs)/se_r studentizes each replication by ITS OWN standard
// error, and the observed statistic t_obs = b_obs/se_obs uses the observed run's
// same-formula (analytic) standard error - the two have to be the same functional or the
// statistic is not asymptotically pivotal and the refinement is lost. The tail count is
// the rank of t_obs among the t_r WITH the add-one Monte-Carlo-null adjustment (F033,
// author ruling 2026-09-07); the critical value stays a plain order statistic, so the p
// and the interval are no longer exact duals (see below).
// Columns with fewer than 2 usable replications, or no usable observed standard error,
// return missing and the caller leaves the Wald numbers in place.
void wooldid_boottp(string scalar accname, string scalar seaccname, string scalar obsbname,
                    string scalar aseobsname, string scalar exclname, real scalar lvl,
                    real scalar ties, string scalar outname)
{
  // Studentized (percentile-t) inference on the SIGNED t-draws t_r = (b_r - b_obs)/se_r,
  // referred against t_obs = b_obs/se_obs (both analytic standard errors, so the two are
  // the same functional - the property the replication-clustering fix restores).
  // Author rulings 2026-09-03/04: the interval is EQUAL-TAILED from the signed quantiles,
  // [b_obs - q_hi*se_obs, b_obs - q_lo*se_obs], NOT the earlier symmetric b +/- q(|t|)*se.
  // The p is the two-sided rank of t_obs among the t_r under the ADD-ONE Monte-Carlo-null
  // convention, 2*(1 + min(m0, m1))/(n + 1) (F033, author ruling 2026-09-07). Unlike the
  // percentile family in wooldid_bootci, these draws ARE a null reference set: t_r is
  // built from (b_r - b_obs), so it is recentred at the null and the observed statistic
  // legitimately joins it - the same argument that already governs every bootstrapped
  // joint statistic in the package. That makes one convention for every recentred
  // rank-based p here, and removes the measured six-percent size at 99 replications.
  // The price is that the p and the interval are no longer exact duals: the critical
  // value stays a plain order statistic (wooldid_bootci's rule), so at the boundary an
  // interval can exclude the null while the p sits just above the level. The p now has a
  // floor of 2/(n + 1) and can no longer print as exactly 0. The former
  // BC-t twin (rows 6-9: z0 from the rank of 0 among the centred t-draws, then Efron's
  // shift on the reflected quantiles) is WITHDRAWN (quality review F003, author ruling
  // 2026-09-06): the shift moved the interval against the median bias, and percentile-t
  // is already second-order accurate; bootdisplay(bc) is refused with boottstat.
  // Output: a 6 x K matrix, rows = n, t_obs, p_pct, qlo_pct, qhi_pct, live. Columns
  // with fewer than 2 usable replications, no usable observed standard error, or zero
  // dispersion post missing with live = 0.
  real matrix A, S, ob, ase, ex, out
  real colvector tr, v
  real scalar K, R, k, r, n, so, to, al, x, m0, m1, sd
  A = st_matrix(accname)
  S = st_matrix(seaccname)
  ob = st_matrix(obsbname)
  ase = st_matrix(aseobsname)
  ex = st_matrix(exclname)
  K = cols(A)
  R = rows(A)
  out = J(6, K, .)
  out[6, .] = J(1, K, 0)
  al = (100 - lvl)/200
  for (k = 1; k <= K; k++) {
    if (ex[1, k] != 0) continue
    so = ase[1, k]
    if (so >= . | so <= 0 | ob[1, k] >= .) continue
    to = ob[1, k] / so
    tr = J(0, 1, .)
    for (r = 1; r <= R; r++) {
      if (A[r, k] >= . | S[r, k] >= . | S[r, k] <= 0) continue
      tr = tr \ ((A[r, k] - ob[1, k]) / S[r, k])
    }
    n = rows(tr)
    out[1, k] = n
    if (n < 2) continue
    v = sort(tr, 1)
    sd = sqrt(variance(v))
    if (sd <= 0 | sd >= .) continue
    out[2, k] = to
    // percentile-t: equal-tailed signed quantiles, wooldid_bootci's rule
    x = n*al
    if (x == floor(x) & x >= 1) out[4, k] = (v[x] + v[min((x+1, n))])/2
    else out[4, k] = v[max((min((ceil(x), n)), 1))]
    x = n*(1-al)
    if (x == floor(x) & x >= 1) out[5, k] = (v[x] + v[min((x+1, n))])/2
    else out[5, k] = v[max((min((ceil(x), n)), 1))]
    m0 = sum(v :<= to)
    m1 = sum(v :>= to)
    out[3, k] = min((1, 2*(min((m0, m1)) + 1)/(n + 1)))
    out[6, k] = 1
  }
  st_matrix(outname, out)
}

void wooldid_bootgrab(string scalar accname, real scalar r, string scalar obsbname,
                      | string scalar seaccname, string scalar pwseaccname)
{
  real matrix A, rb, ob, S, rV, P, rP
  string matrix an, rn, pn, rpn
  real scalar k, j, nb, dose, v, dopw
  transmorphic t, tp
  A = st_matrix(accname)
  an = st_matrixcolstripe(accname)
  rb = st_matrix("e(b)")
  rn = st_matrixcolstripe("e(b)")
  ob = st_matrix(obsbname)
  // boottstat: this rep's own standard errors, off the diagonal of its e(V), aligned by
  // the same colstripe lookup. absent (fastnose reps) or non-positive entries stay
  // missing, and the studentization drops those replications for that column
  dose = (args() >= 4 && seaccname != "")
  if (dose) {
    // ORDERING CONSTRAINT: this read must stay AHEAD of the conformability test below that
    // clears `dose'. The writeback near the end of this function is guarded on
    // `seaccname != ""', NOT on `dose', so it pushes S back unconditionally - and it is
    // only safe because S was already populated by the time dose flipped. Restructure this
    // the natural way (test conformability first, skip the read when the rep's e(V) does
    // not conform) and S stays Mata's default 0 x 0 matrix, which then overwrites the
    // accumulator and wipes every replication's standard errors gathered so far. The
    // failure is silent in the way that matters: later replications refill into whatever
    // the accumulator became and the studentized columns merely report fewer usable draws.
    // The pairwise path below is structured the safe way (its read and its writeback sit
    // inside the same test); this one is not
    S = st_matrix(seaccname)
    rV = st_matrix("e(V)")
    if (rows(rV) != rows(rn)) dose = 0
  }
  t = asarray_create()
  for (j = 1; j <= rows(rn); j++) asarray(t, rn[j,1] + ":" + rn[j,2], j)
  nb = 0
  for (k = 1; k <= rows(an); k++) {
    if (asarray_contains(t, an[k,1] + ":" + an[k,2])) {
      j = asarray(t, an[k,1] + ":" + an[k,2])
      A[r, k] = rb[1, j]
      // A12b, mechanism (review_process/bootstrap_plan.md): in ame-family columns a
      // replication can come back at exact 0.0 not because the effect is zero but because
      // the resample lost dose variation and fastnose skips the amescreen that would have
      // caught it - MEASURED. So the rep is treated as UNAVAILABLE for that column rather
      // than as a draw at zero. Two things to know: the family test is an unanchored
      // SUBSTRING search on the coefficient name, so any future coefficient whose name
      // merely contains "ame" inherits the rule; and a genuinely-zero ame draw is
      // discarded with the structural ones, so the retained distribution is conditional on
      // being away from zero - narrowing that column's interval and shrinking its
      // percentile p. wdb_gblank counts the blanks but cannot separate the two
      if (strpos(an[k,2], "ame") & rb[1, j] == 0 & ob[1, k] != 0 & ob[1, k] < .) {
        A[r, k] = .
        nb++
      }
      if (dose & A[r, k] < .) {
        v = rV[j, j]
        if (v > 0 & v < .) S[r, k] = sqrt(v)
      }
    }
  }
  st_matrix(accname, A)
  st_matrixcolstripe(accname, an)
  if (args() >= 4 && seaccname != "") {
    st_matrix(seaccname, S)
    st_matrixcolstripe(seaccname, an)
  }
  // boottstat: this rep's own pairwise-contrast standard errors, matched to the OBSERVED
  // contrast rows by full row name (eq:JvsK) - never by position, since a resample can
  // lose a subgroup level and shorten the table. Unmatched or non-positive entries stay
  // missing and the studentization drops those replications for that contrast
  dopw = (args() >= 5 && pwseaccname != "")
  if (dopw) {
    rP = st_matrix("e(pwsgs)")
    if (rows(rP) > 0) {
      P = st_matrix(pwseaccname)
      pn = st_matrixcolstripe(pwseaccname)
      rpn = st_matrixrowstripe("e(pwsgs)")
      tp = asarray_create()
      for (j = 1; j <= rows(rpn); j++) asarray(tp, rpn[j,1] + ":" + rpn[j,2], j)
      for (k = 1; k <= rows(pn); k++) {
        if (asarray_contains(tp, pn[k,1] + ":" + pn[k,2])) {
          j = asarray(tp, pn[k,1] + ":" + pn[k,2])
          v = rP[j, 2]
          if (v > 0 & v < .) P[r, k] = v
        }
      }
      st_matrix(pwseaccname, P)
      st_matrixcolstripe(pwseaccname, pn)
    }
  }
  st_local("wdb_gblank", strofreal(nb))
}

// finalize: per-column mean/se over available reps (R-1, deviations from the rep
// mean); LISTWISE complete-case covariance over non-excluded columns (A11);
// sentinel -1e100 diagonals + zero off-diagonals for excluded/<2-available columns;
// low-availability disclosure via locals
// THE TWO MOMENT PASSES RUN ON DIFFERENT ROW SETS, and both are shipped to the user under
// names that imply they agree. se[k] is computed over rows where column k alone is
// available; V[k,k] over rows complete across every non-excluded column (that count is
// reported as e(boot_N_V)). Whenever any column loses a replication - a failed fit, or an
// A12b ame blank from wooldid_bootgrab - the two disagree, so sqrt(diag(e(V_boot))) is not
// the printed se column and a hand-run test or lincom does not exactly reproduce the
// table's Wald numbers. A11 ruled this a documented consequence rather than a defect.
void wooldid_bootfin(string scalar accname, string scalar exclname, string scalar sename, string scalar meanname, string scalar vname, string scalar cntname)
{
  real matrix A, V, CC, dev
  real rowvector excl, se, mn, cnt
  real colvector use
  real colvector ok
  string matrix an
  real scalar K, k, j, ndone, nlw
  string scalar lowav, exn
  A = st_matrix(accname)
  an = st_matrixcolstripe(accname)
  excl = st_matrix(exclname)
  K = cols(A)
  se = J(1, K, .)
  mn = J(1, K, .)
  cnt = J(1, K, 0)
  // completed reps = rows with at least one nonmissing cell
  ok = rowsum(A :< .) :> 0
  ndone = sum(ok)
  lowav = ""
  exn = ""
  for (k = 1; k <= K; k++) {
    use = selectindex(A[., k] :< .)
    cnt[k] = length(use)
    if (excl[k] == 0 & cnt[k] >= 2) {
      mn[k] = mean(A[use, k])
      se[k] = sqrt(variance(A[use, k]))
      if (cnt[k] < 0.5*ndone) lowav = lowav + " " + an[k,2] + "(" + strofreal(cnt[k]) + ")"
    }
    if (excl[k] == 1) exn = exn + " " + an[k,2]
  }
  // listwise complete-case rows over non-excluded columns
  V = J(K, K, 0)
  nlw = 0
  real rowvector keep
  keep = selectindex(excl :== 0)
  if (length(keep) > 0) {
    CC = A[., keep']
    ok = rowmissing(CC) :== 0
    if (sum(ok) >= 2) {
      CC = select(CC, ok)
      nlw = rows(CC)
      dev = CC :- mean(CC)
      V[keep, keep] = (dev' * dev) / (nlw - 1)
    }
  }
  for (k = 1; k <= K; k++) {
    // -1e100 IS THE PACKAGE'S "this variance does not exist; display missing" SENTINEL -
    // the same one the partial-identification rescue writes - and not a number to be read.
    // It survives the trip into e(V) because `ereturn repost V=' does not re-run the
    // positive-definiteness screen that `ereturn post' applies (measured; see the notes at
    // the repost site and in wooldid_vspl). So the matrix leaving this function is
    // deliberately NON-PSD, and a user's hand test or lincom touching a sentinel column
    // builds a quadratic form on a -1e100 variance with no error raised. NOTE the branch
    // also fires on V[k,k] == 0, i.e. when the LISTWISE pass found no variance for column
    // k - the per-column se[k] above may still be finite and displayed, so a column can
    // show a standard error in the table and carry the sentinel in e(V)
    if (excl[k] == 1 | cnt[k] < 2 | V[k,k] == 0) {
      V[k, .] = J(1, K, 0)
      V[., k] = J(K, 1, 0)
      V[k, k] = -1e100
    }
  }
  st_matrix(sename, se)
  st_matrix(meanname, mn)
  st_matrix(vname, V)
  st_matrix(cntname, cnt)
  st_local("wdb_lowav", strtrim(lowav))
  st_local("wdb_exclnames", strtrim(exn))
  st_local("wdb_nlw", strofreal(nlw))
}

// V_es subblock: listwise covariance of the acc columns matching the analytic
// V_es matrix's row names (relyr-anchored); succeeds only if every name maps and
// >= 2 listwise rows exist
void wooldid_bootves(string scalar accname, string scalar exclname, string scalar fk, string scalar oldname, string scalar newname)
{
  real matrix A, CC, dev, V
  real rowvector idx, excl
  real colvector ok
  string matrix an, vn
  string scalar key, nm
  real scalar m, k, j, found
  transmorphic t
  A = st_matrix(accname)
  an = st_matrixcolstripe(accname)
  excl = st_matrix(exclname)
  vn = st_matrixrowstripe(oldname)
  m = rows(vn)
  t = asarray_create()
  for (j = 1; j <= rows(an); j++) asarray(t, an[j,2], j)
  idx = J(1, m, .)
  // V_es row names are pre#/post# labels; the e(b) name is ES_<fk-ish>_<label>.
  // map by suffix match on the estimand family embedded in fk
  for (k = 1; k <= m; k++) {
    nm = vn[k,2]
    // candidate e(b) names for this family
    key = wooldid_bootveskey(fk, nm)
    if (asarray_contains(t, key)) idx[k] = asarray(t, key)
  }
  if (missing(idx)) {
    st_local("wdb_vesok", "0")
    return
  }
  for (k = 1; k <= m; k++) {
    if (excl[idx[k]] == 1) {
      st_local("wdb_vesok", "0")
      return
    }
  }
  CC = A[., idx']
  ok = rowmissing(CC) :== 0
  if (sum(ok) < 2) {
    st_local("wdb_vesok", "0")
    return
  }
  CC = select(CC, ok)
  dev = CC :- mean(CC)
  V = (dev' * dev) / (rows(CC) - 1)
  st_matrix(newname, V)
  st_matrixrowstripe(newname, vn)
  st_matrixcolstripe(newname, st_matrixcolstripe(oldname))
  st_local("wdb_vesok", "1")
}

// e(b) column name for V_es family fk + relyr label nm:
// fk in {att, att_it, ame, ame_it} [+ _s<sg>]: ES_<base>_<label> or
// ES_<base>_SG<sg>_<label> with base = att/att_it/ame/ame_it
string scalar wooldid_bootveskey(string scalar fk, string scalar nm)
{
  if (regexm(fk, "^(att_it|att|ame_it|ame)_s([0-9]+)$")) {
    return("ES_" + regexs(1) + "_SG" + regexs(2) + "_" + nm)
  }
  return("ES_" + fk + "_" + nm)
}

// dense bootdist (completed rows only)
void wooldid_bootdist(string scalar accname, string scalar outname, | string scalar seaccname, string scalar seoutname)
{
  // rows are kept iff the COEFFICIENT row has any non-missing entry, and the same mask is
  // applied to the standard-error accumulator when one is given - so e(bootdist) and
  // e(bootdist_se) stay row-aligned and a user pairing them row by row studentizes the
  // same replication (previously each was compacted on its own mask, so a replication
  // with coefficients but no usable VCE was kept in one and dropped from the other)
  real matrix A, S
  real colvector ok
  A = st_matrix(accname)
  ok = rowsum(A :< .) :> 0
  st_matrix(outname, select(A, ok))
  if (args() >= 4 && seaccname != "") {
    S = st_matrix(seaccname)
    st_matrix(seoutname, select(S, ok))
  }
}


// observed pre-horizon coefficients for every event-study family key that exists, stacked
// into one row vector with colstripe (fk, prename). Built while the OBSERVED e() is live
// (wooldid_bootrun, before any replication), and carried into the rep loop by name
void wooldid_bootjobs(string scalar keys, string scalar outname, | real scalar allh)
{
  real scalar useall
  useall = (args() >= 3 ? allh : 0)
  string rowvector ks
  string matrix cs, out_cs
  real rowvector b, out
  real scalar i, j
  ks = tokens(keys)
  out = J(1, 0, .)
  out_cs = J(0, 2, "")
  for (i = 1; i <= cols(ks); i++) {
    b = st_matrix("e(b_es_" + ks[i] + ")")
    if (cols(b) == 0) continue
    cs = st_matrixcolstripe("e(b_es_" + ks[i] + ")")
    for (j = 1; j <= rows(cs); j++) {
      if (!useall & substr(cs[j, 2], 1, 3) != "pre") continue
      out = out, b[1, j]
      out_cs = out_cs \ (ks[i], cs[j, 2])
    }
  }
  if (cols(out) == 0) {
    st_local("wdb_jobsn", "0")
    return
  }
  st_matrix(outname, out)
  st_matrixcolstripe(outname, out_cs)
  st_local("wdb_jobsn", strofreal(cols(out)))
}

// m x m matrix of e(pwsgs) row indices for the pair (a,b) of the ordered subgroup list,
// for one (eq, item): rows are named "JvsK" (with an optional bn/b marker on either level),
// J = the higher level minus K = the lower, in the eq given. 0 on the diagonal; missing
// where a pair has no row (the caller then declines the polarisation)
void wooldid_bootpairidx(string scalar eqwant, string scalar sglist, string scalar outname)
{
  string matrix rs
  string rowvector sg
  real matrix M
  real scalar m, a, b, j, J, K
  string scalar nm
  rs = st_matrixrowstripe("e(pwsgs)")
  sg = tokens(sglist)
  m = cols(sg)
  M = J(m, m, .)
  for (a = 1; a <= m; a++) M[a, a] = 0
  for (j = 1; j <= rows(rs); j++) {
    if (rs[j, 1] != eqwant) continue
    nm = rs[j, 2]
    if (!regexm(nm, "^([0-9]+)(bn?)?vs([0-9]+)(bn?)?([.].*)?$")) continue
    J = strtoreal(regexs(1))
    K = strtoreal(regexs(3))
    for (a = 1; a <= m; a++) {
      for (b = 1; b <= m; b++) {
        if (strtoreal(sg[a]) == J & strtoreal(sg[b]) == K) {
          M[a, b] = j
          M[b, a] = j
        }
      }
    }
  }
  st_matrix(outname, M)
}

// ---- bootstrapped CELL surfaces: e(gteffects) and e(cfx*) (task 8, ruling 2026-09-04) ----
// key for a gteffects row: "i|onset|relyr" from its first, second and fourth columns; key
// for a cfx row: its row name (the lattice level). An R x ncells accumulator carries the
// cell estimate (or, in the *se twin, the cell's analytic se) per replication, matched by
// key so a replication that lacks a cell leaves it missing rather than shifting the rest
string colvector wooldid_bootcellkeym(real matrix M)
{
  string colvector k
  real scalar i
  k = J(rows(M), 1, "")
  for (i = 1; i <= rows(M); i++) k[i] = strtoname("k" + strofreal(M[i,1], "%16.0g") + "_" + strofreal(M[i,2], "%16.0g") + "_" + strofreal(M[i,4], "%16.0g"))
  return(k)
}
string colvector wooldid_bootcellkey(string scalar mname, real scalar isgt)
{
  real matrix M
  string colvector k
  string matrix rs
  real scalar i
  M = st_matrix(mname)
  if (isgt) k = wooldid_bootcellkeym(M)
  else {
    rs = st_matrixrowstripe(mname)
    k = J(rows(rs), 1, "")
    for (i = 1; i <= rows(rs); i++) k[i] = strtoname("r" + rs[i, 1] + "_" + rs[i, 2])
  }
  return(k)
}
void wooldid_bootcellkeys(string scalar mname, string scalar accname, real scalar reps, real scalar isgt)
{
  string colvector k
  real matrix A
  k = wooldid_bootcellkey(mname, isgt)
  st_local("wdb_cellkeydup", "")
  if (rows(uniqrows(k)) != rows(k)) {
    // sanitised keys collided: two cells would share a column, so this surface is not
    // bootstrapped rather than bootstrapped wrongly. The caller cautions
    st_local("wdb_cellkeydup", "1")
    return
  }
  A = J(reps, rows(k), .)
  st_matrix(accname, A)
  st_matrixcolstripe(accname, (J(rows(k), 1, ""), k))
}
// per-replication harvest: gteffects estimate (col 5) / se (col 6); cfx b (col 1) / se (col 2)
// mapfr/ivar (optional): under the physical re-id leg the replication ran on fresh unit
// codes (__wdbnewi = group(resampled cluster, i)), so its e(gteffects) rows are keyed on
// codes that mean nothing to the observed table - matching them by code harvested the
// wrong cells when the original codes were 1..N and nothing otherwise (2026-09-05 review).
// The frame that ran the replication still holds both columns: translate each row's unit
// code back to the ORIGINAL i before keying. A twice-drawn cluster yields two copies
// with identical cell coefficients; the first copy is kept. Cells of undrawn units stay
// missing for the replication
void wooldid_bootcells(string scalar prefix, real scalar r, real scalar ncfx, real scalar dose, | string scalar mapfr, string scalar ivar)
{
  real matrix A, S, M
  string colvector ak, mk
  string matrix cs
  real scalar j, k, c, remap
  real colvector nv, ov
  transmorphic t, tm
  string scalar nm, cur
  remap = (args() >= 6 ? mapfr != "" : 0)
  // gteffects
  M = st_matrix("e(gteffects)")
  if (rows(M) > 0 & st_matrix(prefix + "_gt") != J(0,0,.)) {
    A = st_matrix(prefix + "_gt")
    cs = st_matrixcolstripe(prefix + "_gt")
    if (remap) {
      cur = st_framecurrent()
      st_framecurrent(mapfr)
      nv = st_data(., "__wdbnewi")
      ov = st_data(., ivar)
      st_framecurrent(cur)
      tm = asarray_create("real")
      for (j = 1; j <= rows(nv); j++) asarray(tm, nv[j], ov[j])
      for (j = 1; j <= rows(M); j++) {
        if (asarray_contains(tm, M[j, 1])) M[j, 1] = asarray(tm, M[j, 1])
      }
    }
    mk = wooldid_bootcellkeym(M)
    t = asarray_create()
    for (j = 1; j <= rows(mk); j++) {
      if (!asarray_contains(t, mk[j])) asarray(t, mk[j], j)
    }
    if (dose) S = st_matrix(prefix + "_gtse")
    for (k = 1; k <= rows(cs); k++) {
      if (!asarray_contains(t, cs[k, 2])) continue
      j = asarray(t, cs[k, 2])
      A[r, k] = M[j, 5]
      if (dose) {
        if (M[j, 6] > 0 & M[j, 6] < .) S[r, k] = M[j, 6]
      }
    }
    st_matrix(prefix + "_gt", A)
    st_matrixcolstripe(prefix + "_gt", cs)
    if (dose) {
      st_matrix(prefix + "_gtse", S)
      st_matrixcolstripe(prefix + "_gtse", cs)
    }
  }
  // cfx matrices, by the names recorded in the observed e(wdid_cfxmats) order
  for (c = 1; c <= ncfx; c++) {
    nm = st_local("wdb_cfxname" + strofreal(c))
    if (nm == "") continue
    M = st_matrix("e(" + nm + ")")
    if (rows(M) == 0) continue
    A = st_matrix(prefix + "_cfx" + strofreal(c))
    if (A == J(0,0,.)) continue
    cs = st_matrixcolstripe(prefix + "_cfx" + strofreal(c))
    mk = wooldid_bootcellkey("e(" + nm + ")", 0)
    t = asarray_create()
    for (j = 1; j <= rows(mk); j++) asarray(t, mk[j], j)
    if (dose) S = st_matrix(prefix + "_cfx" + strofreal(c) + "se")
    for (k = 1; k <= rows(cs); k++) {
      if (!asarray_contains(t, cs[k, 2])) continue
      j = asarray(t, cs[k, 2])
      A[r, k] = M[j, 1]
      if (dose & cols(M) >= 2) {
        if (M[j, 2] > 0 & M[j, 2] < .) S[r, k] = M[j, 2]
      }
    }
    st_matrix(prefix + "_cfx" + strofreal(c), A)
    st_matrixcolstripe(prefix + "_cfx" + strofreal(c), cs)
    if (dose) {
      st_matrix(prefix + "_cfx" + strofreal(c) + "se", S)
      st_matrixcolstripe(prefix + "_cfx" + strofreal(c) + "se", cs)
    }
  }
}
// rewire one cell table in place: given the accumulator (and optional se accumulator +
// observed analytic se), fill columns (secol, pcol, lbcol, ubcol) of matrix mname per
// bootdisplay(): percentile / bc from the cell draws with the bootstrap SD as se; wald;
// or, when studentized (seacc given), the percentile-t family with the analytic se
// (bootdisplay(bc) is refused with boottstat - F003 - so fam is "percentile" there).
// obscol is the column of mname holding the observed estimate. Rows whose key is not in
// the accumulator, or with < 2 usable draws, are blanked. Returns nothing; edits mname.
// bootsave(): fill the (empty) current frame from the accumulators. Row 1 is the observed
// run (_rep = 0), rows 2.. the replications by index (failed ones all missing). Variable
// names are the accumulator column names made valid with strtoname(); a name collision
// or an over-long name makes st_addvar fail, which the caller reports
void wooldid_bootsavemk(string scalar acc, string scalar obsb, string scalar seacc, string scalar aseobs, string scalar pwseacc,
                        | string scalar wmat, real scalar wobs, string scalar jwacc, string scalar jobs, string scalar awacc, string scalar aobs)
{
  real matrix A, B, S, AS, P, W, JW, JO, AW, AO
  string matrix cs, pcs, wcs, jcs
  real scalar R, K, k, j, idx
  string scalar nm
  A = st_matrix(acc); cs = st_matrixcolstripe(acc); B = st_matrix(obsb)
  R = rows(A); K = cols(A)
  st_addobs(R + 1)
  idx = st_addvar("double", "_rep"); st_store(., idx, (0::R))
  for (k = 1; k <= K; k++) {
    nm = strtoname(cs[k, 2])
    idx = st_addvar("double", nm); st_store(., idx, B[1, k] \ A[., k])
  }
  if (seacc != "") {
    S = st_matrix(seacc); AS = st_matrix(aseobs)
    for (k = 1; k <= K; k++) {
      nm = substr("se_" + strtoname(cs[k, 2]), 1, 32)
      idx = st_addvar("double", nm); st_store(., idx, AS[1, k] \ S[., k])
    }
  }
  if (pwseacc != "") {
    P = st_matrix(pwseacc); pcs = st_matrixcolstripe(pwseacc)
    for (j = 1; j <= cols(P); j++) {
      // the pairwise names repeat across the main and pre-treatment rows ("1vs0.sg"
      // twice), so the column index is part of the name: pwse<j>_ matches column j
      // of e(bootdist_pwse)
      // ... and the index makes a truncation to Stata's 32-character limit collision-
      // free (a long subgroup variable name pushed the full name past it: r(3300))
      nm = substr("pwse" + strofreal(j) + "_" + strtoname(pcs[j, 2]), 1, 32)
      idx = st_addvar("double", nm); st_store(., idx, . \ P[., j])
    }
  }
  // F063: the resample - w_<cluster> multiplicities (standard) or Dirichlet multipliers
  // (bayes); the observed row carries wobs (1, or 1/C under bayes where the multipliers
  // sum to 1 per replication). Then the per-replication joint and adversarial statistics
  // under boottstat, with the observed statistic in row 0
  if (args() >= 6) {
    if (wmat != "") {
      W = st_matrix(wmat); wcs = st_matrixcolstripe(wmat)
      for (k = 1; k <= cols(W); k++) {
        nm = substr(strtoname(wcs[k, 2]), 1, 32)
        if (_st_varindex(nm) < .) continue
        idx = st_addvar("double", nm); st_store(., idx, wobs \ W[., k])
      }
    }
    if (jwacc != "") {
      JW = st_matrix(jwacc); jcs = st_matrixcolstripe(jwacc); JO = (jobs != "" ? st_matrix(jobs) : J(1, cols(JW), .))
      for (k = 1; k <= cols(JW); k++) {
        nm = substr("jointW_" + strtoname(jcs[k, 2]), 1, 32)
        if (_st_varindex(nm) < .) continue
        idx = st_addvar("double", nm); st_store(., idx, (cols(JO) >= k ? JO[1, k] : .) \ JW[., k])
      }
    }
    if (awacc != "") {
      AW = st_matrix(awacc); AO = (aobs != "" ? st_matrix(aobs) : J(1, 1, .))
      for (k = 1; k <= cols(AW); k++) {
        nm = substr("advW" + (cols(AW) > 1 ? strofreal(k) : ""), 1, 32)
        if (_st_varindex(nm) < .) continue
        idx = st_addvar("double", nm); st_store(., idx, (cols(AO) >= k ? AO[1, k] : .) \ AW[., k])
      }
    }
  }
}

// F063 helpers. bootwnames: name the resample matrix's columns w_<cluster value> from the
// map frame (called inside it; rows are in map-frame order). bootwrow: row r from the map
// frame's __wdbdraw (count and bayes legs). bootwrowx: row r from an expanded resample
// frame - the multiplicity of each ORIGINAL cluster is the number of distinct resampled
// cluster ids (idcluster) it received; clusters not drawn get 0
void wooldid_bootwnames(string scalar wmat, string scalar clusv)
{
  real colvector cv
  real scalar c
  string matrix cs
  cv = st_data(., clusv)
  cs = J(rows(cv), 2, "")
  for (c = 1; c <= rows(cv); c++) cs[c, 2] = strtoname("w_" + strofreal(cv[c], "%20.0g"))
  st_matrixcolstripe(wmat, cs)
}
void wooldid_bootwrow(string scalar wmat, real scalar r, string scalar drawv)
{
  real matrix W
  string matrix cs
  W = st_matrix(wmat)
  cs = st_matrixcolstripe(wmat)
  W[r, .] = st_data(., drawv)'
  // st_matrix() resets the stripes; put the w_<cluster> names back
  st_matrix(wmat, W)
  st_matrixcolstripe(wmat, cs)
}
void wooldid_bootwrowx(string scalar wmat, real scalar r, string scalar clusv, string scalar nsv)
{
  real matrix W, D, U
  string matrix cs
  real scalar c, k
  transmorphic t
  W = st_matrix(wmat)
  cs = st_matrixcolstripe(wmat)
  t = asarray_create()
  asarray_notfound(t, 0)
  for (k = 1; k <= rows(cs); k++) asarray(t, cs[k, 2], k)
  D = uniqrows(st_data(., (clusv, nsv)))
  W[r, .] = J(1, cols(W), 0)
  for (c = 1; c <= rows(D); c++) {
    k = asarray(t, strtoname("w_" + strofreal(D[c, 1], "%20.0g")))
    if (k > 0) W[r, k] = W[r, k] + 1
  }
  st_matrix(wmat, W)
  st_matrixcolstripe(wmat, cs)
}

void wooldid_bootcellwire(string scalar mname, real scalar isgt, string scalar accname,
                          string scalar seaccname, real scalar obscol, real scalar secol,
                          real scalar pcol, real scalar lbcol, real scalar ubcol,
                          string scalar fam, real scalar lvl, real scalar zc, real scalar dfden, real scalar ties)
{
  real matrix M, A, S
  string colvector mk
  string matrix cs
  real colvector v, use, tr
  real scalar i, k, n, al, x, se, m0, m1, z0, pl, pu, kk2, r0, qbc, b, so, to, lo, hi, p, miss
  transmorphic t
  miss = 0
  M = st_matrix(mname)
  A = st_matrix(accname)
  cs = st_matrixcolstripe(accname)
  t = asarray_create()
  for (k = 1; k <= rows(cs); k++) asarray(t, cs[k, 2], k)
  mk = wooldid_bootcellkey(mname, isgt)
  if (seaccname != "") S = st_matrix(seaccname)
  al = (100 - lvl)/200
  for (i = 1; i <= rows(M); i++) {
    M[i, secol] = .
    M[i, pcol] = .
    M[i, lbcol] = .
    M[i, ubcol] = .
    if (!asarray_contains(t, mk[i])) continue
    k = asarray(t, mk[i])
    b = M[i, obscol]
    if (b >= .) continue
    if (seaccname == "") {
      use = selectindex(A[., k] :< .)
      n = length(use)
      if (n < 2) continue
      v = sort(A[use, k], 1)
      se = sqrt(variance(v))
      if (se <= 0 | se >= .) continue
      M[i, secol] = se
      if (fam == "wald") {
        to = b / se
        M[i, pcol] = (dfden < . ? 2*ttail(dfden, abs(to)) : 2*normal(-abs(to)))
        M[i, lbcol] = b - zc*se
        M[i, ubcol] = b + zc*se
        continue
      }
      x = n*al
      lo = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))])
      x = n*(1-al)
      hi = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))])
      m0 = sum(v :<= 0)
      m1 = sum(v :>= 0)
      p = min((1, (2*min((m0, m1)))/n))
      if (fam == "bc") {
        // z0 ranks the DISPLAYED cell estimate b (the coefficient families rank the rep-0
        // comparator). Deliberate and measured immaterial (F018 research record,
        // 2026-09-06): gteffects cells are stage-1 e(b) entries that fastnose cannot move,
        // so b equals the comparator bitwise; cfx cells differ from it only in margins'
        // last digits (~2e-9 bootstrap SDs); 0 of 4,425 cells had a different z0. The
        // amevoid case is caught by the b >= . guard above
        if (ties) z0 = invnormal((sum(v :< b) + 0.5*sum(v :== b))/n)
        else z0 = invnormal(sum(v :<= b)/n)
        if (z0 >= .) continue
        pl = normal(2*z0 + invnormal(al))
        pu = normal(2*z0 + invnormal(1-al))
        lo = wooldid_bootqs(v, pl)
        hi = wooldid_bootqs(v, pu)
        r0 = m0/n
        if (r0 <= 0 | r0 >= 1) qbc = (r0 <= 0 ? 0 : 1)
        else qbc = normal(invnormal(r0) - 2*z0)
        p = min((1, 2*min((qbc, 1-qbc))))
      }
      M[i, pcol] = p
      M[i, lbcol] = lo
      M[i, ubcol] = hi
    }
    else {
      // studentized: analytic se of the OBSERVED cell is the value already in secol before
      // we blanked it - so it must be passed in via the matrix's own column; re-read it
      so = st_matrix(mname)[i, secol]
      if (so >= . | so <= 0) {
        miss++
        continue
      }
      to = b / so
      use = selectindex((A[., k] :< .) :& (S[., k] :< .) :& (S[., k] :> 0))
      n = length(use)
      if (n < 2) continue
      tr = (A[use, k] :- b) :/ S[use, k]
      v = sort(tr, 1)
      if (sqrt(variance(v)) <= 0) continue
      M[i, secol] = so
      x = n*al
      lo = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))])
      x = n*(1-al)
      hi = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))])
      m0 = sum(v :<= to)
      m1 = sum(v :>= to)
      // F033 (author ruling 2026-09-07, propagated here on the 2026-09-11 ruling): these
      // draws ARE recentred - tr = (b_r - b)/se_r ranked against to = b/se_obs - which is
      // the property F033 said earns the add-one Monte-Carlo-null form, so the studentized
      // CELL p-value uses it like wooldid_boottp and wooldid_bootpwt. One convention for
      // every recentred rank-based p in the package. It used to be the plain 2*min/n share,
      // which put a cell p and the coefficient p for the same quantity roughly 2/(n+1)
      // apart - about 0.02 at reps(99), the size defect F033 was measured to remove - and
      // let a cell p print as exactly 0 where the coefficient p could not. The floor is now
      // 2/(n+1). (The non-studentized cell branch above correctly has NO add-one: those
      // draws are not recentred.)
      p = min((1, 2*(min((m0, m1)) + 1)/(n + 1)))
      M[i, pcol] = p
      M[i, lbcol] = b - hi*so
      M[i, ubcol] = b - lo*so
    }
  }
  st_matrix(mname, M)
  st_local("wdbcellmiss", strofreal(miss))
}

// ---- bootstrapped ADVERSARIAL line-test statistic (task 6, author ruling 2026-09-04) ----
// Under boottstat the adversarial p is the bootstrap of its test STATISTIC: per replication
// W_r = wooldid_advW(b_r - b_obs, V_r, rel) with the replication's own e(V_es_*), against
// W_obs = wooldid_advW(b_obs, V_obs_analytic, rel); p = add-one tail. Recentring by b_obs
// alone imposes the composite null "lies on SOME line" because the minimum-distance
// statistic is invariant to adding any element of col(X) - VERIFIED numerically
// 2026-09-05 (max |W(b) - W(b + Xc)| = 3.6e-14 full-rank, 1.2e-14 rank-deficient; recentre
// by b_obs equals recentre to the fitted line exactly). The forecast (altmethod) statistic
// is linear in b and recentres the same way. rel is derived from the horizon names
// (pre# -> -#, post# -> +#), which is the F2 capture's own convention.
// aobs: 1 x P observed coefficients over ALL horizons with colstripe (fk, name); accname:
// R x (keys x alts) with colnames <fk>_a0 / <fk>_a1
real colvector wooldid_bootadvrel(string colvector names)
{
  real colvector rel
  real scalar j
  rel = J(rows(names), 1, .)
  for (j = 1; j <= rows(names); j++) {
    if (substr(names[j], 1, 3) == "pre") rel[j] = -strtoreal(substr(names[j], 4, .))
    else if (substr(names[j], 1, 4) == "post") rel[j] = strtoreal(substr(names[j], 5, .))
  }
  return(rel)
}
void wooldid_bootadvr(string scalar accname, real scalar r, string scalar aobsn,
                      real scalar esrel, real scalar tt)
{
  real matrix Wacc, V, Vs
  real rowvector bo
  real colvector d, rel, idx
  string matrix wn, on, vn
  string colvector nm
  string scalar fk, key
  real scalar k, j, np, m, alt, W, q
  transmorphic t
  Wacc = st_matrix(accname)
  wn = st_matrixcolstripe(accname)
  on = st_matrixcolstripe(aobsn)
  bo = st_matrix(aobsn)
  for (k = 1; k <= rows(wn); k++) {
    key = wn[k, 2]
    alt = (substr(key, -3, .) == "_a1")
    fk = substr(key, 1, strlen(key) - 3)
    V = st_matrix("e(V_es_" + fk + ")")
    if (rows(V) == 0) continue
    vn = st_matrixcolstripe("e(V_es_" + fk + ")")
    t = asarray_create()
    for (j = 1; j <= rows(vn); j++) asarray(t, vn[j, 2], j)
    np = 0
    for (j = 1; j <= rows(on); j++) if (on[j, 1] == fk) np++
    if (np == 0) continue
    d = J(np, 1, .)
    idx = J(np, 1, .)
    nm = J(np, 1, "")
    m = 0
    for (j = 1; j <= rows(on); j++) {
      if (on[j, 1] != fk) continue
      m++
      if (!asarray_contains(t, on[j, 2])) {
        m = -1
        break
      }
      idx[m] = asarray(t, on[j, 2])
      nm[m] = on[j, 2]
      d[m] = st_matrix("e(b_es_" + fk + ")")[1, idx[m]] - bo[1, j]
    }
    if (m < 0 | hasmissing(d)) continue
    Vs = V[idx, idx]
    rel = wooldid_bootadvrel(nm)
    if (hasmissing(rel)) continue
    q = .
    W = wooldid_advW(d, Vs, rel, esrel, tt, alt, q)
    if (W < .) Wacc[r, k] = W
  }
  st_matrix(accname, Wacc)
  st_matrixcolstripe(accname, wn)
}
void wooldid_bootadvobs(string scalar fk, real scalar alt, string scalar aobsn,
                        real scalar esrel, real scalar tt, string scalar outloc)
{
  real matrix V, Vs
  real rowvector bo
  real colvector b, rel, idx
  string matrix on, vn
  string colvector nm
  string scalar vname
  real scalar j, np, m, W, q
  transmorphic t
  st_local(outloc, ".")
  vname = "e(V_es_analytic_" + fk + ")"
  V = st_matrix(vname)
  if (rows(V) == 0) {
    vname = "e(V_es_" + fk + ")"
    V = st_matrix(vname)
  }
  if (rows(V) == 0) return
  vn = st_matrixcolstripe(vname)
  t = asarray_create()
  for (j = 1; j <= rows(vn); j++) asarray(t, vn[j, 2], j)
  on = st_matrixcolstripe(aobsn)
  bo = st_matrix(aobsn)
  np = 0
  for (j = 1; j <= rows(on); j++) if (on[j, 1] == fk) np++
  if (np == 0) return
  b = J(np, 1, .)
  idx = J(np, 1, .)
  nm = J(np, 1, "")
  m = 0
  for (j = 1; j <= rows(on); j++) {
    if (on[j, 1] != fk) continue
    m++
    if (!asarray_contains(t, on[j, 2])) return
    idx[m] = asarray(t, on[j, 2])
    nm[m] = on[j, 2]
    b[m] = bo[1, j]
  }
  Vs = V[idx, idx]
  rel = wooldid_bootadvrel(nm)
  if (hasmissing(rel)) return
  q = .
  W = wooldid_advW(b, Vs, rel, esrel, tt, alt, q)
  if (W < .) st_local(outloc, strofreal(W, "%21x"))
}

// ---- bootstrapped JOINT statistics (task 5, author ruling 2026-09-03/04) ----
// A joint test under boottstat is the bootstrap of its TEST STATISTIC: per replication
// the recentred Wald W_r = (b_r - b_obs)' V_r^+ (b_r - b_obs) on the ES pre-horizon block,
// using the replication's OWN e(V_es_*) (a genuine full covariance among horizons - the
// posted e(V) is diagonal, but the F2 event-study capture is not), against
// W_obs = b_obs' V_obs^+ b_obs from the observed ANALYTIC e(V_es_*). The add-one tail
// (1 + #{W_r >= W_obs})/(n + 1) is the Monte-Carlo-null convention: unlike the percentile
// p on raw coefficient draws, this IS a reference distribution simulated under the null.
// The scale-free eigen guard is wooldid_bootjw's, verbatim in spirit: correlation-
// normalise, reject on relative rank deficiency, undo the scaling on the inverse.

// recentred per-replication statistic, called INSIDE the rep from wooldid_bootgrab where
// the rep's e() is live. obsn: 1 x P row of observed pre-coefficients with colstripe
// (fk, prename); accname: R x nkeys accumulator with colnames = fk. Name-anchored: a rep
// missing any observed pre-name for a family is dropped for that family (stays missing).
void wooldid_bootjwr(string scalar accname, real scalar r, string scalar obsn)
{
  real matrix Wacc, V, R, Si
  real rowvector bo, ev, sc
  real colvector d
  string matrix wn, on, vn
  string scalar fk
  real scalar k, j, m, np, q, W
  transmorphic t
  Wacc = st_matrix(accname)
  wn = st_matrixcolstripe(accname)
  on = st_matrixcolstripe(obsn)
  bo = st_matrix(obsn)
  for (k = 1; k <= rows(wn); k++) {
    fk = wn[k, 2]
    V = st_matrix("e(V_es_" + fk + ")")
    if (rows(V) == 0) continue
    vn = st_matrixcolstripe("e(V_es_" + fk + ")")
    t = asarray_create()
    for (j = 1; j <= rows(vn); j++) asarray(t, vn[j, 2], j)
    // the observed pre-names of this family, in observed order
    np = 0
    for (j = 1; j <= rows(on); j++) if (on[j, 1] == fk) np++
    if (np == 0) continue
    d = J(np, 1, .)
    R = J(np, np, .)
    m = 0
    for (j = 1; j <= rows(on); j++) {
      if (on[j, 1] != fk) continue
      m++
      if (!asarray_contains(t, on[j, 2])) {
        m = -1
        break
      }
      d[m] = st_matrix("e(b_es_" + fk + ")")[1, asarray(t, on[j, 2])] - bo[1, j]
    }
    if (m < 0 | hasmissing(d)) continue
    // V_r on the same names, same order. A pre-fill of R with missings used to walk these
    // names first and was DELETED as dead (D-5, author ruling 2026-09-11): the block below
    // overwrites R wholesale with V[idx, idx]. A missing name is caught by the
    // `m = -1; break' path above, and the hasmissing() guard below tests V's own entries.
    // fill R by name lookups
    {
      real colvector idx
      idx = J(np, 1, .)
      q = 0
      for (j = 1; j <= rows(on); j++) {
        if (on[j, 1] != fk) continue
        q++
        idx[q] = asarray(t, on[j, 2])
      }
      R = V[idx, idx]
    }
    if (hasmissing(R)) continue
    R = (R + R')/2
    sc = sqrt(diagonal(R))'
    if (min(sc) <= 0 | missing(sc)) continue
    R = R :/ (sc' * sc)
    ev = symeigenvalues(R)
    if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) continue
    Si = invsym(R) :/ (sc' * sc)
    W = d' * Si * d
    if (W < 0 | W >= .) continue
    Wacc[r, k] = W
  }
  st_matrix(accname, Wacc)
  st_matrixcolstripe(accname, wn)
}

// observed comparator for one family: W_obs = b' V^+ b on the pre block, from the OBSERVED
// analytic covariance (V_es_analytic_<fk> when bootwire replaced the live one, else V_es_<fk>).
// Returns W_obs in the local, "." when the guard rejects.
void wooldid_bootjwobs(string scalar fk, string scalar obsn, string scalar outloc)
{
  real matrix V, R, Si
  real rowvector bo, ev, sc
  real colvector b, idx
  string matrix on, vn
  string scalar vname
  real scalar j, np, q, W
  transmorphic t
  st_local(outloc, ".")
  vname = "e(V_es_analytic_" + fk + ")"
  V = st_matrix(vname)
  if (rows(V) == 0) {
    vname = "e(V_es_" + fk + ")"
    V = st_matrix(vname)
  }
  if (rows(V) == 0) return
  vn = st_matrixcolstripe(vname)
  t = asarray_create()
  for (j = 1; j <= rows(vn); j++) asarray(t, vn[j, 2], j)
  on = st_matrixcolstripe(obsn)
  bo = st_matrix(obsn)
  np = 0
  for (j = 1; j <= rows(on); j++) if (on[j, 1] == fk) np++
  if (np == 0) return
  b = J(np, 1, .)
  idx = J(np, 1, .)
  q = 0
  for (j = 1; j <= rows(on); j++) {
    if (on[j, 1] != fk) continue
    q++
    if (!asarray_contains(t, on[j, 2])) return
    idx[q] = asarray(t, on[j, 2])
    b[q] = bo[1, j]
  }
  R = V[idx, idx]
  if (hasmissing(R) | hasmissing(b)) return
  R = (R + R')/2
  sc = sqrt(diagonal(R))'
  if (min(sc) <= 0 | missing(sc)) return
  R = R :/ (sc' * sc)
  ev = symeigenvalues(R)
  if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) return
  Si = invsym(R) :/ (sc' * sc)
  W = b' * Si * b
  if (W < 0 | W >= .) return
  st_local(outloc, strofreal(W, "%21x"))
}

// across-subgroup equality under boottstat, by POLARISATION. For the m subgroup LEVEL
// columns of one estimand (cols) the per-replication covariance among them is rebuilt from
// what the replication actually produced: variances from its level standard errors
// (seacc, the e(V) diagonal) and off-diagonals from its pairwise standard errors
// (pwseacc, margins' pairwise se) via Cab = (Va + Vb - Var(a-b))/2 - exact when all three
// come off one margins pass, which they do (verified 2026-09-03 on an asymmetric fixture:
// three different triples gave one identical covariance and a PSD matrix). W_r is the
// chain-contrast Wald of (b_r - b_obs); W_obs the same from the OBSERVED analytic level se
// (aseobs) and observed pairwise se (pwobs). pairidx: m x m matrix of pwseacc/pwobs column
// indices for the (a,b) pair (0 on the diagonal); symmetric.
// THE INDEX INVARIANT, which the two headers describe in incompatible vocabularies
// (wooldid_bootpairidx says it returns "e(pwsgs) row indices", this one calls its argument
// "pwseacc/pwobs column indices"). Both are true only because those are the same integer by
// construction: column j of pwseacc == column j of pwobs0 == row j of the OBSERVED
// e(pwsgs). Three separate places keep it so - wooldid_bootrun creates the accumulator as
// J(reps, npw, .) with colnames taken from e(pwsgs)'s row names; wooldid_bootjnt snapshots
// e(pwsgs) before anything rewrites it and transposes column 2 into pwobs0; and
// wooldid_bootpairidx re-reads the LIVE e(pwsgs) row stripe. Reorder, append to, or repost
// e(pwsgs) between accumulator creation and this call and the wrong pairwise standard error
// is paired into Cab = (Va + Vb - Var(a-b))/2, returning a plausible and silently wrong p
// for joint_<item>_allsgssame. No guard can catch it: every index stays in range.
void wooldid_bootsgw(string scalar accname, string scalar seaccname, string scalar pwseaccname,
                     string scalar obsbname, string scalar aseobsname, string scalar pwobsname,
                     string scalar colslist, string scalar pairname, string scalar outloc,
                     string scalar nloc)
{
  real matrix A, S, P, Sr, Cm, R, Si, pairidx
  real rowvector idx, ob, ase, pwo, sc, ev
  real colvector d, Wr, b, sr
  real scalar m, r, a, bb, j, W, Wobs, n, nex, ok
  st_local(outloc, ".")
  st_local(nloc, "0")
  idx = strtoreal(tokens(colslist))
  m = cols(idx)
  if (m < 2) return
  A = st_matrix(accname)[., idx']
  S = st_matrix(seaccname)[., idx']
  P = st_matrix(pwseaccname)
  ob = st_matrix(obsbname)[1, idx']
  ase = st_matrix(aseobsname)[1, idx']
  pwo = st_matrix(pwobsname)
  pairidx = st_matrix(pairname)
  Cm = J(m-1, m, 0)
  for (j = 1; j <= m-1; j++) {
    Cm[j, j] = 1
    Cm[j, j+1] = -1
  }
  // W_obs
  Sr = J(m, m, .)
  for (a = 1; a <= m; a++) {
    Sr[a, a] = ase[a]^2
    for (bb = a+1; bb <= m; bb++) {
      j = pairidx[a, bb]
      if (j < 1 | j >= .) return
      Sr[a, bb] = (ase[a]^2 + ase[bb]^2 - pwo[1, j]^2)/2
      Sr[bb, a] = Sr[a, bb]
    }
  }
  if (hasmissing(Sr) | hasmissing(ob) | hasmissing(ase)) return
  R = Cm * Sr * Cm'
  R = (R + R')/2
  sc = sqrt(diagonal(R))'
  if (min(sc) <= 0 | missing(sc)) return
  R = R :/ (sc' * sc)
  ev = symeigenvalues(R)
  if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) return
  Si = invsym(R) :/ (sc' * sc)
  d = Cm * ob'
  Wobs = d' * Si * d
  if (Wobs < 0 | Wobs >= .) return
  // W_r per replication
  Wr = J(0, 1, .)
  for (r = 1; r <= rows(A); r++) {
    b = A[r, .]'
    sr = S[r, .]'
    if (hasmissing(b) | hasmissing(sr)) continue
    if (min(sr) <= 0) continue
    ok = 1
    Sr = J(m, m, .)
    for (a = 1; a <= m; a++) {
      Sr[a, a] = sr[a]^2
      for (bb = a+1; bb <= m; bb++) {
        j = pairidx[a, bb]
        if (P[r, j] >= . | P[r, j] <= 0) {
          ok = 0
          break
        }
        Sr[a, bb] = (sr[a]^2 + sr[bb]^2 - P[r, j]^2)/2
        Sr[bb, a] = Sr[a, bb]
      }
      if (!ok) break
    }
    if (!ok) continue
    R = Cm * Sr * Cm'
    R = (R + R')/2
    sc = sqrt(diagonal(R))'
    if (min(sc) <= 0 | missing(sc)) continue
    R = R :/ (sc' * sc)
    ev = symeigenvalues(R)
    if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) continue
    Si = invsym(R) :/ (sc' * sc)
    d = Cm * (b - ob')
    W = d' * Si * d
    if (W < 0 | W >= .) continue
    Wr = Wr \ W
  }
  n = rows(Wr)
  st_local(nloc, strofreal(n))
  if (n < 2) return
  nex = sum(Wr :>= Wobs)
  st_local(outloc, strofreal((1 + nex)/(n + 1), "%21x"))
}

// studentized inference for ONE pairwise contrast: the contrast draws d_r = b_r[ka] -
// b_r[kb] over the replication's own contrast se (harvested from its e(pwsgs) column 2),
// referred against t_obs = d_obs/se_obs with the OBSERVED analytic pairwise se. Same
// construction as wooldid_boottp - equal-tailed signed quantiles from wooldid_bootci's
// order-statistic rule, add-one Monte-Carlo-null rank p (F033) - so a contrast row and a
// coefficient row are studentized by one rule (no BC-t twin: F003). Output: a 6 x 1 column with the same row
// layout as wooldid_boottp (n, t_obs, p_pct, qlo_pct, qhi_pct, live)
void wooldid_bootpwt(string scalar accname, real scalar ka, real scalar kb, string scalar pwseaccname,
                     real scalar jcol, real scalar dobs, real scalar seobs, real scalar lvl,
                     real scalar ties, string scalar outname)
{
  real matrix A, P
  real colvector ok, d, se, tr, v
  real scalar n, al, x, m0, m1, to, sd
  real colvector out
  out = J(6, 1, .)
  out[6] = 0
  st_matrix(outname, out)
  if (seobs >= . | seobs <= 0 | dobs >= .) return
  A = st_matrix(accname)[., (ka, kb)]
  P = st_matrix(pwseaccname)[., jcol]
  ok = (rowmissing(A) :== 0) :& (P :< .) :& (P :> 0)
  if (sum(ok) < 2) return
  d = select(A[., 1] :- A[., 2], ok)
  se = select(P, ok)
  tr = (d :- dobs) :/ se
  v = sort(tr, 1)
  n = rows(v)
  out[1] = n
  sd = sqrt(variance(v))
  if (sd <= 0 | sd >= .) {
    st_matrix(outname, out)
    return
  }
  to = dobs / seobs
  out[2] = to
  al = (100 - lvl)/200
  x = n*al
  if (x == floor(x) & x >= 1) out[4] = (v[x] + v[min((x+1, n))])/2
  else out[4] = v[max((min((ceil(x), n)), 1))]
  x = n*(1-al)
  if (x == floor(x) & x >= 1) out[5] = (v[x] + v[min((x+1, n))])/2
  else out[5] = v[max((min((ceil(x), n)), 1))]
  m0 = sum(v :<= to)
  m1 = sum(v :>= to)
  out[3] = min((1, 2*(min((m0, m1)) + 1)/(n + 1)))
  out[6] = 1
  st_matrix(outname, out)
}

void wooldid_bootpp(string scalar accname, real scalar ka, real scalar kb, real scalar lvl,
                    string scalar compname, real scalar ties,
                    string scalar lbn, string scalar ubn, string scalar pn, string scalar cnt,
                    string scalar blbn, string scalar bubn, string scalar bpn, string scalar z0n)
{
  real matrix A, cmp
  real colvector ok, d, v
  real scalar n, al, x, m0, m1, se, dcmp, z0, pl, pu, r0, qbc, kk2
  st_local(lbn, ".")
  st_local(ubn, ".")
  st_local(pn, ".")
  st_local(cnt, "0")
  st_local(blbn, ".")
  st_local(bubn, ".")
  st_local(bpn, ".")
  st_local(z0n, ".")
  A = st_matrix(accname)[., (ka, kb)]
  ok = rowmissing(A) :== 0
  if (sum(ok) < 2) return
  A = select(A, ok)
  d = A[., 1] :- A[., 2]
  v = sort(d, 1)
  n = rows(v)
  st_local(cnt, strofreal(n))
  se = sqrt(variance(v))
  if (se <= 0 | se >= .) return
  al = (100 - lvl)/200
  x = n*al
  if (x == floor(x) & x >= 1) st_local(lbn, strofreal((v[x] + v[min((x+1, n))])/2, "%21x"))
  else st_local(lbn, strofreal(v[max((min((ceil(x), n)), 1))], "%21x"))
  x = n*(1-al)
  if (x == floor(x) & x >= 1) st_local(ubn, strofreal((v[x] + v[min((x+1, n))])/2, "%21x"))
  else st_local(ubn, strofreal(v[max((min((ceil(x), n)), 1))], "%21x"))
  m0 = sum(v :<= 0)
  m1 = sum(v :>= 0)
  st_local(pn, strofreal(min((1, (2*min((m0, m1)))/n)), "%21x"))
  // the bias-corrected twin for this contrast. z0 is computed exactly as wooldid_bootci
  // computes it for a coefficient column - the share of contrast draws at or below the
  // REP-0 comparator's contrast, pushed through invnormal, with the same `ties'
  // half-weighting option - so a contrast row and a coefficient row are bias-corrected by
  // the same rule rather than by two rules that merely resemble each other
  cmp = st_matrix(compname)
  dcmp = cmp[1, ka] - cmp[1, kb]
  if (dcmp >= .) return
  if (ties) z0 = invnormal((sum(v :< dcmp) + 0.5*sum(v :== dcmp))/n)
  else z0 = invnormal(sum(v :<= dcmp)/n)
  // guard exactly as wooldid_bootci does: invnormal() of 0 or 1 returns MISSING in Stata,
  // so a single missing test covers both degenerate tails. (An earlier version wrote
  // `z0 <= -.' intending "minus infinity"; in Mata `-.' is missing and `x <= .' is true
  // for every real x, so that guard rejected every column and the BC contrast family came
  // out uniformly blank.)
  if (z0 >= .) return
  st_local(z0n, strofreal(z0, "%21x"))
  pl = normal(2*z0 + invnormal(al))
  pu = normal(2*z0 + invnormal(1-al))
  st_local(blbn, strofreal(wooldid_bootqs(v, pl), "%21x"))
  st_local(bubn, strofreal(wooldid_bootqs(v, pu), "%21x"))
  r0 = m0/n
  if (r0 <= 0 | r0 >= 1) qbc = (r0 <= 0 ? 0 : 1)
  else qbc = normal(invnormal(r0) - 2*z0)
  st_local(bpn, strofreal(min((1, 2*min((qbc, 1-qbc)))), "%21x"))
}

void wooldid_bootpw(string scalar accname, real scalar ka, real scalar kb)
{
  real matrix A, dev, S
  real colvector ok
  A = st_matrix(accname)[., (ka, kb)]
  ok = rowmissing(A) :== 0
  if (sum(ok) < 2) {
    st_local("wdbse", ".")
    return
  }
  A = select(A, ok)
  dev = A :- mean(A)
  S = (dev' * dev) / (rows(A) - 1)
  st_local("wdbse", strofreal(sqrt(S[1,1] + S[2,2] - 2*S[1,2]), "%21x"))
}

// order statistic at probability p under bstat's rule (integer n*p: the average of the
// adjacent order statistics; otherwise the ceil(n*p) order statistic) - the rule the
// percentile interval already used, now shared by every bias-corrected endpoint (the
// v[ceil(n*p)] rule they used disagreed with bstat whenever z0 = 0 and n*alpha was
// integral, i.e. at exactly half the draws below the comparator; 2026-09-05 review)
real scalar wooldid_bootqs(real colvector v, real scalar p)
{
  real scalar n, x
  n = rows(v)
  x = n*p
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1, n))])/2)
  return(v[max((min((ceil(x), n)), 1))])
}
// CI flavors per column: normal (obs +/- zc*se), percentile (classic Stata
// rule, pinned vs bstat at even AND odd rep counts), bias-corrected (order-statistic
// endpoints at p_adj under the same order-statistic rule as the percentile family,
// bstat's, so the two coincide at z0 = 0) with z0 computed
// against the REP-0 comparator; ties adjusts z0 by half-weighting exact ties.
// Columns with fewer than 2 available reps or zero dispersion post missing.
void wooldid_bootci(string scalar accname, string scalar obsbname, string scalar compname, real scalar lvl, real scalar ties, string scalar cinname, string scalar cipname, string scalar cibname, string scalar cippname, string scalar cbpname, string scalar z0name, real scalar dfden)
{
  real matrix A, cin, cip, cib, cipp, cbp
  real rowvector obs, cmp, z0v
  real colvector v, use
  real scalar K, kcol, al, zc, nn, se, z0, x, lo, fr, kk2, pl, pu, m0, m1, r0, qbc
  A = st_matrix(accname)
  obs = st_matrix(obsbname)
  cmp = st_matrix(compname)
  K = cols(A)
  al = (100 - lvl)/200
  zc = (dfden < . ? invt(dfden, 1 - al) : invnormal(1 - al))
  cin = J(2, K, .)
  cip = J(2, K, .)
  cib = J(2, K, .)
  cipp = J(1, K, .)
  cbp = J(1, K, .)
  z0v = J(1, K, .)
  for (kcol = 1; kcol <= K; kcol++) {
    use = selectindex(A[., kcol] :< .)
    nn = length(use)
    if (nn < 2) continue
    v = sort(A[use, kcol], 1)
    se = sqrt(variance(v))
    if (se <= 0 | se >= .) continue
    cin[1, kcol] = obs[kcol] - zc*se
    cin[2, kcol] = obs[kcol] + zc*se
    // percentile: the classic Stata rule (integer n*p: average adjacent order
    // statistics; otherwise the ceil(n*p) order statistic) - pinned vs bstat at
    // both even and odd replication counts
    x = nn*al
    if (x == floor(x) & x >= 1) cip[1, kcol] = (v[x] + v[min((x+1, nn))])/2
    else cip[1, kcol] = v[max((min((ceil(x), nn)), 1))]
    x = nn*(1-al)
    if (x == floor(x) & x >= 1) cip[2, kcol] = (v[x] + v[min((x+1, nn))])/2
    else cip[2, kcol] = v[max((min((ceil(x), nn)), 1))]
    // equal-tailed percentile p at the null theta = 0: the exact inversion of the
    // percentile interval built directly above, and the reason it is computed HERE
    // rather than alongside it - same selection, same nn, same guards, so the two
    // cannot drift apart. WEAK inequalities on both counts, which is what makes
    // "v[j] > 0 <=> j > m0" hold and with it the count-to-order-statistic duality;
    // a draw sitting exactly at the null therefore counts against neither side.
    // NO add-one, on the STATISTICAL argument only: the permutation engine's
    // (1+#)/(n+1) is the Monte-Carlo-null convention, where the observed statistic
    // legitimately joins a reference set drawn UNDER the null. These draws are
    // resampled under the SAMPLE, so there is no such set for it to join.
    //
    // Duality is near-exact but NOT exact, and it is worth being precise about why.
    // The interval rule directly above averages two adjacent order statistics when
    // nn*al is an integer, and no count-based p can be exactly dual to an averaged
    // endpoint, so a boundary disagreement is possible under any convention.
    // Measured 2026-09-03 over 35,646 constructed cases (levels 90/95/99 x nn in
    // {39,40,99,100} x the full enumeration of M draws at -1, Z draws exactly at the
    // null and the rest at +1, which puts the averaged endpoint on the null on
    // purpose): 6 disagreements with the convention below, against 36 for the
    // standard two-sided add-one 2*(1+min(m0,m1))/(nn+1). The weak inequalities are
    // what earn that margin - they are the reason a draw exactly at the null counts
    // against neither side, which is where add-one does its damage.
    //
    // Consequence to document: p can print as exactly 0, meaning no completed
    // replication fell on the far side of the null; the resolution is 2/nn. The p is
    // LEVEL-INDEPENDENT, so it is not recomputed when level() moves. It reads the
    // raw draws, so neither the rep-0 comparator nor ties() touches it; both of
    // those feed only the bias-corrected z0 below.
    m0 = sum(v :<= 0)
    m1 = sum(v :>= 0)
    cipp[1, kcol] = min((1, (2*min((m0, m1)))/nn))
    if (ties) z0 = invnormal((sum(v :< cmp[kcol]) + 0.5*sum(v :== cmp[kcol]))/nn)
    else z0 = invnormal(sum(v :<= cmp[kcol])/nn)
    if (z0 >= .) continue
    z0v[kcol] = z0
    pl = normal(2*z0 + invnormal(al))
    pu = normal(2*z0 + invnormal(1-al))
    cib[1, kcol] = wooldid_bootqs(v, pl)
    cib[2, kcol] = wooldid_bootqs(v, pu)
    // the BC p-value: the dual of the interval just built, obtained by solving for the
    // level at which a BC endpoint lands exactly on the null. The endpoint at one-sided
    // level q is v[ceil(nn*normal(2*z0 + invnormal(q)))], so it reaches 0 when
    // normal(2*z0 + invnormal(q)) equals r0, the rank of 0 among the draws - giving
    // q = normal(invnormal(r0) - 2*z0) and the two-sided p = 2*min(q, 1-q). With z0 = 0
    // this collapses to the plain percentile p, which is the right degenerate behaviour.
    // r0 uses the same weak inequality the percentile p uses, so a draw exactly at the
    // null counts the same way in both. r0 at 0 or 1 sends invnormal to -/+ infinity and
    // the p to 0, matching the percentile p's floor rather than inventing a different one
    r0 = sum(v :<= 0)/nn
    if (r0 <= 0 | r0 >= 1) qbc = (r0 <= 0 ? 0 : 1)
    else qbc = normal(invnormal(r0) - 2*z0)
    cbp[1, kcol] = min((1, 2*min((qbc, 1-qbc))))
  }
  st_matrix(cinname, cin)
  st_matrix(cipname, cip)
  st_matrix(cibname, cib)
  st_matrix(cippname, cipp)
  st_matrix(cbpname, cbp)
  st_matrix(z0name, z0v)
}


// ============================================================================
// F064 (2026-09-06): bootstrap tail record for the percentile family: per column, the
// shares of draws at or below / at or above 0 (P_lo, P_hi), the Monte-Carlo standard
// error of the percentile p (2*min rule: doubled) and its level() Clopper-Pearson interval
void wooldid_bootptails(string scalar accname, string scalar obsbname, real scalar lvl, string scalar outname)
{
  real matrix A, out
  real colvector v
  real scalar K, k, n, m0, m1, cc, qq, al2, cpl, cph
  A = st_matrix(accname)
  K = cols(A)
  out = J(5, K, .)
  al2 = (100 - lvl)/200
  for (k = 1; k <= K; k++) {
    v = select(A[., k], A[., k] :< .)
    n = rows(v)
    if (n < 2) continue
    m0 = sum(v :<= 0)
    m1 = sum(v :>= 0)
    out[1, k] = m0/n
    out[2, k] = m1/n
    cc = min((m0, m1))
    qq = cc/n
    out[3, k] = 2*sqrt(qq*(1 - qq)/n)
    cpl = (cc <= 0 ? 0 : invibeta(cc, n - cc + 1, al2))
    cph = (cc >= n ? 1 : invibeta(cc + 1, n - cc, 1 - al2))
    out[4, k] = min((1, 2*cpl))
    out[5, k] = min((1, 2*cph))
  }
  st_matrix(outname, out)
}

// F067 (2026-09-06): the BCa family. Per column: z0 as wooldid_bootci computes it (the
// rank of the rep-0 comparator among the draws, ties half-weighted on request), the
// acceleration a from the delete-one-cluster jackknife estimates theta_g (Efron 1987:
// a = sum(m - theta_g)^3 / (6 [sum(m - theta_g)^2]^1.5), m the jackknife mean; missing
// rows - failed fits - are dropped), the endpoints at the BCa-adjusted levels
// alpha_j = Phi(z0 + (z0 + z_j)/(1 - a (z0 + z_j))) by the same order-statistic rule
// (wooldid_bootqs), and the p that inverts the interval at the null: with r0 the rank
// of 0 among the draws, w = (Phi^-1(r0) - z0)/(1 + a (Phi^-1(r0) - z0)), q = Phi(w - z0),
// p = 2 min(q, 1 - q) (a = 0 gives the BC p exactly). Columns with fewer than 2 draws,
// zero dispersion, a degenerate z0 or fewer than 3 jackknife fits post missing
void wooldid_bootbca(string scalar accname, string scalar obsbname, string scalar compname,
                     string scalar jkname, real scalar lvl, real scalar ties,
                     string scalar ciname, string scalar pname, string scalar aname)
{
  real matrix A, JK, ci
  real rowvector cmp, p, av
  real colvector v, use, th
  real scalar K, k, nn, al, z0, a, d2, d3, w, q, r0, zr
  A = st_matrix(accname)
  JK = st_matrix(jkname)
  cmp = st_matrix(compname)
  K = cols(A)
  al = (100 - lvl)/200
  ci = J(2, K, .)
  p = J(1, K, .)
  av = J(1, K, .)
  for (k = 1; k <= K; k++) {
    use = selectindex(A[., k] :< .)
    nn = length(use)
    if (nn < 2) continue
    v = sort(A[use, k], 1)
    if (sqrt(variance(v)) <= 0) continue
    if (ties) z0 = invnormal((sum(v :< cmp[k]) + 0.5*sum(v :== cmp[k]))/nn)
    else z0 = invnormal(sum(v :<= cmp[k])/nn)
    if (z0 >= .) continue
    if (cols(JK) < k) continue
    th = select(JK[., k], JK[., k] :< .)
    if (rows(th) < 3) continue
    d2 = sum((mean(th) :- th):^2)
    d3 = sum((mean(th) :- th):^3)
    if (d2 <= 0) continue
    a = d3 / (6 * d2^1.5)
    w = z0 + invnormal(al)
    if (1 - a*w <= 0) continue
    ci[1, k] = wooldid_bootqs(v, normal(z0 + w/(1 - a*w)))
    w = z0 + invnormal(1 - al)
    if (1 - a*w <= 0) continue
    ci[2, k] = wooldid_bootqs(v, normal(z0 + w/(1 - a*w)))
    r0 = sum(v :<= 0)/nn
    if (r0 <= 0 | r0 >= 1) q = (r0 <= 0 ? 0 : 1)
    else {
      zr = invnormal(r0) - z0
      if (1 + a*zr <= 0) continue
      w = zr/(1 + a*zr)
      q = normal(w - z0)
    }
    p[k] = min((1, 2*min((q, 1 - q))))
    // recorded last, so e(boot_accel) is non-missing only where the whole family exists (BS-17)
    av[k] = a
  }
  st_matrix(ciname, ci)
  st_matrix(pname, p)
  st_matrix(aname, av)
}

// the BCa twin for ONE pairwise contrast: contrast draws d_r = b_r[ka] - b_r[kb], the
// comparator's contrast for z0, the jackknife contrasts for the acceleration
void wooldid_bootppa(string scalar accname, real scalar ka, real scalar kb, real scalar lvl,
                     string scalar compname, real scalar ties, string scalar jkname,
                     string scalar lbn, string scalar ubn, string scalar pn)
{
  real matrix A, cmp, JK
  real colvector ok, v, th
  real scalar n, al, z0, a, d2, d3, w, q, r0, zr, dcmp, lo, hi
  st_local(lbn, ".")
  st_local(ubn, ".")
  st_local(pn, ".")
  A = st_matrix(accname)[., (ka, kb)]
  ok = rowmissing(A) :== 0
  if (sum(ok) < 2) return
  A = select(A, ok)
  v = sort(A[., 1] :- A[., 2], 1)
  n = rows(v)
  if (sqrt(variance(v)) <= 0) return
  cmp = st_matrix(compname)
  dcmp = cmp[1, ka] - cmp[1, kb]
  if (dcmp >= .) return
  if (ties) z0 = invnormal((sum(v :< dcmp) + 0.5*sum(v :== dcmp))/n)
  else z0 = invnormal(sum(v :<= dcmp)/n)
  if (z0 >= .) return
  JK = st_matrix(jkname)
  if (cols(JK) < max((ka, kb))) return
  JK = JK[., (ka, kb)]
  ok = rowmissing(JK) :== 0
  if (sum(ok) < 3) return
  th = select(JK[., 1] :- JK[., 2], ok)
  d2 = sum((mean(th) :- th):^2)
  d3 = sum((mean(th) :- th):^3)
  if (d2 <= 0) return
  a = d3 / (6 * d2^1.5)
  al = (100 - lvl)/200
  w = z0 + invnormal(al)
  if (1 - a*w <= 0) return
  lo = wooldid_bootqs(v, normal(z0 + w/(1 - a*w)))
  w = z0 + invnormal(1 - al)
  if (1 - a*w <= 0) return
  hi = wooldid_bootqs(v, normal(z0 + w/(1 - a*w)))
  r0 = sum(v :<= 0)/n
  if (r0 <= 0 | r0 >= 1) q = (r0 <= 0 ? 0 : 1)
  else {
    zr = invnormal(r0) - z0
    if (1 + a*zr <= 0) return
    q = normal(zr/(1 + a*zr) - z0)
  }
  st_local(lbn, strofreal(lo, "%21x"))
  st_local(ubn, strofreal(hi, "%21x"))
  st_local(pn, strofreal(min((1, 2*min((q, 1 - q)))), "%21x"))
}

// Permutation-inference kernels (chunk 1, 2026-09-05). Every quantity here treats the
// replication draws as NULL draws: nothing is recentred by the observed value, and the
// observed statistic is ranked among them. Conventions (author rulings 2026-09-05):
//  - tail p: EQUAL-TAILED, p = min(1, 2*min(P_lo, P_hi)), P_lo = (1 + #{s_r <= s_obs})
//    /(n + 1), P_hi = (1 + #{s_r >= s_obs})/(n + 1) for sampled draws (the add-one
//    Monte-Carlo-null convention the engine has always used); plain shares #{.}/n when
//    the element is present in every one of the M enumerated assignments (the observed
//    assignment is then a member of its own reference set). One convention for every
//    surface, including the legacy e(perm_p_<est>) scalars (formerly |s|-based).
//  - interval: the SHIFT interval [b - q_hi, b - q_lo], q the equal-tailed quantiles
//    of the null draws by wooldid_bootci's order-statistic rule, on the statistic's
//    scale (times se_obs when the statistic is a t-ratio). It inverts the equal-tailed
//    test under the fixed-zero-null pivotality approximation (NOT Conley-Taber, whose
//    reference comes from placebo onsets among never-treated units only) - it is NOT a
//    sharp-null randomization-test inversion, which would re-permute adjusted outcomes
//    at every candidate value; the footer and the help carry that hedge.
//  - quadratic-form statistics (joint pre-trend, across-subgroup equality, adversarial)
//    are one-sided upper-tail on a FIXED metric: the permutation covariance of the
//    draws (replications plus the observed comparator row, listwise complete; under a
//    complete exact enumeration the enumerated rows alone, the observed assignment being
//    one of them), so W_obs and every W_r are exchangeable under the null. No per-replication covariance is
//    consumed, so the tests exist in both the coefficient and the t-ratio modes.
//  - permabsp (absp = 1): the absolute-value convention instead - p = #{|s_r| >= |s_obs|}
//    with the same add-one / plain-share rule - and the interval its symmetric dual,
//    [b - q, b + q] with q the level/100 quantile of |s_r| (the pre-rework convention of
//    the legacy e(perm_p_<est>) scalars). The null band and the off-centre guard stay on
//    the signed equal-tailed quantiles in both conventions.
// Row layout of the 13-row per-element kernel: 1 n, 2 s_obs, 3 p, 4 qlo, 5 qhi, 6 sd,
// 7 lb, 8 ub, 9 se_shift, 10 offcentre flag, 11 bandlo, 12 bandhi (the signed level()
// quantiles of the null draws), 13 live
// PRECONDITION, unstated until now: v must be SORTED ASCENDING and free of missing. This
// function never sorts and never checks, so an unsorted v returns a different element of
// the vector - silently, with no missing and no error. wooldid_permtail's select()+sort is
// what establishes the precondition; it is not a tidiness step. The rule is x = n*pr; if x
// is a positive integer, the average of order statistics x and x+1, else order statistic
// ceil(x), clamped to [1, n] - deliberately the SAME rule wooldid_bootci uses and that
// wooldid_boottp re-implements inline, so the permutation shift interval and the bootstrap
// percentile interval estimate the same quantile the same way. Three copies of one rule:
// they must move together or they decouple silently at small n.
real scalar wooldid_permq(real colvector v, real scalar pr)
{
  real scalar n, x
  n = rows(v)
  x = n*pr
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1, n))])/2)
  return(v[max((min((ceil(x), n)), 1))])
}

// BS-23 (author ruling 2026-09-07): is the OBSERVED assignment's own replication among
// the rows that form this element's reference set, and does it carry the observed
// statistic? Both halves matter: presence alone would let a present-but-unequal row
// report m = 0, i.e. p = 0, where even the add-one form floors at 1/(n+1). M >= . means
// "not an enumeration" (sampled draws keep the documented add-one convention even when a
// shuffle happens to redraw the observed assignment); orow = 0 means the generator had no
// identity column to point at. The bounds are tested in a separate statement ahead of any
// subscript, so this does not depend on whether & evaluates its right operand.
// Three ways this returns 0 for a reason other than non-membership, all deliberately
// fail-SAFE rather than fail-loud, and none signalled except through wdp_addone: (a) the
// `Am[orow, .] == cm' is a whole-matrix comparison, and Mata returns 0 for non-conformable
// operands rather than raising (measured, Stata 19), so a mis-shaped comparator silently
// falls back to the add-one convention; (b) the test is exact floating-point equality
// between numbers produced by two SEPARATE wooldid fits, so one ulp demotes the element -
// which is why wooldid_permrun holds the comparison statistics in tempname SCALARS rather
// than locals; (c) Am is overloaded, and the declared signature hides it: a column vector
// of statistics from wooldid_permstat/permcellwire, an R x 2 difference column from
// wooldid_permpw, or an R x m matrix from wooldid_permqf/permstep - cm is always one row of
// the same width. Every one of these degrades toward add-one, i.e. conservative, never
// anti-conservative.
real scalar wooldid_permmemb(real matrix Am, real rowvector cm, real scalar M,
                             real scalar orow)
{
  if (M >= . | orow < 1 | orow > rows(Am)) return(0)
  if (hasmissing(cm) | hasmissing(Am[orow, .])) return(0)
  return(Am[orow, .] == cm)
}

real rowvector wooldid_permtail(real colvector draws, real scalar sobs, real scalar exact,
                                real scalar memb, real scalar al, real scalar zc,
                                real scalar absp)
{
  // absp is the CONVENTION CODE (F052/F064, 2026-09-06): 0 equal-tailed, 1 absolute add-one
  // (permabsp), 2 raw c/n |s| share with no add-one (permrawp), 3 left one-sided, 4 right
  // one-sided (permside()). Rows 14-18 (appended; rows 1-13 keep their meaning): 14 P_lo
  // and 15 P_hi (the one-sided shares under the add-one / plain-share rule), 16 se_p (the
  // Monte-Carlo standard error of the reported p's underlying share, doubled under the
  // equal-tailed 2*min rule), 17-18 the Clopper-Pearson interval of that share at level()
  // (doubled and capped at 1 under equal-tailed)
  real colvector v
  real scalar n, qlo, qhi, m0, m1, plo, phi, qa, ma, cc, qq, al2, cpl, cph, mult
  real rowvector out
  out = J(1, 18, .)
  out[13] = 0
  v = select(draws, draws :< .)
  n = rows(v)
  out[1] = n
  if (n < 2 | sobs >= .) return(out)
  v = sort(v, 1)
  out[2] = sobs
  out[6] = (exact ? sqrt(variance(v)*(n-1)/n) : sqrt(variance(v)))
  qlo = wooldid_permq(v, al)
  qhi = wooldid_permq(v, 1 - al)
  out[11] = qlo
  out[12] = qhi
  m0 = sum(v :<= sobs)
  m1 = sum(v :>= sobs)
  ma = sum(abs(v) :>= abs(sobs))
  // BS-23 (author ruling 2026-09-07): `exact' is COMPLETENESS and still governs the
  // population-SD convention below. `memb' is MEMBERSHIP - the observed statistic is one
  // of the n reference values - which is what the SHARE convention turns on. The add-one
  // form exists to inject an observed statistic the reference set does not contain; where
  // the set already contains it, adding it again counts it twice and doubles the floor.
  // Callers set memb only under an enumeration, and only after checking the observed row's
  // value equals sobs.
  if (exact | memb) {
    plo = m0/n
    phi = m1/n
  }
  else {
    plo = (1 + m0)/(n + 1)
    phi = (1 + m1)/(n + 1)
  }
  out[14] = plo
  out[15] = phi
  mult = 1
  if (absp == 1 | absp == 2) {
    // absolute-value conventions: the level/100 quantile of |s_r| on each side
    qa = wooldid_permq(sort(abs(v), 1), 1 - 2*al)
    if (absp == 1) out[3] = (exact | memb ? ma/n : (1 + ma)/(n + 1))
    else out[3] = ma/n
    out[4] = -qa
    out[5] = qa
    cc = ma
  }
  else {
    out[4] = qlo
    out[5] = qhi
    if (absp == 3) {
      out[3] = plo
      cc = m0
    }
    else if (absp == 4) {
      out[3] = phi
      cc = m1
    }
    else {
      out[3] = min((1, 2*min((plo, phi))))
      cc = min((m0, m1))
      mult = 2
    }
  }
  // Monte-Carlo precision of the reported p: the standard error and the level()
  // Clopper-Pearson interval of the underlying share cc/n (Stata's permute and ritest
  // report both); under the equal-tailed 2*min rule both are doubled, the upper end
  // capped at 1
  // under an enumeration the shares are exact, not Monte-Carlo estimates, so rows 16-18
  // keep the missing they were initialised to (bug sweep BS-22, extended to the
  // incomplete-but-member case by BS-23, author ruling 2026-09-07: an enumeration that
  // lost assignments is still deterministic, so its shares carry no sampling error either.
  // The only element under permexact that still gets a binomial SE is one whose observed
  // row is genuinely absent - the same element that takes the add-one form.)
  if (!exact & !memb) {
    qq = cc/n
    out[16] = mult*sqrt(qq*(1 - qq)/n)
    al2 = al
    cpl = (cc <= 0 ? 0 : invibeta(cc, n - cc + 1, al2))
    cph = (cc >= n ? 1 : invibeta(cc + 1, n - cc, 1 - al2))
    out[17] = min((1, mult*cpl))
    out[18] = min((1, mult*cph))
  }
  // the shift interval takes the estimate's sampling error to have the NULL permutation
  // distribution, which presumes that distribution is centred at zero. When its central
  // band excludes 0 the presumption fails (a quantity that is non-zero under reassigned
  // onsets - e.g. an ame row carrying a baseline dose slope - or a probabilities design
  // far from the observed assignment) and the interval would not even contain the
  // estimate: it is left blank and flagged (row 10); the rank p-value stands
  if (qlo > 0 | qhi < 0) out[10] = 1
  else {
    out[7] = sobs - out[5]
    out[8] = sobs - out[4]
    // row 9 is the standard error IMPLIED BY INVERTING the reported shift interval,
    // (q_hi - q_lo)/(2*z(level)), on the statistic's scale. It is NOT the dispersion of the
    // null draws (that is row 6, e(perm_sd)) and not an analytic SE. It is published under
    // its own name, e(perm_se_shift), and deliberately not posted into e(V) or e(se)
    // (author rulings 2026-09-03/04, revisited under PS-2 2026-09-10; the bootstrap side
    // publishes no interval-implied SE at all). It inherits the off-centre guard: when row
    // 10 fires, rows 7, 8 and 9 stay missing together, because an SE implied by an interval
    // that was refused is meaningless - the SE column vanishing there is the design, not a
    // bug
    if (zc < . & zc > 0) out[9] = (out[5] - out[4])/(2*zc)
  }
  out[13] = 1
  return(out)
}

// per-column kernel over the accumulator. Non-studentized: s_r = b_r, s_obs = the rep-0
// comparator's b (the engine's observed-is-a-member convention; under the fast-mode
// default it equals the displayed estimate bitwise). Studentized (seaccname != ""):
// s_r = b_r/se_r, s_obs = b_cmp/se_cmp from the comparator run's own e(b)/e(V). The
// interval is always placed about the DISPLAYED estimate ob. M is the number of
// enumerated assignments (. when sampled): an element gets plain shares only when its
// n equals M
void wooldid_permstat(string scalar accname, string scalar seaccname, string scalar obsbname,
                      string scalar cmpname, string scalar cmpsename, string scalar exclname,
                      real scalar lvl, real scalar M, real scalar orow, real scalar absp,
                      string scalar outname)
{
  real matrix A, S, out
  real rowvector ob, cm, cs, ex, row
  real colvector s
  real scalar K, k, r, al, zc, sobs, ts, exk, memb
  string scalar lowav, incom, nonm
  A = st_matrix(accname)
  ob = st_matrix(obsbname)
  cm = st_matrix(cmpname)
  ex = st_matrix(exclname)
  ts = (seaccname != "")
  if (ts) {
    S = st_matrix(seaccname)
    cs = st_matrix(cmpsename)
  }
  K = cols(A)
  al = (100 - lvl)/200
  zc = invnormal(1 - al)
  out = J(18, K, .)
  out[13, .] = J(1, K, 0)
  lowav = ""
  incom = ""
  nonm = ""
  for (k = 1; k <= K; k++) {
    if (ex[1, k] != 0) continue
    if (!ts) {
      s = A[., k]
      sobs = cm[1, k]
    }
    else {
      // s is kept at FULL replication length on purpose, with holes for replications whose
      // SE is unusable, rather than compressed to the usable draws: orow is a positional
      // index into wooldid_permx's assignment columns, so s, A and the accumulator rows
      // must stay in 1:1 replication order all the way to wooldid_permmemb below.
      // wooldid_permtail does its own select() internally, so compressing here - the
      // obvious tidy-up, and exactly what permtail itself does - is not an optimisation but
      // a silent change of WHICH replication membership is tested against, and hence of the
      // share convention (BS-23). The sibling call sites in wooldid_permpw and
      // wooldid_permcellwire carry the same note
      s = J(rows(A), 1, .)
      for (r = 1; r <= rows(A); r++) {
        if (A[r, k] < . & S[r, k] < . & S[r, k] > 0) s[r] = A[r, k]/S[r, k]
      }
      if (cs[1, k] >= . | cs[1, k] <= 0 | cm[1, k] >= .) sobs = .
      else sobs = cm[1, k]/cs[1, k]
    }
    // BS-23 (author ruling 2026-09-07): exk is completeness, memb is membership - the
    // observed assignment's own row survived for this element and carries sobs. An element
    // that lost some assignments keeps the exact share over the ones that completed; the
    // add-one form is for a reference set that does not contain the observed statistic.
    // nonm names those (essentially unreachable: the identity replication is the same fit
    // as the rep-0 comparator), so the add-one caution stays attached to the case that
    // actually uses the add-one form.
    exk = (M < . & sum(s :< .) == M)
    memb = wooldid_permmemb(s, sobs, M, orow)
    if (M < . & sum(s :< .) != M) incom = incom + " " + st_matrixcolstripe(accname)[k, 2]
    if (M < . & !memb & sobs < . & sum(s :< .) >= 2) nonm = nonm + " " + st_matrixcolstripe(accname)[k, 2]
    row = wooldid_permtail(s, sobs, exk, memb, al, zc, absp)
    if (row[13] == 1 & row[10] != 1) {
      if (!ts) {
        row[7] = ob[1, k] - row[5]
        row[8] = ob[1, k] - row[4]
      }
      else {
        row[7] = ob[1, k] - row[5]*cs[1, k]
        row[8] = ob[1, k] - row[4]*cs[1, k]
        row[9] = row[9]*cs[1, k]
      }
    }
    out[., k] = row'
  }
  st_matrix(outname, out)
  st_local("wdp_incomplete", strtrim(incom))
  st_local("wdp_addone", strtrim(nonm))
}

// contrast kernel: d_r = b_r[ka] - b_r[kb] over listwise-complete rows, ranked against
// the comparator's difference; the interval about the displayed difference dobs
void wooldid_permpw(string scalar accname, string scalar cmpname, real scalar ka, real scalar kb,
                    real scalar dobs, real scalar lvl, real scalar M, real scalar orow,
                    real scalar absp, string scalar outname)
{
  real matrix A
  real rowvector cm, row
  real colvector ok, d
  real scalar al, zc, sobs, exk, memb
  A = st_matrix(accname)[., (ka, kb)]
  cm = st_matrix(cmpname)
  ok = rowmissing(A) :== 0
  d = select(A[., 1] :- A[., 2], ok)
  sobs = cm[1, ka] - cm[1, kb]
  al = (100 - lvl)/200
  zc = invnormal(1 - al)
  // BS-23 (author ruling 2026-09-07): membership, not completeness. A is the R x 2 slice
  // with every replication row kept, so the observed row indexes it directly; a row missing
  // on either column gives a missing difference and fails the test.
  exk = (M < . & rows(d) == M)
  memb = wooldid_permmemb(A[., 1] :- A[., 2], sobs, M, orow)
  row = wooldid_permtail(d, sobs, exk, memb, al, zc, absp)
  if (row[13] == 1 & row[10] != 1 & dobs < .) {
    row[7] = dobs - row[5]
    row[8] = dobs - row[4]
  }
  else {
    row[7] = .
    row[8] = .
    row[9] = .
  }
  st_matrix(outname, row')
}

// fixed-metric quadratic form. cols: the accumulator columns of the statistic's
// coefficients; C: the constraint matrix (q x m; identity for a joint zero test, the
// chain contrast for an equality test). Sigma is the covariance of the listwise-complete
// replication rows PLUS the comparator row over those columns; the metric is
// correlation-normalised with wooldid_bootjw's relative-rank guard, and a rejected metric
// posts p = . (never a residual-rank p). W_r = (C a_r)' (C Sigma C')^-1 (C a_r), W_obs the
// same on the comparator; p is the upper tail, add-one unless every M assignments
// contributed. Outputs: locals p, n, wobs; matrix outW = the W_r column (missing rows kept)
void wooldid_permqf(string scalar accname, string scalar cmpname, string scalar colslist,
                    string scalar Cname, real scalar M, real scalar orow, string scalar outp,
                    string scalar outn, string scalar outwobs, string scalar outW)
{
  real matrix A, CC, Cm, R, Si, dev
  real rowvector idx, cm, sc, ev
  real colvector ok, Wr, d
  real scalar n, r, Wobs, nex, exk, memb
  st_local(outp, ".")
  st_local(outn, "0")
  st_local(outwobs, ".")
  idx = strtoreal(tokens(colslist))
  A = st_matrix(accname)[., idx']
  cm = st_matrix(cmpname)[1, idx']
  Cm = st_matrix(Cname)
  Wr = J(rows(A), 1, .)
  st_matrix(outW, Wr)
  if (hasmissing(cm)) return
  ok = rowmissing(A) :== 0
  n = sum(ok)
  st_local(outn, strofreal(n))
  if (n < 2) return
  // the observed assignment is already one of the rows whenever its own replication
  // survived listwise on the tested columns - which is what the appending, and the share,
  // must turn on. Appending it there counts it twice: it shrinks the metric and, since
  // that row's W equals W_obs identically, guarantees nex >= 1. Fixed for the COMPLETE
  // case on 2026-09-05 and extended to the incomplete one by BS-23 (author ruling
  // 2026-09-07). exk stays the completeness flag.
  exk = (M < . & n == M)
  memb = wooldid_permmemb(A, cm, M, orow)
  CC = (exk | memb ? select(A, ok) : select(A, ok) \ cm)
  dev = CC :- mean(CC)
  R = Cm * ((dev' * dev) / (rows(CC) - 1)) * Cm'
  R = (R + R')/2
  sc = sqrt(diagonal(R))'
  if (min(sc) <= 0 | missing(sc)) return
  R = R :/ (sc' * sc)
  ev = symeigenvalues(R)
  if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) return
  Si = invsym(R) :/ (sc' * sc)
  d = Cm * cm'
  Wobs = d' * Si * d
  if (Wobs < 0 | Wobs >= .) return
  for (r = 1; r <= rows(A); r++) {
    if (!ok[r]) continue
    d = Cm * A[r, .]'
    Wr[r] = d' * Si * d
  }
  st_matrix(outW, Wr)
  nex = sum(select(Wr, Wr :< .) :>= Wobs)
  st_local(outwobs, strofreal(Wobs, "%21x"))
  // BS-23 (author ruling 2026-09-07): the share turns on membership too - the comparator
  // was not appended above whenever the observed assignment's own row survived, so adding
  // it back into the tail count here would count it twice
  if (exk | memb) st_local(outp, strofreal(nex/n, "%21x"))
  else st_local(outp, strofreal((1 + nex)/(n + 1), "%21x"))
}

// F065 (2026-09-06): max-|studentized| step-down over a set of columns (the post-treatment
// horizons of one event-study family). Fixed-metric studentization: each column's null
// draws (replications, plus the comparator row unless the enumeration is complete and
// exact) are scaled by their standard deviation; the observed statistic is the comparator.
// Rows of the output: 1 |t_obs|, 2 p_raw (the absolute-value share per column, add-one
// or exact), 3 p_singlestep (max |t_r| over ALL columns), 4 p_stepdown (max over the
// columns not yet rejected, in decreasing |t_obs| order, made monotone). Columns keep
// the accumulator's names.
// ON FAILURE (fewer than 2 complete replications, or a column with zero dispersion) the
// kernel writes J(0,0,.). Stata has no 0 x 0 matrix: MEASURED (2026-09-10, Stata 19 MP)
// st_matrix(name, J(0,0,.)) leaves NO matrix of that name - `confirm matrix' returns
// r(111) - and DROPS any matrix the name already held. So a caller must test with
// `cap confirm matrix' before any rowsof() or matlist, never with rowsof() alone; the
// drop-on-write is also what makes it safe to reuse one tempname across loop iterations.
// NOT absp-AWARE, deliberately, unlike every sibling permutation kernel: the max-|t|
// step-down (Westfall-Young / Romano-Wolf) is a statement about the two-sided |t| maximum
// over a family, so a one-sided permside() or an add-one-free permrawp share has no meaning
// for it. p_raw is reported in the step-down's own studentized metric so it is comparable
// with p_singlestep and p_stepdown beside it, NOT with the table's p. Do not thread absp
// through here: giving this kernel an absp argument "for consistency" would silently move
// every published step-down p-value.
void wooldid_permstep(string scalar accname, string scalar cmpname, string scalar colslist,
                      real scalar M, real scalar orow, string scalar outname)
{
  real matrix A, CC, T, out
  real rowvector idx, cm, sd, tobs, ord
  real colvector ok, mx
  real scalar n, H, j, h, exk, memb, nex, pprev, pj
  string matrix cs
  idx = strtoreal(tokens(colslist))
  A = st_matrix(accname)[., idx']
  cs = st_matrixcolstripe(accname)[idx', .]
  cm = st_matrix(cmpname)[1, idx']
  st_matrix(outname, J(0, 0, .))
  if (hasmissing(cm)) return
  ok = rowmissing(A) :== 0
  n = sum(ok)
  if (n < 2) return
  // BS-23 (author ruling 2026-09-07): the comparator is appended only when the observed
  // assignment's own row is NOT among the survivors. The population-SD rescale stays tied
  // to completeness; it is a common scalar across columns, so it moves the reported t_obs
  // row but cannot move a rank or a p.
  exk = (M < . & n == M)
  memb = wooldid_permmemb(A, cm, M, orow)
  CC = (exk | memb ? select(A, ok) : select(A, ok) \ cm)
  sd = sqrt(diagonal(variance(CC)))'
  if (exk) sd = sd :* sqrt((rows(CC) - 1)/rows(CC))
  if (min(sd) <= 0 | missing(sd)) return
  T = abs(select(A, ok)) :/ sd
  tobs = abs(cm) :/ sd
  H = cols(T)
  out = J(4, H, .)
  out[1, .] = tobs
  for (h = 1; h <= H; h++) {
    nex = sum(T[., h] :>= tobs[h])
    out[2, h] = (exk | memb ? nex/n : (1 + nex)/(n + 1))
  }
  mx = rowmax(T)
  for (h = 1; h <= H; h++) {
    nex = sum(mx :>= tobs[h])
    out[3, h] = (exk | memb ? nex/n : (1 + nex)/(n + 1))
  }
  // step-down: order by |t_obs| descending; at step j the max runs over columns j..H
  // of that order; monotone in j
  ord = order(-tobs', 1)'
  pprev = 0
  for (j = 1; j <= H; j++) {
    mx = rowmax(T[., ord[j..H]])
    nex = sum(mx :>= tobs[ord[j]])
    pj = (exk | memb ? nex/n : (1 + nex)/(n + 1))
    if (pj < pprev) pj = pprev
    out[4, ord[j]] = pj
    pprev = pj
  }
  st_matrix(outname, out)
  st_matrixcolstripe(outname, cs)
  st_matrixrowstripe(outname, (J(4, 1, ""), ("t_obs", "p_raw", "p_singlestep", "p_stepdown")'))
}

// F034 (2026-09-06): the bootstrap twin of wooldid_permqf. Draws A (rows = replications,
// listwise complete on the tested columns) are recentred at the OBSERVED estimates b0
// (the rep-0 comparator), the metric is the bootstrap covariance of the draws about
// their mean (the same listwise covariance e(V_boot) carries for those columns, with
// wooldid_permqf's correlation-normalised eigen guard), W_r = (Cm (b_r - b0))' S^-1
// (Cm (b_r - b0)), W_obs = (Cm b0)' S^-1 (Cm b0), p = add-one share of W_r >= W_obs
void wooldid_bootqf(string scalar accname, string scalar obsbname, string scalar colslist,
                    string scalar Cname, string scalar outp, string scalar outn,
                    string scalar outwobs, string scalar outW)
{
  real matrix A, D, Cm, R, Si, dev
  real rowvector idx, b0, sc, ev
  real colvector ok, Wr, d
  real scalar n, r, Wobs, nex
  st_local(outp, ".")
  st_local(outn, "0")
  st_local(outwobs, ".")
  idx = strtoreal(tokens(colslist))
  A = st_matrix(accname)[., idx']
  b0 = st_matrix(obsbname)[1, idx']
  Cm = st_matrix(Cname)
  Wr = J(rows(A), 1, .)
  st_matrix(outW, Wr)
  if (hasmissing(b0)) return
  ok = rowmissing(A) :== 0
  n = sum(ok)
  st_local(outn, strofreal(n))
  if (n < 2) return
  D = select(A, ok)
  dev = D :- mean(D)
  R = Cm * ((dev' * dev) / (n - 1)) * Cm'
  R = (R + R')/2
  sc = sqrt(diagonal(R))'
  if (min(sc) <= 0 | missing(sc)) return
  R = R :/ (sc' * sc)
  ev = symeigenvalues(R)
  if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) return
  Si = invsym(R) :/ (sc' * sc)
  d = Cm * b0'
  Wobs = d' * Si * d
  if (Wobs < 0 | Wobs >= .) return
  for (r = 1; r <= rows(A); r++) {
    if (!ok[r]) continue
    d = Cm * (A[r, .] - b0)'
    Wr[r] = d' * Si * d
  }
  st_matrix(outW, Wr)
  nex = sum(select(Wr, Wr :< .) :>= Wobs)
  st_local(outwobs, strofreal(Wobs, "%21x"))
  st_local(outp, strofreal((1 + nex)/(n + 1), "%21x"))
}

// adversarial line statistic on the fixed metric: wooldid_advW(b, Sigma, rel, ...) per
// replication and on the comparator, Sigma the permutation covariance of the family's
// horizon columns (replications plus comparator, listwise). rel: the relative-year list
// aligned with cols. One-sided upper tail, same accounting as wooldid_permqf
void wooldid_permadv(string scalar accname, string scalar cmpname, string scalar colslist,
                     string scalar rellist, real scalar esrel, real scalar tt, real scalar alt,
                     real scalar M, real scalar orow, string scalar outp, string scalar outn,
                     string scalar outwobs, string scalar outW)
{
  real matrix A, CC, V, dev
  real rowvector idx, cm
  real colvector ok, Wr, rel
  real scalar n, r, Wobs, nex, exk, memb, mrow, q
  st_local(outp, ".")
  st_local(outn, "0")
  st_local(outwobs, ".")
  idx = strtoreal(tokens(colslist))
  rel = strtoreal(tokens(rellist))'
  A = st_matrix(accname)[., idx']
  cm = st_matrix(cmpname)[1, idx']
  Wr = J(rows(A), 1, .)
  st_matrix(outW, Wr)
  if (hasmissing(cm) | hasmissing(rel)) return
  ok = rowmissing(A) :== 0
  n = sum(ok)
  st_local(outn, strofreal(n))
  if (n < 2) return
  // same rule as wooldid_permqf: the comparator is already one of the rows whenever the
  // observed assignment's own row survived (BS-23, author ruling 2026-09-07). mrow is that
  // row's validity on the COEFFICIENT columns; the tail below re-tests it on the W column,
  // which is the statistic actually ranked.
  exk = (M < . & n == M)
  mrow = wooldid_permmemb(A, cm, M, orow)
  CC = (exk | mrow ? select(A, ok) : select(A, ok) \ cm)
  dev = CC :- mean(CC)
  V = (dev' * dev) / (rows(CC) - 1)
  q = .
  Wobs = wooldid_advW(cm', V, rel, esrel, tt, alt, q)
  if (Wobs >= .) return
  for (r = 1; r <= rows(A); r++) {
    if (!ok[r]) continue
    q = .
    Wr[r] = wooldid_advW(A[r, .]', V, rel, esrel, tt, alt, q)
  }
  st_matrix(outW, Wr)
  ok = Wr :< .
  n = sum(ok)
  st_local(outn, strofreal(n))
  if (n < 2) return
  exk = (M < . & n == M)
  memb = 0
  if (mrow) memb = (Wr[orow] < .)
  nex = sum(select(Wr, ok) :>= Wobs)
  st_local(outwobs, strofreal(Wobs, "%21x"))
  if (exk | memb) st_local(outp, strofreal(nex/n, "%21x"))
  else st_local(outp, strofreal((1 + nex)/(n + 1), "%21x"))
}

// cell tables (the e(cfx*) matrices): per row, the cell's draws ranked against the
// comparator cell (cmp: the rep-0 harvest, 1 row, same keys), the interval about the
// displayed cell estimate; studentized under permtstat by the replication and comparator
// cell standard errors. Column secol is left as the worker posted it (analytic or blank),
// pcol/lbcol/ubcol are rewritten. Rows with no usable draws are blanked
void wooldid_permcellwire(string scalar mname, string scalar accname, string scalar cmpname,
                          string scalar seaccname, string scalar cmpsename, real scalar obscol,
                          real scalar secol, real scalar pcol, real scalar lbcol, real scalar ubcol,
                          real scalar lvl, real scalar M, real scalar orow, real scalar absp)
{
  real matrix Mx, A, S, C, CS
  string colvector mk
  string matrix cs
  real colvector s
  real rowvector row
  real scalar i, k, r, al, zc, b, sobs, ts, exk, memb, so
  transmorphic t
  Mx = st_matrix(mname)
  A = st_matrix(accname)
  C = st_matrix(cmpname)
  cs = st_matrixcolstripe(accname)
  ts = (seaccname != "")
  if (ts) {
    S = st_matrix(seaccname)
    CS = st_matrix(cmpsename)
  }
  t = asarray_create()
  for (k = 1; k <= rows(cs); k++) asarray(t, cs[k, 2], k)
  mk = wooldid_bootcellkey(mname, 0)
  al = (100 - lvl)/200
  zc = invnormal(1 - al)
  for (i = 1; i <= rows(Mx); i++) {
    Mx[i, pcol] = .
    Mx[i, lbcol] = .
    Mx[i, ubcol] = .
    if (!asarray_contains(t, mk[i])) continue
    k = asarray(t, mk[i])
    b = Mx[i, obscol]
    if (b >= .) continue
    if (!ts) {
      s = A[., k]
      sobs = C[1, k]
      so = 1
    }
    else {
      s = J(rows(A), 1, .)
      for (r = 1; r <= rows(A); r++) {
        if (A[r, k] < . & S[r, k] < . & S[r, k] > 0) s[r] = A[r, k]/S[r, k]
      }
      so = CS[1, k]
      if (so >= . | so <= 0 | C[1, k] >= .) continue
      sobs = C[1, k]/so
    }
    // BS-23 (author ruling 2026-09-07): membership, not completeness (the cell accumulator
    // carries one row per replication, so orow indexes it directly)
    exk = (M < . & sum(s :< .) == M)
    memb = wooldid_permmemb(s, sobs, M, orow)
    row = wooldid_permtail(s, sobs, exk, memb, al, zc, absp)
    if (row[13] != 1) continue
    Mx[i, pcol] = row[3]
    if (row[10] == 1) continue
    Mx[i, lbcol] = b - row[5]*so
    Mx[i, ubcol] = b - row[4]*so
  }
  st_matrix(mname, Mx)
}

// permsave(): after wooldid_bootsavemk has laid out _rep and the coefficient (and se_*)
// columns, append the assignment block - one column per cluster, holding the onset that
// cluster carried in each replication (row 1: the observed onset)
void wooldid_permsaveasg(string scalar assignname, string scalar clusvals, string scalar obsons)
{
  real matrix Asg
  real rowvector cv, oo
  real scalar c, idx, nd
  string scalar nm
  Asg = st_matrix(assignname)
  cv = strtoreal(tokens(clusvals))
  oo = strtoreal(tokens(obsons))
  // D-10 (author ruling 2026-09-11): both drops now report, because the output is a FILE
  // the user analyses elsewhere and a silent gap in it is the worst case. (a) the return
  // below writes NO assignment block at all, leaving a dataset indistinguishable from a run
  // that never asked for one; (b) the `_st_varindex' skip further down drops one cluster's
  // column whenever its strtoname()d name is already taken - by a coefficient column laid
  // out by wooldid_bootsavemk, or by an earlier cluster whose id maps to the same string -
  // and because the columns are named rather than positional, nothing in the file marks the
  // gap. wooldid_permsave turns these two locals into cautions.
  st_local("wdpsasgerr", "")
  st_local("wdpsasgdrop", "0")
  if (cols(cv) != rows(Asg) | cols(oo) != rows(Asg)) {
    st_local("wdpsasgerr", "1")
    return
  }
  nd = 0
  for (c = 1; c <= rows(Asg); c++) {
    nm = strtoname("onset_" + strofreal(cv[c], "%20.0g"))
    if (_st_varindex(nm) < .) {
      nd++
      continue
    }
    idx = st_addvar("double", nm)
    st_store(., idx, oo[c] \ Asg[c, .]')
  }
  st_local("wdpsasgdrop", strofreal(nd))
}


// the comparator run's own standard errors by column (the e(V) diagonal; sentinel and
// non-positive entries stay missing, as in wooldid_bootgrab)
void wooldid_permcse(string scalar outname, real scalar K)
{
  real matrix V
  real rowvector out
  real scalar k, v
  out = J(1, K, .)
  V = st_matrix("e(V)")
  // POSITIONAL is admissible here and would not be in wooldid_bootgrab, whose header is
  // emphatic that the same harvest must be name-aligned: this reads the REP-0 COMPARATOR's
  // e(V), and wooldid_permrun has just proven the comparator's colfullnames equal the
  // observed run's (error 459 otherwise), so no resample-shortened table can reach this
  // function. The guard below is therefore not a fallback: rows(V) != K leaves cmpse
  // entirely missing, wooldid_permstat then sets sobs = . for every column, and the whole
  // studentized permutation surface returns live = 0 and blanks - silently, rc 0, no
  // caution. If it is ever reachable, that is a bug rather than a degradation
  if (rows(V) == K) {
    for (k = 1; k <= K; k++) {
      v = V[k, k]
      if (v > 0 & v < .) out[k] = sqrt(v)
    }
  }
  st_matrix(outname, out)
}

// the adversarial DISPLAY LINE, V-free, exactly as wooldid_advstat builds it: the real
// horizon coefficients (missing -> 0) against their relative years, then the planted
// zero-valued reference row(s) - esf: relyr = esrel (preceded by esrel - 1 under
// timetrends); pooled: one horizon below the earliest real horizon - fitted by OLS
// (altmethod: on the rows with relyr <= -1 only); published as [slope, intercept_rt0]
void wooldid_permadvline(string scalar cmpname, string scalar colslist, string scalar rellist,
                         real scalar esrel, real scalar tt, real scalar esf, real scalar alt,
                         string scalar outname)
{
  real rowvector idx, cm
  real colvector rel, b, r, sel
  real matrix X, F
  real scalar k
  idx = strtoreal(tokens(colslist))
  rel = strtoreal(tokens(rellist))'
  cm = st_matrix(cmpname)[1, idx']
  b = cm'
  // CONSEQUENCE of the missing -> 0 substitution, which the note above states but does not
  // price: a horizon with no estimable coefficient enters the OLS fit below as a real
  // observation at exactly zero, so a family that lost horizons reports a slope pulled
  // toward the reference line in proportion to how many dropped out, and nothing in
  // e(adv_line) marks it. Dropping those rows instead would be equally defensible; the
  // zero is chosen for bitwise agreement with the display frame wooldid_advstat plants
  for (k = 1; k <= rows(b); k++) if (b[k] >= .) b[k] = 0
  r = rel
  if (esf) {
    if (tt) {
      b = b \ 0
      r = r \ (esrel - 1)
    }
    b = b \ 0
    r = r \ esrel
  }
  else {
    b = b \ 0
    r = r \ (min(rel) - 1)
  }
  if (alt) {
    sel = selectindex(r :<= -1)
    b = b[sel]
    r = r[sel]
  }
  st_matrix(outname, J(1, 2, .))
  if (rows(b) < 2) return
  X = J(rows(r), 1, 1), r
  if (rank(X'X) < 2) return
  F = invsym(X'X) * X'b
  st_matrix(outname, (F[2], F[1]))
}

// PS-45 (author ruling 2026-09-07): the relative period a horizon label names. The F2
// labels are pre#/post#; e(V)'s own row for relative period 0 reads "contemp" under
// oldsyntax. Anything else - including the "it_post0" tail a shorter family prefix picks
// up from a longer family's row name - returns missing, which is what makes the
// prefix-plus-parse column scan below self-disambiguating
real scalar wooldid_esrelx(string scalar s)
{
  if (substr(s, 1, 4) == "post") return(strtoreal(substr(s, 5, .)))
  if (substr(s, 1, 3) == "pre") return(-strtoreal(substr(s, 4, .)))
  if (s == "contemp") return(0)
  return(.)
}

// PS-45: the SIMULTANEOUS (sup-t / uniform) band for one event-study family, simulated
// from the family's own analytic covariance. The pointwise intervals in the table cover
// each horizon with probability level(); a simultaneous band covers EVERY horizon of the
// family at once with that probability, which is what a reader who scans a whole event
// study for "the effect is zero everywhere before onset" actually needs.
// Construction: studentize the family's covariance to a correlation matrix R, draw
// z ~ N(0, R), and take c = the level() quantile of max_h |z_h| (Montiel Olea and
// Plagborg-Moller 2019). The band is b_h +/- c*se_h, so it is the pointwise interval
// widened by one common factor c/z_{level} >= 1. R is eigen-decomposed with the
// eigenvalues floored at zero rather than factored by cholesky, because a rank-deficient
// family - a structurally collinear horizon set - is ordinary here, and cholesky would
// simply fail on it.
// Horizons with a non-positive or missing variance are dropped from the family before
// the draw: they have no pointwise interval either.
// WHICH COVARIANCE, and why no check: the half-width is c*sqrt(e(V_es_*)[j,j]), while the
// pointwise interval this band is drawn over comes from e(V)/r(table). Those are the same
// quantity computed twice and agree to about 1e-16 only on the ordinary analytic path -
// wooldid_vspl records the measurement and gates every splice on a 1e-8 relative diagonal
// check for exactly that reason. This kernel applies no such gate, and needs none only
// because every divergent path also WITHHOLDS e(V_es_*): the F2 capture that creates it is
// gated on impvalues being empty and on marginsdof not being nose, and e(V_es_ame*) is
// withheld under AMESE containment. PS-45 shows the cost if that invariant ever breaks - on
// the bootstrap side the band came out NARROWER than the pointwise interval beside it.
// RNG: this kernel draws with rnormal() on the CALLER's stream. Reproducibility is entirely
// wooldid_esband's doing (it saves c(rngstate), sets the fixed internal seed, and restores),
// and because that seed is set once for the whole family loop each family consumes
// successive draws - so a family's critical value depends on how many families preceded it,
// and adding subgroups() or an estimand moves an unrelated family's band in its low-order
// digits (simulation noise at nd = 10000). The nd >= 100 floor below is a defensive guard on
// a parameter no user sets.
real scalar wooldid_esbanda(string scalar vn, string scalar bn, real scalar lev,
                            real scalar nd, string scalar outn, string scalar critl)
{
  real matrix V0, V, U, L, Z, out
  real rowvector b0, lam
  string matrix cs, ocs
  real colvector sel, sd, M
  real scalar k, m, j, r, c
  V0 = st_matrix(vn)
  b0 = st_matrix(bn)
  cs = st_matrixcolstripe(vn)
  k = rows(V0)
  if (k < 2 | cols(V0) != k | cols(b0) != k | nd < 100) return(0)
  sel = J(0, 1, .)
  for (j = 1; j <= k; j++) {
    if (V0[j, j] > 0 & V0[j, j] < . & b0[j] < . & wooldid_esrelx(cs[j, 2]) < .) sel = sel \ j
  }
  m = rows(sel)
  if (m < 2) return(0)
  V = J(m, m, 0)
  for (j = 1; j <= m; j++) for (r = 1; r <= m; r++) V[j, r] = V0[sel[j], sel[r]]
  sd = J(m, 1, .)
  for (j = 1; j <= m; j++) sd[j] = sqrt(V[j, j])
  for (j = 1; j <= m; j++) for (r = 1; r <= m; r++) V[j, r] = V[j, r]/(sd[j]*sd[r])
  V = (V + V')/2
  U = .
  lam = .
  symeigensystem(V, U, lam)
  for (j = 1; j <= m; j++) if (lam[j] < 0 | lam[j] >= .) lam[j] = 0
  L = U * diag(sqrt(lam))
  Z = L * rnormal(m, nd, 0, 1)
  M = J(nd, 1, .)
  for (r = 1; r <= nd; r++) M[r] = max(abs(Z[., r]))
  M = sort(M, 1)
  c = wooldid_bootqs(M, lev/100)
  if (c >= . | c <= 0) return(0)
  out = J(3, m, .)
  ocs = J(m, 2, "")
  for (j = 1; j <= m; j++) {
    out[1, j] = b0[sel[j]] - c*sd[j]
    out[2, j] = b0[sel[j]] + c*sd[j]
    out[3, j] = wooldid_esrelx(cs[sel[j], 2])
    ocs[j, 2] = cs[sel[j], 2]
  }
  st_matrix(outn, out)
  st_matrixrowstripe(outn, (J(3, 1, ""), ("lb", "ub", "relyr")'))
  st_matrixcolstripe(outn, ocs)
  st_local(critl, strofreal(c, "%21x"))
  return(m)
}

// PS-45: the same band from the BOOTSTRAP replications - one pass over the accumulator
// that already feeds the joint test, and no distributional assumption. For each
// listwise-complete replication take M_r = max_h |b_rh - b_obs_h| / sd_h, with sd_h the
// standard deviation of that horizon's draws over the same rows (the fixed metric the
// permutation step-down also uses); c is the level() quantile of the M_r under
// wooldid_bootci's order-statistic rule, and the band is b_obs_h +/- c*sd_h.
// LISTWISE is load-bearing: the maximum is a property of a whole replication, so a row
// missing one horizon cannot contribute a maximum over all of them. The sd is computed
// on the same rows, so the band's own scale and its critical value come from one sample -
// which is why sd_h can differ very slightly from the table's per-column standard error
// when some replications failed.
// The family's columns are found by scanning the observed coefficient names for the
// family prefix and keeping those whose remainder parses as a relative period.
// DEGENERATE-HORIZON POLICY, and it is the OPPOSITE of its matched pair wooldid_esbanda's.
// esbanda DROPS a horizon with non-positive/missing variance and bands the rest; this
// kernel REFUSES THE WHOLE FAMILY for one bad horizon (each `return(0)' below: a missing
// observed coefficient, a missing or non-positive observed analytic SE in the studentized
// branch, a zero or missing draw SD). So a run that yields an 8-horizon analytic band can
// yield no bootstrap band at all when horizon 9 is structurally unidentified, and the only
// user-visible trace is the caller's generic "no simultaneous band could be formed"
// caution. The case for it: the critical value is a maximum over the family's horizons
// taken replication by replication, so a horizon that cannot supply a denominator at every
// retained replication changes what the maximum is over. Dropping it would be defensible
// too - this asymmetry is recorded, not ruled.
// LOCALS THAT CHANGE MEANING in the body below, in a 155-line function with two
// near-duplicate construction branches: `mx' is a row-completeness FLAG in the selection
// loop and a running max|t| in the band loops; `d' is a parsed relative year in the column
// scan, then a standardized deviation in the band loops, then a running sum of squares in
// the SD loop. The two branches differ only in the denominator (the replication's own SE
// versus the draw SD) and in which SE sets the endpoints.
real scalar wooldid_bootband(string scalar accn, string scalar obsn, string scalar px,
                             real scalar lev, string scalar outn, string scalar critl,
                             string scalar nl, string scalar sen, string scalar aon,
                             string scalar stl, string scalar whyl)
{
  real matrix A, O, out, S, AO
  string matrix cs, ocs
  real colvector sel, rel, idx, sd, M, p
  real scalar K, R, m, j, r, n, c, d, mx, mu, lp, stud
  A = st_matrix(accn)
  O = st_matrix(obsn)
  cs = st_matrixcolstripe(obsn)
  K = cols(O)
  R = rows(A)
  if (cols(A) != K | R < 2) return(0)
  // PS-45 (2026-09-10): the studentized construction, taken whenever boottstat supplied the
  // replications' own standard errors AND the observed analytic ones. Anything short of both,
  // conformable, falls back to the coefficient-draw SD - which is the only form available
  // without boottstat and stays the default there.
  stud = 0
  st_local(stl, "0")
  // D-8 (author ruling 2026-09-11): every refusal below says WHY, naming the horizon that
  // caused it, so the caller's caution can name it too. An EMPTY reason means "this family
  // is not present here" (fewer than two horizons carry the prefix, or there are no
  // replications) rather than "this family was refused" - which is what lets the caller tell
  // a genuinely absent family apart from one it declined to band.
  st_local(whyl, "")
  if (sen != "" & aon != "") {
    S = st_matrix(sen)
    AO = st_matrix(aon)
    if (rows(S) == R & cols(S) == K & cols(AO) == K & rows(AO) >= 1) stud = 1
  }
  lp = strlen(px)
  sel = J(0, 1, .)
  rel = J(0, 1, .)
  for (j = 1; j <= K; j++) {
    if (substr(cs[j, 2], 1, lp) != px) continue
    d = wooldid_esrelx(substr(cs[j, 2], lp + 1, .))
    if (d >= .) continue
    sel = sel \ j
    rel = rel \ d
  }
  m = rows(sel)
  if (m < 2) return(0)
  p = order(rel, 1)
  sel = sel[p]
  rel = rel[p]
  idx = J(0, 1, .)
  for (r = 1; r <= R; r++) {
    mx = 1
    for (j = 1; j <= m; j++) {
      if (A[r, sel[j]] >= .) mx = 0
      // a studentized replication also needs a usable denominator at every horizon
      if (stud) {
        if (S[r, sel[j]] >= . | S[r, sel[j]] <= 0) mx = 0
      }
    }
    if (mx) idx = idx \ r
  }
  n = rows(idx)
  if (n < 2) {
    st_local(whyl, "fewer than two replications are complete over every horizon of the family")
    return(0)
  }
  sd = J(m, 1, .)
  if (stud) {
    // the endpoints are the OBSERVED analytic standard errors, as the printed percentile-t
    // interval's are; the draws only supply the critical value
    for (j = 1; j <= m; j++) {
      if (O[1, sel[j]] >= .) {
        st_local(whyl, "the observed coefficient at horizon " +
                       substr(cs[sel[j], 2], lp + 1, .) + " is missing")
        return(0)
      }
      if (AO[1, sel[j]] >= . | AO[1, sel[j]] <= 0) {
        st_local(whyl, "the observed analytic standard error at horizon " +
                       substr(cs[sel[j], 2], lp + 1, .) +
                       " is missing or non-positive, so it cannot studentize")
        return(0)
      }
      sd[j] = AO[1, sel[j]]
    }
    M = J(n, 1, .)
    for (r = 1; r <= n; r++) {
      mx = 0
      for (j = 1; j <= m; j++) {
        d = abs(A[idx[r], sel[j]] - O[1, sel[j]])/S[idx[r], sel[j]]
        if (d > mx) mx = d
      }
      M[r] = mx
    }
    M = sort(M, 1)
    c = wooldid_bootqs(M, lev/100)
    if (c >= . | c <= 0) {
      st_local(whyl, "the simultaneous critical value could not be formed from the draws")
      return(0)
    }
    out = J(3, m, .)
    ocs = J(m, 2, "")
    for (j = 1; j <= m; j++) {
      out[1, j] = O[1, sel[j]] - c*sd[j]
      out[2, j] = O[1, sel[j]] + c*sd[j]
      out[3, j] = rel[j]
      ocs[j, 2] = substr(cs[sel[j], 2], lp + 1, .)
    }
    st_matrix(outn, out)
    st_matrixrowstripe(outn, (J(3, 1, ""), ("lb", "ub", "relyr")'))
    st_matrixcolstripe(outn, ocs)
    st_local(critl, strofreal(c, "%21x"))
    st_local(nl, strofreal(n))
    st_local(stl, "1")
    return(m)
  }
  for (j = 1; j <= m; j++) {
    if (O[1, sel[j]] >= .) {
      st_local(whyl, "the observed coefficient at horizon " +
                     substr(cs[sel[j], 2], lp + 1, .) + " is missing")
      return(0)
    }
    mu = 0
    for (r = 1; r <= n; r++) mu = mu + A[idx[r], sel[j]]
    mu = mu/n
    d = 0
    for (r = 1; r <= n; r++) d = d + (A[idx[r], sel[j]] - mu)^2
    sd[j] = sqrt(d/(n - 1))
    if (sd[j] <= 0 | sd[j] >= .) {
      st_local(whyl, "the bootstrap standard deviation at horizon " +
                     substr(cs[sel[j], 2], lp + 1, .) + " is zero or missing")
      return(0)
    }
  }
  M = J(n, 1, .)
  for (r = 1; r <= n; r++) {
    mx = 0
    for (j = 1; j <= m; j++) {
      d = abs(A[idx[r], sel[j]] - O[1, sel[j]])/sd[j]
      if (d > mx) mx = d
    }
    M[r] = mx
  }
  M = sort(M, 1)
  c = wooldid_bootqs(M, lev/100)
  if (c >= . | c <= 0) {
    st_local(whyl, "the simultaneous critical value could not be formed from the draws")
    return(0)
  }
  out = J(3, m, .)
  ocs = J(m, 2, "")
  for (j = 1; j <= m; j++) {
    out[1, j] = O[1, sel[j]] - c*sd[j]
    out[2, j] = O[1, sel[j]] + c*sd[j]
    out[3, j] = rel[j]
    ocs[j, 2] = substr(cs[sel[j], 2], lp + 1, .)
  }
  st_matrix(outn, out)
  st_matrixrowstripe(outn, (J(3, 1, ""), ("lb", "ub", "relyr")'))
  st_matrixcolstripe(outn, ocs)
  st_local(critl, strofreal(c, "%21x"))
  st_local(nl, strofreal(n))
  return(m)
}

// PS-62 (six-lens feasibility review, author rulings 2026-09-07): splice ONE event-study
// family's covariance block into the posted e(V). vn is the working copy of e(V) (a
// Stata tempname, replaced in place so its stripes survive), en the family's e(V_es_*),
// px the e(V) row-name prefix for the family ("ES_att_", "ES_ame_it_SG2_", ...).
// Returns the number of off-diagonal PAIRS written, or 0 when the family was skipped.
//
// The whole safety argument is the per-family gate below. Every horizon that carries a
// positive-finite variance in e(V) must (a) match a row by name and (b) agree with the
// V_es diagonal to 1e-8 relative. The two are the same quantity computed twice - e(V)
// squares the standard error column of margins' r(table), V_es is r(V) itself - so on
// the ordinary analytic path they agree to about 1e-16 (measured), and ANY path where
// the reported scale is not the margins scale, or the variance is not the margins
// variance, fails the check and the family is left alone. That single test is what keeps
// impvalues (per-cluster t-test standard errors), the AMESE-scrubbed rows, a rescaled
// reporting scale and any future divergence out, without naming any of them.
// Rows whose e(V) diagonal is NOT positive-finite - a structural zero, or the -1e100
// sentinel the partial-identification rescue writes - keep a zero row and column: their
// diagonal is left exactly as posted, so a standard error that displays as missing today
// still does. (The ruling's wording is "zero its row and column"; the reason it gives is
// that `ereturn post' can zero a nonsymmetric or highly singular VCE wholesale. That
// hazard does not arise here - the splice reaches e(V) through `ereturn repost', which
// was measured to accept an indefinite matrix and a -1e100 diagonal verbatim - so only
// the off-diagonals are zeroed and the sentinel stands. See review_process/
// held_for_ruling_2026-09-09.md.)
real scalar wooldid_vspl(string scalar vn, string scalar en, string scalar px)
{
  real matrix V, E
  string matrix vs, es
  real colvector idx
  real scalar n, m, i, j, a, b, nsp
  string scalar want, want2
  V = st_matrix(vn)
  E = st_matrix(en)
  if (rows(V) < 2 | rows(E) < 2) return(0)
  vs = st_matrixcolstripe(vn)
  es = st_matrixcolstripe(en)
  n = rows(vs)
  m = rows(es)
  if (rows(E) != m | cols(E) != m) return(0)
  idx = J(m, 1, 0)
  for (j = 1; j <= m; j++) {
    want = px + es[j, 2]
    // under oldsyntax the e(V) row for relative period 0 is named "contemp" while the
    // F2 label is always "post0"; without this alias the family half-splices in silence
    want2 = (es[j, 2] == "post0" ? px + "contemp" : "")
    for (i = 1; i <= n; i++) {
      if (vs[i, 2] == want) idx[j] = i
      else if (want2 != "" & vs[i, 2] == want2) idx[j] = i
    }
  }
  for (j = 1; j <= m; j++) {
    if (idx[j] == 0) {
      // a horizon that HAS a variance but no row is a naming mismatch, not a dropped
      // row: refuse the family rather than splice part of it
      if (E[j, j] > 0 & E[j, j] < .) return(0)
      continue
    }
    a = idx[j]
    if (V[a, a] > 0 & V[a, a] < .) {
      if (E[j, j] <= 0 | E[j, j] >= .) return(0)
      if (abs(V[a, a] - E[j, j]) > 1e-8*abs(E[j, j])) return(0)
    }
  }
  nsp = 0
  for (j = 1; j <= m; j++) {
    if (idx[j] == 0) continue
    a = idx[j]
    if (V[a, a] <= 0 | V[a, a] >= .) continue
    for (i = j + 1; i <= m; i++) {
      if (idx[i] == 0) continue
      b = idx[i]
      if (V[b, b] <= 0 | V[b, b] >= .) continue
      if (E[j, i] >= .) continue
      V[a, b] = E[j, i]
      V[b, a] = E[j, i]
      nsp++
    }
  }
  if (nsp == 0) return(0)
  st_replacematrix(vn, V)
  return(nsp)
}

end
