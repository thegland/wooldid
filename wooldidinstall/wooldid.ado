
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
    * rewire the coefficient tables and post the stripe-aligned surfaces (eclass sub)
    wooldid_permwire , acc(`wdpacc') obsb(`wdpobsb') cmpb(`wdpcmpb') seacc(`wdpsacc') cmpse(`wdpcmpse') level(`level') m(`wdp_M') done(`wdp_r_done') fail(`wdp_r_fail') blank(`wdp_r_blank') ts(`wdp_ts') method(`permute') reps(`wdp_r_reps') ase(`permanalyticse') clus(`cluster') gt(`gteffects') absp(`permabsp') side(`permside') rawp(`permrawp')
    * the legacy per-estimand scalars, read from the posted kernel by name so they are
    * the SAME numbers the tables carry (one convention: equal-tailed p, 5th-95th null
    * band, permutation SD; on the t-ratio scale under permtstat)
    tempname wdpk
    matrix `wdpk' = e(perm_kernel)
    foreach wdpe in `wdp_est' {
      local wdpc = colnumb(`wdpk', "main_`wdpe'")
      if `wdpc' < . {
        ereturn scalar perm_p_`wdpe' = `wdpk'[3, `wdpc']
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
    if `wdp_r_done' <= 10000 cap ereturn matrix permdist = `wdpdist'
    * NFPERM-12 (lean kernel): the estimand comparators actually used for the legacy
    * scalars, on the e(perm_stat) scale (the full vector is e(perm_comparator_b))
    cap ereturn matrix perm_comparator = `wdpcomp'
    * joint tests, subgroup equality, adversarial and contrasts (eclass sub)
    wooldid_permjnt , acc(`wdpacc') obsb(`wdpobsb') cmpb(`wdpcmpb') level(`level') m(`wdp_M') ts(`wdp_ts') jt(`jointtests') adv(`adversarial') advalt(`altmethod_adversarial') advesrel(`esrelativeto') advtt(`timetrends') esf(`esfixedbaseperiod') absp(`permabsp') side(`permside') rawp(`permrawp')
    * the continuous-treatment plot cells (eclass sub)
    if "`makecfxplots'" != "" & `wdp_ncfx' > 0 {
      local wdp_cfxnames
      forval wdk = 1/`wdp_ncfx' {
        local wdp_cfxnames `wdp_cfxnames' `wdb_cfxname`wdk''
      }
      wooldid_permcellsx , cell(`wdpcell') cell0(`wdpcell0') ncfx(`wdp_ncfx') names(`wdp_cfxnames') level(`level') m(`wdp_M') ts(`wdp_ts') lattype(`latticetype') cfxlat(`cfxlattice') absp(`permabsp') side(`permside') rawp(`permrawp')
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
    wooldid_bootwire , mode(`bootstrap') reps(`reps') level(`level') bootalpha(`bootalpha') clus(`cluster') acc(`wdbacc') obsb(`wdbobsb') comp(`wdbcomp') `ties' done(`wdb_r_done') fail(`wdb_r_fail') blank(`wdb_r_blank') rseed(`rseed') bootdf(`bootdf') bdisp(`bootdisplay') nclus(`wdb_r_nclus') seacc(`=cond("`boottstat'" != "", "`wdbsacc'", "")') aseout(`=cond("`boottstat'" != "", "`wdbaseobs'", "")') jkacc(`=cond("`bootdisplay'" == "bca", "`wdbjk'", "")')
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
    if "`bootsave'" != "" wooldid_bootsave , file(`bootsave') acc(`wdbacc') obsb(`wdbobsb') seacc(`=cond("`boottstat'" != "", "`wdbsacc'", "")') aseobs(`=cond("`boottstat'" != "", "`wdbaseobs'", "")') pwseacc(`=cond("`boottstat'" != "", "`wdbpwsacc'", "")') wmat(`=cond("`cluster'" != "", "`wdbW'", "")') wobs(`=cond("`bootstrap'" == "bayes", 1/max(`wdb_r_nclus', 1), 1)') jwacc(`=cond("`boottstat'" != "", "e(boot_jointW)", "")') jobs(`=cond("`boottstat'" != "", "`wdbjobs'", "")') awacc(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "e(boot_advW)", "")') aobs(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbaobs'", "")') seed(`rseed') rng(`wdb_r_rng0') `replace'
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
syntax [varlist(numeric default=none)] [if/]  [aw fw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist numeric) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose fastcollapse fullprecisionslow GTEffects PERMute(string) PERMREPS(integer 0) PERMEXact PERMTStat PERMTci PERMSTRata(varname numeric) PERMPRobs(string asis) PERMSAve(string) PERMANalyticse PERMABsp PERMSide(string) PERMRAWp BOOTstrap(string) reps(integer 0) ties Level(real 0) dots(integer -1) rseed(string) BOOTALPha(string) BOOTDf(string) BOOTDISPlay(string) BOOTTStat BOOTSAve(string)]

* Execute update before running main program
  if "`update'" == "update" {
    net install wooldid, from(https://raw.githubusercontent.com/thegland/wooldid/master/wooldidinstall/) replace
    dis as text "Update complete; exiting wooldid. Other dependencies (reghdfe, ftools, gtools, ppmlhdfe (optional)) must be updated manually if necessary."
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
      if "`safety'" == "off" | `esprelength' <= 1{
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
          dis as text "(Or request estimates for fewer pre-treatment horizons.)"
          dis as text " "
          * If this condition triggers, event study coefficients are very difficult to interpret
        }
      }
    }
    else{
      * if using timetrends, never use the 2nd observation either - dont just always skip the first!
      if "`safety'" == "off" | `esprelength' <= 1{
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
          dis as text "(Or: request estimates for fewer pre-treatment horizons, drop request to include timetrends.)"
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
  * the + 1 is to account for the 0 reference period
  if (`espre' + `espl'+ 1 < 3 & "`timetrends'" != "") | `espl' == 0 | `espre' == 0 local adversarial

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

  if "`subgroups'" != "" {
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
    * initialize these for later
    foreach sgcase in `sglist' {
      local joint_horiz_i_s`sgcase'
      local joint_Choriz_i_s`sgcase'
    }
  }






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

    * initialize
		local contreatcontrolinserterFE
		local contreatcontrolinserter

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
      local contreatcontrolinserterFE  i.`null'##`ringerpoly'
			if "`latticepart2hash'" != ""  local contreatcontrolinserterFE `contreatcontrolinserterFE'  i.`latticevar'##`ringerpoly'
			*using latticepart2hash as a check here is arbitrary and a different lattice byproduct could be used
    }

    * built on top of the initialized contreatcontrolisnerterFE -- add lattice components only if useful
    if strpos("`contreatcontrols'",  "byTreated") != 0  local contreatcontrolinserterFE `contreatcontrolinserterFE' i.`evertreated'`ringerinj'
    if strpos("`contreatcontrols'",  "byTreated") != 0 & "`latticepart2hash'" != "" local contreatcontrolinserterFE `contreatcontrolinserterFE'  i.`evertreated'`latticepart2hash'
    if strpos("`contreatcontrols'",  "byTtre") != 0  local contreatcontrolinserterFE `contreatcontrolinserterFE' i.`modified_ttre'`ringerinj'
    if strpos("`contreatcontrols'",  "byTtre") != 0  & "`latticepart2hash'" != "" local contreatcontrolinserterFE `contreatcontrolinserterFE'  i.`modified_ttre'`latticepart2hash'
    if strpos("`contreatcontrols'",  "byCohort") != 0  local contreatcontrolinserterFE `contreatcontrolinserterFE'  i.`i'`ringerinj'
    if strpos("`contreatcontrols'",  "byCohort") != 0  & "`latticepart2hash'" != "" local contreatcontrolinserterFE `contreatcontrolinserterFE'   i.`i'`latticepart2hash'

    * handle by period - this is in addition to, not in lieu of the others
    if strpos("`contreatcontrols'",  "byTime") != 0  local contreatcontrolinserterFE `contreatcontrolinserterFE' i.`t'`ringerinj'
    if strpos("`contreatcontrols'",  "byTime") != 0  & "`latticepart2hash'" != "" local contreatcontrolinserterFE `contreatcontrolinserterFE'     i.`t'`latticepart2hash'


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
  if "`safety'" != "off" {
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
  if "`jointtests'" != "" & "`bootstrap'`permute'" == "" {
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
      if "`safetyinput'" != "off" {
        dis as error "Error: jointtests is not compatible with wooldid's analytic clustered SEs when clustering at the level of the cohort variable i or above (`i' is nested within the clustering variable `jtcoarse'). The joint test p-values are severely oversized in this configuration: the cluster-robust covariance cannot capture the within-cluster cross-horizon variation that the joint tests depend on."
        dis as error "Consider wooldid's built-in randomization inference (the permute() option) or bootstrap-based clustered SEs for joint-type inference here, or cluster at a finer level. Alternatively, specify safety(off) to proceed anyway; this error then becomes a caution and the joint test p-values are reported as-is (they remain unreliable)."
        error 198
      }
      else {
        dis as text "Caution: jointtests was combined with analytic clustered SEs clustered at the level of the cohort variable i or above (`i' is nested within the clustering variable `jtcoarse'), and safety(off) was specified, so estimation is proceeding - but the joint test p-values are severely oversized in this configuration: the cluster-robust covariance cannot capture the within-cluster cross-horizon variation that the joint tests depend on. Consider wooldid's built-in randomization inference (the permute() option) or bootstrap-based clustered SEs for joint-type inference here, or cluster at a finer level."
      }
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
      if "`verbose'" != "" dis "ppmlhdfe `y' `contreatcontrolinserter' `trendinjector'    i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'   `contreatcontrolinserterFE'     `ccc_absorbinserter'  ) `vcespecifier' tol(1e-`regtol') d "
      qui cap ppmlhdfe `y' `contreatcontrolinserter' `trendinjector'    i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'    `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'    `contreatcontrolinserterFE'  `ccc_absorbinserter' ) `vcespecifier' tol(1e-`regtol') d
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
      if "`verbose'" != "" dis "reghdfe `y' `contreatcontrolinserter' `trendinjector' i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'   `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'     `contreatcontrolinserterFE'   `ccc_absorbinserter'  ) `vcespecifier' tol(1e-`regtol')  "
      qui reghdfe `y' `contreatcontrolinserter' `trendinjector'  i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'   `interinject'     `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' , absorb(i.`i' i.`t' `fe'     `contreatcontrolinserterFE'   `ccc_absorbinserter'  ) `vcespecifier' tol(1e-`regtol')
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
	qui count if `dropcount' & e(sample) == 0
	local droppedobs = r(N)

  * grab some items from the model - sample indicator variable, residual degrees of freedom, etc
  tempvar esamp
  qui gen byte `esamp' = e(sample)
  local df_r = e(df_r)

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
  * quiet estimation call swallows the only warning). Detection + caution are
  * unconditional (outside fast mode); the att-rescue REFIT through the
  * explicit-dummies regress path is deliberately NARROW: it runs only under a firing
  * containment predicate (the blackout-collateral cases this remedy exists for -
  * binary timetrends' designed b-only fallback is NOT overridden), only for a
  * single-variable cluster() (official regress has no multi-way cluster), never for
  * poisson, and behind a maxvar feasibility count.
  if "`fastnose'" == "" {
    wooldid_amescreen , mode(vzero) r1fired(`amescrub')
    if r(vzero) == 1 & `amescrub' == 1 & "`poisson'" != "poisson" & wordcount("`cluster'") == 1 & "`customvce'" == "" & "`fe'" == "" & "`ccc_absorbinserter'" == "" {
      * cap: on huge-cardinality i the levels list can overflow the macro length
      * (rc 920) - treat as refit-infeasible rather than aborting the run
      local amevzi 99999
      local amevzt 99999
      cap qui glevelsof `i' if `esamp'
      if _rc == 0 local amevzi : word count `r(levels)'
      cap qui glevelsof `t' if `esamp'
      if _rc == 0 local amevzt : word count `r(levels)'
      if `amevzi' + `amevzt' <= 450 {
        dis as text "Caution: the estimation command returned an all-zero variance-covariance matrix; wooldid is re-estimating through the explicit-dummies fallback to recover the standard errors that survive (recovered cluster-robust variances for treated-cell quantities remain subject to the cell-containment caution above)."
        qui reg `y' i.`i' i.`t' `fe'   `trendinjector' `contreatcontrolinserter' `contreatcontrolinserterFE'   `ccc_absorbinserter'   i.`_i'#i.`_t'#i.`targetflag'`contreatinj'  `controls'    `interinject'    `coarsecontrolinserter'  `weightinserter'   if   `ifstatement' ,  `vcespecifier'
        cap est drop `wdfm'
        qui est sto `wdfm'
        cap drop `esamp'
        qui gen byte `esamp' = e(sample)
        local df_r = e(df_r)
      }
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
      local wdreal
      qui levelsof `_t' if `_t' > 0 & `_t' < ., local(wdtl)
      foreach v of local wdtl {
        local wdreal `wdreal' `=`v' - `min_t' - 1'
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
  if "`safety'" != "off"  {
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
              local problem_t = `problem_t' - `min_t' - 1
              dis as error  "Error: a treated cohort-period dummy coefficient was dropped for multicollinearity."
              if "`verbose'" != "" dis as error  "The specific problematic coefficient is: ``testbvar''"
              if "`verbose'" != "" dis as error  "For reference, _i is `_i', _t is `_t', and targetflag is `targetflag'"
              dis as error  "The problem seems to be with the treatment dummy for `i' == `problem_i', relative-`t' == `problem_t'"
              dis as error  "Specify contreatwithin if you have a continuous treatment intensity variable and wish to partial out any effect of treatment when the treatment intensity is 0."
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
        local amevoid 1
        local amescrub 1
        local amescrublo 1
      }
    }
  }




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
    local rownamessum
    if "`contreat'" == "" local sumvarlist y
    if "`contreat'" != "" local sumvarlist y contreat
    * gstats sum accepts pweights directly, so the regression weights pass through to these descriptive computations unchanged
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
    tempname summarystats
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
      * the ringer (wooldid's double copy of the dose) is averaged the same way: wherever
      * contreatcontrolinserterFE is non-empty it is an explicit regressor whose level shifts
      * the round-off of every prediction
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
  if "`semielasticity'" == "eydx" local effecttype Semi-elasticities
  if "`contreatelasticitytype'" == "eyex" local Ceffecttype Elasticities
  if "`contreatelasticitytype'" == "eydx" local Ceffecttype Semi-elasticities (eydx)
  if "`contreatelasticitytype'" == "dyex" local Ceffecttype Semi-elasticities (dyex)

  * for histograms of i-level effects, we do the following:
  * ask for a margins estimation to give us treated cohort specific effects without SEs
  * stick them into a matrix
  * use the matrix to populate a dataset that we will use to make the histogram
  if "`histogramcohorteffects'" != "" {
    if "`verbose'" == "verbose" dis as text "Wooldid is now producing histograms."
    tempname ___hist
    * graph window management: graphs are created with name(..., replace) so re-runs overwrite wooldid's own prior graphs; wooldid no longer clears the user's unrelated graphs from memory

    * we pool att with att_i b/c they conceptually should match into the same histograms
    if "`att_i'" == "att_i" |  "`att'" == "att"  {
      local whichweight w_att
      if "`att'" == "" local whichweight w_att_i
      qui margins [`weight'=``whichweight''], `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext' `fcnoes'
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
          if "`marginsdof'" != "nose" {
            tempname f2bmat f2vmat
            matrix `f2bmat' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']
            matrix `f2vmat' = r(V)[`matrixstartgrab'..`=colsof(r(b))', `matrixstartgrab'..`=colsof(r(b))']
            local f2k `horizitem'
            local f2b_`f2k' `f2bmat'
            local f2v_`f2k' `f2vmat'
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
        if "`cfxplotnose'" != "" local cfxmatendpoint 1
        * F9: fast mode keeps the cfx matrices in the established 1-column (point
        * estimates only) shape that cfxplotnose already produces, rather than minting
        * a new 6-column-with-missings shape
        * under permute() the 6-column shape is kept (missing-padded) so the shell's cell
        * rewiring can fill p/ll/ul from the replications
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
            if "`marginsdof'" != "nose" {
              tempname f2bmat f2vmat
              matrix `f2bmat' = r(b)[1, 1..`=colsof(r(b))']
              matrix `f2vmat' = r(V)[1..`=colsof(r(b))', 1..`=colsof(r(b))']
              local f2k = cond("`horizitem'"=="att","ame","ame_it")
              local f2b_`f2k' `f2bmat'
              local f2v_`f2k' `f2vmat'
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
            qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups') `poismarginstext' `marginsdof' `fcnoes'
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

          if "`impvalues'" == "" {
            * EFF3: see the posttreatment twin. the pre-period comparison pass is
            * capture-wrapped (pwcompare hard-errors when the pre-period subpopulation
            * realizes a single subgroup), so its failure branch runs this call as a
            * fallback - the deferral must never be able to lose the level margins.
            if wordcount("`sglist'") > 1 local sgeff3pre_`item' 1
            else {
              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck over(`subgroups') `poismarginstext' `marginsdof' `fcnoes'
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
              matrix `___`horizitem'_horiz_s`sgcase'' =  `f9tab'[1..6, `matrixstartgrab'..`=colsof(r(b))']'
              * F2: capture event-study b/V (relyr-labeled) for honestdid/lincom/test [pure-add]
              if "`marginsdof'" != "nose" {
                tempname f2bmat f2vmat
                matrix `f2bmat' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']
                matrix `f2vmat' = r(V)[`matrixstartgrab'..`=colsof(r(b))', `matrixstartgrab'..`=colsof(r(b))']
                local f2k `horizitem'_s`sgcase'
                local f2b_`f2k' `f2bmat'
                local f2v_`f2k' `f2vmat'
                local f2keys `f2keys' `f2k'
              }
              *the above will be updated w/ rownames / relyrs by the below

              * SGs need not have the same number of relative periods appear for each
              * this code block below figures out what relative periods appear for them
              tempname nameget
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
              local sgkeep
              local sgdroplbl
              local sgndrop 0
              local sgidx 0
              foreach sgnm of local names {
                local ++sgidx
                if strpos("`sgnm'", "o.") != 0 {
                  local ++sgndrop
                  if regexm("`sgnm'", "^([0-9]+)") {
                    local sgdrelyr = real(regexs(1)) - `min_t' - 1
                    if `sgdrelyr' < 0 local sgdroplbl `sgdroplbl' pre`=abs(`sgdrelyr')'
                    else              local sgdroplbl `sgdroplbl' post`sgdrelyr'
                  }
                }
                else local sgkeep `sgkeep' `sgidx'
              }
              if `sgndrop' > 0 {
                dis as text "Caution: the relative period effect(s) `sgdroplbl' are not estimable for subgroup `sgcase' (margins reports them as omitted), so they are dropped from that subgroup's event study; its remaining relative periods, and all other subgroups, are unaffected. This typically arises when the specified fixed effects absorb the subgroup-by-time variation (for example fe() interacting the subgroups() variable with time), or when a subgroup has no usable observations at that relative period."
                if "`sgkeep'" == "" {
                  dis as text "Caution: no relative period effect is estimable for subgroup `sgcase', so its event study is skipped entirely."
                  cap matrix drop `___`horizitem'_horiz_s`sgcase''
                  continue
                }
                local sgnames
                foreach sgi of local sgkeep {
                  local sgw : word `sgi' of `names'
                  local sgnames `sgnames' `sgw'
                }
                local names `sgnames'
                tempname sgkm
                local sgfirst 1
                foreach sgi of local sgkeep {
                  if `sgfirst' {
                    matrix `sgkm' = `___`horizitem'_horiz_s`sgcase''[`sgi', 1...]
                    local sgfirst 0
                  }
                  else matrix `sgkm' = `sgkm' \ `___`horizitem'_horiz_s`sgcase''[`sgi', 1...]
                }
                matrix `___`horizitem'_horiz_s`sgcase'' = `sgkm'
                matrix drop `sgkm'
                if "`marginsdof'" != "nose" {
                  tempname sgkb sgkv sgkrow
                  local sgfirst 1
                  foreach sgi of local sgkeep {
                    if `sgfirst' {
                      matrix `sgkb' = `f2bmat'[1, `sgi']
                      local sgfirst 0
                    }
                    else matrix `sgkb' = (`sgkb', `f2bmat'[1, `sgi'])
                  }
                  local sgfirst 1
                  foreach sgi of local sgkeep {
                    local sgfirst2 1
                    foreach sgj of local sgkeep {
                      if `sgfirst2' {
                        matrix `sgkrow' = `f2vmat'[`sgi', `sgj']
                        local sgfirst2 0
                      }
                      else matrix `sgkrow' = (`sgkrow', `f2vmat'[`sgi', `sgj'])
                    }
                    if `sgfirst' {
                      matrix `sgkv' = `sgkrow'
                      local sgfirst 0
                    }
                    else matrix `sgkv' = `sgkv' \ `sgkrow'
                  }
                  matrix `f2bmat' = `sgkb'
                  matrix `f2vmat' = `sgkv'
                  matrix drop `sgkb' `sgkv' `sgkrow'
                }
                * the adversarial gate counts margins columns; discount the removed ones
                local finaladvtest = `finaladvtest' - `sgndrop'
              }

              local names = subinstr("`names'", "bn","",.)
              local names = subinstr("`names'", "`_t'","",.)
              local names = subinstr("`names'", "1.`targetflag'","",.)
              if wordcount("`names'") == 0 {
                * how do we end up in this if statement?
                * margins will not report the levels of _t examined if there is only 1 level of _t examined
                * instead, it will return 1.`targetflag' as the relevant level (because the if made the over(_t) statement irrelevant by restricting to a sg with only one _t)
                * when we remove 1.`targetflag', there should be nothing left in names, meaning we need to figure out what the missing _t level is and put it back into names
                * why not just by default figure out what the correct _t levels are?
                * this approach requires a levelsof() statement, which can take longer than the text processing for large datasets
                qui  glevelsof(`_t')  if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' & `_t' != 0, clean local(names)
              }

              local relyrbysg
              local relyrlist_s`sgcase'
              foreach relyrsg in `names' {
                local relyrinj = `relyrsg' -`min_t' -1
                if `relyrinj' >= 0 local relyrlist_s`sgcase' `relyrlist_s`sgcase'' post`relyrinj'
                if `relyrinj' == 0 & `oldsyntax' ==1 local relyrlist_s`sgcase' `relyrlist_s`sgcase'' contemp

                if `relyrinj' < 0 {
                  local absinj = abs(`relyrinj')
                  local relyrlist_s`sgcase' `relyrlist_s`sgcase'' pre`absinj'
                }
                * per-estimand copy: att and att_it can keep DIFFERENT horizon sets for a
                * subgroup (SI-7 drops), and the rowname filter below names each ES block
                * from its own estimand's list; the unsuffixed local (last estimand
                * computed) still serves the contreat ES rows
                local relyrlist_s`sgcase'_`horizitem' `relyrlist_s`sgcase''
                if "`relyrbysg'" == "" {
                  local relyrbysg  `relyrinj'
                }
                else{
                  local relyrbysg `relyrbysg' ,  `relyrinj'
                }
              }
              tempname relyrcolumn_s`sgcase'
              matrix `relyrcolumn_s`sgcase'' = (`relyrbysg')'
              matrix `___`horizitem'_horiz_s`sgcase'' = ( `___`horizitem'_horiz_s`sgcase'',`relyrcolumn_s`sgcase'')
              local f2rl_`horizitem'_s`sgcase'
              forval f2ri = 1/`=rowsof(`relyrcolumn_s`sgcase'')' {
                local f2rv = `relyrcolumn_s`sgcase''[`f2ri',1]
                if `f2rv' < 0  local f2rl_`horizitem'_s`sgcase' `f2rl_`horizitem'_s`sgcase'' pre`=abs(`f2rv')'
                else          local f2rl_`horizitem'_s`sgcase' `f2rl_`horizitem'_s`sgcase'' post`f2rv'
              }

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
              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck over(`subgroups') `poismarginstext' `marginsdof' `fcnoes'
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
                matrix `___`horizitem'_Choriz_s`sgcase'' = ( `f9tab'[1..6, 1..`=colsof(r(b))']',`relyrcolumn_s`sgcase'')
                * AMESE R1 scrub: null se/t/p/ci columns of the per-sg ame ES rows
                if `amescrub' {
                  forval amer = 1/`=rowsof(`___`horizitem'_Choriz_s`sgcase'')' {
                    forval amej = `amescrublo'/6 {
                      matrix `___`horizitem'_Choriz_s`sgcase''[`amer', `amej'] = .
                    }
                  }
                }
                * F2: capture event-study b/V (relyr-labeled) for honestdid/lincom/test [pure-add]
                if "`marginsdof'" != "nose" {
                  tempname f2bmat f2vmat
                  matrix `f2bmat' = r(b)[1, 1..`=colsof(r(b))']
                  matrix `f2vmat' = r(V)[1..`=colsof(r(b))', 1..`=colsof(r(b))']
                  local f2k = cond("`horizitem'"=="att","ame","ame_it") + "_s`sgcase'"
                  local f2b_`f2k' `f2bmat'
                  local f2v_`f2k' `f2vmat'
                  local f2keys `f2keys' `f2k'
                  local f2rl_`f2k'
                  forval f2ri = 1/`=rowsof(`relyrcolumn_s`sgcase'')' {
                    local f2rv = `relyrcolumn_s`sgcase''[`f2ri',1]
                    if `f2rv' < 0  local f2rl_`f2k' `f2rl_`f2k'' pre`=abs(`f2rv')'
                    else          local f2rl_`f2k' `f2rl_`f2k'' post`f2rv'
                  }
                }
                * handle the adversarial pre-trends test
                * NFADV-01/02/05: statistics + display line via wooldid_advstat (see the
                * overall binary site); a gap-cancelled fit only skips THIS subgroup
                local finaladvtest
                cap confirm matrix r(b)
                if _rc == 0 {
                  local finaladvtest  `=colsof(r(b))'
                }
                if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
                if "`adversarial'" == "adversarial" & "`marginsdof'" != "nose" & `amescrub' == 0  {
                  if  `finaladvtest' > 1  {
                    wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn_s`sgcase'') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') sg(`sgcase') fm(`wdfm')
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
        local sglistlength = wordcount("`sglist'")
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
                  * include an adjustment in case only 1 sg group exists
                  if strpos("`word'","`subgroups'") == 0 & `sglistlength' == 1 local truerow `wordcounter'
                  local ++wordcounter
                }
                if "`truerow'" == "" local sghaspre no
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
            if "`att'" == "att" {
              matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \  `___att_Choriz_s`sgcase''
              matrix drop `___att_Choriz_s`sgcase''
            }
            if "`att_i'" == "att_i"  |  "`att_it'" == "att_it" |  "`att_itime'" == "att_itime" {
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
            * repeating its last name) or errors r(503). the contreat (ES_ame*) rows keep
            * the unsuffixed list, as before
            local sgeslist `relyrlist_s`sgcase''
            if strpos("`rownametest'","ES_att:") == 1    local sgeslist `relyrlist_s`sgcase'_att'
            if strpos("`rownametest'","ES_att_it:") == 1 local sgeslist `relyrlist_s`sgcase'_att_it'
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


  * return the command line input and estimate_type
  * F2: expose captured event-study b/V (honestdid/lincom/test compatible) [pure-add, after post]
  foreach f2ky of local f2keys {
    local f2bm `f2b_`f2ky''
    local f2vm `f2v_`f2ky''
    local f2rl `f2rl_`f2ky''
    if "`f2bm'" != "" & "`f2rl'" != "" {
      if `=colsof(`f2bm')' == `:list sizeof f2rl' & `=colsof(`f2bm')' > 0 {
        cap mata: st_matrix("`f2bm'", st_matrix("`f2bm'"))
        cap mata: st_matrix("`f2vm'", st_matrix("`f2vm'"))
        matrix colnames `f2bm' = `f2rl'
        matrix colnames `f2vm' = `f2rl'
        matrix rownames `f2vm' = `f2rl'
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
        local f2k = colsof(`f2vm')
        forval f2d = 1/`f2k' {
          if `f2vm'[`f2d',`f2d'] > 0 & `f2vm'[`f2d',`f2d'] < . local f2ok 1
          else local ++f2nbad
        }
        if `f2ok' {
          * SI-A4/R6: under a voided ame block the captured points are structural
          * zeros, so b_es_ame* is withheld along with V_es_ame* (do not trust
          * f2ok=0 to do this - it was observed on the audited d3 runs, not
          * guaranteed across the class)
          if !(`amevoid' & strpos("`f2ky'", "ame") == 1) cap ereturn matrix b_es_`f2ky' = `f2bm'
          * AMESE R1: withhold e(V_es_ame*) under cell containment — the captured V is
          * not the sampling covariance of the ame estimates (honestdid/lincom/test
          * consumers would ingest it silently); the points still post above
          local f2vrc 111
          if !(`amescrub' & strpos("`f2ky'", "ame") == 1) {
            cap ereturn matrix V_es_`f2ky' = `f2vm'
            local f2vrc = _rc
          }
          cap ereturn scalar numpre_`f2ky' = `f2np'
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
          if "`joint_horiz_s`sgcase''" != "" ereturn local joint_horiz_s`sgcase' = `joint_horiz_s`sgcase''
          if "`joint_horiz_i_s`sgcase''" != "" ereturn local joint_horiz_i_s`sgcase' = `joint_horiz_i_s`sgcase''
          if "`joint_Choriz_s`sgcase''" != "" ereturn local  joint_Choriz_s`sgcase' = `joint_Choriz_s`sgcase''
          if "`joint_Choriz_i_s`sgcase''" != "" ereturn local  joint_Choriz_i_s`sgcase' =  `joint_Choriz_i_s`sgcase''
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
  ereturn hidden local wdid_vcespec `"`vcespecifier'"'
  ereturn hidden local wdid_poisson "`poisson'"
  ereturn hidden local wdid_poisexp "`poisexpresults'"
  ereturn hidden local wdid_contreat "`contreat'"
  ereturn hidden local wdid_sgvar "`subgroups'"
  ereturn hidden local wdid_sglist "`sglist'"
  ereturn hidden local wdid_estimands "`att' `att_it' `att_i' `att_itime'"
  ereturn hidden local wdid_cwm "`customweightmultiplier'"
  ereturn hidden local wdid_semielas "`semielasticity'"
  ereturn scalar level = `wdlvl'
  ereturn hidden local wdid_fcol "`fcon'"
  * F061 (2026-09-06): the realised cell count (missing when the collapse did not run, so
  * "no collapse" and "collapsed to zero cells" cannot be confused) and the skip reason(s)
  if "`fconn'" != "" ereturn hidden scalar wdid_fcncell = `fconn'
  else ereturn hidden scalar wdid_fcncell = .
  ereturn hidden local wdid_fcskip "`fcskip'"
  ereturn hidden local wdid_fcguard "`fcguard'"
  ereturn hidden local wdid_cetype "`contreatelasticitytype'"
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
  ereturn local cmd "wooldid"




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
syntax [varlist(numeric default=none)] [if/]  [aw fw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist numeric) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose fastcollapse fullprecisionslow GTEffects PERMute(string) PERMREPS(integer 0) PERMEXact PERMTStat PERMTci PERMSTRata(varname numeric) PERMPRobs(string asis) PERMSAve(string) PERMANalyticse PERMABsp PERMSide(string) PERMRAWp BOOTstrap(string) reps(integer 0) ties Level(real 0) dots(integer -1) rseed(string) BOOTALPha(string) BOOTDf(string) BOOTDISPlay(string) BOOTTStat BOOTSAve(string)]

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
    if "`safety'" != "off" dis as text "Caution: espostlength set to 1 due to use of subgroups in conjunction with esprelength > 0."
    if "`safety'" != "off" dis as text "Without this modification, subgroups may crash if pre-treatment effects are not estimable for them."
    if "`safety'" != "off" dis as text "Specify espostlength(-999) to override this issue. Specify safety(off) to avoid seeing this message."
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
    if "`permstray'" != "" {
      dis as error "Error: `permstray' require(s) the permute() option."
      error 198
    }
  }
  * permtci's critical value comes from the permuted |t| distribution, so it implies
  * permtstat. Resolved HERE (it used to live in the shell alone) so that the fast-mode
  * default further down sees the implied value: permtci alone must not send the observed
  * run down the no-SE path it is about to be studentized against
  if "`permtci'" != "" & "`permtstat'" == "" local permtstat permtstat
  if "`permute'" != "" {
    cap noi wooldid_permctl , permute(`permute') permreps(`permreps') `permexact' clus(`cluster') ttrevar(`ttre_asgiven') absp(`permabsp') side(`permside') rawp(`permrawp') permstrata(`permstrata') permprobs(`permprobs') impv(`impvalues') uncond(`unconditionalse') cw(`contreatwithin') os(`oldsyntax') ts(`permtstat') tc(`permtci') fn(`fastnose') cc(`contreatcontrols') ifr(`if') ptp(`pretrendpower') ssg(`suppresssgprimaryeffects') sg(`subgroups') ivar(`i') tvar(`t') si(`safety') psave(`permsave') level(`level') `replace'
    if _rc != 0 error _rc
  }

  * capture the user's own safety() input before any internal relaxation (contreatwithin sets
  * the working copy to off further down); the identifying-assumption checks gate on this copy
  * so an internal relaxation does not silently disable them
  local safetyinput `safety'

  * F-boot: bootstrap-SE screens; deep checks live in the subprogram (wooldid_bootctl)
  * to conserve this program's compiled-size headroom. called UNCONDITIONALLY: with
  * bootstrap() absent it enforces the stray-companion screen and exits. it sits AFTER
  * the safetyinput capture above deliberately - the bayes cell-weighted gate keys on
  * the user's own safety() value (contreatwithin rewrites the working copy later)
  cap noi wooldid_bootctl , bootstrap(`bootstrap') reps(`reps') `ties' level(`level') dots(`dots') rseed(`rseed') bootalpha(`bootalpha') bootdf(`bootdf') bdisp(`bootdisplay') bts(`boottstat') bsave(`bootsave') pois(`poisson') permute(`permute') impv(`impvalues') uncond(`unconditionalse') cvce(`customvce') fn(`fastnose') os(`oldsyntax') ssg(`suppresssgprimaryeffects') attit(`att_it') atti(`att_i') attitime(`att_itime') si(`safetyinput') clus(`cluster') ivar(`i') ifr(`if') jt(`jointtests') ptp(`pretrendpower') ttrevar(`ttre_asgiven') `replace'
  if _rc != 0 error _rc

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
    if "`cfxlattice'" != "" | "`makecfxplots'" != ""                            local fcbad `fcbad' makecfxplots/cfxlattice()
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
      if "`safety'" == "off" dis as text "Caution: contreatwithin is specified in conjunction with the byTime contreatcontrol type."
      if "`safety'" != "off" dis as error "Error: contreatwithin is specified in conjunction with the byTime contreatcontrol type."
      if "`verbose'" == "verbose" dis as text "Typically, this partials out some treatment effects and yields unidentified estimates."
      if "`safety'" != "off" {
        dis as error "To proceed anyway in a new estimation run, manually set safety(off) before proceeding."
        error 198
      }
    }
    local safety off
    * the above is required since we are intentionally making the base effects irrelevant
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
    if "`adversarial'" != "" {
      dis as error "Error: fastnose is not compatible with adversarial; the adversarial pre-trends test needs a variance-covariance matrix, which fast mode does not estimate."
      error 198
    }
    if "`jointtests'" != "" {
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
  foreach wdv_ret in att att_it att_itime att_i cluster customvce unconditionalse robust impvalues fe coarsecohortcontrols ccc_absorb controls interactivecontrols timetrends subgroups customweightmultiplier esprelength espostlength saveplots poisson poisexpresults oostreatedcontrols anticipation esfixedbaseperiod esrelativeto safety makeplots makecfxplots summarystats histogramcohorteffects regtol jointtests suppresssgprimaryeffects contreat contreatpoly lattice latticeignoreweights ccnoabsorb contreatcontrols verbose cfxplottypes cfxlattice makecfxplotsbysg cfxplotnose semielasticity contreatelasticitytype contreatwithin cleanmatrices emptycellsoverride replace altmethod_adversarial adversarial update oldsyntax pretrendpower fastnose fastcollapse fullprecisionslow fcon fcskip gteffects permute permreps permexact permtstat permtci permstrata permprobs permsave permanalyticse permabsp permside permrawp y i t ttre_asgiven safetyinput suppressjointtests latticenum latticetype esposttrip cfxpost cfxpre cfxrelyrlist bootstrap reps ties level dots rseed bootalpha bootdf bootdisplay boottstat bootsave esprelengthraw {
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
  * convention (1/(permreps + 1) under permabsp), so below 39 replications a 5% rejection
  * is unreachable by construction and every row prints the floor as if it were evidence
  if `permreps' > 0 & "`permexact'" == "" {
    local wdp_floor = cond("`absp'" != "", 1/(`permreps' + 1), 2/(`permreps' + 1))
    if `wdp_floor' > .05 dis as text "Caution: with permreps(`permreps') the smallest attainable permutation p-value is " %5.3f `wdp_floor' " (" cond("`absp'" != "", "1", "2") "/(permreps + 1)), above .05, so no result can reject at the 5% level whatever the data say, and a p-value printed at that floor is the resolution limit rather than evidence. Use at least " cond("`absp'" != "", "20", "40") " replications for a 5% test (roughly 100 for .02, 200 for .01 and 1000 for .002)."
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


*** F6: the permutation engine. builds the cluster->onset map, generates assignments
*** (mata: wooldid_permx), runs the rep loop as nested quiet wooldid calls on a
*** permuted onset column, and hands the summary statistics back to the shell via
*** c_local (the shell restores the observed results and posts them to e())
cap prog drop wooldid_bootctl
prog def wooldid_bootctl
version 16.0
*** F-boot (plan: review_process/bootstrap_plan.md): bootstrap-SE option screens.
*** commit-B state: the full gate set is live, but a passing call still exits r(198)
*** loudly at the tail - the engine lands in a later commit and the option must never
*** silently no-op before then (the F6 commit-A precedent).
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

  * companion-value validation
  if `level' != 0 & (`level' < 10 | `level' > 99.99) {
    dis as error "Error: level() must be between 10 and 99.99."
    error 198
  }
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

  * C-stage refusals. jointtests: the bootstrap Wald joint tests land in a later
  * commit; a silent analytic/bootstrap mixed surface is worse than a loud refusal
  * (they will be exempt from the SI-A1 screen when they land - its own text
  * recommends bootstrap). pretrendpower: permanent v1 refusal (plan A13) - the ptp
  * subprogram computes and posts from the ANALYTIC e(V_es_*) inside the worker,
  * before any rewiring. (A cleanmatrices refusal was shipped here and in permctl and
  * has been REMOVED: the r(3499) it guarded against only reproduces when wooldid is
  * loaded with `do <path>`, as this project's harness does; an adopath-resolved ado
  * recovers because Stata re-reads it when a mata function is missing, so the refusal
  * was blocking a combination that works for installed users - see SI-17/SI-18.)
  * standard mode: engine lands in the next
  * commit; loud refusal until then (F6 precedent).
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
    local wdbsn = _N
  }
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
        local wdbcellmiss 0
        if "`bts'" != "" mata: wooldid_bootcellwire("`gt`fam''", 1, "`cell'_gt", "`cell'_gtse", 5, 6, 7, 8, 9, "`fam'", `level', `zc', `wdbdf', ("`ties'" != ""))
        else            mata: wooldid_bootcellwire("`gt`fam''", 1, "`cell'_gt", "", 5, 6, 7, 8, 9, "`fam'", `level', `zc', `wdbdf', ("`ties'" != ""))
        if "`fam'" == "percentile" local wdbcellmissall = `wdbcellmissall' + `wdbcellmiss'
        matrix colnames `gt`fam'' = `: colnames e(gteffects)'
        if "`fam'" == "`bdisp'" matrix `gtsel' = `gt`fam''
        ereturn hidden matrix boot_`ftag'_gteffects = `gt`fam''
      }
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
      mata: M = st_matrix("`cf`fam''"); M[., 3] = M[., 1] :/ M[., 2]; st_matrix("`cf`fam''", M)
      matrix colnames `cf`fam'' = `: colnames e(`nm')'
      matrix rownames `cf`fam'' = `: rowfullnames e(`nm')'
      if "`fam'" == "`bdisp'" matrix `cfsel' = `cf`fam''
      * e() names are capped at 32 characters; a cfx matrix name too long for the twin
      * prefix keeps only its selected family in the live matrix
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
  local kk 0
  local rnames : rowfullnames `wm'
  foreach rn of local rnames {
    local ++kk
    gettoken eqn nmn : rn, parse(":")
    local nmn = subinstr("`nmn'", ":", "", 1)
    local wdbeq`kk' `eqn'
    local wdbnm`kk' `nmn'
    local wdbsg`kk' .
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
          if "`sgtag'" == "" ereturn scalar `en' = `jp'
          else ereturn local `en' = `jp'
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
      * per subgroup (ereturn LOCALS - the analytic macro convention)
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
        if `jp' < . ereturn local `ename'_s`sg' = `jp'
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
            matrix `pwobs0' = `pw0'[1..., 2]'
            local eqwant0 = cond("`pre'" == "", "`it'", "pre_`it'")
            mata: wooldid_bootpairidx("`eqwant0'", "`sgl'", "`pairidx'")
            mata: wooldid_bootsgw("`acc'", "`seacc0'", "`pwseacc'", "`obsb'", "`aseobs'", "`pwobs0'", "`cols'", "`pairidx'", "wdbjp", "wdbjn")
          }
          if `wdbjp' < . ereturn scalar joint_`pre'`it'_allsgssame = `wdbjp'
          else {
            ereturn scalar joint_`pre'`it'_allsgssame = .
            local ++wdbncaut
            local wdbc`wdbncaut' the bootstrapped `pre'`it' subgroup-equality test is unavailable: forming its studentized statistic needs the observed and per-replication pairwise-contrast standard errors, and at least one of these is missing (the pairwise margins standard error is unavailable in this specification, e.g. under clustering at or above the i level)
          }
        }
        else if "`e(joint_`pre'`it'_allsgssame)'" != "" {
          * boottstat requested but the columns could not be mapped: still never leave
          * the analytic value standing
          ereturn scalar joint_`pre'`it'_allsgssame = .
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
      local ++wdbncaut
      local wdbc`wdbncaut' the across-subgroup equality tests were not computed, for the same reason as the joint tests: bootstrapping them means bootstrapping their test statistic, which needs boottstat. They are reported as absent rather than as analytic tests
    }
  }

  * -------- pairwise subgroup contrasts: keep the analytic point differences,
  * rewire se/z/p/ci from the bootstrap covariance of the matching estimand columns
  cap confirm matrix e(pwsgs)
  if _rc == 0 {
    tempname pw pwsd pwpct pwbc pwbca pwtp pwts pwt1
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
          if `wdbse' < . & `wdbse' > 0 & `wdpalb' < . & `wdpaub' < . {
            matrix `pwbca'[`rr', 2] = `wdbse'
            matrix `pwbca'[`rr', 3] = `pwbca'[`rr', 1]/`wdbse'
            matrix `pwbca'[`rr', 4] = `wdpap'
            matrix `pwbca'[`rr', 5] = `wdpalb'
            matrix `pwbca'[`rr', 6] = `wdpaub'
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

  if `ndone' < max(20, `wdp_R'/2) & "`permexact'" == "" {
    dis as error "Error: only `ndone' of `wdp_R' permutation replications estimated successfully; too few for inference. Check whether reassigned onset configurations are frequently inestimable on these data."
    error 459
  }
  * NFPERM-01 (lean residue): under permexact ANY failure means the enumeration is
  * incomplete, so the caution fires at nfail > 0 with wording that states the
  * substance; sampled mode keeps the 10% threshold (adjudicated accounting design)
  if `nfail' > 0 & "`permexact'" != "" dis as text "Caution: `nfail' of `wdp_R' enumerated assignments failed to estimate, so the exact enumeration is incomplete; the reported exact p-values are shares over completed assignments only (the observed-is-a-member and 2/M-floor readings no longer apply)."
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
  syntax , acc(string) obsb(string) cmpb(string) [seacc(string) cmpse(string) level(real 0) m(string) done(integer 0) fail(integer 0) blank(integer 0) ts(integer 0) method(string) reps(integer 0) ase(string) clus(varlist) gt(string) absp(string) side(string) rawp(string)]
  if `level' == 0 local level = c(level)
  * p-value convention code (F052/F064, 2026-09-06): 0 equal-tailed (default), 1 absolute
  * add-one (permabsp), 2 raw |s| share c/n without add-one (permrawp: the ritest / Stata
  * permute number), 3 left one-sided, 4 right one-sided (permside())
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
  mata: wooldid_permstat("`acc'", "`seacc'", "`obsb'", "`cmpb'", "`cmpse'", "`excl'", `level', `wdpM', `wdpabs', "`kern'")
  local obn : colfullnames `obsb'

  * -------- rewire the results matrices: col 3 = permuted t-ratio (permtstat) or blank,
  * col 4 = permutation p, cols 5-6 = shift interval. col 2 stays as posted (analytic se
  * or blank). Blanked rows are counted for the cautions
  local wdpmiss 0
  local wdplow 0
  local wdpdeg
  local wdptsmiss 0
  local wdpoff
  tempname wma
  matrix `wma' = `wm'
  local kk 0
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
  if `level' != 95 {
    local wmc : subinstr local wmc "lb_95ci" "lb_`wdplvltxt'ci", word
    local wmc : subinstr local wmc "ub_95ci" "ub_`wdplvltxt'ci", word
  }
  matrix colnames `wm' = `wmc'
  ereturn hidden matrix perm_analytic_mainresults = `wma'
  ereturn matrix wdidmainresults = `wm'
  foreach sg of local sgl {
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
    if `level' != 95 {
      local wsc : subinstr local wsc "lb_95ci" "lb_`wdplvltxt'ci", word
      local wsc : subinstr local wsc "ub_95ci" "ub_`wdplvltxt'ci", word
    }
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
      mata: M = st_matrix("`gtm'"); M[., 6..9] = J(rows(M), 4, .); st_matrix("`gtm'", M)
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
  ereturn hidden local wdid_fastnose ""

  * -------- cautions (estimation-time; the footer carries the standing disclosures)
  if `wdpmiss' > 0 dis as text "Caution: `wdpmiss' reported coefficient(s) have no usable permutation distribution (fewer than 2 completed replications carry them`=cond(`ts', ", or no usable observed standard error to form the t-ratio with", "")'); their p-value and interval are left blank rather than filled from any other source."
  if `wdplow' > 0 dis as text "Caution: `wdplow' reported coefficient(s) were present in fewer than half of the completed replications (a reassigned onset can make a horizon or a subgroup level inestimable in a replication); their permutation p-values rest on the counts in e(perm_n)."
  if "`wdpdeg'" != "" dis as text "Caution: the permutation distribution is degenerate (every completed draw equal) for: `wdpdeg'; those permutation p-values are uninformative."
  if "`wdpoff'" != "" dis as text "Caution: the permutation distribution is not centred at zero (its central `level'% band excludes 0) for: `wdpoff'. The shift interval takes the estimate's sampling error to have that null distribution, so it is not interpretable there and is left blank; the permutation p-value stands. This happens when a quantity is non-zero under reassigned onsets - for example the ame rows under lattice() with contreatcontrols(basic), where the baseline dose slope is not partialled out (see the standing issue on the reference cell) - or under a permute(probabilities) design far from the observed assignment."
  if `blank' > 0 dis as text "Caution: `blank' replication cell(s) in ame-family columns returned exact structural zeros (reassignments without dose variation in some cells) and were excluded from those columns' permutation distributions."
  if "`wdp_incomplete'" != "" dis as text "Caution: under permexact the following elements were not estimable in every enumerated assignment, so their p-values use the add-one form over the completed assignments rather than exact shares: `wdp_incomplete'"
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
  syntax , acc(string) obsb(string) cmpb(string) [level(real 0) m(string) ts(integer 0) jt(string) adv(string) advalt(string) advesrel(string) advtt(string) esf(string) absp(string) side(string) rawp(string)]
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
        mata: wooldid_permqf("`acc'", "`cmpb'", "`cols'", "`Cm'", `wdpM', "jp", "jn", "jw", "`Wr'")
        local en = cond("`sgtag'" == "", "`ename'", "`ename'_s`sgtag'")
        if `jp' < . {
          if "`sgtag'" == "" ereturn scalar `en' = `jp'
          else ereturn local `en' = `jp'
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
  * the POST-treatment event-study horizons - the max-|studentized| step-down (Westfall-Young
  * / Romano-Wolf) from the same accumulator, exact under the sharp null: every horizon's
  * draws are studentized on the fixed metric (the SD of its null draws, as the kernel's
  * row 6), the horizons are ordered by |t_obs|, and the j-th adjusted p is the share of
  * replications whose max |t_r| over the horizons not yet rejected reaches |t_obs(j)|,
  * made monotone. Posted per family as e(perm_stepdown[_i][_s<sg>]) with rows t_obs,
  * p_raw (the per-horizon absolute-value share), p_singlestep (max over all post
  * horizons) and p_stepdown, one named column per horizon; computed whenever at least
  * two post horizons exist (jointtests is not required)
  foreach fam in att att_it ame ame_it {
    local sname = cond(strpos("`fam'", "ame"), "perm_Cstepdown", "perm_stepdown")
    if "`fam'" == "att_it" | "`fam'" == "ame_it" local sname `sname'_i
    foreach sgtag in "" `sgl' {
      local sgw = cond("`sgtag'" == "", ".", "`sgtag'")
      local cols
      forval k = 1/`K' {
        if "`wdpeq`k''" == "ES_`fam'" & "`wdpsg`k''" == "`sgw'" & `wdprl`k'' >= 0 & `wdprl`k'' < . local cols `cols' `k'
      }
      if `: word count `cols'' < 2 continue
      tempname SD
      mata: wooldid_permstep("`acc'", "`cmpb'", "`cols'", `wdpM', "`SD'")
      cap confirm matrix `SD'
      if _rc != 0 continue
      if rowsof(`SD') < 4 continue
      local en = cond("`sgtag'" == "", "`sname'", "`sname'_s`sgtag'")
      if strlen("`en'") <= 32 ereturn matrix `en' = `SD'
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
          mata: wooldid_permqf("`acc'", "`cmpb'", "`cols'", "`Cm'", `wdpM', "jp", "jn", "jw", "`Wr'")
          if `jp' < . {
            ereturn scalar joint_`pre'`it'_allsgssame = `jp'
            matrix colnames `Wr' = joint_`pre'`it'_allsgssame
            tempname Wo
            matrix `Wo' = J(1, 1, `jw')
            matrix colnames `Wo' = joint_`pre'`it'_allsgssame
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
          mata: wooldid_permadv("`acc'", "`cmpb'", "`cols'", "`rels'", `advesrel', ("`advtt'" != ""), `alt', `wdpM', "jp", "jn", "jw", "`Wr'")
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
          mata: wooldid_permpw("`acc'", "`cmpb'", `ka', `kb', st_matrix("`pw'")[`rr', 1], `level', `wdpM', `wdpabs', "`o1'")
          matrix `pk'[1, `rr'] = `o1'
          local ok 1
        }
      }
      matrix `pw'[`rr', 3] = .
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
  syntax , cell(string) cell0(string) [ncfx(integer 0) names(string) level(real 0) m(string) ts(integer 0) lattype(string) cfxlat(string) absp(string) side(string) rawp(string)]
  if `level' == 0 local level = c(level)
  * p-value convention code (F052/F064, 2026-09-06): 0 equal-tailed (default), 1 absolute
  * add-one (permabsp), 2 raw |s| share c/n without add-one (permrawp: the ritest / Stata
  * permute number), 3 left one-sided, 4 right one-sided (permside())
  local wdpabs = cond("`rawp'" != "", 2, cond("`side'" == "left", 3, cond("`side'" == "right", 4, ("`absp'" != ""))))
  local wdpM = cond("`m'" == "" | "`m'" == ".", ".", "`m'")
  local drift = (inlist("`lattype'", "baseTreated", "baseTreatedPost", "baseTreatedPre", "baseTreatedRef") | inlist("`lattype'", "byTreated", "byTreatedByTime", "byTtre", "byTtreByTime")) & "`cfxlat'" == ""
  local ndone 0
  forval c = 1/`ncfx' {
    local nm : word `c' of `names'
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    if colsof(e(`nm')) < 6 continue
    cap confirm matrix `cell'_cfx`c'
    if _rc != 0 continue
    tempname cf cfa
    matrix `cf' = e(`nm')
    matrix `cfa' = `cf'
    if `drift' {
      mata: M = st_matrix("`cf'"); M[., 4..6] = J(rows(M), 3, .); if (!`ts') M[., 3] = J(rows(M), 1, .); st_matrix("`cf'", M)
    }
    else {
      mata: wooldid_permcellwire("`cf'", "`cell'_cfx`c'", "`cell0'_cfx`c'", "`=cond(`ts', "`cell'_cfx`c'se", "")'", "`=cond(`ts', "`cell0'_cfx`c'se", "")'", 1, 2, 4, 5, 6, `level', `wdpM', `wdpabs')
      mata: M = st_matrix("`cf'"); M[., 3] = (`ts' ? M[., 1] :/ M[., 2] : J(rows(M), 1, .)); st_matrix("`cf'", M)
    }
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
    dis as text "Caution: the lattice type `lattype' bins the dose on treated units' rows, so its bins are functions of the treatment assignment and differ from replication to replication; the continuous-treatment plot cells therefore have no permutation distribution and their inference columns are blank (points only). latticetype all, byCohort, byTime, byCohortByTime, or a user-supplied cfxlattice() variable, give assignment-invariant bins that are permuted."
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
*** weighted estimand); under bootalpha(size) the mass is split over the cluster's base
*** weight (row multiplier = draw_c / W_c), which matches bayesboot BITWISE on unweighted
*** data. The rwgen priorpowers(-n_c) oracle for equal is its weights times n_c.
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

  * map frame: one row per cluster with n_c (row count) and W_c (base-weight total);
  * frozen generation-time row ids (NFPERM-02: frlink re-sorts the linked frame)
  tempname mbfr
  frame put `clusv' `basew' if `smp', into(`mbfr')
  frame `mbfr' {
    qui bysort `clusv': gen double __wdbn = _N
    qui gegen double __wdbW = total(`basew'), by(`clusv')
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
  * consumption order; row multiplier = (g/colsum)/W_c. standard: draws happen
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
  if "`wtype'" == "" & "`pois'" == "" & "`seaccname'" != "" & "`clus'" == "" & "`mode'" == "standard" local wtype2 fweight
  local repcall `y' `i' `t' `ttrevar'
  tempvar repw lnk drawv

  * rep-0 comparator: fastnose run of the OBSERVED weights (A3/NFPERM-10 precedent;
  * BC/ties in a later commit compare replicates against THESE values, and the
  * master-name sanity check below guards the name-aligned accumulation)
  local wtif0
  if "`wtype'" != "" local wtif0 [`wtype'=`wexp']
  if `"`ifr'"' != "" local wtif0 `wtif0' if (`ifr')
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
  * observed pre-horizon coefficients stacked here while the observed e() is still live.
  * Keys mirror wooldid_bootwire's family walk: att att_it ame ame_it and their _s<sg>
  * twins. A key with no observed e(b_es_*) simply contributes no column
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
      local wdb_proj = r(t100)/5*`reps'
      timer clear 100
      * no "." in the text: the dots line is counted by tests and by eye
      local wdb_ptxt = cond(`wdb_proj' < 90, strofreal(round(`wdb_proj')) + " s", strofreal(round(`wdb_proj'/60)) + " min")
      if `dots' > 0 dis as text " (~`wdb_ptxt' projected for `reps' replications)" _continue
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
      if `dots' > 0 {
        if `badrep' dis as error "x" _continue
        else dis as text "." _continue
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
  syntax , [mode(string) reps(integer 0) acc(string) obsb(string) done(string) fail(string) blank(string) level(real 0) bootalpha(string) clus(varlist) comp(string) TIES rseed(string) bootdf(string) bdisp(string) nclus(integer 0) seacc(string) aseout(string) jkacc(string)]
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
    if `excl'[1, `kk'] == 0 & `bse'[1, `kk'] < . & `bse'[1, `kk'] > 0 & `ciba'[1, `kk'] < . & `ciba'[2, `kk'] < . {
      matrix `wmbca'[`rr', 2] = `bse'[1, `kk']
      matrix `wmbca'[`rr', 3] = `wmbca'[`rr', 1]/`bse'[1, `kk']
      matrix `wmbca'[`rr', 4] = `cbap'[1, `kk']
      matrix `wmbca'[`rr', 5] = `ciba'[1, `kk']
      matrix `wmbca'[`rr', 6] = `ciba'[2, `kk']
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
  ereturn hidden matrix boot_pctl_mainresults = `wmpct'
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
      if `excl'[1, `kkP'] == 0 & `bse'[1, `kkP'] < . & `bse'[1, `kkP'] > 0 & `ciba'[1, `kkP'] < . & `ciba'[2, `kkP'] < . {
        matrix `wsbca`sg''[`rr', 2] = `bse'[1, `kkP']
        matrix `wsbca`sg''[`rr', 3] = `wsbca`sg''[`rr', 1]/`bse'[1, `kkP']
        matrix `wsbca`sg''[`rr', 4] = `cbap'[1, `kkP']
        matrix `wsbca`sg''[`rr', 5] = `ciba'[1, `kkP']
        matrix `wsbca`sg''[`rr', 6] = `ciba'[2, `kkP']
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
  cap ereturn repost V = `vnew'
  local wdb_vposted = _rc == 0
  if !`wdb_vposted' dis as text "Caution: the bootstrap covariance could not be posted to e(V) (rc `=_rc'); it remains available as e(V_boot). Reported tables carry the bootstrap standard errors either way."
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
  }
  if "`ties'" != "" ereturn local boot_ties ties

  ereturn hidden local wdid_vcespec "`wdbmodetext', `wdbcl', reps(`reps')"
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
  local wdbdfr0 = e(df_r)
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
    if `wdbdfr0' < . ereturn hidden scalar df_r_analytic = `wdbdfr0'
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
    * coefficient is individually significant (using wooldid's own critical value)
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
      cap mata: st_local("ptpnissize", strofreal(1 - mvnormalcvqp(-`ptpcrit'*st_matrix("`ptpsd'"), `ptpcrit'*st_matrix("`ptpsd'"), J(1, `ptpk', 0), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
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
            cap mata: st_local("ptppw2", strofreal(1 - mvnormalcvqp(-`ptpcrit'*st_matrix("`ptpsd'"), `ptpcrit'*st_matrix("`ptpsd'"), `ptpghi'*st_matrix("`ptpmrow'"), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
            if _rc != 0 local ptpiter 999
            else local ptppw = real("`ptppw2'")
            local ++ptpiter
          }
          if `ptpiter' < 999 {
            forval ptpbi = 1/80 {
              local ptpgmid = (`ptpglo' + `ptpghi')/2
              cap mata: st_local("ptppw2", strofreal(1 - mvnormalcvqp(-`ptpcrit'*st_matrix("`ptpsd'"), `ptpcrit'*st_matrix("`ptpsd'"), `ptpgmid'*st_matrix("`ptpmrow'"), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
              if _rc == 0 {
                if real("`ptppw2'") < `ptplev'/100 local ptpglo = `ptpgmid'
                else local ptpghi = `ptpgmid'
              }
            }
            local ptpgn = (`ptpglo' + `ptpghi')/2
            cap mata: st_local("ptppw2", strofreal(1 - mvnormalcvqp(-`ptpcrit'*st_matrix("`ptpsd'"), `ptpcrit'*st_matrix("`ptpsd'"), `ptpgn'*st_matrix("`ptpmrow'"), vech(st_matrix("`ptpS'"))', 1000), "%21.0g"))
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

* F038 (quality review 2026-09-05, author ruling 2026-09-06): a rank-based bootstrap p of
* exactly 0 means no completed replication fell on the far side of the null - the
* resolution limit 2/n, not a probability of 0. matlist renders a numeric matrix, so the
* cell keeps the 0 (e() keeps it too, by the ruling); this line under the table says what
* the 0 means. Bootstrap only: the permutation add-one p never prints 0, and the Wald
* display's p is a tail probability
cap prog drop wooldid_pfloor
prog def wooldid_pfloor
  args mat
  if "`e(boot_completed)'" == "" exit
  if "`e(boot_display)'" == "wald" & "`e(boot_tstat)'" != "1" exit
  cap confirm matrix `mat'
  if _rc != 0 exit
  tempname pm
  matrix `pm' = `mat'
  local pc = colnumb(`pm', "p")
  if `pc' >= . local pc 4
  if `pc' > colsof(`pm') exit
  mata: st_local("wdp0", strofreal(sum(st_matrix("`pm'")[., `pc'] :== 0)))
  if `wdp0' > 0 dis as text "  p = 0 above means p < " %5.3f 2/e(boot_completed) " (2/" e(boot_completed) " completed replications, the resolution of the rank-based p): no completed replication fell on the far side of the null. The cell and e() keep the 0."
end

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
  local estimandlist `e(wdid_estimands)'
  local customweightmultiplier `e(wdid_cwm)'
  local semielasticity `e(wdid_semielas)'
  local contreatelasticitytype `e(wdid_cetype)'
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
    matlist e(summarystats)

    if "`poisson'" == "poisson" dis as result  "Note: summary stats above are for untransformed Y, while coefficients are from Y=exp(regression equation) and have a log-linear type interpretation."
  }


  ************* now show main results with matlist
  dis as result  " "
  dis as result  "Main`poistext' Results (Full Estimation Sample): `elasinj'"
  if "`customweightmultiplier'" != "" dis as result  "Note: all Estimates Shown Apply `customweightmultiplier' Weights Post-Estimation"
  matlist e(wdidmainresults)
  wooldid_pfloor e(wdidmainresults)

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
     if "`e(adversarial_pval_`typetext')'" != "" {
       if "`e(perm_method)'" == "" dis as result "p-value for joint test of null that all (non-reference) `typetext' coefficients are on a single line: " %-6.4f e(adversarial_pval_`typetext')
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
        matlist e(wdidsgresults_`sgcase')
        wooldid_pfloor e(wdidsgresults_`sgcase')
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
            if "`e(adversarial_pval_`typetext'_s`sgcase')'" != "" {
              if "`e(perm_method)'" == "" dis as result "p-value for joint test of null that all (non-reference) `typetext' coefficients for sg=`sgcase' are on a single line: " %-6.4f e(adversarial_pval_`typetext'_s`sgcase')
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
      matlist e(pwsgs)
      wooldid_pfloor e(pwsgs)
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
      local wdpex = cond("`e(perm_exact)'" == "1", "exact shares over the enumerated assignments (the observed assignment is a member of its own reference set)", "the add-one Monte-Carlo form (1 + count)/(replications + 1)")
      if "`e(perm_pconv)'" == "absolute" {
        local wdppc two-sided absolute-value permutation p-values (permabsp): the share of replications whose absolute statistic is at or above the observed absolute statistic, `wdpex'
        local wdpiv symmetric shift intervals dual to that p-value - the estimate plus or minus the `wdplv'th-percentile absolute
      }
      else if "`e(perm_pconv)'" == "raw" {
        local wdppc two-sided absolute-value permutation p-values WITHOUT the add-one adjustment (permrawp): c/n, the share of replications whose absolute statistic is at or above the observed absolute statistic - the number ritest and Stata's permute prefix report
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
      dis as result "  Monte-Carlo precision of the p-values: e(perm_p_se) holds the standard error of each p-value's underlying share and e(perm_p_ci) its `wdplv'% Clopper-Pearson interval (doubled under the equal-tailed rule); the one-sided shares are in e(perm_p_lo) and e(perm_p_hi)."
      local wdpqf 0
      cap confirm matrix e(perm_jointW)
      if _rc == 0 local wdpqf 1
      cap confirm matrix e(perm_sgW)
      if _rc == 0 local wdpqf 1
      cap confirm matrix e(perm_advW)
      if _rc == 0 local wdpqf 1
      if `wdpqf' dis as result "  The joint pre-trend tests, across-subgroup equality tests and adversarial line tests shown are one-sided permutation tests of a quadratic-form statistic in the permutation-covariance metric (the covariance of the replications' coefficients), so they need no per-replication covariance; the adversarial statistic is wooldid's line-fit statistic on that metric, and its p-value tests the sharp null of no effect at any horizon rather than the composite on-a-line null of the analytic test. Replication statistics: e(perm_jointW), e(perm_sgW), e(perm_advW)."
      * F065: multiplicity-adjusted p-values across the post-treatment horizons
      foreach wdpsdn in perm_stepdown perm_stepdown_i perm_Cstepdown perm_Cstepdown_i {
        cap confirm matrix e(`wdpsdn')
        if _rc != 0 continue
        local wdpsdf = cond("`wdpsdn'" == "perm_stepdown", "att", cond("`wdpsdn'" == "perm_stepdown_i", "att_it", cond("`wdpsdn'" == "perm_Cstepdown", "ame", "ame_it")))
        dis as result "  Multiplicity across the post-treatment `wdpsdf' horizons (max-|studentized| step-down, Romano-Wolf; the studentization is the SD of each horizon's null draws): horizon: raw p -> step-down p"
        tempname wdpsdm
        matrix `wdpsdm' = e(`wdpsdn')
        local wdpsdc : colnames `wdpsdm'
        local wdpsdk 0
        local wdpsdl
        foreach wdpsdh of local wdpsdc {
          local ++wdpsdk
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
      if "`e(perm_cfx_blank)'" == "1" dis as result "  The continuous-treatment plot cells were not permuted (assignment-dependent lattice bins); their inference columns are blank."
      if "`e(perm_incomplete)'" != "" dis as result "  Elements not estimable in every enumerated assignment (add-one p-values rather than exact shares): `e(perm_incomplete)'"
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
   local wdbmt = cond("`e(boot_method)'" == "bayes", "Bayesian cluster bootstrap (Dirichlet reweighting, alpha = `e(boot_alpha)')", cond(strpos("`e(wdid_vcespec)'", "observation-level") == 0, "cluster bootstrap (resampling with replacement)", "observation-level bootstrap (resampling with replacement)"))
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
     if "`wdbdsp'" == "bca" dis as result "  The BCa acceleration comes from a delete-one-cluster jackknife (" e(boot_jk_completed) " fits; e(boot_accel)); its intervals are in e(boot_ci_bca) and p-values in e(boot_p_bca). The gteffects/cfx cells carry the bias-corrected family."
     dis as result "  A percentile interval is not symmetric about the estimate, so the estimate plus or minus a critical value times the reported standard error will not in general reproduce the reported interval - this is the ordinary reading of a bootstrap table, not a wooldid peculiarity. Use bootdisplay(wald) for the normal/t package in which they do agree."
   }
   dis as result "  Every family is computed whichever is displayed: intervals in e(boot_ci_normal)/e(boot_ci_percentile)/e(boot_ci_bc), p-values in e(boot_p_percentile)/e(boot_p_bc), and the complete tables in the hidden e(boot_wald_mainresults)/e(boot_pctl_mainresults)/e(boot_bc_mainresults)."
   if "`wdbdsp'" != "wald" | "`e(boot_tstat)'" == "1" dis as result "  The rank-based coefficient p-values have resolution 2/" e(boot_completed) " = " %5.3f 2/e(boot_completed) ": a p printed as 0 means p below that (no completed replication on the far side of the null); it is stored as 0 in e(). The joint, equality and adversarial p-values are add-one shares over their own listwise-complete replication counts (their smallest value is 1/(n + 1), never 0)."
   if e(boot_v_posted) == 1 dis as result "  e(V) holds the bootstrap covariance (complete-replication listwise, e(boot_N_V) replications; the analytic covariance, when one existed, remains in the hidden e(V_analytic))."
   else dis as result "  The bootstrap covariance is stored in e(V_boot); no covariance is posted to e(V) for this run."
   * e(V) and e(df_r) carry the sd-family inference whichever family is displayed, so a
   * hand-run test/lincom reproduces the sd package and NOT a displayed percentile p
   if "`wdbdsp'" != "wald" | "`e(boot_tstat)'" == "1" dis as result "  e(V) and e(df_r) carry the standard-deviation family referred to the `wdbrefd' distribution, so test and lincom run by hand reproduce the bootdisplay(wald) numbers rather than the columns printed above."
   if "`e(boot_jt_family)'" == "fixedmetric" dis as result "  The joint pre-trend tests (and the across-subgroup equality tests, when subgroups are present) are rank tests of the replication draws recentred at the observed estimates on the fixed bootstrap-covariance metric, p = add-one share of replications at or beyond the observed statistic (e(boot_jt_family) = fixedmetric; the replication statistics are in e(boot_jointW) and e(boot_sgW)). Under boottstat they would instead be bootstrapped test statistics."
   if "`e(boot_jt_absent)'" != "" | "`e(boot_sgjt_absent)'" != "" dis as result "  Some joint tests and/or across-subgroup equality tests are NOT reported for this run (no usable replication statistic - too few complete replications or a rank-deficient bootstrap covariance of the tested columns). They are absent rather than shown as analytic tests, which would not be bootstrap quantities."
   if "`e(boot_tstat)'" == "1" dis as result "  Subgroup contrasts, when requested, are studentized by their own analytic pairwise standard error in every replication (boottstat) and carry the same package as the coefficient rows."
   else dis as result "  Subgroup contrasts, when requested, are bootstrapped as coefficient differences and carry the same package as the coefficient rows."
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
* wooldid_worker for per-program size headroom. Two modes:
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
      dis as text "Caution: estimation returned an all-zero variance-covariance matrix. The issue is likely related to continuous treatment effect estimation coupled with use of analytic clustered SEs, with clusters occuring at a high level of aggregation. Consider switching to use of bootstrap SEs, IM p-values, or randomization inference if finer clustering levels are non-viable."
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
  syntax [, SAVEplots(name) replace]
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
      if `nc' >= 6 {
        qui replace se = `M'[`r', 2] in `r'
        qui replace p = `M'[`r', 4] in `r'
        qui replace ll = `M'[`r', 5] in `r'
        qui replace ul = `M'[`r', 6] in `r'
      }
    }
    sort lattice
    qui count if ll < . & ul < . & b < .
    local band = r(N) > 0
    local ttl `"`e(wdid_cfxttl`k')'"'
    local ytx `"`e(wdid_cfxytx`k')'"'
    if `band' twoway (rarea ll ul lattice, color(navy%50)) (connected b lattice, color(navy)), `xt' `xlab' ytitle(`"`ytx'"') title(`"`ttl'"') scale(.85) legend(off) name("`nm'", replace)
    else twoway (connected b lattice, color(navy)), `xt' `xlab' ytitle(`"`ytx'"') title(`"`ttl'"') scale(.85) legend(off) name("`nm'", replace)
    if "`saveplots'" != "" {
      qui graph save "`saveplots'_`nm'.gph", `replace'
      tempfile part
      qui save `part'
      frame `af' {
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
    dis as error "Error: wooldid_esplots is internal to wooldid and needs a stored wooldid estimation (replay contract '1') in e(); run wooldid first."
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
      twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr,color(navy%50)) if (rownamesvar =="reference" | rownamesvar=="`graph'")  , `xlineinserter' xtitle("Time Relative to Time of Initial Treatment") ytitle("``graphtext'_text' on `y' by Horizon") legend(off)   name(`graphtext'_fullmodel, replace) `xlabelstarter'
      if "`saveplots'" != "" qui graph save "`saveplots'_`graph'.gph", `replace'
    }


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

        sort relyr
        foreach graph in `graphtypes' {
          local graphtext = subinstr("`graph'","ES_","",.)
          local att_it_text ATT_it `elasinj'
          local att_text ATT `elasinj'
          local ame_it_text AME_it `elasinj2'
          local ame_text AME `elasinj2'
          twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr,color(navy%50)) if (rownamesvar =="reference" | rownamesvar=="`graph'")  , `xlineinserter' xtitle("Time Relative to Time of Initial Treatment") ytitle("``graphtext'_text' on `y' by Horizon") legend(off)   name(`graphtext'_sg`sgcase', replace) `xlabelstarter'   title("`subgroups' == `sgcase'")
          if "`saveplots'" != "" qui graph save "`saveplots'_`graph'_sg`sgcase'.gph", `replace'
        }

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
*** for the event-study levels under `levels', which nothing contrasts) has its largest
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
    real scalar C, si, lnM, M, R, r, k, capexact, capwarn
    string scalar nps
    st_local("wdp_err", "")
    st_local("wdp_caut1", "")
    st_local("wdp_caut2", "")
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
    A = J(C, R, .)
    if (exact) {
        // enumerate per-stratum multiset permutations, then take the cartesian product
        A = wooldid_permx_enum(ons, str)
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
*** alt == 0, primary (NFADV-01): minimum-distance test of the composite null that
*** all horizon effects lie on SOME line, E[b] = X gamma with X = [1, relyr] and
*** gamma unknown. gamma is estimated by GLS under a generalized inverse of V; the
*** minimized distance W = min_g (b - Xg)' V^+ (b - Xg) is computed in the whitened
*** eigenspace of V (z = L^-1/2 U' b, M = L^-1/2 U' X, W = ||z - M ghat||^2), which
*** is exactly chi-squared with q = rank(V) - rank(M) degrees of freedom under the
*** null. The planted zero-variance reference rows receive zero weight under V^+ and
*** so are mathematically inert here (they are constants equal to their expectation);
*** omitting them from the computation is identical to including them.
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
        X = J(rows(rel), 1, 1), rel
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
        M = T*X
        sv = svdsv(M)
        rankM = sum(sv :> max(sv)*1e-12)
        q = k1 - rankM
        if (q <= 0) return(.)
        g = svsolve(M, z)
        r = z - M*g
        W = r'r
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
mata:
// draw matrix: per-rep (outer) per-cluster-ascending (inner) Gamma(alpha,1) draws -
// the exact rwgen bayes stream-consumption order; alpha = n_c (size) or 1 (equal);
// row multiplier = (g/colsum)/W_c. writes a C x R matrix (transposed for el() reads)
void wooldid_bootdw(string scalar fr, string scalar amode, real scalar R, string scalar outname)
{
  string scalar cur
  real matrix N, W, E, P
  real scalar C, r, k, a
  cur = st_framecurrent()
  st_framecurrent(fr)
  N = st_data(., "__wdbn")
  W = st_data(., "__wdbW")
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
  // bootalpha(size) only: split the cluster's mass over its base weight so that every row
  // carries p_c / W_c - the bayesboot / rwgen priorpowers convention. Under bootalpha(equal)
  // (the default) every row carries p_c itself, so the cluster's TOTAL replicate mass is
  // W_c * p_c: Rubin's Bayesian bootstrap with clusters as the exchangeable units, whose
  // replications stay centred on the observation-weighted estimand. Until 2026-09-06 the
  // split applied under equal too, which centred the replications on the equal-cluster-
  // weighted estimand (rep mean .80 vs att 1.69 on a 12-state unbalanced fixture)
  if (amode == "size") P = P :/ W
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
// the rank of t_obs among the t_r WITHOUT add-one (see below: the p and the critical
// value share one order-statistic rule, the one wooldid_bootci already uses).
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
  // [b_obs - q_hi*se_obs, b_obs - q_lo*se_obs], NOT the earlier symmetric b +/- q(|t|)*se;
  // the p is the rank of t_obs among the t_r, two-sided, with NO add-one and the SAME
  // order-statistic rule wooldid_bootci uses for the coefficient families - so the p and
  // the interval are exact duals up to the even-count averaging boundary measured there
  // (this also retires the two defects where the p was add-one while the critical value
  // was a plain order statistic, and where a third quantile rule coexisted). The former
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
    out[3, k] = min((1, (2*min((m0, m1)))/n))
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
      p = min((1, (2*min((m0, m1)))/n))
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
    // V_r on the same names, same order
    q = 0
    for (j = 1; j <= rows(on); j++) {
      if (on[j, 1] != fk) continue
      q++
      R[q, .] = J(1, np, .)
    }
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
// construction as wooldid_boottp - equal-tailed signed quantiles, rank p with no add-one,
// wooldid_bootci's order-statistic rule - so a contrast row and a coefficient row are
// studentized by one rule (no BC-t twin: F003). Output: a 6 x 1 column with the same row
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
  out[3] = min((1, (2*min((m0, m1)))/n))
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
    av[k] = a
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
real scalar wooldid_permq(real colvector v, real scalar pr)
{
  real scalar n, x
  n = rows(v)
  x = n*pr
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1, n))])/2)
  return(v[max((min((ceil(x), n)), 1))])
}

real rowvector wooldid_permtail(real colvector draws, real scalar sobs, real scalar exact,
                                real scalar al, real scalar zc, real scalar absp)
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
  if (exact) {
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
    if (absp == 1) out[3] = (exact ? ma/n : (1 + ma)/(n + 1))
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
  qq = cc/n
  out[16] = mult*sqrt(qq*(1 - qq)/n)
  al2 = al
  cpl = (cc <= 0 ? 0 : invibeta(cc, n - cc + 1, al2))
  cph = (cc >= n ? 1 : invibeta(cc + 1, n - cc, 1 - al2))
  out[17] = min((1, mult*cpl))
  out[18] = min((1, mult*cph))
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
                      real scalar lvl, real scalar M, real scalar absp, string scalar outname)
{
  real matrix A, S, out
  real rowvector ob, cm, cs, ex, row
  real colvector s
  real scalar K, k, r, al, zc, sobs, ts, exk
  string scalar lowav, incom
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
  for (k = 1; k <= K; k++) {
    if (ex[1, k] != 0) continue
    if (!ts) {
      s = A[., k]
      sobs = cm[1, k]
    }
    else {
      s = J(rows(A), 1, .)
      for (r = 1; r <= rows(A); r++) {
        if (A[r, k] < . & S[r, k] < . & S[r, k] > 0) s[r] = A[r, k]/S[r, k]
      }
      if (cs[1, k] >= . | cs[1, k] <= 0 | cm[1, k] >= .) sobs = .
      else sobs = cm[1, k]/cs[1, k]
    }
    exk = (M < . & sum(s :< .) == M)
    if (M < . & sum(s :< .) != M) incom = incom + " " + st_matrixcolstripe(accname)[k, 2]
    row = wooldid_permtail(s, sobs, exk, al, zc, absp)
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
}

// contrast kernel: d_r = b_r[ka] - b_r[kb] over listwise-complete rows, ranked against
// the comparator's difference; the interval about the displayed difference dobs
void wooldid_permpw(string scalar accname, string scalar cmpname, real scalar ka, real scalar kb,
                    real scalar dobs, real scalar lvl, real scalar M, real scalar absp,
                    string scalar outname)
{
  real matrix A
  real rowvector cm, row
  real colvector ok, d
  real scalar al, zc, sobs, exk
  A = st_matrix(accname)[., (ka, kb)]
  cm = st_matrix(cmpname)
  ok = rowmissing(A) :== 0
  d = select(A[., 1] :- A[., 2], ok)
  sobs = cm[1, ka] - cm[1, kb]
  al = (100 - lvl)/200
  zc = invnormal(1 - al)
  exk = (M < . & rows(d) == M)
  row = wooldid_permtail(d, sobs, exk, al, zc, absp)
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
                    string scalar Cname, real scalar M, string scalar outp, string scalar outn,
                    string scalar outwobs, string scalar outW)
{
  real matrix A, CC, Cm, R, Si, dev
  real rowvector idx, cm, sc, ev
  real colvector ok, Wr, d
  real scalar n, r, Wobs, nex, exk
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
  // exact enumeration with every assignment complete: the observed assignment is one of
  // the n rows already, so the comparator is NOT appended (it was counted twice until
  // 2026-09-05, which shrank the metric and moved the exact-mode joint p-values)
  exk = (M < . & n == M)
  CC = (exk ? select(A, ok) : select(A, ok) \ cm)
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
  if (exk) st_local(outp, strofreal(nex/n, "%21x"))
  else st_local(outp, strofreal((1 + nex)/(n + 1), "%21x"))
}

// F065 (2026-09-06): max-|studentized| step-down over a set of columns (the post-treatment
// horizons of one event-study family). Fixed-metric studentization: each column's null
// draws (replications, plus the comparator row unless the enumeration is complete and
// exact) are scaled by their standard deviation; the observed statistic is the comparator.
// Rows of the output: 1 |t_obs|, 2 p_raw (the absolute-value share per column, add-one
// or exact), 3 p_singlestep (max |t_r| over ALL columns), 4 p_stepdown (max over the
// columns not yet rejected, in decreasing |t_obs| order, made monotone). Columns keep
// the accumulator's names; the output is empty (0 x 0) when fewer than 2 complete
// replications exist or a column has zero dispersion
void wooldid_permstep(string scalar accname, string scalar cmpname, string scalar colslist,
                      real scalar M, string scalar outname)
{
  real matrix A, CC, T, out
  real rowvector idx, cm, sd, tobs, ord
  real colvector ok, mx
  real scalar n, H, j, h, exk, nex, pprev, pj
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
  exk = (M < . & n == M)
  CC = (exk ? select(A, ok) : select(A, ok) \ cm)
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
    out[2, h] = (exk ? nex/n : (1 + nex)/(n + 1))
  }
  mx = rowmax(T)
  for (h = 1; h <= H; h++) {
    nex = sum(mx :>= tobs[h])
    out[3, h] = (exk ? nex/n : (1 + nex)/(n + 1))
  }
  // step-down: order by |t_obs| descending; at step j the max runs over columns j..H
  // of that order; monotone in j
  ord = order(-tobs', 1)'
  pprev = 0
  for (j = 1; j <= H; j++) {
    mx = rowmax(T[., ord[j..H]])
    nex = sum(mx :>= tobs[ord[j]])
    pj = (exk ? nex/n : (1 + nex)/(n + 1))
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
                     real scalar M, string scalar outp, string scalar outn, string scalar outwobs,
                     string scalar outW)
{
  real matrix A, CC, V, dev
  real rowvector idx, cm
  real colvector ok, Wr, rel
  real scalar n, r, Wobs, nex, exk, q
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
  // same rule as wooldid_permqf: under complete exact enumeration the comparator is
  // already one of the rows
  exk = (M < . & n == M)
  CC = (exk ? select(A, ok) : select(A, ok) \ cm)
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
  nex = sum(select(Wr, ok) :>= Wobs)
  st_local(outwobs, strofreal(Wobs, "%21x"))
  if (exk) st_local(outp, strofreal(nex/n, "%21x"))
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
                          real scalar lvl, real scalar M, real scalar absp)
{
  real matrix Mx, A, S, C, CS
  string colvector mk
  string matrix cs
  real colvector s
  real rowvector row
  real scalar i, k, r, al, zc, b, sobs, ts, exk, so
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
    exk = (M < . & sum(s :< .) == M)
    row = wooldid_permtail(s, sobs, exk, al, zc, absp)
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
  real scalar c, idx
  string scalar nm
  Asg = st_matrix(assignname)
  cv = strtoreal(tokens(clusvals))
  oo = strtoreal(tokens(obsons))
  if (cols(cv) != rows(Asg) | cols(oo) != rows(Asg)) return
  for (c = 1; c <= rows(Asg); c++) {
    nm = strtoname("onset_" + strofreal(cv[c], "%20.0g"))
    if (_st_varindex(nm) < .) continue
    idx = st_addvar("double", nm)
    st_store(., idx, oo[c] \ Asg[c, .]')
  }
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

end
