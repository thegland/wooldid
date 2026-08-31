
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
      dis as error "Error: the last wooldid call failed with r(${WOOLDID_LASTRC}); the stored results predate that failed call. Refit wooldid before replaying."
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
  * worker set wdp_on; re-running wooldid_validate materializes every normalized
  * option here (the S5 c_local blanket), from which the rep command line is built.
  * the observed results are stored, the rep loop re-runs wooldid quietly on permuted
  * onset assignments, then the observed results are restored, the permutation
  * summaries are appended to e() (this shell is eclass), and part(perm) prints them
  if "`wdp_on'" == "1" {
    qui wooldid_validate `0'
    * permtci's critical value comes from the permuted |t| distribution, so it implies permtstat
    if "`permtci'" != "" & "`permtstat'" == "" local permtstat permtstat
    * rep option string: keep everything that shapes the point estimates; strip
    * display/diagnostic/plot work and the permutation options themselves; reps run
    * under fastnose (only coefficients are consumed). subgroups are stripped: v1
    * permutation statistics cover the pooled estimands only
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
    if "`contreat'" != "" local wdp_ro `wdp_ro' contreat(`contreat') contreatpoly(`contreatpoly')
    if `"`contreatcontrols'"' != "" local wdp_ro `wdp_ro' contreatcontrols(`contreatcontrols')
    if `"`lattice'"' != "" local wdp_ro `wdp_ro' lattice(`lattice')
    if "`semielasticity'" == "eydx" local wdp_ro `wdp_ro' semielasticity
    if "`contreat'" != "" & "`contreatelasticitytype'" != "" local wdp_ro `wdp_ro' contreatelasticitytype(`contreatelasticitytype')
    if `esprelength' != 0 local wdp_ro `wdp_ro' esprelength(`esprelength')
    if "`esposttrip'" == "no" local wdp_ro `wdp_ro' espostlength(-999)
    else if `espostlength' != 0 local wdp_ro `wdp_ro' espostlength(`espostlength')
    local wdp_ro `wdp_ro' esrelativeto(`esrelativeto')
    if "`safety'" != "" local wdp_ro `wdp_ro' safety(`safety')
    if `regtol' != 9 local wdp_ro `wdp_ro' regtol(`regtol')
    if `oostreatedcontrols' != 0 local wdp_ro `wdp_ro' oostreatedcontrols(`oostreatedcontrols')
    if `anticipation' != 0 local wdp_ro `wdp_ro' anticipation(`anticipation')
    local wdp_est `att' `att_it' `att_i' `att_itime'
    if "`wdp_est'" == "" local wdp_est att
    tempname wdpobs wdpdist wdpcomp
    qui est sto `wdpobs'
    cap noi wooldid_permrun , y(`y') i(`i') t(`t') ttrevar(`ttre_asgiven') clus(`cluster') permute(`permute') permreps(`permreps') `permexact' permstrata(`permstrata') permprobs(`permprobs') `permtstat' repopts(`wdp_ro') estimands(`wdp_est') wtype(`wdp_wt') wexp(`wdp_wx') ifr(`wdp_if') distname(`wdpdist') compname(`wdpcomp')
    local wdp_rc = _rc
    qui est restore `wdpobs'
    cap est drop `wdpobs'
    global WOOLDID_LASTRC
    if `wdp_rc' != 0 exit `wdp_rc'
    * append the permutation summaries to the restored observed results
    ereturn local perm_method `permute'
    ereturn hidden local wdid_permest `wdp_est'
    ereturn scalar perm_reps = `wdp_r_reps'
    ereturn scalar perm_completed = `wdp_r_done'
    ereturn scalar perm_failed = `wdp_r_fail'
    ereturn scalar perm_npossible = `wdp_r_npd'
    ereturn local perm_npossible_str `wdp_r_nps'
    ereturn local perm_stat = cond("`permtstat'" != "", "t", "b")
    foreach wdpe in `wdp_est' {
      ereturn scalar perm_p_`wdpe' = `wdp_r_p_`wdpe''
      ereturn scalar perm_lb_`wdpe' = `wdp_r_lb_`wdpe''
      ereturn scalar perm_ub_`wdpe' = `wdp_r_ub_`wdpe''
      ereturn scalar perm_sd_`wdpe' = `wdp_r_sd_`wdpe''
      if "`permtci'" != "" {
        * inverted-t CI: observed b +/- (95th pctile of the permuted |t|) x observed se
        local wdptb = _b[main_`wdpe']
        local wdpts = _se[main_`wdpe']
        if `wdpts' < . & `wdpts' > 0 & `wdp_r_q_`wdpe'' < . {
          ereturn scalar perm_tlb_`wdpe' = `wdptb' - `wdp_r_q_`wdpe''*`wdpts'
          ereturn scalar perm_tub_`wdpe' = `wdptb' + `wdp_r_q_`wdpe''*`wdpts'
        }
        else dis as text "Caution: the inverted-t CI for `wdpe' is unavailable (no usable observed standard error)."
      }
    }
    if `wdp_r_done' <= 10000 cap ereturn matrix permdist = `wdpdist'
    * NFPERM-12 (lean kernel): the comparison statistics actually used for the tail
    * counts, so p-values are reconstructable from e() alone. matrix-handle transport -
    * a c_local decimal round-trip re-parses one ulp off (see the rep-0 comment in
    * wooldid_permrun), which would break bitwise reconstruction
    cap ereturn matrix perm_comparator = `wdpcomp'
    wooldid_display , part(perm)
    dis as result " "
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
    * bootstrap options themselves; reps run under fastnose (only e(b) is consumed)
    local wdb_ro `att' `att_it' `att_itime' `att_i' `esfixedbaseperiod' `timetrends' `poisson' `poisexpresults' `robust' `ccnoabsorb' `latticeignoreweights' `emptycellsoverride' `contreatwithin' fastnose
    if "`cluster'" != "" local wdb_ro `wdb_ro' cluster(`cluster')
    if "`fe'" != "" local wdb_ro `wdb_ro' fe(`fe')
    if `"`coarsecohortcontrols'"' != "" local wdb_ro `wdb_ro' coarsecohortcontrols(`coarsecohortcontrols')
    if "`ccc_absorb'" != "" local wdb_ro `wdb_ro' ccc_absorb(`ccc_absorb')
    if "`controls'" != "" local wdb_ro `wdb_ro' controls(`controls')
    if "`interactivecontrols'" != "" local wdb_ro `wdb_ro' interactivecontrols(`interactivecontrols')
    if "`customweightmultiplier'" != "" local wdb_ro `wdb_ro' customweightmultiplier(`customweightmultiplier')
    if "`contreat'" != "" local wdb_ro `wdb_ro' contreat(`contreat') contreatpoly(`contreatpoly')
    if `"`contreatcontrols'"' != "" local wdb_ro `wdb_ro' contreatcontrols(`contreatcontrols')
    if `"`lattice'"' != "" local wdb_ro `wdb_ro' lattice(`lattice')
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
    tempname wdbobs wdbobsb wdbcomp wdbacc
    qui est sto `wdbobs'
    matrix `wdbobsb' = e(b)
    local wdb_reid = cond("`bootstrap'" == "standard" & "`att_it'`att_i'`att_itime'" != "", "reid", "")
    cap noi wooldid_bootrun , y(`y') i(`i') t(`t') ttrevar(`ttre_asgiven') mode(`bootstrap') `wdb_reid' reps(`reps') dots(`dots') rseed(`rseed') bootalpha(`bootalpha') clus(`cluster') repopts(`wdb_ro') wtype(`wdp_wt') wexp(`wdp_wx') ifr(`wdp_if') pois(`poisson') obsb(`wdbobsb') accname(`wdbacc') compname(`wdbcomp')
    local wdb_rc = _rc
    qui est restore `wdbobs'
    cap est drop `wdbobs'
    global WOOLDID_LASTRC
    if `wdb_rc' != 0 exit `wdb_rc'
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
    * rewire the posted inference surfaces to the bootstrap covariance (eclass sub)
    wooldid_bootwire , mode(`bootstrap') reps(`reps') level(`level') bootalpha(`bootalpha') clus(`cluster') acc(`wdbacc') obsb(`wdbobsb') comp(`wdbcomp') done(`wdb_r_done') fail(`wdb_r_fail') blank(`wdb_r_blank') rseed(`rseed') bootdf(`bootdf') nclus(`wdb_r_nclus')
    * bootstrapped joint tests, subgroup contrasts and CI flavors (eclass sub)
    wooldid_bootjnt , level(`level') acc(`wdbacc') obsb(`wdbobsb') comp(`wdbcomp') jt(`jointtests') `ties' done(`wdb_r_done') bootdf(`bootdf') nclus(`wdb_r_nclus')
    * deferred display: worker skipped these under bootstrap; they read only e(),
    * which now carries the bootstrap inference (identical formatting, free replay)
    wooldid_display , part(1)
    wooldid_display , part(2)
    wooldid_display , part(boot)
    if "`makeplots'`makecfxplots'`saveplots'`gteffects'`histogramcohorteffects'`adversarial'`altmethod_adversarial'" != "" dis as text "Caution: plots, saved plot data, gteffects cell standard errors, histogram outputs and adversarial line-test p-values are produced from the ANALYTIC covariance during estimation; the reported tables and e() inference above are bootstrap-based."
    dis as result " "
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
  * Every Stata replication prefix (bootstrap/jackknife/simulate/permute) runs its
  * replication loop under `set coeftabresults off', which silently removes margins
  * r(table) - and this worker harvests r(table) at 16+ sites, dying r(503)
  * "conformability error" on every replicate. The setting is program-scoped, so this
  * line is automatically undone at worker exit and is a no-op in a default session.
  set coeftabresults on

* S5: this line must remain BYTE-IDENTICAL to the syntax line in wooldid_validate
syntax [varlist(numeric default=none)] [if/]  [aw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist numeric) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose GTEffects PERMute(string) PERMREPS(integer 0) PERMEXact PERMTStat PERMTci PERMSTRata(varname numeric) PERMPRobs(string asis) BOOTstrap(string) reps(integer 0) ties Level(real 0) dots(integer -1) rseed(string) BOOTALPha(string) BOOTDf(string)]

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
      dis as error "Common causes: an onset recorded only from the onset period onward (fill it in on every in-sample row of a treated group, including pre-treatment rows); underlying units with different onsets combined under one `i' value; or IDs stored as float, which pools distinct large IDs (recast double `i'). Alternatively, restrict the estimation sample with if()."
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
    local min_t = abs(r(min))
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
    * interactinos between the de-meaned interactive controls and the dummies for treated units' treatment periods
    * (I think the null untreated dummy X interactive controls object included should be harmless to include / appears to exist in some of Wooldridge's example files)
    local interinject i.`_i'#i.`_t'#i.`targetflag'#c.(`interlist')  i.`t'#c.(`interactivecontrols')
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
    * this text initializes contreatcontrolinserterfe correctly if using basic w/o a lattice
    * and also initializes it this way if there is any chance that this FE object might get kicked outside of an HDFE absorb, in which case this section being included
    * is necessary in order to correctly manipulate the behavior of reg and poisson's colinearity drop rules (nothing in an HDFE is dropped b/c of colinearity with something not in the FE)
		* re: the nature of the manipulation, specifying null##ringer and lattice##ringer (as opposed to ringer alone or lattice#ringer) -- and doing so explicitly, rather than
		*     implicitly like in the other control types -- ensures that the base level of contreat gets partialled out. if not specified, it does not always get partialled out, and
		*     instead one of the contreat control coefficients gets dropped for colinearity instead
    if (strpos("`contreatcontrols'",  "basic") != 0 & "`latticepart2hash'" == "" ) | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
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
      dis as text "Caution: the cohort variable `i' takes `si8ilev' distinct values. wooldid builds cohort-by-period interaction terms from it, so very granular `i' variables can be slow and can exceed Stata's variable limits. You may wish to specify unit-level fixed effects via fe() and swap `i' to a coarser level of aggregation, such as cohorts or treatment onset groups. (safety(off) suppresses this check.)"
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
  if "`jointtests'" != "" & "`bootstrap'" == "" {
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
  qui est sto ___wooldidfullmodel
  global WOOLDID_CLEANUP $WOOLDID_CLEANUP fullmodel
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
        cap est drop ___wooldidfullmodel
        qui est sto ___wooldidfullmodel
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
      mata: wooldid_gtx("`_i'", "`_t'", "`targetflag'", "`gtraw'", "`gtsed'", `min_t', `gthstart', `gthend', "`df_r'", "`gtm'")
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
      qui est restore ___wooldidfullmodel
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
  local poismarginstext
  if "`poisson'" == "poisson" & "`poisexpresults'" == "" local poismarginstext expression(predict(xb))




  ******************************************************************************************
  *****************************         PRODUCE HISTOGRAMS       ***************************
  ******************************************************************************************

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
      qui margins [`weight'=``whichweight''], `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext'
      local matrixstartgrab floor(`=colsof(r(b))'/2)+1
      matrix `___hist' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
      matrix colnames `___hist' = att_i
      if "`contreat'" != "" {
        qui margins [`weight'=``whichweight''], `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext'
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
      qui margins [`weight'=`w_att_it'], `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext'
      local matrixstartgrab floor(`=colsof(r(b))'/2)+1
      tempname ___hist2
      matrix `___hist2' = r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
      matrix colnames `___hist2' = att_itime
      if "`contreat'" != "" {
        tempname ___Chist2
        qui margins [`weight'=`w_att_it'], `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `_i' != 0 & `t' >= `ttre') noestimcheck over(`_i') nose `poismarginstext'
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
      qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' `marginsdof'
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
      qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' nose over(`cluster')
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
      local im_ll = `b' + invt(r(df_t),.025)*`se'
      local im_ul = `b' - invt(r(df_t),.025)*`se'
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
        qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' `marginsdof'
        matrix `f9tab' = r(table)
        if rowsof(`f9tab') < 6 {
          matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
          matrix rownames `f9tab' = b se t pvalue ll ul
        }
        matrix `___`item'_pre' = `f9tab'[1..6,2]'
      }

      * this produces an IM p-values
      if "`impvalues'" == "impvalues" {
        qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' nose over(`cluster')
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
        local im_ll = `b' + invt(r(df_t),.025)*`se'
        local im_ul = `b' - invt(r(df_t),.025)*`se'
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
          qui margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_t' >= `horizonstarter' & `_t' <= `horizonend' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof'
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
          if "`adversarial'" == "adversarial"     {
            wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial')
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
            qui  margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `_t' == `horizon' & `targetflag'==1  & `_t' != 0) noestimcheck `poismarginstext' nose over(`cluster')
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
            local im_ll = `b' + invt(r(df_t),.025)*`se'
            local im_ul = `b' - invt(r(df_t),.025)*`se'
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
  if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'"=="" {
    if "`verbose'" == "verbose" dis as text "Wooldid is now computing pre-treatment horizon effect joint tests."
    if "`att'" == "att" {
      * contrast(overjoint) joint-tests the hypothesis that each estimate by group specified in over() is the same
      * since we have specified our subpop to include _t == 0, which is the omitted case, this is a joint test vs 0
      qui margins [`weight'=`w_att'] ,`semielasticity'(`targetflag') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & (`_t' == 0 | `_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' `jointmarginsdof'
      local joint_horiz = r(p)[1,2]
    }
    if "`joint_horiz_i'" == "" & (  "`att_i'" != "" | "`att_itime'" != "" | "`att_it'" != "" ){
      qui margins [`weight'=`w_att_it'],`semielasticity'(`targetflag') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . &  (`_t' == 0 | `_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' `jointmarginsdof'
      local joint_horiz_i = r(p)[1,2]
    }
  }


  if "`verbose'" == "verbose" dis as text "Wooldid has completed production of all main average treatment effects for the full sample."









  **********************************************************************************************************************************************
  *****************************         MAKE CONTINUOUS TREATMENT EFFECTS PLOTS AND CALCULATE CONTINUOUS TREATMENT EFFECTS      ****************************
  **********************************************************************************************************************************************
  * graph window management: graphs are created with name(..., replace) so re-runs overwrite wooldid's own prior graphs; wooldid no longer clears the user's unrelated graphs from memory
  local cfxmatricestoreturn

  * continuous treatment estimates are constructed borderline identically to the main estimates until we get to cfx plots
  if "`contreat'" != "" {

    foreach item in `att' `att_it' `att_i' `att_itime' {
      if "`verbose'" == "verbose" dis as text "Wooldid is now working on `item'-weighted average marginal effect estimation and production of any requested continuous treatment effects plots, including plots by subgroup."
      tempname ___C`item'
      * posttreatment average marginal effect of the continuous treatment variable generated here - just a matter of specifying contreat vs targetflag to get an AME of some sort
      * we also no longer need the skip-half-the-matrix adjustment, because w/ a continuous variable things are not computed at 0 and 1
      if "`impvalues'" == "" {
        qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' `marginsdof'
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
        qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' nose over(`cluster')
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
        local im_ll = `b' + invt(r(df_t),.025)*`se'
        local im_ul = `b' - invt(r(df_t),.025)*`se'
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
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' `marginsdof'
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
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' nose over(`cluster')
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
          local im_ll = `b' + invt(r(df_t),.025)*`se'
          local im_ul = `b' - invt(r(df_t),.025)*`se'
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
        if "`fastnose'" != "" local cfxmatendpoint 1
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
          qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstext'
          tempname cfx`item'_`nametext'
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `cfx`item'_`nametext'' = `f9tab'[1..`cfxmatendpoint',floor(`=colsof(r(table))'/2)+1..`=colsof(r(table))']'
          qui marginsplot, recastci(rarea) plotopts(color(navy)) ciopts(color(navy%50)) `cfx_xtitle'  scale(.85) ytitle("`ytext'") title("`titletext'") name("cfx`item'_`nametext'", replace)
          if "`saveplots'" != "" qui graph save "`saveplots'_cfx`item'_`nametext'.gph", `replace'

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
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstext'
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
          qui marginsplot, `amecfxci' plotopts(color(navy)) `cfx_xtitle' ytitle("`ytext'")  scale(.85) title("`titletext'") name("cfx`itemtext'_`nametext'", replace)
          if "`saveplots'" != "" qui graph save "`saveplots'_cfx`itemtext'_`nametext'.gph", `replace'
          local cfxmatricestoreturn `cfxmatricestoreturn' cfx`itemtext'_`nametext' cfx`item'_`nametext'
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

              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `subgroups' == `sgcase' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstext'
              tempname cfx`item'_`nametext'_sg`sgcase'
              matrix `f9tab' = r(table)
              if rowsof(`f9tab') < 6 {
                matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                matrix rownames `f9tab' = b se t pvalue ll ul
              }
              matrix `cfx`item'_`nametext'_sg`sgcase''  = `f9tab'[1..`cfxmatendpoint',floor(`=colsof(r(b))'/2)+1..`=colsof(r(b))']'
              qui marginsplot, recastci(rarea) plotopts(color(navy)) ciopts(color(navy%50)) `cfx_xtitle' ytitle("`ytext'") title("`titletext'") scale(.85) name("cfx`item'_`nametext'_sg`sgcase'", replace)
              if "`saveplots'" != ""  qui graph save "`saveplots'_cfx`item'_`nametext'_sg`sgcase'.gph", `replace'
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
              qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `subgroups' == `sgcase' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlatticeplot') `cfxseopt' `poismarginstext'
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
              qui marginsplot, plotopts(color(navy)) `amecfxci' `cfx_xtitle' ytitle("`ytext'") title("`titletext'")  scale(.85) name("cfx`itemtext'_`nametext'_sg`sgcase'", replace)
              if "`saveplots'" != "" qui graph save "`saveplots'_cfx`itemtext'_`nametext'_sg`sgcase'.gph", `replace'
              local cfxmatricestoreturn `cfxmatricestoreturn'   cfx`itemtext'_`nametext'_sg`sgcase' cfx`item'_`nametext'_sg`sgcase'
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
            qui  margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof'
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
            if "`adversarial'" == "adversarial" & `amescrub' == 0   {
              wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial')
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
              qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `_t' == `horizon' & `targetflag'==1  ) noestimcheck   over(`cluster') `poismarginstext' nose
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
              local im_ll = `b' + invt(r(df_t),.025)*`se'
              local im_ul = `b' - invt(r(df_t),.025)*`se'
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
      if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'"=="" & `amescrub' == 0 {
        if "`att'" == "att"  {
          qui margins [`weight'=`w_att'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != .  &  (`_t' == 0 |`_t' >= `horizonstarter')     ) contrast(overjoint) `poismarginstext' `marginsdof'
          local joint_Choriz = r(p)[1,1]
        }

        if "`joint_Choriz_i'" == "" & ( "`att_it'" != "" | "`att_itime'" != "" | "`att_i'" != "" ) {
          qui margins [`weight'=`w_att_it'],`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != .  &  (`_t' == 0 |`_t' >= `horizonstarter')  ) contrast(overjoint) `poismarginstext' `marginsdof'
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
          qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups') `poismarginstext' `marginsdof'
          local matrixstartgrab floor(`=colsof(r(b))'/2)+1
          matrix `f9tab' = r(table)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `___`item'_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
        }


        * pretreatment effects
        if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
          tempname ___`item'_pre_s
          local sgpre 1

          if "`impvalues'" == "" {
            qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck over(`subgroups') `poismarginstext' `marginsdof'
            local matrixstartgrab floor(`=colsof(r(b))'/2)+1
            matrix `f9tab' = r(table)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `___`item'_pre_s' = `f9tab'[1..6,`matrixstartgrab'..`=colsof(r(b))']'
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
              qui margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof'
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
              if "`adversarial'" == "adversarial"   {
                * SI-7 (author-ruled 2026-08-27): if any horizon was dropped as not estimable,
                * refuse the line test for THIS subgroup. The adversarial null is that the
                * coefficients lie on a single line in relative time; with a horizon missing the
                * fit is run on an incomplete - and possibly interior-holed - grid, so the
                * statistic would not mean what its name says. Other subgroups are unaffected.
                if `sgndrop' > 0 {
                  dis as text "Caution: the adversarial pre-trend test was skipped for subgroup `sgcase': `sgndrop' relative period effect(s) were dropped as not estimable, so its coefficients no longer span a complete relative-time grid and a line fit through them would not be interpretable."
                }
                else if  `finaladvtest' > 2  {
                  wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn_s`sgcase'') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') sg(`sgcase')
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
      if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'"=="" {

        foreach sgcase in `sglist' {
          if "`att'" == "att"   {
             qui margins  [`weight'=`w_att'] ,`semielasticity'(`targetflag') noestimcheck  subpop(if `ifstatement'  & `subgroups' == `sgcase' & `t' < `ttre' & `ttre' != .   &  (`_t' == 0 |`_t' >= `horizonstarter')    ) contrast(overjoint) over(`_t') `poismarginstext' `marginsdof'
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
            qui margins  [`weight'=`w_att_it'] ,`semielasticity'(`targetflag') noestimcheck  subpop(if `ifstatement' & `subgroups' == `sgcase' & `t' < `ttre' & `ttre' != .   &  (`_t' == 0 |`_t' >= `horizonstarter')   ) contrast(overjoint) over(`_t') `poismarginstext' `marginsdof'
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
        qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof' `sgpost'
        tempname pwsgs_`item'
        matrix `f9tab' = r(table_vs)
        if rowsof(`f9tab') < 6 {
          matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
          matrix rownames `f9tab' = b se t pvalue ll ul
        }
        matrix `pwsgs_`item''  =  `f9tab'[1..6,floor(`=colsof(r(table_vs))'/2)+1..`=colsof(r(table_vs))']'
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
          qui est restore ___wooldidfullmodel
          if "`sgfallback'" != "" {
            qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof'
            local joint_`item'_sdiff = r(p)[1,2]
          }
        }
        else local joint_`item'_sdiff

        * repeat for pre-treatment
        if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
          cap qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof' `sgpost'
          if _rc != 0 {
            * the pre-ame only exists for 1 sg (or margins otherwise refused): no joint
            * test and no pairwise rows, exactly as under the previous contrast-first order
            local prejoint_`item'_sdiff
            if "`sgpost'" != "" qui est restore ___wooldidfullmodel
          }
          else {
            tempname prepwsgs_`item'
            matrix `f9tab' = r(table_vs)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `prepwsgs_`item''  =  `f9tab'[1..6,floor(`=colsof(r(table_vs))'/2)+1..`=colsof(r(table_vs))']'
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
              qui est restore ___wooldidfullmodel
              if "`sgfallback'" != "" {
                qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof'
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
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre'  ) noestimcheck over(`subgroups') `poismarginstext' `marginsdof'
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
            qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0  ) noestimcheck over(`subgroups') `poismarginstext' `marginsdof'
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
                qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof'
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
                if "`adversarial'" == "adversarial" & `amescrub' == 0  {
                  if  `finaladvtest' > 1  {
                    wooldid_advstat , b(`f2bmat') v(`f2vmat') relyr(`relyrcolumn_s`sgcase'') mint(`min_t') esrel(`esrelativeto') dfr(`df_r') esf(`esfixedbaseperiod') tt(`timetrends') alt(`altmethod_adversarial') sg(`sgcase')
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
        if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'`bootstrap'"=="" & `amescrub' == 0 {
          foreach sgcase in `sglist' {
            if "`att'" == "att" {
              qui margins  [`weight'=`w_att'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & `subgroups' == `sgcase'   &  (`_t' == 0 |`_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' `marginsdof'
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
              qui margins [`weight'=`w_att_it'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & `subgroups' == `sgcase'  &  (`_t' == 0 |`_t' >= `horizonstarter')    ) contrast(overjoint) `poismarginstext' `marginsdof'
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

          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstext' noestimcheck over(`subgroups') pwcompare(effects) `marginsdof' `sgpost'
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
            qui est restore ___wooldidfullmodel
            if "`sgfallback'" != "" {
              qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstext' noestimcheck over(`subgroups') contrast(overjoint) `marginsdof'
              local joint_C`item'_sdiff = r(p)[1,1]
            }
          }
          else local joint_C`item'_sdiff

          * repeat for pre-treatment
          if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
            cap qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof' `sgpost'
            if _rc != 0 {
              local prejoint_C`item'_sdiff
              if "`sgpost'" != "" qui est restore ___wooldidfullmodel
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
                qui est restore ___wooldidfullmodel
                if "`sgfallback'" != "" {
                  qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof'
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
  if "`poisson'" == "" matrix colnames `wdidpooled' =  estimate se t p lb_95ci ub_95ci   relyr
  if "`poisson'" == "poisson" matrix colnames `wdidpooled' =  estimate se z p lb_95ci ub_95ci   relyr


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
          if "`att'" == "att" {
            matrix `wdid_sg`sgcase'' = `wdid_sg`sgcase'' \ `___att_horiz_s`sgcase''
            matrix drop `___att_horiz_s`sgcase''
          }
          if "`att_i'" == "att_i"  |  "`att_it'" == "att_it" |  "`att_itime'" == "att_itime" {
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
            if strpos("`relyrlist_s`sgcase'' ","`testtext' ") != 0 local rownamessg `rownamessg' `rownametest'
            * note: the spaces above are required to avoid triggering "1 is present" when presented with "11"
          }
        }
        matrix rownames `wdid_sg`sgcase'' = `rownamessg'
        local rownamessg`sgcase' `rownamessg'
        if "`poisson'" == "" matrix colnames `wdid_sg`sgcase'' =  estimate se t p lb_95ci ub_95ci   relyr
        if "`poisson'" == "poisson" matrix colnames `wdid_sg`sgcase'' =  estimate se z p lb_95ci ub_95ci   relyr


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
      if "`poisson'" == "" matrix colnames `pwsgs' =  estimate se t p lb_95ci ub_95ci
      if "`poisson'" == "poisson" matrix colnames `pwsgs' =  estimate se z p lb_95ci ub_95ci

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
  foreach matrix in `cfxmatricestoreturn' {
    wooldid_cfxstripe ``matrix''
    ereturn matrix `matrix' = ``matrix''
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
  ereturn hidden local wdid_replayv "1"


  ********* begin actual result presentation *************
  * F11: the presentation code lives in wooldid_display (non-eclass, reads e() only).
  * it is split in two parts so the eclass wooldid_ptp call below stays at worker level
  * (its ereturns would be discarded inside a non-eclass subprogram)
  * F-boot deferred display: under bootstrap the shell rewires e() to the bootstrap
  * covariance first and then calls the (e()-driven) display itself - printing here
  * would show analytic inference that is about to be replaced. pretrendpower is
  * refused with bootstrap at parse (plan A13), so the ptp call needs no gate
  if "`bootstrap'" == "" wooldid_display , part(1)
  * F7: Roth (2022) pre-trend power diagnostics (compute + display + ereturns live in a
  * subprogram; cap noi so a failure here cannot abort an otherwise-complete run)
  if "`pretrendpower'" != "" cap noi wooldid_ptp , levels(`pretrendpower') esrel(`esrelativeto') dfr(`df_r') pois(`poisson') poisexp(`poisexpresults') tvar(`t') ttrevar(`ttre') ivar(`i') wvar(`w') cwm(`customweightmultiplier')
  if "`bootstrap'" == "" wooldid_display , part(2)




  ***************************************************************************************************************************
  ********************************************         PRODUCE EVENT STUDY PLOTS       **************************************
  ***************************************************************************************************************************
  * present event study plots if requested; a little complicated because of naming practices in wdidpooled
  * this is also somewhat awkward placement b/c I have this coded to use the wdidpooled matrix
  * this choice is vestigial to a prior program design and is kept here on the principle of not redoing what already works
  * though I should probably redo this someday

  if `espre' +`espl' > 0 & "`makeplots'" == "makeplots" {
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

  if "`verbose'" != "verbose" qui est drop ___wooldidfullmodel
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
syntax [varlist(numeric default=none)] [if/]  [aw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist numeric) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose GTEffects PERMute(string) PERMREPS(integer 0) PERMEXact PERMTStat PERMTci PERMSTRata(varname numeric) PERMPRobs(string asis) BOOTstrap(string) reps(integer 0) ties Level(real 0) dots(integer -1) rseed(string) BOOTALPha(string) BOOTDf(string)]

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
    if "`permstray'" != "" {
      dis as error "Error: `permstray' require(s) the permute() option."
      error 198
    }
  }
  if "`permute'" != "" {
    cap noi wooldid_permctl , permute(`permute') permreps(`permreps') `permexact' clus(`cluster') ttrevar(`ttre_asgiven') permstrata(`permstrata') permprobs(`permprobs') impv(`impvalues') uncond(`unconditionalse') cw(`contreatwithin') os(`oldsyntax') ts(`permtstat') tc(`permtci') fn(`fastnose') cc(`contreatcontrols') ifr(`if')
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
  cap noi wooldid_bootctl , bootstrap(`bootstrap') reps(`reps') `ties' level(`level') dots(`dots') rseed(`rseed') bootalpha(`bootalpha') bootdf(`bootdf') permute(`permute') impv(`impvalues') uncond(`unconditionalse') cvce(`customvce') fn(`fastnose') os(`oldsyntax') ssg(`suppresssgprimaryeffects') attit(`att_it') atti(`att_i') attitime(`att_itime') si(`safetyinput') clus(`cluster') ivar(`i') ifr(`if') jt(`jointtests') ptp(`pretrendpower')
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
  foreach wdv_ret in att att_it att_itime att_i cluster customvce unconditionalse robust impvalues fe coarsecohortcontrols ccc_absorb controls interactivecontrols timetrends subgroups customweightmultiplier esprelength espostlength saveplots poisson poisexpresults oostreatedcontrols anticipation esfixedbaseperiod esrelativeto safety makeplots makecfxplots summarystats histogramcohorteffects regtol jointtests suppresssgprimaryeffects contreat contreatpoly lattice latticeignoreweights ccnoabsorb contreatcontrols verbose cfxplottypes cfxlattice makecfxplotsbysg cfxplotnose semielasticity contreatelasticitytype contreatwithin cleanmatrices emptycellsoverride replace altmethod_adversarial adversarial update oldsyntax pretrendpower fastnose gteffects permute permreps permexact permtstat permtci permstrata permprobs y i t ttre_asgiven safetyinput suppressjointtests latticenum latticetype esposttrip cfxpost cfxpre cfxrelyrlist bootstrap reps ties level dots rseed bootalpha bootdf esprelengthraw {
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



*** F6: parse-time screens for permutation inference (author's brief: cluster nesting
*** in treatment-onset levels; refuse impvalues/unconditionalse/contreatwithin; reps
*** xor exact; strata subset checks). kept out of wooldid_validate for size headroom.
cap prog drop wooldid_permctl
prog def wooldid_permctl
version 16.0
  syntax , [permute(string) permreps(integer 0) PERMEXact clus(varlist) ttrevar(varname) permstrata(varname numeric) permprobs(string asis) impv(string) uncond(string) cw(string) os(string) ts(string) tc(string) fn(string) cc(string asis) ifr(string asis)]

  local permmethod `permute'
  if "`permmethod'" != "shuffle" & "`permmethod'" != "probabilities" {
    dis as error "Error: permute() must be permute(shuffle) or permute(probabilities). Strata and assignment probabilities have their own options: permstrata() and permprobs()."
    error 198
  }
  if "`impv'" != "" {
    dis as error "Error: permute() is not compatible with impvalues."
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
  syntax , [bootstrap(string) reps(integer 0) ties level(real 0) dots(integer -1) rseed(string) bootalpha(string) bootdf(string) permute(string) impv(string) uncond(string) cvce(string asis) fn(string) os(string) ssg(string) attit(string) atti(string) attitime(string) si(string) clus(varlist) ivar(varname) ifr(string asis) jt(string) ptp(string asis)]

  * stray-companion screen: every companion requires the master option
  if "`bootstrap'" == "" {
    local bootstray
    if `reps' != 0 local bootstray reps()
    if "`ties'" != "" local bootstray `bootstray' ties
    if `level' != 0 local bootstray `bootstray' level()
    if `dots' != -1 local bootstray `bootstray' dots()
    if "`rseed'" != "" local bootstray `bootstray' rseed()
    if "`bootalpha'" != "" local bootstray `bootstray' bootalpha()
    if "`bootdf'" != "" local bootstray `bootstray' bootdf()
    if "`bootstray'" != "" {
      dis as error "Error: `bootstray' require(s) the bootstrap() option."
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
    dis as error "Error: level() must be between 10 and 99.99 for bootstrap confidence intervals."
    error 198
  }
  if `dots' != -1 & `dots' < 0 {
    dis as error "Error: dots() must be a nonnegative integer (0 suppresses the progress display)."
    error 198
  }
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
  if "`bootstrap'" == "standard" & "`cellwtd'" != "" {
    * the resampling copies of a drawn cluster must be re-identifiable as distinct
    * units, which is defined only when every i lies wholly inside one resampling
    * cluster. observation-level resampling (no cluster()) can never satisfy this.
    if "`clus'" == "" {
      dis as error "Error: bootstrap(standard) with `cellwtd' requires cluster(): observation-level resampling cannot re-identify duplicated treatment cells. Specify cluster() at or above the i level, or request att/ame estimands only."
      error 198
    }
    tempvar smp cmin cmax
    if `"`ifr'"' != "" qui gen byte `smp' = ((`ifr')==1) & !missing(`clus')
    else qui gen byte `smp' = !missing(`clus')
    qui gegen double `cmin' = min(`clus') if `smp', by(`ivar')
    qui gegen double `cmax' = max(`clus') if `smp', by(`ivar')
    cap assert `cmin' == `cmax' if `smp'
    if _rc != 0 {
      dis as error "Error: bootstrap(standard) with `cellwtd' requires the cluster() variable to nest i (every i wholly inside one cluster): resampled clusters are duplicated as whole units, which is undefined when a unit straddles clusters. Cluster at or above the i level, or request att/ame estimands only."
      error 459
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
  if "`ptp'" != "" {
    dis as error "Error: bootstrap() may not be combined with pretrendpower(); the Roth diagnostics are computed from the analytic event-study covariance before the bootstrap rewiring. Feed e(b_es_*)/e(V_es_*) from a bootstrap run to the pretrends package directly if bootstrap-based power diagnostics are needed."
    error 198
  }
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
  syntax , [level(real 0) acc(string) obsb(string) comp(string) jt(string) ties done(string) bootdf(string) nclus(integer 0)]
  if `level' == 0 local level 95
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

  * -------- joint pre-trend tests (only when the user requested jointtests)
  if "`jt'" != "" {
    foreach fam in att att_it ame ame_it {
      * overall
      local cols
      forval k = 1/`K' {
        if "`wdbeq`k''" == "ES_`fam'" & "`wdbsg`k''" == "." & `wdbrl`k'' < 0 & `wdbrl`k'' != . local cols `cols' `k'
      }
      * one pre horizon is still a joint test (of a single restriction) - the analytic
      * path posts one, so the bootstrap must too (F8)
      if `: word count `cols'' >= 1 {
        local wdbjp .
        mata: wooldid_bootjw("`acc'", "`obsb'", "`cols'", "", `wdbdf')
        local ename = cond(strpos("`fam'", "ame"), "joint_Choriz", "joint_horiz")
        if "`fam'" == "att_it" | "`fam'" == "ame_it" local ename `ename'_i
        if `wdbjp' != . ereturn scalar `ename' = `wdbjp'
        else {
          local ++wdbncaut
          local wdbc`wdbncaut' the bootstrapped `ename' joint test is unavailable (rank-deficient or thin replication covariance)
        }
      }
      * per subgroup (ereturn LOCALS - the analytic macro convention)
      foreach sg of local sgl {
        local cols
        forval k = 1/`K' {
          if "`wdbeq`k''" == "ES_`fam'" & "`wdbsg`k''" == "`sg'" & `wdbrl`k'' < 0 & `wdbrl`k'' != . local cols `cols' `k'
        }
        if `: word count `cols'' >= 1 {
          local wdbjp .
          mata: wooldid_bootjw("`acc'", "`obsb'", "`cols'", "", `wdbdf')
          local ename = cond(strpos("`fam'", "ame"), "joint_Choriz", "joint_horiz")
          if "`fam'" == "att_it" | "`fam'" == "ame_it" local ename `ename'_i
          if `wdbjp' != . ereturn local `ename'_s`sg' = `wdbjp'
          else {
            * F7: say so - an absent per-subgroup joint test used to be indistinguishable
            * from one that was never requested
            local ++wdbncaut
            local wdbc`wdbncaut' the bootstrapped `ename' joint test for subgroup `sg' is unavailable (rank-deficient or thin replication covariance)
          }
        }
      }
    }
  }

  * -------- allsgssame equality chains (mirror the analytic surface: computed
  * whenever the analytic run posted the scalar, i.e. >1 realized sg; blanked-then-
  * filled is now compute-directly)
  if `: word count `sgl'' > 1 {
    foreach it in att att_it att_i att_itime ame ame_it {
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
        if `ok' & "`e(joint_`pre'`it'_allsgssame)'" != "" {
          local wdbjp .
          mata: wooldid_bootjw("`acc'", "`obsb'", "`cols'", "chain", `wdbdf')
          if `wdbjp' != . ereturn scalar joint_`pre'`it'_allsgssame = `wdbjp'
          else {
            ereturn scalar joint_`pre'`it'_allsgssame = .
            local ++wdbncaut
            local wdbc`wdbncaut' the bootstrapped `pre'`it' subgroup-equality test is unavailable (rank-deficient or thin replication covariance)
          }
        }
      }
    }
  }

  * -------- pairwise subgroup contrasts: keep the analytic point differences,
  * rewire se/z/p/ci from the bootstrap covariance of the matching estimand columns
  cap confirm matrix e(pwsgs)
  if _rc == 0 {
    tempname pw
    matrix `pw' = e(pwsgs)
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
          local wdbse .
          mata: wooldid_bootpw("`acc'", `ka', `kb')
          if `wdbse' < . & `wdbse' > 0 {
            matrix `pw'[`rr', 2] = `wdbse'
            matrix `pw'[`rr', 3] = `pw'[`rr', 1]/`wdbse'
            matrix `pw'[`rr', 4] = cond(`wdbdf' < ., 2*ttail(`wdbdf', abs(`pw'[`rr', 3])), 2*normal(-abs(`pw'[`rr', 3])))
            local zc = cond(`wdbdf' < ., invt(`wdbdf', 1 - (100 - `level')/200), invnormal(1 - (100 - `level')/200))
            matrix `pw'[`rr', 5] = `pw'[`rr', 1] - `zc'*`wdbse'
            matrix `pw'[`rr', 6] = `pw'[`rr', 1] + `zc'*`wdbse'
          }
          else {
            forval cc = 2/6 {
              matrix `pw'[`rr', `cc'] = .
            }
          }
        }
        else {
          * unmapped: blank rather than leave the analytic se/t/p/CI in a table the
          * display attributes to the bootstrap covariance (the commit-C blanking
          * invariant, restored here for every future label surprise as well)
          forval cc = 2/6 {
            matrix `pw'[`rr', `cc'] = .
          }
          local ++pwunmapped
        }
      }
      else {
        forval cc = 2/6 {
          matrix `pw'[`rr', `cc'] = .
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
    matrix colnames `pw' = `pwc'
    ereturn matrix pwsgs = `pw'
    if `pwunmapped' > 0 {
      local ++wdbncaut
      local wdbc`wdbncaut' `pwunmapped' pairwise subgroup contrast row(s) could not be mapped to replication columns; their bootstrap inference is left missing
    }
  }

  * -------- CI flavors: normal / percentile / bias-corrected at level(), plus the
  * BC median-bias z0 (comparison against the REP-0 fastnose comparator per
  * amendment A3; `ties' adjusts z0 for replicate-comparator ties, the built-in
  * bootstrap's ties semantics)
  tempname cin cip cib cz0
  mata: wooldid_bootci("`acc'", "`obsb'", "`comp'", `level', ("`ties'" != ""), "`cin'", "`cip'", "`cib'", "`cz0'", `wdbdf')
  local obn : colfullnames `obsb'
  foreach cm in cin cip cib {
    matrix colnames ``cm'' = `obn'
    matrix rownames ``cm'' = ll ul
  }
  matrix colnames `cz0' = `obn'
  ereturn matrix boot_ci_normal = `cin'
  ereturn matrix boot_ci_percentile = `cip'
  ereturn matrix boot_ci_bc = `cib'
  ereturn matrix boot_z0 = `cz0'
  if "`ties'" != "" ereturn local boot_ties ties

  forval cc = 1/`wdbncaut' {
    dis as text "Caution: `wdbc`cc''."
  }
end

cap prog drop wooldid_permrun
prog def wooldid_permrun
version 16.0
  syntax , y(varname) i(varname) t(varname) ttrevar(varname) clus(varlist) permute(string) estimands(string) distname(string) [permreps(integer 0) PERMEXact permstrata(varname numeric) permprobs(string asis) PERMTStat repopts(string asis) wtype(string) wexp(string) ifr(string asis) compname(string)]

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
  }

  * assignments + counts (mata): fills the C x R matrix, sets wdp_R, the possible-
  * assignment count (double + string), and caution flags
  tempname wdpassign wdpP wdpL
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
  * observed-is-a-member property (measured: b >= local-copy-of-b evaluates FALSE)
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
  forval r = 1/`wdp_R' {
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
    if `badrep' == 0 local ++ndone
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
  }
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
  if `nfail' > 0 & "`permexact'" != "" dis as text "Caution: `nfail' of `wdp_R' enumerated assignments failed to estimate, so the exact enumeration is incomplete; the reported exact p-value is a share over completed assignments only (the observed-is-a-member and 2/M-floor readings no longer apply)."
  else if `nfail' > 0 & `nfail' > 0.1*`wdp_R' dis as text "Caution: `nfail' of `wdp_R' permutation replications failed to estimate and were excluded."

  * summaries per estimand: add-one two-sided p (exact: plain share - the observed
  * assignment is a member of the enumeration), 5th/95th percentile band, SD
  tempname rsfr
  matrix `distname' = J(`ndone', `: word count `estimands'', .)
  local ec 0
  foreach ee in `estimands' {
    local ++ec
    frame create `rsfr' double s
    frame `rsfr' {
      qui set obs `wdp_R'
      qui replace s = el("`res'", _n, `ec')
      qui drop if missing(s)
      qui count if abs(s) >= abs(`wdpo`ec'')
      local nex = r(N)
      qui count
      local nuse = r(N)
      _pctile s, p(5 95)
      local lb = r(r1)
      local ub = r(r2)
      qui sum s
      * NFPERM-11: an exact enumeration IS the complete population of assignments, so
      * report the population SD (denominator M); sampled draws keep the n-1 sample SD.
      * r(Var) is full double precision - no display-rounding hazard
      if "`permexact'" != "" local sdv = sqrt(r(Var)*(r(N)-1)/r(N))
      else local sdv = r(sd)
      if `sdv' == 0 dis as text "Caution: the permuted `ee' distribution is degenerate (every completed draw equals every other); the permutation p-value is uninformative."
      local q95 .
      if "`permtstat'" != "" {
        qui gen double abss = abs(s)
        _pctile abss, p(95)
        local q95 = r(r1)
        qui drop abss
      }
    }
    frame drop `rsfr'
    if "`permexact'" != "" local p_`ee' = `nex'/`nuse'
    else local p_`ee' = (1 + `nex')/(`nuse' + 1)
    c_local wdp_r_p_`ee' `p_`ee''
    c_local wdp_r_lb_`ee' `lb'
    c_local wdp_r_ub_`ee' `ub'
    c_local wdp_r_sd_`ee' `sdv'
    c_local wdp_r_q_`ee' `q95'
  }
  * dense permdist matrix (completed rows only), built once from res
  local rr 0
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
  c_local wdp_r_npd `wdp_npd'
  c_local wdp_r_nps `wdp_nps'
end




cap prog drop wooldid_bootrun
prog def wooldid_bootrun
version 16.0
*** F-boot engine (plan sections 2-4 + amendments): builds the cluster map frame,
*** generates the Dirichlet draw matrix (mata, rwgen-bitwise arithmetic), runs the
*** rep loop as nested fastnose wooldid calls under reweighted [aw/pw], and
*** accumulates the FULL rep e(b) vectors name-aligned into the caller's matrix.
*** bayes semantics (locked 2026-08-30): cluster c receives total Dirichlet mass
*** draw_c; within the cluster the mass is split in proportion to existing weights
*** (row multiplier = draw_c / W_c, W_c = cluster base-weight total). This matches
*** bayesboot BITWISE on unweighted data (bootalpha(size), the default) and rwgen
*** with priorpowers(-n_c) for bootalpha(equal), and is the Rubin generalization
*** for weighted data ("multiply existing weights by the cluster draw").
  syntax , [mode(string) reid reps(integer 0) obsb(string) accname(string) y(varname) i(varname) t(varname) ttrevar(varname) clus(varlist) dots(integer 1) rseed(string) bootalpha(string) repopts(string asis) wtype(string) wexp(string) ifr(string asis) pois(string) compname(string)]

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

  * draws. bayes: C x R Dirichlet matrix up front - per-rep (outer) per-cluster-
  * ascending (inner) Gamma draws from the seeded mt64 stream, the exact rwgen bayes
  * consumption order; row multiplier = (g/colsum)/W_c. standard: draws happen
  * INSIDE the rep loop (one bsample per rep - the official bootstrap prefix's
  * consumption sequence, so same-seed runs are comparable; nested wooldid calls
  * leave c(rngstate) untouched, measured)
  if "`rseed'" != "" set seed `rseed'
  tempname wdbdrawm
  if "`mode'" == "bayes" mata: wooldid_bootdw("`mbfr'", "`bootalpha'", `reps', "`wdbdrawm'")

  * rep command line
  local wtif
  local wtype2 `wtype'
  if "`wtype2'" == "" local wtype2 = cond("`pois'" != "", "pweight", "aweight")
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

  qui frlink m:1 `clusv', frame(`mbfr') generate(`lnk')
  local ndone 0
  local nfail 0
  local nblank 0
  if `dots' > 0 dis as text "Bootstrap replications (`reps'): " _continue
  tempname wdbcf wdbrf
  forval r = 1/`reps' {
    if "`mode'" == "bayes" {
      frame `mbfr': qui replace __wdbdraw = el("`wdbdrawm'", __wdbrow, `r')
    }
    else if "`reid'" == "" {
      * standard count path: cluster multiplicities from bsample on a COPY of the
      * map frame (map-frame draws == full-data draws under the cluster()+
      * idcluster() form, measured; plain bsample diverges - A17)
      cap frame drop `wdbcf'
      frame copy `mbfr' `wdbcf'
      frame `wdbcf' {
        bsample, cluster(`clusv') idcluster(__wdbns)
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
        qui gegen double __wdbnewi = group(__wdbns `i')
        nobreak {
          global WOOLDID_REPMODE 1
        cap qui wooldid `y' __wdbnewi `t' `ttrevar' `wtif0' , `repopts'
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
        mata: wooldid_bootgrab("`accname'", `r', "`obsb'")
        local nblank = `nblank' + `wdb_gblank'
        local ++ndone
      }
      else local ++nfail
      cap frame drop `wdbrf'
      if `dots' > 0 {
        if `badrep' dis as error "x" _continue
        else dis as text "." _continue
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
      mata: wooldid_bootgrab("`accname'", `r', "`obsb'")
      local nblank = `nblank' + `wdb_gblank'
    }
    if `badrep' == 0 local ++ndone
    else local ++nfail
    if `dots' > 0 {
      if `badrep' dis as error "x" _continue
      else dis as text "." _continue
      if mod(`r', 50) == 0 dis as text " " `r' _newline "  " _continue
    }
  }
  if `dots' > 0 dis as text " done"
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
  syntax , [mode(string) reps(integer 0) acc(string) obsb(string) done(string) fail(string) blank(string) level(real 0) bootalpha(string) clus(varlist) comp(string) rseed(string) bootdf(string) nclus(integer 0)]
  if `level' == 0 local level 95
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
  local sgl `e(wdid_sglist)'
  local kk 0
  forval rr = 1/`nmain' {
    local ++kk
    if `wm'[`rr', 1] >= . matrix `excl'[1, `kk'] = 1
  }
  foreach sg of local sgl {
    tempname wsg`sg'
    cap matrix `wsg`sg'' = e(wdidsgresults_`sg')
    if _rc == 0 {
      forval rr = 1/`=rowsof(`wsg`sg'')' {
        local ++kk
        if `wsg`sg''[`rr', 1] >= . matrix `excl'[1, `kk'] = 1
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

  * -------- rewire the results matrices: cols 2-6 = se t/z p lb ub at level()
  local zc = cond(`wdbdf' < ., invt(`wdbdf', 1 - (100 - `level')/200), invnormal(1 - (100 - `level')/200))
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
  ereturn matrix wdidmainresults = `wm'
  foreach sg of local sgl {
    cap confirm matrix `wsg`sg''
    if _rc != 0 continue
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
  }


  * -------- header + metadata
  local wdbmodetext = cond("`mode'" == "bayes", "Bayesian cluster bootstrap", "cluster bootstrap")
  local wdbcl `clus'
  if "`wdbcl'" == "" {
    local wdbmodetext = cond("`mode'" == "bayes", "Bayesian bootstrap", "bootstrap")
    local wdbcl observation-level
  }
  else local wdbcl cluster(`wdbcl')
  ereturn hidden local wdid_vcespec "`wdbmodetext', `wdbcl', reps(`reps')"
  ereturn local boot_method `mode'
  ereturn scalar boot_reps = `reps'
  ereturn scalar boot_completed = `done'
  ereturn scalar boot_failed = `fail'
  ereturn scalar boot_level = `level'
  ereturn local boot_stat `wdblab'
  if `wdbdf' < . {
    ereturn scalar boot_df = `wdbdf'
    * align e(df_r) with the reference distribution actually used, so that test and
    * lincom run by hand on the bootstrap covariance reproduce the reported p-values
    ereturn scalar df_r = `wdbdf'
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
    mata: wooldid_bootdist("`acc'", "`bdist'")
    cap matrix colnames `bdist' = `obn'
    cap ereturn matrix bootdist = `bdist'
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
  syntax , levels(numlist) [esrel(string) dfr(string) pois(string) poisexp(string) tvar(string) ttrevar(string) ivar(string) wvar(string) cwm(string)]

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
  if `ptpdfok' local ptpcrit = invttail(`dfr', .025)
  else local ptpcrit = invnormal(.975)

  * hbar: estimand-weighted mean post-treatment relative period over the FULL post
  * support of the estimation sample - the pooled estimands aggregate over all post
  * cells (subpop t >= onset), not only the reported event-study window, and using the
  * reported window instead understates the implied bias whenever espostlength is
  * smaller than the largest post relative period
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

  local ptpheader 0
  foreach ptpty in att att_it ame ame_it {
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
    dis as result "    slope = " %-9.0g `ptpslope' " (se " %-9.0g `ptpse' "; 95% CI " %-9.0g `ptpll' " to " %-9.0g `ptpul' ")"
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

cap prog drop wooldid_display
prog def wooldid_display
version 16.0
  syntax , part(string)
  * replay contract gate: e(wdid_replayv) is the version of the display contract the
  * store was written under. bump the posted value (in wooldid_worker) whenever the
  * display contract changes; the exact match here refuses cross-version replays
  * instead of misrendering a store this code no longer understands
  if "`e(wdid_replayv)'" != "1" {
    dis as error "Error: the stored results were posted by a wooldid version without replay support or with a different replay contract (e(wdid_replayv)='`e(wdid_replayv)'', this version requires '1'); nothing can be redisplayed. Refit wooldid."
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
       dis as result "p-value for joint test of null that all (non-reference) `typetext' coefficients are on a single line: " %-6.4f e(adversarial_pval_`typetext')
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
    dis as result "    slope = " %-9.0g `ptpslope' " (se " %-9.0g `ptpse' "; 95% CI " %-9.0g `ptpll' " to " %-9.0g `ptpul' ")"
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
              dis as result "p-value for joint test of null that all (non-reference) `typetext' coefficients for sg=`sgcase' are on a single line: " %-6.4f e(adversarial_pval_`typetext'_s`sgcase')
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
    dis as text "Note: standard errors are the stage-1 regression's coefficient-level VCE; rows outside the reported event-study window have inwindow = 0."
    if `gtrows' <= 200 {
      tempname gtds
      matrix `gtds' = `gtd'[1..., 1..9]
      matlist `gtds', format(%9.0g)
    }
    else dis as text "  (`gtrows' rows; table not printed - see e(gteffects))"
  }
  }
  * F6: permutation-inference block - existence-gated on e(perm_method), so replay
  * reproduces it and older stores print nothing
  if "`part'" == "perm" {
    if "`e(perm_method)'" != "" {
      dis as result "Randomization Inference (permute(`e(perm_method)')): `e(perm_completed)' of `e(perm_reps)' replications completed; `e(perm_npossible_str)' possible assignments."
      dis as text "Two-sided permutation p-values; the interval is the 5th-95th percentile band of the permuted-effect distribution (a null band, not a parameter CI); SD is the permuted-effect standard deviation."
      if "`e(perm_stat)'" == "t" dis as text "Statistics are randomized on t-ratios (permtstat)."
      local wdptcihedge 0
      foreach wdpe in `e(wdid_permest)' {
        local wdpsc = cond("`e(perm_stat)'" == "t", " (t-scale)", "")
        dis as result "  `wdpe': estimate " %10.0g e(b)[1, colnumb(e(b), "main_`wdpe'")] "   p " %6.4f e(perm_p_`wdpe') "   null band`wdpsc' [" %10.0g e(perm_lb_`wdpe') ", " %10.0g e(perm_ub_`wdpe') "]   SD`wdpsc' " %10.0g e(perm_sd_`wdpe')
        if e(perm_tlb_`wdpe') < . {
          dis as result "        inverted-t 95% CI [" %10.0g e(perm_tlb_`wdpe') ", " %10.0g e(perm_tub_`wdpe') "]"
          local wdptcihedge 1
        }
      }
      * NFPERM-08 (honesty hedge, once per run): the inverted-t CI is an approximation,
      * not a sharp-null randomization-test inversion
      if `wdptcihedge' dis as text "The inverted-t CI is a fixed-zero-null pivotality approximation (it relies on the permuted t-ratio distribution approximating the t-ratio's distribution under shifted nulls); it is not a sharp-null randomization-test inversion."
      dis as result " "
    }
  }


 if "`part'" == "boot" {
   * F-boot summary block (existence-gated: replay of non-bootstrap stores prints
   * nothing). all content from e() only
   if "`e(boot_method)'" == "" exit
   local wdbmt = cond("`e(boot_method)'" == "bayes", "Bayesian cluster bootstrap (Dirichlet reweighting, alpha = `e(boot_alpha)')", "cluster bootstrap (resampling with replacement)")
   dis as result "Bootstrap Inference (`wdbmt'):"
   dis as result "  " e(boot_completed) " of " e(boot_reps) " replications completed" _continue
   if e(boot_failed) > 0 dis as result " (" e(boot_failed) " failed and excluded)" _continue
   local wdbrefd = cond("`e(boot_stat)'" == "t", "t(" + strofreal(e(boot_df)) + ")", "standard normal")
   dis as result "; reported standard errors, `e(boot_stat)' statistics, p-values and " e(boot_level) "% CIs are bootstrap quantities referred to the `wdbrefd' distribution; percentile and bias-corrected CIs are stored in e(boot_ci_percentile)/e(boot_ci_bc)."
   if e(boot_v_posted) == 1 dis as result "  e(V) holds the bootstrap covariance (complete-replication listwise, e(boot_N_V) replications; the analytic covariance, when one existed, remains in the hidden e(V_analytic))."
   else dis as result "  The bootstrap covariance is stored in e(V_boot); no covariance is posted to e(V) for this run."
   dis as result "  Joint tests and subgroup contrasts, when requested, are " cond("`e(boot_stat)'" == "t", "F", "chi-squared") " Wald tests on the bootstrap covariance."
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
  syntax , b(string) v(string) relyr(string) mint(string) esrel(string) [dfr(string) esf(string) tt(string) alt(string) sg(string)]

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
  * the little reg clobbered e(); put the main model back (stored at worker start)
  qui est restore ___wooldidfullmodel
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
cap prog drop wooldid_esplots
prog def wooldid_esplots
version 16.0
  syntax [, SAVEplots(name) replace]
  * guard the hand-call surface: without a posted wooldid result, the body would die
  * later inside its own tempname frame with a misleading "relyr not found" (an absent
  * e() matrix silently reads as a 1x1 missing scalar), stranding a caller who typed
  * this directly from a non-default frame in the default frame
  * exact-match on the replay contract version, mirroring the wooldid_display gate
  if "`e(wdid_replayv)'" != "1" {
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


    if "`saveplots'" != "" qui save "`saveplots'_ES.dta", `replace'
    clear

    * repeat by SG if appropriate
    if "`subgroups'" != "" & "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        tempname wdidsgmat
        matrix `wdidsgmat' = e(wdidsgresults_`sgcase')
        qui svmat double `wdidsgmat', names(col)
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

        if "`saveplots'" != ""  qui save "`saveplots'_ES_`subgroups'_`sgcase'.dta", `replace'
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
                 string scalar outname)
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
                crit = invttail(df, .025)
            }
            else {
                p = 2 * normal(-abs(z))
                crit = invnormal(.975)
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
    st_matrixcolstripe(outname, (J(15, 1, ""), ("i", "onset", "caltime", "relyr", "estimate", "se", "p", "lb_95ci", "ub_95ci", "inwindow", "nobs", "w_att", "w_att_it", "w_att_i", "w_att_itime")'))
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
mata:
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
    if (alt == 0) {
        X = J(rows(rel), 1, 1), rel
        V = (V + V')/2
        U = .
        lam = .
        symeigensystem(V, U, lam)
        scale = max(abs(lam))
        if (scale <= 0 | scale >= .) return
        tol = scale*1e-12
        sel = selectindex(lam' :> tol)
        k1 = rows(sel)
        if (k1 == 0) return
        T = diag(1:/sqrt(lam[sel]')) * U[., sel]'
        z = T*b
        M = T*X
        sv = svdsv(M)
        rankM = sum(sv :> max(sv)*1e-12)
        q = k1 - rankM
        if (q <= 0) return
        g = svsolve(M, z)
        r = z - M*g
        W = r'r
    }
    else {
        idxP = selectindex(rel :<= -1)
        idxQ = selectindex(rel :>= 0)
        npr = rows(idxP)
        nq = rows(idxQ)
        if (npr == 0 | nq == 0) return
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
        if (ZZi[1,1] == 0 | ZZi[2,2] == 0) return
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
        if (scale <= 0 | scale >= .) return
        tol = scale*1e-12
        sel = selectindex(lamS' :> tol)
        q = rows(sel)
        if (q == 0) return
        c = US[., sel]'d
        resid = d - US[., sel]*c
        if (norm(resid) > 1e-8*max((norm(d), sqrt(scale)))) return
        W = sum(c:^2 :/ lamS[sel]')
    }
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
cap mata: mata drop wooldid_bootgrab()
cap mata: mata drop wooldid_bootfin()
cap mata: mata drop wooldid_bootves()
cap mata: mata drop wooldid_bootdist()
cap mata: mata drop wooldid_bootjw()
cap mata: mata drop wooldid_bootpw()
cap mata: mata drop wooldid_bootci()
cap mata: mata drop wooldid_bootveskey()
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
  P = P :/ W
  st_matrix(outname, P)
}

// name-aligned grab of the current rep e(b) into row r of the accumulation matrix;
// ame-family structural zeros (rep value exactly 0 where the observed b is nonzero)
// are blanked and counted (A12b)
void wooldid_bootgrab(string scalar accname, real scalar r, string scalar obsbname)
{
  real matrix A, rb, ob
  string matrix an, rn
  real scalar k, j, nb
  transmorphic t
  A = st_matrix(accname)
  an = st_matrixcolstripe(accname)
  rb = st_matrix("e(b)")
  rn = st_matrixcolstripe("e(b)")
  ob = st_matrix(obsbname)
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
    }
  }
  st_matrix(accname, A)
  st_matrixcolstripe(accname, an)
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
void wooldid_bootdist(string scalar accname, string scalar outname)
{
  real matrix A
  real colvector ok
  A = st_matrix(accname)
  ok = rowsum(A :< .) :> 0
  st_matrix(outname, select(A, ok))
}


// Wald test of the OBSERVED coefficients in the given columns against zero (or,
// mode=="chain", of their successive differences - the subgroup-equality null),
// using the LISTWISE complete-case bootstrap covariance of those columns. Posts
// st_local wdbjp = p (chi2, df = #constraints) or "." on any guard failure:
// too few listwise rows, non-PSD covariance, or invsym dropping a constraint
// (a dropped constraint must post MISSING, never a residual-rank p).
void wooldid_bootjw(string scalar accname, string scalar obsbname, string scalar colslist, string scalar mode, real scalar dfden)
{
  real matrix A, S, Si, CC, dev, Cm, R
  real rowvector idx
  real colvector ok, b, d, ev, sc
  real scalar k, j, W, df
  idx = strtoreal(tokens(colslist))
  A = st_matrix(accname)[., idx']
  b = st_matrix(obsbname)[1, idx']'
  ok = rowmissing(A) :== 0
  k = cols(A)
  if (mode == "chain") {
    Cm = J(k-1, k, 0)
    for (j = 1; j <= k-1; j++) {
      Cm[j, j] = 1
      Cm[j, j+1] = -1
    }
    d = Cm * b
    df = k - 1
  }
  else {
    Cm = I(k)
    d = b
    df = k
  }
  if (sum(ok) < df + 2 | missing(d)) {
    st_local("wdbjp", ".")
    return
  }
  CC = select(A, ok)
  dev = CC :- mean(CC)
  S = Cm * ((dev' * dev) / (rows(CC) - 1)) * Cm'
  // scale-free guard: invsym's pivot tolerance is ABSOLUTE, so on its own it both
  // misses relative rank deficiency at ordinary scale and refuses well-conditioned
  // covariances at tiny scale. Correlation-normalise, judge the eigenvalues of the
  // normalised matrix, then undo the scaling on the inverse. The 1e-8 relative floor
  // (condition number 1e8) was calibrated on measured conditioning: real event-study
  // pre-blocks sit around 1.5e-2, a near-collinear construction that invsym alone let
  // through at ordinary scale sits at 5.7e-10, and the floor separates them with six
  // orders of margin on the live side
  sc = sqrt(diagonal(S))
  if (min(sc) <= 0 | missing(sc)) {
    st_local("wdbjp", ".")
    return
  }
  R = S :/ (sc * sc')
  ev = symeigenvalues(R)'
  if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) {
    st_local("wdbjp", ".")
    return
  }
  Si = invsym(R) :/ (sc * sc')
  if (min(abs(diagonal(Si))) == 0) {
    st_local("wdbjp", ".")
    return
  }
  W = d' * Si * d
  if (W < 0 | W >= .) {
    st_local("wdbjp", ".")
    return
  }
  // with finite denominator degrees of freedom the joint test is reported as
  // F(df, dfden) = W/df, the small-sample analogue of the chi-squared Wald form
  if (dfden < .) st_local("wdbjp", strofreal(Ftail(df, dfden, W/df), "%21x"))
  else st_local("wdbjp", strofreal(chi2tail(df, W), "%21x"))
}

// bootstrap se of a difference of two columns (listwise over the pair)
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

// CI flavors per column: normal (obs +/- zc*se), percentile (classic Stata
// rule, pinned vs bstat at even AND odd rep counts), bias-corrected (order-statistic
// endpoints v[ceil(n*p_adj)] - the bstat rule, pinned exact) with z0 computed
// against the REP-0 comparator; ties adjusts z0 by half-weighting exact ties.
// Columns with fewer than 2 available reps or zero dispersion post missing.
void wooldid_bootci(string scalar accname, string scalar obsbname, string scalar compname, real scalar lvl, real scalar ties, string scalar cinname, string scalar cipname, string scalar cibname, string scalar z0name, real scalar dfden)
{
  real matrix A, cin, cip, cib
  real rowvector obs, cmp, z0v
  real colvector v, use
  real scalar K, kcol, al, zc, nn, se, z0, x, lo, fr, kk2, pl, pu
  A = st_matrix(accname)
  obs = st_matrix(obsbname)
  cmp = st_matrix(compname)
  K = cols(A)
  al = (100 - lvl)/200
  zc = (dfden < . ? invt(dfden, 1 - al) : invnormal(1 - al))
  cin = J(2, K, .)
  cip = J(2, K, .)
  cib = J(2, K, .)
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
    if (ties) z0 = invnormal((sum(v :< cmp[kcol]) + 0.5*sum(v :== cmp[kcol]))/nn)
    else z0 = invnormal(sum(v :<= cmp[kcol])/nn)
    if (z0 >= .) continue
    z0v[kcol] = z0
    pl = normal(2*z0 + invnormal(al))
    pu = normal(2*z0 + invnormal(1-al))
    kk2 = max((min((ceil(nn*pl), nn)), 1))
    cib[1, kcol] = v[kk2]
    kk2 = max((min((ceil(nn*pu), nn)), 1))
    cib[2, kcol] = v[kk2]
  }
  st_matrix(cinname, cin)
  st_matrix(cipname, cip)
  st_matrix(cibname, cib)
  st_matrix(z0name, z0v)
}

end

