
* Wooldid Version 1.3, 7/31/2024
* Author: Thomas Hegland, Agency for Healthcare Research and Quality
* Contact: thomas.hegland@ahrq.hhs.gov, @thomas_hegland (twitter), thomashegland.com
* This code file is the work of the author. The code comes with no endorsement, guarantee, or warranty
* from the Agency for Healthcare Research and Quality, the Department of Health
* and Human Services, or the United States government. Any views expressed within
* the code file are also those of the author, and not necessarily that of any of
* the aforementioned governmental entities.

cap prog drop wooldid
prog def wooldid, eclass
version 16.0
  * F11: replay - wooldid typed alone after a successful wooldid estimation redisplays
  * the stored results from e() (tables and diagnostics; plots are not redrawn, and no
  * options are accepted on replay). gated on a STRICTLY EMPTY command line, so every
  * other invocation - including "wooldid, update" and an option typo like "wooldid,
  * att" - falls through to the worker and behaves exactly as before
  if `"`0'"' == "" & "`e(cmd)'" == "wooldid" {
    wooldid_display , part(1)
    wooldid_display , part(roth)
    wooldid_display , part(2)
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
    if `wdsh_rc' != 0 {
      * set emptycells is program-scoped (stata restores it automatically at program
      * exit); this restore is belt-and-braces against version drift only
      if "`c(emptycells)'" != "`wdsh_ec'" set emptycells `wdsh_ec'
      * stata auto-drops the worker's tempname frames on error but then leaves the
      * user in the default frame, not the frame wooldid was called from
      qui frame
      if "`r(currentframe)'" != "`wdsh_frame'" cap frame change `wdsh_frame'
      cap est drop __marginsmodel
      if strpos("`wdsh_flags'", "keep") == 0 {
        if strpos("`wdsh_flags'", "fullmodel") != 0 cap est drop ___wooldidfullmodel
        if strpos("`wdsh_flags'", "ppml") != 0 cap drop _ppmlhdfe_d
      }
    }
  }
  if `wdsh_rc' != 0 exit `wdsh_rc'
end

cap prog drop wooldid_worker
prog def wooldid_worker, eclass
version 16.0

* S5: this line must remain BYTE-IDENTICAL to the syntax line in wooldid_validate
syntax [varlist(numeric default=none)] [if/]  [aw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose]

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
  * note: this does not apply to direct pass controls and fe
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
  qui gen `ttre' = `ttre_asgiven' `anticshift' if `ttre_asgiven' `anticshift' <= `maxyearobserved' + `oostreatedcontrols' & `ifstatement'

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
      dis as error "wooldid requires exactly one treatment onset per unit: rows of the same `i' carrying different onsets, or a mix of missing and non-missing onsets, are silently split across cohorts, and part of the unit can be used as its own control."
      dis as error "Common causes: an onset recorded only from the onset period onward (fill it in on every in-sample row of a treated unit, including pre-treatment rows); two different units sharing one `i' value; or unit IDs stored as float, which pools distinct large IDs (recast double `i'). Alternatively, restrict the estimation sample with if()."
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
      qui gen `tp1' = `t' if `t' != `mint'
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
    if "`timetrends'" == "" qui replace `ifvar' = 0 if `ttre' <= `mint'
    if "`timetrends'" == "timetrends" qui replace `ifvar' = 0 if `ttre' <= `mintp1'
    qui replace `ttre' = . if `ifvar' == 0


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
          dis as text "This does not break estimation, but event study estimates can carry counter intuitive meanings and reflect effect levels by cohort."
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
          dis as text "This does not break estimation, but event study estimates can carry counter intuitive meanings and reflect effect levels by cohort."
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
    * by default, estimates are relative to the first time period prior to treatment (-1); we default to this if unacceptable input is given
    * we also modify esprelength to stretch to the relative period
    if `esrelativeto' > -1 local esrelativeto -1
    if `esrelativeto' < -1 & `esprelength' > 0 & `esprelength' < abs(`esrelativeto') {
      local esprelength = abs(`esrelativeto') -1
    }
    if `esprelength' >= abs(`esrelativeto') & "`timetrends'" == "timetrends" local ++esprelength

    local refperiod = `esrelativeto'

    qui gen `_t' = `t' - `ttre' if `ifstatement'


    * here we implement a check to see whether or not the reference period exists for a given cohort
    * if not, we do not use the cohort
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

    * grab information about range of _t
    qui gstats sum `_t' if `ifstatement', meanonly
    local min_t = abs(r(min))
    local max_t = r(max) + 1


    * set _t into all positive terms to allow usge in regressions
    qui replace `_t' = `_t'+ `min_t' + 1

    *redefine the ref period macro for future use and then do a little bit of housekeeping
    * LOAD-BEARING: the four adversarial blocks plant the reference row(s) at t == refperiod
    * via this identity (relyr = t - min_t - 1); if this definition moves, they must move with it
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
      if `failure' == 1 & "`timetrends'" == "timetrends"  dis as error  "Note that inclusion of cohort-specific time trends increases the required number of pre-treatment units."
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

  ****** after the above estimability checks and processing, does the SG variable have surviving members?
  * while we're at it, figure out how many levels the sg variable takes on during the treatment period - will be useful down the line

  if "`subgroups'" != "" {
    qui glevelsof(`subgroups') if `ifstatement' & `t' >= `ttre' , clean local(sglist)
    if wordcount("`sglist'") == 0 {
      dis as error  "Error: subgroups are specified, but contains no values within sample during the treatment period."
      error 459
    }
    * subgroup values reaching the treated estimation sample must be integer-valued: a non-integer
    * value yields an sglist token like 1.5 that later crashes the naming machinery (locals,
    * matrices, ereturns named ..._s1.5) with an unexplained r(198). the check deliberately uses
    * the levelsof domain above, so values appearing only outside it (never-treated units, screened
    * observations) stay tolerated exactly as before, and missing values pass (int(.) == .)
    cap assert `subgroups' == int(`subgroups') if `ifstatement' & `t' >= `ttre'
    if _rc != 0 {
      dis as error "The subgroup variable contains non-integer values within the treated estimation sample; it must contain only non-negative integers."
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
      if "`att_i'" != "" | "`att_itime'" != "" {
        tempvar f13sgimin f13sgimax
        qui gegen double `f13sgimin' = min(cond(`subgroups' >= ., ., `subgroups')) if `ifstatement' & `ttre' != ., by(`i')
        qui gegen double `f13sgimax' = max(cond(`subgroups' >= ., ., `subgroups')) if `ifstatement' & `ttre' != ., by(`i')
        cap assert `f13sgimin' == `f13sgimax' if `ifstatement' & `ttre' != .
        if _rc != 0 {
          dis as text "Caution: the subgroup variable `subgroups' changes over time within some units. For att_i and att_itime, a unit appearing in more than one subgroup enters each subgroup weighted by its share of sample weight there, rather than with the equal unit weighting those estimands normally imply; by-subgroup att and att_it results are unaffected. Consider whether a time-invariant subgroup definition (or a finer unit id) better matches your design."
        }
        qui drop `f13sgimin' `f13sgimax'
      }
    }
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
          qui gen `modified_ttre' = `ttre'
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
  * this does not occur b/c this statement is placed toward the end of the regression statement, meaning multicolinearity forces drops among the controls
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
  if "`timetrends'" == "timetrends"   local trendinjector i.`i'##c.`t'
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
      qui gen `modified_ttre' = `ttre'
      qui replace `modified_ttre' = `minyearobserved' if `ttre' == .
      * minyearobserved is just any value not appearing in ttre; ttre is guaranteed to not include that value
    }

    * initialize
		local contreatcontrolinserterFE
		local contreatcontrolinserter

    * handle basic, byTreated, byTtre, and byCohort in order since these are successively more intensive and supersede eachother

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
  qui est sto ___wooldidfullmodel
  * register this run's scratch objects with the shell's error-path cleanup epilogue
  global WOOLDID_CLEANUP $WOOLDID_CLEANUP fullmodel
  if "`poisson'" == "poisson" & "`unconditionalse'" == "" & "`customvce'" == "" global WOOLDID_CLEANUP $WOOLDID_CLEANUP ppml
  if "`verbose'" == "verbose" dis as text "Wooldid has finished regression estimation and is preparing to begin margins calculations."

	* check on dropped obs
	qui count if `dropcount' & e(sample) == 0
	local droppedobs = r(N)

  * grab some items from the model - sample indicator variable, residual degrees of freedom, etc
  tempvar esamp
  qui gen byte `esamp' = e(sample)
  local df_r = e(df_r)

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
    local relyrlst
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


  * here, we check and see if the main treatment effects have been dropped from the regression for multicolinearity - my tests for this are imperfect, but still useful
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
              dis as error  "Error: a treated cohort-period dummy coefficient was dropped for multicolinearity."
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
              if `catastrophe' == 0 dis as text  "Caution: a treatment period coefficient interacted with something (likely, a term implementing a continuous treatment variable or an interactive control variable) was dropped for multicolinearity."
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
    if `catastrophe' > 0 dis as text "Note: the multicolinearity warning above was triggered `catastrophe' times total."
    if `bigcatastrophe' == 1 error 322
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
    local converged = e(converged)
    local chi2 = e(chi2)
  }
  else {
    local r2 = e(r2)
    local r2adj = e(r2_a)
    local r2_within = e(r2_within)
    local r2_withinadj = e(r2_a_within)
    local F = e(F)
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

  * if doing adversarial checks, inject the "post" statement where needed
  local postinject
  if "`adversarial'" == "adversarial" local postinject post

  * initialize some times
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
    if "`verbose'" == "verbose" dis as text "Wooldid is now computing `item' estimates for full sample, including event study estimates."

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
          qui margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `_t' >= `horizonstarter' & `_t' <= `horizonend' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof' `postinject'
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
          if "`postinject'" == "post" qui est sto __marginsmodel

          * this implements the f test that all points are well explained by a single line fit on the pre-period coefficients
          * the finaladvtest is for the adversarial test being estimable - mainly just an issue w/ SGs later. sometimes we have enough coefs to do the
          *     test overall, but not in a particular sg. since this is the overall section, if it triggers badly here, call off adversarial as a whole
          *     threshold of 2 here b/c that corresponds w/ 1 coef (targetflag==0 will get dummy columns )
          local finaladvtest
          cap confirm matrix r(b)
          if _rc == 0 {
            local finaladvtest  `=colsof(r(b))'
          }
          if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
          if `finaladvtest' <= 2 {
            local adversarial
            qui est restore ___wooldidfullmodel
            local postinject
          }
          if "`adversarial'" == "adversarial"     {
            tempname adversarialframe
            tempname pretrendregmat
            matrix `pretrendregmat' = r(b)[1,`matrixstartgrab'..`=colsof(r(b))']
            matrix `pretrendregmat' = `pretrendregmat''
            frame create `adversarialframe'
            frame change `adversarialframe'
            qui svmat double `pretrendregmat'
            qui rename `pretrendregmat'1 b
            qui gen t = .
            local bignorig = _N
            local bign = _N
            forval interpret_t = 1(1)`bign' {
              local special_pos = strpos(r(by`interpret_t'),".") - 1
              local special_pos = substr(r(by`interpret_t'),1,`special_pos')
              qui replace t = `special_pos' if _n == `interpret_t'
            }
            local bignjump = 1

						*for timetrends, add 2 0 ref periods for esf; not otherwise. for pooled, still count only 1 ref period
						*the best way to do this for pooled is ambiguous, but makes sense since it isnt strictly limited to 2 periods here
            if "`timetrends'" == "timetrends" & "`esfixedbaseperiod'" == "esfixedbaseperiod" local bignjump = 2
            local bign = `bign' + `bignjump'
            qui set obs `bign'
            if "`esfixedbaseperiod'" == "" {
              qui gstats sum t, meanonly
              qui replace b = 0 if b == .
              qui replace t = r(min) - 1 if _n == _N
            }
            else{
              qui replace b = 0 if b == .
              * the esf reference horizon is analytic: relyr = t - min_t - 1 and the
              * reference relyr is esrelativeto, so t == refperiod (worker L456). the
              * timetrends second reference row sits one horizon below (worker L478).
              qui replace t = `refperiod' if _n == _N
              if "`timetrends'" == "timetrends" qui replace t = `refperiod' - 1 if _n == _N - 1
            }


            * an alternate implementation might estimate the regression below using only pre-treatment coefficients, projecting them forward
            local altadv_restrict
            local altadv_restrict_num
            if "`altmethod_adversarial'" != "" local altadv_restrict if t <= `espre'
            if "`altmethod_adversarial'" != "" local altadv_restrict_num `espre'
            if "`altmethod_adversarial'" != "" & "`esfixedbaseperiod'" != "" {
              * esf: relyr = t - min_t - 1, so `t <= min_t' is exactly every pre-period
              * coefficient plus the fixed reference row(s), in either reference regime
              local altadv_restrict if t <= `min_t'
              local altadv_restrict_num = `min_t'
            }

            qui reg b t `altadv_restrict'
            tempname adversarialB_`horizitem'
            matrix `adversarialB_`horizitem''   = e(b)
            qui predict double linearproject, xb
            local testconstructor
            * this ensures that when projecting forward we etst only the post-treatment period effects
            forval advframe_i = 1(1)`bignorig' {
              local horiztocheck = t[`advframe_i']
              local proceed 0 
              if "`altadv_restrict_num'" == "" local proceed 1
              if "`altadv_restrict_num'" != "" {
                if `horiztocheck' > `altadv_restrict_num' local proceed 1
              }
              if `proceed' == 1 {
                local xbline = linearproject[`advframe_i']
                local testconstructor `testconstructor' ([1.`targetflag']`horiztocheck'.`_t' = `xbline')
             }
            }
            frame change `defaultframe'
            frame drop `adversarialframe'
            qui est restore __marginsmodel
            qui test `testconstructor'
            local adversarialP_`horizitem' =  r(p)
            qui est drop __marginsmodel
            qui est restore ___wooldidfullmodel
          }
          * end f test above
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
  if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'"=="" {
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
      if "`verbose'" == "verbose" dis as text "Wooldid is now working on average marginal effect estimation and production of any requested continuous treatment effects plots, including plots by subgroup."
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
            local titletext "Post-Treatment `effecttype' Across the `contreat' distribution"
          }
          if "`cfxplotcase'" == "cfxpre" & `cfxpre' ==1 {
            local cfxif & `t' < `ttre' & `ttre' != .
            local titletext "Pre-Treatment `effecttype' Across the `contreat' distribution"
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
          qui marginsplot, recastci(rarea) ciopts(color(navy%50)) plotopts(color(navy)) `cfx_xtitle' ytitle("`ytext'")  scale(.85) title("`titletext'") name("cfx`itemtext'_`nametext'", replace)
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
              qui marginsplot, plotopts(color(navy)) recastci(rarea) ciopts(color(navy%50)) `cfx_xtitle' ytitle("`ytext'") title("`titletext'")  scale(.85) name("cfx`itemtext'_`nametext'_sg`sgcase'", replace)
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
            qui  margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof' `postinject'
            matrix `f9tab' = r(table)
            if rowsof(`f9tab') < 6 {
              matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
              matrix rownames `f9tab' = b se t pvalue ll ul
            }
            matrix `___C`horizitem'_horiz' = (`f9tab'[1..6, 1..`=colsof(r(b))']',`relyrcolumn')
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
            if "`postinject'" == "post" qui est sto __marginsmodel

            * this implements the f test that all points are well explained by a single line fit on the pre-period, constrained to run through 0
            local finaladvtest
            cap confirm matrix r(b)
            if _rc == 0 {
              local finaladvtest  `=colsof(r(b))'
            }
            if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
            if `finaladvtest' <= 1 {
              local adversarial
              qui est restore ___wooldidfullmodel
              local postinject
            }
            if "`adversarial'" == "adversarial"    {
              tempname adversarialframe
              tempname pretrendregmat
              matrix `pretrendregmat' = r(b)[1,1..`=colsof(r(b))']
              matrix `pretrendregmat' = `pretrendregmat''
              frame create `adversarialframe'
              frame change `adversarialframe'
              qui svmat double `pretrendregmat'
              qui rename `pretrendregmat'1 b
              qui gen t = .
              local bignorig = _N
              local bign = _N
              forval interpret_t = 1(1)`bign' {
                local special_pos = strpos(r(by`interpret_t'),".") - 1
                local special_pos = substr(r(by`interpret_t'),1,`special_pos')
                qui replace t = `special_pos' if _n == `interpret_t'
              }
              local bignjump = 1
              if "`timetrends'" == "timetrends" & "`esfixedbaseperiod'" == "esfixedbaseperiod" local bignjump = 2
              local bign = `bign' + `bignjump'
              qui set obs `bign'
              if "`esfixedbaseperiod'" == "" {
                qui gstats sum t, meanonly
                qui replace b = 0 if b == .
                qui replace t = r(min) - 1 if _n == _N
              }
              else{
                qui replace b = 0 if b == .
                * the esf reference horizon is analytic: relyr = t - min_t - 1 and the
                * reference relyr is esrelativeto, so t == refperiod (worker L456). the
                * timetrends second reference row sits one horizon below (worker L478).
                qui replace t = `refperiod' if _n == _N
                if "`timetrends'" == "timetrends" qui replace t = `refperiod' - 1 if _n == _N - 1
              }

              local altadv_restrict
              if "`altmethod_adversarial'" != "" local altadv_restrict if t <= `espre'
              local altadv_restrict_num
              if "`altmethod_adversarial'" != "" local altadv_restrict_num  `espre'
              if "`altmethod_adversarial'" != "" & "`esfixedbaseperiod'" != "" {
                * esf: relyr = t - min_t - 1, so `t <= min_t' is exactly every pre-period
                * coefficient plus the fixed reference row(s), in either reference regime
                local altadv_restrict if t <= `min_t'
                local altadv_restrict_num = `min_t'
              }

              qui reg b t `altadv_restrict'
              tempname adversarialB_C`horizitem'
              matrix `adversarialB_C`horizitem'' = e(b)
              qui predict double linearproject, xb
              local testconstructor

              forval advframe_i = 1(1)`bignorig' {
                local horiztocheck = t[`advframe_i']
                local proceed 0 
                if "`altadv_restrict_num'" == "" local proceed 1
                if "`altadv_restrict_num'" != "" {
                  if `horiztocheck' > `altadv_restrict_num' local proceed 1
                }
                if `proceed' == 1 {
                  local xbline = linearproject[`advframe_i']
                  local testconstructor `testconstructor' (`horiztocheck'.`_t' = `xbline')
                }
              }
              frame change `defaultframe'
              frame drop `adversarialframe'
              qui est restore __marginsmodel
              qui test `testconstructor'
              local adversarialP_C`horizitem' = r(p)
              qui est drop __marginsmodel
              qui est restore ___wooldidfullmodel
            }
            * end f test above
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
      if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'"=="" {
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
  if "`verbose'" == "verbose" dis as text "Wooldid has finished production of average continuous treatment effects plots and all full treatment sample average marginal effects."








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
              * problem: this doens't guarantee we see the same relative periods, however
              qui margins [`weight'=`w_`horizitem''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof' `postinject'
              local matrixstartgrab = floor(`=colsof(r(b))'/2)+1
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
              *the above will be udpated w/ rownames / relyrs by the below

              * SGs need not have the same number of relative periods appear for each
              * this code block below figures out what relative periods appear for them
              tempname nameget
              matrix `nameget' = r(table)[1, `matrixstartgrab'..`=colsof(r(b))']
              local names: colnames `nameget'
              matrix drop `nameget'
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
              if "`postinject'" == "post" qui est sto __marginsmodel

              * this implements the f test that all points are well explained by a single line fit on the pre-period, constrained to run through 0
              local finaladvtest
              cap confirm matrix r(b)
              if _rc == 0 {
                local finaladvtest  `=colsof(r(b))'
              }
              if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
              if "`adversarial'" == "adversarial"   {
                if  `finaladvtest' > 2  {
                  tempname adversarialframe
                  tempname pretrendregmat
                  matrix `pretrendregmat' = r(b)[1,`matrixstartgrab'..`=colsof(r(b))']
                  matrix `pretrendregmat' = `pretrendregmat''
                  frame create `adversarialframe'
                  frame change `adversarialframe'
                  qui svmat double `pretrendregmat'
                  qui rename `pretrendregmat'1 b
                  qui gen t = .
                  local bignorig = _N
                  local bign = _N
                  forval interpret_t = 1(1)`bign' {
                    local special_pos = strpos(r(by`interpret_t'),".") - 1
                    local special_pos = substr(r(by`interpret_t'),1,`special_pos')
                    qui replace t = `special_pos' if _n == `interpret_t'
                  }
                  local bignjump = 1
                  if "`timetrends'" == "timetrends" & "`esfixedbaseperiod'" == "esfixedbaseperiod" local bignjump = 2
                  local bign = `bign' + `bignjump'
                  qui set obs `bign'
                  if "`esfixedbaseperiod'" == "" {
                    qui gstats sum t, meanonly
                    qui replace b = 0 if b == .
                    qui replace t = r(min) - 1 if _n == _N
                  }
                  else{
                    qui replace b = 0 if b == .
                    * the esf reference horizon is analytic: relyr = t - min_t - 1 and the
                    * reference relyr is esrelativeto, so t == refperiod (worker L456). the
                    * timetrends second reference row sits one horizon below (worker L478).
                    qui replace t = `refperiod' if _n == _N
                    if "`timetrends'" == "timetrends" qui replace t = `refperiod' - 1 if _n == _N - 1
                  }

                  local altadv_restrict
                  local altadv_restrict_num
                  if "`altmethod_adversarial'" != "" local altadv_restrict if t <= `espre'
                  if "`altmethod_adversarial'" != "" local altadv_restrict_num  `espre'
                  if "`altmethod_adversarial'" != "" & "`esfixedbaseperiod'" != "" {
                    * esf: relyr = t - min_t - 1, so `t <= min_t' is exactly every pre-period
                    * coefficient plus the fixed reference row(s), in either reference regime
                    local altadv_restrict if t <= `min_t'
                    local altadv_restrict_num = `min_t'
                  }
                  qui reg b t  `altadv_restrict'
                  tempname adversarialB_`horizitem'_s`sgcase'
                  matrix `adversarialB_`horizitem'_s`sgcase'' = e(b)
                  qui predict double linearproject, xb
                  local testconstructor
                  forval advframe_i = 1(1)`bignorig' {
                    local horiztocheck = t[`advframe_i']
                    local proceed 0 
                    if "`altadv_restrict_num'" == "" local proceed 1
                    if "`altadv_restrict_num'" != "" {
                      if `horiztocheck' > `altadv_restrict_num' local proceed 1
                    }
                    if `proceed' == 1 {
                      local xbline = linearproject[`advframe_i']
                      local testconstructor `testconstructor' ([1.`targetflag']`horiztocheck'.`_t' = `xbline')
                    }
                  }
                  frame change `defaultframe'
                  frame drop `adversarialframe'
                  qui est restore __marginsmodel
                  qui test `testconstructor'
                  local adversarialP_`horizitem'_s`sgcase' = r(p)
                  qui est drop __marginsmodel
                }
              * the above ends f test; this restore is broken out in case the "not enough pre-trend coefs" clause triggers for the sg
              qui est restore ___wooldidfullmodel
              * the below ends the adversarial check, different form the inner one just so as to handle the est restore appropriately
              }
            }
          }
        }
      }

      * produce joint tests of pre-treatment horizons by sg - accomplish by looping through sgs
      if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'"=="" {

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


    * joint test the sg effect estimates against eachother - are they the same or different?
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
          }

          * average marginal effects by horizon
          if `espre' + `espl' > 0  {
            local horizitem att
            if "`item'" != "att" local horizitem att_it

            foreach sgcase in `sglist' {
              cap confirm  matrix `___`horizitem'_Choriz_s`sgcase''
              * this prevents estimating horizons if done already for a non att
              if _rc != 0 {
                qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstext' `marginsdof' `postinject'
                tempname ___`horizitem'_Choriz_s`sgcase'
                matrix `f9tab' = r(table)
                if rowsof(`f9tab') < 6 {
                  matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
                  matrix rownames `f9tab' = b se t pvalue ll ul
                }
                matrix `___`horizitem'_Choriz_s`sgcase'' = ( `f9tab'[1..6, 1..`=colsof(r(b))']',`relyrcolumn_s`sgcase'')
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
                if "`postinject'" == "post" qui est sto __marginsmodel

                * this implements the f test that all points are well explained by a single line fit on the pre-period, constrained to run through 0
                local finaladvtest
                cap confirm matrix r(b)
                if _rc == 0 {
                  local finaladvtest  `=colsof(r(b))'
                }
                if "`finaladvtest'" == ""  | "`finaladvtest'" == "." local finaladvtest 0
                if "`adversarial'" == "adversarial"   {
                  if  `finaladvtest' > 1  {
                    tempname adversarialframe
                    tempname pretrendregmat
                    matrix `pretrendregmat' = r(b)[1,1..`=colsof(r(b))']
                    matrix `pretrendregmat' = `pretrendregmat''
                    frame create `adversarialframe'
                    frame change `adversarialframe'
                    qui svmat double `pretrendregmat'
                    qui rename `pretrendregmat'1 b
                    qui gen t = .
                    local bignorig = _N
                    local bign = _N
                    forval interpret_t = 1(1)`bign' {
                      local special_pos = strpos(r(by`interpret_t'),".") - 1
                      local special_pos = substr(r(by`interpret_t'),1,`special_pos')
                      qui replace t = `special_pos' if _n == `interpret_t'
                    }
                    local bignjump = 1
                    if "`timetrends'" == "timetrends" & "`esfixedbaseperiod'" == "esfixedbaseperiod" local bignjump = 2
                    local bign = `bign' + `bignjump'
                    qui set obs `bign'
                    if "`esfixedbaseperiod'" == "" {
                      qui gstats sum t, meanonly
                      qui replace b = 0 if b == .
                      qui replace t = r(min) - 1 if _n == _N
                    }
                    else{
                      qui replace b = 0 if b == .
                      * the esf reference horizon is analytic: relyr = t - min_t - 1 and the
                      * reference relyr is esrelativeto, so t == refperiod (worker L456). the
                      * timetrends second reference row sits one horizon below (worker L478).
                      qui replace t = `refperiod' if _n == _N
                      if "`timetrends'" == "timetrends" qui replace t = `refperiod' - 1 if _n == _N - 1
                    }
                   local altadv_restrict
                   local altadv_restrict_num
                   if "`altmethod_adversarial'" != "" local altadv_restrict if t <= `espre'
                   if "`altmethod_adversarial'" != "" local altadv_restrict_num  `espre'
                   if "`altmethod_adversarial'" != "" & "`esfixedbaseperiod'" != "" {
                     * esf: relyr = t - min_t - 1, so `t <= min_t' is exactly every pre-period
                     * coefficient plus the fixed reference row(s), in either reference regime
                     local altadv_restrict if t <= `min_t'
                     local altadv_restrict_num = `min_t'
                   }

                    qui reg b t  `altadv_restrict'
                    tempname adversarialB_C`horizitem'_s`sgcase'
                    matrix `adversarialB_C`horizitem'_s`sgcase'' = e(b)
                    qui predict double linearproject, xb
                    local testconstructor
                    forval advframe_i = 1(1)`bignorig' {
                      local horiztocheck = t[`advframe_i']
                      local proceed 0 
                      if "`altadv_restrict_num'" == "" local proceed 1
                      if "`altadv_restrict_num'" != "" {
                        if `horiztocheck' > `altadv_restrict_num' local proceed 1
                      }
                      if `proceed' == 1 {
                        local xbline = linearproject[`advframe_i']
                        local testconstructor `testconstructor' (`horiztocheck'.`_t' = `xbline')
                      }
                    }
                    frame change `defaultframe'
                    frame drop `adversarialframe'
                    qui est restore __marginsmodel
                    qui test `testconstructor'
                    local adversarialP_C`horizitem'_s`sgcase' = r(p)
                    qui est drop __marginsmodel
                  }
                  * the above ends f test; this restore is broken out in case the "not enough pre-trend coefs" clause triggers for the sg
                qui est restore ___wooldidfullmodel
                }
              }
            }
          }
          * end item loop
        }

        * produce joint tests of pre-treatment AMES by horizon by sg
        if (  ("`esfixedbaseperiod'" == "" & `espre' > 1) |  ( `espre' > 1 & "`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) )  & "`suppressjointtests'"=="" {
          foreach sgcase in `sglist' {
            if "`att'" == "att" {
              qui margins  [`weight'=`w_att'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & `subgroups' == `sgcase'   &  (`_t' == 0 |`_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' `marginsdof'
              cap confirm matrix r(p)
              if _rc == 0 {
                local joint_Choriz_s`sgcase' = r(p)[1,2]
              }
              else{
                local joint_Choriz_s`sgcase'
              }
            }

            if ("`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" ) & "`joint_Choriz_i_s`sgcase''" == "" {
              qui margins [`weight'=`w_att_it'] ,`contreatelasticitytype'(`contreat') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & `subgroups' == `sgcase'  &  (`_t' == 0 |`_t' >= `horizonstarter')    ) contrast(overjoint) `poismarginstext' `marginsdof'
              cap confirm matrix r(p)
              if _rc == 0 {
                local joint_Choriz_i_s`sgcase' = r(p)[1,2]
              }
              else{
                local joint_Choriz_i_s`sgcase'
              }
            }
            * these parameters give you the 2nd half of the sg-specific joint test matrix, minus the joint test of all sgs
          }
        }

      }



      *************** joint test the sg AMEs  against eachother -
      * E1: see the binary twin above - one posted pwcompare pass per (estimand x pre/post)
      * serves both the pairwise table and the all-subgroups-equal joint test
      if wordcount("`sglist'") > 1 {
        * F9: see the binary twin - no post, no joint test under fast mode
        local sgpost post
        if "`fastnose'" != "" local sgpost

        foreach item in `att' `att_it' `att_i' `att_itime' {

          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstext' noestimcheck over(`subgroups') pwcompare(effects) `marginsdof' `sgpost'
          tempname pwsgs_C`item'
          matrix `f9tab' = r(table_vs)
          if rowsof(`f9tab') < 6 {
            matrix `f9tab' = `f9tab' \ J(`=6-rowsof(`f9tab')', `=colsof(`f9tab')', .)
            matrix rownames `f9tab' = b se t pvalue ll ul
          }
          matrix `pwsgs_C`item''  = `f9tab'[1..6,1..`=colsof(r(table_vs))']'
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
        foreach rownametest in `rownamespool' {
          if strpos("`rownametest'","ES_")==0 {
            if strpos("`rownametest'","pre")==0 {
              local rownamessg `rownamessg' `rownametest'
            }
            else{
              if "`sghaspre'" == "yes" local rownamessg `rownamessg' `rownametest'
            }
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
      matrix colnames `pwsgs' =  estimate se t p lb_95ci ub_95ci

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
  * composition is a diagnostic aid, not a test: flat columns rule composition out as the cause of
  * a jumpy profile; moving columns do not by themselves establish composition bias.
  * this block runs NO margins call and NO e-class command, so the active e() is untouched and no
  * est restore is needed. it DOES clear r(); nothing between here and the ereturn post consumes
  * r(). gtools by() is hash-based: the user's sort order is not disturbed.
  local f8keys
  if `espre' + `espl' > 0 {
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
  * this means we will get the correct coef and SE when using esttab... but you won't be able to test estimates against eachother, use lincom to combine estimates, etc.
  if "`subgroups'" != "" {
    if "`suppresssgprimaryeffects'" == "" {
      local sgnum 0
      foreach sgcase in `sglist' {
        matrix `ebformer' = `ebformer' \ `wdid_sg`sgcase''
        local ++sgnum
        local rownamespool_sginj = subinstr("`rownamessg`sgcase''",":","_SG`sgcase'_",.)
        local rownamespool `rownamespool' `rownamespool_sginj'
      }
    }
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
        * skip a fully non-identified VCV (e.g. timetrends): no positive-finite diagonal
        local f2ok 0
        forval f2d = 1/`=colsof(`f2vm')' {
          if `f2vm'[`f2d',`f2d'] > 0 & `f2vm'[`f2d',`f2d'] < . local f2ok 1
        }
        if `f2ok' {
          cap ereturn matrix b_es_`f2ky' = `f2bm'
          cap ereturn matrix V_es_`f2ky' = `f2vm'
          cap ereturn scalar numpre_`f2ky' = `f2np'
        }
      }
    }
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
  ereturn local cmdline `"`0'"'
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
       if "`adversarialP`scenario'`type''" != "" {
         local typetext `type'
         if "`scenario'" == "_C" local typetext = subinstr("`typetext'","att","ame",.)
         ereturn scalar adversarial_pval_`typetext' = `adversarialP`scenario'`type''
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
              if "`adversarialP`scenario'`type'_s`sgcase''" != "" {
                local typetext `type'
                if "`scenario'" == "_C" local typetext = subinstr("`typetext'","att","ame",.)
                ereturn scalar adversarial_pval_`typetext'_s`sgcase' = `adversarialP`scenario'`type'_s`sgcase''
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
  foreach matrix in `cfxmatricestoreturn' {
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
  ereturn hidden local wdid_replayv "1"


  ********* begin actual result presentation *************
  * F11: the presentation code lives in wooldid_display (non-eclass, reads e() only).
  * it is split in two parts so the eclass wooldid_ptp call below stays at worker level
  * (its ereturns would be discarded inside a non-eclass subprogram)
  wooldid_display , part(1)
  * F7: Roth (2022) pre-trend power diagnostics (compute + display + ereturns live in a
  * subprogram; cap noi so a failure here cannot abort an otherwise-complete run)
  if "`pretrendpower'" != "" cap noi wooldid_ptp , levels(`pretrendpower') esrel(`esrelativeto') dfr(`df_r') pois(`poisson') poisexp(`poisexpresults') tvar(`t') ttrevar(`ttre') ivar(`i') wvar(`w') cwm(`customweightmultiplier')
  wooldid_display , part(2)




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

  * __marginsmodel is a scratch intermediate; the adversarial branches that store it do
  * not all reach their own est drop (the finaladvtest short-circuits skip it), so drop
  * it unconditionally here (B22; this also removes the _est___marginsmodel marker
  * variable the store planted in the dataset)
  cap est drop __marginsmodel
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
syntax [varlist(numeric default=none)] [if/]  [aw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ANTICipation(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace  altmethod_adversarial ADVersarial update oldsyntax PRETRENDPOWer(numlist) FASTnose]

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

  if "`verbose'" == "verbose" dis as text "Wooldid is processing input in preparation of regression estimation."

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

  *Confirm correct weight usge if using poisson
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
    dis as text "Set option safety(off) to suppress this warning message."
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

  * extract contreatcontrol type and error out if misspecifeid
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
      dis as text "As such, wooldid will not produce these plots unless you specify safety = off."
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

  * capture the user's own safety() input before any internal relaxation (contreatwithin sets
  * the working copy to off further down); the identifying-assumption checks gate on this copy
  * so an internal relaxation does not silently disable them
  local safetyinput `safety'

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

  * if no elasticity type specified, assume regular non-elasticity estimates
  if "`contreatelasticitytype'" != "eyex" & "`contreatelasticitytype'" != "dyex" & "`contreatelasticitytype'" != "eydx" local contreatelasticitytype dydx
  if "`semielasticity'" == "" {
    local semielasticity dydx
  }
  else {
    local semielasticity eydx
  }

  * do some light processing of the "number of pre-treatment coefs" input - it needs some adjustment when esfixedbaseperiod is used since the period normalized to 0 could fall within your range of requested coefs
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
    if strpos("`cfxplottypes'","post") != 0 local cfxpost 1
    local cfxplottypes = subinstr("`cfxplottypes'","post"," ",.)
    if strpos("`cfxplottypes'","pre") != 0 & ("`esfixedbaseperiod'" == "esfixedbaseperiod" | `esprelength' > 0 ) local cfxpre 1
    local cfxplottypes = subinstr("`cfxplottypes'","pre"," ",.)
    foreach item in `cfxplottypes' {
      * entries other than the (already-stripped) post/pre keywords must be numeric: a stray token
      * would otherwise be evaluated as a variable reference below -- crashing with an unexplained
      * r(111) if no such variable exists, or silently injecting the variable into the plot list if
      * one does (stata does not short-circuit |, so the reference-period comparison always runs)
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
      if "`esfixedbaseperiod'" == "" | `item' != `esrelativeto' {
        cap if int(`item')==`item' & `item' >= 0 & (`item' <= `espostlength' | "`esfixedbaseperiod'" == "esfixedbaseperiod") local cfxrelyrlist `cfxrelyrlist' `item'
        cap if int(`item')==`item' & `item' < 0 & (abs(`item') <= `esprelength' | "`esfixedbaseperiod'" == "esfixedbaseperiod") local cfxrelyrlist `cfxrelyrlist' `item'
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

  * subgroups cannot be negative
  if "`subgroups'" != "" {
    cap assert `subgroups' >= 0 `ifassert'
    if _rc != 0 {
      dis as error "Error: the subgroup variable contains negative values; it must contain only non-negative integers."
      error 411
    }
  }

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
      cap assert (`cfxlattice' >= 0 & int(`cfxlattice') == `cfxlattice') | `cfxlattice' >= . `ifassert'
      if _rc != 0 {
        dis as error "Error: the cfxlattice() variable must contain only non-negative integers (margins allows only non-negative integers in over())."
        error 411
      }
      local wdvcfxmiss `cfxlattice' >= .
      if `"`if'"' != "" local wdvcfxmiss (`wdvcfxmiss') & ((`if') == 1)
      qui count if `wdvcfxmiss'
      if r(N) > 0 dis as text "Caution: the cfxlattice() variable is missing on `r(N)' observation(s); margins omits those observations when computing the continuous treatment effect plot points, so the plotted averages reflect only the nonmissing-cfxlattice subsample."
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
      local wdvtlo = r(min) + `anticipation'
      local wdvthi = r(max) + `oostreatedcontrols' + `anticipation'
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

  * ---- return the parsed options and every derived/normalized local to the worker ----
  * blanket-by-construction: every option the syntax line defines is listed (plus the
  * derived locals below), so a future edit that mutates one more local here cannot
  * silently strand the worker on the stale syntax-parsed value. varlist/if/weight/exp
  * are excluded: the worker's own syntax call already parsed them identically and this
  * subprogram never alters them.
  foreach wdv_ret in att att_it att_itime att_i cluster customvce unconditionalse robust impvalues fe coarsecohortcontrols ccc_absorb controls interactivecontrols timetrends subgroups customweightmultiplier esprelength espostlength saveplots poisson poisexpresults oostreatedcontrols anticipation esfixedbaseperiod esrelativeto safety makeplots makecfxplots summarystats histogramcohorteffects regtol jointtests suppresssgprimaryeffects contreat contreatpoly lattice latticeignoreweights ccnoabsorb contreatcontrols verbose cfxplottypes cfxlattice makecfxplotsbysg cfxplotnose semielasticity contreatelasticitytype contreatwithin cleanmatrices emptycellsoverride replace altmethod_adversarial adversarial update oldsyntax pretrendpower fastnose y i t ttre_asgiven safetyinput suppressjointtests latticenum latticetype esposttrip cfxpost cfxpre cfxrelyrlist {
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
    dis as text "Caution: anticipation(`antic') places the effective reference period of `r(N)' unit(s) on a period they are not observed in; those unit(s) and their `anticdropobs' observations are dropped (their reference period would have been observed without the anticipation shift)."
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
    dis as result "    pooled estimate if it continued into the post-treatment period:"
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
cap prog drop wooldid_display
prog def wooldid_display
version 16.0
  syntax , part(string)
  if "`e(wdid_replayv)'" == "" {
    dis as error "Error: the stored results were posted by a wooldid version without replay support; nothing can be redisplayed."
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
  dis as result  "Wooldid Estimation for Outcome: `y'; Standard Errors: `vcespecifier'"
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
      dis as result  "R2 = " %-6.4f `r2' ";  R2adj = " %-6.4f `r2adj' "; R2-within = " %-6.4f `r2_within' ";  R2-withinadj = " %-6.4f `r2_withinadj'
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
  * re-running the quadrature (the NIS computation can cost minutes). estimation-time
  * Cautions from wooldid_ptp are deliberately not replayed. on the live path the worker
  * calls wooldid_ptp itself between part(1) and part(2), so this part is never run live
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
    dis as result "    pooled estimate if it continued into the post-treatment period:"
    matlist e(roth_power_`ptpty'), format(%12.4g)
    dis as text "    note: mdt_wald powers the joint pre-trend test above; mdt_nis powers the test that any"
    dis as text "    pre-period coefficient is individually significant (the convention of Roth's pretrends package)."
    if "`poisson'" == "poisson" & "`e(wdid_poisexp)'" == "" dis as text "    note: poisson event-study coefficients are on the linear-index (xb) scale, and so are these diagnostics."
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
        if "`poisson'" == "poisson" local poistext "Poisson Regression"
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
          if "`e(joint_`itemtext'_allsgssame)'" != "" dis as result  "p-value for joint test that all SGs have the same `itemtext': " %-8.4f e(joint_`itemtext'_allsgssame)
          if "`e(joint_pre_`itemtext'_allsgssame)'" != "" dis as result  "p-value for joint test that all SGs have the same pre-treatment `itemtext': " %-6.4f e(joint_pre_`itemtext'_allsgssame)
        }
        dis as result  " "
      }
    }
  }
  }
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
  if "`e(wdid_replayv)'" == "" {
    dis as error "Error: wooldid_esplots is internal to wooldid and needs a stored wooldid estimation in e(); run wooldid first."
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
      local initialflagt `esrelativeto'
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
          local initialflagt `esrelativeto'
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

        if "`saveplots'" != ""  qui save "`saveplots'_ES_`subgroups'`sgcase'.dta", `replace'
        clear
      }
    }

    frame change `defaultframe'
    frame drop `__wooldidframe'
end

