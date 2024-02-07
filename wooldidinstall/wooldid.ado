
* Wooldid Version 1.1.1, 2/07/2024
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

syntax [varlist(numeric default=none)] [if/]  [aw pw / ]  ,  [att att_it att_itime att_i CLUSter(varlist) customvce(string) UNCONDitionalse  ROBust IMPvalues fe(varlist fv) COARSEcohortcontrols(string) ccc_absorb(varlist numeric)   CONTROLs(varlist fv) INTERactivecontrols(varlist numeric)  TIMETRENDs SUBgroups(varname numeric) CUSTOMWeightmultiplier(varname numeric)  ESPRElength(integer 0) ESPOSTlength(integer 0) SAVEplots(name)  POISson POISEXPresults   OOStreatedcontrols(integer 0) ESFixedbaseperiod ESRelativeto(integer -1) safety(string)   MAKEPlots MAKECfxplots SUMmarystats HISTogramcohorteffects   regtol(integer 9)  JOINTtests SUPPRESSsgprimaryeffects contreat(varname numeric) CONTREATPoly(int 1) lattice(string) LATTICEIGnoreweights  CCNOabsorb CONTREATControls(string) VERbose   CFXPLOTTypes(string) cfxlattice(varname) makecfxplotsbysg cfxplotnose SEMIelasticity CONTREATELASticitytype(string)    CONTREATWithin CLEANmatrices   emptycellsoverride  replace   ADVersarial update oldsyntax]

* Execute update before running main program
  if "`update'" == "update" {
    net install wooldid, from(https://raw.githubusercontent.com/thegland/wooldid/master/wooldidinstall/) replace
    dis as text "Update complete; exiting wooldid. Other dependencies (reghdfe, ftools, ppmlhdfe (optional), gtools (optional)) must be updated manually if necessary."
    exit
  }

  ***************************************************************************************************
  ******************************  Initialize and Process Command Input ******************************
  ***************************************************************************************************
  * This section converts user input into actionable information and checks for validity of the user input.


  if "`verbose'" == "verbose" dis as text "Wooldid is processing input in preparation of regression estimation."

  * Clean the workspace if cleanmatrices is specified - this is helpful b/c sometimes reghdfe enters into a failure state
  * that leaves junk matrices floating around in stata and mata which prevent further estimation. cleanmatrices should fix that.
  if "`cleanmatrices'" != "" {
    clear mata
    matrix drop _all
  }

 * Unless disabled by the user, set emptycells to drop -- greatly improves performance since structure of the reghdfe matrix
 * involves lots of empty cells in the interaction term.
 if "`emptycellsoverride'" == "" {
    set emptycells drop
 }

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

  if "`safety'" != "off" & "`unconditionalse'" == "unconditionalse" & ("`weight'" == "aw" ) {
    dis " "
    dis as text "Error: unconditional standard errors may be unreliable with aweights."
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


  *check conditions for using im p-values - in principle, SG estimates could be computed provided SGs never cross cut clusters, but I have not implemented this
  if "`impvalues'" != "" & "`subgroups'" != "" {
    dis as error "Error: Ibragimov and Muller p-values not compatible with subgroups."
    error 198
  }

  *check conditions for using im p-values - for CFX plots, we have problems similar to with SGs, but even worse
  if "`impvalues'" != "" & "`makecfxplots'" != "" {
    dis as error "Error: Ibragimov and Muller p-values not compatible with continuous treatment effects plots."
    error 198
  }

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
    if "`latticenum'" == "1" {
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
    local numbase 0
    foreach word in `contreatcontrols' {
      if inlist("`word'",  "byTreated", "byTtre", "byCohort", "byTime") == 0 & inlist("`word'","basic","untreatedZero") == 0 {
        dis as error  "Error: continuous treatment control type misspecified; `word' is not a valid continuous treatment control type."
        error 198
      }
      if strpos("`word'","base") != 0 local ++numbase
    }
    if `numbase' > 1  {
      dis as error "Error: cannot specify more than one base[]ByTime control set."
      error 198
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

 * issue warnings pertaining to cfxplot types
  if "`makecfxplots'" != "" & ("`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" )   {
    if "`safety'" == "off" {
      dis as text "Caution: interpretation of continuous treatment effects plots for estimands other than the ATT is difficult unless i-t cells (and sometimes cohorts i) are wholly contained within lattice bins."
    }
    else{
      dis as error "Caution: interpretation of continuous treatment effects plots for estimands other than the ATT is difficult unless i-t cells (and sometimes cohorts i) are wholly contained within lattice bins."
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

  * check that i is specified acceptably
  if "`safety'" != "off" {
    cap assert `i' > 0  & int(`i') == `i'
    if _rc != 0 {
      dis as error  "Error: `i' must consist of positive, non-zero integers, but does not"
      error 411
    }
  }

  * if poisson, check that y is not negative
  if "`poisson'" != "" {
    cap assert `y' >= 0
    if _rc != 0 {
      dis as error "The outcome variable contains negative values; this is not allowed with poisson."
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


  * check if gtools is present and enable program to use it if so
  local gtoolscheck
  local gtoolscheck2
  cap which gtools
  if _rc == 0{
    local gtoolscheck g
    local gtoolscheck2 gstats
  }

  * object to centering the interactive controls without gstats
  if "`interactivecontrols'" != "" & "`gtoolscheck2'" != "gstats" {
    dis as error "Error: the option interactivecontrols() requires that the package gtools be installed. It does not appear to be installed."
    error 198
  }



  * process cfxplottypes
  * note that the tests here on what cfxplottypes are allowed are kinda weak in that you can request impractical relative periods (relyrs) and such, but better than no checks
  if "`makecfxplots'" != "" {
    local cfxpost 0
    local cfxpre 0
    local cfxrelyr 0
    local cfxrelyrlist
    if "`cfxplottypes'" == "" local cfxplottypes post
    if strpos("`cfxplottypes'","post") != 0 local cfxpost 1
    local cfxplottypes = subinstr("`cfxplottypes'","post"," ",.)
    if strpos("`cfxplottypes'","pre") != 0 & ("`esfixedbaseperiod'" == "esfixedbaseperiod" | `esprelength' > 0 ) local cfxpre 1
    local cfxplottypes = subinstr("`cfxplottypes'","pre"," ",.)
    foreach item in `cfxplottypes' {
      if "`esfixedbaseperiod'" == "" | `item' != `esrelativeto' {
        cap if int(`item')==`item' & `item' >= 0 & (`item' <= `espostlength' | "`esfixedbaseperiod'" == "esfixedbaseperiod") local cfxrelyrlist `cfxrelyrlist' `item'
        cap if int(`item')==`item' & `item' < 0 & (abs(`item') <= `esprelength' | "`esfixedbaseperiod'" == "esfixedbaseperiod") local cfxrelyrlist `cfxrelyrlist' `item'
      }
    }
    if "`cfxrelyrlist'" != "" local cfxrelyr 1
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
  if  "`safety'" != "off" {
    cap assert `w' >= 0
    if _rc != 0 {
      dis as error  "Error: weights must not be negative."
      error 411
    }
  }


  * subgroups cannot be negative
  if "`subgroups'" != "" {
    cap assert `subgroups' >= 0
    if _rc != 0 {
      dis as error "The subgroup variable contains negative values; it must contain only non-negative integers."
      error 411
    }
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
  qui replace `ifvar' =  0 if `i' == . | `t' == . | `y' == . | `w' == .  | `w' <= 0

  * sanitize controls passed that may contain factor variable notation
  local loopablecontorls
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
  cap assert `subgroups' >= 0
  foreach mustbenonmissingvar in `checkclus'  `loopablecontrols' `ccc_absorb' `customweightmultiplier' `contreat' `checklattice' `includesg'  `interactivecontrols' {
    qui replace `ifvar' = 0 if `mustbenonmissingvar' == .
  }


  local ifstatement `ifvar' == 1





  *********************************************************************************************************************************
  ******************************  PROCESS EVENT STUDY INFO, ESTIMABILITY CHECKS, FORM TREATMENT DUMMIES ***************************
  *********************************************************************************************************************************
  * In this section, we form treated_cohort X treated_period dummies
  * and also figure out how large the panel is and how many observations we have across specific time periods - imposing
  * estimability constraints where appropriate.
  * We do this differently depending on if we are using a fixed or an pooled reference period.

  `gtoolscheck2' sum `t' if `ifstatement', meanonly
  local maxyearobserved = r(max)
  local minyearobserved = r(min)
  tempvar ttre

  * create a modified version of the period-of-treatment-onset variable ttre
  * the modified version will treat cohorts with treatment onsets a set distance (specified in oostreatedcontrols)
  * from the end of the sample identically as to units that are never treated at all
  qui gen `ttre' = `ttre_asgiven' if `ttre_asgiven' <= `maxyearobserved' + `oostreatedcontrols' & `ifstatement'

  * prepare the _i and _t variables, which are cohort and time variables where non-treated-unit non-treated-times are all pooled
  * noting that "non-treated-times" should be extended here to actually mean "any time period for treated units in its reference period"
  tempvar _i
  tempvar _t

  * we begin by handling the case where the user wants the comparison period to be the beginning of the sample, possibly multiple such periods
  if "`esfixedbaseperiod'" != "esfixedbaseperiod" {

    *  calculate the minimum time unit observed w/in the sample for each unit i
    * make a variable to drop all i where we see no pre period
    tempvar mint

    if "`gtoolscheck'" == "g" {
      qui gegen `mint' = min(`t') if `ifstatement', by(`i')
    }
    else {
      qui bysort `i': egen `mint' = min(`t') if `ifstatement'
    }


    * if requesting time trends, figuring the 2nd from the minimum time unit observed for each i
    if "`timetrends'" == "timetrends" {
      tempvar mintp1
      tempvar tp1
      qui gen `tp1' = `t' if `t' != `mint'
      if "`gtoolscheck'" == "g" {
        qui gegen `mintp1' = min(`tp1' ) if `ifstatement' , by(`i')
      }
      else {
        qui bysort `i': egen `mintp1' = min(`tp1') if `ifstatement'
      }
      qui drop `tp1'
    }

    * drop cohorts which do not have a pre-treatment period (or 2 such periods if using time trends)
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
          dis as result " "
          dis as error "Warning: Reference period for some units forced (for want of panel length) to include time period for which treatment effects were requested."
          dis as result "This does not break estimation, but event study estimates can carry counter intuitive meanings and reflect effect levels by cohort."
          dis as result "Recommendation: use option esfixedbaseperiod to change estimation approach to using a specific time period relative to treatment as your reference period."
          dis as result "(Or request estimates for fewer pre-treatment horizons.)"
          dis as result " "
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
          dis as error "Warning: Reference period for some units forced (for want of panel length) to include time period for which treatment effects were requested."
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
      if "`gtoolscheck'" == "g" qui gegen `controlcheck2' = max(`controlcheck') if `ifstatement' , by(`t')
      if "`gtoolscheck'" != "g" qui bysort `t': egen `controlcheck2' = max(`controlcheck') if `ifstatement'
      qui replace `ifvar'  = 0 if `controlcheck2' == 0

      * update _t and ttre
      qui replace `_t' = . if `ifvar' == 0
      qui replace `ttre'  = . if `ifvar' == 0

      if "`safety'" != "off"  {
        cap assert `controlcheck2' != 0
        if _rc != 0 {
          dis as text "Estimability checks found time periods exist in sample without a control unit, causing them to be dropped."
          if "`contreat'" != "" dis as text "If intending to only use within treated-cohorts variation in your continuous treatment variable, specify option contreatwithin to prevent these time periods from being dropped."
        }
      }
      qui drop `controlcheck' `controlcheck2'
    }

    *get min viable _t and max; set up relative period variable to have no negative values
    qui `gtoolscheck2' sum `_t' if `ifstatement', meanonly
    local min_t = abs(r(min))
    local max_t = r(max) + 1
    qui replace `_t' = `_t' + `min_t' + 1
    qui replace `_t' = 0 if `_t' == . & `ifstatement'

    * check that some cohort has a post period and that some has a pre-period
    cap assert `max_t' > 0 & `max_t' < .
    if _rc != 0 {
      dis as error  "Error: cannot find either a treated time period for any treated units or sufficient untreated periods for treated units."
      if "`timetrends'" != "" dis as error  "Note that including cohort-specific time trends raises the number of untreated periods required."
      error 459
    }
    * produce modified identifiers _i
    * _i is unit identifiers that are positive integers if i is ever treated in sample and 0 otherwise
    qui gen `_i' = `i' if  `_t' > 0 & `_t' < . & `ifstatement'
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
    if "`gtoolscheck'" == "g" qui gegen `refperiodobserved' = max(`refperiodobserved_temp') if `ifstatement' & `ttre' != ., by(`i')
    if "`gtoolscheck'" != "g" qui bysort `i': egen `refperiodobserved' = max(`refperiodobserved_temp') if `ifstatement' & `ttre' != .
    qui drop `refperiodobserved_temp'
    qui replace `ifvar' = 0 if `refperiodobserved' == 0  & `ttre' != .
    * update _t and ttre
    qui replace `_t' = . if `ifvar' == 0
    qui replace `ttre'  = . if `ifvar' == 0


    * if using time trends, you need to see a reference period block that is 2 periods in length; drop units where this is impossible
    if "`timetrends'" == "timetrends" {
      tempvar ttrefperiodobserved_temp
      tempvar ttrefperiodobserved
      qui gen `ttrefperiodobserved_temp' = `_t' == `refperiod' -1 if `ifstatement' & `ttre' != .
      if "`gtoolscheck'" == "g" qui gegen `ttrefperiodobserved' = max(`ttrefperiodobserved_temp') if `ifstatement' & `ttre' != ., by(`i')
      if "`gtoolscheck'" != "g" qui bysort `i': egen `ttrefperiodobserved' = max(`ttrefperiodobserved_temp') if `ifstatement' & `ttre' != .
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
      if "`gtoolscheck'" == "g" qui gegen `controlcheck2' = max(`controlcheck') if `ifstatement' , by(`t')
      if "`gtoolscheck'" != "g" qui bysort `t': egen `controlcheck2' = max(`controlcheck') if `ifstatement'
      qui replace `ifvar'  = 0 if `controlcheck2' == 0
      if "`safety'" != "off"  {
        cap assert `controlcheck2' != 0
        if _rc != 0 {
          dis as text "Estimability checks found time periods exist in sample without a control unit, causing them to be dropped."
          if "`contreat'" != "" dis as text "If intending to only use within treated-cohorts variation in your continuous treatment variable, specify option contreatwithin to avoid these checks."
        }
      }
      qui drop `controlcheck' `controlcheck2'
    }

    * update _t and ttre
    qui replace `_t' = . if `ifvar' == 0
    qui replace `ttre'  = . if `ifvar' == 0

    * grab information about range of _t
    qui `gtoolscheck2' sum `_t' if `ifstatement', meanonly
    local min_t = abs(r(min))
    local max_t = r(max) + 1


    * set _t into all positive terms to allow usge in regressions
    qui replace `_t' = `_t'+ `min_t' + 1

    *redefine the ref period macro for future use and then do a little bit of housekeeping
    local refperiod = `esrelativeto'+`min_t'+1
    qui drop `refperiodobserved'
    if "`timetrends'" == "timetrends" qui drop `ttrefperiodobserved'


    * check that there is a real sample in that someone, somewhere, has a pre-period and a post period
    if "`safety'" != "off"   {
      qui  `gtoolscheck2' sum `_t'  if `ifstatement' & `ttre' != ., meanonly
      local failure 0
      cap assert r(N) != 0
      if _rc != 0 local failure 1
      cap assert r(min) <= `refperiod'
      if _rc != 0 local failure 1
      cap assert r(max) > `refperiod'
      if _rc != 0 local failure 1
      if `failure' == 1 dis as error  "Error: cannot find sufficient pre- and post-treatment periods for the treated units."
      if `failure' == 1 & "`timetrends'" == "timetrends"  dis as error  "Note that inclusion of cohort-specific time trends increases the required number of pre-treatment units."
      if `failure' == 1 error 459
    }
    * send _t to 0  when we want it to be  a comparison period
    qui replace `_t' = 0 if (`_t' == . | `_t' == `refperiod' ) & `ifstatement'
    if "`timetrends'" == "timetrends" qui replace `_t' = 0 if `_t' == `refperiod' - 1 & `ifstatement'

    * produce modified identifiers
    * _i is unit identifiers that are positive integers if i is ever treated in sample and 0 otherwise
    qui gen `_i' = `i' if  `_t' > 0 & `_t' < . & `ifstatement'
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
      error 459
    }
  }

  * if number of requested horizons is too short, don't do adversarial pre-trend testing
  * the + 1 is to account for the 0 reference period
  if (`espre' + `espl'+ 1 < 3 & "`timetrends'" != "") | `espl' == 0 | `espre' == 0 local adversarial

  ****** after the above estimability checks and processing, does the SG variable have surviving members?
  * while we're at it, figure out how many levels the sg variable takes on during the treatment period - will be useful down the line

  if "`subgroups'" != "" {
    qui `gtoolscheck'levelsof(`subgroups') if `ifstatement' & `t' >= `ttre' , clean local(sglist)
    if wordcount("`sglist'") == 0 {
      dis as error  "Error: subgroups are specified, but contains no values within sample during the treatment period."
      error 459
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
          if "`gtoolscheck'" == "" _pctile `contreat' `latticeweightinserter' if `ifstatement' `xtileif', nquantiles(`latticenum')
          if "`gtoolscheck'" == "g" gquantiles `contreat' `latticeweightinserter' if `ifstatement' `xtileif', _pctile   nquantiles(`latticenum')
          gen `latticevar' = 1
          forval latticecase = 2(1)`latticenum'{
            local percentiletograb = `latticecase'-1
            qui replace `latticevar' = `latticevar' + 1 if `contreat' > r(r`percentiletograb')
          }
        }
        else {
          if "`latticetype'" != "all" & "`gtoolscheck'" == "g" qui gquantiles `latticevar' = `contreat' `latticeweightinserter' if `ifstatement' , xtile nquantiles(`latticenum') by(`byable')
          if "`latticetype'" != "all" & "`gtoolscheck'" == "" qui bysort `byable': egen `latticevar' = xtile(`contreat'  ) `latticeweightinserter' if `ifstatement' , nquantiles(`latticenum')
          if "`latticetype'" == "all" & "`gtoolscheck'" == ""  qui egen `latticevar' = xtile(`contreat') `latticeweightinserter' if `ifstatement' ,   nquantiles(`latticenum')
          if "`latticetype'" == "all" & "`gtoolscheck'" == "g"  qui gquantiles `latticevar' = `contreat' `latticeweightinserter' if `ifstatement' ,   nquantiles(`latticenum') xtile
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
    tempvar ringer
    qui gen `ringer' = `contreat'
    local contreatpolypart c.`contreat'
    local ringerpoly c.(`ringer'
      if `contreatpoly' > 1 {
        forval polycount = 2(1)`contreatpoly' {
          local contreatpolypart `contreatpolypart'##c.`contreat'
          tempvar ringer`polycount'
          qui gen `ringer`polycount'' = `ringer'^`polycount'
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
      qui gen `ic`iter'' = `intercon'
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



  *****************************************************************************************************************
  ******************************  PRODUCE THE WEIGHTS FOR CREATING VARIOUS OUTCOMES  ******************************
  *****************************************************************************************************************
  * special weights are required to go beyond a simple ATT; this generates them

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
    qui gen `posttreatment' = `t' >= `ttre' if `ifstatement'
    qui replace `posttreatment' = 2 if `t' < `ttre' & `ttre' != . & `_i' != 0 & `_t' != 0 & `ifstatement'
    qui replace `posttreatment' = 3 if `ttre' == .  & `ifstatement'
    qui replace `posttreatment' = 4 if `posttreatment' == .
  }

  * this variable helps count the number of i-t cells that exist in various subsets - useful for att_itime weights
  if "`att_itime'" == "att_itime" {
    qui `gtoolscheck'egen `tag' = tag(`i' `t') if  `ifstatement'
    qui replace `tag' = 0 if `tag' == .
  }

  * generate some inputs for weights calculation
  * w_per_it contains the amount of weight falling into each i-t cell
  * w per i gives the total amount of sample weight associated for each i falling into the categories specified by posttreatment
  * t per i is the total number of periods appearing for i in each posttreatment-defined period
  * in each case, we only make the input if needed for a weight type the user actually requested
  * that said, we always produce w_att_it whenever a non att weight is requested, since they are needed for the event horizon estimates and sometimes for histogram
  if "`gtoolscheck'" != "g" {
    if "`att_it'" == "att_it" | "`att_itime'" == "att_itime"  | "`att_i'" == "att_i" qui bysort `i' `t': egen `w_per_it' = total(`ifvar' * `w')
    if "`att_i'" == "att_i"  qui bysort  `i' `posttreatment': egen `w_per_i' = total(`w'*`ifvar')
    if "`att_itime'" == "att_itime" qui bysort  `i' `posttreatment': egen `t_per_i' = total(`tag')
  }
  else{
    if "`att_it'" == "att_it" | "`att_itime'" == "att_itime" | "`att_i'" == "att_i" qui gegen `w_per_it' = total(`ifvar' * `w')  , by(`i' `t')
    if "`att_i'" == "att_i"  qui gegen `w_per_i' = total(`w'*`ifvar')      , by(`posttreatment' `i')
    if "`att_itime'" == "att_itime" qui gegen `t_per_i' = total(`tag')    , by(`posttreatment' `i')
  }

  * att requires no special weights - it gets to stay as the given weight (w), albeit with the customweight text multiplier coming through
  * dividing w by w_per_it weights all i-t cells equally, but still lets w be influential within i-t cells (att_it)
  * dividing w by total w in i weights all i equally, but still lets w be influential within i-t cells and across i-t cells within i (att_i)
  * dividing w by total w in i-t cells * number of t in the i, yields w_att_itime. this lets w be influential inside i-t cells, but then
  *        weights every i-t cell equally when contributing to estimates of i-specific effects, which are themselves then averaged

  if "`att_i'" == "att_i" {
    qui gen `w_att_i' = `w' `customweighttext'/`w_per_i'
    qui replace `w_att_i' = 0 if `w_att_i' == .
    qui drop `w_per_i'
    qui drop  `posttreatment'
  }

  if "`att_itime'" == "att_itime" {
    qui gen `w_att_itime' = `w' `customweighttext' /(`w_per_it'*`t_per_i')
    qui replace `w_att_itime' = 0 if `w_att_itime' == .
    qui drop `tag' `t_per_i'
    qui cap drop  `posttreatment'
  }

  if "`att_it'" == "att_it" | "`att_itime'" == "att_itime"  | "`att_i'" == "att_i" {
    qui gen `w_att_it' = `w' `customweighttext'/`w_per_it'
    qui replace `w_att_it' = 0 if `w_att_it' == .
    qui drop `w_per_it'
  }

  * generate w_att whenever att is requested
  if "`att'" == "att" qui gen `w_att' = `w' `customweighttext'

  * some comments on the weights:
  * provided each i is always wholly contained in a given sg, we have no issues using these weights for sg purposes
  * since we don't allow fw type weights, weight scale here actually doesn't matter - no need to multiply weights by scalar adjustment factors
  * for event study estimates, note that att_i and att_itime are not conceptually all that sensible --
  *       we substitute att_it for them since it retains a spirit of 'all same across i'



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
  if "`verbose'" == "verbose" dis as text "Wooldid has finished regression estimation and is preparing to begin margins calculations."

	* check on dropped obs
	qui count if `dropcount' & e(sample) == 0
	local droppedobs = r(N)

  * grab some items from the model - sample indicator variable, residual degrees of freedom, etc
  tempvar esamp
  qui gen byte `esamp' = e(sample)
  local df_r = e(df_r)

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
    tokenize `bvars'
    local maxwordcount = wordcount("`bvars'")
    local catastrophe 0
    local bigcatastrophe 0
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

          * are there observations falling in the _i _t cell associated with the missing coefficient? if not, this isn't a problem
          * note: this regrettably is a bit computationally intensive
          qui count if `_i' == `problem_i' & `_t' == `problem_t' & `targetflag' == 1 & `w' > 0 & `ifstatement'
          if r(N) > 0 {
            * if we proceeded to this point, then the omitted object has associated observations, so this is an actual problem!
            * 2 remaining cases. wordcount = 3, in which case what we have is probably a dropped i-t dummy
            * or wordcount > 3, in which case we have something else, like a contreat interaction

            local wordcountproblem = wordcount("`checkifpresent'")

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
    foreach sumvar in `sumvarlist' {
      local sumvar2 c
      if "`sumvar'" == "y" local sumvar2 y
      qui `gtoolscheck2' sum ``sumvar'' `weightinserter' if e(sample), detail
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
        qui `gtoolscheck2' sum ``sumvar'' `weightinserter' if e(sample) & `ttre' <. & `t' < `ttre', detail
        local `sumvar'pre_mean = r(mean)
        local `sumvar'pre_med = r(p50)
        local `sumvar'pre_sd = r(sd)
        local `sumvar'pre_p25 = r(p25)
        local `sumvar'pre_p75 = r(p75)
        local `sumvar'pre_iqr = ``sumvar'_p75' - ``sumvar'_p25'
        local `sumvar'pre_min = r(min)
        local `sumvar'pre_max = r(max)
        local rownamessum `rownamessum' `sumvar2'_pretreat
      }

      * do the same in the post-treatment treated sample
      qui `gtoolscheck2' sum ``sumvar'' `weightinserter' if e(sample) & `ttre' <. & `t' >= `ttre', detail
      local `sumvar'post_mean = r(mean)
      local `sumvar'post_med = r(p50)
      local `sumvar'post_sd = r(sd)
      local `sumvar'post_p25 = r(p25)
      local `sumvar'post_p75 = r(p75)
      local `sumvar'post_iqr = ``sumvar'_p75' - ``sumvar'_p25'
      local `sumvar'post_min = r(min)
      local `sumvar'post_max = r(max)
      local rownamessum `rownamessum' `sumvar2'_posttreat

      * do the same in the control group
      if "`contreatwithin'" == "" & ("`sumvar'" == "y" | ("`contreatcontrols'" != "untreatedZero" & "`sumvar'" == "contreat"))  {
        qui `gtoolscheck2' sum ``sumvar'' `weightinserter' if e(sample) & `ttre' <. & `t' < `ttre', detail
        local `sumvar'con_mean = r(mean)
        local `sumvar'con_med = r(p50)
        local `sumvar'con_sd = r(sd)
        local `sumvar'con_p25 = r(p25)
        local `sumvar'con_p75 = r(p75)
        local `sumvar'con_iqr = ``sumvar'_p75' - ``sumvar'_p25'
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
    * graph window management
    cap graph drop *

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
      qui svmat `___hist'
      rename `___hist'1 att_i
      if "`contreat'" != "" rename `___hist'2 ame_i
      local bininserter
      if _N < 100 local bininserter bin(20)
      qui hist att_i, freq   xtitle("Cohort-Specific Average Treatment `effecttype'") name(hist_att_i) gap(5) `bininserter' color(navy)
      if "`saveplots'" != "" qui graph save "`saveplots'_hist_att_i.gph", `replace'
      if "`contreat'" != "" qui hist ame_i, freq   xtitle("Cohort-Specific Average Marginal Treatment `Ceffecttype'") name(hist_ame_i) gap(5) `bininserter' color(navy)
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
      qui svmat `___hist2'
      rename `___hist2'1 att_itime
      if "`contreat'" != "" rename `___hist2'2 ame_itime
      local bininserter
      if _N < 100 local bininserter bin(20)
      qui hist att_itime, freq   xtitle("Cohort-Specific i-t Cell Average Treatment `effecttype'") name(hist_att_itime) gap(5) `bininserter' color(navy)
      if "`saveplots'" != "" qui graph save "`saveplots'_hist_att_itime.gph", `replace'
      if "`contreat'" != "" qui hist ame_itime, freq   xtitle("Cohort-Specific i-t Cell Average Marginal Treatment `Ceffecttype'") name(hist_ame_itime) gap(5) `bininserter' color(navy)
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

  * if doing adversarial checks, inject the "post" statement where needed
  local postinject
  if "`adversarial'" == "adversarial" local postinject post

  * initialize some times
  local joint_horiz_i
  local joint_Choriz_i
  local pretlist

  * do not calculate traditional SEs if producing Ibragimov and Muller style p-values
  if "`impvalues'" != "" local marginsdof nose


  * loop through requested estimate types - basically just different weights
  foreach item in `att' `att_it' `att_i' `att_itime' {
    if "`verbose'" == "verbose" dis as text "Wooldid is now computing `item' estimates for full sample, including event study estimates."

    tempname ___`item'
    * use margins to produce the post treatment average effect within the treated group. captured via targetflag = 1 & t >= treatment period
    * do this for the non-impvalues case
    if "`impvalues'" == "" {
      qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' `marginsdof'
      matrix `___`item'' = r(table)[1..6,2]'
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
      qui svmat `ibmucoefs'
      qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
      local b = r(mean)
      qui ttest `ibmucoefs'1 == 0
      local se = r(se)
      local ll = `b' + invt(r(df_t),.025)*`se'
      local ul = `b' - invt(r(df_t),.025)*`se'
      matrix `___`item'' = (`b',`se', r(t), r(p), `ll', `ul')
      frame change `defaultframe'
      frame drop `ibmuframe'
      matrix drop `ibmucoefs'
    }

    * average pre-treatment effect
    if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
      tempname ___`item'_pre

      if "`impvalues'" == "" {
        qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' `marginsdof'
        matrix `___`item'_pre' = r(table)[1..6,2]'
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
        qui svmat `ibmucoefs'
        qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
        local b = r(mean)
        qui ttest `ibmucoefs'1 == 0
        local se = r(se)
        local ll = `b' + invt(r(df_t),.025)*`se'
        local ul = `b' - invt(r(df_t),.025)*`se'
        matrix `___`item'_pre' = (`b',`se', r(t), r(p), `ll', `ul')
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
          matrix `___`horizitem'_horiz' = ( r(table)[1..6, `matrixstartgrab'..`=colsof(r(b))']',`relyrcolumn')
          if "`postinject'" == "post" qui est sto __marginsmodel

          * this implements the f test that all points are well explained by a single line fit on the pre-period, constrained to run through 0
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
            qui svmat `pretrendregmat'
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
              qui `gtoolscheck2' sum t, meanonly
              qui replace b = 0 if b == .
              qui replace t = r(min) - 1 if _n == _N
            }
            else{
              qui replace b = 0 if b == .
              qui gen missingfind = t - t[_n -1]
              qui `gtoolscheck'levelsof(t) if missingfind == `bignjump' + 1, clean local(skiplocation)
              qui drop missingfind
              if "`skiplocation'" != "" {
                qui replace t = `skiplocation' - 1 if _n == _N
                if "`timetrends'" == "timetrends" qui replace t = `skiplocation' - 2 if _n == _N - 1
              }
              else{
                qui `gtoolscheck2' sum t, meanonly
                qui replace t = r(min) -1 if _n == _N
                if "`timetrends'" == "timetrends" qui replace t = r(min) - 2 if _n == _N - 1
              }
            }

            qui reg b t
            tempname adversarialB_`horizitem'
            matrix `adversarialB_`horizitem''   = e(b)
            qui predict linearproject, xb
            qui `gtoolscheck'levelsof(t) if _n <= _N - `bignjump', clean local(adversarialhlist)
            local testconstructor
            forval advframe_i = 1(1)`bignorig' {
              local horiztocheck = t[`advframe_i']
              local xbline = linearproject[`advframe_i']
              local testconstructor `testconstructor' ([1.`targetflag']`horiztocheck'.`_t' = `xbline')
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
            qui svmat `ibmucoefs'
            qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
            local b = r(mean)
            qui ttest `ibmucoefs'1 == 0
            local se = r(se)
            local ll = `b' + invt(r(df_t),.025)*`se'
            local ul = `b' - invt(r(df_t),.025)*`se'
            if `currenthoriz' == 1 matrix `horizhold' = (`b',`se', r(t), r(p), `ll', `ul')
            if `currenthoriz' >1  matrix `horizhold' = `horizhold' \ (`b',`se', r(t), r(p), `ll', `ul')
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
      qui margins [`weight'=`w_att'] ,`semielasticity'(`targetflag') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . & (`_t' == 0 | `_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' df(`df_r')
      local joint_horiz = r(p)[1,2]
    }
    if "`joint_horiz_i'" == "" & (  "`att_i'" != "" | "`att_itime'" != "" | "`att_it'" != "" ){
      qui margins [`weight'=`w_att_it'],`semielasticity'(`targetflag') noestimcheck over(`_t') subpop(if `ifstatement' & `t' < `ttre' & `ttre' != . &  (`_t' == 0 | `_t' >= `horizonstarter')   ) contrast(overjoint) `poismarginstext' df(`df_r')
      local joint_horiz_i = r(p)[1,2]
    }
  }


  if "`verbose'" == "verbose" dis as text "Wooldid has completed production of all main average treatment effects for the full sample."









  **********************************************************************************************************************************************
  *****************************         MAKE CONTINUOUS TREATMENT EFFECTS PLOTS AND CALCULATE CONTINUOUS TREATMENT EFFECTS      ****************************
  **********************************************************************************************************************************************
  * graph window management
  if "`makecfxplots'" != "" & "`histogramcohorteffects'" == "" cap graph drop *
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
        matrix `___C`item'' = r(table)[1..6,1]'
      }

      if "`impvalues'" == "impvalues" {
        qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck `poismarginstext' nose over(`cluster')
        tempname ibmucoefs
        tempname ibmuframe
        matrix `ibmucoefs' =  r(b)[1, 1..`=colsof(r(b))']'
        frame create `ibmuframe'
        frame change `ibmuframe'
        qui svmat `ibmucoefs'
        qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
        local b = r(mean)
        qui ttest `ibmucoefs'1 == 0
        local se = r(se)
        local ll = `b' + invt(r(df_t),.025)*`se'
        local ul = `b' - invt(r(df_t),.025)*`se'
        matrix `___C`item'' = (`b',`se', r(t), r(p), `ll', `ul')
        frame change `defaultframe'
        frame drop `ibmuframe'
        matrix drop `ibmucoefs'
      }

      * pretreatment average marginal effect of the continuous treatment variable
      if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
        tempname ___C`item'_pre

        if "`impvalues'" == "" {
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' `marginsdof'
          matrix `___C`item'_pre' = r(table)[1..6,1]'
        }

        if "`impvalues'" == "impvalues" {
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck `poismarginstext' nose over(`cluster')
          tempname ibmucoefs
          tempname ibmuframe
          matrix `ibmucoefs' =  r(b)[1, 1..`=colsof(r(b))']'
          frame create `ibmuframe'
          frame change `ibmuframe'
          qui svmat `ibmucoefs'
          qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
          local b = r(mean)
          qui ttest `ibmucoefs'1 == 0
          local se = r(se)
          local ll = `b' + invt(r(df_t),.025)*`se'
          local ul = `b' - invt(r(df_t),.025)*`se'
          matrix `___C`item'_pre' = (`b',`se', r(t), r(p), `ll', `ul')
          frame change `defaultframe'
          frame drop `ibmuframe'
          matrix drop `ibmucoefs'
        }
      }



      ************** MAKE CFX PLOTS - effects of treatment (ATT type and AME type) across the contreat distribution (across the lattice)
      if "`makecfxplots'" != "" & ("`item'" == "att" | "`safety'" == "off") {

        local cfxmatendpoint 6
        if "`cfxplotnose'" != "" local cfxmatendpoint 1
        if "`cfxplotnose'" == "" local cfxplotnose `marginsdof'
        if "`cfxplotnose'" == "cfxplotnose" local cfxplotnose nose

        * specify title for cfx plots x axis
        if "`cfxlattice'" != ""   local cfx_xtitle xtitle("estimation lattice: `cfxlattice'")
        if "`cfxlattice'" == "" & wordcount("`lattice'")==2 {
          local cfx_xtitle xtitle("`latticenum' quantiles (type: `latticetype')")
        }

        * use the lattice var for the plots if cfxlattice is not specified
        if "`cfxlattice'" == "" local cfxlattice `latticevar'

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
          qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlattice') `cfxplotnose' `poismarginstext'
          tempname cfx`item'_`nametext'
          matrix `cfx`item'_`nametext'' = r(table)[1..`cfxmatendpoint',floor(`=colsof(r(table))'/2)+1..`=colsof(r(table))']'
          qui marginsplot, recastci(rarea) plotopts(color(navy)) ciopts(color(navy%50)) `cfx_xtitle'  scale(.85) ytitle("`ytext'") title("`titletext'") name("cfx`item'_`nametext'")
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
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlattice') `cfxplotnose' `poismarginstext'
          tempname cfx`itemtext'_`nametext'
          matrix `cfx`itemtext'_`nametext'' = r(table)[1..`cfxmatendpoint',1..`=colsof(r(table))']'
          qui marginsplot, recastci(rarea) ciopts(color(navy%50)) plotopts(color(navy)) `cfx_xtitle' ytitle("`ytext'")  scale(.85) title("`titletext'") name("cfx`itemtext'_`nametext'")
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

              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `subgroups' == `sgcase' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlattice') `cfxplotnose' `poismarginstext'
              tempname cfx`item'_`nametext'_sg`sgcase'
              matrix `cfx`item'_`nametext'_sg`sgcase''  = r(table)[1..`cfxmatendpoint',floor(`=colsof(r(b))'/2)+1..`=colsof(r(b))']'
              qui marginsplot, recastci(rarea) plotopts(color(navy)) ciopts(color(navy%50)) `cfx_xtitle' ytitle("`ytext'") title("`titletext'") scale(.85) name("cfx`item'_`nametext'_sg`sgcase'")
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
              qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `subgroups' == `sgcase' & `targetflag'==1  `cfxif') noestimcheck over(`cfxlattice') `cfxplotnose' `poismarginstext'
              tempname cfx`itemtext'_`nametext'_sg`sgcase'
              matrix `cfx`itemtext'_`nametext'_sg`sgcase'' = r(table)[1..`cfxmatendpoint',1..`=colsof(r(b))']'
              qui marginsplot, plotopts(color(navy)) recastci(rarea) ciopts(color(navy%50)) `cfx_xtitle' ytitle("`ytext'") title("`titletext'")  scale(.85) name("cfx`itemtext'_`nametext'_sg`sgcase'")
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
            matrix `___C`horizitem'_horiz' = (r(table)[1..6, 1..`=colsof(r(b))']',`relyrcolumn')
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
              qui svmat `pretrendregmat'
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
                qui `gtoolscheck2' sum t, meanonly
                qui replace b = 0 if b == .
                qui replace t = r(min) - 1 if _n == _N
              }
              else{
                qui replace b = 0 if b == .
                qui gen missingfind = t - t[_n -1]
                qui `gtoolscheck'levelsof(t) if missingfind == `bignjump' + 1, clean local(skiplocation)
                qui drop missingfind
                if "`skiplocation'" != "" {
                  qui replace t = `skiplocation' - 1 if _n == _N
                  if "`timetrends'" == "timetrends" qui replace t = `skiplocation' - 2 if _n == _N - 1
                }
                else{
                  qui `gtoolscheck2' sum t, meanonly
                  qui replace t = r(min) -1 if _n == _N
                  if "`timetrends'" == "timetrends" qui replace t = r(min) - 2 if _n == _N - 1
                }
              }

              qui reg b t
              tempname adversarialB_C`horizitem'
              matrix `adversarialB_C`horizitem'' = e(b)
              qui predict linearproject, xb
              qui `gtoolscheck'levelsof(t) if _n <= _N - `bignjump', clean local(adversarialhlist)
              local testconstructor
              forval advframe_i = 1(1)`bignorig' {
                local horiztocheck = t[`advframe_i']
                local xbline = linearproject[`advframe_i']
                local testconstructor `testconstructor' (`horiztocheck'.`_t' = `xbline')
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
              qui svmat `ibmucoefs'
              qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
              local b = r(mean)
              qui ttest `ibmucoefs'1 == 0
              local se = r(se)
              local ll = `b' + invt(r(df_t),.025)*`se'
              local ul = `b' - invt(r(df_t),.025)*`se'
              if `currenthoriz' == 1 matrix `horizhold' = (`b',`se', r(t), r(p), `ll', `ul')
              if `currenthoriz' >1  matrix `horizhold' = `horizhold' \ (`b',`se', r(t), r(p), `ll', `ul')
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
          matrix `___`item'_s' = r(table)[1..6,`matrixstartgrab'..`=colsof(r(b))']'
        }
        * this handles the impvalues, though there is reason for suspicion about these depending on how cohorts and clusters interrelate
        if "`impvalues'" == "impvalues" {
          local iter 1
          tempname sgimpmat
          foreach sgcase in `sglist' {
            qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre' & `subgroups' == `sgcase') noestimcheck `poismarginstext' nose over(`cluster')
            tempname ibmucoefs
            tempname ibmuframe
            local matrixstartgrab floor(`=colsof(r(b))'/2)+1
            matrix `ibmucoefs' =  r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
            frame create `ibmuframe'
            frame change `ibmuframe'
            qui svmat `ibmucoefs'
            qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
            local b = r(mean)
            qui ttest `ibmucoefs'1 == 0
            local se = r(se)
            local ll = `b' + invt(r(df_t),.025)*`se'
            local ul = `b' - invt(r(df_t),.025)*`se'
            if `iter'==1 matrix `sgimpmat' = (`b',`se', r(t), r(p), `ll', `ul')
            if `iter'>1 matrix `sgimpmat' = `sgimpmat' \ (`b',`se', r(t), r(p), `ll', `ul')
            frame change `defaultframe'
            frame drop `ibmuframe'
            local ++iter
            matrix drop `ibmucoefs'
          }
          matrix `___`item'_s' =   `sgimpmat'
          matrix drop `sgimpmat'
        }


        * pretreatment effects
        if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
          tempname ___`item'_pre_s
          local sgpre 1

          if "`impvalues'" == "" {
            qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0) noestimcheck over(`subgroups') `poismarginstext' `marginsdof'
            local matrixstartgrab floor(`=colsof(r(b))'/2)+1
            matrix `___`item'_pre_s' = r(table)[1..6,`matrixstartgrab'..`=colsof(r(b))']'
          }

          if "`impvalues'" == "impvalues" {
            local iter 1
            tempname sgimpmat
            foreach sgcase in `sglist' {
              qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `subgroups' == `sgcase') noestimcheck `poismarginstext' nose over(`cluster')
              tempname ibmucoefs
              tempname ibmuframe
              local matrixstartgrab floor(`=colsof(r(b))'/2)+1
              matrix `ibmucoefs' =  r(b)[1, `matrixstartgrab'..`=colsof(r(b))']'
              frame create `ibmuframe'
              frame change `ibmuframe'
              qui svmat `ibmucoefs'
              qui ttest `ibmucoefs'1 == 0
              qui `gtoolscheck2' sum `ibmucoefs'1, meanonly
              local b = r(mean)
              qui ttest `ibmucoefs'1 == 0
              local se = r(se)
              local ll = `b' + invt(r(df_t),.025)*`se'
              local ul = `b' - invt(r(df_t),.025)*`se'
              if `iter'==1 matrix `sgimpmat' = (`b',`se', r(t), r(p), `ll', `ul')
              if `iter'>1 matrix `sgimpmat' = `sgimpmat' \ (`b',`se', r(t), r(p), `ll', `ul')
              frame change `defaultframe'
              frame drop `ibmuframe'
              local ++iter
              matrix drop `ibmucoefs'
            }
            matrix `___`item'_pre_s' =  `sgimpmat'
            matrix drop `sgimpmat'
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
              matrix `___`horizitem'_horiz_s`sgcase'' =  r(table)[1..6, `matrixstartgrab'..`=colsof(r(b))']'
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
                qui  `gtoolscheck'levelsof(`_t')  if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' & `_t' != 0, clean local(names)
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
                  qui svmat `pretrendregmat'
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
                    qui `gtoolscheck2' sum t, meanonly
                    qui replace b = 0 if b == .
                    qui replace t = r(min) - 1 if _n == _N
                  }
                  else{
                    qui replace b = 0 if b == .
                    qui gen missingfind = t - t[_n -1]
                    qui `gtoolscheck'levelsof(t) if missingfind == `bignjump' + 1, clean local(skiplocation)
                    qui drop missingfind
                    if "`skiplocation'" != "" {
                      qui replace t = `skiplocation' - 1 if _n == _N
                      if "`timetrends'" == "timetrends" qui replace t = `skiplocation' - 2 if _n == _N - 1
                    }
                    else{
                      qui `gtoolscheck2' sum t, meanonly
                      qui replace t = r(min) -1 if _n == _N
                      if "`timetrends'" == "timetrends" qui replace t = r(min) - 2 if _n == _N - 1
                    }
                  }

                  qui reg b t
                  tempname adversarialB_`horizitem'_s`sgcase'
                  matrix `adversarialB_`horizitem'_s`sgcase'' = e(b)
                  qui predict linearproject, xb
                  qui `gtoolscheck'levelsof(t) if _n <= _N - `bignjump', clean local(adversarialhlist)
                  local testconstructor
                  forval advframe_i = 1(1)`bignorig' {
                    local horiztocheck = t[`advframe_i']
                    local xbline = linearproject[`advframe_i']
                    local testconstructor `testconstructor' ([1.`targetflag']`horiztocheck'.`_t' = `xbline')
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
          if "`joint_horiz_i_s`sgcase''" == "" & ("`att_it'" == "att_it" | "`att_itime'" == "att_itime'" | "`att_i'" == "att_i" ){
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
    if wordcount("`sglist'") > 1 {
      foreach item in `att' `att_it' `att_i' `att_itime' {
        * this is a joint test that all sgs have the same effect, not vs 0 since no 0 category is included
        qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof'
        local joint_`item'_sdiff = r(p)[1,2]
        * the above gets a p-value for the joint test of all sgs vs the sg baselevel

        * the below grabs all pairwise comparisons among the sgs
        qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof'
        tempname pwsgs_`item'
        matrix `pwsgs_`item''  =  r(table_vs)[1..6,floor(`=colsof(r(table_vs))'/2)+1..`=colsof(r(table_vs))']'


        * repeat for pre-treatment
        if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
          * this is a joint test that all sgs have the same effect, not vs 0 since no 0 category is included
          qui  margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof'
          cap local prejoint_`item'_sdiff = r(p)[1,2]
          if _rc != 0 {
            local prejoint_`item'_sdiff
          * the above gets a p-value for the joint test of all sgs vs the sg baselevel. the if _rc is to handle cases where the pre-ame only exists for 1 sg
          }
          else {
            * the below grabs all pairwise comparisons among the sgs
            qui margins [`weight'=`w_`item''] , `semielasticity'(`targetflag') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof'
            tempname prepwsgs_`item'
            matrix `prepwsgs_`item''  =  r(table_vs)[1..6,floor(`=colsof(r(table_vs))'/2)+1..`=colsof(r(table_vs))']'
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
          matrix `___C`item'_s' = r(table)[1..6,1..`=colsof(r(b))']'

          * pretreatment average marginal effect of the continuous treatment variable
          if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
            qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0  ) noestimcheck over(`subgroups') `poismarginstext' `marginsdof'
            tempname ___C`item'_pre_s
            matrix `___C`item'_pre_s' = r(table)[1..6,1..`=colsof(r(b))']'
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
                matrix `___`horizitem'_Choriz_s`sgcase'' = ( r(table)[1..6, 1..`=colsof(r(b))']',`relyrcolumn_s`sgcase'')
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
                    qui svmat `pretrendregmat'
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
                      qui `gtoolscheck2' sum t, meanonly
                      qui replace b = 0 if b == .
                      qui replace t = r(min) - 1 if _n == _N
                    }
                    else{
                      qui replace b = 0 if b == .
                      qui gen missingfind = t - t[_n -1]
                      qui `gtoolscheck'levelsof(t) if missingfind == `bignjump' + 1, clean local(skiplocation)
                      qui drop missingfind
                      if "`skiplocation'" != "" {
                        qui replace t = `skiplocation' - 1 if _n == _N
                        if "`timetrends'" == "timetrends" qui replace t = `skiplocation' - 2 if _n == _N - 1
                      }
                      else{
                        qui `gtoolscheck2' sum t, meanonly
                        qui replace t = r(min) -1 if _n == _N
                        if "`timetrends'" == "timetrends" qui replace t = r(min) - 2 if _n == _N - 1
                      }
                    }

                    qui reg b t
                    tempname adversarialB_C`horizitem'_s`sgcase'
                    matrix `adversarialB_C`horizitem'_s`sgcase'' = e(b)
                    qui predict linearproject, xb
                    qui `gtoolscheck'levelsof(t) if _n <= _N - `bignjump', clean local(adversarialhlist)
                    local testconstructor
                    forval advframe_i = 1(1)`bignorig' {
                      local horiztocheck = t[`advframe_i']
                      local xbline = linearproject[`advframe_i']
                      local testconstructor `testconstructor' (`horiztocheck'.`_t' = `xbline')
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

            if ("`att_it'" == "att_it" | "`att_itime'" == "att_itime'" | "`att_i'" == "att_i" ) & "`joint_Choriz_i_s`sgcase''" == "" {
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
      if wordcount("`sglist'") > 1 {

        foreach item in `att' `att_it' `att_i' `att_itime' {

          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstext' noestimcheck over(`subgroups') contrast(overjoint) `marginsdof'
          local joint_C`item'_sdiff = r(p)[1,1]
           * the above gets a p-value for the joint test of all sgs vs the sg baselevel

          * the below grabs all pairwise comparisons among the sgs
          qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' >= `ttre') `poismarginstext' noestimcheck over(`subgroups') pwcompare(effects) `marginsdof'
          tempname pwsgs_C`item'
          matrix `pwsgs_C`item''  = r(table_vs)[1..6,1..`=colsof(r(table_vs))']'

          * repeat for pre-treatment
          if (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) {
            * this is a joint test that all sgs have the same effect, not vs 0 since no 0 category is included
            qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups') contrast(overjoint)  `poismarginstext' `marginsdof'
            cap local prejoint_C`item'_sdiff = r(p)[1,1]
            if _rc != 0 {
              local prejoint_C`item'_sdiff
              * the above gets a p-value for the joint test of all sgs vs the sg baselevel
            }
            else{
              * the below grabs all pairwise comparisons among the sgs
              qui margins [`weight'=`w_`item''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1 & `t' < `ttre' & `_t' != 0 & `ttre' != .) noestimcheck over(`subgroups')  pwcompare(effects) `poismarginstext' `marginsdof'
              tempname prepwsgs_C`item'
              matrix `prepwsgs_C`item''  =  r(table_vs)[1..6,1..`=colsof(r(table_vs))']'
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
                local trurerow
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
          local sgnumhoriz_start = (`sgnum'-1)*(`horizoncount'-`horizonstarter'+1)+1
          local sgnumhoriz_end= `sgnum'*(`horizoncount'-`horizonstarter'+1)
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
        if ("`prejoint_`item'_sdiff'" !=  "") & ( (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == "")))) ){
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
          if ("`prejoint_C`item'_sdiff'" != "") & ( (`espre' > 0 & "`esfixedbaseperiod'" == "") | ("`esfixedbaseperiod'" == "esfixedbaseperiod" & (`esrelativeto' < -1 | (`min_t' > 2 | (`min_t' > 1 & "`timetrends'" == ""))))) {
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
  if _rc != 0 {
    cap ereturn post `eb' , esample(`esamp') depname(`y') obs(`N') dof(`df_r')
    if _rc != 0 ereturn post   , esample(`esamp') depname(`y') obs(`N') dof(`df_r')
  }


  * return the command line input and estimate_type
  ereturn local cmdline `"`0'"'
  local estimate_type "pooled_reference_period"
  if "`esfixedbaseperiod'" != "" local estimate_type "fixed_reference_period"
  ereturn local estimate_type "`estimate_type'"



  ********* begin actual result presentation *************
  * the if statements are here to ensure appropriate info is shown depending on if you chose reg vs reghdfe, poisson vs ppmlhdfe, clustered ses vs robust, etc.
  dis as result  " "
  dis as result  "Wooldid Estimation for Outcome: `y'; Standard Errors: `vcespecifier'"
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
    matlist `summarystats'
    ereturn matrix summarystats = `summarystats'

    * most summary stats will live in the summary stats matrix, but these may be worth returning on their own
    ereturn scalar ypre_mean = `ypre_mean'
    ereturn scalar ypre_med =  `ypre_med'
    ereturn scalar ypre_sd =  `ypre_sd'
    ereturn scalar ypre_iqr =  `ypre_iqr'

    if "`contreat'" != "" {
	    ereturn scalar contreatpre_mean = `contreatpre_mean'
	    ereturn scalar contreatpre_med =  `contreatpre_med'
	    ereturn scalar contreatpre_sd =  `contreatpre_sd'
	    ereturn scalar contreatpost_mean = `contreatpost_mean'
	    ereturn scalar contreatpost_med =  `contreatpost_med'
	    ereturn scalar contreatpost_sd =  `contreatpost_sd'
	    ereturn scalar contreat_mean = `contreat_mean'
	    ereturn scalar contreat_med =  `contreat_med'
	    ereturn scalar contreat_sd =  `contreat_sd'
    }

    if "`poisson'" == "poisson" dis as result  "Note: summary stats above are for untransformed Y, while coefficients are from Y=exp(regression equation) and have a log-linear type interpretation."
  }


  * return various scalars and summary stats
  local r2list r2 r2adj r2_within r2_withinadj F
  if "`poisson'" == "poisson" local r2list r2pseud chi2
  if "`regtype'" != "poisson" local r2list `r2list' rmse rss
  foreach scalar in N Nclus ll ll_0 `r2list'   {
    if "``scalar''" != "." ereturn scalar `scalar' = ``scalar''
  }
  ereturn scalar N_dropped = `droppedobs'


  ************* now return main results with matlist `wdidpooled'
  dis as result  " "
  local poistext
  if "`poisson'" == "poisson" local poistext " Poisson Regression"
  local elasinj
  if "`semielasticity'" == "eydx" local elasinj Treatment Effects Presented as Semi-elasticities
  if "`contreatelasticitytype'" != "dydx" & "`elasinj'" != "" local elasinj `elasinj',
  if "`contreatelasticitytype'" == "eyex" local elasinj `elasinj' Marginal Treatment Effects Presented as Elasticities
  if "`contreatelasticitytype'" == "eydx" local elasinj `elasinj' Marginal Treatment Effects Presented as Semi-elasticities (eydx)
  if "`contreatelasticitytype'" == "dyex" local elasinj `elasinj' Marginal Treatment Effects Presented as Semi-elasticities (dyex)



  dis as result  "Main`poistext' Results (Full Estimation Sample): `elasinj'"
  if "`customweightmultiplier'" != "" dis as result  "Note: all Estimates Shown Apply `customweightmultiplier' Weights Post-Estimation"
  matlist `wdidpooled'

  * return the joint  tests
  if "`suppressjointtests'" == "" {
    dis as result  " "
    if "`joint_horiz'" != "" {
      dis as result  "p-value for joint test that all pre-treatment att horizon estimates are 0: " %-6.4f `joint_horiz'
      ereturn scalar joint_horiz = `joint_horiz'
    }
    if "`joint_horiz_i'" != "" {
      dis as result  "p-value for joint test that all pre-treatment att_it horizon estimates are 0: " %-6.4f `joint_horiz_i'
      ereturn scalar joint_horiz_i = `joint_horiz_i'
    }
    if "`joint_Choriz'" != "" {
      dis as result  "p-value for joint test that all pre-treatment ame horizon estimates are 0: " %-6.4f `joint_Choriz'
      ereturn scalar joint_Choriz = `joint_Choriz'
    }
    if "`joint_Choriz_i'" != "" {
      dis as result  "p-value for joint test that all pre-treatment ame_it horizon estimates are 0: " %-6.4f `joint_Choriz_i'
      ereturn scalar joint_Choriz_i = `joint_Choriz_i'
    }
  }
  dis as result  " "

 if "`adversarial'" != "" {
   foreach scenario in _ _C {
     foreach type in att att_it {
       if "`adversarialP`scenario'`type''" != "" {
         local typetext `type'
         if "`scenario'" == "_C" local typetext = subinstr("`typetext'","att","ame",.)
         dis as result "p-value for joint test of null that all (non-reference) `typetext' coefficients are on a single line: " %-6.4f `adversarialP`scenario'`type''
         ereturn scalar adversarial_pval_`typetext' = `adversarialP`scenario'`type''
         ereturn matrix adversarial_line_`typetext' =   `adversarialB`scenario'`type''
       }
     }
   }
 }
  dis as result  " "

  * now return the sg results
  if "`subgroups'" != "" {
    if "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        dis as result  " "
        local poistext
        if "`poisson'" == "poisson" local poistext "Poisson Regression"
        dis as result  "`poistext'Results for Subgroup `subgroups' == `sgcase': `elasinj'"
        if "`customweightmultiplier'" != "" dis as result  "Note: all Estimates Shown Apply `customweightmultiplier' Weights Post-Estimation"
        matlist `wdid_sg`sgcase''
        dis as result  " "
        * return the joint  tests
        if "`suppressjointtests'" == "" {

          if "`joint_horiz_s`sgcase''" != "" {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment att horizon estimates are 0: " %-6.4f `joint_horiz_s`sgcase''
            ereturn local joint_horiz_s`sgcase' = `joint_horiz_s`sgcase''
          }
          if "`joint_horiz_i_s`sgcase''" != ""  {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment att_it horizon estimates are 0: " %-6.4f `joint_horiz_i_s`sgcase''
            ereturn local joint_horiz_i_s`sgcase' = `joint_horiz_i_s`sgcase''
          }


          if "`joint_Choriz_s`sgcase''" != "" {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment ame horizon estimates are 0: " %-6.4f `joint_Choriz_s`sgcase''
            ereturn local  joint_Choriz_s`sgcase' = `joint_Choriz_s`sgcase''
          }

          if "`joint_Choriz_i_s`sgcase''" != ""  {
            dis as result  "p-value for joint test that sg=`sgcase' pre-treatment ame_it horizon estimates are 0: " %-6.4f `joint_Choriz_i_s`sgcase''
            ereturn local  joint_Choriz_i_s`sgcase' =  `joint_Choriz_i_s`sgcase''
          }

          dis as result  " "
        }
        if "`adversarial'" != "" {
          foreach scenario in _ _C {
            foreach type in att att_it {
              if "`adversarialP`scenario'`type'_s`sgcase''" != "" {
                local typetext `type'
                if "`scenario'" == "_C" local typetext = subinstr("`typetext'","att","ame",.)
                dis as result "p-value for joint test of null that all (non-reference) `typetext' coefficients for sg=`sgcase' are on a single line: " %-6.4f `adversarialP`scenario'`type'_s`sgcase''
                ereturn scalar adversarial_pval_`typetext'_s`sgcase' = `adversarialP`scenario'`type'_s`sgcase''
                ereturn matrix adversarial_line_`typetext'_s`sgcase' =   `adversarialB`scenario'`type'_s`sgcase''

              }
            }
          }
        }
        dis as result " "
      }
    }
    * having return the main sg materials, return the pairwise comparisons
    if wordcount("`sglist'") > 1 {
      dis as result  "Pairwise Comparisons of Different Subgroups' Estimates: `elasinj'"
      matlist `pwsgs'
      dis as result " "
      ereturn matrix pwsgs = `pwsgs'
      foreach item in `att'  `att_it' `att_i' `att_itime' {
        if "`joint_`item'_sdiff'" != "" dis as result  "p-value for joint test that all SGs have the same `item': " %-6.4f `joint_`item'_sdiff'
        if "`joint_`item'_sdiff'" != "" ereturn scalar joint_`item'_allsgssame =   `joint_`item'_sdiff'
        if "`prejoint_`item'_sdiff'" != "" dis as result  "p-value for joint test that all SGs have the same pre-treatment `item': " %-6.4f `prejoint_`item'_sdiff'
        if "`prejoint_`item'_sdiff'" != "" ereturn scalar joint_pre_`item'_allsgssame =   `prejoint_`item'_sdiff'
      }
      dis as result  " "
      if "`contreat'" != "" {
        foreach item in `att'  `att_it' `att_i' `att_itime' {
          local itemtext = subinstr("`item'","att","ame",.)
          if "`joint_C`item'_sdiff'" != "" dis as result  "p-value for joint test that all SGs have the same `itemtext': " %-8.4f `joint_C`item'_sdiff'
          if "`joint_C`item'_sdiff'" != "" ereturn scalar joint_`itemtext'_allsgssame =   `joint_C`item'_sdiff'
          if "`prejoint_C`item'_sdiff'" != "" dis as result  "p-value for joint test that all SGs have the same pre-treatment `itemtext': " %-6.4f `prejoint_C`item'_sdiff'
          if "`prejoint_C`item'_sdiff'" != "" ereturn scalar joint_pre_`itemtext'_allsgssame =   `prejoint_C`item'_sdiff'

        }
        dis as result  " "
      }
    }
  }

  *return the histogram matrix
  if "`histogramcohorteffects'" != "" ereturn matrix histogramestimates = `___hist'




  ***************************************************************************************************************************
  ********************************************         PRODUCE EVENT STUDY PLOTS       **************************************
  ***************************************************************************************************************************
  * present event study plots if requested; a little complicated because of naming practices in wdidpooled
  * this is also somewhat awkward placement b/c I have this coded to use the wdidpooled matrix
  * this choice is vestigial to a prior program design and is kept here on the principle of not redoing what already works
  * though I should probably redo this someday

  if `espre' +`espl' > 0 & "`makeplots'" == "makeplots" {

    if "`histogramcohorteffects'" == "" & "`makecfxplots'" == "" cap graph drop *

    * this creates and enters a frame where we will produce graphs using wdidpooled and similar for sgs
    tempname __wooldidframe
    frame create `__wooldidframe'
    frame change `__wooldidframe'
    qui svmat `wdidpooled', names(col)
    qui gen rownamesvar = ""
    local names : roweq `wdidpooled'
    forvalues i=1/`: word count `names'' {
      qui replace rownamesvar=`"`: word `i' of `names''"' in `i'
    }
    qui drop if relyr == .
    qui `gtoolscheck'levelsof(rownamesvar), clean local(graphtypes)
    qui gen postperiod =  relyr >= 0
    keep   postperiod relyr  estimate lb_95ci ub_95ci rownamesvar

    * add the pre treatment observation / reference period
    if "`esfixedbaseperiod'" == "" {
      local newN = _N + 1
      qui `gtoolscheck2' sum relyr, meanonly
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
      qui   `gtoolscheck2' sum relyr, meanonly
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
    qui   `gtoolscheck2' sum relyr, meanonly
    if r(min) < 0 local xlineinserter xline(-0.5,   lpattern(dash) lcolor(black))

    * inject elasticity info
    if "`semielasticity'" == "eydx" local elasinj Treatment Effects Presented as Semi-elasticities
    if "`contreatelasticitytype'" != "dydx" & "`elasinj'" != "" local elasinj `elasinj',
    if "`contreatelasticitytype'" == "eyex" local elasinj `elasinj' Marginal Treatment Effects Presented as Elasticities
    if "`contreatelasticitytype'" == "eydx" local elasinj `elasinj' Marginal Treatment Effects Presented as Semi-elasticities (eydx)
    if "`contreatelasticitytype'" == "dyex" local elasinj `elasinj' Marginal Treatment Effects Presented as Semi-elasticities (dyex)


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
      twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr,color(navy%50)) if (rownamesvar =="reference" | rownamesvar=="`graph'")  , `xlineinserter' xtitle("Time Relative to Time of Initial Treatment") ytitle("``graphtext'_text' on `y' by Horizon") legend(off)   name(`graphtext'_fullmodel) `xlabelstarter'
      if "`saveplots'" != "" qui graph save "`saveplots'_`graph'.gph", `replace'
    }


    if "`saveplots'" != "" qui save "`saveplots'_ES.dta", `replace'
    clear

    * repeat by SG if appropriate
    if "`subgroups'" != "" & "`suppresssgprimaryeffects'" == "" {
      foreach sgcase in `sglist' {
        qui svmat `wdid_sg`sgcase'', names(col)
        qui gen rownamesvar = ""
        local names : roweq `wdid_sg`sgcase''
        forvalues i=1/`: word count `names'' {
          qui replace rownamesvar=`"`: word `i' of `names''"' in `i'
        }
        qui drop if relyr == .
        qui `gtoolscheck'levelsof(rownamesvar), clean local(graphtypes)
        qui gen postperiod =  relyr >= 0
        keep  rownamesvar postperiod relyr  estimate lb_95ci ub_95ci

        * add the pre treatment observation
        if "`esfixedbaseperiod'" == "" {
          local newN = _N + 1
          qui   `gtoolscheck2' sum relyr, meanonly
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
          qui   `gtoolscheck2' sum relyr, meanonly
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
          twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr,color(navy%50)) if (rownamesvar =="reference" | rownamesvar=="`graph'")  , `xlineinserter' xtitle("Time Relative to Time of Initial Treatment") ytitle("``graphtext'_text' on `y' by Horizon") legend(off)   name(`graphtext'_sg`sgcase') `xlabelstarter'   title("`subgroups' == `sgcase'")
          if "`saveplots'" != "" qui graph save "`saveplots'_`graph'_sg`sgcase'.gph", `replace'
        }

        if "`saveplots'" != ""  qui save "`saveplots'_ES_`subgroups'`sgcase'.dta", `replace'
        clear
      }
    }

    frame change `defaultframe'
    frame drop `__wooldidframe'
  }



  ***************************************************************************************************************************
  ***************************************         SHIP REMAINING MATRICES TO ERETURN AND EXIT       **************************************
  ***************************************************************************************************************************

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

  if "`cleanmatrices'" != "" {
    clear mata
    matrix drop _all
  }

  if "`verbose'" != "verbose" qui est drop ___wooldidfullmodel
  * clean up after ppmlhdfe
  if "`poisson'" == "poisson" & ("`unconditionalse'" != "" | "`customvce'" != "" ) qui cap drop _ppmlhdfe_d


  * program ends!
  dis as result  " "
  ereturn local cmd "wooldid"




  end
