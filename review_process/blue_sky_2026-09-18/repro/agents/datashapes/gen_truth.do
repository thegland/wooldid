*==========================================================================
* gen_truth.do -- the DESIGN TRUTH of wooldid's four estimands, computed
*                 directly from a gen_core dataset.
*
* WHY THIS EXISTS.  With `local gysd 0' the gen_core DGP is exactly
*   y = a_i + g_t (+ b1*x1 + b2*x2) + truetau
* so a correctly specified wooldid run has NO sampling error left: every
* estimand must come back equal to its design value to machine precision.
* That converts "is this number right?" from a judgement call into an
* equality test.
*
* THE FOUR WEIGHTS (read off wooldid.ado:1594-1650 of the pinned build):
*   att       w
*   att_it    w / (total w in the i-t cell)            -> every i-t cell equal
*   att_i     w / (total w in i, within post/pre)      -> every i equal
*   att_itime w / (cell w * #post cells of i)          -> cells equal in i, i equal
* and the post-treatment margins subpop is  targetflag==1 & t >= ttre.
*
* USAGE:
*     wooldid ync unit t ttre, att att_i att_it att_itime cluster(unit)
*     qui gen byte insmp = e(sample)
*     local TRUTHSMP insmp      // REQUIRED: wooldid's realized sample
*     local TRUTHI   unit       // the variable you passed as i
*     local TRUTHW   w          // the weight variable (or leave empty for none)
*     include "`GENDIR'/gen_truth.do"
*     di `TRUE_att'  `TRUE_att_i'  `TRUE_att_it'  `TRUE_att_itime'
*
* CAVEATS -- read these before filing anything:
*  1. Pass TRUTHSMP = e(sample).  wooldid drops units with no pre-period,
*     periods with no control unit, singletons, and rows with a zero or
*     missing weight; a truth computed on the raw data will not match and
*     the mismatch is NOT a defect.
*  2. Valid for the BINARY treatment path with the default pooled
*     reference period.  Under esfixedbaseperiod, anticipation(),
*     oostreatedcontrols() or contreat() the post set changes; adjust
*     `__trpost' below to match.
*  3. Rows whose only pre-period observation is their first period are
*     forced into the reference pool but stay POST rows if t >= ttre --
*     that is already handled, since post membership is t >= ttre.
*  4. Requested-horizon truncation: if you pass espost(k) with k smaller
*     than the longest realized horizon, the event-study rows are
*     truncated but the att family is NOT (its subpop is t >= ttre).
*==========================================================================

if "`TRUTHSMP'" == "" local TRUTHSMP 1
if "`TRUTHI'"   == "" local TRUTHI   unit
if "`TRUTHW'"   == "" local TRUTHW   1
if "`TRUTHTAU'" == "" local TRUTHTAU truetau

cap drop __trsmp __trpost __trw __trwit __trwip __trtag __trtpi
cap drop __trw_att __trw_it __trw_i __trw_itime

qui gen byte   __trsmp = (`TRUTHSMP') == 1
qui gen byte   __trpost = __trsmp & ttre < . & t >= ttre
qui gen double __trw   = `TRUTHW'

qui egen double __trwit = total(cond(__trsmp, __trw, 0)), by(`TRUTHI' t)
qui egen double __trwip = total(cond(__trpost, __trw, 0)), by(`TRUTHI')
qui egen byte   __trtag = tag(`TRUTHI' t) if __trpost
qui egen double __trtpi = total(__trtag == 1), by(`TRUTHI')

qui gen double __trw_att   = __trw
qui gen double __trw_it    = __trw / __trwit
qui gen double __trw_i     = __trw / __trwip
qui gen double __trw_itime = __trw / (__trwit * __trtpi)

foreach _e in att it i itime {
  local _nm `_e'
  if "`_e'" == "it"    local _nm att_it
  if "`_e'" == "i"     local _nm att_i
  if "`_e'" == "itime" local _nm att_itime
  qui summarize `TRUTHTAU' [aw = __trw_`_e'] if __trpost == 1, meanonly
  local TRUE_`_nm' = r(mean)
  local TRUEN_`_nm' = r(N)
}

qui count if __trpost == 1
local TRUE_npost = r(N)
qui summarize __trtag if __trpost == 1, meanonly
local TRUE_ncell = r(sum)

di as text "design truth on " as result `TRUE_npost' as text " post rows / " ///
   as result `TRUE_ncell' as text " post cells:"
di as text "  att       = " as result %12.9f `TRUE_att'
di as text "  att_it    = " as result %12.9f `TRUE_att_it'
di as text "  att_i     = " as result %12.9f `TRUE_att_i'
di as text "  att_itime = " as result %12.9f `TRUE_att_itime'

drop __trsmp __trpost __trw __trwit __trwip __trtag __trtpi
drop __trw_att __trw_it __trw_i __trw_itime
