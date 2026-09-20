* R4_px01_sharp.do -- SELF-CONTAINED.
* (i)  PX-01 sharpening: with subgroups() the subgroup-by-horizon margins call keeps its
*      VCE even under fastnose (ado:1981). Does fastnose+permanalyticse then print a
*      PARTIALLY filled "analytic standard error" column?
* (ii) known_bydesign section 7 item 11: does plain fastnose post e(b_es_att)/e(numpre_att)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103

set obs 60
gen int firm = _n
gen double fe_i = rnormal(0, 1.4)
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
gen double size = 3 + rnormal(0, .8)
gen byte sg = mod(_n, 2)
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
gen double y = 2.0 + fe_i + fe_t + 0.4*size + effect + eps
drop rel effect

cap prog drop rpt
prog def rpt
  args tag
  dis "`tag'_ecmd=" "`e(cmd)'"
  dis "`tag'_permase=" cond("`e(perm_analyticse)'"=="", "ABSENT", "`e(perm_analyticse)'")
  foreach m in V V_es_att b_es_att b numpre_att wdidmainresults wdidsgresults {
    cap confirm matrix e(`m')
    local h1 = (_rc==0)
    cap confirm scalar e(`m')
    local h2 = (_rc==0)
    dis "`tag'_has_`m'=" (`h1' | `h2')
  }
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    matrix list WR
  }
  cap confirm matrix e(wdidsgresults)
  if _rc == 0 {
    matrix SG = e(wdidsgresults)
    matrix list SG
  }
end

*========= (i) fastnose + permanalyticse + subgroups()
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   subgroups(sg) permute(shuffle) permreps(30) permanalyticse fastnose rseed(4242)
dis "S1_RC=" _rc
rpt S1

*========= (i-b) the honest comparator: permanalyticse + subgroups(), no fastnose
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   subgroups(sg) permute(shuffle) permreps(30) permanalyticse rseed(4242)
dis "S2_RC=" _rc
rpt S2

*========= (ii) section-7 item 11: plain fastnose, no permute
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) fastnose
dis "S3_RC=" _rc
rpt S3
dis "S3_numpre_att=" e(numpre_att)
ereturn list
dis "---- R4 done ----"
