* SKEPTIC probe of EC-3: is the e() left behind by the 8-digit-subgroup r(198)
* bit-for-bit the same as the e() of the identical run with short labels?
* Also: the contreat() destroying site EC-5 claims (rc 322, e(cmd)=reghdfe).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 8181
set sortseed 8181
which wooldid

clear
set obs 96
gen long id = _n
gen byte clu = ceil(id/4)
gen byte sgs = mod(id,2)
gen double sg8 = 10000000 + mod(id,2)
gen double ui = rnormal()
gen double dose = 0.5 + runiform()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2007 if mod(id,4)==1
replace ttre = 2010 if mod(id,4)==2
expand 13
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 1.5 + 0.4*(t-ttre) if t>=ttre & ttre<.
save p6base, replace

local SPEC "att esfixedbaseperiod esprelength(3) espostlength(3) jointtests adversarial cluster(clu)"

cap noisily wooldid y id t ttre, `SPEC' subgroups(sgs)
dis "P6 SHORT_RC=" _rc
mat bS = e(b)
mat VS = e(V)
local nS : colnames bS
dis "P6 SHORT names: `nS'"
* capture the full e() name set and scalar values
qui ereturn list, all
local eS `"`:  e(scalars)'"'
dis "P6 SHORT scalar count=" wordcount("`eS'")
foreach s of local eS {
  scalar SH_`=substr("`s'",1,26)' = e(`s')
}

cap noisily wooldid y id t ttre, `SPEC' subgroups(sg8)
dis "P6 LONG_RC=" _rc
dis "P6 LONG ecmd=|" e(cmd) "|"
mat bL = e(b)
mat VL = e(V)
local nL : colnames bL
dis "P6 LONG names: `nL'"

* --- bit comparison of the shared coefficient names ---
mata:
  bs = st_matrix("bS"); bl = st_matrix("bL")
  ns = st_matrixcolstripe("bS"); nl = st_matrixcolstripe("bL")
  printf("P6MATA cols short=%f long=%f\n", cols(bs), cols(bl))
  nd = 0
  for (a=1; a<=cols(bs); a++) {
    nm = ns[a,2]
    hit = 0
    for (bq=1; bq<=cols(bl); bq++) {
      if (nl[bq,2] == nm) {
        hit = 1
        if (bs[1,a] != bl[1,bq]) {
          printf("P6MATA DIFF %s  %21x vs %21x\n", nm, bs[1,a], bl[1,bq])
          nd = nd + 1
        }
      }
    }
    if (!hit) printf("P6MATA ONLY-IN-SHORT %s\n", nm)
  }
  printf("P6MATA numeric differences among shared names = %f\n", nd)
end
* the subgroup-keyed names differ only in the label, so line them up positionally too
mata:
  bs = st_matrix("bS"); bl = st_matrix("bL")
  if (cols(bs) == cols(bl)) {
    d = 0
    for (a=1; a<=cols(bs); a++) if (bs[1,a] != bl[1,a]) d = d + 1
    printf("P6MATA positional differences = %f of %f\n", d, cols(bs))
  }
  else printf("P6MATA column counts differ\n")
end

* --- how complete is the surviving e()? ---
qui ereturn list, all
local eL `"`: e(scalars)'"'
dis "P6 LONG scalar count=" wordcount("`eL'")
local eLm `"`: e(macros)'"'
dis "P6 LONG macro count=" wordcount("`eLm'")
local eLx `"`: e(matrices)'"'
dis "P6 LONG matrix count=" wordcount("`eLx'")
cap noisily lincom main_att
dis "P6 LONG lincom rc=" _rc " b=" r(estimate) " se=" r(se)
cap noisily ereturn display
dis "P6 LONG ereturn display rc=" _rc
dis "P6 LASTRC=|${WOOLDID_LASTRC}|"
cap noisily wooldid
dis "P6 bare replay rc=" _rc

* --- is the message TRUE?  compare the surviving e(cmdline) to the failed call ---
dis "P6 e(cmdline)=|" e(cmdline) "|"

*=========== EC-5's contreat destroying site, done properly ===========
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use p6base, clear
cap noisily wooldid y id t ttre, att att_i esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
dis "P6C base rc=" _rc " ecmd=|" e(cmd) "| K=" colsof(e(b)) " N=" e(N)
cap estimates store USERM2
matrix umat = J(3,3,1)
foreach CC in "pooled" "byCohort" "none" {
  cap noisily wooldid y id t ttre, att contreat(dose) contreatcontrols(`CC') esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
  dis "P6C contreat(`CC') rc=" _rc " ecmd=|" e(cmd) "| K=" colsof(e(b)) " N=" e(N)
  cap count if e(sample)
  dis "P6C   e(sample) N=" r(N)
  cap estimates describe USERM2
  dis "P6C   USERM2 rc=" _rc
  cap confirm matrix umat
  dis "P6C   umat rc=" _rc
  cap noisily predict double phat
  dis "P6C   predict rc=" _rc
  cap drop phat
  macro drop WOOLDID_LASTRC
  * re-establish the good baseline for the next site
  cap noisily wooldid y id t ttre, att att_i esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
}
dis "ALLDONE_P6"
