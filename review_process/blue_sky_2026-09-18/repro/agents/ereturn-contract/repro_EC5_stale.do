* EC-5: which failure sites destroy the user's previous wooldid e(), and what the
* consumers do afterwards. Twenty sites; only three land in reghdfe's internal model.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
replace ttre = 11 if grp==4
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + 0.25*t + 1.2*post + rnormal()
gen long sg8 = 10000000 + mod(id,2)
gen double dose = runiform()*post
save _ec5.dta, replace

cap prog drop STATE
prog def STATE
  args tag rc
  local bk = .
  cap confirm matrix e(b)
  if _rc == 0 local bk = colsof(e(b))
  local ns = 0
  cap count if e(sample)
  if _rc == 0 local ns = r(N)
  dis "SITE|`tag'|rc=`rc'|e(cmd)=`e(cmd)'|K=`bk'|e(sample) N=`ns'|userest=" ///
      cond(_rc==0,"",".")
  cap est dir USERM
  dis "   user estimates survive? " (_rc==0)
  cap confirm matrix USERMAT
  dis "   user matrix survives?    " (_rc==0)
end
cap prog drop GOOD
prog def GOOD
  use _ec5.dta, clear
  qui wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2)
  cap est drop USERM
  qui est store USERM
  cap matrix drop USERMAT
  matrix USERMAT = J(2,2,7)
end
GOOD
STATE baseline 0
* --- the three that DESTROY the previous e() ---
GOOD
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(untreatedZero)
STATE F08_contreat_collinear `=_rc'
GOOD
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod subgroups(sg8)
STATE F10_nameoverflow `=_rc'
GOOD
cap noi wooldid y id t ttre, att ccc_absorb(id) esprelength(2) espostlength(2)
STATE F16_ccc_absorb `=_rc'
* --- a representative sample of the seventeen that do NOT ---
GOOD
cap noi wooldid y id t ttre, attt
STATE F01_badoption `=_rc'
GOOD
cap noi wooldid y id t ttre, att permute(shuffle) permreps(10)
STATE F12_permfloor `=_rc'
* --- what the consumers do in the reghdfe state ---
GOOD
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(untreatedZero)
dis "=== consumers after the r(322) ==="
dis "e(cmd) = `e(cmd)'   K = " colsof(e(b))
cap noi predict double yh
dis "predict          -> rc " _rc
cap drop yh
cap noi predict double xb1, xb
dis "predict, xb      -> rc " _rc
cap drop xb1
cap noi estimates store AFTER
dis "estimates store  -> rc " _rc
cap noi wooldid
dis "bare wooldid     -> rc " _rc "  (the replay-staleness message never fires: e(cmd) is not wooldid)"
